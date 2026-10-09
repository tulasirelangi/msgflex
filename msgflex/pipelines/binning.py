#!/usr/bin/env python3
"""
MSGFLEX - Binning/Split database MS-GF+ Pipeline.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import psutil
import subprocess
import time
import uuid
import threading
import numpy as np
import pandas as pd

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

from msgflex.core.process_registry import managed_popen, managed_run
from msgflex.core import ProjectPaths, SparxError, get_logger
from msgflex.modules import masic_merger
from msgflex.modules.dbcurator import remove_redundancy
from msgflex.modules.MASIC_wrapper import run_masic_batch
from msgflex.modules.tic_plot import generate_tic_plot
from msgflex.tools import tool

log = get_logger("pipelines.binning")

# E-VALUE CORRECTION UTILITIES
def parse_msgfplus_conf(conf_path):
    """Parse an MS-GF+ parameter file (key = value lines, # comments)."""
    params = {}
    if not conf_path or not os.path.exists(conf_path):
        return params
    with open(conf_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                params[key.strip()] = val.strip()
    return params
 
# Cleavage regex per MS-GF+ enzyme ID (pattern, terminus-side)
_ENZYME_PATTERNS = {
    0: None,  # unspecific cleavage
    1: (r"(?<=[KR])(?!P)", "C"), # Trypsin
    2: (r"(?<=[FYWL])(?!P)", "C"),  # Chymotrypsin
    3: (r"(?<=[K])(?!P)", "C"),  # Lys-C
    4: (r"(?=[K])", "N"), # Lys-N
    5: (r"(?<=[ED])(?!P)", "C"), # Glu-C
    6: (r"(?<=[R])(?!P)", "C"),  # Arg-C
    7: (r"(?=[DE])", "N"), # Asp-N
    9: None, # no cleavage
}
 
def count_tryptic_peptides(
    fasta_path,
    min_len = 6,
    max_len = 40,
    missed_cleavages = -1,
    ntt = 2,
    enzyme_id = 1,
    ):
    """
    Count candidate peptides in a FASTA matching the MS-GF+ search params.
    """
    import re as _re
 
    pattern_info = _ENZYME_PATTERNS.get(enzyme_id)
    count = 0
    seq_parts = []
 
    def _process(sequence):
        nonlocal count
        if not sequence:
            return
        if pattern_info is None:
            n = len(sequence)
            for length in range(min_len, min(max_len, n) + 1):
                count += n - length + 1
            return
        pattern, _ = pattern_info
        parts = _re.split(pattern, sequence)
        max_mc = missed_cleavages if missed_cleavages >= 0 else len(parts)
        if ntt == 2:
            for start in range(len(parts)):
                for end in range(start + 1, min(start + max_mc + 2, len(parts) + 1)):
                    pep = "".join(parts[start:end])
                    if min_len <= len(pep) <= max_len:
                        count += 1
        else:
            for start in range(len(parts)):
                for end in range(start + 1, min(start + max_mc + 2, len(parts) + 1)):
                    pep = "".join(parts[start:end])
                    plen = len(pep)
                    if plen < min_len:
                        continue
                    for sl in range(min_len, min(max_len, plen) + 1):
                        count += 1
 
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                _process("".join(seq_parts))
                seq_parts = []
            elif line:
                seq_parts.append(line)
    _process("".join(seq_parts))
    return count
 
def get_db_peptide_counts(
    db_name,
    split_dbs,
    database_path,
    conf_path,
    ):
    params = parse_msgfplus_conf(conf_path)
 
    def _int(key, default):
        try:
            return int(params.get(key, default))
        except ValueError:
            return default
 
    min_len = _int("minLength",6)
    max_len = _int("maxLength", 40)
    missed = _int("maxMissedCleavages", -1)
    ntt = _int("ntt",  2)
    enzyme_id = _int("e", 1)
 
    log.info(
        f"Peptide count params from conf: minLen={min_len}, maxLen={max_len}, "
        f"missedCleavages={missed}, ntt={ntt} enzyme={enzyme_id}")
 
    def _count(path):
        return count_tryptic_peptides(
            path, min_len=min_len, max_len=max_len,
            missed_cleavages=missed, ntt=ntt, enzyme_id=enzyme_id,
        )
 
    if split_dbs:
        part_counts = [_count(p) for p in split_dbs]
        full_count  = sum(part_counts)
    else:
        full_db_path = os.path.join(database_path, db_name)
        full_count   = _count(full_db_path)
        part_counts  = [full_count]
 
    return full_count, part_counts
 
def correct_evalue(df, n_full_pep, n_part_pep):
    """
    Multiplicative EValue correction using tryptic peptide counts.
    factor = N_peptides_full_DB / N_peptides_part
    EValue_corrected = EValue_part * factor
    """
    factor = n_full_pep / n_part_pep
    df = df.copy()
    df["MSGFDB_SpecEValue"] = df["SpecEValue"]
    if "EValue" in df.columns:
        df["EValue"] = df["EValue"] * factor
    df["CorrectionFactor"] = factor
    return df

def is_decoy_protein(protein_str, prefix = "XXX_"):
    """
    All proteins for this PSM are decoys (semicolon-separated entries).
    """
    proteins = [p.strip() for p in str(protein_str).split(";")]
    return all(p.startswith(prefix) for p in proteins)

def compute_qvalue(
    df,
    score_col = "MSGFDB_SpecEValue",
    decoy_col = "_IsDecoy",
):
    """
    PSM-level QValue via target-decoy FDR with monotone-non-decreasing fix.
    """
    order = df[score_col].argsort()  # best first (lowest SpecEValue)
    is_dec = df[decoy_col].values[order]

    n_decoy = np.cumsum(is_dec)
    n_target = np.cumsum(~is_dec)
    fdr = np.where(n_target > 0, (n_decoy + 1) / n_target, 1.0)

    qvalue_ordered = np.minimum.accumulate(fdr[::-1])[::-1]

    qvalue = np.empty(len(df))
    qvalue[order] = qvalue_ordered
    return qvalue

def compute_pep_qvalue(
    df,
    score_col = "MSGFDB_SpecEValue",
    decoy_col = "_IsDecoy",
    peptide_col = "Peptide",
):
    """
    Peptide-level QValue. Best PSM per bare-sequence peptide.
    """
    
    df = df.copy()
    df['_coreseq'] = (
        df[peptide_col]
        .str.replace(r'^[A-Z-]\.|\.[ A-Z-]$', '', regex=True)
        .str.replace(r'[^A-Z]', '', regex=True)
    )

    best_idx = df.groupby('_coreseq')[score_col].idxmin()
    best_df = df.loc[best_idx].copy().reset_index(drop=True)
    best_df["_PepQValue"] = compute_qvalue(best_df, score_col, decoy_col)

    pep_qval_map = best_df.set_index("_coreseq")["_PepQValue"].to_dict()
    return df["_coreseq"].map(pep_qval_map).values

def _get_core_config(total_ram_gb: float):
    """
    Detect core count and NUMA topology; derive safe thread/task/job settings.
    """
    total_cores = os.cpu_count() or 4

    numa_nodes = 1
    try:
        numa_path = "/sys/devices/system/node"
        if os.path.isdir(numa_path):
            nodes = [d for d in os.listdir(numa_path) if d.startswith("node") and d[4:].isdigit()]
            numa_nodes = max(1, len(nodes))
    except Exception:
        pass

    cores_per_node = max(1, total_cores // numa_nodes)
    return total_cores, numa_nodes, cores_per_node

NUMACTL_AVAILABLE = shutil.which("numactl") is not None

# PER-SAMPLE MS-GF+ WORKER
def _run_msgfplus_parts(
    mzML,
    db_name,
    datapath,
    database_path,
    conf_loc,
    java_mem,
    msgf_threads=None,
    msgf_tasks=None,
    numa_node=None
    ):
    """
    Run MS-GF+ against every split part for one mzML.
    """
    success_all = True
    db_prefix = os.path.splitext(db_name)[0]
    db_full_path = os.path.join(database_path, db_name)

    split_dbs = sorted(
        glob.glob(os.path.join(database_path, f"{db_prefix}_part*.fasta"))
        + glob.glob(os.path.join(database_path, f"{db_prefix}_part*.faa")),
        key=lambda p: int(m.group(1))
        if (m := re.search(r"_part(\d+)", os.path.basename(p))) else 0,
    )
    is_split = len(split_dbs) > 0

    if not is_split:
        split_dbs = [db_full_path]
    
    base_name = os.path.splitext(mzML)[0]

    # Skip when every part output already exists
    all_exist = all(
        os.path.exists(os.path.join(datapath, f"{base_name}_part{i}.mzid"))
        for i in range(1, len(split_dbs) + 1)
    )
    if all_exist:
        log.info("[SKIP] %s: all mzid parts exist.", mzML)
        return True

    for part_idx, db_part in enumerate(split_dbs, start=1):
        mzid_path = os.path.join(datapath, f"{base_name}_part{part_idx}.mzid")

        if os.path.exists(mzid_path):
            log.info(f"[SKIP] {mzML} part {part_idx}: mzid exists.")
            continue

        temp_fasta = os.path.join(
            database_path, f'temp_database_{uuid.uuid4().hex[:6]}.fasta'
        )
        try:
            shutil.copyfile(db_part, temp_fasta)
            cmd = [
                'java', f'-Xmx{java_mem}',
                '-jar', tool('MSGFPlus.jar'),
                '-s', os.path.join(datapath, mzML),
                '-d', temp_fasta,
                '-o', mzid_path,
                '-conf', str(conf_loc),
                '-addFeatures', '1',
                '-tda', '1',
            ]
            if msgf_threads is not None:
                cmd += ['-thread', str(msgf_threads)]
            if msgf_tasks is not None:
                cmd += ['-tasks', str(msgf_tasks)]
            if numa_node is not None and NUMACTL_AVAILABLE:
                cmd = ["numactl", f"--cpunodebind={numa_node}", f"--membind={numa_node}"] + cmd
            elif numa_node is not None:
                log.debug(f"numactl not found — running {mzML} part {part_idx} without NUMA pinning.")
            
            log.info(
                f"Running MS-GF+: {mzML} [part {part_idx}/{len(split_dbs)}:]"
                f"using:{os.path.basename(db_part)}"
            )

            process = managed_popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1
            )
            
            def _stream(pipe, log_func, label):
                for line in iter(pipe.readline, ''):
                    if line:
                        log_func(f"[{mzML} {label}] {line.strip()}")
                pipe.close()

            t_out = threading.Thread(target=_stream, args=(process.stdout, log.debug, part_idx))
            t_err = threading.Thread(target=_stream, args=(process.stderr, log.error, part_idx))

            t_out.start()
            t_err.start()

            process.wait()
            t_out.join()
            t_err.join()

            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, cmd)

            log.info(f"Success:{mzML} part {part_idx}")

        except Exception:
            log.exception(f"Failed: {mzML}")
            success_all = False
            continue
        
        #temp file clean-up
        finally:
            stem = os.path.splitext(temp_fasta)[0]
            temp_f = set([temp_fasta] + glob.glob(stem + ".*"))
            for t in temp_f:
                try:
                    if os.path.exists(t):
                        os.remove(t)
                except OSError as e:
                    log.warning(f"Could not remove temp file {t}: {e}")
    return success_all

# STEP 1.1: Per-part correction + TSV-level merge for one sample
def _merge_and_correct_parts(
    msfile,
    split_dbs,
    db_name,
    datapath,
    database_path,
    conf_loc,
    parts_tsv_dir,
    resultpath,
    java_mem
    ):
    """
    correction + merge for a single sample.
    """
    base_name = os.path.splitext(msfile)[0]

    part_mzids = sorted(
        glob.glob(os.path.join(datapath, f"{base_name}_part*.mzid")),
        key=lambda p: int(m.group(1))
        if (m := re.search(r"_part(\d+)", os.path.basename(p))) else 0,
    )
    if not part_mzids:
        log.warning(f"No part mzids found for {msfile}")
        return None

    n_full_aa, part_aa_counts = get_db_peptide_counts(
        db_name, split_dbs, database_path, conf_loc
    )
    log.info(f" DB size: {n_full_aa:,} peptides total across {len(split_dbs) or 1} part(s)")

    if len(part_aa_counts) != len(part_mzids):
        log.error(f"{msfile}: {len(part_mzids)} mzid parts but {len(part_aa_counts)} FASTA parts -- cannot correct.")
        return None

    corrected_frames: list[pd.DataFrame] = []

    for part_idx, (mzid_path, n_part_aa) in enumerate(
        zip(part_mzids, part_aa_counts), start=1
    ):
        tsv_part = os.path.join(parts_tsv_dir, f"{base_name}_part{part_idx}.tsv")

        if not os.path.exists(tsv_part):
            log.info(
                f" Converting part {part_idx}/{len(part_mzids)}"
                f"to {os.path.basename(mzid_path)} TSV"
            )
            managed_run(
                [
                    "java", f"-Xmx{java_mem}",
                    "-cp", tool("MSGFPlus.jar"),
                    "edu.ucsd.msjava.ui.MzIDToTsv",
                    "-i", mzid_path,
                    "-o", tsv_part,
                    "-showDecoy", "1",
                    "-unroll", "1",
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )
        else:
            log.info(f" [SKIP] TSV already exists for part {part_idx}")

        df_part = pd.read_csv(tsv_part, sep="\t", low_memory=False)
        df_part = correct_evalue(df_part, n_full_aa, n_part_aa)
        df_part["_PartIdx"] = part_idx

        factor = n_full_aa / n_part_aa
        log.info(f" Part {part_idx}: {n_part_aa:,} AA | factor: {factor:.4f} | {len(df_part):,} PSMs")
        corrected_frames.append(df_part)

    if not corrected_frames:
        return None

    merged = pd.concat(corrected_frames, ignore_index=True)
    log.info(f" Total PSMs after concat: {len(merged):,}")

    # Best PSM per spectrum (use SpecFile + ScanNum as composite key)
    specfile_col = merged.columns[0]  # always '#SpecFile' from MzIDToTsv
    spectrum_key = [specfile_col, "ScanNum"]
    merged = (
        merged
        .sort_values("MSGFDB_SpecEValue", ascending=True)
        .drop_duplicates(subset=spectrum_key, keep="first")
        .reset_index(drop=True)
    )
    log.info(
        f" Unique spectra after best-per-scan: {len(merged):,}"
    )

    merged["_IsDecoy"] = merged["Protein"].apply(
        lambda p: is_decoy_protein(p, "XXX_")
    )
    n_targets = (~merged["_IsDecoy"]).sum()
    n_decoys = merged["_IsDecoy"].sum()
    log.info(f" Target PSMs: {n_targets:,} | Decoy PSMs: {n_decoys:,}")

    merged["QValue"] = compute_qvalue(merged, "MSGFDB_SpecEValue", "_IsDecoy")
    merged["PepQValue"] = compute_pep_qvalue(
        merged, "MSGFDB_SpecEValue", "_IsDecoy", "Peptide"
    )

    out_df = merged.drop(
        columns=["_IsDecoy", "_PartIdx", "CorrectionFactor", "_coreseq"],
        errors="ignore",
    )
    out_tsv = os.path.join(resultpath, f"{base_name}.tsv")
    out_df.to_csv(out_tsv, sep="\t", index=False)
    log.info(f"Written: {out_tsv}")
    return out_tsv

# ----helpers---
def _normalize(p):
    return os.path.abspath(os.path.expanduser(p))

# Call functions to run pipeline
def run(
    paths,
    conf_loc,
    java_mem = "4G",
    *,
    max_workers: Optional[int] = None,
    msgf_tasks: Optional[int] = None,
    stop_event: Optional[threading.Event] = None
    ):
    """
    Run the binning MS-GF+ pipeline end-to-end.
    """
    log.info(f"Binning pipeline => base: {paths.base}")
    log.info(f"Params file: {conf_loc}")
    log.info(f"Java heap: {java_mem}")

    paths.require_sparx_inputs()
    if not paths.inputfile.exists():
        raise SparxError(
            f"Binning pipeline requires a pre-built inputfile.tsv at {paths.inputfile}."
        )
    paths.ensure_dirs()

    datapath = str(paths.data) + os.sep
    database_path = str(paths.database) + os.sep
    QCdir = str(paths.qcdir) + os.sep
    SICdir = str(paths.sicdir)
    resultpath = str(paths.results) + os.sep
    PHRPOut = str(paths.phrp_out)
    MasicOut = str(paths.masic_out)

    parts_tsv_dir  = str(paths.base / "parts_tsv")
    os.makedirs(parts_tsv_dir, exist_ok=True)

    decoder = pd.read_csv(paths.inputfile, sep="\t", header=0)
    log.info(f"Loaded {len(decoder)} rows from inputfile.tsv")

    start_time = time.time()

    # STEP 0: RAW -> Spectrum Files + QC plots
    VALID_EXTENSIONS = (".mzML", ".mgf", ".mzXML", ".ms2")
    # 1. Check if ALL samples have AT LEAST ONE valid spectrum file
    _step0_spec_exist = all(
        any(
            os.path.exists(os.path.join(datapath, str(row["msfilename"]) + ext))
            for ext in VALID_EXTENSIONS
        )
        for _, row in decoder.iterrows()
    )
    # 2. Check if QC plots exist, but ONLY for samples that have an .mzML
    _step0_qc_exist = all(
        os.path.exists(os.path.join(QCdir, str(row["QCplot"]) + ".png"))
        for _, row in decoder.iterrows()
        if os.path.exists(os.path.join(datapath, str(row["msfilename"]) + ".mzML"))
    )

    if _step0_spec_exist and _step0_qc_exist:
        log.info("[SKIP] Step-0: all valid spectrum files and QC plots already exist.")
    else:
        processed_db: set[str] = set()
        for _, row in decoder.iterrows():
            msfile  = str(row["msfilename"])
            msraw   = msfile + ".raw"
            msmzml  = msfile + ".mzML"
            db_name = str(row["database"])
            outname = str(row["QCplot"]) + ".png"

            log.info("-" * 40)
            log.info(f"{msraw}")

            if db_name not in processed_db:
                input_fasta = os.path.join(database_path, db_name)
                try:
                    log.info(f"Removing redundancy: {db_name}...")
                    remove_redundancy(input_fasta, input_fasta)
                    processed_db.add(db_name)
                except FileNotFoundError:
                    log.warning(f"FASTA not found: {input_fasta}. Skipping redundancy removal.")

            # Check if ANY valid spectrum file exists for this specific sample
            spec_exists = any(
                os.path.exists(os.path.join(datapath, msfile + ext)) 
                for ext in VALID_EXTENSIONS
            )

            # Only convert from RAW if NO valid spectrum files were found
            if not spec_exists:
                log.info("No compatible spectrum files detected. Converting .raw → .mzML ...")
                managed_run(
                    [
                        "mono",
                        tool("ThermoRawFileParser/ThermoRawFileParser.exe"),
                        "-i", os.path.join(datapath, msraw),
                        "-f", "2",
                        "-L", "1-",
                    ],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT
                )

            mzml_full = os.path.join(datapath, msmzml)
            if os.path.exists(mzml_full):
                try:
                    spectrum_count = 0
                    # errors='ignore' prevents crashes from weird vendor-encoded characters
                    with open(mzml_full, 'r', errors='ignore') as f:
                        for line in f:
                            if "<spectrum " in line:  # Note the space to avoid matching <spectrumList>
                                spectrum_count += 1
                    log.info(f"Conversion complete. Detected {spectrum_count} spectra in {msmzml}")
                except Exception as e:
                    log.warning(f"Could not count spectra in {msmzml}: {e}")
            
            qc_full = os.path.join(QCdir, outname)
            
            if os.path.exists(mzml_full) and not os.path.exists(qc_full):
                generate_tic_plot(mzml_full, qc_full)
            elif os.path.exists(qc_full):
                log.info(f"[SKIP] QC plot exists: {outname}")

        log.info("Step-0 complete.")

    # STEP 1: MS-GF+ search (multi-format aware)
    VALID_EXTENSIONS = (".mzML", ".mgf", ".mzXML", ".ms2")
    
    peak_to_db = {}
    for _, row in decoder.iterrows():
        msfile = str(row['msfilename'])
        db_name = str(row['database'])
        
        # Check which format actually exists on disk for this sample
        for ext in VALID_EXTENSIONS:
            if os.path.exists(os.path.join(datapath, msfile + ext)):
                peak_to_db[msfile + ext] = db_name
                break
        else:
            log.warning(f"No compatible spectrum file found on disk for entry: {msfile}")

    assert isinstance(peak_to_db, dict)
    
    # Look for all valid formats in the directory
    peak_files = [f for f in os.listdir(datapath) if f.endswith(VALID_EXTENSIONS)]
    num_files  = len(peak_files)
    
    if not peak_files:
        raise SparxError(f"No compatible spectrum files ({', '.join(VALID_EXTENSIONS)}) found in {datapath}")
    
    # --- Universal Resource Guard ---
    try:
        total_ram_gb = psutil.virtual_memory().total / (1024**3)
    except ImportError:
        log.warning("psutil not installed. Assuming 16GB RAM for safety limits.")
        total_ram_gb = 32.0
        
    # binning pipeline is recommended to run on dedicated workstations with minimum of 64GB RAM
    mem_per_job_gb = int(java_mem.upper().replace("G", ""))
    if mem_per_job_gb < 4:
        log.warning(
            f"java_mem={java_mem} is very low for multi-threading MS-GF+ search."
            "Consider at least 4GB per job to avoid heap OOM")
    safe_ram_gb = total_ram_gb * 0.90
    max_jobs_by_ram = max(1, int(safe_ram_gb / mem_per_job_gb))
    max_jobs_by_cpu = 8

    total_cores, numa_nodes, cores_per_node = _get_core_config(total_ram_gb)

    MIN_THREADS_PER_JOB = 4
    max_jobs_by_cpu = max(1, total_cores // MIN_THREADS_PER_JOB)
    num_jobs = max(1, min(max_jobs_by_ram, max_jobs_by_cpu, num_files))
    jobs_per_node = max(1, -(-num_jobs // numa_nodes))

    threads_per_job = max(1, cores_per_node // jobs_per_node)
    TASK_MULTIPLIER = 2
    if msgf_tasks is not None:
        tasks_per_job = max(1, min(msgf_tasks, threads_per_job))
    else:
        tasks_per_job = threads_per_job * TASK_MULTIPLIER
    
    log.info(f"System Resources -> Total RAM: {total_ram_gb}GB | Java/job: {mem_per_job_gb}GB")
    log.info(f"Concurrency Limits -> By RAM: {max_jobs_by_ram} | By CPU Cap: {max_jobs_by_cpu} | Final Jobs: {num_jobs}")
    log.info(f"Spectrum files: {len(peak_files)}  |  parallel jobs: {num_jobs}")
    log.info(f"Thread Allocation -> Total cores: {total_cores} | jobs: {num_jobs} | Threasds/job: {cores_per_node} | Tasks/job: {tasks_per_job}")
    log.info(f"Spectrum files: {len(peak_files)}  |  Parallel jobs: {num_jobs}")

    # Filter out any orphan files that exist on disk but aren't in inputfile.tsv
    valid_peak_files = [f for f in peak_files if f in peak_to_db]
    
    def _step1_all_outputs_exist() -> bool:
        for peak_file in valid_peak_files:
            db_name   = peak_to_db[peak_file]
            db_prefix = os.path.splitext(db_name)[0]
            base_name = os.path.splitext(peak_file)[0] # Works perfectly for .mgf, .mzML, etc.
            
            split_dbs = sorted(
                glob.glob(os.path.join(database_path, f"{db_prefix}_part*.fasta"))
                + glob.glob(os.path.join(database_path, f"{db_prefix}_part*.faa")),
                key=lambda p: int(m.group(1))
                if (m := re.search(r"_part(\d+)", os.path.basename(p))) else 0,
            )
            n_parts = len(split_dbs) if split_dbs else 1
            
            for idx in range(1, n_parts + 1):
                if not os.path.exists(
                    os.path.join(datapath, f"{base_name}_part{idx}.mzid")
                ):
                    return False
        return True

    if valid_peak_files and _step1_all_outputs_exist():
        log.info("[SKIP] Step-1: all mzid part files already exist.")
    else:
        failed_jobs: list[str] = []
        with ThreadPoolExecutor(max_workers=num_jobs) as executor:
            future_to_peak = {}
            for idx, peak_file in enumerate(valid_peak_files):
                numa_node = idx % numa_nodes if numa_nodes > 1 else None
                future = executor.submit(
                    _run_msgfplus_parts,
                    peak_file, peak_to_db[peak_file],
                    datapath, database_path, conf_loc, java_mem,
                    msgf_threads=threads_per_job,
                    msgf_tasks=tasks_per_job,
                    numa_node=numa_node,
                )
                future_to_peak[future] = peak_file

            for future in as_completed(future_to_peak):
                peak_file = future_to_peak[future]
                try:
                    if not future.result():
                        failed_jobs.append(peak_file)
                except Exception as e:
                    log.error(f"Unhandled error for {peak_file}: {e}")
                    failed_jobs.append(peak_file)

        if failed_jobs:
            log.warning(f"{len(failed_jobs)} job(s) failed: {failed_jobs}")

    # STEP 1.2 -- Per-part correction + TSV-level merge
    _step1_1_all_tsvs_exist = all(
        os.path.exists(os.path.join(
            resultpath, f"{os.path.splitext(str(row['msfilename']).strip())[0]}.tsv"
        ))
        for _, row in decoder.iterrows()
    )
    if _step1_1_all_tsvs_exist:
        log.info("[SKIP] Step-1.5: all corrected TSVs already exist.")
    else:
        correction_failed: list[str] = []
        for _, row in decoder.iterrows():
            msfile = str(row["msfilename"]).strip()
            db_name = str(row["database"])

            log.info("=" * 55)
            log.info(f"Step-1.1: correcting {msfile}")

            out_tsv = os.path.join(
                resultpath, f"{os.path.splitext(msfile)[0]}.tsv"
            )
            if os.path.exists(out_tsv):
                log.info(f"[SKIP] Corrected TSV already exists: {out_tsv}")
                continue

            db_prefix = os.path.splitext(db_name)[0]
            split_dbs = sorted(
                glob.glob(os.path.join(database_path, f"{db_prefix}_part*.fasta"))
                + glob.glob(os.path.join(database_path, f"{db_prefix}_part*.faa")),
                key=lambda p: int(m.group(1))
                if (m := re.search(r"_part(\d+)", os.path.basename(p))) else 0,
            )

            result = _merge_and_correct_parts(
                msfile, split_dbs, db_name,
                datapath, database_path, conf_loc, parts_tsv_dir, resultpath, java_mem,
            )
            if result is None:
                correction_failed.append(msfile)

        if correction_failed:
            log.warning(f"Correction failed for: {correction_failed}")
        log.info("Step-1.1 complete.")
 
    # STEP 1.2: Per-part mzid merging for external 
    _step1_2_all_mzids_exist = all(
        os.path.exists(os.path.join(
            datapath,
            f"{os.path.splitext(str(row['msfilename']).strip())[0]}.mzid",
        ))
        for _, row in decoder.iterrows()
    )
    if _step1_2_all_mzids_exist:
        log.info("[SKIP] Step-1.6: all merged mzid files already exist.")
    else:
        merge_failed: list[str] = []

        for _, row in decoder.iterrows():
            msfile   = os.path.splitext(str(row["msfilename"]).strip())[0]
            mzid_out = os.path.join(datapath, f"{msfile}.mzid")

            if os.path.exists(mzid_out):
                log.info("[SKIP] Merged mzid already exists: %s.mzid", msfile)
                continue

            pattern = os.path.join(datapath, f"{msfile}_part*.mzid")
            mzid_parts = sorted(
                glob.glob(pattern),
                key=lambda p: int(m.group(1))
                if (m := re.search(r"_part(\d+)", os.path.basename(p))) else 0,
            )

            log.info(f"Searching for part mzids: {pattern}")

            if not mzid_parts:
                log.warning(f"No mzid parts found for {msfile}")
                merge_failed.append(msfile)
                continue

            tmp_dir = os.path.join(datapath, f"_tmp_mzidmerge_{uuid.uuid4().hex[:6]}")
            os.makedirs(tmp_dir, exist_ok=True)

            try:
                for f in mzid_parts:
                    shutil.copy2(f, os.path.join(tmp_dir, os.path.basename(f)))

                managed_run(
                    [
                        "mono",
                        tool("MzidMerger/net48/MzidMerger.exe"),
                        "-InDir", tmp_dir,
                        "-Out", mzid_out,
                        "-KeepOnlyBestResults",
                    ],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT
                )
                log.info(f"Merged mzid ==> {mzid_out}")

                if os.path.exists(mzid_out):
                    for _part_mzid in mzid_parts:
                        try:
                            os.remove(_part_mzid)
                        except OSError as _e:
                            log.warning(f"Could not remove part mzid {_part_mzid}: {_e}")
                    log.info(f"Removed {len(mzid_parts)} part mzid(s) for {msfile}")

            except Exception as e:
                log.error(f"MzidMerger failed for {msfile}: {e}")
                merge_failed.append(msfile)

            finally:
                shutil.rmtree(tmp_dir, ignore_errors=True)

        if merge_failed:
            log.warning(f"MzidMerger failed for: {merge_failed}")
        log.info("Step-1.2 complete.")

    if os.path.exists(parts_tsv_dir):
        shutil.rmtree(parts_tsv_dir)
        log.info(f"Removed parts_tsv directory: {parts_tsv_dir}")

    # STEP 3a. pre-pass: Reconstruct split FASTAs for PHRP
    _all_dbs_exist = all(
        os.path.exists(os.path.join(database_path, str(row["database"])))
        for _, row in decoder.iterrows()
    )
    if _all_dbs_exist:
        log.info("[SKIP] DB reconstruction: all full DB files already exist.")
    else:
        reconstructed: set[str] = set()
        for _, row in decoder.iterrows():
            db_name = str(row["database"])
            if db_name in reconstructed:
                continue

            db_prefix = os.path.splitext(db_name)[0]
            expected_db = os.path.join(database_path, db_name)

            split_dbs = sorted(
                glob.glob(os.path.join(database_path, f"{db_prefix}_part*.fasta"))
                + glob.glob(os.path.join(database_path, f"{db_prefix}_part*.faa")),
                key=lambda p: int(m.group(1))
                if (m := re.search(r"_part(\d+)", os.path.basename(p))) else 0,
            )

            if split_dbs:
                log.info(f"Reconstructing {db_name} from {len(split_dbs)} parts ...")
                with open(expected_db, "wb") as outfile:
                    for part in split_dbs:
                        with open(part, "rb") as infile:
                            shutil.copyfileobj(infile, outfile)
                size_mb = os.path.getsize(expected_db) / 1e6
                log.info(f"{db_name} reconstructed ({size_mb:.2f}MB).")
                for part in split_dbs:
                    os.remove(part)
                log.info(f"Removed {len(split_dbs)} part files.")
            else:
                log.info(f"Using existing DB: {db_name}")

            reconstructed.add(db_name)
    
    # STEP 3: PHRP
    decoder_db_lookup: dict[str, str] = {
        os.path.splitext(str(row["msfilename"]))[0]: str(row["database"])
        for _, row in decoder.iterrows()
    }

    _step3_tsv_files = [f for f in os.listdir(resultpath) if f.endswith(".tsv")]
    _step3_phrp_done = (
        os.path.isdir(PHRPOut) and bool(os.listdir(PHRPOut))
        and all(
            any(f.startswith(os.path.splitext(tsv)[0]) for f in os.listdir(PHRPOut))
            for tsv in _step3_tsv_files
        )
    )
    if _step3_phrp_done and _step3_tsv_files:
        log.info("[SKIP] Step-3: PHRP outputs already exist.")
    else:
        os.makedirs(PHRPOut, exist_ok=True)
        tsv_files = [f for f in os.listdir(resultpath) if f.endswith(".tsv")]
        log.info(f"Starting PHRP for {len(tsv_files)} TSV file(s) ...")

        for tsv_file in tsv_files:
            tsv_path = os.path.join(resultpath, tsv_file)
            base_name = os.path.splitext(tsv_file)[0]
            sample_out_dir = os.path.join(PHRPOut, base_name)
            os.makedirs(sample_out_dir, exist_ok=True)

            base_stripped = os.path.splitext(base_name)[0]
            db_name = (
                decoder_db_lookup.get(base_name)
                or decoder_db_lookup.get(base_stripped)
            )
            if db_name is None:
                log.warning(f"No DB match for {tsv_file}; falling back to first DB in manifest.")
                db_name = str(decoder.iloc[0]["database"])

            log.info(f"PHRP: {tsv_file}  DB: {db_name}")

            managed_run(
                [
                    "mono",
                    _normalize(tool("PHRP/PeptideHitResultsProcRunner.exe")),
                    f"/I:{_normalize(tsv_path)}",
                    f"/M:{_normalize(tool('PHRP/MSGFDB_Mods.txt'))}",
                    f"/N:{_normalize(conf_loc)}",
                    f"/T:{_normalize(tool('PHRP/Mass_Correction_Tags.txt'))}",
                    f"/F:{_normalize(os.path.join(database_path, db_name))}",
                    f"/O:{_normalize(sample_out_dir)}",
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            log.info(f"PHRP done: {tsv_file}")

        # Flatten PHRPOut subdirectories
        for sample_dir in os.listdir(PHRPOut):
            sample_dir_path = os.path.join(PHRPOut, sample_dir)
            if os.path.isdir(sample_dir_path):
                for fname in os.listdir(sample_dir_path):
                    shutil.move(
                        os.path.join(sample_dir_path, fname),
                        os.path.join(PHRPOut, fname),
                    )
                os.rmdir(sample_dir_path)

        log.info(f"All PHRP outputs consolidated into PHRPOut/")

    # STEP 4: MASIC SIC generation
    os.makedirs(MasicOut, exist_ok=True)
    msraw_files = [f for f in os.listdir(datapath) if f.endswith(".raw")]
    if not msraw_files:
        log.warning(f"No .raw files found in {datapath}; MASIC SIC generation will be skipped.")

    _step4_masic_done = (
        os.path.isdir(MasicOut) and bool(os.listdir(MasicOut))
        and all(
            any(f.startswith(os.path.splitext(raw)[0]) for f in os.listdir(MasicOut))
            for raw in msraw_files
        )
    )
    if _step4_masic_done:
        log.info("[SKIP] Step-4: MASIC outputs already exist.")
    else:
        log.info("Starting MASIC SIC generation ...")
        successful, failed, skipped = run_masic_batch(
            datapath, MasicOut,
            workers=4,
            masic_exe=tool("MASIC/MASIC_Console.exe"),
            masic_params=tool("MASICParameters.xml"),
            stop_event=stop_event,
        )
        log.info(f"MASIC: {successful} succeeded, {failed} failed, {skipped} skipped.")

    # STEP 5: MASIC merger
    _step5_sic_done = os.path.isdir(SICdir) and bool(os.listdir(SICdir))
    if _step5_sic_done:
        log.info("[SKIP] Step-5: merged SIC outputs already exist.")
    else:
        try:
            temp_merger = os.path.join(resultpath, "temp_merger")
            masic_merger.copy_files_to_temp(PHRPOut, MasicOut, temp_merger)
            masic_merger.merge_files(temp_merger, SICdir)
            log.info("All intensity merging completed.")
        except Exception as e:
            log.error(f"Error in MASIC merger: {e}")
            raise SparxError(f"MASIC merger failed: {e}") from e

    elapsed_hours = (time.time() - start_time) / 3600
    log.info("*" * 60)
    log.info(f"Binning pipeline completed! Time: {elapsed_hours:.2f} hours")
    log.info("*" * 60)
