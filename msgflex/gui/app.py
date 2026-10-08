#!/usr/bin/env python3
"""
MSGFLEX GUI application shell.
"""

from __future__ import annotations
import threading
import os
from pathlib import Path

import dearpygui.dearpygui as dpg

from msgflex import __version__
from msgflex.gui.fonts import register_fonts
from msgflex.gui.logo import _draw_logo_on
from msgflex.gui.log_pane import LogPane
from msgflex.gui.state import AppState
from msgflex.gui.tabs import check, full_run, project, quantix, sparx, xpectra
from msgflex.gui.theme import TOKENS, apply_global_theme
from msgflex.core.process_registry import ProcessRegistry

_ASSETS = Path(__file__).parent / "assets"
_ICON_LARGE = _ASSETS / "msgflex.png"
_ICON_SMALL = _ASSETS / "msgflex_64.png"

def _on_quit_clicked():
    threading.Thread(target=ProcessRegistry.kill_all, daemon=True).start()
    dpg.stop_dearpygui()

def _build_header():
    """Top header bar: logo + wordmark + tagline + version """
    with dpg.group(horizontal=True):
        # Logo
        logo_dl = dpg.add_drawlist(width=28, height=28)
        _draw_logo_on(logo_dl, size=28)

        dpg.add_spacer(width=6)

        # Wordmark + tagline stacked
        with dpg.group():
            dpg.add_text("MSGFLEX", color=TOKENS["text"])
            dpg.add_text(
                "PROTEOMICS · UNIFIED",
                color=TOKENS["accent"],
            )

        # Pushes everything after this to the right
        dpg.add_spacer(width=-1)

        # Right side — version + project path + Quit
        dpg.add_text(f"v{__version__}", color=TOKENS["accent"], tag="header_version")
        dpg.add_text(" · ", color=TOKENS["text_muted"])
        dpg.add_text(
            "(no project)",
            tag="header_project_path",
            color=TOKENS["text_muted"],
        )

def _build_status_bar(themes: dict) -> None:
    """Bottom status bar: running state + summary + quit."""
    with dpg.group(horizontal=True):
        dpg.add_text("●", color=TOKENS["text_faint"], tag="statusbar_dot")
        dpg.add_text(
            "Ready — select a project to begin",
            tag="statusbar_text",
            color=TOKENS["text"],
        )
        dpg.add_spacer(width=-1)
        dpg.add_text(
            "",
            tag="statusbar_summary",
            color=TOKENS["text_muted"],
        )
        dpg.add_spacer(width=12)
        quit_btn = dpg.add_button(
            label="Quit",
            width=70, height=24,
            callback=_on_quit_clicked,
            tag="header_quit_btn",
        )
        dpg.bind_item_theme(quit_btn, themes["danger"])

def _refresh_header_and_status(state: AppState) -> None:
    """Subscribed to state changes — updates the header path + status bar."""
    if dpg.does_item_exist("header_project_path"):
        dpg.set_value(
            "header_project_path",
            str(state.paths.base) if state.paths else "(no project)",
        )
    if dpg.does_item_exist("statusbar_text"):
        dpg.set_value(
            "statusbar_text",
            "Ready" if state.paths else "Ready — select a project to begin",
        )

def build_app(state: AppState, themes: dict, fonts: dict | None = None) -> None:
    """Assemble the full app inside the primary window."""
    fonts = fonts or {}
    state.subscribe(_refresh_header_and_status)

    with dpg.window(
        tag="primary",
        no_title_bar=True,
        no_resize=False,
        no_move=True,
        no_close=True,
        no_collapse=True,
    ):
        #  Header
        _build_header()
        dpg.add_separator()

        # Tab + log split 
        with dpg.group(horizontal=True):
            # Left: tab bar and contents (flex)
            with dpg.child_window(width=-300, height=-40, border=False):
                with dpg.tab_bar(tag="main_tabs"):
                    with dpg.tab(label="Project"):
                        with dpg.group() as g:
                            pass
                        project.build(g, state, _log_pane, themes)

                    with dpg.tab(label="SPARX"):
                        with dpg.group() as g:
                            pass
                        sparx.build(g, state, _log_pane, themes)

                    with dpg.tab(label="XPECTRA"):
                        with dpg.group() as g:
                            pass
                        xpectra.build(g, state, _log_pane, themes)

                    with dpg.tab(label="QUANTIX"):
                        with dpg.group() as g:
                            pass
                        quantix.build(g, state, _log_pane, themes)

                    with dpg.tab(label="Full run"):
                        with dpg.group() as g:
                            pass
                        full_run.build(g, state, _log_pane, themes)

                    with dpg.tab(label="Check"):
                        with dpg.group() as g:
                            pass
                        check.build(g, state, _log_pane, themes)

            # Right: live log pane (fixed 290px)
            with dpg.child_window(width=290, height=-40, border=True) as log_container:
                _log_pane.parent = log_container
                _log_pane.build()
                # Bind monospace font to the log content for code-like readability
                if fonts.get("mono") and _log_pane._content_id:
                    dpg.bind_item_font(_log_pane._content_id, fonts["mono"])

        # Status bar
        dpg.add_separator()
        _build_status_bar(themes)

    dpg.set_primary_window("primary", True)

    # Welcome log line
    _log_pane.append(f"MSGFLEX v{__version__} ready", "success")
    _log_pane.append("Open the Project tab to select a base directory.", "info")

def run_gui() -> None:
    """Main entry point — called by `msgflex gui`."""
    global _log_pane, _fonts

    dpg.create_context()
    from msgflex.core.logging import setup_logging
    setup_logging(mode="gui", verbosity="normal")

    _fonts = register_fonts(base_size=16, mono_size=14, global_scale=1.0)

    state = AppState()
    themes = apply_global_theme()
    _log_pane = LogPane(parent=0, state=state)  # parent assigned in build_app

    # additonal icon kwargs for Windows and Linux
    icon_kwargs = {}
    if _ICON_LARGE.is_file() and _ICON_SMALL.is_file():
        icon_kwargs = dict(small_icon=str(_ICON_SMALL), large_icon=str(_ICON_LARGE))
    
    # Viewport — resizes naturally
    dpg.create_viewport(
        title=f"MSGFLEX v{__version__}",
        width=1240,
        height=720,
        min_width=980,
        min_height=600,
        **icon_kwargs
    )

    build_app(state, themes, fonts=_fonts)

    dpg.setup_dearpygui()
    dpg.show_viewport()
    while dpg.is_dearpygui_running():
        _log_pane._drain()
        dpg.render_dearpygui_frame()

    # so any in-flight worker thread DPG calls become no-ops.
    _log_pane._content_id = 0

    # Kill child subprocesses (Java / Mono) synchronously
    ProcessRegistry.kill_all(grace_seconds=3.0)

    import time
    for _ in range(10):
        active = [t for t in threading.enumerate()
                  if t.name.startswith("Thread-") and t.is_alive()
                  and t is not threading.current_thread()]
        if not active:
            break
        time.sleep(0.1)

    try:
        dpg.destroy_context()
    except Exception:
        pass

    os._exit(0)

if __name__ == "__main__":
    run_gui()
