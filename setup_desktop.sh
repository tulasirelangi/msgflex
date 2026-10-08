#!/usr/bin/env bash

set -euo pipefail

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m';  GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; BOLD='\033[1m';     RESET='\033[0m'

info()    { echo -e "${BLUE}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; }
die()     { error "$*"; exit 1; }
header()  { echo -e "\n${BOLD}── $* ──${RESET}"; }

# ── Resolve bundle directory ──────────────────────────────────────────────────
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_SCRIPT="$BUNDLE_DIR/install.sh"

# ── Paths ─────────────────────────────────────────────────────────────────────
ICON_SRC="$BUNDLE_DIR/msgflex/gui/assets/msgflex.png" # optional
ICON_DIR="$HOME/.local/share/icons/hicolor/256x256/apps"
ICON_DEST="$ICON_DIR/MSGFLEX.png"
APP_DIR="$HOME/.local/share/applications"
DESKTOP_DEST="$APP_DIR/MSGFLEX.desktop"

# =============================================================================
# Banner
# =============================================================================
echo -e "${BOLD}"
echo "  ╔══════════════════════════════════════════╗"
echo "  ║   MSGFLEX Pipeline v1.1.0 Desktop Setup  ║"
echo "  ╚══════════════════════════════════════════╝"
echo -e "${RESET}"
echo "  Bundle directory : $BUNDLE_DIR"
echo "  install.sh       : $INSTALL_SCRIPT"
echo ""

# Preflight
header "Checking requirements"

if [[ "$(uname -s)" != "Linux" ]]; then
    die "Desktop setup is Linux-only (.desktop files are a Linux standard)."
fi
success "Linux detected"

if [[ ! -f "$INSTALL_SCRIPT" ]]; then
    die "install.sh not found at $BUNDLE_DIR — run setup_desktop.sh from inside the MSGFLEX bundle."
fi
success "install.sh found"

# Detect a usable terminal emulator (tries common ones in priority order)
detect_terminal() {
    local terminals=(gnome-terminal konsole xfce4-terminal mate-terminal xterm)
    for t in "${terminals[@]}"; do
        if command -v "$t" &>/dev/null; then
            echo "$t"
            return 0
        fi
    done
    return 1
}

TERMINAL="$(detect_terminal || true)"
if [[ -z "$TERMINAL" ]]; then
    warn "No supported terminal emulator found."
    warn "Supported: gnome-terminal, konsole, xfce4-terminal, mate-terminal, xterm"
    warn "The .desktop file will still be created but may not open automatically."
else
    success "Terminal emulator detected: $TERMINAL"
fi

# =============================================================================
# Step 1 — Install icon
# =============================================================================
header "Step 1/4  Installing icon"

mkdir -p "$ICON_DIR"

if [[ -f "$ICON_SRC" ]]; then
    cp "$ICON_SRC" "$ICON_DEST"
    ICON_NAME="MSGFLEX"
    success "Icon installed: $ICON_DEST"
else
    warn "msgflex.png not found in bundle — using system fallback icon (application-x-executable)."
    warn "To add a custom icon, place msgflex.png (256×256) in: $BUNDLE_DIR"
    ICON_NAME="application-x-executable"
fi

# Refresh icon cache if gtk-update-icon-cache is available
if command -v gtk-update-icon-cache &>/dev/null; then
    gtk-update-icon-cache -f -t "$HOME/.local/share/icons/hicolor" 2>/dev/null || true
    success "Icon cache refreshed"
fi

# =============================================================================
# Step 2 — Build the Exec line based on detected terminal
# =============================================================================
header "Step 2/4  Building launcher command"

# Each terminal has a different syntax for "run this command and hold open"
case "$TERMINAL" in
    gnome-terminal)
        EXEC_CMD="gnome-terminal -- bash -c 'bash \"${INSTALL_SCRIPT}\"; echo \"\"; echo \"Press Enter to close...\"; read'"
        ;;
    konsole)
        EXEC_CMD="konsole --hold -e bash -c 'bash \"${INSTALL_SCRIPT}\"'"
        ;;
    xfce4-terminal)
        EXEC_CMD="xfce4-terminal --hold -e 'bash -c \"bash \\\"${INSTALL_SCRIPT}\\\"; echo \\\"\\\"; echo \\\"Press Enter to close...\\\"; read\"'"
        ;;
    mate-terminal)
        EXEC_CMD="mate-terminal -- bash -c 'bash \"${INSTALL_SCRIPT}\"; echo \"\"; echo \"Press Enter to close...\"; read'"
        ;;
    xterm)
        EXEC_CMD="xterm -hold -e bash -c 'bash \"${INSTALL_SCRIPT}\"'"
        ;;
    *)
        # Fallback: try xterm which is nearly universal
        EXEC_CMD="xterm -hold -e bash -c 'bash \"${INSTALL_SCRIPT}\"'"
        ;;
esac

success "Exec line built for: ${TERMINAL:-xterm (fallback)}"

# =============================================================================
# Step 3 — Write .desktop file
# =============================================================================
header "Step 3/4  Writing .desktop file"

mkdir -p "$APP_DIR"

cat > "$DESKTOP_DEST" << EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=MSGFLEX Installer
GenericName=Proteomics Pipeline Installer
Comment=Install the MSGFLEX Pipeline v0.5.0 (SPARX + XPECTRA + QUANTIX)
Icon=${ICON_NAME}
Exec=bash -c 'bash "${INSTALL_SCRIPT}"; echo; read -rp "Press Enter to close..."'
Terminal=true
Categories=Science;Biology;
Keywords=proteomics;mass-spectrometry;MSGFLEX;pipeline;
StartupNotify=true
Path=${BUNDLE_DIR}
EOF

chmod +x "$DESKTOP_DEST"
success "Desktop entry written: $DESKTOP_DEST"

# Mark as trusted (required by GNOME 3.x+ to allow execution on double-click)
if command -v gio &>/dev/null; then
    gio set "$DESKTOP_DEST" metadata::trusted true 2>/dev/null || true
    success "Marked as trusted (gio)"
elif command -v dbus-launch &>/dev/null; then
    # Older GNOME fallback
    dbus-launch gio set "$DESKTOP_DEST" metadata::trusted true 2>/dev/null || true
fi

# Refresh application database
if command -v update-desktop-database &>/dev/null; then
    update-desktop-database "$APP_DIR" 2>/dev/null || true
    success "Application database refreshed"
fi

# =============================================================================
# Step 4 — Optional Desktop shortcut
# =============================================================================
header "Step 4/4  Desktop shortcut"

DESKTOP_HOME="$HOME/Desktop"

if [[ -d "$DESKTOP_HOME" ]]; then
    read -rp "  Place shortcut on ~/Desktop? [Y/n]: " SHORTCUT_CONFIRM
    if [[ ! "$SHORTCUT_CONFIRM" =~ ^[Nn]$ ]]; then
        cp "$DESKTOP_DEST" "$DESKTOP_HOME/MSGFLEX.desktop"
        chmod +x "$DESKTOP_HOME/MSGFLEX.desktop"

        # Trust the Desktop copy too
        if command -v gio &>/dev/null; then
            gio set "$DESKTOP_HOME/MSGFLEX.desktop" metadata::trusted true 2>/dev/null || true
        fi

        success "Shortcut placed: $DESKTOP_HOME/MSGFLEX.desktop"
    else
        info "Skipping Desktop shortcut."
    fi
else
    warn "~/Desktop not found — skipping Desktop shortcut."
    warn "You can still launch the installer from your application menu."
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo -e "${BOLD}══════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Desktop setup complete                          ${RESET}"
echo -e "${BOLD}══════════════════════════════════════════════════${RESET}"
echo ""
echo "  Desktop entry : $DESKTOP_DEST"
echo "  Icon ${ICON_NAME}"
echo "  Terminal : ${TERMINAL:-xterm (fallback)}"
echo "  Bundle : $BUNDLE_DIR"
echo ""
echo -e "  ${BOLD}What to do next:${RESET}"
echo "    • Look for 'MSGFLEX Installer' in your application menu, or"
echo "    • Double-click the icon on ~/Desktop (if created above), or"
echo "    • Run the installer directly: bash $INSTALL_SCRIPT"
echo ""
echo -e "  ${BOLD}To unregister the desktop entry:${RESET}"
echo "    rm $DESKTOP_DEST"
if [[ -f "$DESKTOP_HOME/MSGFLEX.desktop" ]] 2>/dev/null; then
echo "    rm $DESKTOP_HOME/MSGFLEX.desktop"
fi
echo ""
