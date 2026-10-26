#!/bin/sh
# Installs the 3D Feature Viewer plugin for VMD (Linux / macOS). Run from the unzipped folder:
#   sh install-vmd.sh              install / update
#   sh install-vmd.sh --uninstall  remove
# Copies the plugin to ~/3d-feature-viewer and adds a block to ~/.vmdrc (backed up first).
set -e
src=$(cd "$(dirname "$0")" && pwd)
dst="$HOME/3d-feature-viewer"
rc="$HOME/.vmdrc"
begin="# >>> 3D Feature Viewer >>>"
end="# <<< 3D Feature Viewer <<<"

strip() {  # print stdin without this plugin's block (also the block written by 2.1.0)
    awk -v b1="$begin" -v e1="$end" -v b2="# >>> 3D Feature Viewer for UniProt >>>" -v e2="# <<< 3D Feature Viewer for UniProt <<<" '
        $0 == b1 || $0 == b2 { skip = 1; next }
        $0 == e1 || $0 == e2 { skip = 0; next }
        !skip { print }'
}

[ -f "$rc" ] && { [ -f "$rc.3dfv-backup" ] || cp "$rc" "$rc.3dfv-backup"; }
old=""; [ -f "$rc" ] && old=$(strip < "$rc")

rm -rf "$dst" "$HOME/ufv-vmd"
if [ "$1" = "--uninstall" ]; then
    if [ -n "$old" ]; then printf '%s\n' "$old" > "$rc"; else rm -f "$rc"; fi
    echo "Removed. Restart VMD."
    exit 0
fi

[ -f "$src/vmd/ufv_vmd.tcl" ] || { echo "Run this script from the unzipped vmd-plugin folder." >&2; exit 1; }
mkdir -p "$dst/vmd"
cp "$src/ufv_core.py" "$dst/"
cp "$src"/vmd/*.tcl "$dst/vmd/"
{
    [ -n "$old" ] && printf '%s\n' "$old"
    cat <<BLOCK
$begin
set ::ufv_dir {$dst}
if {![info exists env(UFV_NO_AUTOLOAD)] && [file exists \$::ufv_dir/vmd/ufv_vmd.tcl] && [catch { source \$::ufv_dir/vmd/ufv_vmd.tcl } e]} { puts "3D Feature Viewer: \$e" }
catch { vmd_install_extension ufv ufv_gui "Analysis/3D Feature Viewer" }
$end
BLOCK
} > "$rc.new"
mv "$rc.new" "$rc"
echo "Done. Restart VMD, then open Extensions > Analysis > 3D Feature Viewer."
