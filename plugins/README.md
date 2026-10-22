# PyMOL and VMD plugins

The extension's annotations, predictors and analyses on structures and trajectories you open in
PyMOL or VMD: PTMs, variants, AlphaMissense, sites, domains, ligands, pockets and more.

Download both from the [Releases page](https://github.com/aminkvh/3D-Feature-Viewer-for-UniProt/releases): `pymol-plugin-v<version>.zip` and `vmd-plugin-v<version>.zip`.

## PyMOL

1. *Plugin → Plugin Manager → Install New Plugin → Choose file…* and pick `pymol-plugin-v<version>.zip` (don't unzip it).
2. Restart PyMOL and open the panel: *Plugin → 3D Feature Viewer for UniProt*, or type `ufv_gui`.

If you installed an older single-file `ufv_pymol.py`, uninstall it first (*Installed Plugins* tab).

## VMD

Requires VMD 1.9.4+ and Python 3.8+ on the system (numpy optional).

**Windows:** unzip `vmd-plugin-v<version>.zip`, then in that folder run

```powershell
powershell -ExecutionPolicy Bypass -File install-vmd.ps1
```

Start VMD from the new Start-menu shortcut **VMD (3D Feature Viewer)** and open
*Extensions → Analysis → 3D Feature Viewer for UniProt*.

**macOS / Linux:** unzip to `~/ufv-vmd` and add this line to `~/.vmdrc`:

```tcl
source ~/ufv-vmd/vmd/ufv_vmd.tcl
```

If Python isn't found, set it in the VMD console: `ufv_python /path/to/python3`.

## Quick start

1. Enter a UniProt accession → **Fetch** → pick a structure → **Load** (or annotate a structure you already have open).
2. **Annotations:** pick a category, tick rows to show them in 3D. Type in the search box to filter (`pathogenic`, `clinvar`, `1300-1400`, a disease name). ↺ unticks everything.
3. Click a row or residue for its report; double-click to zoom.
4. **Colour** the cartoon by domains, AlphaMissense, hotspots, pockets, and more.
5. **Export** a per-residue CSV or a PDB with scores in the B-factor column.

Every action is also a console command, and both plugins use the same commands and give the
same results. See [REFERENCE.md](REFERENCE.md).
