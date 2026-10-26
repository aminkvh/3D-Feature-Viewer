# Plugin command reference

Install and first steps: [README.md](README.md). Open the panel with `ufv_gui`, or script it:


```
ufv_load P35498            # download AlphaFold model + annotate
ufv_ptms                   # PTM spheres        (also: ufv_variants, ufv_sites, ufv_ligands)
ufv_alphamissense          # colour cartoon     (also: ufv_domains, ufv_hotspots, ufv_pockets…)
ufv_report 1378            # residue report + zoom in
ufv_csv out.csv            # per-residue matrix
ufv_clear
```

More commands (all take an optional object / accession; the panel has the same controls):

| Command | What it does |
|---|---|
| `ufv_nearby 8` | Nearby radius, 3–15 Å (default 5): focus sticks, the report's Nearby list (+ other chains) and ligand pockets |
| `ufv_pockets_all` | Every ProtVar (AutoSite) pocket: scans all residues once (cached), colours each pocket, lists them |
| `ufv_show_pocket obj, pocket=67` | One pocket's lining as sticks + surface (also `ufv_show_pocket obj, 12+15+40`) |
| `ufv_motif_search obj, pocket=67` | RCSB Structure Motif Search for a pocket (PDB entries and AlphaFold/AlphaFill models) |
| `ufv_binding 108` / `ufv_pdbekb` | PDBe-KB experimental ligands / PPI interfaces at a residue / as a sphere layer |
| `ufv_tractability` | Open Targets tractability (small molecule, antibody, PROTAC) |
| `ufv_ligand_identity 50` | AlphaFill: show transplants with ≥ 50 % donor identity (presets 25–70) |
| `ufv_save_bfactor out.pdb, obj, hotspots` | PDB with an analysis in the B-factor column (`alphamissense`, `hotspots`, `hubs`, `burden`, `pocket`, `plddt`, `selected`) |
| `ufv_select` | PyMOL selection of the residues shown in 3D |
| `ufv_show_annotations obj, variants, search=pathogenic` | Tick annotation rows from a script (`all`, row ids, or `search=`; add `add` to keep existing ticks) |
| `ufv_reset_annotations` | Untick every annotation in every category (the panel's ↺) |
| `ufv_load_annotations my.csv` | Your own residue annotations (CSV / TSV: position, label, value, color, note — the extension's format) as the *Custom* category; the panel's **+**. `ufv_unload_annotations` removes them |
| `ufv_align` | Superpose every loaded structure of the protein on the first one, fitted on Cα atoms matched by UniProt position (falls back to sequence alignment below 10 shared residues) |

The residue report adds CADD to the per-substitution predictor table, one block per ProtVar
pocket, PDBe-KB experimental binding, analysis flags (● hotspot / burden / hub / pocket, once
computed) and the Nearby list. The ligand panel shows molecular weight, H-bond donors/acceptors,
PubChem / DrugBank / 2D-3D similarity / ChEMBL links, the pocket's mean pLDDT or B-factor,
functional sites, residue character and the closest-contact clash band (AlphaFill: donor identity
and clash score). The CSV's `nearby_ligands` column lists ligands within 5 Å of each residue.

## Residue numbering

Loaded residue numbers rarely equal UniProt positions. Pick one (panel → *Numbering*, or
`ufv_map` / `ufv_chain`):

- **Identity** — `resi == UniProt position` (AlphaFold / UniProt-numbered models).
- **SIFTS** — map through PDBe/SIFTS for a PDB id (handles gaps and author numbering).
- **Manual** — anchor each chain by hand; ideal for trajectories with custom numbering.

If SIFTS has no residue mapping for an entry, only the chains whose author-numbered sequence
matches the UniProt entry are mapped (with the offset found); otherwise nothing is projected and
a warning tells you to set the numbering by hand.

## Trajectories (MD, NMR ensembles)

Load a multi-state object (`load top.pdb, traj` + `load_traj run.dcd, traj`, or a multi-model
PDB), map it (`ufv_map traj, P35498, identity` / `ufv_chainmap …`) and use the layers as usual:
they follow every frame, and focus / Nearby / the report use the frame shown. Hotspots, hubs and
pockets are computed on the current frame and cached per frame (`ufv_traj_tiers traj, pockets, 250`).
For big trajectories the plugin switches PyMOL to `defer_builds_mode 3` once, so recolouring stays
fast (revert with `set defer_builds_mode, 0`).

| Command | What it does |
|---|---|
| `ufv_traj_frames 1, 0, 5, 0` | analysis window: start, stop (0 = last), stride, reference (0 = first) |
| `ufv_traj_rmsf` | Cα RMSD per frame + per-residue RMSF; cartoon blue → red; mean RMSF of PTM / pathogenic / site / pocket / hotspot residues vs the rest |
| `ufv_traj_dist traj, auto` | distance time series: `120~245` (Cα), `120~lig:ATP` (closest atoms), `10+12+15~200` (centroids); `auto` = PTM ↔ nearest pathogenic variant + ligand ↔ lining; pairs changing > 2 Å flagged and drawn in ONE measurement object |
| `ufv_traj_pockets traj, auto` | pocket size / Rg / buriedness per frame, opening–closing events, lining persistence (constraint-pocket clusters, ProtVar pockets, ligand pockets, or `120+124+156`) |
| `ufv_traj_contacts traj, ligand` / `partner` | contact occupancy (heavy atoms ≤ 4.5 Å) per residue, cartoon coloured by occupancy |
| `ufv_traj_annot` / `ufv_traj_nearby traj, 1378` | annotated residues' ligand / partner contact occupancy; Nearby occupancy of one residue |
| `ufv_traj_csv prefix` / `ufv_save_bfactor out.pdb, traj, rmsf` | per-frame + per-residue CSV; RMSF / `ligocc` / `partnerocc` in the B-factor column |

The panel's *Trajectory* box (shown for multi-state objects) runs the same analyses off the UI
thread, plots them (click or drag on the plot to jump to a frame; a red cursor follows the frame)
and lists residues / pairs / pockets (click a row to focus it or draw that distance). Timings for
500 frames of a 16 k-atom model: RMSF 0.7 s, distances 0.1 s + drawing, pockets 0.3 s, Nearby
occupancy 0.5 s.

## Cache

Annotation and metadata downloads, structure files (≤ 64 MB each) and ProtVar scores are cached in
`~/.ufv_cache` for 7 days, which makes reloads 5–12× faster; delete the folder any time to free the
disk space (override the folder
with `UFV_CACHE_DIR`, the lifetime with `UFV_CACHE_TTL_DAYS`, or disable it with
`UFV_NO_CACHE=1`). `ufv_fetch <uid>` always downloads a fresh copy. When the network fails, an
expired cached copy is used instead of nothing; `UFV_OFFLINE=1` never touches the network (cache only).
If 3D-Beacons cannot be reached, the structure list says so in one line (console and panel status):
computed models (SWISS-MODEL, ModelArchive, AlphaFill) are then missing until the next fetch.

## Commands (PyMOL and VMD)

Both plugins have the same command names, and they give the same numbers. The backend
(`ufv_core.py`) is shared. The suites check parity for groups, tiers, CSV, B-factor PDBs and
trajectory series.

The argument syntax differs. PyMOL takes comma-separated arguments with the object first
(`ufv_ptms myobj, P35498`). VMD takes Tcl words, and the molecule id is optional: it defaults to
the current molecule, which is the last one loaded or the one chosen with `ufv_use`.

| Group | Commands |
|---|---|
| Data, loading | `ufv_fetch` `ufv_structures` `ufv_load` `ufv_use` `ufv_align` `ufv_info` (VMD also `ufv_load_structure <key>`) |
| Numbering | `ufv_map` `ufv_chain` `ufv_chainmap` `ufv_custom` `ufv_chimera` |
| Marker layers | `ufv_load_annotations` `ufv_unload_annotations` `ufv_ptms` `ufv_variants` `ufv_diseases` `ufv_sites` `ufv_mutagenesis` `ufv_ligands` `ufv_ligands_hide` `ufv_partners` `ufv_partner_layers` `ufv_pdbekb` `ufv_prox` `ufv_hide` `ufv_select` `ufv_clear` (VMD also `ufv_show`, `ufv_layers`) |
| Cartoon colourings | `ufv_domains` `ufv_topology` `ufv_alphamissense` `ufv_plddt` `ufv_bfactor` `ufv_burden` `ufv_hotspots` `ufv_contacthubs` `ufv_pockets` `ufv_pocket_q` `ufv_pockets_all` (VMD also `ufv_cartoon plain`) |
| Residue, ligand | `ufv_report` `ufv_focus` `ufv_resetview` `ufv_nearby` `ufv_ligand_focus` `ufv_show_pocket` `ufv_motif_search` `ufv_binding` `ufv_tractability` `ufv_ligand_identity` (VMD also `ufv_residue`, `ufv_ligand_info`, `ufv_open_uniprot`) |
| Export | `ufv_csv` `ufv_save_bfactor` |
| Trajectory | `ufv_traj` `ufv_traj_frames` `ufv_traj_rmsf` `ufv_traj_dist` `ufv_traj_pockets` `ufv_traj_contacts` `ufv_traj_annot` `ufv_traj_nearby` `ufv_traj_tiers` `ufv_traj_csv` (VMD also `ufv_traj_bfactor`, `ufv_traj_show`, `ufv_traj_engine`, `ufv_traj_cancel`, `ufv_traj_gui`) |
| Panel, backend | `ufv_gui` (VMD also `ufv_python`, `ufv_core`, `ufv_backend`) |

### Annotations in VMD

- **Layers.** Each annotation is a layer (`ufv_ptms`, `ufv_variants`, and so on). `ufv_hide` and
  `ufv_show` toggle a layer instantly.
- **Cartoon colourings** are exclusive, and switching between them is instant. `ufv_cartoon plain`
  removes them.
- **`ufv_clear`** removes only the plugin's representations.
- **Your own files.** To annotate a molecule you loaded yourself, use `ufv_use <molid> <uid>`. If
  the file is a local copy (or an MD topology) of a listed structure, add its key, for example
  `ufv_use 0 P35498 AF-P35498`. That structure's entry then applies, including the AlphaFold PAE
  gate for constraint pockets and the motif search.
- **Colours.** The plugin defines VMD colour ids 17-32 once and maps other colours to the nearest
  defined one.

### Trajectories in VMD

Load a topology and its frames (`mol new top.pdb` + `mol addfile run.dcd waitfor all`), run
`ufv_use <molid> <uid> ?key?`, then use the commands above. The analyses match the PyMOL plugin's
definitions and numbers:

- `ufv_traj_frames 0 end 5 first`: sets the window (start, stop, stride, reference). Frames are
  0-based.
- `ufv_traj_rmsf`, `ufv_traj_dist auto`, `ufv_traj_pockets auto`, `ufv_traj_contacts ligand|partner`,
  `ufv_traj_annot`, `ufv_traj_nearby 1378`: the analyses.
- `ufv_traj_tiers pockets 250`: hotspots, hubs or constraint pockets on one frame, cached per frame.
- `ufv_traj_csv <prefix>` and `ufv_traj_bfactor out.pdb rmsf`: export.

Analyses run as background jobs, so VMD stays responsive, and `ufv_traj_cancel` stops one.
`ufv_traj_engine auto|native|backend` chooses the engine. `native` uses VMD's `measure`; `backend`
uses the shared numpy code. `auto` picks the faster engine for each analysis.

Per-frame values are written to the `user` field. RMSF and occupancy colour the cartoon, and
distance pairs are drawn as labels that follow the frames. The panel's *Trajectory* section
(or `ufv_traj_gui`) appears for multi-frame molecules. It plots the series (click to jump to a
frame) and lists residues, pairs and pockets.

Timings for 500 frames of a 16 k-atom model:

| Analysis | Time |
|---|---|
| RMSF | 0.25 s |
| Distances | 0.04-0.07 s |
| Pockets | 0.2-0.3 s |
| Nearby occupancy | 0.34 s |

---

`ufv_core.py` also works on its own:

- **Request server.** `python ufv_core.py --serve` reads JSON lines (`{"id", "op", "args"}`) and
  replies with `{"id", "ok", "result"}`. The VMD plugin uses it.
- **One-shot CLI.** For example `--emit-tcl <uid>`, `--structures`, `--residue`, `--csv` and
  `--download <uid>`.

It runs with any Python ≥ 3.8 and never imports PyMOL.
