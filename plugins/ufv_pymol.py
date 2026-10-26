"""
3D Feature Viewer — PyMOL plugin
============================================

Brings UniProt residue-level annotations onto a structure or trajectory already
open in PyMOL, or downloads the AlphaFold model directly, and projects PTMs,
disease variants, ClinVar/AlphaMissense, functional sites, domains/regions, and
membrane topology onto it — the desktop counterpart of the browser extension.

It pulls from the same public resources (UniProtKB, EBI Proteins API, AlphaFold DB,
PDBe/SIFTS) with no API keys and only the Python standard library.

Quick start (PyMOL command line)
--------------------------------
    run ufv_pymol.py                 # load the plugin (ufv_core.py must sit next to it;
                                     # or install the ufv plugin .zip via Plugin Manager)

    ufv_gui                          # open the control panel (Fetch -> pick structure -> layers)
    ufv_load P35498                  # or just quick-load the AlphaFold model (no auto-projection)

Annotate a structure / trajectory YOU already loaded
----------------------------------------------------
    load mytraj.pdb, traj            # your object (single PDB or a trajectory)
    ufv_fetch P35498                 # fetch the annotations once

    # Tell the plugin how `traj` residue numbers relate to UniProt positions:
    ufv_map traj, P35498, identity              # resi == UniProt position (AlphaFold-style)
    ufv_map traj, P35498, sifts, 7dtd           # map through PDBe/SIFTS for PDB 7DTD
    # ...or define each chain by hand (great for trajectories with custom numbering):
    ufv_chain traj, A, 1                          # chain A resi 1  == UniProt 1
    ufv_chain traj, A, 200, 5, 480                # chain A resi 5  == UniProt 200, valid 5..480
    # ...or per chain, like the extension's "Custom PDB" (auto | offset | adopt | ignore + accession):
    ufv_custom traj, P35498                       # detect every chain (sequence match), print table
    ufv_chainmap traj, B, offset, 28              # chain B's first residue == UniProt 28
    ufv_chainmap traj, C, adopt, 7DTD_A/A         # copy a loaded structure's mapping
    ufv_chainmap traj, D, ignore                  # no annotations on chain D
    ufv_chainmap traj, E, auto, , P18505          # chain E is another protein (partner)
    ufv_structures P35498 / ufv_use <key>         # AlphaFold, isoforms (VSP-mapped), PDB entries
                                                  # (chains merged, SIFTS numbering), 3D-Beacons models
    ufv_partners obj / ufv_chimera obj            # other proteins in the entry / unmapped residues

    ufv_ptms traj            # PTM Ca spheres, coloured by category
    ufv_variants traj        # variant spheres, coloured by clinical consequence
    ufv_sites traj           # active / binding / metal sites (amber)
    ufv_domains traj         # colour the cartoon by domain / region / repeat
    ufv_topology traj        # colour by membrane topology
    ufv_alphamissense traj   # colour by mean AlphaMissense pathogenicity
    ufv_variants traj, , all, clinvar, AD1;CAA, disease   # filters + colour profile (consequence|disease|provenance|frequency|substitution)
    ufv_diseases             # disease summary (palette colours, counts)
    ufv_hotspots / ufv_contacthubs / ufv_pockets traj      # extension-exact analyses (mulberry32, all-source
    ufv_pocket_q 0.2         #   Brandes, PAE-gated pockets); pocket sensitivity recolours from cache
    ufv_prox 8, 12           # PTM / pathogenic-variant proximity radii for the residue report + CSV
    ufv_partner_layers traj  # partner proteins' variants / PTMs / sites on their own chains
    ufv_nearby 8             # Nearby radius 3-15 Å (focus sticks, report Nearby list, ligand pocket)
    ufv_pockets_all traj     # every ProtVar pocket (one crawl, cached), palette-coloured; then
    ufv_show_pocket traj, pocket=67 / ufv_motif_search traj, pocket=67   # surface / RCSB motif search
    ufv_binding 108 / ufv_pdbekb traj / ufv_tractability   # PDBe-KB binding, interfaces; Open Targets
    ufv_ligand_identity 50   # AlphaFill: transplants with >= 50 % donor identity only
    ufv_save_bfactor out.pdb, traj, hotspots   # analysis in the B-factor column; ufv_select = selection
    ufv_show_annotations traj, variants, search=cons:pathogenic src:clinvar 1300-1500
                             # check (= show) exactly these annotation-list rows; ids = all | none | row ids
    ufv_reset_annotations traj   # uncheck every annotation row (colouring / camera untouched)
    ufv_clear traj           # remove all UFV colouring/representations

All commands default to the active object when the object name is omitted.
"""

from __future__ import annotations

import contextlib
import functools
import json
import math
import os
import re
import sys
import threading
import time
import types

# PyMOL / Qt layer only. Everything PyMOL-independent (fetching, parsing, numbering maps, analyses,
# CSV / report formatting, the VMD command-line backend) lives in ufv_core.py next to this file.
try:
    from pymol import cmd
    _HAS_PYMOL = True
except Exception:  # not inside PyMOL: the module still imports, commands are not registered
    cmd = None
    _HAS_PYMOL = False

# ---- backend -------------------------------------------------------------------------------------
# ufv_core is imported ONCE and its names are re-exported here, so the code below (and `ufv.<name>`
# from scripts / tests) sees one shared backend state: caches, UFV_MAPS, _STATE and the analysis
# settings are the very same objects in both modules. Assigning a backend name on this module
# (ufv.HOTSPOT_RNG = "python", a monkeypatched fetcher) is forwarded to ufv_core (see the end).
try:
    from . import ufv_core as _core         # installed as the `ufv` package (Plugin Manager .zip)
except ImportError:                         # plain file next to ufv_core.py (`run ufv_pymol.py`, tests)
    _here = os.path.dirname(os.path.abspath(globals().get("__script__") or __file__))
    if _here not in sys.path:
        sys.path.insert(0, _here)
    import ufv_core as _core
for _n, _v in vars(_core).items():
    if not (_n.startswith("__") and _n.endswith("__")):
        globals()[_n] = _v
del _n, _v


# ----------------------------------------------------------------------------------------------
# Endpoints (mirror api.js)
# ----------------------------------------------------------------------------------------------
CHIMERA_COLOR = "#8fa5c0"     # residues of OUR chain with no UniProt mapping (chimeric partner layer)
# Ions / small inorganics drawn as spheres in the ligand layer (viewer-molstar.js ION_CODES).
ION_CODES = ("NA", "K", "LI", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "FE", "FE2", "MN", "MN3", "CU", "CU1",
             "CO", "NI", "CD", "HG", "CL", "BR", "IOD", "FLO", "F", "AL", "PT", "AU", "AG", "PB", "SO4", "PO4")


# ----------------------------------------------------------------------------------------------
# Residue-number mapping (UniProt position -> object chain+resi)
# ----------------------------------------------------------------------------------------------


_AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
        "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C"}


def _chain_residues(obj, chain):
    """[(author resv, one-letter)] of `obj`'s modelled polymer residues on `chain` (Cα present;
    insertion-code residues skipped)."""
    rows = []
    sel = "(%s) and polymer and name CA" % obj + (" and chain %s" % chain if chain else "")
    try:
        cmd.iterate(sel, "rows.append((resi, resv, resn))", space={"rows": rows})
    except Exception:
        return []
    out, seen = [], set()
    for resi, resv, resn in rows:
        if _resi_int(resi) is None or resv in seen:
            continue
        seen.add(resv)
        out.append((resv, _AA3.get(str(resn).upper(), "X")))
    return out


def _verify_segments(obj, segs, seq):
    """Check the segments that can't be trusted blindly against the model's own sequence: "label"
    (SIFTS gave only SEQRES numbers and no residue_listing was available) and "unverified" (no SIFTS entry
    for this accession; best_structures SEQRES ranges translated through residue_listing). A chain whose
    mapped residues match the UniProt sequence (>= 90 %) is kept as is; otherwise its author offset is
    re-derived from the sequence (_seq_offset); failing that, label segments are kept with a warning and
    unverified ones are dropped unless at least half of their mapped residues match."""
    for ch in list(segs):
        sl = segs[ch]
        if not any(s.get("label") or s.get("unverified") for s in sl):
            continue
        res = _chain_residues(obj, ch)
        tot, match = _segs_identity(res, seq, sl)
        if tot >= 5 and match >= 0.9 * tot:
            for s in sl:
                s.pop("label", None)
                s.pop("unverified", None)
            continue
        hit = _seq_offset(res, seq)
        if hit:
            off = hit[0]
            for s in sl:
                for k in ("label", "unverified", "s2a", "a2s", "seqres_start"):
                    s.pop(k, None)
                s["pdb_start"], s["pdb_end"] = s["u_start"] + off, s["u_end"] + off
                s["via"] = "sequence"
            print("[UFV] %s chain %s: PDBe numbering did not match the model; re-anchored from the sequence "
                  "(author resi = UniProt %+d, %.0f%% identity)." % (obj, ch, off, hit[1] * 100))
            continue
        if any(s.get("unverified") for s in sl) and not (tot and match >= 0.5 * tot):
            segs.pop(ch)
            print("[UFV] WARNING: %s chain %s: residue_listing mapping does not match the UniProt sequence "
                  "(%d/%d identical) — chain left unmapped." % (obj, ch, match, tot))
        else:
            for s in sl:
                s.pop("unverified", None)
            print("[UFV] WARNING: PDBe gives only label numbering for chain %s of %s and it could not be "
                  "verified against the sequence; positions on that chain may be offset." % (ch, obj))
    return segs


def _set_sifts_map(obj, pdb_id, uid, segs=None, hints=None, meta=None):
    """Map `obj` (PDB entry `pdb_id`) to `uid` through PDBe/SIFTS (+ residue_listing). `segs` / `hints` /
    `meta` may be passed in when already fetched off the UI thread (segs={} means "PDBe has nothing").
    When PDBe has no residue mapping, the accession is NOT sprayed over every chain with identity numbering
    (that projected e.g. GABA-A α1 annotations onto β/γ subunits and nanobodies): only chains named by PDBe
    best_structures / molecules — or, failing that, any chain — that match the UniProt sequence by their
    author numbering are mapped, with the offset found; otherwise nothing is mapped + warning. The entry's
    chimeric flag and partner chains (other accessions) are stored on the map."""
    uid = uid.upper()
    if segs is None:
        segs = _sifts_segments(pdb_id, uid)
    ann = _CACHE.get(uid) or fetch_annotations(uid)
    seq = ann.get("sequence", "")
    obj_chains = list(cmd.get_chains(obj) or [])
    ok = True
    if segs:
        segs = _verify_segments(obj, segs, seq)
    if segs:
        mp = {"mode": "sifts", "segments": segs, "status": _segments_status(segs)}
        print("[UFV] %s mapped for %s (%s): chains %s" % (obj, pdb_id.upper(), mp["status"],
                                                          ", ".join(sorted(segs))))
    else:
        if hints is None:
            hints = _pdb_chain_hints(pdb_id, uid, (ann.get("protnlm") or {}).get("name"))
        cands = [c for c in hints if c in obj_chains] if hints else obj_chains
        found, offs = {}, []
        for ch in cands:
            res = _chain_residues(obj, ch)
            hit = _seq_offset(res, seq)
            if hit:
                found[ch] = [_offset_segment(res, seq, hit[0])]
                offs.append("%s%+d" % (ch, hit[0]) if hit[0] else ch)
        mp = {"mode": "sifts", "segments": found,
              "status": "Mapped by sequence (no PDBe residue mapping)." if found else "Residue mapping unavailable."}
        if found:
            print("[UFV] no SIFTS residue mapping for %s / %s — mapped chain(s) %s by matching their "
                  "author-numbered sequence to %s (other chains left unmapped)."
                  % (pdb_id.upper(), uid, ", ".join(offs), uid))
        else:
            ok = False
            print("[UFV] WARNING: no SIFTS mapping for %s / %s and no chain%s matched the UniProt sequence — "
                  "nothing will be projected. Set numbering by hand with ufv_chainmap / ufv_chain / ufv_map."
                  % (pdb_id.upper(), uid, (" among " + ", ".join(cands)) if hints else ""))
    if meta is None:
        try:
            meta = _pdb_entry_meta(pdb_id, uid, sorted(mp["segments"]))
        except Exception:
            meta = None
    mp["pdbId"] = pdb_id.upper()
    mp["uid"] = uid
    if meta:
        for k in ("chimeric", "chimericChains", "organisms", "partners", "chainNames"):
            mp[k] = meta.get(k)
        if meta.get("chimeric"):
            print("[UFV] ⚛ %s: chain %s is chimeric (%s) — residues outside the %s mapping stay unannotated; "
                  "'ufv_chimera %s' highlights them." % (pdb_id.upper(), ",".join(meta.get("chimericChains") or []),
                                                         " / ".join(meta.get("organisms") or []) or "fusion",
                                                         uid, obj))
    UFV_MAPS[obj] = mp
    return ok


def _mapped_chains(obj):
    mp = UFV_MAPS.get(obj)
    if not mp:
        return cmd.get_chains(obj)
    if mp["mode"] in ("sifts", "custom"):
        return list(mp["segments"].keys())
    if mp["mode"] == "manual":
        return list(mp["chains"].keys())
    ex = mp.get("exclude")
    if ex:                      # AlphaFill: transplant chains are not our protein
        return [c for c in (cmd.get_chains(obj) or [""]) if c not in ex] or [""]
    return cmd.get_chains(obj) or [""]


# ----------------------------------------------------------------------------------------------
# Custom-structure mapping (modal.js Custom PDB, L3311-3600): per-chain auto | offset | adopt | ignore,
# per-chain accession. Works on any loaded object (your own PDB, a trajectory, a listed structure).
# ----------------------------------------------------------------------------------------------
CUSTOM_MODES = ("auto", "offset", "adopt", "ignore")


def _chain_extents(obj):
    """{chain: [first resv, last resv, n residues]} over the object's polymer Cα atoms (parsePdbChains)."""
    rows = []
    try:
        cmd.iterate("(%s) and polymer and name CA" % obj, "rows.append((chain, resv))", space={"rows": rows})
    except Exception:
        return {}
    ext = {}
    for ch, rv in rows:
        e = ext.get(ch)
        if e is None:
            ext[ch] = [rv, rv, 1]
        else:
            if rv < e[0]:
                e[0] = rv
            if rv > e[1]:
                e[1] = rv
            e[2] += 1
    return ext


def _resolve_chain_config(obj, ch, cfg, uid, seq, base, extent):
    """(segments or None, note) for one chain of a Custom mapping (modal.js applyBtn handler)."""
    mode = (cfg.get("mode") or "auto").lower()
    acc = (cfg.get("accession") or uid).upper()
    if mode == "ignore":
        return None, "ignored"
    if not _same_acc(acc, uid):
        return None, "partner chain (%s)" % acc
    if mode == "offset":
        start = int(cfg.get("start") or 1)
        first, last = extent[0], extent[1]
        return [{"u_start": start, "u_end": start + (last - first), "pdb_start": first, "pdb_end": last}], \
            "offset: resi %d == UniProt %d" % (first, start)
    if mode == "adopt":
        src = cfg.get("src")
        if not src:
            return None, "adopt: no source structure given"
        smp = UFV_MAPS.get(src) or {"mode": "identity"}
        src_chain = cfg.get("src_chain") or _map_primary_chain(smp) or ch
        sl = _map_chain_segments(smp, src_chain, len(seq))
        if not sl:
            return None, "adopt: %s chain %s has no mapping" % (src, src_chain)
        return sl, "adopted from %s/%s" % (src, src_chain)
    # auto: keep a known mapping (SIFTS / isoform / manual), else detect the offset from the sequence
    if base and base.get("mode") not in (None, "identity", "custom"):
        sl = _map_chain_segments(base, ch, len(seq))
        if sl:
            return sl, "auto: kept %s mapping" % base.get("mode")
    res = _chain_residues(obj, ch)
    hit = _seq_offset(res, seq)
    if hit:
        return [_offset_segment(res, seq, hit[0])], \
            "auto: sequence match (resi = UniProt %+d, %.0f%% identity)" % (hit[0], hit[1] * 100)
    return None, "auto: no sequence match — set an offset or adopt a mapping"


def _apply_custom_map(obj, uid, configs, reset=False, quiet=False):
    """Apply per-chain Custom mapping `configs` {chain: {mode, start, src, src_chain, accession}} to `obj`
    for accession `uid`. Chains not named keep their current Custom config (or 'auto'). The previous
    non-custom map (SIFTS / isoform / manual / identity) is kept as the 'auto' base."""
    uid = uid.upper()
    ann = _CACHE.get(uid) or fetch_annotations(uid)
    seq = ann.get("sequence", "")
    prev = UFV_MAPS.get(obj) or {"mode": "identity"}
    base = prev.get("base") if prev.get("mode") == "custom" else prev
    old_cfg = {} if (reset or prev.get("mode") != "custom") else (prev.get("config") or {})
    ext = _chain_extents(obj)
    cfg_all = {ch: dict(old_cfg.get(ch) or {"mode": "auto"}) for ch in ext}
    for ch, c in (configs or {}).items():
        if ch not in ext:
            print("[UFV] %s has no polymer chain '%s' (chains: %s)." % (obj, ch, ", ".join(sorted(ext)) or "none"))
            continue
        cfg_all[ch].update({k: v for k, v in c.items() if v not in (None, "")})
        if c.get("mode") and c["mode"] != "adopt":
            cfg_all[ch].pop("src", None)
            cfg_all[ch].pop("src_chain", None)
    segs, notes, partners, ignored = {}, {}, {}, []
    for ch in sorted(ext):
        cfg = cfg_all[ch]
        sl, note = _resolve_chain_config(obj, ch, cfg, uid, seq, base, ext[ch])
        notes[ch] = note
        if sl:
            segs[ch] = sl
        elif (cfg.get("mode") or "auto") == "ignore":
            ignored.append(ch)
        elif cfg.get("accession") and not _same_acc(cfg["accession"], uid):
            partners[ch] = cfg["accession"].upper()
    mp = {"mode": "custom", "segments": segs, "config": cfg_all, "notes": notes, "ignored": ignored,
          "accessions": {ch: (cfg_all[ch].get("accession") or uid).upper() for ch in ext},
          "partnersCustom": partners, "base": base, "uid": uid,
          "status": "Custom mapping: %d of %d chain(s) mapped." % (len(segs), len(ext))}
    for k in ("pdbId", "chimeric", "chimericChains", "organisms", "partners", "chainNames", "isoform"):
        if k in (base or {}):
            mp[k] = base[k]
    mp["chainId"] = (base or {}).get("chainId") if (base or {}).get("chainId") in segs else \
        (sorted(segs)[0] if segs else None)
    UFV_MAPS[obj] = mp
    _STATE.setdefault("obj_uid", {})[obj] = uid
    _invalidate_obj(obj)
    if not quiet:
        print("[UFV] %s — per-chain mapping for %s:" % (obj, uid))
        print("    %-5s %-13s %-11s %-7s %s" % ("chain", "resi", "accession", "mode", "mapping"))
        for ch in sorted(ext):
            e = ext[ch]
            print("    %-5s %-13s %-11s %-7s %s" % (
                ch or "''", "%d-%d (%d)" % (e[0], e[1], e[2]), mp["accessions"][ch],
                cfg_all[ch].get("mode", "auto"),
                (_seg_text(segs[ch]) + "  [" + notes[ch] + "]") if ch in segs else notes[ch]))
    return mp


# ----------------------------------------------------------------------------------------------
# Projection helpers
# ----------------------------------------------------------------------------------------------
_COLOR_REGISTERED = set()
_COLOR_INDEX = {}         # colour name -> PyMOL colour index (stable once defined)


def _color_index(name):
    """cmd.get_color_index, memoised (an API-locked round trip per call otherwise)."""
    i = _COLOR_INDEX.get(name)
    if i is None:
        i = _COLOR_INDEX[name] = cmd.get_color_index(name)
    return i


def _color_name(hex_str):
    h = hex_str.lstrip("#")
    name = "ufv_%s" % h
    if name not in _COLOR_REGISTERED:   # set_color is cheap but this is called per colour per redraw
        cmd.set_color(name, [int(h[0:2], 16) / 255.0, int(h[2:4], 16) / 255.0, int(h[4:6], 16) / 255.0])
        _COLOR_REGISTERED.add(name)
    return name


def _shift_obj(obj, val):
    """For commands documented as `cmd [object,] value`: when only one positional was given and it
    isn't a loaded object, it is the value — return (None, it) so the current object is used."""
    if obj and val in (None, "") and str(obj) not in cmd.get_names("objects"):
        return None, obj
    return obj, val


def _resolve_object(obj):
    if obj:
        return _traj_sync(obj) and obj      # trajectories: overlay copies / build mode follow the state count
    objs = cmd.get_object_list()
    for o in objs:
        if o in UFV_MAPS:
            return _traj_sync(o) and o
    return objs[-1] if objs else None


def _collapse_resi(resis):
    """Collapse a sorted int list into PyMOL resi tokens with ranges: [1,2,3,5] -> '1-3+5'. Keeps
    selection strings short so PyMOL doesn't choke parsing thousands of '+'-separated residues
    (a 260-residue domain becomes one '100-360' token instead of 260)."""
    def tok(r):                     # negative residue numbers need the minus escaped: \-5-\-3
        return ("\\%d" % r) if r < 0 else str(r)
    out, i, n = [], 0, len(resis)
    while i < n:
        j = i
        while j + 1 < n and resis[j + 1] == resis[j] + 1:
            j += 1
        out.append(tok(resis[i]) if i == j else "%s-%s" % (tok(resis[i]), tok(resis[j])))
        i = j + 1
    return "+".join(out)


def _sel_for_positions(obj, positions, ca_only=False, target=None):
    """Build a compact PyMOL selection covering the given UniProt positions. Numbering is resolved
    with `obj`'s map; the selection is rooted at `target` (defaults to `obj`) so it can address a
    copied sphere layer object that shares the same chain/resi numbering."""
    per_chain = {}
    for ch in _mapped_chains(obj):
        resis = []
        for p in positions:
            r = _uni_to_resi(obj, ch, p)
            if r is not None:
                resis.append(int(r))
        if resis:
            per_chain[ch] = sorted(set(resis))
    if not per_chain:
        return None
    clauses = []
    for ch, resis in per_chain.items():
        rsel = "resi " + _collapse_resi(resis)
        clauses.append("(chain {} and {})".format(ch, rsel) if ch else "({})".format(rsel))
    sel = "({}) and ({})".format(target or obj, " or ".join(clauses))
    if ca_only:
        sel += " and name CA"
    return sel
_CARTOON_LAYER = {}   # obj -> active cartoon-colouring mode (informational)
POINT_TAGS = {"ptms": "ptm", "variants": "var", "sites": "site"}
CARTOON_TAGS = ("domains", "topology", "alphamissense", "burden", "plddt", "bfactor",
                "hotspots", "contacthubs", "pockets", "allpockets")


_UNDO_TLS = threading.local()
_NAME_HEALS = []      # (object, how) for every object selection re-registered by _heal_object_names
_HEAL_WARNED = set()  # objects whose lost name could not be re-registered (warned once)


_HEAL_TYPES = {}      # object name -> "object:molecule" (names last seen as molecules; others re-queried)
_HEAL = {"names": None, "full": 0.0, "suspect": set()}   # see _heal_check
HEAL_FULL_EVERY = 2.0  # s: a plugin command re-checks EVERY object at most this often (bursts pay once)


def _heal_check():
    """Which objects _heal_object_names must look at on this plugin command (undo active). The
    selector check costs ~1 ms per object in a 200 k-atom session (it scans every atom), which made it
    the bulk of each command's overhead there, so instead of every object every time:
      * every object, when nothing is recorded yet or the last full check is > HEAL_FULL_EVERY s old;
      * otherwise the objects the bug can hit: names whose deletion we SAW between two plugin commands
        (the trigger: a delete recorded by the undo manager, then a new object of that name — kept as
        'suspect' for the session), objects that are new since the last command, and every object when
        the order of the known ones changed (a delete + re-load between two commands moves the name).
    Our own deletes happen with the undo manager paused (safe), and the name list is recorded when each
    plugin command ends (_heal_note), so a difference at the next start is the user's doing."""
    try:
        names = cmd.get_names("objects") or []
    except Exception:
        return []
    last = _HEAL["names"]
    now = time.time()
    cur = set(names)
    lastset = set(last or ())
    for n in lastset - cur:
        _HEAL["suspect"].add(n)
        _HEAL_TYPES.pop(n, None)
    if last is None or now - _HEAL["full"] > HEAL_FULL_EVERY:
        _HEAL["full"] = now
        return names
    if [n for n in last if n in cur] != [n for n in names if n in lastset]:
        _HEAL["full"] = now
        return names
    return [n for n in names if n not in lastset or n in _HEAL["suspect"]]


def _heal_note():
    """End of a plugin command: remember the object list (the baseline _heal_check compares against)."""
    try:
        _HEAL["names"] = tuple(cmd.get_names("objects") or ())
    except Exception:
        _HEAL["names"] = None


def _heal_object_names(names=None):
    """Work around a PyMOL 3.x undo-manager bug (seen in Incentive 3.1.4; repro in the session notes):
    deleting an object while undo is recording keeps the old object in the undo history; when that
    history entry later ages out, PyMOL drops the selector entry BY NAME, so a newer, live object with
    the same name (a structure the user deleted and re-loaded, a layer object we rebuilt after the user
    deleted it) silently stops resolving: "Invalid selection name '4HHB_B_D'", although get_names()
    still lists it and "model 4HHB_B_D" still selects its atoms. Not name-specific (digit-leading or
    not), and only possible while the undo manager is active.

    Check: cmd.get_object_list(<bare object name>) looks the name up in the selector directly (no
    parse, no whole-session atom table; ~1 ms per 16 k-atom object) and returns None when the entry
    is gone. On a miss, re-register the object's own selection via cmd.sort (ObjectMoleculeSort ->
    invalidate(cRepInvAtoms) -> SelectorUpdateObjectSele; already-sorted atoms are not moved; if it
    did reorder, the object's derived caches are dropped). Call with the undo manager paused (_no_undo
    does). `names`: only these objects (default: every object). Non-molecular objects (measurements,
    maps, CGOs) are skipped without querying the selector (that printed a Selector-Error for each).
    Returns the healed names."""
    if cmd is None:
        return []
    if names is None:
        try:
            names = cmd.get_names("objects")
        except Exception:
            return []
    healed = []
    for o in names:
        try:
            if _HEAL_TYPES.get(o) != "object:molecule":
                t = cmd.get_type(o)          # executive lookup (cheap); cached for molecules only
                if t != "object:molecule":
                    continue
                _HEAL_TYPES[o] = t
            if cmd.get_object_list(o) is not None:
                continue
            if cmd.get_type(o) != "object:molecule":
                _HEAL_TYPES.pop(o, None)
                continue
        except Exception:
            continue
        before = []
        try:
            cmd.iterate("model %s" % o, "before.append(ID)", space={"before": before})
            if not before:
                continue                    # no atoms: nothing a selection could find anyway
            cmd.sort(o)
            ok = cmd.get_object_list(o) is not None
        except Exception:
            ok = False
        if not ok:
            if o not in _HEAL_WARNED:
                _HEAL_WARNED.add(o)
                print("[UFV] warning: PyMOL lost the selection name of object '%s' (PyMOL undo bug) and it "
                      "could not be re-registered; use 'model %s' or re-load it." % (o, o))
            continue
        after = []
        cmd.iterate("model %s" % o, "after.append(ID)", space={"after": after})
        how = "sort"
        if after != before:
            _invalidate_obj(o)
            how = "sort (atoms reordered; caches dropped)"
        _NAME_HEALS.append((o, how))
        healed.append(o)
    return healed


@contextlib.contextmanager
def _no_undo():
    """Pause PyMOL 3.x's undo manager for the duration. In PyMOL 3.x every cmd.set / cmd.unset is
    wrapped in an undo session that snapshots the WHOLE session forward and reverse (60-470 ms per
    call, growing with session size) — that, not the object count, was the stutter on every toggle
    and residue click. Pausing records nothing for our (rebuildable) overlay edits; the previous mode
    is restored afterwards so the user's own undo keeps working. Nesting-safe: an inner block sees the
    manager already paused (or disabled) and leaves it alone. No-op on PyMOL 2.x (no undo manager)."""
    depth = getattr(_UNDO_TLS, "depth", 0)
    prev = None
    if depth == 0 and cmd is not None and hasattr(cmd, "undo_get_mode"):
        try:
            mode = cmd.undo_get_mode()
            if mode in (1, 2):              # UndoMode.Enable / Unpause -> Pause (3)
                cmd.undo_mode(3)
                prev = mode
        except Exception:
            prev = None
        if prev is not None:                # undo active: an aged-out history entry may have dropped a
            chk = _heal_check()             # live object's selection name (PyMOL bug, see the helper)
            if chk:
                _heal_object_names(chk)
    _UNDO_TLS.depth = depth + 1             # inner blocks skip the (API-locked) mode query entirely
    try:
        if depth == 0:
            with _keep_view():              # plugin commands never move the camera unless they say so
                yield
        else:
            yield
    finally:
        _UNDO_TLS.depth = depth
        if depth == 0 and cmd is not None:
            _heal_note()
        if prev is not None:
            try:
                cmd.undo_mode(prev)
            except Exception:
                pass


_VIEW_TLS = threading.local()   # .saved = the camera view the outermost plugin command must end on
_FOCUS_VIEW = {"pre": None, "post": None}   # camera before the first focus of a chain / right after the last


@contextlib.contextmanager
def _keep_view():
    """Camera guard for plugin commands: the view (cmd.get_view()) is recorded when the OUTERMOST
    plugin command starts and put back exactly when it ends, so loading, layers, colourings, numbering,
    alignment, trajectory analyses, frame jumps, enable / disable ... can never nudge the camera (no
    auto_zoom, no clipping drift). Only explicit view actions move it: they call _view_commit() after
    their zoom / orient, which makes the new view the one to end on. Nesting-safe (inner blocks are
    no-ops). Main thread only: a worker thread must never put back a view the user has since changed."""
    if (cmd is None or getattr(_VIEW_TLS, "saved", None) is not None
            or threading.current_thread() is not threading.main_thread()):
        yield
        return
    try:
        v0 = tuple(cmd.get_view())
    except Exception:
        yield
        return
    _VIEW_TLS.saved = v0
    try:
        yield
    finally:
        target, _VIEW_TLS.saved = _VIEW_TLS.saved, None
        try:
            if target is not None and tuple(cmd.get_view()) != target:
                cmd.set_view(target)
        except Exception:
            pass


def _view_commit():
    """An explicit view action (focus / zoom / reset / first-load orient) just set the camera: keep it."""
    if getattr(_VIEW_TLS, "saved", None) is not None:
        try:
            _VIEW_TLS.saved = tuple(cmd.get_view())
        except Exception:
            pass


def _focus_zoom(sel, buffer=3):
    """Explicit focus zoom. Remembers the view before the FIRST focus of a chain (a later focus made
    without the user moving the camera keeps it), so clearing the focus (ufv_resetview) returns to
    exactly that view."""
    cur = tuple(cmd.get_view())
    if _FOCUS_VIEW["post"] is None or cur != _FOCUS_VIEW["post"]:
        _FOCUS_VIEW["pre"] = cur                # the user moved the camera (or first focus): start here
    cmd.zoom(sel, buffer, animate=0)
    _FOCUS_VIEW["post"] = tuple(cmd.get_view())
    _view_commit()


@contextlib.contextmanager
def _auto_zoom_off():
    """auto_zoom 0 around an internal load (zoom=0 is passed as well; this also covers PyMOL builds that
    zoom on the first object regardless). Call with the undo manager paused (cmd.set is undo-tracked)."""
    prev = 0
    try:
        prev = int(cmd.get_setting_int("auto_zoom"))
        if prev:
            cmd.set("auto_zoom", 0, quiet=1)
    except Exception:
        prev = 0
    try:
        yield
    finally:
        if prev:
            try:
                cmd.set("auto_zoom", prev, quiet=1)
            except Exception:
                pass


def _undo_paused(fn):
    """Decorator: run `fn` with the undo manager paused (see _no_undo). Used on every ufv_* command
    (they are both cmd.extend commands and what the GUI calls) and the layer helpers the GUI uses."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with _no_undo():
            return fn(*args, **kwargs)
    return wrapper


_LAYER_LISTENERS = []   # callables the open panel(s) register; run after any layer / cartoon change


def _notify_layers():
    """Tell open panels the layer bookkeeping changed (cheap: panels coalesce and only read dicts)."""
    for fn in list(_LAYER_LISTENERS):
        try:
            fn()
        except Exception:
            try:
                _LAYER_LISTENERS.remove(fn)       # panel was deleted
            except ValueError:
                pass


def _notifies(fn):
    """Decorator (apply under @_undo_paused): notify layer listeners after `fn` returns."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        finally:
            _notify_layers()
    return wrapper


@contextlib.contextmanager
def _batch():
    """Group many cmd ops: runs them with the undo manager paused. (The old global auto_zoom /
    suspend_updates toggles are gone: they were themselves undo-tracked cmd.set calls costing ~0.8 s
    per use under PyMOL 3.x. New objects are created with zoom=0 instead, so a layer never yanks
    the camera, and commands issued from one Qt slot never repaint in between anyway.)"""
    with _no_undo():
        yield


def _ufv_layer_obj(obj, tag):
    # Characters other than [A-Za-z0-9_] are hex-escaped (not all folded to "_") so two objects such
    # as "6HUO-A" and "6HUO_A" can never share — and delete — each other's layer objects.
    safe = re.sub(r"[^A-Za-z0-9_]", lambda m: "_%02x" % ord(m.group(0)), obj or "x")
    return "ufv_%s_%s" % (safe, tag)


# Every per-structure overlay object is ufv_<obj>_<tag> with one of these tags (the hidden Cα
# scaffold, the sphere layers, and the two ligand-layer objects). Deleting by this exact list instead
# of by name prefix keeps `ufv_clear X` from also deleting the layers of a structure named X_B.
_LAYER_TAGS = ("ca", "ptm", "var", "site", "mutag", "dom", "domr", "filt", "lig", "ion", "prt", "kb")


def _object_set():
    return set(cmd.get_object_list() or [])


def _enabled_set():
    try:
        return set(cmd.get_names("public_objects", enabled_only=1) or [])
    except Exception:
        return _object_set()


def _ca_scaffold(obj):
    """A persistent, hidden CA-only copy of `obj` (ufv_<obj>_ca), built once. Marker layers are then
    carved out of THIS small object instead of re-scanning the full structure (which is slow for big
    models, e.g. AlphaFill) on every toggle. Returns the scaffold name, or `obj` if it can't be built."""
    name = _ufv_layer_obj(obj, "ca")
    if name not in _object_set():
        try:
            cmd.create(name, "(%s) and name CA and polymer" % obj, zoom=0)
            cmd.hide("everything", name)
            cmd.disable(name)
        except Exception:
            return obj
    return name


def _tune_sphere_obj(lobj):
    # GPU impostor spheres render thousands of points cheaply.
    for k, v in (("sphere_mode", 9), ("sphere_quality", 1), ("sphere_scale", 0.8)):
        try:
            cmd.set(k, v, lobj)
        except Exception:
            pass


def _res_key_map(obj, pos_value):
    """{uniprot_pos: value} -> {(chain, resv): value} over the object's mapped chains, so an `alter`
    can colour atoms with a plain dict lookup on (chain, resv) instead of calling back into Python
    for every atom (27 ms vs 77 ms on 16 k atoms)."""
    chains = _mapped_chains(obj)
    m = {}
    for pos, val in pos_value.items():
        pos = int(pos)
        for ch in chains:
            r = _uni_to_resi(obj, ch, pos)
            if r is not None:
                m[(ch, int(r))] = val
    return m


@_undo_paused
@_notifies
def _set_sphere_layer(obj, tag, groups, universe=None):
    """Render a sphere layer on its OWN lightweight object — a copy of just the annotated Cα atoms,
    with no cartoon — so showing thousands of points never rebuilds the main structure's cartoon.
    groups = {color: [positions]} to show, or None to hide the layer.

    The layer object is built once and then only re-shown: hiding DISABLES it (instant), showing the
    same groups again just enables it, and showing a different subset of the positions it was built
    with (e.g. a variant consequence filter) is one `alter` (colour + a visibility flag in q) plus a
    show/hide on that small object — no delete/create per click. `universe` (optional positions) is
    included in the build so later subsets (other filters) stay on the fast path."""
    lobj = _ufv_layer_obj(obj, tag)
    st = _SPHERE_STATE.setdefault(obj, {})
    built = _LAYER_BUILT.setdefault(obj, {})
    exists = tag in built and lobj in _object_set()
    if not groups:
        st.pop(tag, None)
        if exists:
            cmd.disable(lobj)
        return 0
    allpos = {int(p) for positions in groups.values() for p in positions}
    if exists and built[tag]["groups"] == groups:      # unchanged: just show it again
        cmd.enable(lobj)
        st[tag] = groups
        return len(allpos)
    with _batch():
        if not (exists and allpos <= built[tag]["pos"]):
            want = allpos | {int(p) for p in (universe or ())}
            scaffold = _ca_scaffold(obj)                 # carve the layer out of the small CA-only copy
            sel = _sel_for_positions(obj, sorted(want), ca_only=True, target=scaffold)
            if not sel:
                st.pop(tag, None)
                if exists:
                    cmd.disable(lobj)
                return 0
            cmd.delete(lobj)
            cmd.create(lobj, sel, zoom=0)    # copy only the Cα atoms into a separate, cartoon-free object
            cmd.hide("everything", lobj)
            _tune_sphere_obj(lobj)
            built[tag] = {"pos": frozenset(want), "groups": None}
        # Colour + visibility for this subset in ONE alter pass over the (small) layer object: each
        # atom looks up (colour index, shown-flag) by (chain, resv); the flag goes in q (occupancy is
        # meaningless on this Cα copy), so the show/hide below needs no resi-list selection parsing.
        pos_val = {}
        for color, positions in groups.items():
            idx = _color_index(_color_name(color))
            for p in positions:
                pos_val[int(p)] = (idx, 1.0)
        cmd.alter(lobj, "color, q = _m.get((chain, resv), (_b, 0.0))",
                  space={"_m": _res_key_map(obj, pos_val), "_b": _color_index("gray80")})
        cmd.hide("spheres", lobj)
        cmd.show("spheres", "(%s) and q > 0.5" % lobj)
        cmd.recolor(lobj)
        cmd.enable(lobj)
    built[tag]["groups"] = groups
    st[tag] = groups
    return len(allpos)


@_undo_paused
def _hide_layer(obj, tag):
    if tag in _ANN_TAG_CAT:
        _ann_forget(obj, (_ANN_TAG_CAT[tag],))     # annotation rows of that category are unchecked too
    _set_sphere_layer(obj, tag, None)


def _layer_wanted(obj, lobj):
    """Should overlay object `lobj` of `obj` be visible (an active sphere layer, or the ligand layer
    while it is switched on)? Used when re-enabling objects a focus temporarily hid."""
    for tag in (_SPHERE_STATE.get(obj) or {}):
        if _ufv_layer_obj(obj, tag) == lobj:
            return True
    return obj in _LIG_SHOWN and lobj in (_ufv_layer_obj(obj, "lig"), _ufv_layer_obj(obj, "ion"))


@_undo_paused
@_notifies
def _delete_layers(obj):
    """Remove every overlay object belonging to `obj` (exact names, never by prefix)."""
    names = _object_set()
    for tag in _LAYER_TAGS:
        lobj = _ufv_layer_obj(obj, tag)
        if lobj in names:
            try:
                cmd.delete(lobj)
            except Exception:
                pass
    for tag in ("focusnb", "focusset", "focusnlig", "focuslig", "dist"):   # dist: not in get_object_list
        try:
            cmd.delete(_ufv_layer_obj(obj, tag))
        except Exception:
            pass
    _SPHERE_STATE.pop(obj, None)
    _LAYER_BUILT.pop(obj, None)
    _LIG_SHOWN.discard(obj)
    _LIG_ROWS.pop(obj, None)
    _ANN_LIGS.pop(obj, None)
    _ann_forget(obj)


@_undo_paused
@_notifies
def _reset_cartoon(obj, grey=True):
    """Drop any cartoon-colour layer: recolour neutral grey and remove domain single-residue spheres.
    Only the polymer is recoloured so bound ligands/cofactors keep their own colours (this is what made
    cartoon colouring grey out AlphaFill ligands). grey=False: the caller recolours EVERY polymer atom
    right after (its alter has a base colour), so the grey pass would only be one more full-object colour
    change for PyMOL to rebuild."""
    if grey:
        cmd.color("gray80", "(%s) and polymer" % obj)
    _set_sphere_layer(obj, "dom", None)
    _CARTOON_LAYER.pop(obj, None)


@_undo_paused
@_notifies
def _set_cartoon_layer(obj, name, grey=True):
    _reset_cartoon(obj, grey)
    cmd.show("cartoon", obj)
    _CARTOON_LAYER[obj] = name


# ----------------------------------------------------------------------------------------------
# Annotations v3: ONE checked-row state per (object, accession). A checked row of the annotation list
# (ufv_core.ann_rows: PTMs, variants, sites, mutagenesis, domains, ligands, partners, PDBe-KB) is drawn in
# 3D, an unchecked one is not. Every category is rebuilt from its checked set in ONE batched call on its
# one layer object (never per-row objects). Console commands (ufv_ptms, ufv_variants, ufv_hide ...) set
# the same state, so the panel list and the console always agree.
# ----------------------------------------------------------------------------------------------
_ANN_ON = {}        # (obj, uid) -> {cat: set(row ids)}
_ANN_STYLE = {}     # obj -> {"color": variant colour profile, "dis": disease set for 'disease' colours,
                    #         "prt": partner layer set, "kb": PDBe-KB kind}
_ANN_LIGS = {}      # obj -> {resn: copies}   (ligand rows; dropped when the object is (re)loaded / cleared)
_LIG_ROWS = {}      # obj -> set of ligand codes shown by the ligand layer (missing = every ligand)
ANN_TAGS = {"ptm": "ptm", "var": "var", "site": "site", "mutag": "mutag", "dom": "domr", "prt": "prt", "kb": "kb",
            "cus": "cus"}
_ANN_TAG_CAT = {v: k for k, v in ANN_TAGS.items()}


def _ann_key(obj, uid=None):
    return (obj, (_resolve_uid(uid, obj) or "").upper())


def _ann_state(obj, uid=None):
    """{cat: set(ids)} of (obj, uid). Layers are per object, so another accession's state of the same object
    is dropped."""
    key = _ann_key(obj, uid)
    st = _ANN_ON.get(key)
    if st is None:
        for k in [k for k in _ANN_ON if k[0] == obj]:
            _ANN_ON.pop(k, None)
        st = _ANN_ON[key] = {c: set() for c in ANN_CATS}
    return st


def _ann_forget(obj, cats=None):
    """Uncheck `cats` (default all) of every state of obj (layer hidden / deleted elsewhere)."""
    for k, st in _ANN_ON.items():
        if k[0] == obj:
            for c in (cats or ANN_CATS):
                st[c] = set()


def _ann_ligands(obj):
    m = _ANN_LIGS.get(obj)
    if m is None and obj:
        try:
            m = enumerate_ligands(obj)
        except Exception:
            m = {}
        _ANN_LIGS[obj] = m
    return m or {}


def _ann_extra(obj, uid, cats=("lig", "prt", "kb")):
    """Row data that is not in the annotation dict (no network: partner / PDBe-KB rows exist once fetched)."""
    ex = {}
    if "lig" in cats and obj:
        ex["lig"] = [{"resn": r, "count": n, "name": (_LIG_CACHE.get(r) or {}).get("name") or ""}
                     for r, n in sorted(_ann_ligands(obj).items())]
    if "prt" in cats and obj:
        pa = (UFV_MAPS.get(obj) or {}).get("_partnerAnn")
        if pa is not None:
            cols = _partner_sphere_colors(pa, set(PARTNER_LAYERS))
            ex["prt"] = [{"chain": k[0], "resi": k[1], "key": k, "color": c, "text": t} for k, (c, t) in cols.items()]
    if "kb" in cats and uid:
        kb = _PDBEKB_CACHE.get(_canon_acc(uid))
        if kb is not None or _canon_acc(uid) in _PDBEKB_CACHE:
            ex["kb"] = kb or {}
    return ex


def annotation_rows(obj, uid=None, cat="all"):
    """The annotation-list rows of `obj` (ufv_core.ann_rows + this object's ligands / partners / PDBe-KB)."""
    uid = (_resolve_uid(uid, obj) or "").upper()
    return ann_rows(_CACHE.get(uid) or {}, cat, _ann_extra(obj, uid))


def ann_available(obj, uid=None):
    """{cat: row count, or None = not fetched yet (partners / PDBe-KB: network on first use)}."""
    uid = (_resolve_uid(uid, obj) or "").upper()
    ann = _CACHE.get(uid) or {}
    ex = _ann_extra(obj, uid)
    out = {}
    for c in ANN_CATS:
        if c in ("prt", "kb") and c not in ex:
            out[c] = None
        else:
            out[c] = len(ann_rows(ann, c, ex))
    return out


def _ann_style(obj):
    return _ANN_STYLE.setdefault(obj, {"color": "consequence", "dis": None, "prt": set(PARTNER_LAYERS), "kb": "both"})


@_undo_paused
@_notifies
def _ann_apply(obj, uid=None, cats=None):
    """Draw categories `cats` (default all) of (obj, uid) from their checked rows — one layer object per
    category, one batched rebuild / alter each. Returns {cat: rows shown}."""
    uid = (_resolve_uid(uid, obj) or "").upper()
    st = _ann_state(obj, uid)
    ann = _CACHE.get(uid) or {}
    sty = _ann_style(obj)
    out = {}
    with _batch():
        for c in (cats or ANN_CATS):
            on = st.get(c) or set()
            out[c] = len(on)
            if c == "lig":
                _ann_apply_ligands(obj, {i.split(":", 1)[1] for i in on})
                continue
            if c == "prt":
                if not on:
                    _show_partner_layer(obj, None)
                    continue
                ex = _ann_extra(obj, uid, ("prt",)).get("prt") or []
                pa = (UFV_MAPS.get(obj) or {}).get("_partnerAnn") or []
                cols = _partner_sphere_colors(pa, sty["prt"] or set(PARTNER_LAYERS)) if pa else {}
                sel = {}
                for r in ex:
                    if "prt:%s:%s" % (r["chain"], r["resi"]) in on:
                        sel[r["key"]] = (cols.get(r["key"]) or (r["color"],))[0]
                _show_partner_layer(obj, sel or None)
                continue
            tag = ANN_TAGS[c]
            if not on:
                _set_sphere_layer(obj, tag, None)
                continue
            if c == "var":
                vs = ann.get("variants") or []
                chosen = [vs[int(i.split(":", 1)[1])] for i in sorted(on, key=lambda x: int(x.split(":", 1)[1]))
                          if int(i.split(":", 1)[1]) < len(vs)]
                groups = ann_variant_groups(ann, chosen, sty["color"], sty["dis"])
                _set_sphere_layer(obj, tag, groups or None, universe=_variant_positions(ann) if ann else None)
                continue
            rows = ann_rows(ann, c, _ann_extra(obj, uid, ("kb",)) if c == "kb" else None)
            groups, universe = {}, set()
            kb_kind = sty.get("kb") or "both"
            for r in rows:
                universe.update(r["pts"])
                if r["id"] in on:
                    col = r["color"]
                    if c == "kb" and kb_kind != "both":
                        col = KB_COLORS[kb_kind]
                    groups.setdefault(col, []).extend(r["pts"])
            if c == "kb":
                groups = {k: sorted(v) for k, v in groups.items()}
            _set_sphere_layer(obj, tag, groups or None, universe=sorted(universe))
    return out


def _ann_apply_ligands(obj, codes):
    """Ligand rows: build the ligand layer once (ufv_ligands), then show only the checked codes on it."""
    if not codes:
        _LIG_ROWS.pop(obj, None)
        if obj in _LIG_SHOWN or _lig_objs(obj):
            ufv_ligands_hide(obj, _quiet=True)
        return
    allc = set(_ann_ligands(obj))
    if codes >= allc:
        _LIG_ROWS.pop(obj, None)
    else:
        _LIG_ROWS[obj] = set(codes)
    ufv_ligands(obj, _quiet=True)


def _ann_set(obj, uid, cat, ids, mode="set"):
    st = _ann_state(obj, uid)
    ids = set(ids)
    if mode == "add":
        st[cat] |= ids
    elif mode == "remove":
        st[cat] -= ids
    else:
        st[cat] = ids


@_undo_paused
def ann_set_checked(obj, uid, changes):
    """changes = {cat: (ids, on)} -> one batched update of every touched category (the panel's check boxes,
    header check-all, Space on a multi-row selection). Returns {cat: rows shown}."""
    touched = []
    for cat, (ids, on) in changes.items():
        _ann_set(obj, uid, cat, ids, "add" if on else "remove")
        touched.append(cat)
    return _ann_apply(obj, uid, touched)


def ann_counts(obj, uid=None):
    """{cat: checked rows} of (obj, uid) — what is shown in 3D, per category."""
    st = _ANN_ON.get(_ann_key(obj, uid)) or {}
    return {c: len(st.get(c) or ()) for c in ANN_CATS}


def ann_checked(obj, uid=None):
    """{cat: set(ids)} (read-only view of the checked state)."""
    return _ANN_ON.get(_ann_key(obj, uid)) or {c: set() for c in ANN_CATS}


@_undo_paused
def ufv_reset_annotations(obj=None):
    """ufv_reset_annotations [object] — uncheck every annotation row of every category (PTMs, variants, sites,
    mutagenesis, domains, ligands, partners, PDBe-KB) in one batched update. The cartoon colouring and the
    camera are not touched (ufv_clear also resets the colouring and the focus)."""
    obj = _resolve_object(obj)
    if not obj:
        return 0
    n = sum(ann_counts(obj).values())
    _ann_forget(obj)
    _ann_state(obj)
    _ann_apply(obj, None, ANN_CATS)
    print("[UFV] %s: %d annotation row(s) unchecked in every category." % (obj, n))
    return n


@_undo_paused
def ufv_show_annotations(obj=None, category=None, ids="all", mode="set", search=None, uid=None):
    """ufv_show_annotations [object,] category [, ids [, mode [, search]]] — script the annotation list:
    check exactly (mode=set), add or remove rows of one category (ptms | variants | sites | mutagenesis |
    domains | ligands | partners | pdbekb | all). ids = all | none | row ids / positions separated by
    space, '+' or ';' (var:12 or 12 = row 12 of that category; ligand codes for ligands), or use
    search=<query> (the panel's search syntax, e.g. search=cons:pathogenic src:clinvar 1300-1500)."""
    if obj and category is None and ann_category(obj) and obj not in cmd.get_names("objects"):
        obj, category = None, obj
    obj = _resolve_object(obj)
    uid = (_resolve_uid(uid, obj) or "").upper()
    cat = "all" if str(category or "all").lower() == "all" else ann_category(category)
    if cat is None:
        print("[UFV] unknown annotation category '%s' (ptms variants sites mutagenesis domains ligands partners "
              "pdbekb all)." % category)
        return None
    rows = annotation_rows(obj, uid, cat)
    spec = str(ids if ids is not None else "all").strip()
    if search not in (None, ""):
        want = set(ann_search(rows, search))
    elif spec.lower() in ("", "all"):
        want = {r["id"] for r in rows}
    elif spec.lower() == "none":
        want = set()
    else:
        known = {r["id"] for r in rows}
        want = set()
        for t in re.split(r"[\s+;,]+", spec):
            if not t:
                continue
            if ":" not in t and cat != "all":
                t = "%s:%s" % (cat, t.upper() if cat == "lig" else t)
            if t in known:
                want.add(t)
            else:
                print("[UFV] no row '%s' in %s." % (t, cat))
    mode = (mode or "set").lower()
    cats = ANN_CATS if cat == "all" else (cat,)
    for c in cats:
        sub = {i for i in want if i.split(":", 1)[0] == c}
        if mode == "set":
            _ann_set(obj, uid, c, sub)
        else:
            _ann_set(obj, uid, c, sub, "remove" if mode.startswith("rem") else "add")
    res = _ann_apply(obj, uid, cats)
    print("[UFV] %s: %s." % (obj, ", ".join("%s %d" % (ANN_CAT_SHORT[c], res[c]) for c in cats if c in res)))
    return res


# ----------------------------------------------------------------------------------------------
# Projection commands  (each is independently toggleable; pair with ufv_hide to remove)
# ----------------------------------------------------------------------------------------------
@_undo_paused
def ufv_ptms(obj=None, uid=None):
    """ufv_ptms [object [, uniprot_id]] — PTM Ca spheres coloured by category."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    ann = fetch_annotations(uid)
    _ann_set(obj, uid, "ptm", [r["id"] for r in ann_rows(ann, "ptm")])     # = check every PTM row
    _ann_apply(obj, uid, ("ptm",))
    print("[UFV] %s: PTM layer on (%d categories)." % (obj, len({p["color"] for p in ann["ptms"]})))
VARIANT_COLOR_PROFILES = ("consequence", "disease", "provenance", "frequency", "substitution")


@_undo_paused
def ufv_variants(obj=None, uid=None, only="pathogenic", provenance=None, disease=None, color="consequence"):
    """ufv_variants [object [, uniprot_id [, only [, provenance [, disease [, color]]]]]]
    Variant Ca spheres, one per position (its most severe variant). `only` filters by consequence and
    may list several (space separated): pathogenic | benign | uncertain | deleterious | all (default
    pathogenic). `provenance`: any of reviewed / clinvar / cancer / population / other / large (default
    all). `disease`: disease labels separated by ';' (default all; 'none' hides all). `color` =
    consequence | disease | provenance | frequency | substitution (the extension's colour profiles;
    disease = first active disease's palette colour; frequency = gnomAD allele-frequency bin)."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    ann = fetch_annotations(uid)
    only = (only or "").lower()
    tokens = None if only in ("", "all") else set(re.split(r"[ ,]+", only.strip()))
    provs = _parse_provenances(provenance)
    dis = _parse_diseases(ann, disease)
    color = (color or "consequence").lower()
    if not any(color.startswith(c[:3]) for c in VARIANT_COLOR_PROFILES):
        color = "consequence"
    shown = {id(v) for v in _filter_variants(ann["variants"], tokens, provs, dis, False)}
    sty = _ann_style(obj)
    sty["color"], sty["dis"] = color, dis
    # = check exactly the matching variant rows (row id = index in the variant list)
    _ann_set(obj, uid, "var", ["var:%d" % i for i, v in enumerate(ann["variants"]) if id(v) in shown])
    _ann_apply(obj, uid, ("var",))
    groups = _variant_groups(ann, tokens, False, provs, dis, color)
    print("[UFV] %s: variant layer on (%s%s%s; colour by %s, %d positions)."
          % (obj, only or "all", ("; " + ", ".join(sorted(provs))) if provs is not None else "",
             ("; %d disease(s)" % len(dis)) if dis is not None else "", color,
             sum(len(v) for v in groups.values())))


def ufv_diseases(uid=None):
    """ufv_diseases [uniprot_id] — the disease summary (label, variant count, palette colour) used by
    the 'disease' variant colour profile and filter (data.js getDiseaseSummary)."""
    ann = fetch_annotations(_resolve_uid(uid))
    summ = _disease_summary(ann["variants"], ann.get("diseaseColors"))
    for name, d in summ.items():
        print("  %-48s %5d  %s" % (name[:48], d["count"], d["color"]))
    return summ


@_undo_paused
def ufv_sites(obj=None, uid=None):
    """ufv_sites [object [, uniprot_id]] — active / binding / metal sites as amber spheres."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    ann = fetch_annotations(uid)
    _ann_set(obj, uid, "site", [r["id"] for r in ann_rows(ann, "site")])
    _ann_apply(obj, uid, ("site",))
    print("[UFV] %s: site layer on (%d sites)." % (obj, len(ann["sites"])))


@_undo_paused
def ufv_mutagenesis(obj=None, uid=None):
    """ufv_mutagenesis [object [, uniprot_id]] — experimentally mutated residues (UniProt Mutagenesis)
    as brown spheres; hover the report for the observed effect."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    ann = fetch_annotations(uid)
    _ann_set(obj, uid, "mutag", [r["id"] for r in ann_rows(ann, "mutag")])
    _ann_apply(obj, uid, ("mutag",))
    print("[UFV] %s: mutagenesis layer on (%d features)." % (obj, len(ann.get("mutagenesis", []))))


@_undo_paused
@_notifies
def _show_partner_layer(obj, colors):
    """Render {(chain, resi): colour} as ONE sphere object (ufv_<obj>_prt) carved from the Cα scaffold
    and coloured with one alter. Empty -> layer hidden. Returns the number of spheres."""
    lobj = _ufv_layer_obj(obj, "prt")
    st = _SPHERE_STATE.setdefault(obj, {})
    built = _LAYER_BUILT.setdefault(obj, {})
    if not colors:
        st.pop("prt", None)
        if lobj in _object_set():
            cmd.disable(lobj)
        return 0
    by_chain = {}
    for (ch, r) in colors:
        by_chain.setdefault(ch, []).append(r)
    scaffold = _ca_scaffold(obj)
    sel = " or ".join("(chain %s and resi %s)" % (ch, _collapse_resi(sorted(rs))) for ch, rs in by_chain.items())
    with _batch():
        cmd.delete(lobj)
        cmd.create(lobj, "(%s) and (%s)" % (scaffold, sel), zoom=0)
        cmd.hide("everything", lobj)
        _tune_sphere_obj(lobj)
        idx = {}
        m = {}
        for k, c in colors.items():
            if c not in idx:
                idx[c] = cmd.get_color_index(_color_name(c))
            m[k] = idx[c]
        cmd.alter(lobj, "color = _m.get((chain, resv), color)", space={"_m": m})
        cmd.show("spheres", lobj)
        cmd.recolor(lobj)
        cmd.enable(lobj)
    built["prt"] = {"pos": frozenset(), "groups": None}
    st["prt"] = {"partner": sorted(colors)}
    return cmd.count_atoms(lobj)


@_undo_paused
def ufv_partner_layers(obj=None, layers="all"):
    """ufv_partner_layers [object [, layers]] — annotate the OTHER proteins of a complex (partner chains)
    with their own pathogenic variants (disease colours), PTMs and sites, one sphere per residue in the
    priority variants > ptms > sites (the extension's Partners layers). layers = all | off | any of
    variants ptms sites. Fetches each partner's annotations once (first 8 accessions)."""
    if obj and layers == "all" and obj not in cmd.get_names("objects") and all(w in ("all", "off", "none", "0", "variants", "ptms", "sites") for w in str(obj).replace(",", " ").lower().split()):
        obj, layers = None, obj
    obj = _resolve_object(obj)
    want = _parse_partner_layers(layers)
    if not want:
        _ann_forget(obj, ("prt",))
        _show_partner_layer(obj, None)
        print("[UFV] %s: partner layer off." % obj)
        return 0
    if not _partner_maps(obj):
        print("[UFV] %s has no partner proteins with a UniProt mapping." % obj)
        return 0
    cols = _partner_sphere_colors(_partner_annotations(obj), want)
    _ann_style(obj)["prt"] = set(want)
    _ann_set(obj, None, "prt", ["prt:%s:%s" % k for k in cols])             # = check those partner rows
    n = _show_partner_layer(obj, {k: c for k, (c, _t) in cols.items()})
    print("[UFV] %s: partner layer on (%s): %d residue(s) on %d partner chain(s)."
          % (obj, ", ".join(l for l in PARTNER_LAYERS if l in want), n, len({k[0] for k in cols})))
    return n


def _ligand_sel(obj):
    return "(%s) and not polymer and not solvent" % obj


def _has_ligands(obj):
    c = _LIG_PRESENT.get(obj)
    if c is None:
        try:
            c = cmd.count_atoms(_ligand_sel(obj)) > 0
        except Exception:
            c = False
        _LIG_PRESENT[obj] = c
    return c


_LIG_SHOWN = set()  # objects whose ligand layer is switched on


def _lig_objs(obj):
    """The ligand-layer objects that exist for `obj`: ufv_<obj>_lig (every organic ligand/cofactor
    copy, sticks) and ufv_<obj>_ion (ions/inorganics, small spheres)."""
    names = _object_set()
    return [o for o in (_ufv_layer_obj(obj, "lig"), _ufv_layer_obj(obj, "ion")) if o in names]


_ION_SEL = "resn " + "+".join(ION_CODES)
AFF_SETTINGS = {"min_identity": 0.25}          # AlphaFill donor-identity cutoff (0.25 = show every transplant)
AFF_IDENTITY_STEPS = (25, 30, 40, 50, 60, 70)  # the extension's preset buttons (%)


def _apply_aff_filter(obj):
    """Hide AlphaFill transplants whose donor identity is below AFF_SETTINGS['min_identity'] (modal.js
    applyLigandFilters): a show/hide on the two ligand-layer objects, never a rebuild. Returns
    (shown, total) transplant counts, or None when the structure carries no transplant metadata."""
    meta = _AFF_META.get(obj)
    if not meta:
        return None
    cut = float(AFF_SETTINGS["min_identity"])
    low = sorted(ch for ch, m in meta.items() if m.get("identity") is not None and m["identity"] < cut)
    total = sum(1 for m in meta.values() if m.get("identity") is not None)
    _show_lig_rows(obj, low)
    return total - len(low), total


def _show_lig_rows(obj, low=()):
    """Ligand-layer visibility: only the checked ligand codes (_LIG_ROWS; missing = all), minus AlphaFill
    transplants below the identity cut-off (chains `low`). show / hide on the two layer objects only."""
    codes = _LIG_ROWS.get(obj)
    names = _object_set()
    with _batch():
        for o, rep in ((_ufv_layer_obj(obj, "lig"), "sticks"), (_ufv_layer_obj(obj, "ion"), "spheres")):
            if o not in names:
                continue
            if codes is None:
                cmd.show(rep, o)
            else:
                cmd.hide(rep, o)
                if codes:
                    cmd.show(rep, "(%s) and resn %s" % (o, "+".join(sorted(codes))))
            if low:
                cmd.hide(rep, "(%s) and chain %s" % (o, "+".join(low)))


@_undo_paused
def ufv_ligand_identity(min_identity=25, obj=None):
    """ufv_ligand_identity pct [, object] — AlphaFill models: show only transplanted ligands whose donor
    has >= pct % sequence identity (presets 25/30/40/50/60/70; 25 = AlphaFill's floor = all)."""
    v = float(min_identity)
    AFF_SETTINGS["min_identity"] = v / 100.0 if v > 1 else v
    obj = _resolve_object(obj)
    r = _apply_aff_filter(obj)
    if r is None:
        print("[UFV] %s: no AlphaFill transplant metadata (load the AlphaFill model and run ufv_ligands)." % obj)
    else:
        print("[UFV] %s: %d/%d transplants with donor identity >= %d%%." % (obj, r[0], r[1],
                                                                          round(AFF_SETTINGS["min_identity"] * 100)))
    return r


@_undo_paused
@_notifies
def ufv_ligands(obj=None, uid=None, min_identity=None, _quiet=False):
    """ufv_ligands [object [, uid [, min_identity]]] — show bound ligands/cofactors (non-polymer,
    non-water). All organic ligand/cofactor copies go on ONE object (ufv_<obj>_lig, sticks) and ions
    (the extension's ION_CODES list) on another (ufv_<obj>_ion, spheres) — two objects instead of one per
    copy (a glycosylated receptor made 35). Single copies stay reachable with ufv_ligand_focus. Built
    once; toggling afterwards is an instant enable/disable. On an AlphaFill model the transplant metadata
    (donor identity, clash) is fetched once and transplants below `min_identity` % are hidden."""
    obj = _resolve_object(obj)
    if not _quiet:                  # the command = check every ligand row (the panel passes its own subset)
        _LIG_ROWS.pop(obj, None)
        _ann_set(obj, uid, "lig", ["lig:%s" % r for r in _ann_ligands(obj)])
    if min_identity not in (None, ""):
        v = float(min_identity)
        AFF_SETTINGS["min_identity"] = v / 100.0 if v > 1 else v
    if _is_alphafill(obj) and obj not in _AFF_META:
        meta = alphafill_transplants(_resolve_uid(uid, obj))     # network once (cached); GUI prefetches
        _AFF_META[obj] = meta or {}
    sel = _ligand_sel(obj)
    existing = _lig_objs(obj)
    if existing:
        _LIG_SHOWN.add(obj)
        with _batch():
            for o in existing:
                cmd.enable(o)
            if not (_FOCUS_SHOWN.get(obj) or {}).get("lig"):
                cmd.hide("sticks", sel)       # a past ligand focus may have left the main copy drawn
            if _apply_aff_filter(obj) is None:
                _show_lig_rows(obj)
        return sum(cmd.count_atoms(o) for o in existing)
    n = cmd.count_atoms(sel)
    if not n:
        print("[UFV] %s: no ligands/cofactors present." % obj)
        return 0
    lig, ion = _ufv_layer_obj(obj, "lig"), _ufv_layer_obj(obj, "ion")
    cps = ligand_copies(obj)
    copies = len(cps)
    ions_here = sorted({r for r, _c, _i in cps} & set(ION_CODES))   # short resn clause = fast selection
    org_sel = ("(%s) and not resn %s" % (sel, "+".join(ions_here))) if ions_here else sel
    ions_sel = ("(%s) and resn %s" % (sel, "+".join(ions_here))) if ions_here else "none"
    with _batch():
        if cmd.count_atoms(org_sel):
            cmd.create(lig, org_sel, zoom=0)
            cmd.hide("everything", lig)
            cmd.show("sticks", lig)
            try:
                cmd.color("green", "(%s) and elem C" % lig)
                cmd.color("atomic", "(%s) and (not elem C)" % lig)
            except Exception:
                pass
        if cmd.count_atoms(ions_sel):       # ions/metals: one grouped object
            cmd.create(ion, ions_sel, zoom=0)
            cmd.hide("everything", ion)
            cmd.show("spheres", ion)
            cmd.set("sphere_scale", 0.3, ion)
            try:
                cmd.color("atomic", ion)
            except Exception:
                pass
        cmd.hide("sticks", sel)               # ensure the main object's copies aren't drawn twice
        flt = _apply_aff_filter(obj)
        if flt is None and _LIG_ROWS.get(obj) is not None:
            _show_lig_rows(obj)
    _LIG_SHOWN.add(obj)
    print("[UFV] %s: ligands shown (%d copies, %d atoms)%s." % (
        obj, copies, n, (", AlphaFill transplants %d/%d at identity >= %d%%" %
                         (flt[0], flt[1], round(AFF_SETTINGS["min_identity"] * 100))) if flt else ""))
    return n


@_undo_paused
@_notifies
def ufv_ligands_hide(obj=None, _quiet=False):
    """ufv_ligands_hide [object] — hide the ligand layer (instant: disables ufv_<obj>_lig / _ion)."""
    obj = _resolve_object(obj)
    _LIG_SHOWN.discard(obj)
    if not _quiet:
        _ann_forget(obj, ("lig",))
    for o in _lig_objs(obj):
        try:
            cmd.disable(o)
        except Exception:
            pass


def ligand_copies(obj):
    """Distinct ligand/cofactor copies in the structure -> sorted [(resn, chain, resi)] (no water)."""
    rows = []
    try:
        cmd.iterate(_ligand_sel(obj), "rows.append((resn, chain, resi))", space={"rows": rows})
    except Exception:
        pass
    return sorted(set(rows))


def enumerate_ligands(obj):
    """Ligand/cofactor components in the structure -> {resn: instance_count} (distinct chain+resi
    copies, no water) so multiple copies of the same ligand are counted."""
    counts = {}
    for resn, _ch, _ri in ligand_copies(obj):
        counts[resn] = counts.get(resn, 0) + 1
    return counts


def ligand_chemistry(obj=None, codes=None):
    """For each distinct ligand, fetch CCD chemistry + PubChem fingerprint and rank the others by
    Tanimoto similarity. Returns {resn: {info..., 'similar': [(resn, score)]}}. Network-heavy — call
    off the UI thread, passing `codes` (enumerate the ligands with enumerate_ligands on the main
    thread first: PyMOL's API must not be driven from a worker). `obj` alone still works from the
    command line / main thread."""
    if codes is None:
        codes = ligand_codes_ordered(obj)
    codes = list(dict.fromkeys(codes))      # structure order (modal.js renderSimilarLigands), not sorted

    def one(c):                              # CCD record -> PubChem fingerprint + CID, per code
        inf = ligand_info(c)
        ik = inf.get("inchikey")
        return inf, ligand_fingerprint(ik), (pubchem_cid(ik) if ik else None)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=8) as ex:     # independent per code: all codes in flight at once
        got = dict(zip(codes, ex.map(one, codes)))
    out = {}
    for c in codes:
        inf, fp, cid = got[c]
        sims = []
        if fp:
            # the extension scores at most 20 OTHER ligand types (structure order) and keeps score > 0
            for d in [x for x in codes if x != c][:LIGAND_SIM_CAP]:
                t = tanimoto(fp, got[d][1])
                if t is not None and t > 0:
                    sims.append((d, round(t, 3)))
            sims.sort(key=lambda x: -x[1])
        out[c] = dict(inf, similar=sims, cid=cid, links=ligand_links(inf, cid))
    return out


def ligand_codes_ordered(obj):
    """Distinct ligand CCD codes in model (first-appearance) order, water excluded."""
    rows = []
    try:
        cmd.iterate(_ligand_sel(obj), "rows.append(resn)", space={"rows": rows})
    except Exception:
        pass
    return list(dict.fromkeys(rows))


def ligand_instances(obj, resn):
    """Distinct copies of a ligand by name -> sorted [(chain, resi), ...] so the GUI can step through
    each occurrence of the same component (e.g. several HEM/ATP in one structure)."""
    rows = []
    try:
        cmd.iterate("(%s) and resn %s and not solvent" % (obj, resn),
                    "rows.append((chain, resi))", space={"rows": rows})
    except Exception:
        pass
    seen, out = set(), []
    for ch, resi in rows:
        key = (ch, resi)
        if key not in seen:
            seen.add(key); out.append(key)
    out.sort(key=lambda x: (x[0], _resi_int(x[1]) or 0))
    return out


def _resi_tok(resi):
    """A resi value as a PyMOL selection token (negative numbers need the minus escaped)."""
    resi = str(resi)
    return "\\" + resi if resi.startswith("-") else resi


def _color_sticks_by_annotation(obj, nb_sel, cmap):
    """Set stick_color of the residues in `nb_sel` from their annotation colour (variant/PTM/site)."""
    atoms = []
    cmd.iterate("(%s) and name CA" % nb_sel, "atoms.append((chain, resi))", space={"atoms": atoms})
    by_color = {}
    for ch, resi in set(atoms):
        c = cmap.get(_resi_to_uni(obj, ch, resi))
        if c:
            by_color.setdefault(c, {}).setdefault(ch, []).append(_resi_tok(resi))
    for c, chres in by_color.items():
        clauses = " or ".join("(chain %s and resi %s)" % (ch, "+".join(rs)) if ch else "(resi %s)" % "+".join(rs)
                              for ch, rs in chres.items())
        cmd.set("stick_color", _color_name(c), "(%s) and (%s)" % (obj, clauses))


@_undo_paused
def ufv_ligand_focus(obj=None, resn=None, chain=None, resi=None, uid=None, zoom=1):
    """ufv_ligand_focus [object,] resn [, chain, resi] — zoom to a ligand and show it + its 5 Å pocket
    as sticks. With chain+resi, focus that single copy; otherwise all copies of the named component.
    zoom=0 rebuilds the focus without moving the camera (e.g. after the Nearby radius changed)."""
    obj, resn = _shift_obj(obj, resn)
    obj = _resolve_object(obj)
    if chain is not None and resi is not None:
        lig = ("(%s) and resn %s and chain %s and resi %s" % (obj, resn, chain, _resi_tok(resi)) if chain
               else "(%s) and resn %s and resi %s" % (obj, resn, _resi_tok(resi)))
    else:
        lig = "(%s) and resn %s" % (obj, resn)
    if not cmd.count_atoms(lig):
        return
    cmap = _annotation_color_map(_CACHE.get(_resolve_uid(uid, obj)) or {})
    # Named selections: PyMOL stores them as fixed atom sets, so the `around 5` search runs once here
    # instead of on every later show / hide / unset of a string selection (e.g. in _clear_focus).
    lsel = _ufv_layer_obj(obj, "focuslig")
    pocket = _ufv_layer_obj(obj, "focusnb")
    fset = _ufv_layer_obj(obj, "focusset")
    disabled = []
    with _batch():
        _clear_focus(obj)                       # drop any previous residue/ligand focus sticks
        enabled = _enabled_set()
        for lig_layer in _lig_objs(obj):        # hide the ligand-layer objects so this copy isn't drawn twice
            if lig_layer in enabled:
                cmd.disable(lig_layer); disabled.append(lig_layer)
        cmd.select(lsel, lig, enable=0)
        # scoped to the ligand's own object: overlapping models (aligned AF / PDB) must not add residues
        cmd.select(pocket, "(%s) and byres ((%s) around %g) and polymer%s" % (obj, lsel, NEARBY_SETTINGS["dist"],
                                                                              _not_transplant(obj)), enable=0,
                   state=_obj_cur_state(obj))
        cmd.select(fset, "%s or %s" % (lsel, pocket), enable=0)
        cmd.show("sticks", lsel)
        cmd.color("green", "%s and elem C" % lsel)
        cmd.set("stick_radius", 0.18, lsel)     # fatten the focused copy so it reads as 'the selected one'
        cmd.show("sticks", pocket)
        cmd.set("stick_color", "grey70", pocket)
        # Colour the pocket residues by their disease annotation (variant/PTM/site), like residue focus,
        # instead of leaving them all grey.
        _color_sticks_by_annotation(obj, pocket, cmap)
        near, partner_near = _near_from_selection(obj, pocket)
    _FOCUS_SHOWN[obj] = {"sticks": fset, "sel": None, "disabled": disabled, "lig": lsel,
                         "selnames": [fset, pocket, lsel], "near": near, "partnerNear": partner_near,
                         "ligand": (resn, chain, resi), "uid": uid}
    if int(zoom):
        _focus_zoom(lsel, 6)
    return near


@_undo_paused
def ufv_show_pocket(obj=None, positions=None, pocket=None, uid=None):
    """ufv_show_pocket [object,] positions — highlight a predicted (druggable) binding pocket: show
    its lining residues as orange sticks + a translucent surface and zoom in. `pocket=<pocketId>` takes
    the residues of a ProtVar pocket already fetched (residue report or ufv_pockets_all)."""
    obj, positions = _shift_obj(obj, positions)
    obj = _resolve_object(obj)
    if pocket not in (None, ""):
        pk = _find_pocket(_resolve_uid(uid, obj), pocket)
        if not pk:
            print("[UFV] ProtVar pocket %s not fetched yet (ufv_pockets_all or a residue report)." % pocket)
            return
        positions = pk["resid"]
    if isinstance(positions, str):          # from the command line: "12+15+40" / "12,15,40"
        positions = [p for p in re.split(r"[\s,+]+", positions) if p]
    if not positions:
        return
    sel = _sel_for_positions(obj, [int(p) for p in positions])
    if not sel:
        return
    res = "byres (%s)" % sel
    with _batch():
        _clear_focus(obj)
        cmd.show("sticks", res)
        cmd.set("stick_color", _color_name("#fb8c00"), res)
        try:
            cmd.show("surface", res)
            cmd.set("transparency", 0.45, res)
            cmd.set("surface_color", _color_name("#ffd180"), res)
        except Exception:
            pass
    _FOCUS_SHOWN[obj] = {"sticks": res, "sel": None, "disabled": [], "surface": res}
    _focus_zoom(res, 2)


@_undo_paused
def ufv_domains(obj=None, uid=None):
    """ufv_domains [object [, uniprot_id]] — colour the cartoon by domain / region / repeat."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    _set_cartoon_layer(obj, "domains", grey=False)
    dom_spheres = {}
    pos_color = {}
    n = 0
    for d in ann["domains"]:
        if d["isRange"]:
            for p in range(d["position"], d["endPosition"] + 1):
                pos_color[p] = d["color"]
        else:
            dom_spheres.setdefault(d["color"], []).append(d["position"])
        n += 1
    _recolor_positions(obj, pos_color, base_name="gray80")
    _set_sphere_layer(obj, "dom", dom_spheres or None)
    print("[UFV] %s: coloured %d domain/region features." % (obj, n))


@_undo_paused
def ufv_topology(obj=None, uid=None):
    """ufv_topology [object [, uniprot_id]] — colour the cartoon by membrane topology."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    _set_cartoon_layer(obj, "topology", grey=False)
    pos_color = {}
    for t in ann["topology"]:
        for p in range(t["start"], t["end"] + 1):
            pos_color[p] = t["color"]
    _recolor_positions(obj, pos_color, base_name="gray80")
    print("[UFV] %s: coloured %d topology segments." % (obj, len(ann["topology"])))


@_undo_paused
def ufv_alphamissense(obj=None, uid=None):
    """ufv_alphamissense [object [, uniprot_id]] — colour by mean AlphaMissense pathogenicity."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    am = ann["amMean"]
    if not am:
        print("[UFV] no AlphaMissense scores available for %s." % ann["uid"])
        return
    _set_cartoon_layer(obj, "alphamissense", grey=False)
    pos_color = {}
    for pos, avg in am.items():
        pos_color[pos] = ("#b71c1c" if avg >= 0.78 else "#e06666" if avg >= 0.564
                          else "#b9c2cf" if avg >= 0.34 else "#3d85c8")
    _recolor_positions(obj, pos_color, base_hex="#b9c2cf")
    print("[UFV] %s: coloured by AlphaMissense (%d scored positions)." % (obj, len(am)))


@_undo_paused
def ufv_burden(obj=None, uid=None):
    """ufv_burden [object [, uniprot_id]] — colour mutation/phenotype burden-positive residues."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    _set_cartoon_layer(obj, "burden", grey=False)
    _recolor_positions(obj, {p: "#e65100" for p in ann["burden"]}, base_name="gray80")
    print("[UFV] %s: %d burden-positive residues." % (obj, len(ann["burden"])))


@_undo_paused
def ufv_plddt(obj=None, uid=None):
    """ufv_plddt [object] — colour by AlphaFold pLDDT (the B-factor column), exact AFDB bins.
    Only meaningful for predicted models (the B column is pLDDT there)."""
    obj = _resolve_object(obj)
    if not _structure_is_predicted(obj):
        print("[UFV] pLDDT applies to predicted models only (%s is experimental)." % obj)
        return
    _set_cartoon_layer(obj, "plddt")
    pol = "(%s) and polymer" % obj      # leave bound ligands their own colours
    with _batch():
        # AlphaFold DB bins: >=90 dark blue, 70-90 light blue, 50-70 yellow, <50 orange.
        cmd.color(_color_name("#ff7d45"), "(%s) and b<50" % pol)
        cmd.color(_color_name("#ffdb13"), "(%s) and (b>50 or b=50) and b<70" % pol)
        cmd.color(_color_name("#65cbf3"), "(%s) and (b>70 or b=70) and b<90" % pol)
        cmd.color(_color_name("#0053d6"), "(%s) and (b>90 or b=90)" % pol)
    print("[UFV] %s: coloured by pLDDT." % obj)


@_undo_paused
def ufv_bfactor(obj=None, uid=None):
    """ufv_bfactor [object] — colour by crystallographic/EM B-factor with the extension's exact
    blue-white-red gradient (#313695 -> #f7f7f7 -> #d73027), clamped 0-100. Experimental only
    (for predicted models the B column is pLDDT — use that mode instead)."""
    obj = _resolve_object(obj)
    if _structure_is_predicted(obj):
        print("[UFV] B-factor applies to experimental structures only (%s is a predicted model)." % obj)
        return
    _set_cartoon_layer(obj, "bfactor")
    try:
        cmd.spectrum("b", "0x313695 0xf7f7f7 0xd73027", "(%s) and polymer" % obj, minimum=0, maximum=100)
    except Exception:
        cmd.spectrum("b", "blue_white_red", "(%s) and polymer" % obj)
    print("[UFV] %s: coloured by B-factor." % obj)


# ----------------------------------------------------------------------------------------------
# Geometry + structure-dependent analyses (hotspots, contact hubs) + residue report / focus
# ----------------------------------------------------------------------------------------------
_FOCUS_SHOWN = {}     # obj -> stick selection currently shown by a focus


def _ca_geometry(obj, refresh=False, state=None):
    """Modelled Cα atoms of `obj` as [{uniPos, chain, resi, x, y, z, b}] in model order (cached per
    object), one per (chain, resv) — first wins, like the extension's residueGeometry (alt-loc copies are
    dropped). `resi` is the integer author residue number (resv): the pocket analysis does resi±1
    arithmetic on it. Insertion-code residues are kept for geometry but have no UniProt position.
    Multi-state objects: the CURRENT state (or `state`); caches are kept per state (_traj_geom_state)."""
    st = _traj_geom_state(obj, state)
    if not refresh and obj in _GEOM_CACHE:
        return _GEOM_CACHE[obj]
    rows = []
    try:
        cmd.iterate_state(st, "(%s) and name CA and polymer" % obj,
                          "rows.append((chain, resi, resv, x, y, z, b))", space={"rows": rows})
    except Exception:
        rows = []
    geom, seen = [], set()
    for (ch, resi, resv, x, y, z, b) in rows:
        if (ch, resv) in seen:
            continue
        seen.add((ch, resv))
        geom.append({"uniPos": _resi_to_uni(obj, ch, resi), "chain": ch, "resi": resv,
                     "x": x, "y": y, "z": z, "b": b})
    _GEOM_CACHE[obj] = geom
    return geom


HOTSPOT_PARTNERS = True       # fold partner-protein classified residues into the hotspot test


def _analysis(obj, uid, kind, geom=None, ann=None, **kw):
    """Cached structure analysis for (obj, uid): kind = 'hotspots' | 'hubs' (merged {uniPos: tier}),
    'pockets' (full computePockets result: byPos with q for every candidate), 'ptmprox' (PTM-variant
    proximity at the current radii). Per-chain hotspot / hub maps are cached as '<kind>_bychain'.
    Re-picking a colouring mode or exporting the CSV never recomputes. Safe to call from a worker
    thread when `geom` (read on the main thread) and `ann` are passed in (network allowed there)."""
    if kind == "ptmprox":
        d8 = float(kw.get("d8", PROX_SETTINGS["ptm"]))
        d12 = float(kw.get("d12", PROX_SETTINGS["var"]))
        key = (obj, uid, "ptmprox", d8, d12)
    else:
        key = (obj, None if kind.startswith("hubs") else uid, kind)
    if geom is None and _HAS_PYMOL:
        _traj_geom_state(obj)               # trajectories: the cache must belong to the state shown
    hit = _ANALYSIS_CACHE.get(key)
    if hit is not None:
        return hit
    if geom is None:
        geom = _ca_geometry(obj)
    if kind in ("hubs", "hubs_bychain"):
        bych = {}
        res = _per_chain_merge(geom, _betweenness_hubs, _HUB_RANK, bych)
        _ANALYSIS_CACHE[(obj, None, "hubs_bychain")] = bych
        _ANALYSIS_CACHE[(obj, None, "hubs")] = res
        return _ANALYSIS_CACHE[key]
    ann = ann or _CACHE.get(uid) or fetch_annotations(uid)
    if kind in ("hotspots", "hotspots_bychain"):
        variants = ann.get("variants", [])
        classes = _hotspot_classes(variants)
        partners = []
        if HOTSPOT_PARTNERS and len(classes[0]) >= HOTSPOT_MIN_PATH and _partner_maps(obj):
            try:
                partners = _resolve_partner_points(geom, _partner_points(obj))
            except Exception:
                partners = []
        bych = {}
        res = _per_chain_merge(geom, lambda gs: _compute_hotspots(gs, variants, partners=partners, classes=classes),
                               _TIER_RANK, bych)
        _ANALYSIS_CACHE[(obj, uid, "hotspots_bychain")] = bych
        _ANALYSIS_CACHE[(obj, uid, "hotspots")] = res
        return _ANALYSIS_CACHE[key]
    if kind == "pockets":
        pae = None
        info = _af_model_info(obj)
        if info and ann.get("amMap"):
            try:
                pae = _pae_matrix(*info)
            except Exception:
                pae = None
        res = _compute_pockets(geom, ann.get("amMap") or {}, ann.get("sequence") or "", pae)
    elif kind == "ptmprox":
        res = _ptm_variant_proximity(ann.get("ptms", []), ann.get("variants", []), geom, d8, d12)
    else:
        raise ValueError(kind)
    _ANALYSIS_CACHE[key] = res
    return res


def _analysis_cached(obj, uid, kind):
    if obj and _HAS_PYMOL:
        _traj_geom_state(obj)               # trajectories: per-state caches (main thread)
    if kind == "ptmprox":
        return _ANALYSIS_CACHE.get((obj, uid, "ptmprox", float(PROX_SETTINGS["ptm"]), float(PROX_SETTINGS["var"])))
    return _ANALYSIS_CACHE.get((obj, None if kind.startswith("hubs") else uid, kind))


@_undo_paused
def _recolor_positions(obj, pos_color, base_hex=None, base_name=None):
    """Colour `obj`'s polymer atoms from a {uniProtPos: hex} map in ONE `alter` pass. This replaces
    the old 'one cmd.color per colour group with a big resi-list selection' pattern, whose repeated
    selection parsing over the whole structure was the slow part of every cartoon-colouring mode.
    The map is pre-resolved to {(chain, resv): colour index}, so `alter` only does a dict lookup per
    atom (no Python function call per atom). Atoms not in the map keep their colour, or take
    `base_hex` / the named PyMOL colour `base_name` if given. Returns how many positions were coloured."""
    idx_of, idx_map = {}, {}
    for pos, hexc in pos_color.items():
        if hexc not in idx_of:                 # one colour lookup per colour, not per residue
            idx_of[hexc] = _color_index(_color_name(hexc))
        idx_map[int(pos)] = idx_of[hexc]
    if not idx_map and base_hex is None and base_name is None:
        return 0
    m = _res_key_map(obj, idx_map)
    with _batch():
        if base_hex or base_name:
            cmd.alter("(%s) and polymer" % obj, "color = _m.get((chain, resv), _b)",
                      space={"_m": m, "_b": _color_index(_color_name(base_hex) if base_hex else base_name)})
        else:
            cmd.alter("(%s) and polymer" % obj, "color = _m.get((chain, resv), color)", space={"_m": m})
        cmd.recolor(obj)
    return len(idx_map)


@_undo_paused
def _color_tiers(obj, tiers, colors, base_name=None):
    """Colour tiers {uniPos: tier} with `colors` {tier: hex}. base_name (e.g. "gray80"): every other
    polymer atom takes that colour in the same pass (use with _set_cartoon_layer(..., grey=False))."""
    pos_color = {}
    for uni, tier in tiers.items():
        c = colors.get(tier)
        if c:
            pos_color[uni] = c
    _recolor_positions(obj, pos_color, base_name=base_name)
    return len(tiers)


HOTSPOT_COLORS = {"strong": "#b71c1c", "moderate": "#e64a19", "weak": "#ffa726"}
HUB_COLORS = {"strong": "#6a1b9a", "moderate": "#ab47bc"}
POCKET_COLORS = {"pocket": "#00897b", "exposed": "#8e24aa"}


@_undo_paused
def ufv_hotspots(obj=None, uid=None):
    """ufv_hotspots [object [, uniprot_id]] — colour pathogenic-enrichment hotspot tiers (extension-exact
    mulberry32 permutation null; partner-protein disease residues of a complex are folded in)."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    tiers = _analysis(obj, ann["uid"], "hotspots", ann=ann)
    _set_cartoon_layer(obj, "hotspots", grey=False)
    n = _color_tiers(obj, tiers, HOTSPOT_COLORS, "gray80")
    print("[UFV] %s: %d hotspot residues." % (obj, n))


@_undo_paused
def ufv_contacthubs(obj=None, uid=None):
    """ufv_contacthubs [object] — colour long-range contact-hub tiers (Cα betweenness, all sources)."""
    obj = _resolve_object(obj)
    tiers = _analysis(obj, None, "hubs")
    _set_cartoon_layer(obj, "contacthubs", grey=False)
    n = _color_tiers(obj, tiers, HUB_COLORS, "gray80")
    print("[UFV] %s: %d contact-hub residues." % (obj, n))


@_undo_paused
def ufv_pockets(obj=None, uid=None, q=None):
    """ufv_pockets [object [, uniprot_id [, q]]] — constraint-pocket prioritisation: residues more
    evolutionarily constrained (AlphaMissense) than their structural burial explains, clustered in
    3D (AlphaFold models: spatial weights gated by the PAE). Buried/concave candidates are teal
    ('pocket'), exposed ones purple. `q` = FDR sensitivity (0.01-0.40, default 0.10); every candidate's
    q is kept, so 'ufv_pocket_q 0.2' recolours without recomputing."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    if not ann.get("amMap"):
        print("[UFV] no AlphaMissense data — constraint-pocket needs it.")
        return
    if q not in (None, ""):
        _set_pocket_q(q)
    res = _analysis(obj, ann["uid"], "pockets", ann=ann)
    _set_cartoon_layer(obj, "pockets", grey=False)
    n = _color_tiers(obj, _pocket_classes(res), POCKET_COLORS, "gray80")
    print("[UFV] %s: %d constraint-pocket residues at q <= %.2f (%d candidates%s)%s."
          % (obj, n, POCKET_SETTINGS["q"], len(res.get("byPos") or {}), ", PAE-gated" if res.get("hasPae") else "",
             (" — " + res["reason"]) if res.get("reason") else ""))


@_undo_paused
def ufv_pocket_q(q=None, obj=None, uid=None):
    """ufv_pocket_q [q [, object]] — set the constraint-pocket sensitivity (FDR q threshold, 0.01-0.40)
    and recolour objects showing the constraint-pocket mode from the cached result (no recompute).
    Without q, prints the current value."""
    if q in (None, ""):
        print("[UFV] constraint-pocket sensitivity q <= %.2f" % POCKET_SETTINGS["q"])
        return POCKET_SETTINGS["q"]
    q = _set_pocket_q(q)
    objs = [_resolve_object(obj)] if obj else [o for o, m in _CARTOON_LAYER.items() if m in ("pockets", "constraint pocket")]
    for o in objs:
        u = _resolve_uid(uid, o)
        res = _analysis_cached(o, u, "pockets")
        if res is None:
            continue
        _set_cartoon_layer(o, _CARTOON_LAYER.get(o) or "pockets", grey=False)
        n = _color_tiers(o, _pocket_classes(res, q), POCKET_COLORS, "gray80")
        print("[UFV] %s: %d constraint-pocket residues at q <= %.2f." % (o, n, q))
    return q


def ufv_prox(ptm_radius=None, var_radius=None):
    """ufv_prox [ptm_radius [, var_radius]] — PTM–variant proximity radii (Å, 2-30; defaults 8 / 12):
    PTMs within ptm_radius and pathogenic variants within var_radius of a residue are listed in its
    report, and the PTM-centric tiers (CSV) use them as the 8 / 12 Å cut-offs."""
    for k, v in (("ptm", ptm_radius), ("var", var_radius)):
        if v not in (None, ""):
            try:
                PROX_SETTINGS[k] = float(min(30.0, max(2.0, float(v))))
            except ValueError:
                pass
    print("[UFV] proximity radii: PTMs <= %.1f Å, pathogenic variants <= %.1f Å"
          % (PROX_SETTINGS["ptm"], PROX_SETTINGS["var"]))
    return dict(PROX_SETTINGS)


# ----------------------------------------------------------------------------------------------
# Services (WP5-WP7): ProtVar pockets crawl, PDBe-KB binding context, Open Targets tractability,
# AlphaFill transplants, RCSB strucmotif search, pocket evidence / ligand clash, Nearby, flags,
# B-factor export. Network helpers are pure (safe in worker threads); cmd.* only in the commands.
# ----------------------------------------------------------------------------------------------


@_undo_paused
def ufv_pockets_all(obj=None, uid=None):
    """ufv_pockets_all [object [, uid]] — every ProtVar (AutoSite) pocket of the protein: crawls all
    positions once (slow the first time: one request per residue; cached), lists the pockets, and colours
    each pocket's lining residues with the extension's 8-colour pocket palette (best-scoring pocket wins
    where pockets overlap) in ONE alter pass on the structure — no per-pocket objects. Show one pocket's
    surface with: ufv_show_pocket obj, pocket=<id>."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    ann = fetch_annotations(uid)
    n = len(ann.get("sequence") or "")
    res = fetch_all_pockets(uid, n, lambda d, t: print("[UFV] pockets: %d/%d residues scanned" % (d, t))
                            if d % 250 == 0 else None)
    if not res:
        print("[UFV] %s: no ProtVar pockets detected across the protein." % uid)
        return None
    _show_all_pockets(obj, res["pockets"])
    for i, p in enumerate(res["pockets"]):
        print("  pocket %-6s %s  score %s  res %d  pLDDT %s  Rg %s" % (
            p["pocketId"], POCKET_PALETTE[i % len(POCKET_PALETTE)],
            "%.0f" % p["score"] if _finite(p["score"]) else "-", len(p["resid"]),
            "%.0f" % p["meanPlddt"] if _finite(p["meanPlddt"]) else "-",
            "%.1f" % p["radGyration"] if _finite(p["radGyration"]) else "-"))
    print("[UFV] %s: %d ProtVar pocket(s) coloured on %s." % (uid, res["count"], obj))
    return res


@_undo_paused
def _show_all_pockets(obj, pockets):
    """Colour every pocket's lining residues (palette by pocket rank) as a cartoon layer: one alter."""
    pos_color = {}
    for i, p in enumerate(pockets):
        c = POCKET_PALETTE[i % len(POCKET_PALETTE)]
        for r in p.get("resid") or []:
            pos_color.setdefault(int(r), c)
    _set_cartoon_layer(obj, "allpockets", grey=False)
    return _recolor_positions(obj, pos_color, base_name="gray80")
KB_COLORS = {"ligand": "#2e7d32", "interface": "#1565c0", "both": "#00838f"}


@_undo_paused
def ufv_binding(pos, uid=None):
    """ufv_binding uniprot_pos [, uid] — PDBe-KB experimental binding at a residue: ligands bound there
    (top 8 by number of PDB entries) and protein-protein interfaces (top 4), across every PDB structure
    of the protein."""
    uid = _resolve_uid(uid)
    kb = fetch_pdbekb(uid)
    lt, nl, it, ni = _kb_residue(kb, pos)
    print("[UFV] %s %s — PDBe-KB experimental binding (https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/%s)"
          % (uid, pos, _canon_acc(uid)))
    if not (lt or it):
        print("  no ligand / interface evidence at this residue.")
    if lt:
        print("  Ligands: " + ", ".join("%s (%d PDB%s)" % (l["code"], l["pdbCount"], (", e.g. %s" % str(l["pdb"]).upper())
                                                        if l.get("pdb") else "") for l in lt)
              + ((" +%d" % (nl - 8)) if nl > 8 else ""))
    if it:
        print("  PPI interface: " + "; ".join("%s (%d)" % (p.get("name") or p.get("accession"), p["pdbCount"]) for p in it)
              + ((" +%d" % (ni - 4)) if ni > 4 else ""))
    return {"ligands": lt, "interfaces": it}


@_undo_paused
def ufv_pdbekb(obj=None, kind="both", uid=None):
    """ufv_pdbekb [object [, ligand|interface|both [, uid]]] — PDBe-KB layer: Cα spheres on residues seen
    binding a ligand (green) and/or in a protein interface (blue; both = teal) in ANY experimental PDB
    structure of the protein. One sphere object (ufv_<obj>_kb); ufv_hide obj, binding removes it."""
    if obj and str(obj).lower() in ("both", "ligand", "ligands", "interface", "interfaces", "ppi") and obj not in cmd.get_names("objects"):
        obj, kind = None, obj
    obj = _resolve_object(obj)
    kind = (kind or "both").strip().lower()
    kb = fetch_pdbekb(_resolve_uid(uid, obj))
    if not kb:
        print("[UFV] no PDBe-KB ligand / interface data.")
        _hide_layer(obj, "kb")
        return 0
    _ann_style(obj)["kb"] = "both" if kind not in ("ligand", "ligands", "interface", "interfaces", "ppi") else         ("ligand" if kind.startswith("lig") else "interface")
    lig = set(kb["ligandByPos"]) if kind in ("both", "ligand", "ligands") else set()
    ifc = set(kb["interfaceByPos"]) if kind in ("both", "interface", "interfaces", "ppi") else set()
    groups = {}
    for p in lig | ifc:
        c = KB_COLORS["both"] if (p in lig and p in ifc) else KB_COLORS["ligand"] if p in lig else KB_COLORS["interface"]
        groups.setdefault(c, []).append(p)
    _ann_set(obj, None, "kb", ["kb:%d" % p for p in lig | ifc])               # = check those PDBe-KB rows
    n = _set_sphere_layer(obj, "kb", {c: sorted(v) for c, v in groups.items()} or None)
    print("[UFV] %s: PDBe-KB %s layer — %d residue(s) (%d ligand-site, %d interface)." % (obj, kind, n, len(lig), len(ifc)))
    return n


def ufv_tractability(uid=None):
    """ufv_tractability [uid] — Open Targets target tractability (positive buckets per modality)."""
    uid = _resolve_uid(uid)
    ot = fetch_open_targets(uid)
    if not ot or not ot["groups"]:
        print("[UFV] %s: no Open Targets tractability evidence." % uid)
        return ot
    print("[UFV] %s (%s) tractability — %s" % (ot["symbol"], ot["ensemblId"], ot["url"]))
    for g in ot["groups"]:
        print("  %-15s %s" % (g["modality"], ", ".join(g["labels"])))
    return ot


def _transplant_chains(obj):
    """viewer-molstar.js _focusNeighbourhood: chains carrying a non-water HETATM that are NOT the main
    (largest) chain = AlphaFill transplant chains (incl. donor binding-site PEPTIDES). AlphaFold-based
    models ONLY: in an experimental complex every glycosylated subunit carries HETATMs, and excluding
    them emptied Nearby. Cached per object (one iterate)."""
    if obj in _TRANSPLANT:
        return _TRANSPLANT[obj]
    out = frozenset()
    if _is_af_based(obj):
        try:
            if cmd.count_atoms("(%s) and hetatm and not solvent" % obj):
                res = {}
                het = set()
                rows = []
                cmd.iterate("(%s) and not solvent" % obj, "rows.append((chain, resi, type))", space={"rows": rows})
                for ch, resi, typ in rows:
                    res.setdefault(ch, set()).add(resi)
                    if typ == "HETATM":
                        het.add(ch)
                main = max(res, key=lambda c: len(res[c])) if res else None
                out = frozenset(c for c in het if c != main)
        except Exception:
            out = frozenset()
    _TRANSPLANT[obj] = out
    return out


def _not_transplant(obj):
    """' and not chain X+Y' for the AlphaFill transplant chains of `obj` ('' when there are none)."""
    tc = _transplant_chains(obj)
    return (" and not chain %s" % "+".join(sorted(tc))) if tc else ""


def motif_search_url(obj, positions, uid=None, geom=None):
    """modal.js openPocketSearch: an RCSB strucmotif query URL for a pocket — the (up to) 8 residues
    closest to the pocket's Cα centroid. PDB entries: mmCIF label ids from residue_listing; AlphaFold /
    AlphaFill models: entry AF_AF<acc>F1, chain A, label_seq_id = UniProt position. None if < 2 residues
    resolve or the structure is neither."""
    st = _obj_struct(obj)
    pdb_id = st.get("pdbId") or (UFV_MAPS.get(obj) or {}).get("pdbId")
    is_af = not pdb_id and _is_af_based(obj)
    if not pdb_id and not is_af:
        return None
    geo = {}
    for g in (geom if geom is not None else _ca_geometry(obj)):
        if g["uniPos"] is not None:
            geo[g["uniPos"]] = g                  # last chain wins (Map.set), as the extension
    pts = [geo[int(p)] for p in dict.fromkeys(int(p) for p in positions or []) if int(p) in geo]
    if len(pts) < 2:
        return None
    cx = sum(g["x"] for g in pts) / len(pts)
    cy = sum(g["y"] for g in pts) / len(pts)
    cz = sum(g["z"] for g in pts) / len(pts)
    picked = sorted(pts, key=lambda g: (g["x"] - cx) ** 2 + (g["y"] - cy) ** 2 + (g["z"] - cz) ** 2)[:8]
    if pdb_id:
        lab = _pocket_label_ids(pdb_id)
        rid = [{"label_asym_id": lab[k][0], "label_seq_id": lab[k][1]}
               for k in ("%s|%s" % (g["chain"], g["resi"]) for g in picked) if k in lab]
        entry = pdb_id.upper()
    else:
        acc = _canon_acc(_resolve_uid(uid, obj) or "")
        if not acc:
            return None
        entry = "AF_AF%sF1" % acc
        rid = [{"label_asym_id": "A", "label_seq_id": g["uniPos"]} for g in picked]
    if len(rid) < 2:
        return None
    q = {"query": {"type": "terminal", "service": "strucmotif",
                   "parameters": {"value": {"entry_id": entry, "residue_ids": rid}}}, "return_type": "assembly"}
    return "https://www.rcsb.org/search?request=" + _uri(json.dumps(q, separators=(",", ":")))


def ufv_motif_search(obj=None, positions=None, uid=None, pocket=None, browser=1):
    """ufv_motif_search [object,] positions | pocket=<id> — RCSB Structure Motif Search for a pocket
    (similar 3-D arrangements of its residues across the PDB); opens the browser and returns the URL."""
    obj, positions = _shift_obj(obj, positions)
    obj = _resolve_object(obj)
    if pocket not in (None, ""):
        pk = _find_pocket(_resolve_uid(uid, obj), pocket)
        positions = pk["resid"] if pk else None
    if isinstance(positions, str):
        positions = [p for p in re.split(r"[\s,+]+", positions) if p]
    url = motif_search_url(obj, positions or [], uid)
    if not url:
        print("[UFV] motif search unavailable (needs >= 2 mapped residues on a PDB entry or AlphaFold/AlphaFill model).")
        return None
    if int(browser):
        _open_url(url)
    print("[UFV] RCSB strucmotif: %s" % url)
    return url


def ligand_clash(obj, resn, chain, resi):
    """analysis.js ligandClash: closest heavy-atom contact between a ligand copy (HETATM) and the protein
    (non-HETATM) — severe < 2.0 Å, suspicious < 2.6 Å, else normal. -> {minDist, band} or None."""
    try:
        lsel = "(%s) and hetatm and resn %s and resi %s%s and not elem H+D" % (
            obj, resn, _resi_tok(resi), (" and chain %s" % chain) if chain else "")
        st = _obj_cur_state(obj)
        lc = cmd.get_coords(lsel, st)
        pc = cmd.get_coords("(%s) and not hetatm and not elem H+D" % obj, st)
    except Exception:
        return None
    if lc is None or pc is None or not len(lc) or not len(pc):
        return None
    np = _numpy()
    if np is not None:
        best = float("inf")
        for i in range(0, len(lc), 64):
            d = ((pc[None, :, :] - lc[i:i + 64, None, :]) ** 2).sum(-1).min()
            best = min(best, float(d))
    else:
        best = min((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 for a in lc for b in pc)
    d = best ** 0.5
    return {"minDist": d, "band": "severe" if d < 2.0 else "suspicious" if d < 2.6 else "normal"}


def ligand_panel_data(obj, resn, chain=None, resi=None, uid=None):
    """Main-thread inputs of the ligand panel beyond chemistry: the pocket (residues within the Nearby
    distance of this copy, from the current ligand focus when it matches), its evidence, the closest-
    contact clash band, AlphaFill transplant meta of the copy's chain, and whether a motif search works.
    -> kwargs for format_ligand_html."""
    uid = _resolve_uid(uid, obj)
    ann = _CACHE.get(uid) or {}
    f = _FOCUS_SHOWN.get(obj) or {}
    near = f.get("near") if f.get("ligand") == (resn, chain, resi) else None
    if near is None:
        lig = ("(%s) and resn %s%s%s" % (obj, resn, (" and chain %s" % chain) if chain else "",
                                         (" and resi %s" % _resi_tok(resi)) if resi is not None else ""))
        near, _p = _near_from_selection(obj, "(%s) and byres ((%s) around %g) and polymer%s" % (
            obj, lig, NEARBY_SETTINGS["dist"], _not_transplant(obj)))
    ev = pocket_evidence(near, _ca_geometry(obj), ann.get("sites") or [], ann.get("sequence") or "") if near else None
    clash = ligand_clash(obj, resn, chain, resi) if resi is not None else None
    aff = (_AFF_META.get(obj) or {}).get(chain) if chain is not None else None
    return {"aff": aff, "evidence": ev, "clash": clash, "near": near, "nearby_dist": NEARBY_SETTINGS["dist"],
            "motif": bool(near) and len(near) >= 2 and _structure_motif_ok(obj),
            "experimental": not _structure_is_predicted(obj), "sequence": ann.get("sequence") or ""}


# ---- Nearby (viewer-molstar.js _focusNeighbourhood) --------------------------------------------
NEARBY_SETTINGS = {"dist": 5.0}      # Å, the residue-panel slider (3-15)
NEARBY_MIN, NEARBY_MAX = 3.0, 15.0


def _near_from_selection(obj, nb_sel):
    """UniProt positions of the residues in a neighbourhood selection + residues on OTHER-protein chains
    [(chain, resi)] (not mapped to our UniProt; only for structures whose own chains are known)."""
    rows = []
    try:      # only `obj`'s own atoms count (another overlapping object's residues are not this structure's)
        cmd.iterate("(%s) and (%s) and name CA" % (obj, _state_sel(obj, nb_sel, "_ufv_nbsel")),
                    "rows.append((chain, resi))", space={"rows": rows})
    except Exception:
        return [], []
    mp = UFV_MAPS.get(obj) or {}
    ours = set(_mapped_chains(obj)) if mp.get("mode") in ("sifts", "custom", "manual") else None
    near, partner = set(), []
    seen = set()
    for ch, resi in rows:
        if (ch, resi) in seen:
            continue
        seen.add((ch, resi))
        u = _resi_to_uni(obj, ch, resi)
        if u is not None:
            near.add(u)
        elif ours and ch not in ours:
            partner.append((ch, resi))
    partner.sort(key=lambda x: (x[0], _resi_int(x[1], False) or 0))
    return sorted(near), partner


@_undo_paused       # read-only; at top level this also heals lost object names first (_heal_object_names)
def residue_nearby(obj, uni_pos, dist=None):
    """Residues with any (non-HETATM) atom within `dist` Å of the residue (AlphaFill transplant chains
    excluded) -> {"dist", "near": [uniPos], "partners": [(chain, resi, annotation or None)]}."""
    d = float(dist if dist is not None else NEARBY_SETTINGS["dist"])
    f = _FOCUS_SHOWN.get(obj) or {}
    if f.get("pos") == int(uni_pos) and f.get("dist") == d and f.get("nb"):
        near, partner = _near_from_selection(obj, f["nb"])     # reuse the focus' neighbourhood (no search)
    else:
        rsel = _sel_for_positions(obj, [int(uni_pos)])
        if not rsel:
            return {"dist": d, "near": [], "partners": []}
        nb = "byres (((%s) and polymer%s) within %g of (%s))" % (obj, _not_transplant(obj), d, rsel)
        near, partner = _near_from_selection(obj, nb)
    near = [p for p in near if p != int(uni_pos)]
    pann = (UFV_MAPS.get(obj) or {}).get("_partnerAnn")
    plist = []
    for ch, resi in partner:
        info = None
        if pann:
            r = _resi_int(resi)
            p = next((x for x in pann if x["chainId"] == ch and r in x["byResi"]), None)
            if p is None:
                continue                         # extension: annotated partner neighbours only
            info = dict(p["byResi"][r], accession=p["accession"])
        plist.append((ch, resi, info))
    return {"dist": d, "near": near, "partners": plist}


@_undo_paused
def _refocus_all(objs=None):
    """Re-run the current residue / ligand focus of each object (e.g. after the Nearby radius changed)."""
    names = _object_set()
    for o in (objs if objs is not None else list(_FOCUS_SHOWN)):
        f = _FOCUS_SHOWN.get(o) or {}
        if o not in names:
            continue
        try:
            if f.get("ligand"):
                rn, ch, ri = f["ligand"]
                ufv_ligand_focus(o, rn, ch, ri, f.get("uid"), zoom=0)
            elif f.get("pos") is not None:
                ufv_focus(o, f["pos"], f.get("uid"), zoom=0)
        except Exception:
            pass


@_undo_paused
def ufv_nearby(dist=None, obj=None):
    """ufv_nearby [distance [, object]] — the Nearby radius (3-15 Å, default 5) used by residue / ligand
    focus, the report's Nearby list and the ligand pocket evidence; re-focuses the current focus."""
    if dist in (None, ""):
        print("[UFV] Nearby distance: %g Å" % NEARBY_SETTINGS["dist"])
        return NEARBY_SETTINGS["dist"]
    NEARBY_SETTINGS["dist"] = max(NEARBY_MIN, min(NEARBY_MAX, float(dist)))
    _refocus_all([_resolve_object(obj)] if obj else None)
    print("[UFV] Nearby distance set to %g Å." % NEARBY_SETTINGS["dist"])
    return NEARBY_SETTINGS["dist"]


# ---- per-residue ligand contacts (CSV nearby_ligands) -------------------------------------------
LIGAND_CONTACT_CUTOFF = 5.0


def _ligand_contacts(obj, cutoff=LIGAND_CONTACT_CUTOFF):
    """{UniProt position: sorted CCD codes of ligands (non-polymer, non-water) with any atom within
    `cutoff` Å of the residue}. The extension's CSV column is a stub (always empty: viewer-molstar.js
    ligandContactsByResidue returns an empty Map); this fills it. Main thread (cmd); cached per object."""
    key = (obj, None, "ligcontacts", float(cutoff))
    hit = _ANALYSIS_CACHE.get(key)
    if hit is not None:
        return hit
    out = {}
    if _has_ligands(obj):
        lsel = _ligand_sel(obj)
        lig, prot = [], []
        try:
            st = _obj_cur_state(obj)
            cmd.iterate_state(st, lsel, "lig.append((resn, x, y, z))", space={"lig": lig})
            cmd.iterate_state(st, _state_sel(obj, "(%s) and polymer%s and ((%s) within %g of (%s))" % (
                obj, _not_transplant(obj), obj, cutoff, lsel), "_ufv_lcsel"),
                "prot.append((chain, resi, x, y, z))", space={"prot": prot})
        except Exception:
            lig, prot = [], []
        if lig and prot:
            c2 = cutoff * cutoff
            np = _numpy()
            names = [r[0] for r in lig]
            res_of = {}
            if np is not None:
                L = np.array([r[1:] for r in lig], dtype=float)
                P = np.array([r[2:] for r in prot], dtype=float)
                for i in range(0, len(P), 256):
                    d2 = ((P[i:i + 256, None, :] - L[None, :, :]) ** 2).sum(-1)
                    for j, k in zip(*np.nonzero(d2 <= c2)):
                        res_of.setdefault((prot[i + j][0], prot[i + j][1]), set()).add(names[k])
            else:
                for ch, resi, x, y, z in prot:
                    for nm, a, b, c in lig:
                        if (x - a) ** 2 + (y - b) ** 2 + (z - c) ** 2 <= c2:
                            res_of.setdefault((ch, resi), set()).add(nm)
            for (ch, resi), names_ in res_of.items():
                u = _resi_to_uni(obj, ch, resi)
                if u is not None:
                    out.setdefault(u, set()).update(names_)
    out = {u: sorted(v) for u, v in out.items()}
    _ANALYSIS_CACHE[key] = out
    return out


# ---- residue flags ("bulbs", modal.js renderDetails title) -------------------------------------
FLAG_DEFS = (("hotspot", "Pathogenic variant hotspot", "#e53935"),
             ("burden", "Recurrent phenotype residue", "#e65100"),
             ("hub", "Contact-network centrality", "#6a1b9a"),
             ("pocket", "Burial-adjusted constraint cluster", "#00897b"))


def _residue_flags(obj, uid, pos, ann=None):
    """[(key, label, colour, flagged)] for the analyses already computed (memoised — nothing is computed
    here); an analysis not run yet is left out. Constraint pocket = class at the current q threshold."""
    ann = ann or _CACHE.get(uid) or {}
    out = []
    hot = _analysis_cached(obj, uid, "hotspots") if obj else None
    hub = _analysis_cached(obj, uid, "hubs") if obj else None
    pk = _analysis_cached(obj, uid, "pockets") if obj else None
    states = {"hotspot": None if hot is None else pos in hot,
              "burden": None if ann.get("burden") is None else pos in ann["burden"],
              "hub": None if hub is None else pos in hub,
              "pocket": None if pk is None else pos in _pocket_classes(pk)}
    for key, label, color in FLAG_DEFS:
        if states[key] is not None:
            out.append((key, label, color, states[key]))
    return out


# ---- exports: B-factor-encoded PDB (export.js rewritePdbBeta), PyMOL selection (formatSelection) -----
BFACTOR_MODES = ("selected", "alphamissense", "hotspots", "hubs", "burden", "pocket", "plddt", "bfactor")


def _bfactor_values(obj, uid, mode):
    """{uniPos: value} for a B-factor export mode, or None for keep-the-file-value modes (pLDDT/B)."""
    ann = _CACHE.get(uid) or fetch_annotations(uid)
    if mode in ("plddt", "bfactor"):
        return None
    if mode in ("rmsf", "ligocc", "partnerocc"):        # trajectory results (ufv_traj_rmsf / _contacts)
        return _traj_bfactor_values(obj, mode)
    if mode in ("alphamissense", "am"):
        return {int(p): float(v) for p, v in (ann.get("amMean") or {}).items()}
    if mode == "hotspots":
        return {p: {"strong": 1.0, "moderate": 0.67, "weak": 0.33}.get(t, 0.0)
                for p, t in (_analysis(obj, uid, "hotspots") or {}).items()}
    if mode in ("hubs", "contacthubs"):
        return {p: {"strong": 1.0, "moderate": 0.5}.get(t, 0.0) for p, t in (_analysis(obj, uid, "hubs") or {}).items()}
    if mode == "burden":
        return {int(p): 1.0 for p in ann.get("burden") or ()}
    if mode in ("pocket", "pockets", "prism"):
        return {p: {"pocket": 1.0, "exposed": 0.66}.get(c, 1.0)
                for p, c in _pocket_classes(_analysis(obj, uid, "pockets")).items()}
    # "selected": residues currently drawn by the active sphere layers
    return {int(p): 1.0 for p in _displayed_positions(obj)}


@_undo_paused
def ufv_save_bfactor(filename=None, obj=None, mode="selected", uid=None):
    """ufv_save_bfactor [filename [, object [, mode [, uid]]]] — save the structure as PDB with an
    analysis encoded in the B-factor column (export.js rewritePdbBeta): alphamissense (mean score 0-1),
    hotspots (1 / 0.67 / 0.33), hubs (1 / 0.5), burden (1), pocket (1 buried / 0.66 exposed, current q),
    plddt / bfactor (file values kept), selected (1 = residues in the active marker layers). Atoms not
    mapped to the UniProt entry get 0. The object's own B-factors are restored afterwards."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    mode = (mode or "selected").strip().lower()
    if mode == "am":
        mode = "alphamissense"
    vals = _bfactor_values(obj, uid, mode)
    filename = filename or os.path.join(os.getcwd(), "%s_%s_bfactor.pdb" % (re.sub(r"[^A-Za-z0-9_.-]", "_", obj), mode))
    keys, saved = set(), {}
    cmd.iterate(obj, "keys.add((chain, resi)); saved[index] = b", space={"keys": keys, "saved": saved})
    ours = set(_mapped_chains(obj))
    m = {}
    for ch, resi in keys:
        u = _resi_to_uni(obj, ch, resi) if ch in ours else None
        if u is None:
            m[(ch, resi)] = 0.0
        elif vals is not None:
            m[(ch, resi)] = round(float(vals.get(u, 0.0)), 2)
    try:
        if vals is None:      # keep the file's pLDDT / B for mapped residues
            cmd.alter(obj, "b = _m.get((chain, resi), b)", space={"_m": m})
        else:
            cmd.alter(obj, "b = _m.get((chain, resi), 0.0)", space={"_m": m})
        cmd.save(filename, obj)
    finally:
        cmd.alter(obj, "b = _s.get(index, b)", space={"_s": saved})
    print("[UFV] wrote %s (B-factor = %s)." % (filename, mode))
    return filename


@_undo_paused
def ufv_select(name="ufv_sel", obj=None, layer=None):
    """ufv_select [name [, object [, layer]]] — a PyMOL selection of the residues the marker layers show
    (all active layers, or one: ptms | variants | sites | mutagenesis | filtered | binding), like the
    extension's "copy selection"; prints the `resi ...` string."""
    obj = _resolve_object(obj)
    pos = _displayed_positions(obj, (layer or "").strip().lower() or None)
    sel = _sel_for_positions(obj, pos)
    if not sel:
        print("[UFV] nothing displayed to select.")
        return None
    cmd.select(name, sel, enable=0)
    rs = sorted({int(r) for ch in _mapped_chains(obj) for r in (_uni_to_resi(obj, ch, p) for p in pos) if r is not None})
    print("[UFV] %s: %d residue(s) — resi %s" % (name, len(rs), "+".join(str(r) for r in rs)))
    return name


def _residue_csv(ann, obj=None, protvar_by_pos=None, geom=None, chains=None, lig_contacts=None):
    """Per-residue annotation matrix, the extension's buildResidueMatrix (export.js) column for column:
    position, aa, pdb_residue, one-hot ptm_<category> / disease_<label> (diseasePairs labels), site flags,
    mutagenesis, n_variants, n_pathogenic, gnomad_max_af, [ProtVar], am_avg/max/n (19 WT substitutions),
    hotspot_tier, contact_hub_tier, residue_burden, constraint_pocket_q / _class (EVERY candidate),
    PTM-variant proximity (tier, nearest variant, distance, counts at the current radii), nearby_ligands
    (CCD codes within 5 Å, ';'-joined — the extension's column is an always-empty stub). A structure mapped on several chains of OUR protein gets the multi-chain
    layout (pdb_residue_<C>, hotspot_tier_<C>, contact_hub_tier_<C> at the end). From a worker thread
    pass `geom`, `chains` and `lig_contacts` (_ligand_contacts) read on the main thread; tiers come from
    the shared analysis cache."""
    if not ann.get("sequence", ""):
        return ""
    uid = ann.get("uid")
    am_map = ann.get("amMap") or {}
    hot, hub, pocket, prox, hot_ch, hub_ch = {}, {}, {}, {}, {}, {}
    if chains is None:
        chains = _mapped_chains(obj) if obj else []
    if obj:
        if geom is None:
            geom = _ca_geometry(obj)
        for kind in ("hotspots", "hubs"):
            try:
                r = _analysis(obj, uid, kind, geom, ann)
                rc = _analysis(obj, uid, kind + "_bychain", geom, ann)
            except Exception:
                r, rc = {}, {}
            if kind == "hotspots":
                hot, hot_ch = r, rc
            else:
                hub, hub_ch = r, rc
        if am_map:
            try:
                pocket = (_analysis(obj, uid, "pockets", geom, ann) or {}).get("byPos") or {}
            except Exception:
                pass
        if ann["ptms"] and ann["variants"]:
            try:
                prox = _analysis(obj, uid, "ptmprox", geom, ann)
            except Exception:
                pass
    if lig_contacts is None and obj:
        lig_contacts = _ANALYSIS_CACHE.get((obj, None, "ligcontacts", LIGAND_CONTACT_CUTOFF))
        if lig_contacts is None and threading.current_thread() is threading.main_thread():
            try:
                lig_contacts = _ligand_contacts(obj)
            except Exception:
                lig_contacts = None
    return _residue_csv_text(ann, obj, protvar_by_pos, chains, hot, hub, pocket, prox, hot_ch, hub_ch,
                             lig_contacts)


@_undo_paused
def ufv_csv(filename=None, obj=None, uid=None, protvar=0):
    """ufv_csv [filename [, object [, uniprot_id [, protvar]]]] — write the per-residue annotation
    matrix CSV (one-hot PTM/disease/site/mutagenesis + variant counts, gnomAD, AlphaMissense, burden,
    and structure tiers). protvar=1 also adds ProtVar conservation/EVE/ESM1b/M3D (one request per
    residue — slow)."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid, obj))
    pv = None
    if int(protvar):
        print("[UFV] fetching ProtVar for %d residues ..." % len(ann.get("sequence", "")))
        pv = _fetch_protvar_all(ann["uid"], len(ann["sequence"]),
                                lambda d, t: print("[UFV] ProtVar %d/%d" % (d, t)) if d % 200 == 0 else None)
    text = _residue_csv(ann, obj, pv)
    if not text:
        print("[UFV] no sequence to export.")
        return None
    filename = filename or os.path.join(os.getcwd(), "%s_residue_annotations.csv" % ann["uid"])
    with open(filename, "w") as fh:
        fh.write(text)
    print("[UFV] wrote %s" % filename)
    return filename


@_undo_paused       # read-only; at top level this also heals lost object names first (_heal_object_names)
def residue_report(obj, uid, uni_pos):
    """Everything the extension's detail panel shows for one residue, as a plain dict."""
    ann = fetch_annotations(uid)
    geom = _ca_geometry(obj)
    seq = ann["sequence"]
    uni_pos = int(uni_pos)
    aa = seq[uni_pos - 1] if 0 < uni_pos <= len(seq) else "?"
    var_pos = ann["varByPos"]
    # PTM–variant proximity, residue-centric (analysis.js computeResidueProximity): PTMs within the PTM
    # radius and pathogenic variants within the variant radius of THIS residue, nearest first.
    rp = _residue_proximity(uni_pos, ann["ptms"], ann["variants"], geom)
    prox_ptms = [{"pos": p["position"], "dist": d, "category": p.get("category"), "color": p.get("color"),
                  "description": p.get("description")} for p, d in rp["ptms"]]
    prox_vars = [{"pos": v["position"], "dist": d, "label": _vlabel(v), "consequence": v.get("consequence"),
                  "color": v.get("consequenceColor"), "diseases": v.get("diseases") or []} for v, d in rp["variants"]]
    # ligands/cofactors within 5 Å of this residue — skip the spatial query entirely when the
    # structure has no ligands (the common AlphaFold case), so each click stays cheap.
    ligs = []
    if _has_ligands(obj):
        rsel = _sel_for_positions(obj, [uni_pos])
        if rsel:
            try:
                cmd.iterate(_state_sel(obj, "(%s) and not polymer and not solvent and ((%s) around 5)" % (obj, rsel),
                                       "_ufv_rlsel"), "ligs.append(resn)", space={"ligs": ligs})
            except Exception:
                pass
    # PTM-centric tier (only when this residue is itself a PTM — cheap, one feature), current radii.
    proximity = None
    ptms_here = [p for p in ann["ptms"] if p["position"] == uni_pos]
    if ptms_here:
        proximity = _ptm_variant_proximity(ptms_here, ann["variants"], geom).get(uni_pos)
    nb = residue_nearby(obj, uni_pos)
    want = set(nb["near"])
    var_color = {}
    for v in ann["variants"]:                  # modal.js buildAnnotationMap: last variant at a position wins
        if v["position"] in want:
            var_color[v["position"]] = v.get("consequenceColor")
    return {
        "uid": ann["uid"], "pos": uni_pos, "aa": aa,
        "proximity": proximity,
        "proxPtms": prox_ptms, "proxVariants": prox_vars,
        "proxRadii": (PROX_SETTINGS["ptm"], PROX_SETTINGS["var"]),
        "ptms": [p for p in ann["ptms"] if p["position"] <= uni_pos <= p["endPosition"]],
        "variants": var_pos.get(uni_pos, []),
        "sites": [s for s in ann["sites"] if s["position"] <= uni_pos <= s["endPosition"]],
        "domains": [d for d in ann["domains"] if d["position"] <= uni_pos <= d["endPosition"]],
        "mutagenesis": [mtg for mtg in ann.get("mutagenesis", []) if mtg["position"] <= uni_pos <= mtg["endPosition"]],
        "amMean": ann["amMean"].get(uni_pos),
        "amProfile": ann.get("amByPos", {}).get(uni_pos, []),  # full 19-substitution profile
        "nearbyLigands": sorted(set(ligs)),
        "obj": obj,
        "nearby": nb,
        "flags": _residue_flags(obj, ann["uid"], uni_pos, ann),
        "pdbekb": _PDBEKB_CACHE.get(_canon_acc(ann["uid"]), "pending"),   # filled by fetch_pdbekb (worker)
        "varColor": var_color,
    }


_DIMMED = set()  # objects whose cartoon is currently dimmed for focus


@_undo_paused
def _clear_focus(obj, restore_spheres=True):
    """Remove the current focus sticks, restore the focused residue's own atom colours, and re-enable
    the overlay objects the focus hid (only those that are still switched on)."""
    prev = _FOCUS_SHOWN.pop(obj, None)
    if not prev:
        return
    try:
        cmd.hide("sticks", prev["sticks"])
        cmd.unset("stick_color", prev["sticks"])
        if prev.get("surface"):
            cmd.hide("surface", prev["surface"])
        if prev.get("lig"):
            cmd.unset("stick_radius", prev["lig"])  # un-fatten the previously focused ligand copy
        if prev.get("sel"):
            cmd.hide("spheres", prev["sel"])
            cmd.unset("sphere_color", "(%s) and name CA" % prev["sel"])
            if prev.get("colors"):
                # Restore the residue's own atom colours (it used to be painted gray80, punching a
                # hole in the active cartoon colouring). Only non-CA atoms were recoloured and their
                # sticks/spheres are hidden now, so no recolor() of the whole object is needed.
                cmd.alter(prev["sel"], "color = _m.get(index, color)", space={"_m": prev["colors"]})
    except Exception:
        pass
    for nm in prev.get("selnames", []):       # drop the named focus selections we created
        try:
            cmd.delete(nm)
        except Exception:
            pass
    if restore_spheres and prev.get("disabled"):
        if obj not in _enabled_set():         # the structure itself was hidden meanwhile
            return
        for lobj in prev["disabled"]:
            if _layer_wanted(obj, lobj):      # never re-show a layer the user switched off meanwhile
                try:
                    cmd.enable(lobj)
                except Exception:
                    pass


@_undo_paused
def ufv_focus(obj=None, uni_pos=None, uid=None, zoom=1, declutter=1):
    """ufv_focus [object,] uniprot_pos — zoom into a residue: show it and its 5 Å neighbourhood as
    sticks coloured by annotation, and hide the sphere markers there. (The cartoon is NOT made
    transparent: transparency forces a full cartoon rebuild and is expensive to render every frame,
    which made zooming slow and the zoomed-in view laggy to rotate/scroll.) Returns the name of the
    focus selection (zoom=0 lets a caller frame several focused structures at once). declutter=0 keeps the
    annotation sphere layers shown (the panel's single-click highlight: what is checked stays drawn)."""
    obj, uni_pos = _shift_obj(obj, uni_pos)
    obj = _resolve_object(obj)
    ann = _CACHE.get(_resolve_uid(uid, obj)) or {}
    sel = _sel_for_positions(obj, [int(uni_pos)])
    if not sel:
        print("[UFV] position %s not modelled in %s." % (uni_pos, obj))
        return None
    # Materialise the residue + its 5 Å neighbourhood as NAMED selections (per object, so the
    # multi-structure focus loop doesn't collide on one global name). PyMOL stores a named
    # selection as a fixed atom set, so the `around 5` distance search — costly on big models —
    # runs ONCE here instead of being re-evaluated by every show / colour / hide of the string.
    fnb = _ufv_layer_obj(obj, "focusnb")
    fset = _ufv_layer_obj(obj, "focusset")
    cmap = _annotation_color_map(ann)
    nd = NEARBY_SETTINGS["dist"]
    selnames = [fset, fnb]
    with _batch():
        _clear_focus(obj)                       # clear previous focus sticks
        fst = _obj_cur_state(obj)               # trajectories: the frame shown, not "near in any state"
        # scoped to the residue's own object (partner chains of it still count): an overlapping model loaded
        # alongside (aligned AF / PDB) must not add residues — same fix as ufv_ligand_focus' pocket
        cmd.select(fnb, "(%s) and byres ((%s) around %g) and polymer%s" % (obj, sel, nd, _not_transplant(obj)),
                   enable=0, state=fst)
        if _has_ligands(obj) and obj not in _LIG_SHOWN:
            # pocket ligands within nearby + 3 Å stay visible while zoomed in (viewer-molstar.js)
            flig = _ufv_layer_obj(obj, "focusnlig")
            cmd.select(flig, "byres ((%s) and not polymer and not solvent and ((%s) within %g of (%s)))"
                       % (obj, obj, nd + 3, sel), enable=0, state=fst)
            selnames.append(flig)
            cmd.select(fset, "(%s) or (%s) or (%s)" % (sel, fnb, flig), enable=0)
        else:
            cmd.select(fset, "(%s) or (%s)" % (sel, fnb), enable=0)
        cmd.show("sticks", fset)
        cmd.set("stick_color", "grey70", fset)
        _color_sticks_by_annotation(obj, fnb, cmap)
        if len(selnames) > 2:
            cmd.unset("stick_color", selnames[2])      # pocket ligands keep element colours
        # The selected residue: ball-and-stick with atom-type colours. Carbons take the residue's
        # DISEASE annotation colour (variant/PTM/site) or neutral grey — never the active cartoon
        # colouring (hub/burden/pocket), which we explicitly override here so it can't bleed onto the
        # sidechain. Heteroatoms keep element colours. The CA itself is NOT recoloured (the cartoon
        # takes its colour from CA, so recolouring it punched a hole in the cartoon layer); it gets
        # the colour through atom-level stick/sphere colour settings instead. Original atom colours
        # are saved and restored by _clear_focus.
        dcol = cmap.get(int(uni_pos))
        ccol = _color_name(dcol) if dcol else "gray85"
        saved = {}
        cmd.iterate(sel, "saved[index] = color", space={"saved": saved})
        cmd.show("spheres", sel)
        cmd.set("sphere_scale", 0.28, sel)
        cmd.color("atomic", "(%s) and (not elem C)" % sel)
        cmd.color(ccol, "(%s) and elem C and not name CA" % sel)
        cmd.unset("stick_color", sel)  # let element/disease colours show on the selected residue's sticks
        cmd.set("stick_color", ccol, "(%s) and name CA" % sel)
        cmd.set("sphere_color", ccol, "(%s) and name CA" % sel)
        # Hide the annotation SPHERE layers while zoomed in (declutter). Keep the ligands layer shown —
        # a residue near a ligand should still see it as context. Only layers that are currently
        # visible are hidden (and therefore re-shown later).
        disabled = []
        enabled = _enabled_set() if int(declutter) else ()
        for tag in ("ptm", "var", "site", "dom", "domr", "filt", "mutag"):
            lobj = _ufv_layer_obj(obj, tag)
            if lobj in enabled:
                try:
                    cmd.disable(lobj)
                    disabled.append(lobj)
                except Exception:
                    pass
    _FOCUS_SHOWN[obj] = {"sticks": fset, "sel": sel, "disabled": disabled, "selnames": selnames,
                         "colors": saved, "pos": int(uni_pos), "uid": uid, "dist": float(nd), "nb": fnb}
    if int(zoom):
        _focus_zoom(fset, 3)
    return fset


@_undo_paused
def ufv_resetview(obj=None):
    """ufv_resetview [object] — clear focus sticks, restore spheres + cartoon, and go back to the view
    from before the focus (exactly, when the camera was not moved since the focus zoom); otherwise, or
    with no focus, zoom out to the whole object."""
    obj = _resolve_object(obj)
    with _batch():
        _clear_focus(obj)
        if obj in _DIMMED:
            try:
                cmd.set("cartoon_transparency", 0.0, obj)
            except Exception:
                pass
            _DIMMED.discard(obj)
    pre, post = _FOCUS_VIEW["pre"], _FOCUS_VIEW["post"]
    _FOCUS_VIEW["pre"] = _FOCUS_VIEW["post"] = None
    if pre is not None and post is not None and tuple(cmd.get_view()) == post:
        cmd.set_view(pre)
    elif obj:
        cmd.zoom(obj, animate=0)
    _view_commit()


@_undo_paused
def ufv_report(uni_pos, obj=None, uid=None):
    """ufv_report uniprot_pos [, object [, uid]] — print the residue report + focus the residue."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    ufv_focus(obj, uni_pos, uid)             # first: the report reuses the focus' neighbourhood selection
    rep = residue_report(obj, uid, uni_pos)
    print(format_report(rep))
    return rep


def _loaded_structures():
    return [o for o in (cmd.get_object_list() or []) if not o.startswith("ufv_")]


def _uni_ca(obj):
    """{UniProt position: (x, y, z)} of obj's modelled Cα atoms (first chain in model order wins)."""
    out = {}
    for g in _ca_geometry(obj, refresh=True):
        if g["uniPos"] is not None and g["uniPos"] not in out:
            out[g["uniPos"]] = (g["x"], g["y"], g["z"])
    return out


@_undo_paused
def ufv_align(reference=None):
    """ufv_align [reference] — superpose all loaded protein structures onto a reference (default the
    first loaded), so several structures of the protein overlay. Like the extension's Overlay, each is
    fitted on the Cα atoms of the residues both share, matched by UniProt position (no sequence
    alignment); objects of another protein or with < 10 shared residues fall back to cmd.align."""
    objs = _loaded_structures()
    if len(objs) < 2:
        print("[UFV] need at least 2 loaded structures to align.")
        return
    ref = reference or objs[0]
    ref_uid = (_resolve_uid(None, ref) or "").upper() if ref in UFV_MAPS else ""
    ref_ca = _uni_ca(ref) if ref_uid else {}
    for o in objs:
        if o == ref:
            continue
        try:
            same = ref_uid and o in UFV_MAPS and (_resolve_uid(None, o) or "").upper() == ref_uid
            ca = _uni_ca(o) if same else {}
            shared = sorted(set(ref_ca) & set(ca))
            if len(shared) >= 10:
                m, rmsd = fit_rigid([ref_ca[p] for p in shared], [ca[p] for p in shared])
                cmd.transform_object(o, [x for row in m for x in row], homogenous=1)
                _GEOM_CACHE.pop(o, None)
                print("[UFV] aligned %s onto %s: Cα RMSD %.2f Å over %d shared UniProt residues." % (o, ref, rmsd, len(shared)))
            else:
                r = cmd.align(o, ref)
                print("[UFV] aligned %s onto %s (sequence alignment): RMSD %.2f Å over %d atoms." % (o, ref, r[0], r[1]))
        except Exception as e:
            print("[UFV] align %s failed: %s" % (o, e))


@_undo_paused
def ufv_load_annotations(path, obj=None, uid=None):
    """ufv_load_annotations file [, object [, uniprot_id]] — load your own residue annotations (CSV / TSV) as the
    Custom annotation category and show them as spheres. Columns: position (123, R123W, p.Arg123Trp or
    120-130), optional label, value, color, note; without a header row: position, label, value. Kept for
    the session; added to ufv_csv as the custom_annotation column. Same format as the browser extension."""
    obj = _resolve_object(obj)
    uid = (_resolve_uid(uid, obj) or "").upper()
    s = custom_load(uid, path)
    print("[UFV] %s: %d custom residue%s from %s%s; groups: %s." % (
        uid, s["residues"], "" if s["residues"] == 1 else "s", s["file"],
        (" (%d row%s skipped: no valid position)" % (s["skipped"], "" if s["skipped"] == 1 else "s")) if s["skipped"] else "",
        ", ".join("%s %d" % (g["label"], g["count"]) for g in s["groups"]) or "none"))
    if obj:
        _ann_set(obj, uid, "cus", ["cus:%d" % i for i in range(s["residues"])])
        _ann_apply(obj, uid, ("cus",))


@_undo_paused
def ufv_unload_annotations(obj=None, uid=None):
    """ufv_unload_annotations [object [, uniprot_id]] — remove the custom annotations."""
    obj = _resolve_object(obj)
    uid = (_resolve_uid(uid, obj) or "").upper()
    if obj:
        _ann_set(obj, uid, "cus", [])
        _ann_apply(obj, uid, ("cus",))
    print("[UFV] %s: custom annotations %s." % (uid, "removed" if custom_clear(uid) else "were not loaded"))


@_undo_paused
def ufv_structures(uid=None):
    """ufv_structures [uniprot_id] — list the structures available for an accession."""
    uid = (_resolve_uid(uid) or "").upper()
    fetch_annotations(uid)
    structs = list_structures(uid)
    print("[UFV] %d structures for %s:" % (len(structs), uid))
    for s in structs:
        print("    %-26s %s" % (s["key"], s["label"]))
    _STATE["structures"] = structs
    _STATE["structures_uid"] = uid
    return structs


@_undo_paused
def ufv_use(key=None, uid=None):
    """ufv_use key [, uniprot_id] — download and load a structure listed by ufv_structures,
    setting its residue numbering automatically (AlphaFold/computed = identity, PDB = SIFTS)."""
    uid = (_resolve_uid(uid) or "").upper()
    structs = _STATE.get("structures") if _STATE.get("structures_uid") == uid else None
    structs = structs or ufv_structures(uid)       # never pick from another accession's list
    match = next((s for s in structs if s["key"].lower() == str(key).lower()), None)
    if not match:
        print("[UFV] no structure '%s' — run ufv_structures to list them." % key)
        return
    return _load_structure(match, uid)


def _struct_obj_name(struct):
    return re.sub(r"[^A-Za-z0-9_]", "_", struct["key"])


@_undo_paused
def _load_structure(struct, uid, name=None):
    print("[UFV] downloading %s ..." % struct["url"])
    prep = _prepare_structure(struct, uid)       # model + numbering metadata fetched in parallel
    if not prep:
        print("[UFV] download failed.")
        return None
    name = _load_from_path(struct, prep["path"], uid, prep["segments"], prep["hints"], name=name,
                           meta=prep["meta"])
    print("[UFV] loaded %s (%s). Numbering: %s." % (name, struct.get("label", name), struct.get("numbering")))
    return name


@_undo_paused
def ufv_hide(obj=None, layer=None):
    """ufv_hide [object [, layer]] — hide one layer: ptms | variants | sites | mutagenesis | ligands |
    filtered, or a cartoon colouring: domains | topology | alphamissense | burden | plddt | bfactor |
    hotspots | contacthubs | pockets | cartoon. With no layer, hides every UFV overlay (ufv_clear).
    `ufv_hide ptms` (a layer name alone) hides that layer on the current object."""
    if obj and not layer and obj not in cmd.get_names("objects"):
        key = obj.strip().lower()
        if (key in _HIDE_POINT_ALIASES or key in CARTOON_TAGS
                or key in ("ligands", "ligand", "lig", "cartoon", "hubs", "pocket")):
            obj, layer = None, key
    obj = _resolve_object(obj)
    layer = (layer or "").strip().lower()
    if not layer:
        ufv_clear(obj)
        return
    if layer in ("ligands", "ligand", "lig"):
        ufv_ligands_hide(obj)
    elif layer in CARTOON_TAGS or layer in ("cartoon", "hubs", "pocket"):
        _reset_cartoon(obj)
    elif layer in _HIDE_POINT_ALIASES:
        _hide_layer(obj, _HIDE_POINT_ALIASES[layer])
    else:
        print("[UFV] unknown layer '%s'." % layer)
        return
    print("[UFV] %s: hid %s layer." % (obj, layer))


@_undo_paused
def ufv_clear(obj=None):
    """ufv_clear [object] — remove all UFV overlays and reset colour."""
    obj = _resolve_object(obj)
    _clear_focus(obj, restore_spheres=False)
    _delete_layers(obj)
    _CARTOON_LAYER.pop(obj, None)
    if obj in _DIMMED:
        try:
            cmd.set("cartoon_transparency", 0.0, obj)
        except Exception:
            pass
        _DIMMED.discard(obj)
    if obj:
        cmd.color("gray80", "(%s) and polymer" % obj)   # keep ligand colours intact
    print("[UFV] cleared UFV overlays on %s." % obj)


def actual_coverage(obj, seqlen):
    """Fraction of the UniProt sequence actually modelled in `obj` (CA atoms with a mapped UniProt
    position) — the real coverage, like the extension's modelled-residue count (not the SIFTS range)."""
    if not seqlen:
        return None
    covered = {g["uniPos"] for g in _ca_geometry(obj) if g["uniPos"] is not None}
    return round(len(covered) / seqlen * 100, 1)


# ----------------------------------------------------------------------------------------------
# Mapping / loading commands
# ----------------------------------------------------------------------------------------------
@_undo_paused
def ufv_fetch(uid):
    """ufv_fetch uniprot_id — download and cache all annotation layers for an accession."""
    fetch_annotations(uid, force=True)


@_undo_paused
def ufv_map(obj, uid, mode="identity", pdb_id=None):
    """ufv_map object, uniprot_id [, mode [, pdb_id]]
    Establish how `object` residue numbers relate to UniProt positions.
      mode = identity        resi == UniProt position (AlphaFold / UniProt-numbered models)
      mode = sifts, <pdb>    map through PDBe/SIFTS (+ residue_listing) for the given PDB id
      mode = auto            per-chain Custom mapping, every chain auto-detected (= ufv_custom);
                             refine single chains with ufv_chainmap (offset / adopt / ignore / accession)"""
    obj = _resolve_object(obj)
    uid = uid.strip().upper()
    fetch_annotations(uid)
    if mode.lower() in ("auto", "custom"):
        _apply_custom_map(obj, uid, {}, reset=True)
        return
    if mode.lower() == "sifts":
        pid = pdb_id or uid  # allow `ufv_map obj, P35498, sifts, 7dtd`
        _set_sifts_map(obj, pid, uid)
    else:
        _set_identity_map(obj)
        print("[UFV] %s mapped with identity numbering (resi == UniProt position)." % obj)
    _STATE.setdefault("obj_uid", {})[obj] = uid
    _invalidate_obj(obj)


@_undo_paused
def ufv_chain(obj, chain, uniprot_start, resi_start=None, resi_end=None):
    """ufv_chain object, chain, uniprot_start [, resi_start [, resi_end]]
    Manually anchor one chain's numbering: residue `resi_start` (default = uniprot_start)
    corresponds to UniProt position `uniprot_start`, linearly. `resi_end` bounds the chain.
    Ideal for trajectories with custom numbering or to define each chain's start/finish."""
    obj = _resolve_object(obj)
    mp = _ensure_manual_map(obj)
    u_start = int(uniprot_start)
    r_start = int(resi_start) if resi_start not in (None, "", "None") else u_start
    r_end = int(resi_end) if resi_end not in (None, "", "None") else None
    mp["chains"][chain] = {"u_start": u_start, "resi_start": r_start, "resi_end": r_end}
    _invalidate_obj(obj)
    print("[UFV] %s chain %s: resi %d == UniProt %d%s"
          % (obj, chain, r_start, u_start, (" (..%d)" % r_end) if r_end else ""))


@_undo_paused
def ufv_custom(obj=None, uid=None):
    """ufv_custom [object [, uniprot_id]] — per-chain mapping for a structure you loaded yourself (the
    extension's "Custom PDB"): lists every polymer chain of `object` with its first/last residue, maps
    each chain automatically — a known SIFTS / isoform mapping is kept, otherwise the chain's offset is
    found by matching its sequence to the UniProt sequence; chains that don't match stay unmapped — and
    prints the per-chain table. Refine single chains with ufv_chainmap."""
    obj = _resolve_object(obj)
    uid = (_resolve_uid(uid, obj) or "").strip().upper()
    if not obj or not uid:
        print("[UFV] usage: ufv_custom object, uniprot_id")
        return None
    return _apply_custom_map(obj, uid, {}, reset=True)


@_undo_paused
def ufv_chainmap(obj=None, chain=None, mode="auto", value=None, accession=None, uid=None):
    """ufv_chainmap object, chain, mode [, value [, accession [, uniprot_id]]]
    Set one chain's residue mapping (Custom PDB). chain may be * for every chain. mode =
      auto     keep the chain's known mapping (SIFTS / isoform) or detect its offset from the sequence
      offset   value = UniProt position of the chain's FIRST residue (linear over first..last residue)
      adopt    value = a loaded object [/chain] whose mapping this chain copies, e.g. 6HUO_A_D/A
      ignore   no annotations on this chain
    accession assigns the chain to another UniProt entry: it becomes a partner chain and is not
    annotated with this protein's features (ufv_partners lists it).
    Examples:  ufv_chainmap traj, B, offset, 28      ufv_chainmap traj, C, ignore
               ufv_chainmap traj, A, adopt, 6HUO_A_D/A     ufv_chainmap traj, E, auto, , P18505"""
    obj = _resolve_object(obj)
    uid = (_resolve_uid(uid, obj) or "").strip().upper()
    mode = (mode or "auto").strip().lower()
    if not obj or not uid or chain in (None, ""):
        print("[UFV] usage: ufv_chainmap object, chain, auto|offset|adopt|ignore [, value [, accession]]")
        return None
    if mode not in CUSTOM_MODES:
        print("[UFV] unknown mode '%s' — use one of %s." % (mode, ", ".join(CUSTOM_MODES)))
        return None
    cfg = {"mode": mode}
    value = None if value in (None, "", "None") else str(value).strip()
    if mode == "offset":
        try:
            cfg["start"] = int(value)
        except (TypeError, ValueError):
            print("[UFV] offset needs the UniProt position of the chain's first residue.")
            return None
    elif mode == "adopt":
        if not value:
            print("[UFV] adopt needs a loaded object (optionally object/chain) to copy the mapping from.")
            return None
        src, _, src_chain = value.partition("/")
        if src not in _object_set():
            print("[UFV] no loaded object '%s'." % src)
            return None
        cfg["src"], cfg["src_chain"] = src, src_chain or None
    if accession not in (None, "", "None"):
        cfg["accession"] = str(accession).strip().upper()
    chains = sorted(_chain_extents(obj)) if str(chain).strip() == "*" else [str(chain).strip()]
    return _apply_custom_map(obj, uid, {c: dict(cfg) for c in chains})


@_undo_paused
def ufv_chimera(obj=None):
    """ufv_chimera [object] — highlight residues of OUR chains that have no UniProt mapping (chimeric /
    fusion partner segments, tags, linkers) in blue-grey #8fa5c0 — the extension's "Chimeric partner"
    layer. One alter pass; ufv_clear / another colouring removes it."""
    obj = _resolve_object(obj)
    mp = UFV_MAPS.get(obj) or {}
    if not obj or mp.get("mode") in (None, "identity"):
        print("[UFV] %s uses identity numbering: every residue maps to UniProt." % obj)
        return 0
    chains = set(_mapped_chains(obj))
    keys = {(g["chain"], g["resi"]) for g in _ca_geometry(obj) if g["uniPos"] is None and g["chain"] in chains}
    if not keys:
        print("[UFV] %s: every modelled residue of chain(s) %s maps to %s." % (obj, ",".join(sorted(chains)),
                                                                          mp.get("uid", "UniProt")))
        return 0
    idx = cmd.get_color_index(_color_name(CHIMERA_COLOR))
    with _batch():
        cmd.alter("(%s) and polymer" % obj, "color = _m.get((chain, resv), color)",
                  space={"_m": {k: idx for k in keys}})
        cmd.recolor(obj)
    _CARTOON_LAYER[obj] = "chimera"
    orgs = mp.get("organisms") or []
    print("[UFV] %s: %d residue(s) of chain(s) %s have no %s mapping (blue-grey)%s."
          % (obj, len(keys), ",".join(sorted({k[0] for k in keys})), mp.get("uid", "UniProt"),
             (" — chimeric: " + " / ".join(orgs)) if mp.get("chimeric") and orgs else ""))
    return len(keys)


def ufv_partners(obj=None):
    """ufv_partners [object] — other proteins in the loaded PDB entry (UniProt accession -> chains and
    mapped range, from SIFTS) plus chains assigned to another accession with ufv_chainmap."""
    obj = _resolve_object(obj)
    mp = UFV_MAPS.get(obj) or {}
    by = {}
    for p in mp.get("partners") or []:
        by.setdefault(p["accession"], []).append(p)
    custom = mp.get("partnersCustom") or {}
    if not by and not custom:
        print("[UFV] %s: no partner chains (other UniProt entries) known." % obj)
        return []
    print("[UFV] %s partner chains%s:" % (obj, (" (%d accessions)" % len(by)) if by else ""))
    for acc, ms in by.items():
        rng = sorted({"%s-%s" % (m["uniprotStart"], m["uniprotEnd"]) for m in ms if m.get("uniprotStart")})
        print("    %-12s chains %-10s %s%s" % (
            acc, ",".join(sorted({str(m["chainId"]) for m in ms})),
            ("UniProt " + ", ".join(rng)) if rng else "", ("  " + ms[0]["name"]) if ms[0].get("name") else ""))
    for ch, acc in sorted(custom.items()):
        print("    %-12s chain  %-10s (assigned with ufv_chainmap)" % (acc, ch))
    return [dict(p) for p in mp.get("partners") or []] + [{"accession": a, "chainId": c} for c, a in custom.items()]


@_undo_paused
def ufv_load(uid, name=None):
    """ufv_load uniprot_id [, name] — download the AlphaFold model, load it as a plain grey
    cartoon (identity numbering), and fetch the annotations. Nothing is projected automatically:
    choose layers from ufv_gui or the ufv_* commands so the view stays legible."""
    uid = uid.strip().upper()
    name = name or uid
    url = _alphafold_url(uid)
    if not url:
        print("[UFV] no AlphaFold model found for %s." % uid)
        return
    print("[UFV] downloading %s ..." % url)
    data = _get_structure_text(url)
    if not data:
        print("[UFV] failed to download model.")
        return
    fetch_annotations(uid)
    struct = {"key": name, "label": "AlphaFold %s" % uid, "source": "AlphaFold", "numbering": "identity", "url": url}
    _load_from_path(struct, _write_temp(data, "pdb"), uid, name=name)
    print("[UFV] loaded %s. Open the panel (ufv_gui) or add layers: ufv_ptms / "
          "ufv_variants / ufv_sites / ufv_domains / ufv_topology / ufv_alphamissense %s" % (name, name))


@_undo_paused
def ufv_info(uid=None):
    """ufv_info [uniprot_id] — print a summary of the cached annotations."""
    ann = fetch_annotations(_resolve_uid(uid))
    print("[UFV] %s  seq=%daa  PTMs=%d  variants=%d  sites=%d  domains=%d  topology=%d"
          % (ann["uid"], len(ann["sequence"]), len(ann["ptms"]), len(ann["variants"]),
             len(ann["sites"]), len(ann["domains"]), len(ann["topology"])))
    npath = sum(1 for v in ann["variants"] if "pathogenic" in v["consequence"].lower())
    print("[UFV] pathogenic/likely-pathogenic variants: %d" % npath)


# ----------------------------------------------------------------------------------------------
# Register commands
# ----------------------------------------------------------------------------------------------
if _HAS_PYMOL:
    for _name, _fn in [
        ("ufv_load", ufv_load), ("ufv_fetch", ufv_fetch), ("ufv_map", ufv_map),
        ("ufv_custom", ufv_custom), ("ufv_chainmap", ufv_chainmap), ("ufv_chimera", ufv_chimera),
        ("ufv_partners", ufv_partners),
        ("ufv_chain", ufv_chain), ("ufv_ptms", ufv_ptms), ("ufv_variants", ufv_variants),
        ("ufv_sites", ufv_sites), ("ufv_mutagenesis", ufv_mutagenesis),
        ("ufv_domains", ufv_domains), ("ufv_topology", ufv_topology),
        ("ufv_alphamissense", ufv_alphamissense), ("ufv_burden", ufv_burden),
        ("ufv_plddt", ufv_plddt), ("ufv_bfactor", ufv_bfactor),
        ("ufv_hotspots", ufv_hotspots), ("ufv_contacthubs", ufv_contacthubs), ("ufv_pockets", ufv_pockets),
        ("ufv_pocket_q", ufv_pocket_q), ("ufv_prox", ufv_prox), ("ufv_diseases", ufv_diseases),
        ("ufv_partner_layers", ufv_partner_layers),
        ("ufv_structures", ufv_structures), ("ufv_use", ufv_use),
        ("ufv_report", ufv_report), ("ufv_focus", ufv_focus), ("ufv_resetview", ufv_resetview),
        ("ufv_align", ufv_align), ("ufv_load_annotations", ufv_load_annotations),
        ("ufv_unload_annotations", ufv_unload_annotations),
        ("ufv_ligands", ufv_ligands), ("ufv_ligands_hide", ufv_ligands_hide),
        ("ufv_ligand_focus", ufv_ligand_focus), ("ufv_show_pocket", ufv_show_pocket),
        ("ufv_hide", ufv_hide), ("ufv_clear", ufv_clear), ("ufv_info", ufv_info),
        ("ufv_csv", ufv_csv),
        ("ufv_ligand_identity", ufv_ligand_identity), ("ufv_pockets_all", ufv_pockets_all),
        ("ufv_binding", ufv_binding), ("ufv_pdbekb", ufv_pdbekb), ("ufv_tractability", ufv_tractability),
        ("ufv_motif_search", ufv_motif_search), ("ufv_nearby", ufv_nearby),
        ("ufv_save_bfactor", ufv_save_bfactor), ("ufv_select", ufv_select),
        ("ufv_reset_annotations", ufv_reset_annotations), ("ufv_show_annotations", ufv_show_annotations),
    ]:
        cmd.extend(_name, _fn)


# ----------------------------------------------------------------------------------------------
# Trajectory: multi-state objects (MD trajectories, NMR ensembles)
# ----------------------------------------------------------------------------------------------
# Layers already follow the frame: every overlay object is a cmd.create copy of ALL states (the hidden
# Cα scaffold, the sphere layers carved from it, the ligand copies), cartoon colours are per atom, and
# focus / Nearby / report spatial queries run in the object's current state. Structure analyses
# (hotspots / hubs / constraint pocket / proximity / ligand contacts) are computed on the CURRENT state
# and cached per state (_traj_geom_state swaps the per-object caches when the frame changes).
#
# Trajectory analyses read coordinates on the main thread in ONE cmd.get_coords(sel, state=0) call per
# atom set (all states at once, ~30 ms for 500 frames x 2 k atoms; cached per object until its state or
# atom count changes); the numpy work (_traj_calc) is pure and thread-safe, so the panel runs it in a
# QThread. Shapes: X = (frames, atoms, 3). No per-frame or per-pair PyMOL objects are made: distances go
# to ONE measurement object (ufv_<obj>_dist, all states), RMSF / occupancy to ONE cartoon alter.
TRAJ_SETTINGS = {"start": 1, "stop": 0, "stride": 1, "ref": 0, "contact": 4.5, "dthr": 2.0}
_TRAJ_STASH_MAX = 32
_TRAJ_NSTATES = {}      # obj -> state count the overlay copies were last checked against
_TRAJ_FORCE = {}        # obj -> state forced by ufv_traj_tiers
_TRAJ_COORDS_MAX_BYTES = 384 * 1024 * 1024    # all objects together (float32 coordinate sets, re-read on demand)
TRAJ_KINDS = ("rmsf", "dist", "pockets", "ligand", "partner", "annot", "nearby")


def _obj_nstates(obj):
    try:
        return max(1, int(cmd.count_states(obj)))
    except Exception:
        return 1


def _obj_cur_state(obj):
    """The state `obj` currently shows (1 for single-state objects; a ufv_traj_tiers override wins)."""
    n = _obj_nstates(obj)
    if n <= 1:
        return 1
    f = _TRAJ_FORCE.get(obj)
    if f:
        return min(max(1, int(f)), n)
    try:
        s = int(cmd.get_object_state(obj))
    except Exception:
        s = int(cmd.get_state())
    return min(max(1, s), n)


def _state_sel(obj, sel, name="_ufv_cursel"):
    """`sel` evaluated in `obj`'s CURRENT state. PyMOL evaluates distance operators (around / within)
    over ALL states by default (an atom matches if it is near in ANY frame) — wrong and slow on a
    trajectory. Single-state objects get `sel` back unchanged (no extra call)."""
    if _obj_nstates(obj) <= 1:
        return sel
    cmd.select(name, sel, enable=0, state=_obj_cur_state(obj))
    return name


TRAJ_DEFER_ATOM_STATES = 500000    # atoms x states above which PyMOL's large-trajectory build mode is used
_TRAJ_DEFER_DONE = []


def _traj_defer_builds(obj, n):
    """With PyMOL's default defer_builds_mode 0, EVERY recolour of a multi-state object rebuilds the
    cartoon of ALL its states (measured: 2.6 s per colour change for 120 x 16 k atoms, ~11 s at 500
    states). Mode 3 builds only the state shown (scrubbing: ~70 ms/frame with cartoon + layers at 16 k
    atoms; recolour ~0.1 s), so it is switched on once for big trajectories unless the user set a mode."""
    if _TRAJ_DEFER_DONE or not TRAJ_SETTINGS.get("defer", True):
        return
    try:
        if n * cmd.count_atoms(obj) < TRAJ_DEFER_ATOM_STATES or int(float(cmd.get("defer_builds_mode"))) != 0:
            return
        with _no_undo():
            cmd.set("defer_builds_mode", 3)
        _TRAJ_DEFER_DONE.append(obj)
        print("[UFV] %s: %d states — defer_builds_mode set to 3 (only the shown state is built, so colouring "
              "stays fast; 'set defer_builds_mode, 0' to revert)." % (obj, n))
    except Exception:
        pass


def _traj_sync(obj, n=None):
    """When `obj` gained/lost states after its overlay copies were made (e.g. load_traj after mapping),
    rebuild the Cα scaffold, the sphere layers and the ligand layer so they follow every frame."""
    n = n or _obj_nstates(obj)
    if _TRAJ_NSTATES.get(obj) == n:
        return n
    _TRAJ_NSTATES[obj] = n
    if n > 1:
        _traj_defer_builds(obj, n)
    names = _object_set()
    sc = _ufv_layer_obj(obj, "ca")
    if sc not in names or _obj_nstates(sc) == n:
        return n
    with _no_undo():
        _TRAJ_COORDS.pop(obj, None)
        cmd.delete(sc)
        _LAYER_BUILT.pop(obj, None)
        for tag, groups in list((_SPHERE_STATE.get(obj) or {}).items()):
            lobj = _ufv_layer_obj(obj, tag)
            if lobj in names:
                cmd.delete(lobj)
            if tag == "prt":                 # partner layer: re-run ufv_partner_layers to rebuild
                _SPHERE_STATE[obj].pop(tag, None)
                continue
            try:
                _set_sphere_layer(obj, tag, groups)
            except Exception:
                pass
        lig_on = obj in _LIG_SHOWN
        for tag in ("lig", "ion"):
            if _ufv_layer_obj(obj, tag) in names:
                cmd.delete(_ufv_layer_obj(obj, tag))
        if lig_on:
            _LIG_SHOWN.discard(obj)
            try:
                ufv_ligands(obj)
            except Exception:
                pass
    print("[UFV] %s now has %d states: overlay layers rebuilt to follow every frame." % (obj, n))
    return n


def _traj_geom_state(obj, state=None):
    """State whose Cα geometry / analyses the caches hold for `obj`: the current one (or `state`). On a
    frame change the per-object caches are stashed and the new state's restored, so every analysis is
    computed on the frame shown and cached PER STATE. Single-state objects: always 1, no bookkeeping."""
    n = _obj_nstates(obj)
    if n <= 1:
        return 1
    _traj_sync(obj, n)
    st = min(max(1, int(state)), n) if state else _obj_cur_state(obj)
    prev = _TRAJ_GEOM_STATE.get(obj)
    if prev != st:
        stash = _TRAJ_STASH.setdefault(obj, {})
        keys = [k for k in list(_ANALYSIS_CACHE) if k[0] == obj]
        cur = (_GEOM_CACHE.pop(obj, None), {k: _ANALYSIS_CACHE.pop(k) for k in keys})
        if prev is not None:
            stash[prev] = cur
            while len(stash) > _TRAJ_STASH_MAX:
                stash.pop(next(iter(stash)))
        back = stash.pop(st, None)
        if back is not None:
            if back[0] is not None:
                _GEOM_CACHE[obj] = back[0]
            _ANALYSIS_CACHE.update(back[1])
        _TRAJ_GEOM_STATE[obj] = st
    return st


def _traj_frames(obj, start=None, stop=None, stride=None):
    """[state numbers] of the analysis window: TRAJ_SETTINGS start..stop (0 = last) every stride."""
    n = _obj_nstates(obj)
    s = int(start or TRAJ_SETTINGS["start"] or 1)
    e = int(stop or TRAJ_SETTINGS["stop"] or n)
    k = max(1, int(stride or TRAJ_SETTINGS["stride"] or 1))
    s = min(max(1, s), n)
    e = min(max(s, e), n)
    return list(range(s, e + 1, k))


def _traj_ref(obj, frames, ref=None):
    """Reference state (TRAJ_SETTINGS ref; 0 = first frame of the window)."""
    r = int(ref if ref not in (None, "") else TRAJ_SETTINGS["ref"] or 0)
    return min(max(1, r), _obj_nstates(obj)) if r else frames[0]


def _traj_coords(obj, sel):
    """(rows, X): rows = [(index, chain, resi, resv, resn, name)] of `sel` in `obj` (object atom order),
    X = float32 (nstates, n, 3) for EVERY state, read in one get_coords(state=0) call. Cached per object
    until its state or atom count changes. Main thread."""
    np = _numpy()
    n = _obj_nstates(obj)
    sig = (n, cmd.count_atoms(obj))
    ent = _TRAJ_COORDS.get(obj)
    if ent is None or ent["sig"] != sig:
        ent = _TRAJ_COORDS[obj] = {"sig": sig, "sets": {}}
    full = "(%s) and (%s)" % (obj, sel)
    rows = []
    cmd.iterate(full, "rows.append((index, chain, resi, resv, resn, name))", space={"rows": rows})
    key = (sel, hash(tuple(r[0] for r in rows)))      # named selections: same text, other atoms
    hit = ent["sets"].get(key)
    if hit is not None:
        return hit
    na = len(rows)
    X = cmd.get_coords(full, 0) if na else None
    if na and X is not None and X.shape[0] == n * na:
        X = X.reshape(n, na, 3).astype(np.float32, copy=False)
    else:          # atoms missing in some states (e.g. NMR models with different atoms): per state
        X = np.full((n, na, 3), np.nan, np.float32)
        for s in range(1, n + 1) if na else ():
            c = cmd.get_coords(full, s)
            if c is not None and len(c) == na:
                X[s - 1] = c
    ent["sets"][key] = hit = (rows, X)
    # one budget over ALL objects' coordinate sets (oldest first; the set just read always stays)
    pool = [(o, k) for o, e in list(_TRAJ_COORDS.items()) for k in list(e["sets"]) if (o, k) != (obj, key)]
    tot = sum(e["sets"][k][1].nbytes for e in _TRAJ_COORDS.values() for k in e["sets"])
    for o, k in pool:
        if tot <= _TRAJ_COORDS_MAX_BYTES:
            break
        tot -= _TRAJ_COORDS[o]["sets"].pop(k)[1].nbytes
    return hit


def _traj_ca(obj):
    """Deduplicated Cα (first per chain/resv, like _ca_geometry): (res, keep, X) with res =
    [(chain, resi, resv, resn, uniPos)] and keep = column indices into X."""
    rows, X = _traj_coords(obj, "polymer and name CA")
    keep, res, seen = [], [], set()
    for i, (_idx, ch, resi, resv, resn, _nm) in enumerate(rows):
        if (ch, resv) in seen:
            continue
        seen.add((ch, resv))
        keep.append(i)
        res.append((ch, resi, resv, resn, _resi_to_uni(obj, ch, resi)))
    return res, keep, X


def _traj_tags(obj, uid):
    """{uniPos: 'PTM,Path,Site,Hot,CP,PV'} annotation tags from cached data only (no network)."""
    ann = _CACHE.get(uid) or {}
    groups = _traj_groups(obj, uid, ann, with_lig=False)
    short = {"PTM sites": "PTM", "Pathogenic variants": "Path", "Functional sites": "Site", "Hotspots": "Hot",
             "Constraint pocket": "CP", "ProtVar pocket": "PV"}
    out = {}
    for g, ps in groups.items():
        for p in ps:
            out.setdefault(p, []).append(short.get(g, g))
    return {p: ",".join(v) for p, v in out.items()}


def _traj_groups(obj, uid, ann=None, with_lig=True):
    """Annotation residue groups {label: set(uniPos)} for the annotation-aware summaries: PTMs,
    pathogenic variants, functional sites, cached hotspots / constraint pocket, ProtVar pockets
    (best 3 of the crawl), ligand pocket (residues within 5 Å of a ligand in the current state)."""
    ann = ann if ann is not None else (_CACHE.get(uid) or {})
    g = {}
    g["PTM sites"] = {int(p["position"]) for p in ann.get("ptms") or [] if p.get("position")}
    g["Pathogenic variants"] = {int(v["position"]) for v in ann.get("variants") or []
                                if v.get("position") and _is_pathogenic(v)}
    s = set()
    for f in ann.get("sites") or []:
        a, b = int(f.get("position") or 0), int(f.get("endPosition") or f.get("position") or 0)
        if a and b - a <= 20:
            s.update(range(a, b + 1))
    g["Functional sites"] = s
    if uid:
        hs = _analysis_cached(obj, uid, "hotspots")
        if hs:
            g["Hotspots"] = {int(p) for p in hs}
        cp = _analysis_cached(obj, uid, "pockets")
        if cp:
            g["Constraint pocket"] = {int(p) for p in _pocket_classes(cp)}
        acc = _canon_acc(uid)
        pv = [v for k, v in _ALL_POCKETS.items() if k[0] == acc and v]
        if pv:
            g["ProtVar pocket"] = {int(r) for pk in pv[0]["pockets"][:3] for r in pk.get("resid") or []}
    if with_lig and _has_ligands(obj):
        try:
            nm = _state_sel(obj, "byres ((%s) and polymer and not hydro and ((%s) within 5 of ((%s) and not hydro)))"
                            % (obj, obj, _ligand_sel(obj)), "_ufv_tlig")
            near, _p = _near_from_selection(obj, nm)
            g["Ligand pocket"] = set(near)
        except Exception:
            pass
    return {k: v for k, v in g.items() if v}


def _ligand_ids(obj, min_heavy=1):
    """Ligand residues (non-polymer, non-solvent) -> [(resn, chain, resi, n_heavy)] largest first."""
    rows = []
    try:
        cmd.iterate("(%s) and not hydro" % _ligand_sel(obj), "rows.append((resn, chain, resi))", space={"rows": rows})
    except Exception:
        return []
    from collections import Counter
    c = Counter(rows)
    return sorted([(r, ch, ri, k) for (r, ch, ri), k in c.items() if k >= min_heavy], key=lambda x: -x[3])


def _lig_tok_sel(tok):
    """'lig:RESN[/chain[/resi]]' -> selection of that ligand's heavy atoms."""
    parts = tok[4:].split("/")
    sel = "not polymer and not solvent and not hydro and resn %s" % parts[0]
    if len(parts) > 1 and parts[1]:
        sel += " and chain %s" % parts[1]
    if len(parts) > 2 and parts[2]:
        sel += " and resi %s" % parts[2]
    return sel


def _res_tok_sel(obj, tok):
    """'120' / 'A/120' (UniProt position; raw resi when the object is unmapped) -> (selection, label)."""
    m = re.match(r"^\s*(?:([A-Za-z0-9]+)[/:])?(-?\d+)\s*$", str(tok))
    if not m:
        return None, None
    ch, p = m.group(1), int(m.group(2))
    if obj in UFV_MAPS:
        for c in ([ch] if ch else _mapped_chains(obj)):
            r = _uni_to_resi(obj, c, p)
            if r is not None:
                return ("chain %s and resi %s" % (c, _resi_tok(r)) if c else "resi %s" % _resi_tok(r)), tok
        return None, None
    return ("chain %s and resi %d" % (ch, p) if ch else "resi %d" % p), tok


# ---- prep (main thread: cmd + caches) -> calc (pure numpy, any thread) -> apply (main thread) ------
def _traj_prep(kind, obj=None, uid=None, frames=None, ref=None, spec=None, cutoff=None, threshold=None):
    """Main-thread half of a trajectory analysis: resolve the window, read the coordinate arrays and the
    annotation context. Returns a dict for _traj_calc (numpy arrays are copies: thread-safe)."""
    if _numpy() is None:
        raise RuntimeError("trajectory analyses need numpy")
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid, obj)
    n = _traj_sync(obj)
    frames = list(frames) if frames else _traj_frames(obj)
    r = _traj_ref(obj, frames, ref)
    ann = _CACHE.get(uid) or {}
    p = {"kind": kind, "obj": obj, "uid": uid, "nstates": n, "frames": frames, "ref": r, "t0": time.time(),
         "cutoff": float(cutoff or TRAJ_SETTINGS["contact"]),
         "threshold": float(threshold or TRAJ_SETTINGS["dthr"]), "spec": spec}
    res, keep, Xca = _traj_ca(obj)
    if not keep:
        raise RuntimeError("%s has no polymer Cα atoms" % obj)
    p["res"] = res
    p["CA"] = _traj_take(Xca, frames, keep)
    p["CAref"] = _traj_take(Xca, [r], keep)[0]
    p["tags"] = _traj_tags(obj, uid)
    if kind == "rmsf":
        p["groups"] = _traj_groups(obj, uid, ann)
    elif kind == "dist":
        _prep_dist(p, obj, uid, ann, spec)
    elif kind == "pockets":
        _prep_pockets(p, obj, uid, ann, spec)
    elif kind in ("ligand", "partner", "nearby"):
        _prep_contacts(p, kind, obj, spec)
    elif kind == "annot":
        for k in ("ligand", "partner"):
            last = (_TRAJ_RESULTS.get(obj) or {}).get(k)
            if last and last["frames"] == frames and last["ref"] == r and last["cutoff"] == p["cutoff"]:
                p[k] = {"reuse": last}           # same window: reuse the last contact run
                continue
            try:
                q = dict(p)
                _prep_contacts(q, k, obj, None)
                p[k] = q
            except Exception as e:
                p[k] = {"error": str(e)}
    p["prep_ms"] = (time.time() - p["t0"]) * 1000
    return p


def _traj_auto_pairs(obj, uid, ann, maxn=10):
    """Suggested pairs: each PTM with its nearest pathogenic-variant residue (Cα <= 15 Å, current
    state; closest 6), and each ligand (largest 2 with >= 5 heavy atoms) with up to 4 lining residues
    (annotated first)."""
    np = _numpy()
    geom = _ca_geometry(obj)
    xyz = {}
    for g in geom:
        if g["uniPos"] is not None and g["uniPos"] not in xyz:
            xyz[g["uniPos"]] = (g["x"], g["y"], g["z"])
    ptms = sorted({int(p["position"]) for p in ann.get("ptms") or [] if int(p.get("position") or 0) in xyz})
    path = sorted({int(v["position"]) for v in ann.get("variants") or []
                   if _is_pathogenic(v) and int(v.get("position") or 0) in xyz})
    out = []
    if ptms and path:
        P = np.array([xyz[q] for q in path])
        cand = []
        for t in ptms:
            d = np.sqrt(((P - np.array(xyz[t])) ** 2).sum(1))
            d[np.array(path) == t] = np.inf
            j = int(d.argmin())
            if d[j] <= 15.0:
                cand.append((float(d[j]), t, path[j]))
        seen = set()
        for d, a, b in sorted(cand):
            if (a, b) in seen:
                continue
            seen.add((a, b))
            out.append(("%d~%d" % (a, b), "PTM %d ↔ pathogenic %d" % (a, b)))
            if len(out) >= 6:
                break
    tags = _traj_tags(obj, uid)
    for resn, ch, ri, k in _ligand_ids(obj, 5)[:2]:
        lt = "lig:%s/%s/%s" % (resn, ch, ri)
        nm = _state_sel(obj, "byres ((%s) and polymer and not hydro and ((%s) within 4.5 of ((%s) and %s)))"
                        % (obj, obj, obj, _lig_tok_sel(lt)), "_ufv_tlin")
        near, _p = _near_from_selection(obj, nm)
        near.sort(key=lambda q: (q not in tags, q))
        for q in near[:4]:
            out.append(("%d~%s" % (q, lt), "%s ↔ lining %d" % (resn, q)))
    return out[:maxn]


def _prep_dist(p, obj, uid, ann, spec):
    spec = (spec or "auto").strip()
    toks = []
    for part in re.split(r"[;,\n]+", spec):
        part = part.strip()
        if not part:
            continue
        if part.lower() == "auto":
            toks.extend(_traj_auto_pairs(obj, uid, ann))
        else:
            toks.append((part, part))
    pairs, sels = [], []
    for tok, label in toks:
        if "~" in tok:
            a, b = tok.split("~", 1)
        else:
            m = re.match(r"^(.+?)\s*-\s*(.+)$", tok)
            if not m:
                print("[UFV] distance pair '%s': use A~B (e.g. 120~245, A/120~lig:ATP, 10+12+15~200)." % tok)
                continue
            a, b = m.groups()
        eps = []
        for e in (a.strip(), b.strip()):
            el = e.lower()
            if el.startswith("lig:"):
                eps.append(("lig", _lig_tok_sel(e), e))
            elif el.startswith("sel:"):
                eps.append(("cen", e[4:], e))
            elif "+" in e or e.startswith("{"):
                ps = re.findall(r"(?:[A-Za-z0-9]+/)?\d+", e)
                ss = [s for s in (_res_tok_sel(obj, x)[0] for x in ps) if s]
                eps.append(("cen", " or ".join("(%s and name CA)" % s for s in ss) if ss else None, e))
            else:
                s, _l = _res_tok_sel(obj, e)
                eps.append(("res", s and ("(%s) and polymer and not hydro" % s), e))
        if any(s is None for _k, s, _e in eps):
            print("[UFV] distance pair '%s': endpoint not found in %s." % (tok, obj))
            continue
        pairs.append({"label": label, "tok": tok, "eps": eps})
        sels.extend(s for _k, s, _e in eps)
    if not pairs:
        raise RuntimeError("no distance pairs resolved (spec: %s)" % spec)
    union = " or ".join("(%s)" % s for s in dict.fromkeys(sels))
    rows, X = _traj_coords(obj, union)
    np = _numpy()
    index = np.array([r[0] for r in rows])
    pos_of = {ix: i for i, ix in enumerate(index)}
    for pr in pairs:
        atoms = []
        for k, s, e in pr["eps"]:
            ids = []
            cmd.iterate("(%s) and (%s)" % (obj, s), "ids.append(index)", space={"ids": ids})
            c = [pos_of[i] for i in ids if i in pos_of]
            # residue: its Cα vs a residue / centroid; all its heavy atoms vs a ligand (closest contact)
            ca = [j for j in c if rows[j][5] == "CA"] if k == "res" else None
            atoms.append((k, c, ca))
        (ka, ca_, caa), (kb, cb_, cab) = atoms
        if ka == "res" and kb == "res" and caa and cab:
            mode, A, B = "ca", caa[:1], cab[:1]
        elif "cen" in (ka, kb):
            mode = "cen"
            A = caa[:1] if ka == "res" and caa else ca_
            B = cab[:1] if kb == "res" and cab else cb_
        else:
            mode, A, B = "min", ca_, cb_
        pr["mode"] = mode
        pr["A"] = _traj_take(X, p["frames"], A)
        pr["B"] = _traj_take(X, p["frames"], B)
        pr["idxA"] = [int(index[i]) for i in A]
        pr["idxB"] = [int(index[i]) for i in B]
        pr.pop("eps", None)
    p["pairs"] = pairs


def _prep_pockets(p, obj, uid, ann, spec):
    """Pocket residue sets -> Cα columns: 'auto' = constraint pocket (if computed), ProtVar pockets (best
    3, if crawled), ligand pockets (largest 4 ligands, residues within 5 Å in the reference state);
    'pos:1+2+3' / '1+2+3' = UniProt positions; 'pv:<id>' = a fetched ProtVar pocket; 'cp'; 'lig'."""
    col_of = {}
    for i, (ch, resi, resv, resn, u) in enumerate(p["res"]):
        if u is not None:
            col_of.setdefault(u, i)
    key_col = {(r[0], r[2]): i for i, r in enumerate(p["res"])}
    pockets = []
    spec = (spec or "auto").strip()
    for part in re.split(r"[;\n]+", spec):
        part = part.strip()
        pl = part.lower()
        if not part:
            continue
        if pl in ("auto", "cp"):
            cp = _analysis_cached(obj, uid, "pockets") if uid else None
            cls = _pocket_classes(cp) if cp else {}
            pos = sorted(u for u, c in cls.items() if c == "pocket") or sorted(cls)
            for i, cl in enumerate(_traj_clusters(p["CAref"], [col_of[u] for u in pos if u in col_of])[:3]):
                pockets.append(("Constraint pocket %d" % (i + 1), [("k", p["res"][c][0], p["res"][c][2]) for c in cl]))
        if pl in ("auto", "pv"):
            acc = _canon_acc(uid) if uid else None
            pv = [v for k, v in _ALL_POCKETS.items() if k[0] == acc and v]
            for pk in (pv[0]["pockets"][:3] if pv else []):
                pockets.append(("ProtVar pocket %s" % pk.get("pocketId"), [int(r) for r in pk.get("resid") or []]))
        if pl.startswith("pv:"):
            pk = _find_pocket(uid, part[3:].strip()) if uid else None
            if pk:
                pockets.append(("ProtVar pocket %s" % pk.get("pocketId"), [int(r) for r in pk.get("resid") or []]))
        if pl in ("auto", "lig") and _has_ligands(obj):
            for resn, ch, ri, k in _ligand_ids(obj, 5)[:4]:
                lt = "lig:%s/%s/%s" % (resn, ch, ri)
                cmd.select("_ufv_tpk", "byres ((%s) and polymer and ((%s) within 5 of ((%s) and %s)))"
                           % (obj, obj, obj, _lig_tok_sel(lt)), enable=0, state=p["ref"])
                keys = []
                cmd.iterate("_ufv_tpk and name CA", "keys.append((chain, resv))", space={"keys": keys})
                pockets.append(("%s %s/%s pocket" % (resn, ch, ri), [("k",) + k_ for k_ in keys]))
        if re.match(r"^(pos:)?[\d+ ]+$", pl) and "+" in pl:
            pockets.append(("Residues " + part.replace("pos:", ""), [int(x) for x in re.findall(r"\d+", part)]))
    out = []
    for label, items in pockets:
        cols = sorted({key_col.get((it[1], it[2])) if isinstance(it, tuple) else col_of.get(it) for it in items} - {None})
        if len(cols) >= 3:
            out.append({"label": label, "cols": cols})
    if not out:
        raise RuntimeError("no pocket with >= 3 modelled residues (run Constraint pocket / All pockets, load "
                           "a ligand, or give positions like 120+124+156)")
    p["pockets"] = out


def _prep_contacts(p, kind, obj, spec):
    """Atom sets for contact occupancy. ligand: protein heavy atoms of residues near a ligand vs the
    ligands' heavy atoms; partner: interface heavy atoms of every polymer chain; nearby: one residue's
    heavy atoms vs the residues/ligands around it. Candidates come from ONE PyMOL search in the reference
    state within cutoff + pad, pad = 2 x the 95th-percentile Cα displacement (after fitting) + 6 Å
    (<= 25 Å): a pair that ends up in contact in some frame has moved at most that much relative to the
    reference (an any-state search costs seconds at 500 frames). _contact_pairs is then exact."""
    np = _numpy()
    c = p["cutoff"]
    tr = _not_transplant(obj)
    if "pad" not in p:
        Y, _T = _fit(p["CA"], p["CAref"])
        ri = _ref_index(p)
        Yr = Y[ri] if ri is not None else p["CAref"]
        dmax = float(np.percentile(np.sqrt(((Y - Yr) ** 2).sum(-1)).max(0), 95)) if len(Y) > 1 else 0.0
        p["pad"] = min(25.0, 2.0 * dmax + 6.0)
    pad, st = p["pad"], p["ref"]
    if kind == "ligand":
        if not _has_ligands(obj):
            raise RuntimeError("%s has no ligands" % obj)
        lsel = "(%s) and not hydro" % _ligand_sel(obj)
        cmd.select("_ufv_tct", "byres ((%s) and polymer%s and not hydro and ((%s) within %g of (%s)))"
                   % (obj, tr, obj, c + pad, lsel), enable=0, state=st)
        A, B = "_ufv_tct and not hydro", "not polymer and not solvent and not hydro"
    elif kind == "partner":
        chains = [ch for ch in (cmd.get_chains("(%s) and polymer" % obj) or [])]
        if "" in chains and len(chains) > 1:
            raise RuntimeError("%s has polymer atoms without a chain ID (inter-chain contacts need chain IDs)" % obj)
        if len(chains) < 2:
            raise RuntimeError("%s has a single polymer chain (no partner / inter-chain contacts)" % obj)
        ours = [x for x in (_mapped_chains(obj) if obj in UFV_MAPS else []) if x in chains]
        sub = bool(ours) and not set(ours) >= set(chains)
        # our protein's chains vs every other chain (contacts among partner chains are not ours to report):
        # interface atoms on both sides of each of OUR chains only
        parts = []
        for ch in (ours if sub else chains):
            own, other = "(%s) and polymer and chain %s" % (obj, ch), "(%s) and polymer and not chain %s" % (obj, ch)
            parts.append("((%s) and ((%s) within %g of (%s)))" % (own, obj, c + pad, other))
            if sub:
                parts.append("((%s) and ((%s) within %g of (%s)))" % (other, obj, c + pad, own))
        cmd.select("_ufv_tct", "byres (%s)" % " or ".join(parts), enable=0, state=st)
        B = "_ufv_tct and not hydro"
        A = "%s and chain %s" % (B, "+".join(ours)) if sub else B
    else:   # nearby: spec = UniProt position (default: the focused residue)
        pos = spec if spec not in (None, "") else (_FOCUS_SHOWN.get(obj) or {}).get("pos")
        s, _l = _res_tok_sel(obj, str(pos)) if pos is not None else (None, None)
        if not s:
            raise RuntimeError("nearby: give a residue (UniProt position) or focus one first")
        d = float(NEARBY_SETTINGS["dist"])
        p["cutoff"] = c = d
        p["center"] = str(pos)
        cmd.select("_ufv_tcr", "(%s) and (%s)" % (obj, s), enable=0)
        cmd.select("_ufv_tct", "byres ((%s) and not hydro and not solvent and not _ufv_tcr%s and ((%s) within %g of "
                   "(_ufv_tcr and not hydro)))" % (obj, tr, obj, d + pad), enable=0, state=st)
        A, B = "_ufv_tcr and not hydro", "_ufv_tct and not hydro"
    rowsA, XA = _traj_coords(obj, A)
    rowsB, XB = (rowsA, XA) if B == A else _traj_coords(obj, B)
    if not rowsA or not rowsB:
        raise RuntimeError("no %s contacts possible in %s" % (kind, obj))
    p["same"] = B == A
    p["A"] = _traj_take(XA, p["frames"], f32=True)
    p["B"] = p["A"] if p["same"] else _traj_take(XB, p["frames"], f32=True)
    p["Aref"] = _traj_take(XA, [p["ref"]], f32=True)[0]
    p["Bref"] = p["Aref"] if p["same"] else _traj_take(XB, [p["ref"]], f32=True)[0]
    p["rowsA"] = [r[1:] for r in rowsA]
    p["rowsB"] = [r[1:] for r in rowsB]
    p["uniA"] = [_resi_to_uni(obj, r[1], r[2]) for r in rowsA]
    p["uniB"] = [_resi_to_uni(obj, r[1], r[2]) for r in rowsB]
    p["ours"] = set(_mapped_chains(obj)) if obj in UFV_MAPS else None


# ---- apply (main thread) ---------------------------------------------------------------------------
_RAMP_BWR = ["#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"]
_RAMP_OCC = ["#fff5eb", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04"]


def _traj_color_keys(obj, keyval, ramp, lo, hi, base="#c8c8c8"):
    """Colour residues {(chain, resv): value} on a ramp in ONE alter (others -> base)."""
    idx, m = {}, {}
    span = (hi - lo) or 1.0
    for k, v in keyval.items():
        b = int(round((len(ramp) - 1) * min(1.0, max(0.0, (v - lo) / span))))
        h = ramp[b]
        if h not in idx:
            idx[h] = _color_index(_color_name(h))
        m[k] = idx[h]
    with _batch():
        cmd.alter("(%s) and polymer" % obj, "color = _m.get((chain, resv), _b)",
                  space={"_m": m, "_b": _color_index(_color_name(base))})
        cmd.recolor(obj)
    return len(m)


def _traj_apply(r, color=True, show=True):
    """Main-thread half: store the result, colour the cartoon (RMSF / occupancy), draw distances."""
    np = _numpy()
    obj = r["obj"]
    _res = _TRAJ_RESULTS.setdefault(obj, {})
    _res.pop(r["kind"], None)             # re-insert: dict order = recency (the panel shows the last one)
    _res[r["kind"]] = r
    with _no_undo():
        if r["kind"] == "rmsf" and color:
            keys = {(c, rv): float(v) for (c, _ri, rv, _rn, _u), v in zip(r["res"], r["rmsf"])}
            _set_cartoon_layer(obj, "rmsf", grey=False)       # _traj_color_keys colours every polymer atom
            _traj_color_keys(obj, keys, _RAMP_BWR, 0.0, float(np.percentile(r["rmsf"], 98)) or 1.0)
        elif r["kind"] in ("ligand", "partner") and color:
            _set_cartoon_layer(obj, "%s occupancy" % r["kind"], grey=False)
            _traj_color_keys(obj, r["occupancy"], _RAMP_OCC, 0.0, 1.0)
        elif r["kind"] == "dist" and show:
            # PyMOL measures every state of a distance (~0.3-0.9 ms per state per pair), so only the
            # pairs that change most are drawn by default; any other pair: click its row / _traj_show_pairs
            rng = [float(s.max() - s.min()) for _n, s in r["series"]]
            order = sorted(range(len(rng)), key=lambda i: -rng[i])
            flagged = [i for i in order if r["rowmeta"][i]["flag"]]
            _traj_show_pairs(r, (flagged or order)[:int(TRAJ_SETTINGS.get("show_pairs", 3))])
    return r


def _traj_show_pairs(r, which):
    """ONE measurement object ufv_<obj>_dist for the chosen pairs, measured in EVERY state (PyMOL keeps
    one distance per state, so it follows frame scrubbing with no per-frame objects)."""
    obj = r["obj"]
    name = _ufv_layer_obj(obj, "dist")
    with _batch():
        cmd.delete(name)
        for i in which[:12]:
            pr = r["pairs"][i]
            if pr.get("disp"):
                a, b = pr["disp"]
                cmd.distance(name, "(%s) and index %d" % (obj, a), "(%s) and index %d" % (obj, b),
                             mode=0, cutoff=999, state=0, zoom=0)
            else:
                sa = "(%s) and index %s" % (obj, "+".join(str(x) for x in pr["idxA"]))
                sb = "(%s) and index %s" % (obj, "+".join(str(x) for x in pr["idxB"]))
                cmd.distance(name, sa, sb, mode=4, state=0, zoom=0)
        if name in _object_set():
            cmd.color("orange", name)
            cmd.set("dash_gap", 0.25, name)
            cmd.set("label_size", 12, name)
    return name


def _traj_run(kind, obj=None, uid=None, spec=None, color=True, show=True, **kw):
    p = _traj_prep(kind, obj, uid, spec=spec, **kw)
    r = _traj_calc(p)
    return _traj_apply(r, color, show)


def _traj_shift(obj, val):
    """`cmd [object,] value`: a first argument that is not a loaded object is the value."""
    if obj not in (None, "") and str(obj) not in (cmd.get_names("objects") or []):
        return None, obj
    return obj, val


def _traj_cmd(kind, obj, uid, spec=None, **kw):
    try:
        r = _traj_run(kind, obj, uid, spec=spec, **kw)
    except RuntimeError as e:
        print("[UFV] %s: %s" % (kind, e))
        return None
    _traj_print(r)
    return r


# ---- commands ------------------------------------------------------------------------------------
@_undo_paused
def ufv_traj(obj=None):
    """ufv_traj [object] — trajectory info: states, current state, analysis window and results."""
    obj = _resolve_object(obj)
    n = _traj_sync(obj)
    fr = _traj_frames(obj)
    print("[UFV] %s: %d state(s), showing %d; window %d-%d step %d (%d frames), ref %d; results: %s"
          % (obj, n, _obj_cur_state(obj), fr[0], fr[-1], TRAJ_SETTINGS["stride"], len(fr), _traj_ref(obj, fr),
             ", ".join(sorted(_TRAJ_RESULTS.get(obj) or {})) or "none"))
    return n


def ufv_traj_frames(start=None, stop=None, stride=None, ref=None):
    """ufv_traj_frames [start [, stop [, stride [, ref]]]] — analysis window (stop 0 = last state;
    ref 0 = first frame of the window). Without arguments prints the window."""
    for k, v in (("start", start), ("stop", stop), ("stride", stride), ("ref", ref)):
        if v not in (None, ""):
            TRAJ_SETTINGS[k] = max(0, int(float(v)))
    TRAJ_SETTINGS["stride"] = max(1, TRAJ_SETTINGS["stride"])
    TRAJ_SETTINGS["start"] = max(1, TRAJ_SETTINGS["start"])
    print("[UFV] trajectory window: start %(start)d, stop %(stop)d (0 = last), stride %(stride)d, ref %(ref)d"
          % TRAJ_SETTINGS)
    return dict(TRAJ_SETTINGS)


@_undo_paused
def ufv_traj_rmsf(obj=None, uid=None, color=1):
    """ufv_traj_rmsf [object [, uid [, color]]] — Cα RMSD (fit to the reference frame) per frame and
    per-residue RMSF (refit to the average); colours the cartoon blue (rigid) -> red (mobile) and
    compares mean RMSF of PTM sites / pathogenic variants / sites / pockets / hotspots with the rest."""
    return _traj_cmd("rmsf", obj, uid, color=str(color) not in ("0", "false", "off"))


@_undo_paused
def ufv_traj_dist(obj=None, pairs="auto", threshold=None, show=1, uid=None):
    """ufv_traj_dist [object [, pairs [, threshold [, show]]]] — distance time series. pairs separated by
    ';': 120~245 (Cα–Cα, UniProt positions; A/120 picks a chain), 120~lig:ATP[/chain/resi] (closest heavy
    atoms), 10+12+15~200 or sel:<selection> (centroids), auto = PTM ↔ nearest pathogenic variant + ligand
    ↔ lining residues. Pairs whose distance ranges more than `threshold` Å (default 2) are flagged and
    drawn as ONE measurement object (ufv_<obj>_dist) that follows the frame."""
    obj, pairs = _traj_shift(obj, pairs)
    return _traj_cmd("dist", obj, uid, spec=pairs, threshold=threshold, show=str(show) not in ("0", "false", "off"))


@_undo_paused
def ufv_traj_pockets(obj=None, pockets="auto", uid=None):
    """ufv_traj_pockets [object [, pockets]] — pocket dynamics per frame: Rg and mean pairwise distance
    (size) of the lining Cα, buriedness (Cα within 10 Å), opening/closing events (size beyond median ±
    2 robust sd) and lining persistence. pockets: auto (constraint pocket + ProtVar pockets + ligand
    pockets that are available) | cp | pv | pv:<id> | lig | 120+124+156 ; several separated by ';'."""
    obj, pockets = _traj_shift(obj, pockets)
    return _traj_cmd("pockets", obj, uid, spec=pockets)


@_undo_paused
def ufv_traj_contacts(obj=None, kind="ligand", cutoff=None, color=1, uid=None):
    """ufv_traj_contacts [object [, ligand|partner [, cutoff [, color]]]] — contact occupancy (heavy atoms
    within cutoff, default 4.5 Å) across the window: residue–ligand or inter-chain (partner); cartoon
    coloured by occupancy (white -> dark orange), annotated residues tagged in the table."""
    if obj in ("ligand", "partner") and kind in (None, "", "ligand", "partner"):
        obj, kind = None, obj
    return _traj_cmd(kind or "ligand", obj, uid, cutoff=cutoff, color=str(color) not in ("0", "false", "off"))


@_undo_paused
def ufv_traj_annot(obj=None, uid=None):
    """ufv_traj_annot [object] — annotation dynamics: for every annotated residue (PTM, pathogenic
    variant, site, hotspot, pocket) the fraction of frames it contacts a ligand / another chain and the
    number of on/off switches."""
    return _traj_cmd("annot", obj, uid)


@_undo_paused
def ufv_traj_nearby(obj=None, pos=None, uid=None):
    """ufv_traj_nearby [object [, position]] — Nearby occupancy: residues/ligands with any heavy atom
    within the Nearby distance (ufv_nearby) of the residue, as a fraction of frames (default residue =
    the focused one)."""
    obj, pos = _traj_shift(obj, pos)
    return _traj_cmd("nearby", obj, uid, spec=pos)


@_undo_paused
def ufv_traj_tiers(obj=None, kind="hotspots", state=None, uid=None):
    """ufv_traj_tiers [object [, hotspots|hubs|pockets [, state]]] — run a structure analysis on one
    state (default: the current one). Results are cached per state, so revisiting a frame is instant."""
    obj = _resolve_object(obj)
    st = int(state) if state not in (None, "") else _obj_cur_state(obj)
    fn = {"hotspots": ufv_hotspots, "hubs": ufv_contacthubs, "contacthubs": ufv_contacthubs,
          "pockets": ufv_pockets}.get((kind or "hotspots").lower())
    if fn is None:
        print("[UFV] ufv_traj_tiers: kind = hotspots | hubs | pockets")
        return None
    _TRAJ_FORCE[obj] = st
    try:
        fn(obj, uid)
    finally:
        _TRAJ_FORCE.pop(obj, None)
    print("[UFV] (computed on state %d of %s)" % (st, obj))
    return st


@_undo_paused
def ufv_traj_csv(prefix=None, obj=None):
    """ufv_traj_csv [prefix [, object]] — write <prefix>_frames.csv (one row per state: RMSD, Rg, every
    distance pair, pocket metrics, contact counts of the last runs) and <prefix>_residues.csv (chain,
    resi, UniProt, RMSF, ligand / partner occupancy, pocket persistence, annotation tags)."""
    import csv
    obj = _resolve_object(obj)
    res = _TRAJ_RESULTS.get(obj) or {}
    if not res:
        print("[UFV] no trajectory results for %s yet." % obj)
        return None
    prefix = prefix or os.path.join(os.getcwd(), "%s_traj" % re.sub(r"[^A-Za-z0-9_.-]", "_", obj))
    cols, data, frames = [], {}, set()
    for kind in TRAJ_KINDS:
        r = res.get(kind)
        if not r:
            continue
        for name, arr in r.get("series") or []:
            c = "%s:%s" % (kind, name)
            cols.append(c)
            for f, v in zip(r["frames"], arr):
                data[(f, c)] = float(v)
                frames.add(f)
    f1 = prefix + "_frames.csv"
    with open(f1, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["state"] + cols)
        for f in sorted(frames):
            w.writerow([f] + [("%.4f" % data[(f, c)]) if (f, c) in data else "" for c in cols])
    info, rmsf, pers = {}, {}, {}
    rm = res.get("rmsf")
    if rm:
        for (c, ri, rv, rn, u), v in zip(rm["res"], rm["rmsf"]):
            rmsf[(c, rv)] = float(v)
            info[(c, rv)] = (c, ri, rn, u)
    occ = {k: (res.get(k) or {}).get("occupancy") or {} for k in ("ligand", "partner")}
    pk = res.get("pockets")
    if pk:
        for lab, d in pk["persist"].items():
            for col, f in d.items():
                c, ri, rv, rn, u = pk["res"][col]
                info.setdefault((c, rv), (c, ri, rn, u))
                pers.setdefault((c, rv), []).append("%s=%.2f" % (lab, f))
    tags = _traj_tags(obj, _resolve_uid(None, obj))
    for src in (occ["ligand"], occ["partner"]):
        for k in src:
            if k not in info:
                info[k] = (k[0], str(k[1]), "", _resi_to_uni(obj, k[0], k[1]))
    f2 = prefix + "_residues.csv"
    with open(f2, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "resi", "resn", "uniprot_pos", "rmsf", "ligand_occupancy", "partner_occupancy",
                    "pocket_persistence", "tags"])
        for k in sorted(info, key=lambda k: (k[0], k[1])):
            c, ri, rn, u = info[k]
            w.writerow([c, ri, rn, "" if u is None else u,
                        ("%.3f" % rmsf[k]) if k in rmsf else "",
                        ("%.3f" % occ["ligand"][k]) if k in occ["ligand"] else "",
                        ("%.3f" % occ["partner"][k]) if k in occ["partner"] else "",
                        ";".join(pers.get(k, [])), tags.get(u, "") if u is not None else ""])
    print("[UFV] wrote %s and %s." % (f1, f2))
    return f1, f2


if _HAS_PYMOL:
    for _name, _fn in [("ufv_traj", ufv_traj), ("ufv_traj_frames", ufv_traj_frames), ("ufv_traj_rmsf", ufv_traj_rmsf),
                       ("ufv_traj_dist", ufv_traj_dist), ("ufv_traj_pockets", ufv_traj_pockets),
                       ("ufv_traj_contacts", ufv_traj_contacts), ("ufv_traj_annot", ufv_traj_annot),
                       ("ufv_traj_nearby", ufv_traj_nearby), ("ufv_traj_tiers", ufv_traj_tiers),
                       ("ufv_traj_csv", ufv_traj_csv)]:
        cmd.extend(_name, _fn)


# ---- Qt: TrajectoryPanel (self-contained; embedded by the main panel via _ufv_make_traj_panel) --------
TrajectoryPanel = None
_TRAJ_LABELS = [("RMSD / RMSF", "rmsf", ""), ("Distances", "dist", "auto  |  120~245; A/45~lig:ATP; 10+12~200"),
                ("Pocket dynamics", "pockets", "auto  |  cp; pv; lig; 120+124+156"),
                ("Ligand contacts", "ligand", ""), ("Partner contacts", "partner", ""),
                ("Annotation dynamics", "annot", ""), ("Nearby occupancy", "nearby", "position (default: focused)")]


_ELIDED_CLS = None      # module-level elided label (the trajectory panel is built outside the main panel)
_WHEEL_FILTER = None    # the one wheel-guard event filter shared by every guarded widget


def _ufv_elided_label_class():
    """Single-line QLabel that elides to its width ('…'); text() keeps the full text, shown on hover.
    Its size hint does not depend on the text, so new text never re-lays out (or widens) the panel."""
    global _ELIDED_CLS
    if _ELIDED_CLS is None:
        from pymol.Qt import QtWidgets, QtCore, QtGui

        class _UfvElided(QtWidgets.QLabel):
            def __init__(self, text=""):
                super(_UfvElided, self).__init__(text)
                self.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Fixed)
                self.setToolTip(text or "")

            def setText(self, text):
                super(_UfvElided, self).setText(text)
                self.setToolTip(text or "")
                self.update()

            def minimumSizeHint(self):
                return QtCore.QSize(10, self.fontMetrics().height())

            def sizeHint(self):
                return QtCore.QSize(10, self.fontMetrics().height())

            def paintEvent(self, _ev):
                p = QtGui.QPainter(self)
                r = self.contentsRect()
                txt = self.fontMetrics().elidedText((self.text() or "").replace("\n", " · "), QtCore.Qt.ElideRight,
                                                    r.width())
                self.style().drawItemText(p, r, int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter), self.palette(),
                                          self.isEnabled(), txt, self.foregroundRole())
                p.end()
        _ELIDED_CLS = _UfvElided
    return _ELIDED_CLS


def _ufv_wheel_guard(root):
    """Mouse wheel over a combo box / spin box / slider anywhere under `root` never changes its value:
    the wheel scrolls the enclosing scroll area instead (the panel, the numbering table), so scrolling
    the panel with the pointer passing over a control is safe. An OPEN combo popup still scrolls its
    own list (the popup is a separate window, not filtered). Controls get Qt.StrongFocus (no focus
    by wheel). Idempotent: call again after adding widgets."""
    global _WHEEL_FILTER
    from pymol.Qt import QtWidgets, QtCore
    if _WHEEL_FILTER is None:

        class _UfvWheelGuard(QtCore.QObject):
            def eventFilter(self, w, ev):
                if ev.type() != QtCore.QEvent.Wheel:
                    return False
                try:
                    if isinstance(w, QtWidgets.QComboBox) and w.view().isVisible():
                        return False
                    sa = w.parentWidget()
                    while sa is not None and not isinstance(sa, QtWidgets.QAbstractScrollArea):
                        sa = sa.parentWidget()
                    if sa is not None:
                        sb = sa.verticalScrollBar()
                        dy = ev.angleDelta().y()
                        if dy and sb.maximum() > sb.minimum():
                            lines = QtWidgets.QApplication.wheelScrollLines()
                            sb.setValue(sb.value() - int(round(dy / 120.0 * lines * sb.singleStep())))
                except Exception:
                    pass
                ev.accept()
                return True
        _WHEEL_FILTER = _UfvWheelGuard()
    if root is None:
        return 0
    n = 0
    for cls in (QtWidgets.QComboBox, QtWidgets.QAbstractSpinBox, QtWidgets.QSlider):
        ws = list(root.findChildren(cls))
        if isinstance(root, cls):
            ws.append(root)
        for w in ws:
            if w.property("ufvWheelGuard"):
                continue
            w.setFocusPolicy(QtCore.Qt.StrongFocus)
            w.installEventFilter(_WHEEL_FILTER)
            w.setProperty("ufvWheelGuard", True)
            n += 1
    return n


def _build_traj_panel_class():
    global TrajectoryPanel
    if TrajectoryPanel is not None:
        return TrajectoryPanel
    from pymol.Qt import QtWidgets, QtCore, QtGui
    np = _numpy()

    class _TrajWorker(QtCore.QThread):
        done = QtCore.Signal(object)

        def __init__(self, prep):
            super(_TrajWorker, self).__init__()
            self.prep = prep

        def run(self):
            try:
                self.done.emit(_traj_calc(self.prep))
            except Exception as exc:
                self.done.emit(exc)

    class _TrajPlot(QtWidgets.QWidget):
        """Compact QPainter line plot: series vs state, cursor at the shown state, click/drag = jump."""
        frameClicked = QtCore.Signal(int)
        COLORS = ["#00897b", "#d55e00", "#0072b2", "#cc79a7", "#e6a817", "#6a1b9a", "#009e73", "#56b4e9"]

        def __init__(self, parent=None):
            super(_TrajPlot, self).__init__(parent)
            self.setFixedHeight(120)        # fixed: results never change the panel's height
            self.frames, self.series, self.cursor, self.hl, self._last = [], [], None, None, None

        def set_data(self, frames, series, hl=None):
            self.frames = list(frames or [])
            self.series = [(n, np.asarray(a, float)) for n, a in (series or [])]
            self.hl = hl
            self.update()

        def set_cursor(self, f):
            if f != self.cursor:
                self.cursor = f
                self.update()

        def _area(self):
            return self.rect().adjusted(34, 4, -4, -13)

        def _shown(self):
            if self.hl is not None and 0 <= self.hl < len(self.series):
                return [self.series[self.hl]], [self.hl]
            idx = list(range(min(8, len(self.series))))
            return [self.series[i] for i in idx], idx

        def paintEvent(self, ev):
            p = QtGui.QPainter(self)
            p.setRenderHint(QtGui.QPainter.Antialiasing, True)
            p.fillRect(self.rect(), QtGui.QColor("#ffffff"))
            a = self._area()
            p.setPen(QtGui.QColor("#d4d9de"))
            p.drawRect(a)
            f = p.font()
            f.setPointSizeF(7.5)
            p.setFont(f)
            if not self.frames or not self.series:
                p.setPen(QtGui.QColor("#8a929b"))
                p.drawText(a, QtCore.Qt.AlignCenter, "Run an analysis — click the plot to jump to a frame")
                return
            shown, idx = self._shown()
            ys = np.concatenate([s for _n, s in shown]) if shown else np.zeros(1)
            lo, hi = float(np.nanmin(ys)), float(np.nanmax(ys))
            if hi - lo < 1e-9:
                lo, hi = lo - 0.5, hi + 0.5
            f0, f1 = self.frames[0], self.frames[-1]
            fx = lambda fr: a.left() + (fr - f0) / float(max(1, f1 - f0)) * a.width()
            fy = lambda y: a.bottom() - (y - lo) / (hi - lo) * a.height()
            p.setPen(QtGui.QColor("#5a626c"))
            p.drawText(QtCore.QRectF(0, a.top() - 2, 32, 12), QtCore.Qt.AlignRight, "%.3g" % hi)
            p.drawText(QtCore.QRectF(0, a.bottom() - 10, 32, 12), QtCore.Qt.AlignRight, "%.3g" % lo)
            p.drawText(QtCore.QRectF(a.left(), a.bottom() + 1, 60, 12), QtCore.Qt.AlignLeft, str(f0))
            p.drawText(QtCore.QRectF(a.right() - 60, a.bottom() + 1, 60, 12), QtCore.Qt.AlignRight, str(f1))
            fr = np.asarray(self.frames, float)
            step = max(1, int(len(fr) // max(1, a.width())))
            for (name, s), i in zip(shown, idx):
                pen = QtGui.QPen(QtGui.QColor(self.COLORS[i % len(self.COLORS)]))
                pen.setWidthF(1.6 if self.hl is not None else 1.1)
                p.setPen(pen)
                p.drawPolyline(QtGui.QPolygonF([QtCore.QPointF(fx(x), fy(y)) for x, y in zip(fr[::step], s[::step])]))
            for j, ((name, _s), i) in enumerate(list(zip(shown, idx))[:4]):
                p.setPen(QtGui.QColor(self.COLORS[i % len(self.COLORS)]))
                p.drawText(QtCore.QRectF(a.left() + 4, a.top() + 1 + 10 * j, a.width() - 8, 11),
                           QtCore.Qt.AlignLeft, name[:60])
            if self.cursor is not None and f0 <= self.cursor <= f1:
                p.setPen(QtGui.QPen(QtGui.QColor("#e53935"), 1.2))
                x = fx(self.cursor)
                p.drawLine(QtCore.QPointF(x, a.top()), QtCore.QPointF(x, a.bottom()))
                p.drawText(QtCore.QRectF(x + 2, a.bottom() - 11, 50, 11), QtCore.Qt.AlignLeft, str(self.cursor))

        def _emit(self, ev):
            if not self.frames:
                return
            x = ev.position().x() if hasattr(ev, "position") else ev.pos().x()
            a = self._area()
            t = min(1.0, max(0.0, (x - a.left()) / float(max(1, a.width()))))
            target = self.frames[0] + t * (self.frames[-1] - self.frames[0])
            f = min(self.frames, key=lambda fr: abs(fr - target))
            if f != self._last:
                self._last = f
                self.frameClicked.emit(int(f))

        def mousePressEvent(self, ev):
            self._last = None
            self._emit(ev)

        def mouseMoveEvent(self, ev):
            self._emit(ev)

    class _TrajectoryPanel(QtWidgets.QWidget):
        """Compact trajectory tools: window, analysis picker, plot with frame cursor, result table.
        `panel` = the main UFV panel (obj() / cur_uid() give the current object / accession)."""

        def __init__(self, panel=None):
            super(_TrajectoryPanel, self).__init__(panel if isinstance(panel, QtWidgets.QWidget) else None)
            self.main = panel
            self.result = None
            self._workers = set()
            self._obj = None
            self._tick_n = 0
            self.setStyleSheet("QWidget { font-size: 11px; } QPushButton, QToolButton { padding: 1px 6px; } "
                               "QSpinBox, QComboBox, QLineEdit { padding: 0px 2px; } "
                               "QTableWidget { font-size: 10px; }")
            v = QtWidgets.QVBoxLayout(self)
            v.setContentsMargins(0, 0, 0, 0)
            v.setSpacing(3)
            self.head = QtWidgets.QLabel("Trajectory: needs a multi-state object")
            self.head.setStyleSheet("color: #5a626c;")
            v.addWidget(self.head)
            self.body = QtWidgets.QWidget()
            b = QtWidgets.QVBoxLayout(self.body)
            b.setContentsMargins(0, 0, 0, 0)
            b.setSpacing(3)
            v.addWidget(self.body)
            r1 = QtWidgets.QHBoxLayout()
            r1.setSpacing(3)
            self.sp = {}
            for key, lbl, tip in (("start", "", "First state"), ("stop", "–", "Last state (end = last state)"),
                                  ("stride", "step", "Use every n-th state"),
                                  ("ref", "ref", "Reference state (first = the first analysed state)")):
                if lbl:
                    r1.addWidget(QtWidgets.QLabel(lbl))
                s = QtWidgets.QSpinBox()
                s.setRange(0 if key in ("stop", "ref") else 1, 10 ** 7)
                if key in ("stop", "ref"):
                    s.setSpecialValueText("end" if key == "stop" else "first")   # 0 = auto
                s.setValue(int(TRAJ_SETTINGS[key]))
                s.setToolTip(tip)
                s.setMaximumWidth(62)
                s.valueChanged.connect(lambda val, k=key: TRAJ_SETTINGS.__setitem__(k, int(val)))
                self.sp[key] = s
                r1.addWidget(s)
            r1.addStretch(1)
            b.addLayout(r1)
            r2 = QtWidgets.QHBoxLayout()
            r2.setSpacing(3)
            self.kind = QtWidgets.QComboBox()
            for lbl, k, _ph in _TRAJ_LABELS:
                self.kind.addItem(lbl, k)
            self.kind.currentIndexChanged.connect(self._kind_changed)
            r2.addWidget(self.kind)
            self.param = QtWidgets.QLineEdit()
            r2.addWidget(self.param, 1)
            self.run_btn = QtWidgets.QPushButton("Run")
            self.run_btn.setObjectName("primary")
            self.run_btn.clicked.connect(self.run)
            r2.addWidget(self.run_btn)
            self.csv_btn = QtWidgets.QToolButton()
            self.csv_btn.setText("CSV")
            self.csv_btn.setToolTip("Per-frame and per-residue CSV of all trajectory results")
            self.csv_btn.clicked.connect(self.save_csv)
            r2.addWidget(self.csv_btn)
            self.bf_btn = QtWidgets.QToolButton()
            self.bf_btn.setText("B-PDB")
            self.bf_btn.setToolTip("Save a PDB with the last RMSF / occupancy result in the B-factor column")
            self.bf_btn.clicked.connect(self.save_bfactor)
            r2.addWidget(self.bf_btn)
            b.addLayout(r2)
            self.plot = _TrajPlot()
            self.plot.frameClicked.connect(self.goto)
            b.addWidget(self.plot)
            self.table = QtWidgets.QTableWidget(0, 0)
            self.table.verticalHeader().setVisible(False)
            self.table.verticalHeader().setDefaultSectionSize(17)
            self.table.setFixedHeight(5 * 17 + 26)  # header + 5 rows; more rows scroll inside
            self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.table.cellClicked.connect(self.on_row)
            b.addWidget(self.table)
            self.status = _ufv_elided_label_class()()  # one line, full text on hover: never re-flows the panel
            self.status.setStyleSheet("color: #5a626c;")
            b.addWidget(self.status)
            self._kind_changed()
            _ufv_wheel_guard(self)
            self.timer = QtCore.QTimer(self)
            self.timer.timeout.connect(self._tick)
            self.timer.start(200)
            self.refresh()

        # ---- context ----
        def obj(self):
            try:
                o = self.main.obj() if self.main is not None and hasattr(self.main, "obj") else None
            except Exception:
                o = None
            return o or _resolve_object(None)

        def uid(self):
            try:
                if self.main is not None and hasattr(self.main, "cur_uid"):
                    return self.main.cur_uid()
            except Exception:
                pass
            return None

        def refresh(self):
            o = self.obj()
            n = _obj_nstates(o) if o and o in (cmd.get_names("objects") or []) else 1
            self._obj = o
            multi = n > 1
            self.body.setVisible(multi)
            self.head.setVisible(not multi)       # the section header already shows "<obj> · N states"
            self.head.setText("Trajectory tools: the current object has a single state")
            for k in ("start", "stop", "ref"):
                self.sp[k].setMaximum(max(1, n))
            # an analysis run from the command line (ufv_traj_*) is shown too: latest result of this object
            latest = list((_TRAJ_RESULTS.get(o) or {}).values())[-1:] if multi else []
            if latest and latest[0] is not self.result and latest[0].get("table"):
                self.show_result(latest[0])
            return n

        # The 200 ms timer (frame cursor + 1 s state-count refresh) only runs while the panel is on screen:
        # a hidden section (single-state object, collapsed section, closed main panel) polls nothing.
        def hideEvent(self, ev):
            self.timer.stop()
            super(_TrajectoryPanel, self).hideEvent(ev)

        def showEvent(self, ev):
            if not self.timer.isActive():
                self.timer.start(200)
                self.refresh()
            super(_TrajectoryPanel, self).showEvent(ev)

        def _tick(self):
            self._tick_n += 1
            if self._tick_n % 5 == 0 or self.obj() != self._obj:
                self.refresh()
            if self.body.isVisible():
                try:
                    self.plot.set_cursor(int(cmd.get_state()))
                except Exception:
                    pass

        def _kind_changed(self, *_a):
            ph = _TRAJ_LABELS[self.kind.currentIndex()][2]
            self.param.setPlaceholderText(ph)
            self.param.setEnabled(bool(ph))

        # ---- actions ----
        def goto(self, f):
            with _no_undo():
                cmd.frame(int(f))
            self.plot.set_cursor(int(f))

        def run(self):
            k = self.kind.currentData()
            o = self.obj()
            spec = self.param.text().strip() or None
            if k == "nearby" and spec is None:
                spec = (_FOCUS_SHOWN.get(o) or {}).get("pos")
            self.status.setText("Reading coordinates ...")
            try:
                with _no_undo():
                    prep = _traj_prep(k, o, self.uid(), spec=spec)
            except Exception as e:
                self.status.setText("Error: %s" % e)
                return
            self.run_btn.setEnabled(False)
            self.status.setText("Computing %s over %d frames ..." % (self.kind.currentText(), len(prep["frames"])))
            wk = _TrajWorker(prep)
            self._workers.add(wk)
            wk.finished.connect(lambda: self._workers.discard(wk))

            def done(res):
                self.run_btn.setEnabled(True)
                if isinstance(res, Exception):
                    self.status.setText("Error: %s" % res)
                    return
                with _no_undo():
                    _traj_apply(res)
                self.show_result(res)
            wk.done.connect(done)
            wk.start()

        def show_result(self, r):
            self.result = r
            self.plot.set_data(r["frames"], r.get("series"))
            hd, rows = r["table"]
            self.table.clear()
            self.table.setColumnCount(len(hd))
            self.table.setRowCount(len(rows))
            self.table.setHorizontalHeaderLabels([str(h) for h in hd])
            for i, row in enumerate(rows):
                for j, x in enumerate(row):
                    self.table.setItem(i, j, QtWidgets.QTableWidgetItem(str(x)))
            self.table.resizeColumnsToContents()
            self.status.setText("%s  [read %.0f ms, compute %.0f ms]" % (r["summary"], r["prep_ms"], r["calc_ms"]))

        def on_row(self, row, _col):
            r = self.result
            if not r or row >= len(r.get("rowmeta") or []):
                return
            m = r["rowmeta"][row]
            o = r["obj"]
            with _no_undo():
                if "pair" in m:
                    _traj_show_pairs(r, [m["pair"]])
                    self.plot.set_data(r["frames"], r["series"], hl=m["pair"])
                elif "pocket" in m:
                    i = [x["pocket"] for x in r["rowmeta"]].index(m["pocket"])
                    self.plot.set_data(r["frames"], r["series"][3 * i:3 * i + 3])
                elif m.get("pos") is not None:
                    try:
                        ufv_focus(o, m["pos"], r.get("uid"))
                    except Exception as e:
                        self.status.setText("focus: %s" % e)
                elif m.get("key"):
                    ch, rv = m["key"]
                    _focus_zoom("(%s) and chain %s and resi %s" % (o, ch, rv), 6)

        def save_csv(self):
            o = self.obj()
            fn, _f = QtWidgets.QFileDialog.getSaveFileName(self, "Trajectory CSV prefix", "%s_traj" % o, "CSV (*.csv)")
            if fn:
                out = ufv_traj_csv(re.sub(r"\.csv$", "", fn), o)
                self.status.setText("Wrote %s" % ", ".join(out or []))

        def save_bfactor(self):
            r = self.result
            mode = {"rmsf": "rmsf", "ligand": "ligocc", "partner": "partnerocc"}.get(r["kind"]) if r else None
            if not mode:
                self.status.setText("Run RMSF or contacts first.")
                return
            o = r["obj"]
            fn, _f = QtWidgets.QFileDialog.getSaveFileName(self, "B-factor PDB", "%s_%s.pdb" % (o, mode), "PDB (*.pdb)")
            if fn:
                ufv_save_bfactor(fn, o, mode, r.get("uid"))
                self.status.setText("Wrote %s" % fn)

    TrajectoryPanel = _TrajectoryPanel
    return TrajectoryPanel


def _ufv_make_traj_panel(panel):
    """Hook for the main panel: a TrajectoryPanel bound to it (placed in its _traj_slot)."""
    return _build_traj_panel_class()(panel)


# ----------------------------------------------------------------------------------------------
# Qt control panel (PyMOL is Qt-based) + Plugin Manager hook
# ----------------------------------------------------------------------------------------------
_gui_ref = None  # keep a reference so the window isn't garbage-collected


def ufv_gui():
    """ufv_gui - open the Qt control panel (visualization driver + report generator)."""
    global _gui_ref
    try:
        from pymol.Qt import QtWidgets, QtCore, QtGui
    except Exception as e:
        print("[UFV] Qt GUI unavailable (%s). Use the ufv_* commands instead." % e)
        return
    cls = _build_panel_class(QtWidgets, QtCore, QtGui)
    _gui_ref = cls()
    _gui_ref.show()
    return _gui_ref


_UFV_QSS = """
* { font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif; font-size: 12px; color: #2c3038; }
QScrollArea, QScrollArea > QWidget > QWidget { background: #f6f7f9; }
#ufvHeader { background: #00695c; }
#ufvTitle { color: #ffffff; font-size: 13px; font-weight: 700; }
#ufvSub { color: #9ed7d0; font-size: 10px; }
#ufvSecHead { background: #eaf1f0; border-top: 1px solid #d5dfdd; border-bottom: 1px solid #e3e9e8; }
#ufvSecHead QToolButton { background: transparent; border: none; padding: 2px 2px; color: #00796b;
                          font-weight: 700; font-size: 11px; }
#ufvSecHead QToolButton:focus { color: #004d40; text-decoration: underline; }
#ufvSecSum { color: #6b747d; font-size: 10px; }
#ufvSecBody { background: #ffffff; }
#ufvFooter { background: #f0f3f4; border-top: 1px solid #d5dbe0; }
QToolButton[chip="true"] { background: #f0f2f4; border: 1px solid #d4d9de; border-radius: 9px;
                           padding: 1px 6px; font-size: 11px; color: #4a525b; }
QToolButton[chip="true"]:hover { border-color: #00897b; }
QToolButton[chip="true"]:checked { background: #00897b; border-color: #00897b; color: #ffffff; }
QToolButton[chip="true"]:disabled { color: #b4bac1; }
QGroupBox {
    background: #ffffff; border: 1px solid #e4e7eb; border-radius: 9px;
    margin-top: 15px; padding: 10px 9px 9px 9px; font-weight: 600;
}
QGroupBox::title {
    subcontrol-origin: margin; subcontrol-position: top left; left: 11px; padding: 1px 6px;
    color: #00897b; font-weight: 700; font-size: 10px;
}
QPushButton {
    background: #f0f2f4; border: 1px solid #d4d9de; border-radius: 5px; padding: 2px 8px;
}
QPushButton:hover { background: #e7ebef; border-color: #b6bdc5; }
QPushButton:pressed { background: #dbe0e5; }
QPushButton:disabled { color: #aab0b8; background: #f4f5f7; }
QPushButton#primary { background: #00897b; color: #ffffff; border: none; font-weight: 600; }
QPushButton#primary:hover { background: #00796b; }
QPushButton#primary:pressed { background: #00695c; }
QToolButton { background: #f0f2f4; border: 1px solid #d4d9de; border-radius: 5px; padding: 2px 6px; }
QToolButton:hover { background: #e7ebef; }
QToolButton:checked { background: #00897b; color: #ffffff; border-color: #00897b; }
QLineEdit, QComboBox {
    border: 1px solid #d4d9de; border-radius: 5px; padding: 1px 6px; background: #ffffff;
    selection-background-color: #00897b; selection-color: #ffffff;
}
QLineEdit:focus, QComboBox:focus { border: 1px solid #00897b; }
QComboBox::drop-down { border: none; width: 18px; }
QComboBox QAbstractItemView { border: 1px solid #d4d9de; background: #ffffff; selection-background-color: #e0f2f1; selection-color: #2c3038; }
QComboBox QAbstractItemView::item:disabled { color: #b4bac1; }   /* e.g. pLDDT on experimental, B-factor on predicted */
QComboBox:disabled, QLineEdit:disabled, QCheckBox:disabled, QLabel:disabled { color: #b4bac1; }
QCheckBox { spacing: 6px; }
QTableWidget, QTableView {
    border: 1px solid #e4e7eb; border-radius: 6px; gridline-color: #f0f1f3; background: #ffffff;
    selection-background-color: #e0f2f1; selection-color: #2c3038;
}
QHeaderView::section {
    background: #f3f5f7; border: none; border-bottom: 1px solid #e4e7eb; padding: 2px 5px;
    font-weight: 600; color: #5a626c;
}
QTextBrowser { border: 1px solid #e4e7eb; border-radius: 6px; background: #ffffff; padding: 4px; }
QListWidget { border: 1px solid #e4e7eb; border-radius: 6px; background: #ffffff; }
QScrollBar:vertical { width: 11px; background: transparent; margin: 0; }
QScrollBar::handle:vertical { background: #cbd2d9; border-radius: 5px; min-height: 26px; }
QScrollBar::handle:vertical:hover { background: #aeb7c0; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
"""


def _build_panel_class(QtWidgets, QtCore, QtGui):
    from collections import Counter

    class _Worker(QtCore.QThread):
        done = QtCore.Signal(object)

        def __init__(self, fn):
            super(_Worker, self).__init__()
            self._fn = fn

        def run(self):
            try:
                self.done.emit(self._fn())
            except Exception as exc:
                self.done.emit(exc)

    class _CsvWorker(QtCore.QThread):
        """Builds the residue CSV off the UI thread, optionally bulk-fetching ProtVar with progress."""
        done = QtCore.Signal(object)
        progress = QtCore.Signal(int, int)

        def __init__(self, ann, obj, with_pv, geom=None, chains=None, lig_contacts=None):
            super(_CsvWorker, self).__init__()
            self.ann, self.obj, self.with_pv = ann, obj, with_pv
            self.geom, self.chains, self.lig_contacts = geom, chains, lig_contacts

        def run(self):
            try:
                pv = None
                if self.with_pv:
                    pv = _fetch_protvar_all(self.ann["uid"], len(self.ann["sequence"]), self.progress.emit)
                self.done.emit(_residue_csv(self.ann, self.obj, pv, self.geom, self.chains, self.lig_contacts))
            except Exception as exc:
                self.done.emit(exc)

    class _PocketWorker(QtCore.QThread):
        """ProtVar all-pocket crawl (fetch_all_pockets) off the UI thread, with progress."""
        done = QtCore.Signal(object)
        progress = QtCore.Signal(int, int)

        def __init__(self, uid, length):
            super(_PocketWorker, self).__init__()
            self.uid, self.length = uid, length

        def run(self):
            try:
                self.done.emit(fetch_all_pockets(self.uid, self.length, self.progress.emit))
            except Exception as exc:
                self.done.emit(exc)

    # ---- compact-layout helpers -------------------------------------------------------------
    def _ufv_settings():
        """Persisted panel state (collapsed sections, window size). Offscreen (test) runs get the defaults
        unless UFV_SETTINGS_INI names an .ini file to use instead."""
        ini = os.environ.get("UFV_SETTINGS_INI")
        if ini:
            return QtCore.QSettings(ini, QtCore.QSettings.IniFormat)
        if os.environ.get("QT_QPA_PLATFORM") == "offscreen":
            return None
        try:
            return QtCore.QSettings("UFV", "UniProt3DFeatureViewer-PyMOL")
        except Exception:
            return None

    class _FlowLayout(QtWidgets.QLayout):
        """Left-to-right layout that wraps onto further lines (chips, tool rows): the panel never
        needs a horizontal scroll bar, whatever the dock width."""

        def __init__(self, parent=None, hspace=4, vspace=3):
            super(_FlowLayout, self).__init__(parent)
            self._items = []
            self._hs, self._vs = hspace, vspace
            self.setContentsMargins(0, 0, 0, 0)

        def addItem(self, item):
            self._items.append(item)

        def count(self):
            return len(self._items)

        def itemAt(self, i):
            return self._items[i] if 0 <= i < len(self._items) else None

        def takeAt(self, i):
            return self._items.pop(i) if 0 <= i < len(self._items) else None

        def expandingDirections(self):
            return QtCore.Qt.Orientations(0)

        def hasHeightForWidth(self):
            return True

        def heightForWidth(self, w):
            return self._place(QtCore.QRect(0, 0, w, 0), True)

        def setGeometry(self, rect):
            super(_FlowLayout, self).setGeometry(rect)
            self._place(rect, False)

        def sizeHint(self):
            w, h = 0, 0
            for it in self._items:
                if not it.isEmpty():
                    s = it.sizeHint(); w += s.width() + self._hs; h = max(h, s.height())
            return QtCore.QSize(max(0, w - self._hs), h)

        def minimumSize(self):
            s = QtCore.QSize(0, 0)
            for it in self._items:
                if not it.isEmpty():
                    s = s.expandedTo(it.minimumSize())
            return s

        def _place(self, rect, test):
            x, y, line_h, right = rect.x(), rect.y(), 0, rect.x() + rect.width()
            for it in self._items:
                if it.isEmpty():
                    continue
                sz = it.sizeHint()
                if x > rect.x() and x + sz.width() > right:
                    x = rect.x(); y += line_h + self._vs; line_h = 0
                if not test:
                    it.setGeometry(QtCore.QRect(QtCore.QPoint(x, y), sz))
                x += sz.width() + self._hs; line_h = max(line_h, sz.height())
            return y + line_h - rect.y()

    class _ElidedLabel(QtWidgets.QLabel):
        """Single-line label that elides to the available width ('…'). text() stays the full text;
        the tooltip shows the full text plus any extra tooltip set on it."""

        def __init__(self, text=""):
            super(_ElidedLabel, self).__init__(text)
            self._extra_tip = ""
            self.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Fixed)
            self._sync_tip()

        def setText(self, text):
            super(_ElidedLabel, self).setText(text)
            self._sync_tip(); self.update()

        def setToolTip(self, tip):
            self._extra_tip = tip or ""
            self._sync_tip()

        def _sync_tip(self):
            full = self.text() or ""
            QtWidgets.QLabel.setToolTip(self, "\n\n".join(t for t in (full, self._extra_tip) if t))

        def minimumSizeHint(self):
            return QtCore.QSize(10, self.fontMetrics().height())

        def sizeHint(self):
            fm = self.fontMetrics()
            return QtCore.QSize(fm.boundingRect(self.text() or " ").width() + 4, fm.height())

        def paintEvent(self, _ev):
            p = QtGui.QPainter(self)
            r = self.contentsRect()
            txt = self.fontMetrics().elidedText((self.text() or "").replace("\n", " · "), QtCore.Qt.ElideRight, r.width())
            self.style().drawItemText(p, r, int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter), self.palette(),
                                      self.isEnabled(), txt, self.foregroundRole())
            p.end()

    class _Section(QtWidgets.QWidget):
        """Collapsible panel section: slim header (▸/▾ title · summary) over a body. The collapsed
        state is remembered per section (QSettings)."""

        def __init__(self, key, title, tip="", settings=None, expanded=True):
            super(_Section, self).__init__()
            self.key, self._settings = key, settings
            v = QtWidgets.QVBoxLayout(self); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(0)
            self.head = QtWidgets.QFrame(); self.head.setObjectName("ufvSecHead")
            hl = QtWidgets.QHBoxLayout(self.head); hl.setContentsMargins(4, 0, 8, 0); hl.setSpacing(6)
            self.btn = QtWidgets.QToolButton(); self.btn.setText(title.replace("&", "&&"))   # no mnemonic
            self.btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
            self.btn.setToolTip((tip + "\n\n" if tip else "") + "Click (or Space) to collapse / expand")
            self.btn.clicked.connect(lambda _=False: self.setExpanded(not self.isExpanded()))
            hl.addWidget(self.btn)
            self.summary = _ElidedLabel(""); self.summary.setObjectName("ufvSecSum")
            hl.addWidget(self.summary, 1)
            v.addWidget(self.head)
            self.body = QtWidgets.QWidget(); self.body.setObjectName("ufvSecBody")
            self.body.setAttribute(QtCore.Qt.WA_StyledBackground, True)
            self.head.setAttribute(QtCore.Qt.WA_StyledBackground, True)
            self.body_lay = QtWidgets.QVBoxLayout(self.body)
            self.body_lay.setContentsMargins(8, 4, 8, 6); self.body_lay.setSpacing(4)
            v.addWidget(self.body, 1)
            if settings is not None:
                val = settings.value("sections/" + key, None)
                if val is not None:
                    expanded = val in (True, 1, "true", "1")
            self._expanded = None
            self.setExpanded(expanded, save=False)

        def mousePressEvent(self, ev):      # clicking anywhere on the header bar toggles too
            if self.head.geometry().contains(ev.pos()):
                self.setExpanded(not self.isExpanded())
            super(_Section, self).mousePressEvent(ev)

        def isExpanded(self):
            return bool(self._expanded)

        def setExpanded(self, on, save=True):
            on = bool(on)
            if on == self._expanded:
                return
            self._expanded = on
            self.body.setVisible(on)
            self.btn.setArrowType(QtCore.Qt.DownArrow if on else QtCore.Qt.RightArrow)
            # a collapsed section must not soak up the stretch given to it (the report section)
            self.setMaximumHeight(16777215 if on else self.head.sizeHint().height())
            if save and self._settings is not None:
                try:
                    self._settings.setValue("sections/" + self.key, on)
                except Exception:
                    pass

        def setSummary(self, text):
            self.summary.setText(text or "")

    class _ReportBrowser(QtWidgets.QTextBrowser):
        """Residue report view whose preferred height is its minimum, so it takes exactly the height
        left over by the other sections (QTextBrowser's default 192 px hint would force scrolling)."""

        def sizeHint(self):
            return QtCore.QSize(256, max(self.minimumHeight(), 60))

        def loadResource(self, rtype, url):
            """ufvspider:<px>:<h,a,ac,b,p> -> the pocket-character radar (same plot as the extension)."""
            if url.scheme() == "ufvspider":
                try:
                    px, vals = url.toString()[len("ufvspider:"):].split(":", 1)
                    return _spider_image(int(px), [float(v) for v in vals.split(",")])
                except Exception:
                    return None
            return super(_ReportBrowser, self).loadResource(rtype, url)

    def _spider_image(size, fracs):
        """Radar over Hydrophobic / Aromatic / Acidic / Basic / Polar (axes are independent fractions)."""
        k = 2  # drawn at 2x for sharp edges on the 1x-sized <img>
        img = QtGui.QImage(size * k, size * k, QtGui.QImage.Format_ARGB32)
        img.fill(QtCore.Qt.transparent)
        p = QtGui.QPainter(img); p.setRenderHint(QtGui.QPainter.Antialiasing)
        p.scale(k, k)
        cx = cy = size / 2.0; R = size * 48.0 / 168; n = len(fracs)
        pt = lambda i, r: QtCore.QPointF(cx + r * math.cos(-math.pi / 2 + i * 2 * math.pi / n),
                                         cy + r * math.sin(-math.pi / 2 + i * 2 * math.pi / n))
        grid = QtGui.QPen(QtGui.QColor("#cfd6dd")); grid.setWidthF(1)
        p.setPen(grid); p.setBrush(QtCore.Qt.NoBrush)
        for f in (0.5, 1.0):
            p.drawPolygon(QtGui.QPolygonF([pt(i, R * f) for i in range(n)]))
        font = p.font(); font.setPixelSize(max(8, int(size * 9 / 168))); font.setBold(True); p.setFont(font)
        for i in range(n):
            p.setPen(grid); p.drawLine(QtCore.QPointF(cx, cy), pt(i, R))
            p.setPen(QtGui.QColor("#5f6b76"))
            c = pt(i, R + size * 15.0 / 168)
            p.drawText(QtCore.QRectF(c.x() - 30, c.y() - 7, 60, 14), QtCore.Qt.AlignCenter, SPIDER_LABELS[i])
        poly = QtGui.QPolygonF([pt(i, R * max(0.0, min(1.0, f))) for i, f in enumerate(fracs)])
        pen = QtGui.QPen(QtGui.QColor("#00bfa5")); pen.setWidthF(2.2)
        p.setPen(pen); p.setBrush(QtGui.QColor(0, 137, 123, 77)); p.drawPolygon(poly)
        p.setBrush(QtGui.QColor("#00bfa5"))
        for q in poly:
            p.drawEllipse(q, 2.4, 2.4)
        p.end()
        return img

    def _keep_space(*widgets):
        """Hidden but still laid out: showing / hiding the widget never moves anything else."""
        for w in widgets:
            sp = w.sizePolicy(); sp.setRetainSizeWhenHidden(True); w.setSizePolicy(sp)

    def _chip(text, tip="", checked=False, color=None):
        """Small checkable pill button (layer toggles, filters). Same API as a QCheckBox for the slots:
        isChecked / setChecked / clicked(bool) / toggled(bool)."""
        b = QtWidgets.QToolButton(); b.setText(text); b.setCheckable(True); b.setChecked(checked)
        b.setProperty("chip", True); b.setFocusPolicy(QtCore.Qt.StrongFocus)
        if tip:
            b.setToolTip(tip)
        if color:
            b.setStyleSheet("QToolButton { color: #8a929b; padding: 1px 4px; font-size: 10px; } QToolButton:checked { color: %s;"
                            " background: #ffffff; border: 1px solid %s; }" % (color, color))
        return b

    # ---- Annotations v3 list: model (all rows of a category) -> proxy (search filter) -> view --------------
    _CHECKED, _UNCHECKED, _PARTIAL = QtCore.Qt.Checked, QtCore.Qt.Unchecked, QtCore.Qt.PartiallyChecked

    class _AnnModel(QtCore.QAbstractTableModel):
        """Rows of the current category (dicts from ufv_core.ann_rows). Column 0 = check box (checked = the
        row's id is in the panel's checked set of its category), 1 = category (All only), 2 = position,
        3 = annotation. No per-row Qt items: 4,000 rows cost one list."""

        def __init__(self, panel):
            super(_AnnModel, self).__init__(panel)
            self.rows, self.on, self.unmod, self.head = [], {}, frozenset(), "Annotation"
            self._brushes = {}

        def set_rows(self, rows, on, unmod):
            self.beginResetModel()
            self.rows, self.on, self.unmod = rows, on, unmod
            self.endResetModel()

        def rowCount(self, parent=QtCore.QModelIndex()):
            return 0 if parent.isValid() else len(self.rows)

        def columnCount(self, parent=QtCore.QModelIndex()):
            return 0 if parent.isValid() else 4

        def is_on(self, r):
            return r["id"] in (self.on.get(r["cat"]) or ())

        def _brush(self, color):
            b = self._brushes.get(color)
            if b is None:
                b = self._brushes[color] = QtGui.QBrush(QtGui.QColor(color or "#5a626c"))
            return b

        def data(self, idx, role=QtCore.Qt.DisplayRole):
            if not idx.isValid():
                return None
            r, c = self.rows[idx.row()], idx.column()
            if role == QtCore.Qt.DisplayRole:
                return (None, ANN_CAT_SHORT.get(r["cat"], r["cat"]), r["poslabel"], r["label"])[c]
            if role == QtCore.Qt.CheckStateRole and c == 0:
                return _CHECKED if self.is_on(r) else _UNCHECKED
            if role == QtCore.Qt.ForegroundRole and c:
                return self._brush("#b4bac1" if r["id"] in self.unmod else r["color"])
            if role == QtCore.Qt.ToolTipRole:
                return r["label"] + ("  (not modelled in this structure)" if r["id"] in self.unmod else "")
            return None

        def headerData(self, sec, orient, role=QtCore.Qt.DisplayRole):
            if orient == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
                return ("", "Cat", "Pos", self.head)[sec]
            if orient == QtCore.Qt.Horizontal and role == QtCore.Qt.ToolTipRole and sec == 0:
                return "Tick / untick every listed (searched) row"
            return None

        def flags(self, idx):
            # the check box is painted from CheckStateRole and toggled by _AnnView (whole column 0 is the
            # hit area; a click on a row of a multi-row selection toggles the whole selection)
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

        def refresh_checks(self):
            if self.rows:
                self.dataChanged.emit(self.index(0, 0), self.index(len(self.rows) - 1, 0))

    class _AnnProxy(QtCore.QSortFilterProxyModel):
        """Search filter: rows whose id is in `ok` (None = every row). The ids come from ONE ann_search call."""

        def __init__(self, parent=None):
            super(_AnnProxy, self).__init__(parent)
            self.ok = None

        def set_ok(self, ok):
            self.ok = ok
            self.invalidateFilter()

        def filterAcceptsRow(self, r, _parent):
            ok = self.ok
            return ok is None or self.sourceModel().rows[r]["id"] in ok

    class _CheckHeader(QtWidgets.QHeaderView):
        """Horizontal header whose first section is a tri-state check box (all / some / none of the listed rows)."""
        toggled = QtCore.Signal()

        def __init__(self, parent=None):
            super(_CheckHeader, self).__init__(QtCore.Qt.Horizontal, parent)
            self.state = _UNCHECKED
            self.setSectionsClickable(True)

        def set_state(self, st):
            if st != self.state:
                self.state = st
                self.viewport().update()

        def paintSection(self, painter, rect, idx):
            painter.save()
            super(_CheckHeader, self).paintSection(painter, rect, idx)
            painter.restore()
            if idx != 0:
                return
            opt = QtWidgets.QStyleOptionButton()
            s = self.style().pixelMetric(QtWidgets.QStyle.PM_IndicatorWidth, None, self)
            opt.rect = QtCore.QRect(rect.x() + (rect.width() - s) // 2, rect.y() + (rect.height() - s) // 2, s, s)
            opt.state = QtWidgets.QStyle.State_Enabled | (
                QtWidgets.QStyle.State_On if self.state == _CHECKED else
                QtWidgets.QStyle.State_NoChange if self.state == _PARTIAL else QtWidgets.QStyle.State_Off)
            self.style().drawPrimitive(QtWidgets.QStyle.PE_IndicatorCheckBox, opt, painter, self)

        def mousePressEvent(self, ev):
            if self.logicalIndexAt(ev.pos()) == 0:
                self.toggled.emit(); ev.accept(); return
            super(_CheckHeader, self).mousePressEvent(ev)

    class _AnnView(QtWidgets.QTableView):
        """Space = toggle the selected rows; Enter = zoom; a press in column 0 toggles that row (or the whole
        selection when the row is part of a multi-row selection) without changing the selection."""
        space = QtCore.Signal()
        enter = QtCore.Signal()
        toggle_at = QtCore.Signal(object)

        def keyPressEvent(self, ev):
            k = ev.key()
            if k == QtCore.Qt.Key_Space:
                self.space.emit(); ev.accept(); return
            if k in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
                self.enter.emit(); ev.accept(); return
            super(_AnnView, self).keyPressEvent(ev)

        def _col0(self, ev):
            idx = self.indexAt(ev.pos())
            return idx if idx.isValid() and idx.column() == 0 else None

        def mousePressEvent(self, ev):
            idx = self._col0(ev) if ev.button() == QtCore.Qt.LeftButton else None
            if idx is not None:
                self.toggle_at.emit(idx); ev.accept(); return
            super(_AnnView, self).mousePressEvent(ev)

        def mouseDoubleClickEvent(self, ev):
            idx = self._col0(ev) if ev.button() == QtCore.Qt.LeftButton else None
            if idx is not None:          # a fast second click on a check box toggles again (never zooms)
                self.toggle_at.emit(idx); ev.accept(); return
            super(_AnnView, self).mouseDoubleClickEvent(ev)

    _SHORT_ROLE = QtCore.Qt.UserRole + 1

    class _CatCombo(QtWidgets.QComboBox):
        """Category dropdown: the popup lists 'Variants (12 / 811 shown)'; the closed box paints the item's
        short text ('Variants 12/811') so it stays narrow (the search box gets the width)."""

        def paintEvent(self, _ev):
            p = QtWidgets.QStylePainter(self)
            opt = QtWidgets.QStyleOptionComboBox()
            self.initStyleOption(opt)
            short = self.itemData(self.currentIndex(), _SHORT_ROLE)
            if short:
                opt.currentText = short
            p.drawComplexControl(QtWidgets.QStyle.CC_ComboBox, opt)
            p.drawControl(QtWidgets.QStyle.CE_ComboBoxLabel, opt)

    class _UFVPanel(QtWidgets.QWidget):
        _layers_changed = QtCore.Signal()      # emitted by the module's layer hook (queued across threads)

        def __init__(self):
            super(_UFVPanel, self).__init__()
            self._workers = set()
            self._sync_pending = False
            self._pick_timer = None
            self._picker = None
            self.setWindowTitle("3D Feature Viewer")
            self.setStyleSheet(_UFV_QSS)
            self._settings = _ufv_settings()
            self._sections = {}
            # The window size belongs to the user: it is set once (remembered size, else 440 x 900) and the
            # content never drives it (no layout size constraint; fixed minimum). Everything that changes
            # with the data elides or scrolls inside its own fixed-height box; the report takes the rest.
            self.setMinimumSize(340, 360)   # narrow PyMOL docks: rows wrap instead of scrolling sideways
            self.resize(440, 900)
            if self._settings is not None:
                try:
                    g = self._settings.value("window/geometry")
                    if g is not None:
                        self.restoreGeometry(g)
                except Exception:
                    pass

            outer = QtWidgets.QVBoxLayout(self); outer.setContentsMargins(0, 0, 0, 0); outer.setSpacing(0)
            outer.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
            # Slim header bar: title + subtitle on one line.
            header = QtWidgets.QWidget(); header.setObjectName("ufvHeader")
            header.setAttribute(QtCore.Qt.WA_StyledBackground, True)
            hl = QtWidgets.QHBoxLayout(header); hl.setContentsMargins(10, 4, 10, 4); hl.setSpacing(8)
            _title = QtWidgets.QLabel("3D Feature Viewer"); _title.setObjectName("ufvTitle")
            _sub = _ElidedLabel("UniProt annotations & predictors on structures and trajectories")
            _sub.setObjectName("ufvSub")
            hl.addWidget(_title); hl.addWidget(_sub, 1)
            outer.addWidget(header)

            scroll = QtWidgets.QScrollArea(); scroll.setWidgetResizable(True); scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
            scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)   # rows wrap; never scroll sideways
            self._scroll = scroll
            outer.addWidget(scroll, 1)
            bodyw = QtWidgets.QWidget(); scroll.setWidget(bodyw)
            lay = QtWidgets.QVBoxLayout(bodyw); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)
            self._body_lay = lay

            def hrow(parent):
                r = QtWidgets.QHBoxLayout(); r.setContentsMargins(0, 0, 0, 0); r.setSpacing(4)
                parent.addLayout(r)
                return r

            def flow(parent):
                f = _FlowLayout(); parent.addLayout(f)
                return f

            # ── Protein: accession, entry summary, AI name, tractability ──
            pl = self._section("protein", "Protein", "UniProt accession, entry summary, ProtNLM name, "
                                                     "Open Targets tractability").body_lay
            r1 = hrow(pl)
            self.uid_edit = QtWidgets.QLineEdit(_STATE.get("uid") or "")
            self.uid_edit.setPlaceholderText("UniProt accession, e.g. P35498")
            self.uid_edit.setToolTip("UniProt accession — press Enter to fetch")
            self.uid_edit.returnPressed.connect(lambda: self.fetch_btn.isEnabled() and self.on_fetch())
            r1.addWidget(self.uid_edit, 1)
            self.fetch_btn = QtWidgets.QPushButton("Fetch"); self.fetch_btn.setObjectName("primary"); r1.addWidget(self.fetch_btn)
            self.fetch_btn.setToolTip("Fetch UniProt annotations and the list of structures (Enter)")
            self.fetch_btn.clicked.connect(self.on_fetch)
            self.uniprot_btn = QtWidgets.QToolButton(); self.uniprot_btn.setText("UniProt ↗")
            self.uniprot_btn.setToolTip("Open this accession's UniProtKB entry page in your browser")
            self.uniprot_btn.clicked.connect(self.open_uniprot); r1.addWidget(self.uniprot_btn)
            self.cb_protnlm = _chip("AI name", "Show UniProt's ProtNLM AI-generated protein name/function "
                                               "(predicted from sequence)")
            self.cb_protnlm.clicked.connect(self.toggle_protnlm); r1.addWidget(self.cb_protnlm)
            # one line: the entry summary, or the ProtNLM name while "AI name" is on (full text on hover)
            ir = hrow(pl)
            self.info = _ElidedLabel("")     # one line; full summary + function / location on hover
            self.info.setStyleSheet("font-size:11px; color:#555;")
            ir.addWidget(self.info, 1)
            self.protnlm_label = _ElidedLabel("")
            self.protnlm_label.setStyleSheet("font-size:11px; color:#5b3a8e;"); self.protnlm_label.setVisible(False)
            ir.addWidget(self.protnlm_label, 1)
            # Open Targets target tractability (filled in the background after Fetch) — ufv_tractability.
            # One reserved line (elided, bucket list on hover) + the Open Targets link: never re-flows.
            tr0 = hrow(pl)
            self.tract_label = _ElidedLabel("")
            self.tract_label.setStyleSheet("font-size:10px; color:#37474f;"); self.tract_label.setVisible(False)
            self.tract_label.setToolTip("Open Targets tractability — how amenable this whole protein is to a small "
                                        "molecule, antibody or PROTAC (positive evidence only)")
            tr0.addWidget(self.tract_label, 1)
            self.tract_link = QtWidgets.QLabel(""); self.tract_link.setOpenExternalLinks(True)
            self.tract_link.setTextFormat(QtCore.Qt.RichText); self.tract_link.setStyleSheet("font-size:10px;")
            self.tract_link.setVisible(False)
            tr0.addWidget(self.tract_link)
            _keep_space(self.tract_label, self.tract_link)

            # ── Structures: selector, loaded objects, numbering status ──
            sl = self._section("structures", "Structures", "Choose / load / align structures; the loaded "
                                                           "list sets the active one").body_lay
            r2 = hrow(sl)
            self.struct_combo = QtWidgets.QComboBox()
            self.struct_combo.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
            self.struct_combo.setMinimumContentsLength(6)
            self.struct_combo.setToolTip("Structures of this accession (AlphaFold, PDB, AlphaFill, isoforms)")
            r2.addWidget(self.struct_combo, 1)
            self.load_btn = QtWidgets.QPushButton("Load"); r2.addWidget(self.load_btn)
            self.load_btn.setToolTip("Download and load the selected structure")
            self.load_btn.clicked.connect(self.on_load_selected)
            self.loadall_btn = QtWidgets.QPushButton("All")
            self.loadall_btn.setToolTip("Load all — download and load every listed structure")
            self.loadall_btn.clicked.connect(self.on_load_all); r2.addWidget(self.loadall_btn)
            self.align_btn = QtWidgets.QPushButton("Align"); self.align_btn.setToolTip("Superpose all loaded structures")
            self.align_btn.clicked.connect(lambda: ufv_align()); r2.addWidget(self.align_btn)
            self.num_btn = QtWidgets.QToolButton(); self.num_btn.setText("№")
            self.num_btn.setToolTip("Residue numbering (for structures/trajectories you loaded yourself)")
            self.num_btn.clicked.connect(self.open_numbering_dialog); r2.addWidget(self.num_btn)
            # active structure + % modelled: shown in this section's header (always visible, even collapsed)
            self.obj_label = _ElidedLabel(""); self.obj_label.setVisible(False)
            self.obj_list = QtWidgets.QListWidget()     # fixed 2-row box (scrolls), present even when empty
            self.obj_list.setToolTip("Loaded structures: tick = visible, click a name to make it active, ✕ removes")
            sl.addWidget(self.obj_list)
            # numbering / chimera / partner status of the active structure (one line; full text on hover)
            mr = hrow(sl)
            self.map_label = _ElidedLabel("")
            self.map_label.setStyleSheet("font-size:11px; color:#455;")
            mr.addWidget(self.map_label, 1)
            self.chim_btn = QtWidgets.QToolButton(); self.chim_btn.setText("⚛ unmapped")
            self.chim_btn.setToolTip("Highlight residues of this protein's chains that have no UniProt mapping "
                                     "(chimeric / fusion partner, tags) in blue-grey — ufv_chimera")
            self.chim_btn.clicked.connect(self.show_chimera); self.chim_btn.setVisible(False)
            mr.addWidget(self.chim_btn)
            _keep_space(self.chim_btn)
            self._map_row = mr

            # ── Colouring: ONE cartoon-colour dropdown (+ pocket q, all pockets, apply-to-all) ──
            lgl = self._section("colour", "Colouring", "Cartoon colouring of the active structure (one mode at "
                                                       "a time); 'all' applies colouring and annotations to every "
                                                       "loaded structure").body_lay
            cr0 = flow(lgl)                 # wraps when the q spin box appears on a narrow dock
            _cl = QtWidgets.QLabel("Cartoon"); _cl.setStyleSheet("font-size:11px; color:#6b747d;"); cr0.addWidget(_cl)
            self.cart = QtWidgets.QComboBox()
            self.cart.addItems(["None", "Domains", "Topology", "pLDDT", "B-factor", "AlphaMissense",
                                "Burden", "Hotspots", "Contact hubs", "Constraint pocket"])
            self.cart.setMaximumWidth(140)  # fits "Constraint pocket"; keeps the q spin on the same row
            self.cart.setToolTip("Cartoon colouring of the active structure (analyses run in the background, "
                                 "then are cached)")
            cr0.addWidget(self.cart)
            self.cart.currentIndexChanged.connect(lambda _i: self.apply_cartoon())
            self.q_label = QtWidgets.QLabel("q ≤"); cr0.addWidget(self.q_label)
            self.q_label.setToolTip("Constraint-pocket sensitivity (BH q)")
            self.q_spin = QtWidgets.QDoubleSpinBox(); self.q_spin.setRange(0.01, 0.40)
            self.q_spin.setSingleStep(0.01); self.q_spin.setDecimals(2); self.q_spin.setValue(POCKET_SETTINGS["q"])
            self.q_spin.setToolTip("Constraint-pocket FDR threshold (BH q). Every candidate's q is cached, so "
                                   "changing it recolours instantly — ufv_pocket_q")
            cr0.addWidget(self.q_spin)
            self.q_spin.valueChanged.connect(self.on_pocket_q)
            self.q_label.setVisible(False); self.q_spin.setVisible(False)
            _keep_space(self.q_label, self.q_spin)
            self.allpk_btn = QtWidgets.QPushButton("All pockets")
            self.allpk_btn.setToolTip("Every ProtVar (AutoSite) pocket of the protein: scans all residues once (slow "
                                      "the first time, cached), colours each pocket and lists them — ufv_pockets_all")
            self.allpk_btn.clicked.connect(self.on_all_pockets); cr0.addWidget(self.allpk_btn)
            self.cb_apply_all = QtWidgets.QCheckBox("all")
            self.cb_apply_all.setToolTip("Apply to all — colouring, annotation checks and zoom act on every loaded "
                                         "structure (e.g. aligned copies), not only the active one")
            cr0.addWidget(self.cb_apply_all)
            self.aff_combo = QtWidgets.QComboBox()
            self.aff_combo.addItems(["id ≥ %d%%" % p for p in AFF_IDENTITY_STEPS])
            self.aff_combo.setToolTip("AlphaFill: show transplants whose donor has at least this sequence identity "
                                      "(25% = AlphaFill's floor = all) — ufv_ligand_identity")
            self.aff_combo.setVisible(False)
            self.aff_combo.currentIndexChanged.connect(self.on_aff_identity)
            # AlphaFill-only identity preset: in the Structures section's numbering row, next to the elided
            # numbering line (which gives up the width), so showing it never re-flows anything
            self._map_row.addWidget(self.aff_combo)
            _keep_space(self.aff_combo)
            self.cart.currentIndexChanged.connect(lambda _i: self._layers_summary())

            # ── Trajectory (only shown when the active object has more than one state) ──
            self._traj_sec = self._section("trajectory", "Trajectory", "Multi-state object / MD trajectory tools")
            self._traj_slot = self._traj_sec.body_lay          # QVBoxLayout the trajectory panel goes into
            self._traj_widget = None; self._traj_failed = False
            self._traj_sec.setVisible(False)

            # ── Annotations v3: category dropdown + search + Filters ▾ + ↺, ONE checkable list ──
            # Checked row == drawn in 3D (per object + accession); checks are kept per category, so switching
            # the dropdown never changes the 3D scene. The search text IS the filter (ufv_core.ann_search, the
            # same matcher as the VMD panel); Filters ▾ only writes search terms.
            agl = self._section("annotations", "Annotations", "Tick = shown in 3D · click a row = report + highlight "
                                                              "(camera stays) · double-click / Enter = zoom · Space "
                                                              "toggles the selected rows").body_lay
            tr = hrow(agl)
            self.ann_cat = _CatCombo()
            self.ann_cat.addItem("All", "all")
            for _c in ANN_CATS:
                self.ann_cat.addItem(ANN_CAT_LABELS[_c], _c)
            self.ann_cat.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
            self.ann_cat.setMinimumContentsLength(11)
            self.ann_cat.setMinimumWidth(92); self.ann_cat.setMaximumWidth(150)   # never drives the body width
            self.ann_cat.view().setMinimumWidth(190)          # the popup shows the full 'x / y shown' texts
            self.ann_cat.setToolTip("Annotation category (checked / total rows shown in 3D). Every category keeps "
                                    "its own checks; switching here never changes the 3D view")
            self.ann_cat.currentIndexChanged.connect(lambda _i: self.on_ann_category())
            tr.addWidget(self.ann_cat)
            self.ann_search = QtWidgets.QLineEdit()
            self.ann_search.setPlaceholderText("search: 1378, 100-200, N1378H, pathogenic, \"phrase\"")
            self.ann_search.setClearButtonEnabled(True)
            self.ann_search.setMinimumWidth(60)
            self.ann_search.setToolTip(
                "Filter the list (all terms must match, case-insensitive):\n"
                "  1378 · 100-200 · N1378H · R2K — position / range / mutation\n"
                "  pathogenic · likely · benign · uncertain · deleterious — consequence\n"
                "  uniprot · clinvar · cancer · population · reviewed — provenance\n"
                "  disease name or id, description text, PTM type, ligand code, category\n"
                "  \"quoted phrase\" · a|b = either · facets cons: src: af: sub: dis: cat: pos:\n"
                "Filters ▾ writes these terms for you. The header box ticks every listed row.")
            self._search_timer = QtCore.QTimer(self); self._search_timer.setSingleShot(True)
            self._search_timer.setInterval(120)              # debounce: one filter pass per typing pause
            self._search_timer.timeout.connect(self._apply_search)
            self.ann_search.textChanged.connect(lambda _t: self._search_timer.start())
            self.ann_search.textEdited.connect(self._update_completer)
            tr.addWidget(self.ann_search, 1)
            self._sugg_model = QtCore.QStringListModel(self)
            self._completer = QtWidgets.QCompleter(self._sugg_model, self)
            self._completer.setCompletionMode(QtWidgets.QCompleter.UnfilteredPopupCompletion)
            self._completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
            self._completer.setWidget(self.ann_search)
            self._completer.activated[str].connect(self._take_suggestion)
            self.ann_filters_btn = QtWidgets.QToolButton(); self.ann_filters_btn.setText("Filters ▾")
            self.ann_filters_btn.setToolTip("Consequence / source / allele frequency / substitution / disease filters "
                                            "and the variant colour profile — they only write search terms")
            self.ann_filters_btn.setPopupMode(QtWidgets.QToolButton.InstantPopup)
            self.ann_filters_btn.setStyleSheet("QToolButton::menu-indicator { image: none; width: 0px; }")
            self.ann_filters_menu = QtWidgets.QMenu(self.ann_filters_btn)
            self.ann_filters_menu.aboutToShow.connect(self._build_filters_menu)
            self.ann_filters_btn.setMenu(self.ann_filters_menu)
            tr.addWidget(self.ann_filters_btn)
            self.ann_reset_btn = QtWidgets.QToolButton(); self.ann_reset_btn.setText("↺")
            self.ann_reset_btn.setToolTip("Uncheck all categories — nothing shown in 3D")
            self.ann_reset_btn.clicked.connect(self.on_ann_reset)
            tr.addWidget(self.ann_reset_btn)
            self.ann_custom_btn = QtWidgets.QToolButton(); self.ann_custom_btn.setText("+")
            self.ann_custom_btn.setToolTip("Load custom annotations (CSV / TSV: position, label, value, color, note) "
                                           "— shown as the Custom category")
            self.ann_custom_btn.clicked.connect(self.on_ann_custom)
            tr.addWidget(self.ann_custom_btn)

            self.ann_model = _AnnModel(self)
            self.ann_proxy = _AnnProxy(self); self.ann_proxy.setSourceModel(self.ann_model)
            self.ann_view = _AnnView(self)
            self.ann_head = _CheckHeader(self.ann_view)
            self.ann_view.setHorizontalHeader(self.ann_head)
            self.ann_view.setModel(self.ann_proxy)
            hh = self.ann_view.horizontalHeader()
            hh.setStretchLastSection(True); hh.setHighlightSections(False)
            hh.setSectionResizeMode(0, QtWidgets.QHeaderView.Fixed); hh.resizeSection(0, 24)
            hh.resizeSection(1, 34); hh.resizeSection(2, 66)
            self.ann_view.setColumnHidden(1, True)            # category column: "All" only
            self.ann_view.verticalHeader().setVisible(False)
            self.ann_view.verticalHeader().setMinimumSectionSize(16)
            self.ann_view.verticalHeader().setDefaultSectionSize(19)
            self.ann_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.ann_view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.ann_view.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
            self.ann_view.setWordWrap(False); self.ann_view.setShowGrid(False)
            self.ann_view.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            self.ann_view.setFixedHeight(232)    # rows scroll inside; the section never changes height
            self.ann_view.clicked.connect(self.on_ann_click)
            self.ann_view.doubleClicked.connect(self.on_ann_zoom)
            self.ann_view.space.connect(self.toggle_selected)
            self.ann_view.enter.connect(lambda: self.on_ann_zoom(self.ann_view.currentIndex()))
            self.ann_view.toggle_at.connect(self.toggle_row_at)
            self.ann_head.toggled.connect(self.toggle_all_listed)
            agl.addWidget(self.ann_view, 0)
            self._ann_ctx_key = None          # (obj, uid, category) the list was built for
            self._ann_rows = []               # rows of the current category (ufv_core.ann_rows)
            self._ann_vis = []                # rows passing the search
            self._ann_avail = {}              # {cat: total rows | None (not fetched)}
            self._ann_fetching = set()

            # ── Residue report (takes the remaining height) ──
            self._report_sec = self._section("report", "Residue report", "Report of the clicked / picked residue "
                                                                         "or ligand", stretch=1)
            rgl = self._report_sec.body_lay
            br = flow(rgl)
            self.prox_spin = {}
            for key, lbl, tip in (("ptm", "PTM Å", "Residue report: list PTMs within this Cα distance (also the "
                                                   "PTM-variant tier-2 cut-off)"),
                                  ("var", "Var Å", "Residue report: list pathogenic variants within this Cα distance "
                                                   "(also the PTM-variant tier-3 cut-off)")):
                sp = QtWidgets.QSpinBox(); sp.setRange(2, 30); sp.setValue(int(PROX_SETTINGS[key]))
                sp.setToolTip(tip + " — ufv_prox"); sp.setPrefix(lbl + " ")
                sp.valueChanged.connect(self.on_prox_radius); self.prox_spin[key] = sp
                br.addWidget(sp)
            self.nb_spin = QtWidgets.QSpinBox(); self.nb_spin.setRange(int(NEARBY_MIN), int(NEARBY_MAX))
            self.nb_spin.setValue(int(NEARBY_SETTINGS["dist"])); self.nb_spin.setPrefix("Nearby Å ")
            self.nb_spin.setToolTip("Residue / ligand focus neighbourhood (any-atom distance) and the report's "
                                    "Nearby list — ufv_nearby")
            self.nb_spin.valueChanged.connect(self.on_nearby); br.addWidget(self.nb_spin)
            self.cb_pick = _chip("Pick 3D"); br.addWidget(self.cb_pick)
            self.cb_pick.setToolTip("Click any atom in the PyMOL viewer to zoom to its residue/ligand")
            self.cb_pick.clicked.connect(self.toggle_pick)
            self.reset_btn = QtWidgets.QPushButton("Reset view"); br.addWidget(self.reset_btn)
            self.reset_btn.setToolTip("Clear the focus and return to the view from before it (whole structure "
                                     "when there was no focus or the camera was moved meanwhile)")
            self.reset_btn.clicked.connect(lambda: ufv_resetview(self.obj()))

            self.detail = _ReportBrowser(); self.detail.setMinimumHeight(140)
            self.detail.setOpenLinks(False); self.detail.setOpenExternalLinks(False)
            self.detail.anchorClicked.connect(self.on_anchor)
            rgl.addWidget(self.detail, 1)   # the residue report takes the spare vertical space
            # a new report / ligand / pocket listing re-opens the report section if it was collapsed
            self.detail.textChanged.connect(lambda: self._report_sec.setExpanded(True, save=False))
            lay.addStretch(0)               # takes the slack only while the report section is collapsed

            # ── Footer (outside the scroll area, always visible): status + exports ──
            foot = QtWidgets.QWidget(); foot.setObjectName("ufvFooter")
            foot.setAttribute(QtCore.Qt.WA_StyledBackground, True)
            fl = QtWidgets.QVBoxLayout(foot); fl.setContentsMargins(8, 2, 8, 4); fl.setSpacing(2)
            self.status = _ElidedLabel(""); self.status.setStyleSheet("color:#a33; font-size:11px;")
            fl.addWidget(self.status)
            bot = hrow(fl)
            self.csv_btn = QtWidgets.QPushButton("CSV")
            self.csv_btn.setToolTip("Export CSV — save the per-residue annotation matrix (one-hot PTM/disease/site, "
                                    "variant counts, gnomAD, AlphaMissense, burden, structure tiers)")
            self.csv_btn.clicked.connect(lambda: self.on_export_csv(False)); bot.addWidget(self.csv_btn)
            self.csv_pv_btn = QtWidgets.QPushButton("CSV+ProtVar")
            self.csv_pv_btn.setToolTip("CSV including ProtVar conservation/EVE/ESM1b/M3D — fetches one "
                                       "request per residue, so it is slow")
            self.csv_pv_btn.clicked.connect(lambda: self.on_export_csv(True)); bot.addWidget(self.csv_pv_btn)
            self.bfac_btn = QtWidgets.QPushButton("B-factor PDB")
            self.bfac_btn.setToolTip("Save the active structure as PDB with the current cartoon colouring encoded "
                                     "in the B-factor column (AM / hotspots / hubs / burden / pocket; otherwise the "
                                     "residues shown by the marker layers) — ufv_save_bfactor")
            self.bfac_btn.clicked.connect(self.on_export_bfactor); bot.addWidget(self.bfac_btn)
            self.clear_btn = QtWidgets.QPushButton("Clear"); bot.addWidget(self.clear_btn, 1)
            self.clear_btn.setToolTip("Clear all overlays (layers, colouring, focus) — ufv_clear")
            self.clear_btn.clicked.connect(self.on_clear)
            for _b in (self.csv_btn, self.csv_pv_btn, self.bfac_btn, self.clear_btn):
                _b.setMinimumWidth(10)
            outer.addWidget(foot)
            self._layers_summary()

            self._QtGui = QtGui
            self._last_pos = None
            self._sel = None               # unified current selection: {"kind":"res","pos":n} or
                                           # {"kind":"lig","resn":r} — drives 'Zoom selected'
            self._last_rep = None          # cached residue_report dict for re-rendering (evidence toggle)
            self._last_protvar = None      # cached ProtVar fetch for the current residue ('off'/'loading'/dict)
            self._show_evidence = False    # variant-evidence fold state in the detail panel
            self._lig_resn = None          # ligand component currently shown
            self._lig_copies = []          # [(chain, resi), ...] occurrences of that component
            self._lig_idx = 0              # focused occurrence
            self._lig_chem = {}            # cached chemistry for the current ligand
            self._all_pockets = None       # last ufv_pockets_all result shown in the panel
            self.ensurePolished()          # panel style sheet applied to every child: real font metrics below
            self._fit_obj_list()           # fixed 2-row loaded-structures box
            self._install_traj_panel()     # every widget exists now; no-op until the factory is defined
            if self.cur_uid() and self.cur_uid() in _CACHE:
                self.after_fetch()
            self.update_traj_section()
            _ufv_wheel_guard(self)          # wheel over a combo / spin box scrolls the panel, not the value
            self._layers_changed.connect(self._sched_sync)
            self._layer_listener = self._layers_changed.emit
            _LAYER_LISTENERS.append(self._layer_listener)
            self.destroyed.connect(lambda *_a, fn=self._layer_listener:
                                   fn in _LAYER_LISTENERS and _LAYER_LISTENERS.remove(fn))
            self.sync_layers()

        # ---- compact layout: sections, summaries, trajectory slot ----
        def _sched_sync(self):
            if not self._sync_pending:
                self._sync_pending = True
                QtCore.QTimer.singleShot(0, self.sync_layers)

        # chips <- the layer bookkeeping tag they show (the module's _SPHERE_STATE / _LIG_SHOWN)
        _CART_FROM_LAYER = {"domains": "Domains", "topology": "Topology", "plddt": "pLDDT", "bfactor": "B-factor",
                            "alphamissense": "AlphaMissense", "burden": "Burden", "hotspots": "Hotspots",
                            "contacthubs": "Contact hubs", "pockets": "Constraint pocket",
                            "constraint pocket": "Constraint pocket"}

        def sync_layers(self):
            """Make the cartoon combo, the annotation list's check boxes / counts and the section summaries
            reflect the real state (set by console commands as well as by this panel). Reads dicts only;
            signals are blocked so no slot re-runs."""
            self._sync_pending = False
            try:
                o = self.obj()
                if self.cart.isEnabled():
                    txt = self._CART_FROM_LAYER.get(_CARTOON_LAYER.get(o), "None")
                    i = self.cart.findText(txt)
                    if i >= 0 and i != self.cart.currentIndex():
                        was = self.cart.blockSignals(True)
                        self.cart.setCurrentIndex(i)
                        self.cart.blockSignals(was)
                    pk = self.cart.currentText() == "Constraint pocket"
                    self.q_label.setVisible(pk); self.q_spin.setVisible(pk)
                self._layers_summary()
                self._ann_sync()
            except RuntimeError:
                pass

        def _section(self, key, title, tip="", stretch=0):
            sec = _Section(key, title, tip, self._settings)
            self._sections[key] = sec
            self._body_lay.addWidget(sec, stretch)
            return sec

        def _layers_summary(self):
            cart = self.cart.currentText()
            self._sections["colour"].setSummary(("cartoon: " + cart) if cart and cart != "None" else "no colouring")

        # ---- Annotations v3 ----
        def _ann_ctx(self):
            return self.obj(), (self.cur_uid() or "").upper()

        def _ann_modelled(self, o):
            if not o or o not in _object_set():
                return None
            try:
                return {g["uniPos"] for g in _ca_geometry(o) if g["uniPos"] is not None}
            except Exception:
                return None

        def ann_reload(self, force=False):
            """(Re)build the list for the active object / accession / category; the search is re-applied."""
            o, uid = self._ann_ctx()
            cat = self.ann_cat.currentData() or "all"
            key = (o, uid, cat)
            if not force and key == self._ann_ctx_key:
                return self._ann_sync()
            self._ann_ctx_key = key
            self._ann_avail = ann_available(o, uid) if (o or uid in _CACHE) else {}
            rows = annotation_rows(o, uid, cat) if (uid in _CACHE or o) else []
            mod = self._ann_modelled(o)
            unmod = frozenset(r["id"] for r in rows if mod is not None and r["pos"] is not None
                              and r["cat"] not in ("kb",) and not any(p in mod for p in r["pts"][:1] or [r["pos"]]))
            self._ann_rows = rows
            self.ann_model.set_rows(rows, ann_checked(o, uid) if o else {}, unmod)
            self.ann_view.setColumnHidden(1, cat != "all")
            self._apply_search()
            self._ann_sync(counts_only=True)

        def _apply_search(self):
            """One ann_search pass over the category's rows -> proxy filter (< 30 ms on 4,000 rows)."""
            self._search_timer.stop()
            q = self.ann_search.text()
            if q.strip():
                ids = ann_search(self._ann_rows, q)
                ok = set(ids)
                self._ann_vis = [r for r in self._ann_rows if r["id"] in ok]
                self.ann_proxy.set_ok(ok)
            else:
                self._ann_vis = self._ann_rows
                self.ann_proxy.set_ok(None)
            self._ann_header()

        def _ann_header(self):
            on = self.ann_model.on
            n = sum(1 for r in self._ann_vis if r["id"] in (on.get(r["cat"]) or ()))
            v = len(self._ann_vis)
            self.ann_head.set_state(_UNCHECKED if not n else _CHECKED if n == v else _PARTIAL)
            head = "%d listed · %d shown" % (v, n) if self._ann_rows else "Annotation"
            if head != self.ann_model.head:
                self.ann_model.head = head
                self.ann_model.headerDataChanged.emit(QtCore.Qt.Horizontal, 3, 3)

        def _ann_sync(self, counts_only=False):
            """Check boxes + dropdown counts + header + section summary <- module state (cheap)."""
            o, uid = self._ann_ctx()
            cat = self.ann_cat.currentData() or "all"
            if not counts_only and (o, uid, cat) != self._ann_ctx_key:
                return self.ann_reload()
            self.ann_model.on = ann_checked(o, uid) if o else {}
            if not counts_only:
                self.ann_model.refresh_checks()
                self._ann_header()
            counts = ann_counts(o, uid) if o else {c: 0 for c in ANN_CATS}
            avail = self._ann_avail or {}
            was = self.ann_cat.blockSignals(True)
            m = self.ann_cat.model()
            tot_all = sum(v for v in avail.values() if v)
            for i in range(self.ann_cat.count()):
                c = self.ann_cat.itemData(i)
                if c == "all":
                    txt, ok = ("All (%d / %d shown)" % (sum(counts.values()), tot_all)) if tot_all else "All", True
                else:
                    tot = avail.get(c)
                    lbl = ANN_CAT_LABELS[c]
                    if tot is None:          # partners / PDBe-KB: fetched when picked
                        ok = bool(o) and (c != "prt" or bool(_partner_maps(o))) and (c != "kb" or bool(uid))
                        txt = lbl + (" (…)" if ok else "")
                    else:
                        ok = tot > 0
                        txt = "%s (%d / %d shown)" % (lbl, counts.get(c, 0), tot) if tot else lbl
                short = txt.replace(" shown)", "").replace(" / ", "/").replace(" (", " ").rstrip(")")
                if self.ann_cat.itemText(i) != txt:
                    self.ann_cat.setItemText(i, txt)
                    self.ann_cat.setItemData(i, short, _SHORT_ROLE)
                it = m.item(i) if hasattr(m, "item") else None
                if it is not None and it.isEnabled() != ok:
                    it.setEnabled(ok)
            self.ann_cat.blockSignals(was)
            bits = ["%s %d" % (ANN_CAT_SHORT[c], counts[c]) for c in ANN_CATS if counts.get(c)]
            n = sum(counts.values())
            self._sections["annotations"].setSummary(" · ".join(bits) if bits else "nothing shown in 3D")
            self.ann_reset_btn.setToolTip("Uncheck all categories — %d shown in 3D (colouring and camera are kept) "
                                          "— ufv_reset_annotations" % n)

        def on_ann_category(self):
            cat = self.ann_cat.currentData() or "all"
            o, uid = self._ann_ctx()
            if cat in ("prt", "kb") and (self._ann_avail or {}).get(cat, 0) is None and o:
                self._ann_fetch(cat, o, uid)
            self.ann_reload()

        def _ann_fetch(self, cat, o, uid):
            """Partner / PDBe-KB rows need one network fetch (worker); the list refreshes when it lands."""
            if (cat, o, uid) in self._ann_fetching:
                return
            self._ann_fetching.add((cat, o, uid))
            work = (lambda: _partner_annotations(o)) if cat == "prt" else (lambda: fetch_pdbekb(uid))

            def done(_res):
                self._ann_fetching.discard((cat, o, uid))
                self.ann_reload(force=True)
            self._async(work, done, "Fetching %s rows ..." % ANN_CAT_LABELS[cat])

        def _row_of(self, pidx):
            if pidx is None or not pidx.isValid():
                return None
            s = self.ann_proxy.mapToSource(pidx)
            return self.ann_model.rows[s.row()] if 0 <= s.row() < len(self.ann_model.rows) else None

        def _selected_rows(self):
            sm = self.ann_view.selectionModel()
            return [self._row_of(i) for i in (sm.selectedRows(0) if sm else [])]

        def ann_commit(self, rows, on):
            """Check (on=True) / uncheck rows: ONE batched update per object (all touched categories)."""
            rows = [r for r in rows if r]
            o, uid = self._ann_ctx()
            if not rows or not o:
                return
            changes = {}
            for r in rows:
                changes.setdefault(r["cat"], [[], on])[0].append(r["id"])
            objs = self.target_objs() or [o]
            if on and "lig" in changes:
                need = [x for x in objs if _is_alphafill(x) and x not in _AFF_META]
                if need:                     # AlphaFill transplant metadata: network in a worker, then draw
                    def apply(meta):
                        for x in need:
                            _AFF_META[x] = meta or {}
                        self.ann_commit(rows, on)
                    self._async(lambda: alphafill_transplants(uid), apply, "Fetching AlphaFill transplant data ...")
                    return
            if "var" in changes:            # 'disease' colour profile: first disease named in the search
                fac = ann_query_facets(self.ann_search.text()).get("dis")
                known = {d.lower(): d for d in ((_CACHE.get(uid) or {}).get("diseaseColors") or {})}
                dis = {known.get(d.lower(), d) for d in fac} if fac else None
                for x in objs:
                    _ann_style(x)["dis"] = dis
            with _no_undo():
                for x in objs:
                    ann_set_checked(x, uid, {c: (ids, on) for c, (ids, on) in changes.items()})
            self._ann_sync()

        def toggle_row_at(self, pidx):
            r = self._row_of(pidx)
            if r is None:
                return
            sel = [x for x in self._selected_rows() if x]
            rows = sel if (len(sel) > 1 and any(x["id"] == r["id"] for x in sel)) else [r]
            self.ann_commit(rows, not self.ann_model.is_on(r))

        def toggle_selected(self):
            rows = [x for x in self._selected_rows() if x]
            if not rows:
                r = self._row_of(self.ann_view.currentIndex())
                rows = [r] if r else []
            if rows:
                self.ann_commit(rows, not all(self.ann_model.is_on(r) for r in rows))

        def toggle_all_listed(self):
            """Header box: tick every listed (searched) row of this category, or untick them all if all are on."""
            vis = self._ann_vis
            if vis:
                self.ann_commit(vis, not all(self.ann_model.is_on(r) for r in vis))

        def on_ann_custom(self):
            o, uid = self._ann_ctx()
            if not uid:
                return
            path, _f = QtWidgets.QFileDialog.getOpenFileName(self, "Custom annotations for %s" % uid, "",
                                                            "Annotation tables (*.csv *.tsv *.txt);;All files (*)")
            if not path:
                return
            with _no_undo():
                ufv_load_annotations(path, o, uid)
            i = self.ann_cat.findData("cus")
            if i >= 0:
                self.ann_cat.setCurrentIndex(i)
            self.ann_reload(force=True)

        def on_ann_reset(self):
            o = self.obj()
            if not o:
                return
            with _no_undo():
                for x in (self.target_objs() or [o]):
                    ufv_reset_annotations(x)
            self._ann_sync()

        def on_ann_click(self, pidx):
            """Single click: report + transient highlight of the row's residue / ligand. Never moves the camera."""
            if pidx.column() == 0:
                return
            self._ann_open(self._row_of(pidx), zoom=False)

        def on_ann_zoom(self, pidx):
            """Double click / Enter: zoom to the row (explicit view action)."""
            self._ann_open(self._row_of(pidx), zoom=True)

        def _ann_open(self, r, zoom):
            if r is None:
                return
            if r["cat"] == "lig":
                self.on_ligand(r["f"]["resn"], zoom=zoom); return
            if r["cat"] == "prt":
                o = self.obj()
                if o and zoom:
                    with _no_undo():
                        sel = "(%s) and chain %s and resi %s" % (o, r["f"]["chain"], _resi_tok(r["f"]["resi"]))
                        if cmd.count_atoms(sel):
                            _focus_zoom(sel, 6)
                self.status.setText(r["label"]); return
            if r["pos"] is not None:
                self.report_residue(int(r["pos"]), focus=True, zoom=zoom)

        def _update_completer(self, text):
            """Suggestions for the term being typed (consequence / source words, diseases, PTM types ...)."""
            prefix, partial = ann_last_term(text)
            sug = ann_suggest(self._ann_rows, partial) if partial else []
            if not sug:
                self._sugg_model.setStringList([])     # no stale entries from the previous term
                self._completer.popup().hide(); return
            self._sugg_model.setStringList([prefix + s + " " for s in sug])
            self._completer.complete()

        def _take_suggestion(self, text):
            self.ann_search.setText(text)
            self.ann_search.setFocus()
            self._apply_search()

        def _build_filters_menu(self):
            """Filters ▾: checkable shortcuts that rewrite the search terms (the search text IS the filter)."""
            menu = self.ann_filters_menu
            menu.clear()
            q = self.ann_search.text()
            fac = ann_query_facets(q)

            def toggler(key, val, allvals):
                def f(_c=False):
                    cur = fac.get(key)
                    cur = set(allvals) if cur is None else set(cur)
                    cur.symmetric_difference_update({val})
                    self.ann_search.setText(ann_query_set(self.ann_search.text(), key,
                                                          None if cur >= set(allvals) else cur))
                    self._apply_search()
                return f
            menu.addSection("Consequence")
            for tok, lbl in (("pathogenic", "Pathogenic / likely"), ("deleterious", "Predicted deleterious"),
                             ("benign", "Benign / likely"), ("uncertain", "Uncertain significance")):
                a = menu.addAction(lbl); a.setCheckable(True)
                a.setChecked(fac["cons"] is None or tok in fac["cons"])
                a.triggered.connect(toggler("cons", tok, ANN_CONS_TOKENS))
            menu.addSection("Source")
            for tok, lbl in (("uniprot", "UniProt reviewed"), ("clinvar", "ClinVar / ClinGen"), ("cancer", "Cancer somatic (COSMIC/TCGA)"),
                             ("population", "Population (gnomAD/ExAC/TOPMed)"), ("other", "Other large-scale studies")):
                a = menu.addAction(lbl); a.setCheckable(True)
                a.setChecked(fac["src"] is None or tok in fac["src"])
                a.triggered.connect(toggler("src", tok, ANN_SRC_TOKENS))
            a = menu.addAction("Reviewed only"); a.setCheckable(True); a.setChecked(fac["reviewed"])

            def rev(_c=False, on=not fac["reviewed"]):
                self.ann_search.setText(ann_query_set(self.ann_search.text(), "reviewed", on)); self._apply_search()
            a.triggered.connect(rev)
            for key, title, toks, cats in (("af", "Allele frequency (gnomAD)", ANN_AF_TOKENS, AF_COLORS),
                                           ("sub", "Substitution type", ANN_SUB_TOKENS, SUBST_COLORS)):
                sm = menu.addMenu(title if fac[key] is None else "%s (%d of %d)" % (title, len(fac[key]), len(toks)))
                for tok, (lbl, hexc) in zip(toks, cats.items()):
                    pm = QtGui.QPixmap(10, 10); pm.fill(QtGui.QColor(hexc))
                    a = sm.addAction(QtGui.QIcon(pm), lbl); a.setCheckable(True)
                    a.setChecked(fac[key] is None or tok in fac[key])
                    a.triggered.connect(toggler(key, tok, toks))
            ann = _CACHE.get(self.cur_uid() or "") or {}
            summ = _disease_summary(ann.get("variants") or [], ann.get("diseaseColors")) if ann else {}
            names = [n for n in summ if n != "Unclassified"]
            dm = menu.addMenu("Diseases (%d)" % len(names) if fac["dis"] is None else "Diseases (%d of %d)"
                              % (len(fac["dis"]), len(names)))
            dm.setEnabled(bool(names))
            if names:
                sel = None if fac["dis"] is None else {d.lower() for d in fac["dis"]}

                def setdis(vals):
                    def f(_c=False):
                        self.ann_search.setText(ann_query_set(self.ann_search.text(), "dis", vals))
                        self._apply_search()
                    return f
                dm.addAction("All diseases").triggered.connect(setdis(None))
                dm.addAction("None").triggered.connect(setdis(set()))
                dm.addSeparator()
                for n in names:
                    pm = QtGui.QPixmap(10, 10); pm.fill(QtGui.QColor(summ[n]["color"]))
                    a = dm.addAction(QtGui.QIcon(pm), "%s (%d)" % (n, summ[n]["count"])); a.setCheckable(True)
                    on = sel is None or n.lower() in sel
                    a.setChecked(on)
                    cur = set(names) if sel is None else {x for x in names if x.lower() in sel}
                    a.triggered.connect(setdis(None if (cur ^ {n}) >= set(names) else (cur ^ {n})))
            cm = menu.addMenu("Colour variants by")
            o = self.obj()
            prof = _ann_style(o)["color"] if o else "consequence"
            grp = QtWidgets.QActionGroup(cm)
            for nm in VARIANT_COLOR_PROFILES:
                a = cm.addAction(nm); a.setCheckable(True); a.setChecked(prof.startswith(nm[:3])); grp.addAction(a)
                a.triggered.connect(lambda _c=False, nm=nm: self.set_variant_colour(nm))
            menu.addSeparator()
            menu.addAction("Clear search").triggered.connect(lambda: (self.ann_search.clear(), self._apply_search()))

        def set_variant_colour(self, profile):
            """Variant sphere colour profile (consequence | disease | provenance): recolours checked variants."""
            uid = self.cur_uid()
            with _no_undo():
                for x in self.target_objs():
                    _ann_style(x)["color"] = profile
                    if ann_counts(x, uid).get("var"):
                        _ann_apply(x, uid, ("var",))

        def _install_traj_panel(self):
            """Insert the trajectory panel (_ufv_make_traj_panel(panel) -> QWidget) into self._traj_slot,
            once, if that factory exists. Retried whenever the section is about to be shown."""
            if self._traj_widget is not None or self._traj_failed:
                return self._traj_widget is not None
            fn = globals().get("_ufv_make_traj_panel")
            if fn is None:
                return False
            try:
                w = fn(self)
            except Exception as exc:
                self._traj_failed = True
                print("[UFV] trajectory panel unavailable: %s" % exc)
                return False
            if w is not None:
                self._traj_slot.addWidget(w); self._traj_widget = w
                _ufv_wheel_guard(w)
            return w is not None

        def update_traj_section(self, name=None):
            """Show the Trajectory section only when the active object has more than one state."""
            name = name or self.obj()
            n = 0
            try:
                if name and name in _object_set():
                    n = cmd.count_states(name)
            except Exception:
                n = 0
            multi = n > 1
            if multi:
                self._install_traj_panel()
                self._traj_sec.setSummary("%s · %d states" % (name, n))
            self._traj_sec.setVisible(multi)

        def _fit_obj_list(self):
            """Fixed 2-row height whatever the number of loaded structures (more rows scroll inside)."""
            if not getattr(self, "_obj_list_h", None):
                probe = self._obj_row("X", None)
                rh = max(18, probe.sizeHint().height()); probe.deleteLater()
                self._obj_list_h = 2 * rh + 2 * self.obj_list.frameWidth() + 2
            self.obj_list.setFixedHeight(self._obj_list_h)

        # ---- UniProt link + ProtNLM AI name ----
        def open_uniprot(self):
            uid = self.cur_uid()
            if not uid:
                self.status.setText("Enter an accession first."); return
            try:
                import webbrowser
                webbrowser.open("https://www.uniprot.org/uniprotkb/%s/entry" % uid)
            except Exception as exc:
                self.status.setText("Could not open browser: %s" % exc)

        def _show_protnlm(self, on):
            self.protnlm_label.setVisible(on); self.info.setVisible(not on)

        def toggle_protnlm(self):
            if not self.cb_protnlm.isChecked():
                self._show_protnlm(False); return
            uid = self.cur_uid()
            ann = _CACHE.get(uid or "")
            if not ann:
                self.cb_protnlm.setChecked(False)
                self.status.setText("Fetch the accession first."); return
            self._render_protnlm(ann.get("protnlm"))

        def _render_protnlm(self, p):
            if not self.cb_protnlm.isChecked():
                self._show_protnlm(False); return
            tip = ""
            if not p or not p.get("name"):
                txt = "No protein name in the UniProt entry."
            elif p.get("isAI"):
                txt = "ProtNLM (AI-predicted): %s" % p["name"]
                if p.get("source"):
                    txt += " · %s" % p["source"]
                tip = p.get("caution") or ""
            else:
                tag = "curated" if p.get("reviewed") else "not AI-generated"
                txt = "Protein name: %s (%s — no ProtNLM prediction)" % (p["name"], tag)
            self.protnlm_label.setText(txt); self.protnlm_label.setToolTip(tip)
            self._show_protnlm(True)

        # ---- helpers ----
        def cur_uid(self):
            return self.uid_edit.text().strip().upper() or _STATE.get("uid")

        def obj(self):
            return _STATE.get("obj") or _resolve_object(None)

        def _async(self, work, apply, status, ctrl=None):
            """Run `work` off the UI thread, then `apply(result)` on it. Only the control that
            triggered the work (`ctrl`) is disabled meanwhile — not the whole panel, which used to
            freeze every other widget for the 2-8 s of a network fetch."""
            self.status.setText(status)
            if ctrl is not None:
                ctrl.setEnabled(False)
            wk = _Worker(work)

            def on_done(res):
                if ctrl is not None:
                    ctrl.setEnabled(True)
                self.status.setText(""); self._workers.discard(wk)
                if isinstance(res, Exception):
                    self.status.setText("Error: %s" % res); return
                if apply:
                    with _no_undo():
                        apply(res)
            wk.done.connect(on_done); self._workers.add(wk); wk.start()

        def update_info(self):
            ann = _CACHE.get(self.cur_uid() or "")
            if not ann:
                self.info.setText("Not fetched."); self.info.setToolTip("")
                self._sections["protein"].setSummary(""); return
            self._sections["protein"].setSummary("%s · %d aa%s" % (ann["uid"], len(ann["sequence"]), (
                " · " + (ann.get("protnlm") or {}).get("name")) if (ann.get("protnlm") or {}).get("name") else ""))
            cons = Counter(_cons_token(v["consequence"]) for v in ann["variants"])
            txt = "%s · %d aa · PTM %d · Var %d (path %d) · Site %d · Mutag %d · Dom %d · AM %s" % (
                ann["uid"], len(ann["sequence"]), len(ann["ptms"]), len(ann["variants"]),
                cons.get("pathogenic", 0), len(ann["sites"]), len(ann.get("mutagenesis", [])),
                len(ann["domains"]), "y" if ann["amMean"] else "n")
            fn = ann.get("function") or {}
            if fn.get("locations"):
                txt += "\nLocation: " + ", ".join(fn["locations"][:4])
            self.info.setText(txt)
            tip = []
            if fn.get("summary"):
                tip.append("Function: " + fn["summary"])
            if fn.get("catalytic"):
                tip.append("Catalytic: " + "; ".join(fn["catalytic"][:3]))
            if fn.get("locations"):
                tip.append("Subcellular location: " + ", ".join(fn["locations"]))
            self.info.setToolTip("\n\n".join(tip))

        def loaded_objs(self):
            return [o for o in (cmd.get_object_list() or []) if not o.startswith("ufv_")]

        def target_objs(self):
            """Objects layer/colour/zoom operations act on: all loaded if 'apply to all', else active."""
            if self.cb_apply_all.isChecked():
                return self.loaded_objs()
            o = self.obj()
            return [o] if o else []

        def refresh_objs(self):
            active = self.obj()
            self.obj_list.clear()
            for o in self.loaded_objs():
                item = QtWidgets.QListWidgetItem()
                w = self._obj_row(o, active)
                item.setSizeHint(w.sizeHint())
                self.obj_list.addItem(item)
                self.obj_list.setItemWidget(item, w)
            self._fit_obj_list()
            self.update_traj_section(active)
        obj_label_refresh = refresh_objs

        def _obj_row(self, name, active):
            w = QtWidgets.QWidget(); h = QtWidgets.QHBoxLayout(w); h.setContentsMargins(2, 0, 2, 0); h.setSpacing(4)
            vis = QtWidgets.QCheckBox(); vis.setChecked(name in _enabled_set()); vis.setToolTip("Show / hide")
            vis.clicked.connect(lambda checked, n=name: self.toggle_visible(n, checked))
            h.addWidget(vis)
            btn = QtWidgets.QPushButton(name); btn.setFlat(True)
            btn.setStyleSheet("text-align:left; border:none; padding:0 2px; background:transparent;"
                              + ("font-weight:bold; color:#1565c0;" if name == active else ""))
            btn.setToolTip("Click to make active")
            btn.clicked.connect(lambda _=False, n=name: self.set_active(n))
            h.addWidget(btn, 1)
            x = QtWidgets.QToolButton(); x.setText("✕"); x.setToolTip("Remove structure")
            x.setStyleSheet("QToolButton { padding:0 4px; border:none; background:transparent; color:#8a929b; }"
                            " QToolButton:hover { color:#c62828; }")
            x.clicked.connect(lambda _=False, n=name: self.remove_structure(n))
            h.addWidget(x)
            return w

        def set_active(self, name):
            _STATE["obj"] = name
            self.aff_combo.setVisible(bool(name) and _is_alphafill(name))
            self.update_colour_modes()
            ann = _CACHE.get(self.cur_uid()) or {}
            cov = actual_coverage(name, len(ann.get("sequence", "")))
            self.obj_label.setText(("%s · %.0f%% modelled" % (name, cov)) if cov is not None else name)
            self._sections["structures"].setSummary(self.obj_label.text())
            self.update_map_status(name)
            self.refresh_objs()
            self.sync_layers()

        def update_map_status(self, name=None):
            """One-line numbering summary of the active structure: mapping source, chimeric flag, partners."""
            name = name or self.obj()
            mp = UFV_MAPS.get(name) or {}
            mode = mp.get("mode")
            bits = []
            if mode in (None, "identity"):
                bits.append("Numbering: identity (resi == UniProt)")
            elif mode == "manual":
                bits.append("Numbering: manual (%s)" % ",".join(sorted(mp.get("chains") or {})))
            else:
                chains = sorted(mp.get("segments") or {}) if mode != "ranges" else []
                bits.append("Numbering: %s%s" % (mp.get("status") or mode,
                                                 (" · chains %s" % ",".join(chains)) if chains else ""))
            tip = []
            if mp.get("chimeric"):
                orgs = mp.get("organisms") or []
                bits.append("⚛ chimeric" + ((" (" + " / ".join(orgs) + ")") if len(orgs) > 1 else ""))
                tip.append(_chimeric_tip(mp))
            parts = {}
            for p in mp.get("partners") or []:
                parts.setdefault(p["accession"], set()).add(str(p.get("chainId")))
            for c, a in (mp.get("partnersCustom") or {}).items():
                parts.setdefault(a, set()).add(c)
            if parts:
                bits.append("partners: " + ", ".join("%s (%s)" % (a, ",".join(sorted(c))) for a, c in
                                                     list(parts.items())[:PARTNER_ACCESSION_CAP]))
            if mode == "custom":
                tip += ["%s: %s" % (c or "''", n) for c, n in sorted((mp.get("notes") or {}).items())]
            self.map_label.setText(" · ".join(bits) if name else "")
            self.map_label.setToolTip("\n".join(t for t in tip if t))
            self.chim_btn.setVisible(bool(name) and mode in ("sifts", "custom", "ranges", "manual"))

        def show_chimera(self):
            o = self.obj()
            if not o:
                return
            with _no_undo():
                n = ufv_chimera(o)
            self.status.setText("%d unmapped residue(s) highlighted on %s." % (n, o) if n else
                                "Every modelled residue of %s maps to UniProt." % o)

        def toggle_visible(self, name, checked):
            with _no_undo():
                try:
                    if checked:
                        cmd.enable(name)
                        # re-show only the overlays that are switched on (hidden layer objects persist
                        # now, so blindly enabling every ufv_<name>_* would resurrect hidden layers)
                        for tag in (_SPHERE_STATE.get(name) or {}):
                            cmd.enable(_ufv_layer_obj(name, tag))
                        if name in _LIG_SHOWN and not (_FOCUS_SHOWN.get(name) or {}).get("lig"):
                            for lobj in _lig_objs(name):
                                cmd.enable(lobj)
                    else:
                        cmd.disable(name)
                        names = _object_set()
                        for tag in _LAYER_TAGS:
                            lobj = _ufv_layer_obj(name, tag)
                            if lobj in names:
                                cmd.disable(lobj)
                except Exception:
                    pass

        def remove_structure(self, name):
            if not name:
                return
            with _no_undo():
                ufv_clear(name)
                try:
                    cmd.delete(name)
                except Exception:
                    pass
            (_STATE.get("sources") or {}).pop(name, None)
            (_STATE.get("obj_uid") or {}).pop(name, None)
            pdb_objs = _STATE.get("pdb_objs") or {}
            for k in [k for k, v in pdb_objs.items() if v == name]:
                pdb_objs.pop(k, None)
            _invalidate_obj(name); UFV_MAPS.pop(name, None)
            if _STATE.get("obj") == name:
                _STATE["obj"] = (self.loaded_objs() or [None])[0]
            self.refresh_objs()

        # ---- fetch / structures ----
        def on_fetch(self):
            uid = self.cur_uid()
            if not uid:
                self.status.setText("Enter an accession."); return

            def work():
                # annotations and the structure list are independent: fetch them concurrently (the
                # shared UniProt entry JSON is downloaded once and both are disk-cached).
                from concurrent.futures import ThreadPoolExecutor
                with ThreadPoolExecutor(max_workers=2) as ex:
                    fa = ex.submit(fetch_annotations, uid)
                    fl = ex.submit(list_structures, uid)
                    return fa.result(), fl.result()
            self._async(work, lambda res: self.after_fetch(res[1], uid), "Fetching %s ..." % uid, self.fetch_btn)

        def after_fetch(self, structs=None, uid=None):
            structs = structs if structs is not None else (_STATE.get("structures") or [])
            _STATE["structures"] = structs
            _STATE["structures_uid"] = (uid or self.cur_uid() or "").upper()
            if uid:
                _STATE["uid"] = uid.upper()
            self.struct_combo.clear()
            # item data = index into _STATE["structures"] (a dict would be deep-converted by Qt on every read)
            for i, s in enumerate(structs):
                self.struct_combo.addItem(s["label"], i)
                tip = "\n".join(t for t in (s.get("mappingStatus"), _chimeric_tip(s)) if t)
                if tip:
                    self.struct_combo.setItemData(i, tip, QtCore.Qt.ToolTipRole)
            self.update_info(); self.obj_label_refresh()
            self.ann_reload(force=True)
            note = structure_notice(_STATE.get("structures_uid"))
            if note:                             # e.g. 3D-Beacons down: say so instead of a silently short list
                self.status.setText(note)
            self._fetch_tractability(_STATE.get("uid") or self.cur_uid())
            if self.cb_protnlm.isChecked():
                self._render_protnlm((_CACHE.get(self.cur_uid() or "") or {}).get("protnlm"))

        def cur_struct(self):
            i = self.struct_combo.currentData()
            if isinstance(i, dict):              # item added with the struct itself as data
                return i
            structs = _STATE.get("structures") or []
            return structs[i] if isinstance(i, int) and 0 <= i < len(structs) else None

        def on_load_selected(self):
            s = self.cur_struct()
            if not s:
                self.status.setText("Fetch, then choose a structure."); return
            uid = self.cur_uid()
            # A PDB entry is loaded whole (all chains, SIFTS-mapped), so another chain row of an entry
            # that is already loaded just re-selects that object instead of downloading it again.
            if s.get("pdbId"):
                prev = (_STATE.get("pdb_objs") or {}).get((s["pdbId"].upper(), (uid or "").upper()))
                if prev and prev in _object_set():
                    self.set_active(prev)
                    self.status.setText("%s is already loaded as %s." % (s["pdbId"], prev)); return
            self._async(lambda: _prepare_structure(s, uid),
                        lambda prep: self._finish_load(s, prep, uid), "Downloading %s ..." % s["label"],
                        self.load_btn)

        def _finish_load(self, struct, prep, uid):
            if not prep or not prep.get("path"):
                self.status.setText("Download failed."); return
            name = _load_from_path(struct, prep["path"], uid, prep.get("segments"), prep.get("hints"),
                                   meta=prep.get("meta"))
            self.cart.setCurrentIndex(0)
            self.set_active(name)

        def on_load_all(self):
            structs = _STATE.get("structures") or []
            if not structs:
                self.status.setText("Fetch, then load."); return
            uid = self.cur_uid()
            # One load per PDB ENTRY, not per (entry, chain) row: each row downloads and loads the
            # whole entry, so a receptor listed under 5 chains used to be loaded 5 times.
            todo, seen = [], set()
            for s in structs:
                key = ("PDB", s["pdbId"].upper()) if s.get("pdbId") else ("URL", s.get("url"))
                if key in seen:
                    continue
                seen.add(key)
                prev = (_STATE.get("pdb_objs") or {}).get((s["pdbId"].upper(), (uid or "").upper())) \
                    if s.get("pdbId") else None
                if prev and prev in _object_set():
                    continue
                todo.append(s)
            if not todo:
                self.status.setText("Everything listed is already loaded."); return

            def work():
                # downloads (+ SIFTS) in parallel, off the UI thread
                from concurrent.futures import ThreadPoolExecutor

                def prep(s):
                    try:
                        return s, _prepare_structure(s, uid)
                    except Exception:
                        return s, None
                with ThreadPoolExecutor(max_workers=6) as ex:
                    return list(ex.map(prep, todo))

            def apply(results):
                # Load one structure per event-loop tick so the viewer and panel stay responsive
                # (and repaint) between loads; orient once at the end instead of after every load.
                queue = [(s, p) for s, p in results if p and p.get("path")]
                loaded, errors = [], []
                was_empty = not self.loaded_objs()

                def step():
                    if queue:
                        s, prep = queue.pop(0)
                        try:
                            with _no_undo():
                                loaded.append(_load_from_path(s, prep["path"], uid, prep.get("segments"),
                                                              prep.get("hints"), orient=False,
                                                              meta=prep.get("meta")))
                        except Exception as exc:
                            errors.append("%s: %s" % (s.get("key"), exc))
                        self.status.setText("Loading %d / %d ..." % (len(loaded), len(results)))
                        QtCore.QTimer.singleShot(0, step)
                        return
                    self.loadall_btn.setEnabled(True)
                    self.cart.setCurrentIndex(0)
                    if loaded:
                        if was_empty:           # only the first load into an empty session frames it
                            with _no_undo():
                                cmd.orient(loaded[0], animate=0)
                                _view_commit()
                        self.set_active(loaded[0])
                    self.refresh_objs()
                    self.status.setText("Loaded %d of %d structures (%d listed rows).%s"
                                        % (len(loaded), len(results), len(structs),
                                           (" First error: " + errors[0]) if errors else ""))
                self.loadall_btn.setEnabled(False)
                step()
            self._async(work, apply, "Downloading %d structures ..." % len(todo), self.loadall_btn)

        def update_colour_modes(self):
            predicted = _structure_is_predicted(self.obj())
            model = self.cart.model()
            for i in range(self.cart.count()):
                txt = self.cart.itemText(i)
                ok = True
                if txt == "pLDDT":
                    ok = predicted
                elif txt == "B-factor":
                    ok = not predicted
                model.item(i).setEnabled(ok)

        # ---- layers ----
        def on_aff_identity(self, i):
            AFF_SETTINGS["min_identity"] = AFF_IDENTITY_STEPS[i] / 100.0
            with _no_undo():
                for o in self.target_objs():
                    r = _apply_aff_filter(o)
                    if r:
                        self.status.setText("AlphaFill: %d/%d transplants with donor identity ≥ %d%%."
                                            % (r[0], r[1], AFF_IDENTITY_STEPS[i]))

        def _fetch_tractability(self, uid):
            self.tract_label.setVisible(False); self.tract_link.setVisible(False)
            if not uid:
                return

            def done(ot):
                if isinstance(ot, Exception) or not ot or not ot.get("groups") or (self.cur_uid() or "") != uid:
                    return
                # compact one-liner (modality + bucket count); the full bucket list is the tooltip
                self.tract_link.setText('<a href="%s">%s ↗</a>' % (_h(ot["url"]), _h(ot.get("symbol") or "Open Targets")))
                self.tract_link.setToolTip("Open Targets target page")
                self.tract_label.setText("Tractability — %s" % " · ".join(
                    "%s (%d)" % (g["modality"], len(g["labels"])) for g in ot["groups"]))
                self.tract_label.setToolTip(
                    "Open Targets tractability — how amenable this whole protein is to a small molecule, antibody or "
                    "PROTAC (positive evidence only):\n\n" +
                    "\n".join("%s: %s" % (g["modality"], ", ".join(g["labels"])) for g in ot["groups"]))
                self.tract_label.setVisible(True); self.tract_link.setVisible(True)
            self._bg(lambda: fetch_open_targets(uid), done)

        def on_nearby(self, v):
            NEARBY_SETTINGS["dist"] = float(v)
            _refocus_all()
            if self._lig_resn:
                self._render_ligand()
            elif self._last_rep and self._last_rep.get("obj"):
                try:
                    self._last_rep["nearby"] = residue_nearby(self._last_rep["obj"], self._last_rep["pos"])
                    self._last_rep["varColor"] = residue_report_colors(self._last_rep)
                except Exception:
                    pass
                self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence, self._last_protvar or "off"))

        def on_all_pockets(self):
            o, uid = self.obj(), self.cur_uid()
            ann = _CACHE.get(uid or "")
            if not o or not ann:
                self.status.setText("Fetch + load a structure first."); return
            n = len(ann.get("sequence") or "")
            cached = _ALL_POCKETS.get((_canon_acc(uid), n), False)
            if cached is not False:
                self._apply_all_pockets(o, uid, cached); return
            self.allpk_btn.setEnabled(False)
            self.status.setText("Scanning all residues for ProtVar pockets … 0%")
            wk = _PocketWorker(uid, n)

            def on_done(res):
                self.allpk_btn.setEnabled(True); self._workers.discard(wk)
                if isinstance(res, Exception):
                    self.status.setText("Pocket scan error: %s" % res); return
                self.status.setText("")
                self._apply_all_pockets(o, uid, res)
            wk.progress.connect(lambda d, t: self.status.setText(
                "Scanning all residues for ProtVar pockets … %d%%" % round(100.0 * d / max(t, 1))))
            wk.done.connect(on_done); self._workers.add(wk); wk.start()

        def _apply_all_pockets(self, o, uid, res):
            self._all_pockets = res
            if not res:
                self.status.setText("No ProtVar pockets detected across this protein."); return
            with _no_undo():
                _show_all_pockets(o, res["pockets"])
            self._lig_resn = None; self._last_rep = None
            self.detail.setHtml(all_pockets_html(res, o, uid))
            self.status.setText("%d ProtVar pocket(s) coloured." % res["count"])

        def on_export_bfactor(self):
            o, uid = self.obj(), self.cur_uid()
            if not o or not uid or uid not in _CACHE:
                self.status.setText("Fetch + load a structure first."); return
            mode = {"AlphaMissense": "alphamissense", "Hotspots": "hotspots", "Contact hubs": "hubs",
                    "Burden": "burden", "Constraint pocket": "pocket", "pLDDT": "plddt",
                    "B-factor": "bfactor"}.get(self.cart.currentText(), "selected")
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save B-factor-encoded PDB", "%s_%s_bfactor.pdb" % (o, mode), "PDB (*.pdb)")
            if not path:
                return
            try:
                ufv_save_bfactor(path, o, mode, uid)
                self.status.setText("Saved %s (B-factor = %s)." % (os.path.basename(path), mode))
            except Exception as exc:
                self.status.setText("B-factor export error: %s" % exc)

        def on_pocket_q(self, val):
            _set_pocket_q(val)
            if self.cart.currentText() != "Constraint pocket":
                return
            uid = self.cur_uid()
            with _no_undo():
                for o in self.target_objs():
                    res = _analysis_cached(o, uid, "pockets")
                    if res is not None:
                        _set_cartoon_layer(o, "constraint pocket", grey=False)
                        n = _color_tiers(o, _pocket_classes(res), POCKET_COLORS, "gray80")
                        self.status.setText("Constraint pocket: %d residue(s) at q ≤ %.2f." % (n, POCKET_SETTINGS["q"]))

        def on_prox_radius(self, _v=None):
            ufv_prox(self.prox_spin["ptm"].value(), self.prox_spin["var"].value())
            if self._last_pos is not None and (self._sel or {}).get("kind") == "res":
                objs, uid = self.target_objs(), self.cur_uid()
                if objs and uid:
                    try:
                        self._last_rep = residue_report(objs[0], uid, self._last_pos)
                        self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence,
                                                               self._last_protvar or "off"))
                    except Exception:
                        pass

        def apply_cartoon(self):
            objs = self.target_objs(); choice = self.cart.currentText()
            self.q_label.setVisible(choice == "Constraint pocket"); self.q_spin.setVisible(choice == "Constraint pocket")
            if not objs:
                return
            if choice == "None":
                with _no_undo():
                    for o in objs:
                        _reset_cartoon(o)
                return
            uid = self.cur_uid()
            if choice not in ("pLDDT", "B-factor") and not _CACHE.get(uid):
                self.status.setText("Fetch first."); self.cart.setCurrentIndex(0); return
            cheap = {"Domains": ufv_domains, "Topology": ufv_topology, "pLDDT": ufv_plddt,
                     "B-factor": ufv_bfactor, "AlphaMissense": ufv_alphamissense, "Burden": ufv_burden}
            if choice in cheap:
                with _no_undo():
                    for o in objs:
                        cheap[choice](o, uid)
                return
            ann = _CACHE.get(uid) or {}
            if choice == "Constraint pocket" and not ann.get("amMap"):
                self.status.setText("Constraint pocket needs AlphaMissense data."); self.cart.setCurrentIndex(0); return
            kind, colors = {"Hotspots": ("hotspots", HOTSPOT_COLORS),
                            "Contact hubs": ("hubs", HUB_COLORS)}.get(choice, ("pockets", POCKET_COLORS))

            def apply(res):
                if self.cart.currentText() != choice:
                    return                               # user picked another mode meanwhile
                with _no_undo():
                    for o, r in res.items():
                        _set_cartoon_layer(o, choice.lower(), grey=False)
                        _color_tiers(o, _tiers_for(kind, r), colors, "gray80")
                extra = ""
                if kind == "pockets":
                    r0 = next(iter(res.values()), None) or {}
                    extra = " (q ≤ %.2f%s%s)" % (POCKET_SETTINGS["q"], ", PAE-gated" if r0.get("hasPae") else "",
                                                 ("; " + r0["reason"]) if r0.get("reason") else "")
                self.status.setText("%s applied to %d structure(s)%s." % (choice, len(res), extra))
            cached = {o: _analysis_cached(o, uid, kind) for o in objs}
            if all(v is not None for v in cached.values()):
                apply(cached)                            # memoised: instant, no worker
                return
            geoms = {o: _ca_geometry(o) for o in objs}    # read coords on the UI thread
            work = lambda: {o: _analysis(o, uid, kind, g, ann) for o, g in geoms.items()}
            self._async(work, apply, "Computing %s ..." % choice, self.cart)

        # ---- report ----
        def report_residue(self, uni_pos, focus=False, objs=None, zoom=True):
            """Report a residue; focus=True also shows its focus sticks, and zoom=True frames them (explicit
            view action). zoom=False is the list's single-click highlight: sticks only, layers stay shown,
            the camera does not move."""
            objs, uid = objs or self.target_objs(), self.cur_uid()
            if not objs or not uid:
                return
            self._last_pos = uni_pos
            self._sel = {"kind": "res", "pos": uni_pos}
            if focus:  # zoom FIRST so the view moves before any compute
                fsets = []
                with _no_undo():
                    for o in objs:
                        try:
                            f = ufv_focus(o, uni_pos, uid, zoom=0, declutter=1 if zoom else 0)
                            if f:
                                fsets.append(f)
                        except Exception as exc:
                            self.status.setText("Focus error: %s" % exc)
                    if fsets and zoom:          # frame every focused copy at once (one zoom, not N)
                        _focus_zoom(" or ".join(fsets), 3)
            # Build the report on the NEXT event-loop tick so the zoom repaints immediately instead of
            # waiting for residue_report (neighbourhood + ligand-proximity can be heavy on big models).
            self._lig_resn = None
            QtCore.QTimer.singleShot(0, lambda p=uni_pos, o=objs[0], u=uid: self._build_report(o, u, p))

        def _build_report(self, obj, uid, pos):
            if self._last_pos != pos:
                return                      # user already clicked elsewhere
            try:
                self._last_rep = residue_report(obj, uid, pos)
                self._show_evidence = False
                self._last_protvar = "loading"   # ProtVar always shown; fetched in the background
                self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence, "loading"))
                self._fetch_protvar_async(uid, pos)
            except Exception as exc:
                self._last_rep = None
                self.detail.setHtml("<i>Report error: %s</i>" % exc)

        def _fetch_protvar_async(self, uid, pos):
            pos = int(pos)

            def work():
                from concurrent.futures import ThreadPoolExecutor
                with ThreadPoolExecutor(max_workers=1) as ex:
                    fk = ex.submit(fetch_pdbekb, uid)
                    pv = fetch_protvar(uid, pos)
                    try:
                        kb = fk.result()
                    except Exception:
                        kb = None
                return pv, kb

            def done(res):
                if self._last_pos != pos:
                    return                                   # stale — user already clicked elsewhere
                pv, kb = res if not isinstance(res, Exception) else (None, None)
                self._last_protvar = pv
                if self._last_rep:
                    self._last_rep["pdbekb"] = kb
                    self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence, self._last_protvar))
            self._bg(work, done)

        def _bg(self, work, done):
            """Run `work` off the UI thread WITHOUT disabling the panel (unlike _async), for background
            enrichment like ProtVar that shouldn't freeze interaction."""
            wk = _Worker(work)

            def on_done(res):
                self._workers.discard(wk)
                done(res)
            wk.done.connect(on_done); self._workers.add(wk); wk.start()

        def on_anchor(self, url):
            """Handle in-report links: ufv:res:<n> (jump to residue), ufv:evi (toggle evidence),
            ufv:lig:prev / ufv:lig:next (circulate same-named ligand copies)."""
            s = url.toString() if hasattr(url, "toString") else str(url)
            if s.startswith("http://") or s.startswith("https://"):
                try:
                    import webbrowser
                    webbrowser.open(s)
                except Exception:
                    pass
                return
            if s.startswith("ufv:res:"):
                try:
                    self.report_residue(int(s.rsplit(":", 1)[1]), focus=True)
                except ValueError:
                    pass
            elif s == "ufv:pocket":
                pk = self._last_protvar.get("pocket") if isinstance(self._last_protvar, dict) else None
                if pk and pk.get("lining"):
                    ufv_show_pocket(self.obj(), pk["lining"])
            elif s.startswith("ufv:pocket:") or s.startswith("ufv:motif:"):
                o = self.obj()
                if s == "ufv:motif:lig":
                    f = _FOCUS_SHOWN.get(o) or {}
                    resids = f.get("near") or []
                else:
                    pk = self._last_protvar.get("pocket") if isinstance(self._last_protvar, dict) else None
                    pl = (pk or {}).get("pockets") or []
                    i = int(s.rsplit(":", 1)[1])
                    resids = pl[i]["resid"] if 0 <= i < len(pl) else []
                if not o or not resids:
                    return
                if s.startswith("ufv:pocket:"):
                    ufv_show_pocket(o, resids)
                else:
                    url = motif_search_url(o, resids, self.cur_uid())
                    if url:
                        _open_url(url); self.status.setText("RCSB Structure Motif Search opened in the browser.")
                    else:
                        self.status.setText("Motif search unavailable for this pocket / structure.")
            elif s == "ufv:allpockets":
                self.on_all_pockets()
            elif s.startswith("ufv:apk:"):
                pl = (self._all_pockets or {}).get("pockets") or []
                key = s.split(":", 2)[2]
                o = self.obj()
                if not o:
                    return
                if key == "overlay":
                    with _no_undo():
                        _show_all_pockets(o, pl)
                elif key.startswith("motif:"):
                    i = int(key.split(":")[1])
                    url = motif_search_url(o, pl[i]["resid"], self.cur_uid()) if 0 <= i < len(pl) else None
                    if url:
                        _open_url(url)
                else:
                    i = int(key)
                    if 0 <= i < len(pl):
                        ufv_show_pocket(o, pl[i]["resid"])
            elif s.startswith("ufv:lig:show:"):
                self.on_ligand(s.rsplit(":", 1)[1])    # jump to a Tanimoto-similar ligand
            elif s == "ufv:evi":
                self._show_evidence = not self._show_evidence
                if self._last_rep:
                    self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence,
                                                           self._last_protvar or "off"))
            elif s in ("ufv:lig:prev", "ufv:lig:next") and self._lig_copies:
                self._lig_idx = (self._lig_idx + (1 if s.endswith("next") else -1)) % len(self._lig_copies)
                ch, resi = self._lig_copies[self._lig_idx]
                o = self.obj()
                if o:
                    ufv_ligand_focus(o, self._lig_resn, ch, resi, self.cur_uid())
                self._render_ligand()

        def on_ligand(self, resn, pick=None, zoom=True):
            o = self.obj()
            if not o:
                return
            self._lig_resn = resn
            self._sel = {"kind": "lig", "resn": resn}
            self._lig_copies = ligand_instances(o, resn)
            self._lig_idx = 0
            if pick is not None and pick in self._lig_copies:
                self._lig_idx = self._lig_copies.index(pick)
            uid = self.cur_uid()
            if self._lig_copies:
                ch, resi = self._lig_copies[self._lig_idx]
                ufv_ligand_focus(o, resn, ch, resi, uid, zoom=1 if zoom else 0)
            else:
                ufv_ligand_focus(o, resn, uid=uid, zoom=1 if zoom else 0)
            codes = ligand_codes_ordered(o)          # PyMOL read on the UI thread; network in the worker
            if _is_alphafill(o) and o not in _AFF_META:
                def _aff_done(meta, o=o, resn=resn):
                    _AFF_META[o] = meta if isinstance(meta, dict) else {}
                    if self._lig_resn == resn:
                        self._render_ligand()
                self._bg(lambda: alphafill_transplants(uid), _aff_done)
            if self._lig_chem and all(c in self._lig_chem for c in codes):
                self._render_ligand(); return        # chemistry for this structure already fetched
            self.detail.setHtml("<b>%s</b> — fetching chemistry…" % resn)
            self.status.setText("Fetching ligand chemistry ...")
            self._bg(lambda: ligand_chemistry(codes=codes),
                     lambda chem: self.show_ligand(resn, chem))

        def show_ligand(self, resn, chem):
            if self.status.text() == "Fetching ligand chemistry ...":
                self.status.setText("")
            if isinstance(chem, Exception):
                self.status.setText("Ligand chemistry error: %s" % chem); return
            self._lig_chem = chem
            if self._lig_resn == resn:               # still showing this ligand
                self._render_ligand()
            if (self.ann_cat.currentData() or "all") in ("all", "lig"):
                self.ann_reload(force=True)          # fill in the names

        def _render_ligand(self):
            if not self._lig_resn:
                return
            info = (self._lig_chem or {}).get(self._lig_resn, {})
            o = self.obj()
            extra = {}
            if o:
                ch, resi = (self._lig_copies[self._lig_idx] if self._lig_copies else (None, None))
                try:
                    extra = ligand_panel_data(o, self._lig_resn, ch, resi, self.cur_uid())
                except Exception:
                    extra = {}
            self.detail.setHtml(format_ligand_html(self._lig_resn, info,
                                                   instances=len(self._lig_copies) or 1,
                                                   copies=self._lig_copies, cur=self._lig_idx, **extra))

        # ---- 3D picking ----
        # Two capture paths feed one handler so picking works regardless of PyMOL build:
        #   1. a native Wizard (do_select / do_pick) — fires the instant you click an atom;
        #   2. a QTimer that watches the "sele" PyMOL creates on a click — a fallback if (1) is quiet.
        # _handle_pick is debounced so the two paths never double-act, ALWAYS zooms (it does not depend
        # on an accession being fetched), and writes a status line so it's obvious the click registered.
        def toggle_pick(self):
            with _no_undo():
                self._toggle_pick()

        def _toggle_pick(self):
            if self.cb_pick.isChecked():
                self._install_picker()
                t = getattr(self, "_pick_timer", None)
                if t is None:
                    t = QtCore.QTimer(self); t.setInterval(200)
                    t.timeout.connect(self._poll_pick); self._pick_timer = t
                try:
                    cmd.deselect()
                except Exception:
                    pass
                t.start()
                self.status.setText("Pick 3D on — click any atom (protein or ligand) to zoom in.")
            else:
                self._remove_picker()
                t = getattr(self, "_pick_timer", None)
                if t is not None:
                    t.stop()
                self.status.setText("")

        def _install_picker(self):
            try:
                from pymol.wizard import Wizard
            except Exception:
                self._picker = None
                return
            panel = self

            class _Picker(Wizard):
                def get_prompt(self):
                    return ["Click any atom to zoom to its residue / ligand."]

                def get_panel(self):
                    # Panel command strings (3rd field) are run through PyMOL's command parser. An
                    # EMPTY command makes some PyMOL builds dump the whole command list to the console
                    # ("parser: matching commands: abort accept …"), so the banner uses a harmless,
                    # valid no-op call instead of "".
                    return [[1, "UFV pick", "cmd.get_wizard()"],
                            [2, "Stop picking", "cmd.set_wizard()"]]

                def _grab(self, selection):
                    atoms = []
                    try:
                        m = self.cmd.get_model(selection)
                        atoms = list(getattr(m, "atom", []))
                    except Exception:
                        pass
                    try:
                        self.cmd.deselect()
                    except Exception:
                        pass
                    if atoms:
                        a = atoms[0]
                        panel._handle_pick(a.chain, a.resi, a.resn, int(getattr(a, "hetatm", 0) or 0))

                def do_select(self, selection):
                    self._grab(selection)
                    try:
                        self.cmd.refresh_wizard()
                    except Exception:
                        pass
                    return None

                def do_pick(self, bondFlag):
                    self._grab("pk1")
                    try:
                        self.cmd.unpick()
                    except Exception:
                        pass
                    return None

            try:
                self._picker = _Picker()
                cmd.set_wizard(self._picker)
            except Exception:
                self._picker = None

        def _remove_picker(self):
            if getattr(self, "_picker", None) is not None:
                try:
                    cmd.set_wizard()
                except Exception:
                    pass
                self._picker = None

        def _poll_pick(self):
            try:
                # name lookup first (~0.05 ms): count_atoms evaluates a selection over EVERY atom of the
                # session (several ms per 200 ms tick in a large session) and raises when there is no sele
                if "sele" not in (cmd.get_names("selections") or ()) or cmd.count_atoms("sele") == 0:
                    return
            except Exception:
                return
            info = []
            try:
                with _no_undo():
                    cmd.iterate("sele", "info.append((chain, resi, resn, hetatm, model))", space={"info": info})
                    cmd.deselect()      # clear the click so it isn't re-read and the pink dots disappear
            except Exception:
                return
            if info:
                ch, resi, resn, het, model = info[0]
                self._handle_pick(ch, resi, resn, het, model)

        def _handle_pick(self, ch, resi, resn, het, model=None):
            import time
            key = (ch, resi, resn, het)
            now = time.time()
            if key == getattr(self, "_last_pick", None) and now - getattr(self, "_last_pick_t", 0) < 0.6:
                return                                  # debounce: the two capture paths firing together
            self._last_pick, self._last_pick_t = key, now
            o = self.obj()
            if model and not str(model).startswith("ufv_") and model in _object_set():
                o = model                               # a click on another loaded structure maps through IT
            if not o:
                self.status.setText("Pick: load a structure first."); return
            with _no_undo():
                if het:
                    self.on_ligand(resn, pick=(ch, resi))
                    self.status.setText("Picked ligand %s (%s/%s)" % (resn, ch, resi)); return
                uni = _resi_to_uni(o, ch, resi)
                if uni is None:
                    self.status.setText("Picked %s/%s — not mapped to a UniProt position." % (ch, resi)); return
                self._sel = {"kind": "res", "pos": uni}; self._last_pos = uni
                objs = self.target_objs() or [o]
                if o not in objs:
                    objs = [o]
                if self.cur_uid():
                    # same path as a row click: focus + deferred report + background ProtVar
                    self.report_residue(uni, focus=True, objs=objs)
                else:
                    try:                                # ALWAYS zoom — independent of fetched annotations
                        for ob in objs:
                            ufv_focus(ob, uni)
                    except Exception as exc:
                        self.status.setText("Focus error: %s" % exc); return
            self.status.setText("Picked residue %s%d" % (resn[:1] if resn else "", uni))

        # ---- numbering dialog ----
        def open_numbering_dialog(self):
            dlg = self.numbering_dialog()
            if dlg is not None:
                dlg.exec_()

        def numbering_dialog(self):
            """Residue-numbering dialog: identity / SIFTS / one manual chain, or the per-chain table of the
            extension's Custom PDB (accession · auto|offset|adopt|ignore · value) for every chain."""
            o, uid = self.obj(), (self.cur_uid() or "").upper()
            if not o or not uid:
                self.status.setText("Load/select an object first."); return None
            mp = UFV_MAPS.get(o) or {"mode": "identity"}
            seqlen = len((_CACHE.get(uid) or {}).get("sequence") or "")
            dlg = QtWidgets.QDialog(self); dlg.setWindowTitle("Residue numbering — %s" % o)
            dlg.resize(640, 420)
            v = QtWidgets.QVBoxLayout(dlg)
            g = QtWidgets.QGridLayout(); v.addLayout(g)
            rb_id = QtWidgets.QRadioButton("Identity (resi == UniProt)")
            rb_si = QtWidgets.QRadioButton("SIFTS, PDB:"); pdb = QtWidgets.QLineEdit(mp.get("pdbId") or "")
            pdb.setMaximumWidth(80)
            rb_man = QtWidgets.QRadioButton("Manual chain:")
            ch = QtWidgets.QLineEdit("A"); ch.setMaximumWidth(34)
            us = QtWidgets.QLineEdit(); us.setPlaceholderText("UniProt@"); us.setMaximumWidth(72)
            rs = QtWidgets.QLineEdit(); rs.setPlaceholderText("resi@"); rs.setMaximumWidth(58)
            re_ = QtWidgets.QLineEdit(); re_.setPlaceholderText("end"); re_.setMaximumWidth(52)
            rb_cu = QtWidgets.QRadioButton("Per chain (Custom PDB) — accession · mapping mode · value:")
            g.addWidget(rb_id, 0, 0, 1, 5)
            g.addWidget(rb_si, 1, 0); g.addWidget(pdb, 1, 1)
            g.addWidget(rb_man, 2, 0); g.addWidget(ch, 2, 1); g.addWidget(us, 2, 2); g.addWidget(rs, 2, 3); g.addWidget(re_, 2, 4)
            g.addWidget(rb_cu, 3, 0, 1, 5)
            {"sifts": rb_si, "manual": rb_man, "custom": rb_cu}.get(mp.get("mode"), rb_id).setChecked(True)

            ext = _chain_extents(o)
            cfgs = (mp.get("config") or {}) if mp.get("mode") == "custom" else {}
            sources = []
            for x in self.loaded_objs():
                if x == o:
                    continue
                sources.append(x)
                sm = UFV_MAPS.get(x) or {}
                sources += ["%s/%s" % (x, c) for c in sorted(sm.get("segments") or {})]
            tbl = QtWidgets.QTableWidget(len(ext), 6)
            tbl.setHorizontalHeaderLabels(["Chain", "Residues", "Accession", "Mode", "Value", "Current mapping"])
            tbl.verticalHeader().setVisible(False)
            tbl.horizontalHeader().setStretchLastSection(True)
            rows = []
            for i, c in enumerate(sorted(ext)):
                e = ext[c]
                cfg = cfgs.get(c) or {"mode": "auto"}
                it = QtWidgets.QTableWidgetItem(c or "''"); it.setFlags(QtCore.Qt.ItemIsEnabled)
                tbl.setItem(i, 0, it)
                it = QtWidgets.QTableWidgetItem("%d–%d (%d)" % (e[0], e[1], e[2])); it.setFlags(QtCore.Qt.ItemIsEnabled)
                tbl.setItem(i, 1, it)
                acc = QtWidgets.QLineEdit(cfg.get("accession") or uid)
                acc.setToolTip("UniProt accession this chain represents (another accession = partner chain)")
                mode = QtWidgets.QComboBox(); mode.addItems(list(CUSTOM_MODES))
                mode.setCurrentIndex(CUSTOM_MODES.index(cfg.get("mode", "auto")))
                val = QtWidgets.QComboBox(); val.setEditable(True); val.addItem(""); val.addItems(sources)
                val.lineEdit().setPlaceholderText("UniProt@%d  or  object[/chain]" % e[0])
                val.setToolTip("offset: UniProt position of residue %d · adopt: loaded structure to copy" % e[0])
                if cfg.get("mode") == "offset" and cfg.get("start") is not None:
                    val.setEditText(str(cfg["start"]))
                elif cfg.get("mode") == "adopt" and cfg.get("src"):
                    val.setEditText(cfg["src"] + ("/" + cfg["src_chain"] if cfg.get("src_chain") else ""))
                cur = _map_chain_segments(mp, c, seqlen) if mp.get("mode") not in (None, "identity") else None
                note = (mp.get("notes") or {}).get(c, "")
                it = QtWidgets.QTableWidgetItem(_seg_text(cur) if cur else (note or ("identity" if mp.get("mode") in (None, "identity") else "unmapped")))
                it.setFlags(QtCore.Qt.ItemIsEnabled); it.setToolTip(note)
                tbl.setItem(i, 5, it)
                tbl.setCellWidget(i, 2, acc); tbl.setCellWidget(i, 3, mode); tbl.setCellWidget(i, 4, val)

                def on_val(text, mode=mode):
                    t = text.strip()
                    if re.match(r"^-?\d+$", t):
                        mode.setCurrentIndex(CUSTOM_MODES.index("offset"))
                    elif t and t.split("/")[0] in _object_set():
                        mode.setCurrentIndex(CUSTOM_MODES.index("adopt"))
                    rb_cu.setChecked(True)
                val.editTextChanged.connect(on_val)
                mode.currentIndexChanged.connect(lambda _i: rb_cu.setChecked(True))
                acc.textEdited.connect(lambda _t: rb_cu.setChecked(True))
                rows.append((c, acc, mode, val))
            tbl.resizeColumnsToContents()
            v.addWidget(tbl, 1)
            err = QtWidgets.QLabel(""); err.setStyleSheet("color:#a33;"); v.addWidget(err)
            bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
            v.addWidget(bb)

            def custom_configs():
                out = {}
                for c, acc, mode, val in rows:
                    m = mode.currentText()
                    cfg = {"mode": m, "accession": acc.text().strip().upper() or uid}
                    t = val.currentText().strip()
                    if m == "offset":
                        if not re.match(r"^-?\d+$", t):
                            raise ValueError("Chain %s: offset needs the UniProt position of residue %d." % (c, ext[c][0]))
                        cfg["start"] = int(t)
                    elif m == "adopt":
                        src, _, sc = t.partition("/")
                        if src not in _object_set():
                            raise ValueError("Chain %s: choose a loaded structure to adopt the mapping from." % c)
                        cfg["src"], cfg["src_chain"] = src, sc or None
                    out[c] = cfg
                return out

            def apply_and_close():
                try:
                    with _no_undo():
                        if rb_si.isChecked():
                            ufv_map(o, uid, "sifts", pdb.text().strip())
                        elif rb_man.isChecked():
                            ufv_chain(o, ch.text().strip(), us.text().strip(), rs.text().strip() or None,
                                      re_.text().strip() or None)
                        elif rb_cu.isChecked():
                            _apply_custom_map(o, uid, custom_configs(), reset=True)
                        else:
                            ufv_map(o, uid, "identity")
                except ValueError as exc:
                    err.setText(str(exc)); return
                _invalidate_obj(o)
                self.set_active(o)
                dlg.accept()
            bb.accepted.connect(apply_and_close); bb.rejected.connect(dlg.reject)
            dlg.ufv_rows, dlg.ufv_apply, dlg.ufv_radios = rows, apply_and_close, (rb_id, rb_si, rb_man, rb_cu)
            dlg.ufv_table = tbl
            _ufv_wheel_guard(dlg)
            return dlg

        def on_clear(self):
            ufv_clear(self.obj())
            _LIG_SHOWN.discard(self.obj())
            self.cart.setCurrentIndex(0); self.detail.setHtml("")
            self._ann_sync()

        def on_export_csv(self, with_protvar=False):
            uid = self.cur_uid()
            if not uid or uid not in _CACHE:
                self.status.setText("Fetch an accession first."); return
            suffix = "_protvar" if with_protvar else ""
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Export residue annotation CSV", "%s_residue_annotations%s.csv" % (uid, suffix), "CSV (*.csv)")
            if not path:
                return
            for b in (self.csv_btn, self.csv_pv_btn):
                b.setEnabled(False)
            self.status.setText("Building CSV …" + (" (fetching ProtVar — this is slow)" if with_protvar else ""))
            o = self.obj()
            # PyMOL reads (Cα geometry, mapped chains) happen here on the UI thread; the worker only
            # computes (and reuses cached hotspot / hub / pocket tiers).
            geom = _ca_geometry(o) if o else None
            chains = _mapped_chains(o) if o else None
            ligc = _ligand_contacts(o) if o else None
            wk = _CsvWorker(_CACHE[uid], o, with_protvar, geom, chains, ligc)

            def on_done(res):
                for b in (self.csv_btn, self.csv_pv_btn):
                    b.setEnabled(True)
                self._workers.discard(wk)
                if isinstance(res, Exception):
                    self.status.setText("CSV error: %s" % res); return
                self._write_csv(path, res)
            wk.progress.connect(lambda d, t: self.status.setText("ProtVar: %d / %d residues …" % (d, t)))
            wk.done.connect(on_done); self._workers.add(wk); wk.start()

        def _write_csv(self, path, text):
            try:
                with open(path, "w") as fh:
                    fh.write(text or "")
                self.status.setText("Saved %s" % path)
            except Exception as exc:
                self.status.setText("Save error: %s" % exc)

        def _save_geometry(self):
            if self._settings is not None:
                try:
                    self._settings.setValue("window/geometry", self.saveGeometry())
                except Exception:
                    pass

        def closeEvent(self, ev):
            self._save_geometry()
            t = getattr(self, "_pick_timer", None)
            if t is not None:
                t.stop()
            self._remove_picker()
            super(_UFVPanel, self).closeEvent(ev)

        def hideEvent(self, ev):
            self._save_geometry()
            t = getattr(self, "_pick_timer", None)   # don't poll PyMOL while the panel is hidden
            if t is not None:
                t.stop()
            super(_UFVPanel, self).hideEvent(ev)

        def showEvent(self, ev):
            t = getattr(self, "_pick_timer", None)
            if t is not None and self.cb_pick.isChecked():
                t.start()
            self.update_traj_section()
            super(_UFVPanel, self).showEvent(ev)

    return _UFVPanel


@_undo_paused
def _load_from_path(struct, path, uid, segments=None, hints=None, name=None, orient=True, meta=None):
    """Main-thread half of loading: load `path` as a fresh object, grey cartoon, residue mapping."""
    name = name or _struct_obj_name(struct)
    uid = (uid or "").upper()
    first = not _loaded_structures()     # empty session: the one load that may frame the new structure
    if name in _object_set():
        # Re-loading the same structure used to APPEND a second state to the existing object;
        # replace it instead (and drop its overlays, which belong to the old coordinates).
        _clear_focus(name, restore_spheres=False)
        _delete_layers(name)
        _CARTOON_LAYER.pop(name, None)
        cmd.delete(name)
    _ANN_LIGS.pop(name, None)
    _LIG_ROWS.pop(name, None)
    _ann_forget(name)
    with _auto_zoom_off():
        cmd.load(path, name, zoom=0)
    try:
        os.remove(path)
    except OSError:
        pass
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.color("gray80", "(%s) and polymer" % name)
    if struct.get("numbering") == "sifts" and struct.get("pdbId"):
        _set_sifts_map(name, struct["pdbId"], uid, segs=segments, hints=hints, meta=meta)
        mp = UFV_MAPS.get(name)
        if mp is not None and struct.get("chainId") in (mp.get("segments") or {}):
            mp["chainId"] = struct["chainId"]           # primary (highest-coverage) chain
        _STATE.setdefault("pdb_objs", {})[(struct["pdbId"].upper(), uid)] = name
    elif struct.get("numbering") == "ranges" and struct.get("ranges"):
        # isoform model: canonical UniProt positions -> isoform residues through the VSP-derived segments
        UFV_MAPS[name] = {"mode": "ranges", "ranges": [dict(r) for r in struct["ranges"]], "uid": uid,
                          "isoform": struct.get("isoform"), "status": struct.get("mappingStatus")}
        print("[UFV] %s: isoform %s — canonical annotations mapped across %d segment(s)."
              % (name, struct.get("isoform"), len(struct["ranges"])))
    else:
        _set_identity_map(name)
    _invalidate_obj(name)
    _STATE["obj"] = name
    _STATE.setdefault("sources", {})[name] = struct.get("source", "PDB")
    _STATE.setdefault("obj_struct", {})[name] = struct
    if uid:
        _STATE.setdefault("obj_uid", {})[name] = uid
    mp = UFV_MAPS.get(name)
    if mp is not None and mp.get("mode") == "identity" and _is_af_based(name):
        tc = _transplant_chains(name)   # AlphaFill: ligand / donor-peptide chains are not this protein
        if tc:
            mp["exclude"] = frozenset(tc)
    if orient and first:
        cmd.orient(name, animate=0)
        _view_commit()
    return name


def __init_plugin__(app=None):
    """PyMOL Plugin Manager entry point."""
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt("3D Feature Viewer", ufv_gui)
    except Exception:
        pass


if _HAS_PYMOL:
    cmd.extend("ufv_gui", ufv_gui)
if _HAS_PYMOL and __name__ != "__main__":      # stdout of the CLI backend is Tcl: no banner there
    print("[UFV] 3D Feature Viewer loaded — type 'ufv_gui' to open the panel "
          "(or 'help ufv_gui'; all commands are 'ufv_*').")


# ---- one shared backend state ----------------------------------------------------------------------
# Names that still refer to the ufv_core object are backend names: reading them on this module returns
# ufv_core's current binding, assigning / deleting them updates ufv_core as well as this module.
_CORE_NAMES = frozenset(_n for _n, _v in vars(_core).items()
                        if not (_n.startswith("__") and _n.endswith("__")) and globals().get(_n) is _v)


class _UFVModule(types.ModuleType):
    def __getattribute__(self, name):
        d = types.ModuleType.__getattribute__(self, "__dict__")
        if name in d.get("_CORE_NAMES", ()):
            return getattr(d["_core"], name)
        return types.ModuleType.__getattribute__(self, name)

    def __setattr__(self, name, value):
        d = types.ModuleType.__getattribute__(self, "__dict__")
        if name in d.get("_CORE_NAMES", ()):
            setattr(d["_core"], name, value)
        types.ModuleType.__setattr__(self, name, value)

    def __delattr__(self, name):
        d = types.ModuleType.__getattribute__(self, "__dict__")
        if name in d.get("_CORE_NAMES", ()) and hasattr(d["_core"], name):
            delattr(d["_core"], name)
        types.ModuleType.__delattr__(self, name)


def _this_module():
    g = globals()
    m = sys.modules.get(__name__)
    if m is not None and getattr(m, "__dict__", None) is g:
        return m
    import gc                                   # loaded without sys.modules (importlib spec by path)
    for r in gc.get_referrers(g):
        if isinstance(r, types.ModuleType) and r.__dict__ is g:
            return r
    return None


if "__script__" not in globals() and __name__ not in ("__main__", "pymol"):   # not `run` into a shared ns
    _m = _this_module()
    if _m is not None and type(_m) is types.ModuleType:
        _m.__class__ = _UFVModule
    del _m
