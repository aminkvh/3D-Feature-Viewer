"""
3D Feature Viewer — shared backend (ufv_core.py)
=============================================================

Everything that does not need PyMOL: the HTTP layer and on-disk cache, the UniProtKB / EBI Proteins /
AlphaFold DB / PDBe-SIFTS / RCSB / 3D-Beacons / ProtVar / AlphaFill / PDBe-KB / Open Targets / PubChem
fetchers and parsers, UniProt <-> structure numbering (SIFTS, isoform segments, manual / custom maps),
the structure analyses (mutation hotspots, contact hubs, constraint pockets, burden, PTM-variant
proximity, pocket evidence, trajectory maths), the per-residue CSV and the residue / ligand report
formatting. Standard library only (numpy is used when present).

Used two ways:
  * imported by the PyMOL plugin (ufv_pymol.py / the `ufv` plugin package), which adds the PyMOL / Qt
    layer on top and shares this module's state (caches, maps, settings);
  * run as a command-line backend by the VMD plugin (ufv_vmd.tcl). stdout carries only the requested
    payload (Tcl variable assignments or a file path):

    python ufv_core.py --emit-tcl       <uniprot_id>
    python ufv_core.py --emit-sifts     <pdb_id> <uniprot_id>
    python ufv_core.py --structures     <uniprot_id>
    python ufv_core.py --download       <uniprot_id>
    python ufv_core.py --get            <url> <fmt>
    python ufv_core.py --residue        <uniprot_id> <pos>
    python ufv_core.py --compute-tiers  <uniprot_id> <coords.json>
    python ufv_core.py --csv            <uniprot_id> <outfile.csv>
    python ufv_core.py --serve          (persistent JSON-lines request server; see _serve)

This module must never import pymol or Qt.
"""

from __future__ import annotations

import contextlib
import functools
import gzip
import io
import json
import os
import re
import sys
import threading
import time
import urllib.parse
# urllib.request (with http.client, email, ssl: ~75 ms) is imported on the first network call, not at
# import time: plugin start-up and the VMD backend's start-up do not pay for it (_urlreq()).


def _urlreq():
    import urllib.request
    return urllib.request


_SSL_DEFAULT = {"ctx": None}
_SSL_DEFAULT_LOCK = threading.Lock()


def _ssl_default():
    """ONE default HTTPS context for every request. urlopen() without a context builds a new one per
    connection (ssl._create_default_https_context(): ~15 ms here, it re-reads the Windows certificate
    stores each time). Same factory, so the same verification (and a user's override is honoured)."""
    ctx = _SSL_DEFAULT["ctx"]
    if ctx is None:
        import ssl
        with _SSL_DEFAULT_LOCK:
            ctx = _SSL_DEFAULT["ctx"]
            if ctx is None:
                ctx = _SSL_DEFAULT["ctx"] = ssl._create_default_https_context()
    return ctx


def _sha1hex(b):
    import hashlib                          # lazy (OpenSSL load: ~15 ms at import time otherwise)
    return hashlib.sha1(b).hexdigest()


_TMP = {"dir": None}


def _ufv_tmpdir():
    """The directory every temporary file of the backend goes to: <system temp>/ufv. Never the working
    directory: tempfile.gettempdir() silently falls back to os.getcwd() when TEMP / TMP are not set and
    no standard temp folder is usable (e.g. a child process started with a stripped environment) — that
    once wrote an AlphaFill .cif into the repository root. Last resort: ~/.ufv_cache/tmp."""
    d = _TMP["dir"]
    if d and os.path.isdir(d):
        return d
    import tempfile
    cands = [os.environ.get(k) for k in ("TMPDIR", "TEMP", "TMP")]
    if os.environ.get("LOCALAPPDATA"):
        cands.append(os.path.join(os.environ["LOCALAPPDATA"], "Temp"))
    cands += [tempfile.gettempdir(), os.path.join(os.path.expanduser("~"), ".ufv_cache", "tmp")]
    cwd = os.path.normcase(os.path.abspath(os.getcwd()))
    for c in cands:
        if not c or not os.path.isabs(c) or os.path.normcase(os.path.abspath(c)) == cwd:
            continue
        d = os.path.join(c, "ufv")
        try:
            os.makedirs(d, exist_ok=True)
            if os.access(d, os.W_OK):
                _TMP["dir"] = d
                return d
        except OSError:
            continue
    raise OSError("no writable temporary directory for UFV")

# ----------------------------------------------------------------------------------------------
# Endpoints (mirror api.js)
# ----------------------------------------------------------------------------------------------
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
AF_MODEL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v{}.pdb"
FEATURES_URL = "https://www.ebi.ac.uk/proteins/api/features/{}"
VARIATION_URL = "https://www.ebi.ac.uk/proteins/api/variation/{}"
PROTEOMICS_PTM_URL = "https://www.ebi.ac.uk/proteins/api/proteomics-ptm/{}"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
AM_CSV_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-aa-substitutions.csv"
PDBe_SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{}"
PDBe_BEST = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{}"
PDBe_MOLECULES = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{}"
PDBe_PDB = "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{}.ent"
PDBe_CIF = "https://www.ebi.ac.uk/pdbe/entry-files/download/{}.cif"   # large entries have no .ent
BEACONS_SUMMARY = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/{}.json"
PDBe_RESIDUES = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{}"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql?query="
# 3D-Beacons model hosts the extension accepts (api.js BEACON_ALLOWED_HOSTS).
BEACON_ALLOWED_HOSTS = {"alphafold.ebi.ac.uk", "www.ebi.ac.uk", "files.rcsb.org", "swissmodel.expasy.org",
                        "www.modelarchive.org", "modelarchive.org", "proteinensemble.org", "pdb-ihm.org",
                        "alphafill.eu"}
BEACONS_CAP = 25              # computed models kept in the selector (api.js get3DBeaconsModels)
DECORATE_CONCURRENCY = 8      # parallel per-entry metadata fetches (api.js DECORATE_CONCURRENCY)
PARTNER_ACCESSION_CAP = 8     # distinct partner proteins fetched for a complex (api.js)
LIGAND_CCD = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"
PUBCHEM_FP = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/Fingerprint2D/JSON"
# ProtVar (EMBL-EBI) per-residue predictors that complement AlphaMissense: conservation, EVE, ESM1b,
# M3D structural effect, and FoldX ΔΔG stability. Keyed by UniProt accession + position.
PROTVAR_SCORE = "https://www.ebi.ac.uk/ProtVar/api/score/{}/{}"
PROTVAR_FOLDX = "https://www.ebi.ac.uk/ProtVar/api/prediction/foldx/{}/{}"
PROTVAR_POCKET = "https://www.ebi.ac.uk/ProtVar/api/prediction/pocket/{}/{}"
# CADD is genomic (not in /score): /mapping maps "ACC POS REF *" to the codon's single-nucleotide variants,
# each with a caddScore (api.js fetchCaddByMut). alt '*' = return the whole codon's SNV set.
PROTVAR_MAPPING = "https://www.ebi.ac.uk/ProtVar/api/mapping?q={}&assembly=AUTO"
AA3_TO_1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
            "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
            "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PUBCHEM_CID = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/cids/JSON"
ALPHAFILL_JSON = "https://alphafill.eu/v1/aff/{}/json"
PDBEKB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/{}/{}"      # ligand_sites | interface_residues
OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

_USER_AGENT = "UFV-PyMOL-plugin/1.0 (https://github.com/aminkvh/3D-Feature-Viewer-for-UniProt)"

# ----------------------------------------------------------------------------------------------
# Colour tables (mirror data.js)
# ----------------------------------------------------------------------------------------------
PTM_COLORS = {   # data.js PTM_COLORS (colour-blind-safe; no hex shared with the consequence colours)
    "Phosphorylation": "#8b5cf6", "N6-acetyllysine": "#06b6d4", "N6-succinyllysine": "#0e7490",
    "N6-methyllysine": "#16a34a", "N6,N6-dimethyllysine": "#15803d", "Omega-N-methylarginine": "#65a30d",
    "Symmetric dimethylarginine": "#4ade80", "Asymmetric dimethylarginine": "#86efac",
    "Glycyl lysine isopeptide": "#c026d3", "Ubiquitination": "#e879f9", "Ubiquitinated lysine": "#e879f9",
    "SUMOylation": "#9333ea", "ADP-ribosylation": "#6366f1", "Glycosylation": "#be185d",
    "Disulfide bond": "#0d9488", "Lipidation": "#f472b6", "Cross-link": "#a78bfa",
    "Acetylation": "#22d3ee", "Methylation": "#22c55e", "Citrulline": "#2dd4bf",
    "Hydroxyproline": "#38bdf8", "Hydroxylation": "#a3e635", "Nitration": "#facc15",
    "Pyroglutamic acid": "#fbbf24", "S-nitrosocysteine": "#4f46e5", "Deamidation": "#818cf8",
    "Modified residue": "#546e7a", "default": "#ff8fab",
}
# Categories with no PTM_COLORS entry (they would all share the default pink) get these, in order of
# appearance, skipping colours already used (data.js PTM_FALLBACK_PALETTE / groupPTMsByCategory).
PTM_FALLBACK_PALETTE = ["#8b5cf6", "#06b6d4", "#f59e0b", "#ec4899", "#10b981", "#3b82f6",
                        "#ef4444", "#a3e635", "#f97316", "#14b8a6", "#d946ef", "#eab308"]
CONSEQUENCE_COLORS = {   # data.js CONSEQUENCE_CATEGORIES (vermillion / amber / blue / grey)
    "Likely pathogenic or pathogenic": "#d55e00",
    "Predicted deleterious": "#e6a817",
    "Likely benign or benign": "#0072b2",
    "Uncertain significance": "#9e9e9e",
}
PROVENANCE_COLORS = {   # data.js PROVENANCE_CATEGORIES
    "UniProt reviewed": "#42a5f5",
    "ClinVar": "#ab47bc",
    "Large scale studies": "#26c6da",
}
# 50-colour disease palette (data.js DISEASE_PALETTE): diseases are coloured in order of first
# appearance across the variant list; consecutive entries come from different hue families.
DISEASE_PALETTE = [
    "#1565c0", "#bf360c", "#880e4f", "#33691e", "#7b1fa2", "#f9a825", "#00897b", "#6d4c41", "#3949ab", "#039be5",
    "#0288d1", "#e65100", "#c2185b", "#558b2f", "#9c27b0", "#fbc02d", "#004d40", "#795548", "#5c6bc0", "#29b6f6",
    "#0097a7", "#ef6c00", "#e91e63", "#689f38", "#4527a0", "#ff6f00", "#4db6ac", "#a1887f", "#7986cb", "#80deea",
    "#00838f", "#ffa000", "#ec407a", "#7cb342", "#ba68c8", "#ffe082", "#26a69a", "#bcaaa4", "#9fa8da", "#006064",
    "#311b92", "#f06292", "#9ccc65", "#7c4dff", "#01579b", "#0277bd", "#ff8f00", "#00695c", "#4e342e", "#fdd835",
]
UNCLASSIFIED_COLOR = "#9e9e9e"
TOPOLOGY_COLORS = {
    "transmembrane": "#f9a825", "intramembrane": "#ef6c00", "cytoplasmic": "#1e88e5",
    "extracellular": "#e53935", "lumenal": "#43a047", "periplasmic": "#26a69a",
    "nuclear": "#8e24aa", "mitochondrial": "#00897b", "other": "#7e57c2",
}
DOMAIN_PALETTE = ["#1e88e5", "#43a047", "#fb8c00", "#8e24aa", "#00897b", "#e53935", "#3949ab",
                  "#c0ca33", "#6d4c41", "#00acc1", "#d81b60", "#5e35b1", "#7cb342", "#f4511e",
                  "#039be5", "#fdd835"]
SITE_COLOR = "#fbc02d"
DOMAIN_TYPES = ["Domain", "Region", "Repeat", "Compositional bias", "Zinc finger",
                "Coiled coil", "Motif", "DNA binding"]
PTM_TYPES = {"MOD_RES", "CROSSLNK", "LIPID", "CARBOHYD", "DISULFID"}
STRUCTURE_LIST_CAP = 400  # safety cap on experimental PDB ENTRIES listed (the extension lists all)

# ----------------------------------------------------------------------------------------------
# Module state
# ----------------------------------------------------------------------------------------------
_CACHE = {}          # uniprot_id -> annotations dict
UFV_MAPS = {}        # object name -> mapping descriptor
_STATE = {"uid": None}


# ----------------------------------------------------------------------------------------------
# HTTP helpers
# ----------------------------------------------------------------------------------------------
_SSL_FALLBACK = {"ctx": None, "hosts": set()}
_SSL_FALLBACK_LOCK = threading.Lock()


def _ssl_fallback_context():
    """A still-verifying TLS context for hosts the default one rejects. On Windows Python's default context
    also loads the machine's intermediate ("CA") store, which can hold an EXPIRED cross-signed root (e.g.
    "ISRG Root X2" issued by X1, expired 2025-09); OpenSSL then builds the chain through it and fails with
    "certificate has expired" although the server's own chain is valid (swissmodel.expasy.org, Let's
    Encrypt's 2026 "Root YE" chain). The fallback trusts certifi's bundle when present, else only the
    Windows ROOT store."""
    import ssl
    with _SSL_FALLBACK_LOCK:
        if _SSL_FALLBACK["ctx"] is not None:
            return _SSL_FALLBACK["ctx"]
        ctx = None
        try:
            import certifi
            ctx = ssl.create_default_context(cafile=certifi.where())
        except Exception:
            ctx = None
        if ctx is None and hasattr(ssl, "enum_certificates"):
            ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
            for der, enc, _trust in ssl.enum_certificates("ROOT"):
                if enc == "x509_asn":
                    try:
                        ctx.load_verify_locations(cadata=der)
                    except Exception:
                        pass
        _SSL_FALLBACK["ctx"] = ctx or False
        return _SSL_FALLBACK["ctx"]


def _is_cert_error(exc):
    import ssl
    r = getattr(exc, "reason", exc)
    return isinstance(r, ssl.SSLCertVerificationError) or "CERTIFICATE_VERIFY_FAILED" in str(r)


def _http_get(url):
    if _OFFLINE:
        raise OSError("offline (UFV_OFFLINE=1)")
    ureq = _urlreq()
    req = ureq.Request(url, headers={"User-Agent": _USER_AGENT, "Accept": "*/*"})
    host = urllib.parse.urlparse(url).hostname
    if host not in _SSL_FALLBACK["hosts"]:
        try:
            with ureq.urlopen(req, timeout=60, context=_ssl_default()) as r:
                return r.read()
        except Exception as exc:
            if not _is_cert_error(exc) or not _ssl_fallback_context():
                raise
    ctx = _ssl_fallback_context()
    with ureq.urlopen(req, timeout=60, context=ctx) as r:
        data = r.read()
    _SSL_FALLBACK["hosts"].add(host)      # verified with the fallback: use it straight away next time
    return data


# Small on-disk cache of annotation/metadata responses (UniProt, Proteins API, AlphaMissense CSV,
# PDBe, 3D-Beacons, RCSB, PubChem) and of downloaded structure files, so a new PyMOL session doesn't
# re-download megabytes for the same accession / entry. One gzip file per URL in ~/.ufv_cache (override:
# UFV_CACHE_DIR; disable: UFV_NO_CACHE=1), expiring after UFV_CACHE_TTL_DAYS (default 7); responses
# above DISK_CACHE_MAX_ITEM are not stored. Every cache operation is best-effort: any failure
# (read-only home, corrupt file, full disk) silently falls back to the network.
DISK_CACHE_MAX_ITEM = 64 * 1024 * 1024
_DISK_CACHE_DIR = os.environ.get("UFV_CACHE_DIR") or os.path.join(os.path.expanduser("~"), ".ufv_cache")
_DISK_CACHE_ON = not os.environ.get("UFV_NO_CACHE")
# UFV_OFFLINE=1: never touch the network; every request is served from the disk cache (any age) or fails.
# Used by the tests with UFV_CACHE_DIR pointing at a recorded fixture cache.
_OFFLINE = bool(os.environ.get("UFV_OFFLINE"))
try:
    _DISK_CACHE_TTL = float(os.environ.get("UFV_CACHE_TTL_DAYS", "7")) * 86400
except ValueError:
    _DISK_CACHE_TTL = 7 * 86400
_URL_LOCKS = {}
_URL_LOCKS_GUARD = threading.Lock()
_RECENT = []          # [(url, bytes)] last few cached responses, so concurrent callers share one
_RECENT_MAX, _RECENT_MAX_BYTES = 6, 8 * 1024 * 1024    # download even when the disk cache is off


def _cache_file(url):
    return os.path.join(_DISK_CACHE_DIR, _sha1hex(url.encode("utf-8")) + ".gz")


def _cache_read(url, stale=False):
    """Cached bytes of `url` or None. `stale`: ignore the TTL (the network just failed — an old copy is
    better than nothing; also every read in UFV_OFFLINE mode)."""
    if not _DISK_CACHE_ON:
        return None
    try:
        path = _cache_file(url)
        if not (stale or _OFFLINE) and time.time() - os.path.getmtime(path) > _DISK_CACHE_TTL:
            return None
        with gzip.open(path, "rb") as fh:
            return fh.read()
    except Exception:
        return None


def _cache_write(url, data):
    if not _DISK_CACHE_ON or not data:
        return
    tmp = None
    try:
        os.makedirs(_DISK_CACHE_DIR, exist_ok=True)
        path = _cache_file(url)
        tmp = "%s.%d.%d.tmp" % (path, os.getpid(), threading.get_ident())
        with gzip.open(tmp, "wb", compresslevel=5) as fh:
            fh.write(data)
        os.replace(tmp, path)                 # atomic: readers never see a half-written file
        tmp = None
    except Exception:
        pass
    finally:
        if tmp:
            try:
                os.remove(tmp)
            except OSError:
                pass


def _url_lock(url):
    with _URL_LOCKS_GUARD:
        lk = _URL_LOCKS.get(url)
        if lk is None:
            lk = _URL_LOCKS[url] = threading.Lock()
        return lk


def _fetch_bytes(url, cache=False, refresh=False, validate=None, cache_404=None):
    """GET `url` -> bytes or None. With `cache`, serve from / store to the disk cache; concurrent
    requests for the same URL wait for the first (so e.g. the UniProt entry JSON that both
    fetch_annotations and list_structures need is downloaded once). `validate(bytes)` must not raise
    for a response to be cached. `refresh` skips the cache read (still writes). `cache_404`: bytes to
    store and return for a definitive HTTP 404 (PDBe answers 404 for "no SIFTS mapping"; without this
    every new session asked again). Responses above DISK_CACHE_MAX_ITEM are never written to disk."""
    if not cache:
        try:
            return _http_get(url)
        except Exception:
            return _cache_read(url, stale=True)    # network down: a recorded copy, if any (never written here)
    with _url_lock(url):
        if not refresh:
            for u, raw in list(_RECENT):
                if u == url:
                    return raw
            raw = _cache_read(url)
            if raw is not None:
                try:
                    if validate:
                        validate(raw)
                    return raw
                except Exception:
                    pass
        try:
            raw = _http_get(url)
        except Exception as exc:
            if cache_404 is not None and getattr(exc, "code", None) == 404:
                _cache_write(url, cache_404)
                return cache_404
            return _cache_read(url, stale=True)    # network down: serve an expired copy rather than nothing
        try:
            if validate:
                validate(raw)
            if len(raw) <= DISK_CACHE_MAX_ITEM:
                _cache_write(url, raw)
            _RECENT[:] = [e for e in _RECENT if e[0] != url]   # a refresh replaces the old copy
            if len(raw) <= _RECENT_MAX_BYTES:
                _RECENT.append((url, raw))
                del _RECENT[:-_RECENT_MAX]
        except Exception:
            pass
        return raw


_FETCH_BYTES_REAL = _fetch_bytes      # (tests substitute _fetch_bytes with recorded responses)


def _get_json_raw(url, cache=True, refresh=False):
    """(raw bytes or None when nothing could be fetched, parsed JSON or None). The validator's parse is
    kept (tied to the exact bytes object it parsed) and returned, so a cached or freshly downloaded
    response is parsed ONCE, not once to validate and once more to use."""
    box = []

    def validate(b):
        box[:] = [b, json.loads(b.decode("utf-8"))]
    raw = _fetch_bytes(url, cache, refresh, validate=validate, cache_404=b"null")
    if raw is None:
        return None, None
    if box and box[0] is raw:
        return raw, box[1]
    try:
        return raw, json.loads(raw.decode("utf-8"))
    except Exception:
        return raw, None


def _get_json(url, cache=True, refresh=False):
    return _get_json_raw(url, cache, refresh)[1]


def _get_text(url, cache=False, refresh=False):
    raw = _fetch_bytes(url, cache, refresh)
    if raw is None:
        return None
    try:
        return raw.decode("utf-8", "replace")
    except Exception:
        return None


def _get_structure_text(url):
    """A structure file (AlphaFold / PDB / mmCIF / computed model) as text. Disk-cached like the
    metadata (a re-load, also in a new session, reads ~/.ufv_cache instead of re-downloading megabytes;
    files above DISK_CACHE_MAX_ITEM are not stored), but not kept in the in-memory _RECENT list (read
    once; it would only pin megabytes). Same fallbacks as _get_text."""
    raw = _fetch_bytes(url, True)
    if _RECENT:
        with _url_lock(url):
            _RECENT[:] = [e for e in _RECENT if e[0] != url]
    if raw is None:
        return None
    try:
        return raw.decode("utf-8", "replace")
    except Exception:
        return None


# ----------------------------------------------------------------------------------------------
# Data extraction (mirror data.js)
# ----------------------------------------------------------------------------------------------
def _categorize_ptm(desc, type_):
    d = (desc or "").lower()
    if type_ == "DISULFID":
        return "Disulfide bond"
    if type_ == "LIPID":
        return "Lipidation"
    if type_ == "CARBOHYD":
        return "Glycosylation"
    if type_ == "CROSSLNK":
        if "ubiquitin" in d:
            return "Ubiquitination"
        if "sumo" in d:
            return "SUMOylation"
        return "Cross-link"
    if "phospho" in d:
        return "Phosphorylation"
    if "n6-acetyl" in d:
        return "N6-acetyllysine"
    if "n6-succinyl" in d:
        return "N6-succinyllysine"
    if "n6,n6-dimethyl" in d:
        return "N6,N6-dimethyllysine"
    if "n6-methyl" in d:
        return "N6-methyllysine"
    if "omega-n-methyl" in d:
        return "Omega-N-methylarginine"
    if "symmetric dimethyl" in d:
        return "Symmetric dimethylarginine"
    if "asymmetric dimethyl" in d:
        return "Asymmetric dimethylarginine"
    if "adp-ribos" in d:
        return "ADP-ribosylation"
    if "ubiquitin" in d:
        return "Ubiquitinated lysine"
    if "glycyl lysine isopeptide" in d:
        return "Glycyl lysine isopeptide"
    if "citrulline" in d:
        return "Citrulline"
    if "hydroxyproline" in d:
        return "Hydroxyproline"
    if "hydroxylysine" in d or "hydroxylation" in d:
        return "Hydroxylation"
    if "nitrated" in d:
        return "Nitration"
    if "pyroglutam" in d:
        return "Pyroglutamic acid"
    if "s-nitrosocysteine" in d:
        return "S-nitrosocysteine"
    if "deamidated" in d:
        return "Deamidation"
    label = (desc or "").split(";")[0].strip()
    return label or "Modified residue"


def _extract_ptms(features):
    out = []
    for f in (features or {}).get("features", []):
        if f.get("type") not in PTM_TYPES:
            continue
        try:
            pos = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        cat = _categorize_ptm(f.get("description", ""), f.get("type"))
        out.append({
            "position": pos,
            "endPosition": _js_int(f.get("end")) or pos,
            "type": f.get("type"),
            "description": f.get("description") or "Unknown modification",
            "category": cat,
            "color": PTM_COLORS.get(cat, PTM_COLORS["default"]),
        })
    return out


def _js_int(v):
    """JavaScript parseInt for the leading integer of a value (None when JS would give NaN)."""
    if isinstance(v, bool):
        return None
    if isinstance(v, int):
        return v
    if isinstance(v, float):
        return int(v) if v == v and abs(v) != float("inf") else None
    m = re.match(r"\s*([+-]?\d+)", str(v if v is not None else ""))
    return int(m.group(1)) if m else None


def _extract_proteomics_ptms(data):
    """Large-scale PTMs from the EBI proteomics-ptm endpoint (data.js extractProteomicsPTMs): peptide-
    relative positions made absolute, one entry per (position, PTM name), categorised like MOD_RES."""
    out, seen = [], set()
    for f in (data or {}).get("features", []) or []:
        if f.get("type") != "PROTEOMICS_PTM" or not isinstance(f.get("ptms"), list):
            continue
        start = _js_int(f.get("begin"))
        if start is None:
            continue
        for ptm in f["ptms"]:
            try:
                abs_pos = start + (int(ptm.get("position")) - 1)
            except (TypeError, ValueError):
                continue
            key = (abs_pos, ptm.get("name"))
            if key in seen:
                continue
            seen.add(key)
            desc = ptm.get("name") or "Unknown modification"
            cat = _categorize_ptm(desc.lower(), "MOD_RES")
            out.append({
                "position": abs_pos, "endPosition": abs_pos, "type": "MOD_RES",
                "description": "%s (large scale data, %s)" % (desc, ", ".join(str(s) for s in (ptm.get("sources") or []))),
                "category": cat, "color": PTM_COLORS.get(cat, PTM_COLORS["default"]),
                "isLargeScale": True,
            })
    return out


def _merge_ptms(curated, large_scale):
    """api.js loadFeatureData: curated EBI PTMs first, then large-scale ones at (position, category)
    pairs not already present (curated entries win)."""
    seen = {(p["position"], p["category"]) for p in curated}
    out = list(curated)
    for p in large_scale:
        k = (p["position"], p["category"])
        if k not in seen:
            seen.add(k)
            out.append(p)
    return out


def _group_ptm_colors(ptms):
    """data.js groupPTMsByCategory: categories whose colour is the shared default pink get distinct hues
    from PTM_FALLBACK_PALETTE (in order of first appearance, skipping colours already in use); the
    colour is written back onto the items. Returns {category: colour} in first-appearance order."""
    groups = {}
    for p in ptms:
        groups.setdefault(p["category"], p["color"])
    default = PTM_COLORS["default"]
    used = {c for c in groups.values() if c != default}
    pal = [c for c in PTM_FALLBACK_PALETTE if c not in used]
    fi = 0
    for cat, col in list(groups.items()):
        if col != default:
            continue
        groups[cat] = (pal[fi % len(pal)] if pal else col)
        fi += 1
    for p in ptms:
        p["color"] = groups[p["category"]]
    return groups


def _classify_consequence(v):
    sigs = v.get("clinicalSignificances") or []
    for s in sigs:
        t = (s.get("type") or "").lower()
        if "pathogenic" in t and "benign" not in t and "uncertain" not in t:
            return "Likely pathogenic or pathogenic"
    for s in sigs:
        t = (s.get("type") or "").lower()
        if "benign" in t and "pathogenic" not in t:
            return "Likely benign or benign"
    for s in sigs:
        if "uncertain" in (s.get("type") or "").lower():
            return "Uncertain significance"
    for p in (v.get("predictions") or []):
        label = " ".join(str(p.get(k, "")) for k in
                         ("predictionValType", "predictionVal", "value", "name")).lower()
        if re.search(r"\b(deleterious|damaging|pathogenic)\b", label) and \
           not re.search(r"\b(benign|tolerated)\b", label):
            return "Predicted deleterious"
    return "Uncertain significance"


def _classify_provenance(v):
    src = (v.get("sourceType") or "").lower()
    if src in ("uniprot", "mixed"):
        return "UniProt reviewed"
    for x in (v.get("xrefs") or []):
        if "clinvar" in (x.get("name") or "").lower():
            return "ClinVar"
    return "Large scale studies"


_DISEASE_ABBR_EXCLUDE = ("ECO", "MIM", "NCI", "VAR", "PRO")


def _disease_data(v):
    """data.js _extractDiseaseData: (labels, ids, pairs). Labels/ids keep first-appearance order (JS Set);
    pairs [{id, label}] pair each association's identifier with ITS label. Fallbacks when the association
    array has no disease: an abbreviation leading a descriptions[] value, then 'in ABC [and DEF]' in the
    flat description text."""
    labels, ids, pairs = [], [], []

    def add(lst, x):
        if x not in lst:
            lst.append(x)
    for a in (v.get("association") or []) if isinstance(v.get("association"), list) else []:
        aid = a.get("id") or a.get("accession") or ""
        if aid:
            add(ids, aid)
        if a.get("disease") is True and a.get("name"):
            m = re.search(r"\(([A-Z][A-Z0-9]+)\)\s*$", a["name"])
            label = m.group(1) if m else a["name"]
            add(labels, label)
            pairs.append({"id": aid or "", "label": label})
    for d in (v.get("descriptions") or []) if isinstance(v.get("descriptions"), list) else []:
        val = d.get("value") if isinstance(d, dict) else None
        if not val:
            continue
        m = re.match(r"^([A-Z][A-Z0-9]{1,15})(?:;|$)", val)
        if m and m.group(1) not in _DISEASE_ABBR_EXCLUDE and len(m.group(1)) > 1:
            add(labels, m.group(1))
    if not labels:
        desc = v.get("description") or ""
        for mm in re.finditer(r"\bin\s+([A-Z][A-Z0-9]{1,10}(?:\s+and\s+[A-Z][A-Z0-9]{1,10})*)\b", str(desc)):
            for part in re.split(r"\s+and\s+", mm.group(1)):
                t = part.strip()
                if t and t != "dbSNP" and len(t) > 1:
                    add(labels, t)
    return labels, ids, pairs


def _disease_labels(v):
    return _disease_data(v)[0]


def _gnomad_af(f):
    """Largest gnomAD allele frequency from a variation feature's populationFrequencies (the EBI
    variation API already carries gnomAD v4 exome/genome AFs — we just surface them). None if absent."""
    best = None
    for pf in (f.get("populationFrequencies") or []):
        if "gnomad" in (pf.get("source") or "").lower():
            fr = pf.get("frequency")
            if isinstance(fr, (int, float)):
                best = fr if best is None else max(best, fr)
    return best


def _genomic_hgvs(f):
    g = f.get("genomicLocation")
    if isinstance(g, list):
        return g[0] if g else None
    return g or None


_RE_RSID = re.compile(r"^rs\d+$", re.I)


def _extract_variants(variation, am_map):
    out = []
    for f in (variation or {}).get("features", []):
        if f.get("type") != "VARIANT":
            continue
        try:
            pos = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        if not f.get("alternativeSequence"):
            continue
        cons = _classify_consequence(f)
        am = None
        if am_map and f.get("wildType") and len(f.get("alternativeSequence", "")) == 1:
            am = am_map.get("{}{}{}".format(f["wildType"], f["begin"], f["alternativeSequence"]))
        if am is None:      # fallback: an AlphaMissense prediction carried by the EBI variation record
            for p in (f.get("predictions") or []):
                lab = ("%s %s" % (p.get("predAlgorithmNameType") or p.get("predictionName") or "",
                                  p.get("name") or "")).lower()
                sc = p.get("score")
                if sc is None:
                    sc = p.get("predictionVal")
                if sc is None:
                    sc = p.get("value")
                try:
                    sc = float(sc)
                except (TypeError, ValueError):
                    continue
                if "alphamissense" in lab and sc == sc and abs(sc) != float("inf"):
                    am = sc
        sigs = f.get("clinicalSignificances") or []
        rs = []
        for x in (f.get("xrefs") or []):
            val = x.get("id") or x.get("name") or x.get("value") or ""
            if _RE_RSID.match(str(val)) and val not in rs:
                rs.append(val)
            if "dbsnp" in (x.get("database") or x.get("name") or "").lower() and val and val not in rs:
                rs.append(val)
        prov = _classify_provenance(f)
        labels, ids, pairs = _disease_data(f)
        descs = [d.get("value") for d in (f.get("descriptions") or []) if isinstance(d, dict)]
        out.append({
            "position": pos,
            "wildType": f.get("wildType") or "?",
            "mutant": f.get("alternativeSequence"),
            "consequence": cons,
            "consequenceColor": CONSEQUENCE_COLORS.get(cons, "#9e9e9e"),
            "provenance": prov,
            "provenanceColor": PROVENANCE_COLORS.get(prov, "#26c6da"),
            "reviewed": prov == "UniProt reviewed",
            "sourceType": f.get("sourceType") or "unknown",
            "description": "; ".join(str(d) for d in descs) if descs else
                           ", ".join(str(s.get("type")) for s in sigs),
            "diseases": labels,
            "diseaseIds": ids,
            "diseasePairs": pairs,
            "clinVar": ", ".join(s.get("type", "") for s in sigs if s.get("type")),
            "clinVarReview": ", ".join(filter(None, (s.get("reviewStatus") or s.get("review_status") or "" for s in sigs))),
            "rsIds": rs,
            "alphaMissense": am,
            "gnomad": _gnomad_af(f),
            "genomic": _genomic_hgvs(f),
        })
    return out


def _compute_disease_colors(variants):
    """data.js computeDiseaseColors / getDiseaseSummary: one DISEASE_PALETTE colour per disease label in
    order of first appearance; each variant's diseaseColor = its FIRST disease's colour (grey if none).
    Returns the ordered {label: colour} map."""
    cmap = {}
    for v in variants:
        for d in v.get("diseases") or []:
            if d not in cmap:
                cmap[d] = DISEASE_PALETTE[len(cmap) % len(DISEASE_PALETTE)]
    for v in variants:
        ds = v.get("diseases") or []
        v["diseaseColor"] = cmap[ds[0]] if ds else UNCLASSIFIED_COLOR
    return cmap


def _disease_summary(variants, cmap=None):
    """{label: {count, color}} in first-appearance order, plus 'Unclassified' for variants with no disease
    (data.js getDiseaseSummary)."""
    cmap = cmap if cmap is not None else _compute_disease_colors(variants)
    s = {}
    for v in variants:
        for d in v.get("diseases") or []:
            e = s.setdefault(d, {"count": 0, "color": cmap.get(d, UNCLASSIFIED_COLOR)})
            e["count"] += 1
    nd = sum(1 for v in variants if not v.get("diseases"))
    if nd:
        s["Unclassified"] = {"count": nd, "color": UNCLASSIFIED_COLOR}
    return s


def _fmt_af(af):
    """Human-readable gnomAD allele frequency."""
    if af is None:
        return None
    if af == 0:
        return "0 (not observed)"
    if af >= 0.01:
        return "%.2f%%" % (af * 100)
    if af >= 1e-4:
        return "%.4f%%" % (af * 100)
    return "%.1e" % af


def _extract_sites(features):
    labels = {"SITE": "Site", "ACT_SITE": "Active site", "BINDING": "Binding site", "METAL": "Metal binding"}
    out = []
    for f in (features or {}).get("features", []):
        lab = labels.get(f.get("type"))
        if not lab:
            continue
        try:
            pos = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        lig = (f.get("ligand") or {}).get("name") or ""
        if f.get("type") == "SITE":
            desc = f.get("description") or "Site"
        else:
            parts = [p for p in (lig, f.get("description")) if p]
            desc = "{}: {}".format(lab, " — ".join(parts)) if parts else lab
        out.append({"position": pos, "endPosition": int(f.get("end") or pos),
                    "description": desc, "color": SITE_COLOR})
    return out


def _extract_sites_uniprot(uniprot):
    labels = {"Site": "Site", "Active site": "Active site", "Binding site": "Binding site",
              "Metal binding": "Metal binding"}
    out = []
    for f in (uniprot or {}).get("features", []):
        lab = labels.get(f.get("type"))
        if not lab:
            continue
        pos = ((f.get("location") or {}).get("start") or {}).get("value")
        if pos is None:
            continue
        end = ((f.get("location") or {}).get("end") or {}).get("value") or pos
        lig = (f.get("ligand") or {}).get("name") or ""
        note = (f.get("ligand") or {}).get("note") or ""
        if f.get("type") == "Site":
            desc = f.get("description") or "Site"
        else:
            parts = [p for p in (lig, note, f.get("description")) if p]
            desc = "{}: {}".format(lab, "; ".join(parts)) if parts else lab
        out.append({"position": pos, "endPosition": end, "description": desc, "color": SITE_COLOR})
    return out


def _extract_topology(features):
    out = []
    for f in (features or {}).get("features", []):
        if f.get("type") not in ("TOPO_DOM", "TRANSMEM", "INTRAMEM"):
            continue
        try:
            start = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        end = int(f.get("end") or start)
        C = TOPOLOGY_COLORS
        if f.get("type") == "TRANSMEM":
            label, color = "Transmembrane", C["transmembrane"]
        elif f.get("type") == "INTRAMEM":
            label, color = "Intramembrane", C["intramembrane"]
        else:
            d = (f.get("description") or "").lower()
            label = f.get("description") or "Topological domain"
            color = (C["cytoplasmic"] if "cytoplasm" in d else
                     C["extracellular"] if "extracellular" in d else
                     C["lumenal"] if "lumen" in d else
                     C["periplasmic"] if "periplasm" in d else
                     C["nuclear"] if "nuclear" in d else
                     C["mitochondrial"] if "mitochond" in d else C["other"])
        out.append({"start": start, "end": end, "label": label, "color": color})
    return out


def _extract_domains(uniprot):
    out = []
    feats = [f for f in (uniprot or {}).get("features", []) if f.get("type") in DOMAIN_TYPES]
    for i, f in enumerate(feats):
        start = ((f.get("location") or {}).get("start") or {}).get("value")
        if start is None:
            continue
        end = ((f.get("location") or {}).get("end") or {}).get("value") or start
        desc = f.get("description") or f.get("type")
        out.append({"position": start, "endPosition": end, "type": f.get("type"),
                    "description": desc, "color": DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)],
                    "isRange": end > start})
    return sorted(out, key=lambda d: d["position"])


MUTAG_COLOR = "#6d4c41"   # experimentally mutated residues (UniProt Mutagenesis features)


def _extract_mutagenesis(uniprot):
    """UniProt 'Mutagenesis' features: residues experimentally mutated in the lab, with the observed
    effect. The experimental counterpart to the variant-effect predictors. -> [{position, wildType,
    mutants, effect, pubmed, color}]."""
    out = []
    for f in (uniprot or {}).get("features", []):
        if f.get("type") != "Mutagenesis":
            continue
        loc = f.get("location") or {}
        start = (loc.get("start") or {}).get("value")
        if start is None:
            continue
        end = (loc.get("end") or {}).get("value") or start
        alt = f.get("alternativeSequence") or {}
        pubmed = [e.get("id") for e in (f.get("evidences") or [])
                  if e.get("source") == "PubMed" and e.get("id")]
        out.append({"position": start, "endPosition": end,
                    "wildType": alt.get("originalSequence") or "",
                    "mutants": alt.get("alternativeSequences") or [],
                    "effect": f.get("description") or "", "pubmed": pubmed, "color": MUTAG_COLOR})
    return sorted(out, key=lambda d: d["position"])


def _strip_refs(text):
    """Trim UniProt inline citations for a tidy one-line display."""
    text = re.sub(r"\s*\(PubMed:[^)]*\)", "", text or "")
    text = re.sub(r"\s*\(By similarity\)", "", text)
    text = re.sub(r"\s*\(Ref\.[^)]*\)", "", text)
    return re.sub(r"\s+", " ", text).strip()


def _extract_function(uniprot):
    """Protein-level function context from the entry comments: a short function summary, subcellular
    locations, and catalytic activity. None if the entry carries nothing useful."""
    if not uniprot:
        return None
    summary, locations, catalytic = "", [], []
    for c in uniprot.get("comments") or []:
        t = c.get("commentType")
        if t == "FUNCTION" and not summary:
            texts = c.get("texts") or []
            if texts:
                summary = _strip_refs(texts[0].get("value") or "")
        elif t == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations") or []:
                v = (sl.get("location") or {}).get("value")
                if v and v not in locations:
                    locations.append(v)
        elif t == "CATALYTIC ACTIVITY":
            rxn = (c.get("reaction") or {}).get("name")
            if rxn and rxn not in catalytic:
                catalytic.append(rxn)
    if not (summary or locations or catalytic):
        return None
    return {"summary": summary, "locations": locations, "catalytic": catalytic}


def _name_with_evidence(desc):
    """(fullName value, evidences) from a proteinDescription's recommendedName or first submissionName."""
    if not desc:
        return None, []
    rec = desc.get("recommendedName") or {}
    fn = rec.get("fullName")
    if fn and fn.get("value"):
        return fn["value"], fn.get("evidences", []) or []
    for sub in desc.get("submissionNames") or []:
        fn = sub.get("fullName")
        if fn and fn.get("value"):
            return fn["value"], fn.get("evidences", []) or []
    return None, []


def _extract_protnlm(uniprot):
    """UniProt's ProtNLM (and ProtNLM2) AI model names proteins that would otherwise be
    'uncharacterized', from sequence alone. Those names appear in the entry's proteinDescription with
    computational evidence (and a CAUTION note). Returns {name, isAI, source, caution, reviewed} or
    None — `isAI` flags an AI-generated name so the UI can show it as a prediction."""
    if not uniprot:
        return None
    name, evs = _name_with_evidence(uniprot.get("proteinDescription") or {})
    if not name:
        return None
    is_ai, source = False, ""
    for e in evs:
        tag = ("%s %s" % (e.get("source", ""), e.get("id", ""))).strip()
        # ProtNLM names are attributed to source 'Google' / id 'ProtNLM' (so 'ProtNLM2' matches too).
        # We deliberately do NOT flag rule-based automatic names (ARBA/RuleBase/HAMAP) as AI.
        if "protnlm" in tag.lower() or e.get("source", "").lower() == "google":
            is_ai, source = True, tag
            break
    caution = ""
    for c in uniprot.get("comments") or []:
        if c.get("commentType") == "CAUTION":
            for t in c.get("texts") or []:
                val = t.get("value") or ""
                if "protnlm" in val.lower():
                    caution, is_ai = val, True
    reviewed = "reviewed" in (uniprot.get("entryType") or "").lower() and \
               "unreviewed" not in (uniprot.get("entryType") or "").lower()
    return {"name": name, "isAI": is_ai, "source": source or ("ProtNLM" if is_ai else ""),
            "caution": caution, "reviewed": reviewed}


def _parse_am_csv(text):
    """Parse the AlphaFold AlphaMissense substitution CSV into {'<wt><pos><mut>': score}."""
    out = {}
    if not text:
        return out
    for line in text.splitlines()[1:]:
        line = line.strip()
        if not line or "," not in line:
            continue
        var, rest = line.split(",", 1)
        try:
            out[var] = float(rest.split(",")[0])
        except ValueError:
            continue
    return out


_RE_AM_POS = re.compile(r"^[A-Z](\d+)[A-Z]$")          # precompiled: tens of thousands of keys per protein
_RE_AM_KEY = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def _am_mean_by_position(am_map):
    """Mean AlphaMissense score per residue position (across substitutions)."""
    acc = {}
    match = _RE_AM_POS.match
    for key, score in am_map.items():
        m = match(key)
        if not m:
            continue
        pos = int(m.group(1))
        s = acc.setdefault(pos, [0.0, 0])
        s[0] += score
        s[1] += 1
    return {pos: tot / n for pos, (tot, n) in acc.items() if n}


BURDEN_COMPOSITE_MIN = 0.90   # top-decile composite percentile rank (analysis.js COMPOSITE_MIN)


def _am_mean_from_variants(variants):
    """Fallback when there is no AlphaMissense CSV: mean of the annotated variants' AM scores per position."""
    acc = {}
    for v in variants:
        sc = v.get("alphaMissense")
        if isinstance(sc, (int, float)) and sc == sc and abs(sc) != float("inf"):
            acc.setdefault(v["position"], []).append(sc)
    return {p: sum(s) / len(s) for p, s in acc.items()}


def _compute_burden(variants):
    """Recurrent phenotype residues ("residue burden", analysis.js computeResidueBurden; Akbari Ahangar
    et al., iScience 2024). Per position: mutCount = variant records, phenCount = distinct phenotypes
    (diseasePairs keyed by id, else label; a label already paired with an id is not counted again; plus
    bare disease labels). Each metric becomes a percentile rank (fraction of positions with a value <=
    it); flagged when mutCount >= 2, phenCount >= 2 and 0.5*rank(mut) + 0.5*rank(phen) >= 0.90."""
    import bisect
    mut, phen, ident = {}, {}, {}
    for v in variants:
        pos = v["position"]
        mut[pos] = mut.get(pos, 0) + 1
        ph = phen.setdefault(pos, set())
        idl = ident.setdefault(pos, set())
        for p in v.get("diseasePairs") or []:
            key = p.get("id") or p.get("label")
            if p.get("id") and p.get("label"):
                idl.add(p["label"])
                ph.discard(p["label"])
            if key:
                ph.add(key)
        for label in v.get("diseases") or []:
            if label and label not in idl:
                ph.add(label)
    if not mut:
        return set()
    positions = list(mut)
    pc = {p: len(phen.get(p) or ()) for p in positions}
    sm = sorted(mut.values())
    sp = sorted(pc.values())
    n = float(len(positions))
    out = set()
    for p in positions:
        c, k = mut[p], pc[p]
        if c < 2 or k < 2:
            continue
        comp = 0.5 * (bisect.bisect_right(sm, c) / n) + 0.5 * (bisect.bisect_right(sp, k) / n)
        if comp >= BURDEN_COMPOSITE_MIN:
            out.add(p)
    return out


# ----------------------------------------------------------------------------------------------
# Geometry-based analyses (ported from analysis.js). All operate on a `geom` list of modelled Ca
# atoms — each {uniPos, chain, resi, x, y, z} — so they are unit-testable without PyMOL. The
# permutation nulls use the extension's own PRNG (mulberry32, seed 0x9e3779b9) with the same draw
# order, so hotspot tiers and constraint-pocket p/q values are the extension's numbers exactly
# (HOTSPOT_RNG = "python" restores the old random.Random null).
# ----------------------------------------------------------------------------------------------
import math as _math
import random as _random

HOTSPOT_PERMUTATIONS = 1000
HOTSPOT_SEED = 0x9e3779b9 & 0xFFFFFFFF
HOTSPOT_MIN_PATH = 3
HOTSPOT_MIN_BENIGN = 2
HUB_Z_STRONG = 3.0
HUB_Z_MODERATE = 2.0
HUB_MAX_CA = 6000
_TIER_RANK = {"strong": 3, "moderate": 2, "weak": 1}


def _d2(a, b):
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return dx * dx + dy * dy + dz * dz


def _numpy():
    """numpy if importable (PyMOL ships it), else None — every numpy path has a pure-Python twin."""
    try:
        import numpy as np
        return np
    except Exception:
        return None


def _adjacency(coords, threshold):
    """Neighbour-index lists (ascending, self excluded) for points within `threshold` Å. Vectorised
    with numpy (block-wise distance scan, one nonzero() per block) so the O(N²) scan that dominated
    the analyses runs in milliseconds; falls back to a pure-Python loop if numpy is unavailable."""
    n = len(coords)
    th2 = threshold * threshold
    adj = [[] for _ in range(n)]
    np = _numpy()
    if np is not None and n:
        try:
            pts = np.asarray(coords, dtype=np.float64)
            block = 512
            for i0 in range(0, n, block):
                i1 = min(i0 + block, n)
                within = ((pts[i0:i1, None, :] - pts[None, :, :]) ** 2).sum(-1) <= th2
                within[np.arange(i1 - i0), np.arange(i0, i1)] = False      # no self-loops
                rows, cols = np.nonzero(within)                             # row-major: cols ascend
                bounds = np.searchsorted(rows, np.arange(i1 - i0 + 1))
                cols = cols.tolist()
                for ii in range(i1 - i0):
                    adj[i0 + ii] = cols[bounds[ii]:bounds[ii + 1]]
            return adj
        except Exception:
            adj = [[] for _ in range(n)]
    for i in range(n):
        ax, ay, az = coords[i]
        for j in range(i + 1, n):
            bx, by, bz = coords[j]
            if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 <= th2:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _ca_by_uni(geom):
    m = {}
    for g in geom:
        if g["uniPos"] is not None and g["uniPos"] not in m:
            m[g["uniPos"]] = g
    return m


def _is_pathogenic(v):
    s = (v.get("consequence") or v.get("clinVar") or "").lower()
    return "pathogenic" in s or "deleterious" in s


def _residue_neighborhood(geom, center_uni, threshold=8.0):
    """[(uniPos, distance)] for residues within `threshold` Å of center_uni's Cα, nearest first."""
    cab = _ca_by_uni(geom)
    c = cab.get(center_uni)
    if not c:
        return []
    th2 = threshold * threshold
    out = [(uni, _math.sqrt(_d2(c, g))) for uni, g in cab.items()
           if uni != center_uni and _d2(c, g) <= th2]
    out.sort(key=lambda t: t[1])
    return out


def _dist(a, b):
    """Euclidean Cα distance exactly as the JS computes it (sqrt of dx²+dy²+dz², in that order)."""
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return _math.sqrt(dx * dx + dy * dy + dz * dz)


def _vlabel(v, pos=None):
    return "%s%s%s" % (v.get("wildType") or "", v["position"] if pos is None else pos, v.get("mutant") or "")


# Proximity radii (state.js proxPtmRadius / proxVarRadius, sliders 2-30 Å). They drive the residue
# report's PTM–variant proximity lists AND the PTM-centric tiers (modal.js passes them as d8 / d12).
PROX_SETTINGS = {"ptm": 8.0, "var": 12.0}


def _ptm_variant_proximity(ptms, variants, geom, d8=None, d12=None):
    """analysis.js computePtmVariantProximity. For each modelled PTM position: tier1 = every variant at
    that residue (any class); then over the other modelled positions within d12 Å (Cα), ONLY pathogenic /
    deleterious variants count: dist <= d8 -> tier2, else tier3. nearest = closest tier2/3 variant (0 Å /
    the first tier1 variant when tier1 exists). nearbyCount8A = |tier1| + |tier2|, pathCount8A = |tier2|.
    tier = 1 / 2 / 3 by the first non-empty list; PTMs with none are omitted."""
    d8 = PROX_SETTINGS["ptm"] if d8 is None else float(d8)
    d12 = PROX_SETTINGS["var"] if d12 is None else float(d12)
    cab = _ca_by_uni(geom)
    vbypos = {}
    for v in variants:
        vbypos.setdefault(v["position"], []).append(v)
    vpos = [(u, cab[u]) for u in cab if u in vbypos]      # only positions carrying variants matter
    np = _numpy()
    VX = np.array([(g["x"], g["y"], g["z"]) for _, g in vpos], dtype=np.float64).reshape(-1, 3) if np else None
    result, seen = {}, set()
    for ptm in ptms:
        pp = ptm["position"]
        if pp in seen:
            continue
        seen.add(pp)
        pc = cab.get(pp)
        if not pc:
            continue
        tier1 = list(vbypos.get(pp, []))
        tier2, tier3 = [], []
        nearest, nearest_var = float("inf"), None
        if VX is not None and len(vpos):        # coarse numpy pre-filter (+0.01 Å margin); exact test below
            dd = ((VX - (pc["x"], pc["y"], pc["z"])) ** 2).sum(1)
            cand = [vpos[k] for k in np.flatnonzero(dd <= (d12 + 0.01) ** 2).tolist()]
        else:
            cand = vpos
        for uni, g in cand:
            if uni == pp:
                continue
            dist = _dist(g, pc)
            if dist > d12:
                continue
            for v in vbypos[uni]:
                if not _is_pathogenic(v):
                    continue
                (tier2 if dist <= d8 else tier3).append((v, dist))
                if dist < nearest:
                    nearest, nearest_var = dist, _vlabel(v)
        if tier1 and nearest > 0:
            nearest, nearest_var = 0.0, _vlabel(tier1[0], pp)
        tier = 1 if tier1 else 2 if tier2 else 3 if tier3 else None
        if tier is not None:
            result[pp] = {
                "tier": tier, "tier1": tier1, "tier2": tier2, "tier3": tier3,
                "pathCount8A": len(tier2), "nearbyCount8A": len(tier1) + len(tier2),
                "d8": d8, "d12": d12,
                "nearestDist": None if nearest == float("inf") else nearest,
                "nearestVariant": nearest_var,
            }
    return result


def _residue_proximity(pos, ptms, variants, geom, ptm_radius=None, var_radius=None):
    """analysis.js computeResidueProximity, distances FROM residue `pos`: PTM positions (other than pos)
    within ptm_radius, and pathogenic/deleterious variants (other positions) within var_radius, each
    sorted by distance. -> {"ptms": [(ptm, dist)], "variants": [(variant, dist)]}"""
    ptm_radius = PROX_SETTINGS["ptm"] if ptm_radius is None else float(ptm_radius)
    var_radius = PROX_SETTINGS["var"] if var_radius is None else float(var_radius)
    cab = _ca_by_uni(geom)
    rc = cab.get(pos)
    if not rc:
        return {"ptms": [], "variants": []}
    near_p, seen = [], set()
    for p in ptms:
        q = p["position"]
        if q == pos or q in seen:
            continue
        seen.add(q)
        g = cab.get(q)
        if g is None:
            continue
        d = _dist(g, rc)
        if d <= ptm_radius:
            near_p.append((p, d))
    near_p.sort(key=lambda t: t[1])
    near_v = []
    for v in variants:
        if v["position"] == pos or not _is_pathogenic(v):
            continue
        g = cab.get(v["position"])
        if g is None:
            continue
        d = _dist(g, rc)
        if d <= var_radius:
            near_v.append((v, d))
    near_v.sort(key=lambda t: t[1])
    return {"ptms": near_p, "variants": near_v}


# ---- Deterministic PRNG shared with the extension: mulberry32(0x9e3779b9) (analysis.js / algorithms.js).
# Its state advances by a constant per draw, so draw i is a pure function of i: whole blocks of the
# stream are generated with numpy uint32 arithmetic (wrap-around == Math.imul / |0 semantics) and the
# results are bit-identical to the JS. The Fisher-Yates swaps themselves stay sequential (cheap lists).
_MB_INC = 0x6D2B79F5
HOTSPOT_RNG = "js"          # "js" = mulberry32, exact extension numbers; "python" = legacy random.Random


class _Mulberry32(object):
    """mulberry32 stream: next() -> float in [0,1); block(count) -> the next `count` draws (numpy array
    when numpy is available, else a list)."""

    def __init__(self, seed, np=None):
        self.a = seed & 0xFFFFFFFF
        self.n = 0                         # draws consumed
        self.np = np

    def next(self):
        self.n += 1
        a = (self.a + self.n * _MB_INC) & 0xFFFFFFFF
        t = ((a ^ (a >> 15)) * (a | 1)) & 0xFFFFFFFF
        t = ((t + (((t ^ (t >> 7)) * (t | 61)) & 0xFFFFFFFF)) & 0xFFFFFFFF) ^ t
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0

    def block(self, count):
        np = self.np
        if np is None:
            return [self.next() for _ in range(count)]
        i = np.arange(self.n + 1, self.n + count + 1, dtype=np.uint64)
        self.n += count
        a = ((np.uint64(self.a) + i * np.uint64(_MB_INC)) & np.uint64(0xFFFFFFFF)).astype(np.uint32)
        u15, u7, u14 = np.uint32(15), np.uint32(7), np.uint32(14)
        t = (a ^ (a >> u15)) * (a | np.uint32(1))                  # uint32 arrays wrap: Math.imul
        t = (t + (t ^ (t >> u7)) * (t | np.uint32(61))) ^ t
        return (t ^ (t >> u14)).astype(np.float64) / 4294967296.0


def _fy_draws(rng, k, perms):
    """r-indices for `perms` successive JS Fisher-Yates shuffles of a length-k array (for i = k-1..1:
    r = floor(rand() * (i+1))) -> list of per-shuffle lists."""
    if k < 2 or perms <= 0:
        return [[] for _ in range(max(perms, 0))]
    np = rng.np
    u = rng.block(perms * (k - 1))
    if np is not None:
        mult = np.arange(k, 1, -1, dtype=np.float64)
        return np.floor(u.reshape(perms, k - 1) * mult).astype(np.int64).tolist()
    mult = list(range(k, 1, -1))
    return [[int(u[p * (k - 1) + j] * mult[j]) for j in range(k - 1)] for p in range(perms)]


def _apply_fy(arr, rs):
    """One JS shuffleInPlace with precomputed r-indices (in place)."""
    i = len(arr) - 1
    for r in rs:
        arr[i], arr[r] = arr[r], arr[i]
        i -= 1


def _fy_positions(rng, k, perms):
    """numpy form of _fy_draws + _apply_fy: the same RNG draws (same stream consumption), returned as a
    (perms x k) int32 array whose row p is arange(k) after JS shuffle p — a Fisher-Yates pass with fixed
    draws is a fixed position permutation, so shuffling ANY array `a` with those draws gives a[rows[p]].
    The k-1 swap steps are applied to all `perms` rows at once (one vectorised swap per step) instead
    of perms x k Python swaps; the result is identical. Requires rng.np."""
    np = rng.np
    if k < 2 or perms <= 0:
        pos = np.empty((max(perms, 0), k), dtype=np.int32)
        pos[:] = np.arange(k, dtype=np.int32)
        return pos
    u = rng.block(perms * (k - 1))
    R = np.floor(u.reshape(perms, k - 1) * np.arange(k, 1, -1, dtype=np.float64)).astype(np.int64)
    # transposed (k x perms) working copy: position i of every shuffle is one contiguous row
    T = np.empty((k, perms), dtype=np.int32)
    T[:] = np.arange(k, dtype=np.int32)[:, None]
    flat = T.reshape(-1)
    cols = np.arange(perms, dtype=np.int64)
    RT = np.ascontiguousarray(R.T)
    for j in range(k - 1):
        i = k - 1 - j
        fi = RT[j] * perms + cols                  # element (r_p, p) of every shuffle p
        a = T[i].copy()
        T[i] = flat[fi]
        flat[fi] = a
    return np.ascontiguousarray(T.T)


# ---- Contact hubs (analysis.js betweennessHubs) ----
HUB_EXACT = True            # exact Brandes over ALL sources (numpy, level-synchronous); else sampled
HUB_SAMPLE_SOURCES = 400    # sampled fallback (no numpy, or HUB_EXACT off): sources drawn per chain


def _brandes_numpy(adj, np):
    """Unnormalised Brandes betweenness (unweighted, undirected; bc[w] += delta[w] for w != s) over ALL
    sources, exactly what analysis.js computes (to ~1e-14 relative: only the summation order differs).
    A batch of sources is searched level-synchronously with vectorised gathers over a padded neighbour
    table (padding -> a per-source sentinel slot that is never 'unvisited'). Each level keeps its BFS-DAG
    edges (frontier row -> next-level node), so sigma (forward) and delta (backward) are plain bincounts
    over those edges — no per-level re-expansion. ~1 s for a 2000-residue chain."""
    n = len(adj)
    D = max(1, max(len(a) for a in adj))
    m1 = n + 1                                    # per-source block: n nodes + 1 sentinel slot
    NB = np.full((n, D), n, dtype=np.int64)
    for i, a in enumerate(adj):
        NB[i, :len(a)] = a
    bc = np.zeros(n)
    B = max(1, min(n, 800000 // m1))
    big = np.int32(1 << 30)
    for s0 in range(0, n, B):
        src = np.arange(s0, min(n, s0 + B), dtype=np.int64)
        b = len(src)
        N = b * m1
        dist = np.full(N, -1, dtype=np.int32)
        dist[np.arange(b, dtype=np.int64) * m1 + n] = big
        sigma = np.zeros(N)
        owner = np.empty(N, dtype=np.int64)
        start = np.arange(b, dtype=np.int64) * m1 + src
        dist[start] = 0
        sigma[start] = 1.0
        levels, edges = [start], []
        front, d = start, 0
        while True:
            u = front % m1
            nb = (front - u)[:, None] + NB[u]
            rows, cols = np.nonzero(dist[nb] < 0)        # DAG edges: frontier row -> unvisited node
            if not rows.size:
                break
            cand = nb[rows, cols]
            ar = np.arange(cand.size, dtype=np.int64)
            owner[cand] = ar                             # one representative entry per next-level node
            rep = owner[cand]
            is_rep = rep == ar
            nxt = cand[is_rep]
            gid = (np.cumsum(is_rep) - 1)[rep]           # entry -> index in nxt
            dist[nxt] = d + 1
            sigma[nxt] = np.bincount(gid, weights=sigma[front][rows], minlength=nxt.size)
            edges.append((rows, gid))
            levels.append(nxt)
            front, d = nxt, d + 1
        delta = np.zeros(N)
        for lv in range(len(levels) - 1, 0, -1):
            W, V = levels[lv], levels[lv - 1]
            rows, gid = edges[lv - 1]
            coef = (1.0 + delta[W]) / sigma[W]
            delta[V] = sigma[V] * np.bincount(rows, weights=coef[gid], minlength=V.size)
        delta[start] = 0.0
        bc += delta.reshape(b, m1)[:, :n].sum(0)
    return bc


def _brandes_python(adj, sources):
    n = len(adj)
    bc = [0.0] * n
    for s in sources:
        stack, pred = [], [[] for _ in range(n)]
        sigma = [0.0] * n
        sigma[s] = 1.0
        dist = [-1] * n
        dist[s] = 0
        queue, qi = [s], 0
        while qi < len(queue):
            v = queue[qi]
            qi += 1
            stack.append(v)
            for w in adj[v]:
                if dist[w] < 0:
                    dist[w] = dist[v] + 1
                    queue.append(w)
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = [0.0] * n
        for idx in range(len(stack) - 1, -1, -1):
            w = stack[idx]
            for v in pred[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
            if w != s:
                bc[w] += delta[w]
    return bc


def _betweenness_hubs(geom, threshold=8.0, exact=None):
    """analysis.js betweennessHubs for ONE chain: nodes = every Cα of the chain (mapped or not), edges =
    Cα-Cα <= threshold; Brandes betweenness * 1/((n-1)(n-2)); population z-score; strong z >= 3,
    moderate z >= 2 (bc > 0 only), reported at UniProt positions. Exact (all sources) by default via
    numpy; without numpy (or exact=False) betweenness is estimated from HUB_SAMPLE_SOURCES sources.
    Chains with no mapped residue (partner proteins) are skipped. -> {uniPos: tier}"""
    if not any(g["uniPos"] is not None for g in geom):
        return {}
    n = len(geom)
    if n < 8 or n > HUB_MAX_CA:
        return {}
    adj = _adjacency([(g["x"], g["y"], g["z"]) for g in geom], threshold)
    exact = HUB_EXACT if exact is None else exact
    np = _numpy()
    if exact and np is not None:
        bc = _brandes_numpy(adj, np).tolist()
    else:
        sources = range(n)
        if n > HUB_SAMPLE_SOURCES:
            sources = _random.Random(HOTSPOT_SEED).sample(range(n), HUB_SAMPLE_SOURCES)
        bc = _brandes_python(adj, sources)
    scale = 1.0 / ((n - 1) * (n - 2))
    bc = [b * scale for b in bc]
    mean = sum(bc) / n
    std = _math.sqrt(sum((b - mean) * (b - mean) for b in bc) / n)
    if not std > 0:
        return {}
    out = {}
    for i in range(n):
        if bc[i] <= 0:
            continue
        z = (bc[i] - mean) / std
        tier = "strong" if z >= HUB_Z_STRONG else "moderate" if z >= HUB_Z_MODERATE else None
        if not tier:
            continue
        uni = geom[i]["uniPos"]
        if uni is None:
            continue
        if out.get(uni) != "strong":
            out[uni] = tier
    return out


def _bh_fdr(centres):
    """In-place Benjamini-Hochberg q-values; `centres` sorted ascending by 'p'."""
    n = len(centres)
    for i, c in enumerate(centres):
        c["q"] = min(c["p"] * n / (i + 1), 1.0)
    for i in range(n - 2, -1, -1):
        centres[i]["q"] = min(centres[i]["q"], centres[i + 1]["q"])


HOTSPOT_CHUNK = 128          # permutations scored per numpy matrix product
_HOTSPOT_PURE = False        # tests: force the pure-Python permutation loops


def _perm_ge_casecontrol(neigh, labels, local, rng, np=None):
    """Case-control null: per permutation, shuffle the path/benign labels (ONE array shuffled
    cumulatively, as the JS does) and count, for every item, how often its permuted pathogenic
    neighbour count reaches the observed `local`. `rng` is a _Mulberry32 (exact JS draws) or a
    random.Random (legacy). numpy and pure Python give IDENTICAL counts; numpy scores a block of
    permutations at once as a (perms x k) @ (k x k) neighbourhood product."""
    k = len(labels)
    perm = labels[:]
    js = isinstance(rng, _Mulberry32)
    if js and np is not None and rng.np is not None:
        # vectorised: all shuffles as position permutations (_fy_positions), applied cumulatively to ONE
        # label vector exactly as the sequential JS shuffles do -> identical rows, identical counts
        pos = _fy_positions(rng, k, HOTSPOT_PERMUTATIONS)
        A = np.zeros((k, k), dtype=np.float32)
        for i in range(k):
            A[neigh[i], i] = 1.0
        loc = np.asarray(local, dtype=np.float32)
        ge = np.zeros(k, dtype=np.int64)
        cur = np.asarray(perm, dtype=np.float32)
        done = 0
        while done < HOTSPOT_PERMUTATIONS:
            b = min(HOTSPOT_CHUNK, HOTSPOT_PERMUTATIONS - done)
            rows = np.empty((b, k), dtype=np.float32)
            for r in range(b):
                cur = cur[pos[done + r]]
                rows[r] = cur
            ge += ((rows @ A) >= loc).sum(0)
            done += b
        return ge.tolist()
    draws = _fy_draws(rng, k, HOTSPOT_PERMUTATIONS) if js else None

    def shuffle(p):
        if js:
            _apply_fy(perm, draws[p])
        else:
            rng.shuffle(perm)
    if np is None:
        ge = [0] * k
        for p in range(HOTSPOT_PERMUTATIONS):
            shuffle(p)
            for i in range(k):
                if sum(1 for j in neigh[i] if perm[j]) >= local[i]:
                    ge[i] += 1
        return ge
    A = np.zeros((k, k), dtype=np.float32)         # A[j, i] = 1 if j in neigh[i] (self included)
    for i in range(k):
        A[neigh[i], i] = 1.0
    loc = np.asarray(local, dtype=np.float32)
    ge = np.zeros(k, dtype=np.int64)
    done = 0
    while done < HOTSPOT_PERMUTATIONS:
        b = min(HOTSPOT_CHUNK, HOTSPOT_PERMUTATIONS - done)
        rows = np.empty((b, k), dtype=np.float32)
        for r in range(b):
            shuffle(done + r)
            rows[r] = perm
        ge += ((rows @ A) >= loc).sum(0)             # small integer counts: exact in float32
        done += b
    return ge.tolist()


def _perm_ge_spatial(neigh, path_idx, obs, n, rng, np=None, partner_den=None):
    """Spatial-placement null: per permutation, place len(path_idx) pathogenic labels uniformly (partial
    Fisher-Yates on a persistent pool: r = i + floor(rand()*(n-i))) and count, per observed centre, how
    often partnerDen + its simulated pathogenic neighbour count reaches the observed density. Same RNG
    draw sequence in both paths -> identical results; numpy scores blocks of permutations as one product."""
    m = len(path_idx)
    pool = list(range(n))
    pden = partner_den or [0] * n
    js = isinstance(rng, _Mulberry32)
    if js:
        u = rng.block(HOTSPOT_PERMUTATIONS * m)
        span = [n - i for i in range(m)]
        if np is not None:
            R = (np.floor(u.reshape(HOTSPOT_PERMUTATIONS, m) * np.asarray(span, dtype=np.float64))
                 .astype(np.int64) + np.arange(m, dtype=np.int64)).tolist() if m else [[]] * HOTSPOT_PERMUTATIONS
        else:
            R = [[i + int(u[p * m + i] * span[i]) for i in range(m)] for p in range(HOTSPOT_PERMUTATIONS)]

    def place(p):
        if js:
            rs = R[p]
            for i in range(m):
                r = rs[i]
                pool[i], pool[r] = pool[r], pool[i]
        else:
            for i in range(m):
                r = i + rng.randrange(n - i)
                pool[i], pool[r] = pool[r], pool[i]
    if np is None:
        ge = {ci: 0 for ci in path_idx}
        for p in range(HOTSPOT_PERMUTATIONS):
            place(p)
            sim = [0] * n
            for i in range(m):
                sim[pool[i]] = 1
            for ci in path_idx:
                if pden[ci] + sum(1 for j in neigh[ci] if sim[j]) >= obs[ci]:
                    ge[ci] += 1
        return ge
    A = np.zeros((n, m), dtype=np.float32)          # column c = neighbourhood of centre path_idx[c]
    for c, ci in enumerate(path_idx):
        A[neigh[ci], c] = 1.0
    ob = np.asarray([obs[ci] - pden[ci] for ci in path_idx], dtype=np.float32)
    ge = np.zeros(m, dtype=np.int64)
    done = 0
    while done < HOTSPOT_PERMUTATIONS:
        b = min(HOTSPOT_CHUNK, HOTSPOT_PERMUTATIONS - done)
        sim = np.zeros((b, n), dtype=np.float32)
        for r in range(b):
            place(done + r)
            sim[r, pool[:m]] = 1.0
        ge += ((sim @ A) >= ob).sum(0)
        done += b
    return {ci: int(ge[c]) for c, ci in enumerate(path_idx)}


def _hotspot_rng(mode=None):
    mode = mode or HOTSPOT_RNG
    if mode == "python":
        return _random.Random(HOTSPOT_SEED)
    return _Mulberry32(HOTSPOT_SEED, None if _HOTSPOT_PURE else _numpy())


_RE_PATH_CONS = re.compile(r"pathogenic|deleterious", re.I)
_RE_BENIGN_CONS = re.compile(r"benign", re.I)


def _hotspot_classes(variants):
    """Ordered (first-appearance, like the JS Sets) pathogenic positions and benign-only positions."""
    path = {}
    for v in variants:
        if _RE_PATH_CONS.search(v.get("consequence") or ""):
            path.setdefault(v["position"], None)
    benign = {}
    for v in variants:
        if v["position"] not in path and _RE_BENIGN_CONS.search(v.get("consequence") or ""):
            benign.setdefault(v["position"], None)
    return list(path), list(benign)


def _compute_hotspots(geom, variants, threshold=8.0, partners=None, rng=None, classes=None):
    """analysis.js computeHotspots for ONE chain. Case-control label permutation when >= 2 benign-only
    positions exist (else a spatial-placement null over every Cα of the chain). `partners` = resolved
    partner-protein classified residues [(is_pathogenic, geom_entry)] (other chains of a complex): in
    case-control they join the pool when within `threshold` of one of our items; in the spatial null the
    pathogenic ones add a fixed density. They are never reported. Exact JS numbers with rng 'js'
    (mulberry32, default). Returns {uniPos: 'strong'|'moderate'|'weak'}."""
    path, benign = classes or _hotspot_classes(variants)
    if len(path) < HOTSPOT_MIN_PATH:
        return {}
    if not any(g["uniPos"] is not None for g in geom):
        return {}
    partners = partners or []
    np = None if _HOTSPOT_PURE else _numpy()
    th2 = threshold * threshold

    if len(benign) >= HOTSPOT_MIN_BENIGN:  # case-control
        cab = _ca_by_uni(geom)
        items = [(p, True, cab[p], True) for p in path if p in cab] + \
                [(p, False, cab[p], True) for p in benign if p in cab]
        ours = len(items)
        if ours == 0:
            return {}
        for is_path, g in partners:
            if any(_d2(g, it[2]) <= th2 for it in items[:ours]):
                items.append((None, bool(is_path), g, False))
        k = len(items)
        if k < HOTSPOT_MIN_PATH + 1:
            return {}
        adj = _adjacency([(it[2]["x"], it[2]["y"], it[2]["z"]) for it in items], threshold)
        neigh = [[i] + adj[i] for i in range(k)]  # include self, as in the JS
        labels = [it[1] for it in items]
        local = [sum(1 for j in neigh[i] if labels[j]) for i in range(k)]
        base = sum(1 for l in labels if l) / k
        ge = _perm_ge_casecontrol(neigh, labels, local, rng or _hotspot_rng(), np)
        centres = []
        for i in range(k):
            if not items[i][1] or not items[i][3]:
                continue
            tot = len(neigh[i])
            frac = local[i] / tot if tot else 0
            centres.append({"pos": items[i][0], "p": (ge[i] + 1) / (HOTSPOT_PERMUTATIONS + 1),
                            "local": local[i], "frac": frac, "er": frac / max(base, 1e-9)})
        centres.sort(key=lambda c: c["p"])
        _bh_fdr(centres)
        out = {}
        for c in centres:
            if c["q"] <= 0.10 and c["local"] >= 3 and c["er"] >= 2.0:
                out[c["pos"]] = "strong"
            elif c["q"] <= 0.25 and c["local"] >= 2 and c["er"] >= 1.5:
                out[c["pos"]] = "moderate"
            elif c["q"] <= 0.25 and c["local"] >= 2 and c["frac"] >= 0.40:
                out[c["pos"]] = "weak"
        return out

    # spatial-placement null over EVERY Cα of the chain (no benign controls needed)
    universe = geom
    n = len(universe)
    if n < 10 or n > 6000:
        return {}
    idx_by_uni = {}
    for i, u in enumerate(universe):
        if u["uniPos"] is not None and u["uniPos"] not in idx_by_uni:
            idx_by_uni[u["uniPos"]] = i
    path_idx = [idx_by_uni[p] for p in path if p in idx_by_uni]
    m = len(path_idx)
    if m < HOTSPOT_MIN_PATH:
        return {}
    neigh = _adjacency([(u["x"], u["y"], u["z"]) for u in universe], threshold)
    pden = [0] * n
    for is_path, g in partners:
        if not is_path:
            continue
        for i in range(n):
            if _d2(g, universe[i]) <= th2:
                pden[i] += 1
    is_path = [0] * n
    for i in path_idx:
        is_path[i] = 1
    obs = [pden[i] + sum(1 for j in neigh[i] if is_path[j]) for i in range(n)]
    ge = _perm_ge_spatial(neigh, path_idx, obs, n, rng or _hotspot_rng(), np, pden)
    centres = [{"ci": ci, "p": (ge[ci] + 1) / (HOTSPOT_PERMUTATIONS + 1), "den": obs[ci]} for ci in path_idx]
    centres.sort(key=lambda c: c["p"])
    _bh_fdr(centres)
    out = {}
    for c in centres:
        uni = universe[c["ci"]]["uniPos"]
        if uni is None:
            continue
        if c["q"] <= 0.10 and c["den"] >= 3:
            out[uni] = "strong"
        elif c["q"] <= 0.25 and c["den"] >= 2:
            out[uni] = "moderate"
    return out


def _chains_of(geom):
    chains = {}
    for g in geom:
        chains.setdefault(g["chain"], []).append(g)
    return chains


def _per_chain_merge(geom, fn, rank=None, by_chain=None):
    """Run a per-chain analysis `fn(geom_of_chain)` and merge tiers across chains (strongest wins, by
    `rank`). Chains without any mapped residue are skipped. `by_chain` (dict) receives each chain's map."""
    rank = rank or _TIER_RANK
    merged = {}
    for ch, gs in _chains_of(geom).items():
        if not any(g["uniPos"] is not None for g in gs):
            continue
        res = fn(gs)
        if by_chain is not None and res:
            by_chain[ch] = res
        for uni, tier in res.items():
            if uni is None:
                continue
            if uni not in merged or rank.get(tier, 0) > rank.get(merged[uni], 0):
                merged[uni] = tier
    return merged


_HUB_RANK = {"strong": 2, "moderate": 1}


def _am_by_pos(am_map):
    """Index AlphaMissense by position once: {pos: [(wt, mut, score)] highest-first}. Avoids
    rescanning the whole substitution map (tens of thousands of keys) on every residue report."""
    out = {}
    match = _RE_AM_KEY.match
    for key, score in am_map.items():
        mt = match(key)
        if mt:
            out.setdefault(int(mt.group(2)), []).append((mt.group(1), mt.group(3), score))
    for p in out:
        out[p].sort(key=lambda t: -t[2])
    return out


def _am_profile(am_map, pos):
    """Per-substitution AlphaMissense scores at a position: [(wt, mut, score)], highest first."""
    return _am_by_pos(am_map).get(pos, [])


# Constraint-pocket (port of algorithms.js UFVPocket.computePockets). Getis-Ord Gi* on the
# residual of per-residue AlphaMissense mean after regressing out structural burial (coordination
# number) — residues more evolutionarily constrained than their burial explains, clustered in 3D.
# Spatial weights are gated by the AlphaFold PAE (exp(-PAE/10)) for AlphaFold models.
PRISM_HSE = 13.0     # Å coordination-sphere radius (burial)
PRISM_CUT = 13.0     # Å spatial-weight neighbourhood
PRISM_SIGMA = 5.0    # Å Gaussian bandwidth
PRISM_TAU = 10.0     # Å PAE gate scale
PRISM_MAXN = 3000
AM_AAS = "ACDEFGHIKLMNPQRSTVWY"
POCKET_Q_DEFAULT = 0.10                 # display threshold (modal.js sensitivity slider default)
POCKET_SETTINGS = {"q": POCKET_Q_DEFAULT}   # current sensitivity (0.01-0.40); recolours without recompute
PAE_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-predicted_aligned_error_v{}.json"
PAE_MAXN = 2500
_PAE_CACHE = {}      # (ID, version) -> {"n", "data" (numpy float32 n*n)} or None (LRU, _PAE_CACHE_MAX matrices)
_PAE_CACHE_MAX = 3


def _perm_count(n):
    return 999 if n <= 600 else 599 if n <= 1200 else 399


def _parse_pae(raw, np):
    """api.js parsePae: [{predicted_aligned_error: [[...]]}] (or the object), or legacy triplets
    {residue1, residue2, distance}. -> {"n": n, "data": float32 n x n} or None (n == 0 or > 2500).
    The common list-of-rows layout is parsed straight from the text with numpy (several times faster
    than json.loads on the ~20 MB file of a 2000-residue model)."""
    if not raw or np is None:
        return None
    txt = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else str(raw)
    k = txt.find('"predicted_aligned_error"')
    if k >= 0:
        a = txt.find("[[", k)
        b = txt.find("]]", a + 2) if a >= 0 else -1
        if a >= 0 and b > a:
            n = txt.count("]", a, b) + 1
            if n > PAE_MAXN:
                return None
            try:
                arr = np.array(json.loads(txt[a:b + 2]), dtype=np.float32)
            except (ValueError, TypeError):
                arr = None
            if arr is not None and arr.shape == (n, n):
                return {"n": n, "data": np.nan_to_num(arr, nan=0.0)}
    try:
        j = json.loads(txt)
    except Exception:
        return None
    obj = j[0] if isinstance(j, list) and j else j
    if not isinstance(obj, dict):
        return None
    m = obj.get("predicted_aligned_error")
    if isinstance(m, list) and m and isinstance(m[0], list):
        n = len(m)
        if n == 0 or n > PAE_MAXN:
            return None
        data = np.zeros((n, n), dtype=np.float32)
        for i, row in enumerate(m):
            r = np.asarray([(float(x) if isinstance(x, (int, float)) else 0.0) for x in row[:n]], dtype=np.float32)
            data[i, :len(r)] = r
        return {"n": n, "data": data}
    r1, r2, d = obj.get("residue1"), obj.get("residue2"), obj.get("distance")
    if isinstance(r1, list) and isinstance(r2, list) and isinstance(d, list):
        n = max([0] + [int(x) for x in r1] + [int(x) for x in r2])
        if n == 0 or n > PAE_MAXN:
            return None
        data = np.zeros((n, n), dtype=np.float32)
        for a_, b_, v in zip(r1, r2, d):
            data[int(a_) - 1, int(b_) - 1] = float(v) if isinstance(v, (int, float)) else 0.0
        return {"n": n, "data": data}
    return None


def _pae_ok(b):
    if b.lstrip()[:1] not in (b"[", b"{", b"n"):
        raise ValueError("not JSON")


def _pae_matrix(af_id, version=None):
    """AlphaFold PAE for model AF-<af_id>-F1 (api.js getPaeMatrix): tries [version, 4, 3] and keeps the
    first that parses. Disk-cached + memoised per id. Network: call off the UI thread. None if absent."""
    np = _numpy()
    if not af_id or np is None:
        return None
    key = af_id.upper()
    vers = []
    for v in (version, 4, 3):
        try:
            v = int(v) if v else None
        except (TypeError, ValueError):
            v = None
        if v and v not in vers:
            vers.append(v)
    for v in vers:
        if (key, v) in _PAE_CACHE:
            if _PAE_CACHE[(key, v)] is not None:
                _PAE_CACHE[(key, v)] = _PAE_CACHE.pop((key, v))     # most recently used last
                return _PAE_CACHE[(key, v)]
            continue
        raw = _fetch_bytes(PAE_URL.format(af_id, v), cache=True, cache_404=b"null", validate=_pae_ok)
        if raw is None:
            continue                       # network error: not memoised, retried next time
        pae = _parse_pae(raw, np) if raw.strip() != b"null" else None
        del raw
        _PAE_CACHE[(key, v)] = pae
        if pae is not None:
            # bounded: an n x n float32 matrix is up to 25 MB (2500 residues); keep the most recent few
            # (re-parsed from the disk cache if needed again). Misses (None) are tiny and all kept.
            held = [k for k, m in _PAE_CACHE.items() if m is not None]
            for k in held[:max(0, len(held) - _PAE_CACHE_MAX)]:
                del _PAE_CACHE[k]
            return pae
    return None


def _af_model_info(obj):
    """(AlphaFold accession id, version) for an AlphaFold model object loaded through the plugin — the
    isoform id for isoform models — else None (PAE gating applies to AlphaFold models only)."""
    st = (_STATE.get("obj_struct") or {}).get(obj) or {}
    mp = UFV_MAPS.get(obj) or {}
    if st.get("source") != "AlphaFold" and not mp.get("afId"):
        return None
    af_id = mp.get("afId") or st.get("isoform") or mp.get("isoform") or \
        (_STATE.get("obj_uid") or {}).get(obj) or mp.get("uid")
    if not af_id:
        return None
    ver = mp.get("afVersion")
    if ver is None:
        m = re.search(r"_v(\d+)\.(?:pdb|cif|bcif)", str(st.get("url") or ""), re.I)
        ver = int(m.group(1)) if m else None
    return af_id, ver


def _am_mean_wt(am_map, sequence, pos):
    """algorithms.js amMean: mean of the finite AM scores of the 19 substitutions WT->X at `pos`, WT from
    the canonical sequence (summed in AM_AAS order, as the JS does). None if none."""
    if not sequence or pos < 1 or pos > len(sequence):
        return None
    wt = sequence[pos - 1]
    s, c = 0.0, 0
    for mt in AM_AAS:
        if mt == wt:
            continue
        sc = am_map.get("%s%d%s" % (wt, pos, mt))
        if isinstance(sc, (int, float)) and sc == sc and abs(sc) != float("inf"):
            s += sc
            c += 1
    return s / c if c else None


def _loess(X, Y, np, f=0.3):
    """algorithms.js loess: 1-D local linear fit, tricube weights, bandwidth = k-th smallest |x - x_i|."""
    n = len(X)
    k = max(2, min(n, int(_math.floor(f * n + 0.5))))       # Math.round
    pred = np.empty(n)
    kk = min(k - 1, n - 1)
    if LOESS_BLOCK and n > 1:
        return _loess_block(X, Y, np, n, kk, pred)
    for i in range(n):
        d = np.abs(X - X[i])
        h = np.partition(d, kk)[kk] or 1e-9
        u = d / h
        w = np.where(u < 1, (1 - u * u * u) ** 3, 0.0)
        sw = w.sum(); swx = (w * X).sum(); swy = (w * Y).sum()
        swxx = (w * X * X).sum(); swxy = (w * X * Y).sum()
        denom = sw * swxx - swx * swx
        if abs(denom) < 1e-12:
            pred[i] = swy / sw if sw else Y[i]
        else:
            b = (sw * swxy - swx * swy) / denom
            pred[i] = (swy - b * swx) / sw + b * X[i]
    return pred


LOESS_BLOCK = True      # False: the original one-row-at-a-time loop (reference for the tests)


def _loess_block(X, Y, np, n, kk, pred, blk=256):
    """_loess for blocks of rows at once, bit-identical to the per-row loop: the elementwise steps are the
    same IEEE operations on the same operands, np.partition's k-th value is exact, the power is still
    evaluated on one contiguous row at a time (same code path / tail split as the 1-D call), the row sums
    reduce each contiguous row exactly as the 1-D sum does, and the final per-row formula is the
    original scalar code."""
    for i0 in range(0, n, blk):
        i1 = min(i0 + blk, n)
        D = np.abs(X[None, :] - X[i0:i1, None])
        h = np.partition(D, kk, axis=1)[:, kk]
        h = np.where(h != 0, h, 1e-9)                # `... or 1e-9`
        U = D / h[:, None]
        T = 1 - U * U * U
        for r in range(i1 - i0):
            T[r] = T[r].copy() ** 3                  # a fresh 1-D array, exactly like the per-row code's
        Wt = np.where(U < 1, T, 0.0)
        sw = Wt.sum(1)
        swx = (Wt * X).sum(1)
        swy = (Wt * Y).sum(1)
        swxx = (Wt * X * X).sum(1)
        swxy = (Wt * X * Y).sum(1)
        for r in range(i1 - i0):
            i = i0 + r
            s0, s1, s2, s3, s4 = sw[r], swx[r], swy[r], swxx[r], swxy[r]
            denom = s0 * s3 - s1 * s1
            if abs(denom) < 1e-12:
                pred[i] = s2 / s0 if s0 else Y[i]
            else:
                b = (s0 * s4 - s1 * s2) / denom
                pred[i] = (s2 - b * s1) / s0 + b * X[i]
    return pred


def _compute_pockets(geom, am_map, sequence, pae=None, max_n=PRISM_MAXN):
    """algorithms.js computePockets. geom = every Cα of the model in model order (partner chains with
    uniPos None are burial context only); scored residues = mapped residues with an AlphaMissense mean.
    PAE (AlphaFold, {"n","data"}) gates the spatial weights by the residue's absolute model-order index.
    Returns {"byPos": {uniPos: {"cat","p","q","score"}} for EVERY positive-residual candidate (the
    display threshold is applied when colouring), "reason", "hasPae", "n"} — identical to the extension
    (mulberry32 permutations)."""
    empty = {"byPos": {}, "reason": None, "hasPae": False, "n": 0}
    np = _numpy()
    if np is None:
        return dict(empty, reason="numpy is required for the constraint-pocket analysis.")
    if not am_map:
        return dict(empty, reason="No AlphaMissense data for this protein.")
    if not geom:
        return dict(empty, reason="No modelled residues.")
    am_cache = {}
    recs = []
    for i, g in enumerate(geom):
        u = g["uniPos"]
        if u is None:
            continue
        if u not in am_cache:
            am_cache[u] = _am_mean_wt(am_map, sequence, u)
        if am_cache[u] is None:
            continue
        recs.append((i, g, am_cache[u]))
    n = len(recs)
    if n < 12:
        return dict(empty, reason="Too few residues with AlphaMissense scores.")
    if n > max_n:
        return dict(empty, reason="Protein too large for the analysis (%d > %d residues)." % (n, max_n))
    C = np.array([(g["x"], g["y"], g["z"]) for g in geom], float)
    ctx_idx = np.array([i for i, _, _ in recs])
    P = C[ctx_idx]
    am = np.array([a for _, _, a in recs], float)

    by_cr = {}
    for k, (_, g, _) in enumerate(recs):
        by_cr[(g["chain"], g["resi"])] = k                     # later copies win, like Map.set
    dirs = np.zeros((n, 3))
    has_dir = np.zeros(n, bool)
    for k, (_, g, _) in enumerate(recs):
        pv = by_cr.get((g["chain"], g["resi"] - 1))
        nx = by_cr.get((g["chain"], g["resi"] + 1))
        v = None
        if pv is not None and nx is not None:
            v = 2 * P[k] - P[pv] - P[nx]
        elif pv is not None:
            v = P[k] - P[pv]
        elif nx is not None:
            v = P[k] - P[nx]
        if v is not None:
            L = _math.hypot(float(v[0]), float(v[1]), float(v[2])) or 1.0
            dirs[k] = v / L
            has_dir[k] = True

    hse2 = PRISM_HSE ** 2
    cn = np.zeros(n); hse_up = np.zeros(n); hse_down = np.zeros(n)
    blk = 256
    # per-coordinate (b, nctx) difference planes instead of one strided (b, nctx, 3) array: the same
    # subtractions / products / sums in the same order (identical values), on contiguous memory
    Cx, Cy, Cz = (np.ascontiguousarray(C[:, a]) for a in range(3))
    Px, Py, Pz = (np.ascontiguousarray(P[:, a]) for a in range(3))
    for i0 in range(0, n, blk):
        i1 = min(i0 + blk, n)
        dx = Cx[None, :] - Px[i0:i1, None]               # (b, nctx)
        dy = Cy[None, :] - Py[i0:i1, None]
        dz = Cz[None, :] - Pz[i0:i1, None]
        d2 = dx * dx + dy * dy + dz * dz
        within = d2 <= hse2
        within[np.arange(i1 - i0), ctx_idx[i0:i1]] = False
        cn[i0:i1] = within.sum(1)
        proj = dx * dirs[i0:i1, None, 0] + dy * dirs[i0:i1, None, 1] + dz * dirs[i0:i1, None, 2]
        up = (within & (proj > 0)).sum(1)
        hd = has_dir[i0:i1]
        hse_up[i0:i1] = np.where(hd, up, 0)
        hse_down[i0:i1] = np.where(hd, within.sum(1) - up, 0)

    cut2 = PRISM_CUT ** 2
    two_sig2 = 2 * PRISM_SIGMA * PRISM_SIGMA
    paeN = int(pae["n"]) if pae else 0
    W = np.zeros((n, n))
    # the usual AlphaFold case: scored residues == model residues 0..n-1, all inside the PAE -> plain
    # slices of the PAE matrix instead of np.ix_ gathers (same elements, same shapes, same values)
    contig = bool(paeN) and int(ctx_idx[0]) == 0 and int(ctx_idx[-1]) == n - 1 and n <= paeN
    for i0 in range(0, n, blk):
        i1 = min(i0 + blk, n)
        dx = Px[i0:i1, None] - Px[None, :]
        dy = Py[i0:i1, None] - Py[None, :]
        dz = Pz[i0:i1, None] - Pz[None, :]
        d2 = dx * dx + dy * dy + dz * dz
        w = np.where(d2 <= cut2, np.exp(-d2 / two_sig2), 0.0)
        if paeN:
            if contig:
                D = pae["data"]
                sub = D[i0:i1, :n].astype(np.float64)
                subT = D[:n, i0:i1].astype(np.float64).T
                w = w * np.exp(-(0.5 * (sub + subT)) / PRISM_TAU)
                W[i0:i1] = w
                continue
            ci = ctx_idx[i0:i1]
            gate = np.ones_like(w)
            okr = ci < paeN
            okc = ctx_idx < paeN
            if okr.any() and okc.any():
                D = pae["data"]
                sub = D[np.ix_(ci[okr], ctx_idx[okc])].astype(np.float64)
                subT = D[np.ix_(ctx_idx[okc], ci[okr])].astype(np.float64).T
                gate[np.ix_(np.flatnonzero(okr), np.flatnonzero(okc))] = np.exp(-(0.5 * (sub + subT)) / PRISM_TAU)
            w = w * gate
        W[i0:i1] = w
    np.fill_diagonal(W, 0.0)
    W[W <= 1e-6] = 0.0

    rp = am - _loess(cn, am, np, 0.3)
    np.fill_diagonal(W, 1.0)          # M = W + I in place (diagonal 0 + 1, off-diagonal w + 0 == w): no 2nd n x n
    M = W                             # Gi*_i = rp[i] + Σ w_ij rp[j]
    gi_obs = M.dot(rp)

    perm = _perm_count(n)
    rng = _Mulberry32(HOTSPOT_SEED, np)
    pos = _fy_positions(rng, n, perm)  # == _fy_draws + sequential _apply_fy on one index list (identical)
    idx = np.arange(n)
    ge = np.zeros(n, dtype=np.int64)
    chunk = 64
    for p0 in range(0, perm, chunk):
        p1 = min(perm, p0 + chunk)
        cols = np.empty((n, p1 - p0), dtype=np.int64)
        for p in range(p0, p1):
            idx = idx[pos[p]]
            cols[:, p - p0] = idx
        ge += (M.dot(rp[cols]) >= gi_obs[:, None]).sum(1)

    cn_sorted = np.sort(cn)
    med_cn = float(cn_sorted[n >> 1]) if n else 0.0
    best = {}  # uniPos -> (idx, p)
    for k, (_, g, _) in enumerate(recs):
        if rp[k] <= 0:
            continue
        p = (int(ge[k]) + 1) / (perm + 1)
        u = g["uniPos"]
        if u not in best or p < best[u][1]:
            best[u] = (k, p)
    cands = [{"u": u, "k": k, "p": p} for u, (k, p) in best.items()]
    cands.sort(key=lambda c: c["p"])
    _bh_fdr(cands)
    by_pos = {}
    for c in cands:
        k = c["k"]
        concavity = hse_up[k] / (hse_down[k] + 1)
        by_pos[c["u"]] = {"cat": "pocket" if (cn[k] >= med_cn or concavity >= 1.3) else "exposed",
                          "p": c["p"], "q": c["q"], "score": float(rp[k])}
    return {"byPos": by_pos, "reason": None, "hasPae": bool(pae), "n": n}


def _pocket_classes(res, q=None):
    """{uniPos: 'pocket'|'exposed'} for candidates with q <= the sensitivity threshold."""
    q = POCKET_SETTINGS["q"] if q is None else float(q)
    by = (res or {}).get("byPos") if isinstance(res, dict) and "byPos" in res else (res or {})
    return {u: v["cat"] for u, v in by.items() if isinstance(v, dict) and v.get("q", 1.0) <= q}


# ----------------------------------------------------------------------------------------------
# Ligand chemistry + Tanimoto similarity (ported from api.js getLigandInfo / getLigandFingerprint)
# ----------------------------------------------------------------------------------------------
_LIG_CACHE = {}
_FP_CACHE = {}


def _num_or_none(x):
    return x if isinstance(x, (int, float)) and not isinstance(x, bool) else None


def _uri(s):
    """JavaScript encodeURIComponent."""
    return urllib.parse.quote(str(s), safe="-_.!~*'()")


def ligand_info(ccd):
    """RCSB Chemical Component Dictionary entry for a CCD code (api.js getLigandInfo): name, formula,
    molecular weight (Da), H-bond donor/acceptor counts, SMILES, InChIKey, DrugBank id."""
    key = str(ccd).upper()
    if key in _LIG_CACHE:
        return _LIG_CACHE[key]
    j = _get_json(LIGAND_CCD.format(_uri(key))) or {}
    d = j.get("rcsb_chem_comp_descriptor") or {}
    cc = j.get("chem_comp") or {}
    ci = j.get("rcsb_chem_comp_info") or {}
    db = next((x for x in (j.get("rcsb_chem_comp_related") or []) if x.get("resource_name") == "DrugBank"), None)
    info = {"id": key, "name": cc.get("name"), "formula": cc.get("formula"),
            "weight": _num_or_none(cc.get("formula_weight")),
            "hbondDonors": _num_or_none(ci.get("atom_count_hydrogen_bond_donors")),
            "hbondAcceptors": _num_or_none(ci.get("atom_count_hydrogen_bond_acceptors")),
            "smiles": d.get("SMILES_stereo") or d.get("SMILES"), "inchikey": d.get("InChIKey"),
            "drugbank": (db or {}).get("resource_accession_code")}
    _LIG_CACHE[key] = info
    return info


_CID_CACHE = {}


def pubchem_cid(inchikey):
    """InChIKey -> PubChem CID (api.js pubchemCid), or None. Cached."""
    if not inchikey:
        return None
    if inchikey in _CID_CACHE:
        return _CID_CACHE[inchikey]
    j = _get_json(PUBCHEM_CID.format(_uri(inchikey)))
    cid = ((((j or {}).get("IdentifierList") or {}).get("CID")) or [None])[0]
    _CID_CACHE[inchikey] = cid
    return cid


def ligand_links(info, cid=None):
    """External compound links for the ligand panel (api.js pubchemUrl / pubchemSimilarity*): PubChem,
    DrugBank, PubChem 2D / 3D similarity (keyed by InChIKey, upgraded to the CID query once resolved), and
    a ChEMBL search by SMILES. -> [(label, url)]"""
    info = info or {}
    ik, out = info.get("inchikey"), []
    if ik:
        out.append(("PubChem", "https://pubchem.ncbi.nlm.nih.gov/compound/%s" % _uri(ik)))
    if info.get("drugbank"):
        out.append(("DrugBank", "https://go.drugbank.com/drugs/%s" % info["drugbank"]))
    if ik:
        q = _uri("CID%s structure" % cid) if cid else _uri(ik)
        out.append(("2D similarity", "https://pubchem.ncbi.nlm.nih.gov/#query=%s&tab=similarity" % q))
        out.append(("3D similarity", "https://pubchem.ncbi.nlm.nih.gov/#query=%s&tab=similarity_3d" % q))
    if info.get("smiles"):
        out.append(("ChEMBL", "https://www.ebi.ac.uk/chembl/g/#search_results/all/query=%s" % _uri(info["smiles"])))
    return out


def ligand_fingerprint(inchikey):
    """PubChem 881-bit CACTVS 2D substructure fingerprint (decoded) for an InChIKey, or None."""
    if not inchikey:
        return None
    if inchikey in _FP_CACHE:
        return _FP_CACHE[inchikey]
    import base64
    import urllib.parse
    j = _get_json(PUBCHEM_FP.format(urllib.parse.quote(inchikey)))
    b64 = (((j or {}).get("PropertyTable") or {}).get("Properties") or [{}])[0].get("Fingerprint2D")
    fp = None
    if b64:
        raw = base64.b64decode(b64)  # first 4 bytes are a length prefix, then 881 bits
        fp = bytes(((raw[4 + (i >> 3)] if 4 + (i >> 3) < len(raw) else 0) >> (7 - (i & 7))) & 1 for i in range(881))
    _FP_CACHE[inchikey] = fp
    return fp


def tanimoto(a, b):
    """Tanimoto coefficient between two bit vectors (bytes of 0/1), or None if either is missing."""
    if not a or not b:
        return None
    inter = union = 0
    for x, y in zip(a, b):
        if x and y:
            inter += 1
        if x or y:
            union += 1
    return inter / union if union else 0.0


# ----------------------------------------------------------------------------------------------
# Fetch orchestration
# ----------------------------------------------------------------------------------------------
def fetch_annotations(uid, force=False):
    """Fetch and cache every annotation layer for a UniProt accession."""
    uid = uid.strip().upper()
    if not force and uid in _CACHE:
        _STATE["uid"] = uid          # a cache hit is still "the current accession" for omitted uids
        return _CACHE[uid]
    print("[UFV] fetching annotations for %s ..." % uid)
    # Fetch the independent endpoints in parallel — the AlphaMissense CSV is large, so doing them
    # concurrently cuts wall time from the sum to the slowest single request. All four go through
    # the on-disk cache (force=True bypasses it for a fresh copy).
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=5) as ex:
        f_features = ex.submit(_get_json, FEATURES_URL.format(uid), True, force)
        f_variation = ex.submit(_get_json, VARIATION_URL.format(uid), True, force)
        f_uniprot = ex.submit(_get_json, UNIPROT_URL.format(uid), True, force)
        f_am = ex.submit(_get_text, AM_CSV_URL.format(uid), True, force)
        f_prot = ex.submit(_get_json, PROTEOMICS_PTM_URL.format(uid), True, force)
        features = f_features.result()
        variation = f_variation.result()
        uniprot = f_uniprot.result()
        am_map = _parse_am_csv(f_am.result())
        proteomics = f_prot.result()

    sites_uni = _extract_sites_uniprot(uniprot)
    sites_api = _extract_sites(features)
    seen = {}
    for s in sites_uni + sites_api:
        seen.setdefault((s["position"], s["endPosition"]), s)
    sites = sorted(seen.values(), key=lambda s: s["position"])

    ptms = _merge_ptms(_extract_ptms(features), _extract_proteomics_ptms(proteomics))
    variants = _extract_variants(variation, am_map)
    ann = {
        "uid": uid,
        "sequence": ((uniprot or {}).get("sequence") or {}).get("value", "")
                    or (features or {}).get("sequence") or (variation or {}).get("sequence") or "",
        "ptms": ptms,
        "ptmGroups": _group_ptm_colors(ptms),          # fallback-palette colours written onto the PTMs
        "variants": variants,
        "diseaseColors": _compute_disease_colors(variants),
        "sites": sites,
        "topology": _extract_topology(features),
        "domains": _extract_domains(uniprot),
        # analysis.js aggregateAlphaMissense: the full substitution CSV, else the variants' own AM scores
        "amMean": _am_mean_by_position(am_map) or _am_mean_from_variants(variants),
        "amMap": am_map,
        "amByPos": _am_by_pos(am_map),
        "protnlm": _extract_protnlm(uniprot),
        "mutagenesis": _extract_mutagenesis(uniprot),
        "function": _extract_function(uniprot),
    }
    ann["burden"] = _compute_burden(ann["variants"])
    # Precompute per-position lookups once (residue_report reuses them on every click instead of
    # rebuilding a dict over all variants each time).
    var_by_pos = {}
    for v in ann["variants"]:
        var_by_pos.setdefault(v["position"], []).append(v)
    ann["varByPos"] = var_by_pos
    ann["ptmPosSet"] = {p["position"] for p in ann["ptms"]}
    _CACHE[uid] = ann
    _STATE["uid"] = uid
    if force:
        _invalidate_uid(uid)
        _STRUCT_CACHE.pop(uid, None)
    print("[UFV] %s: %d PTMs, %d variants, %d sites, %d domains, %d topology segments, AM=%s"
          % (uid, len(ann["ptms"]), len(ann["variants"]), len(ann["sites"]),
             len(ann["domains"]), len(ann["topology"]), "yes" if ann["amMean"] else "no"))
    return ann


def _alphafold_url(uid):
    data = _get_json(ALPHAFOLD_API.format(uid))
    if isinstance(data, list):
        for e in data:
            if e.get("pdbUrl"):
                return e["pdbUrl"]
    for ver in (6, 5, 4, 3):
        url = AF_MODEL.format(uid, ver)
        try:
            req = _urlreq().Request(url, method="HEAD", headers={"User-Agent": _USER_AGENT})
            with _urlreq().urlopen(req, timeout=30, context=_ssl_default()) as r:
                if r.status == 200:
                    return url
        except Exception:
            continue
    return None


# ----------------------------------------------------------------------------------------------
# Structure metadata (api.js decoratePdbStructures / getExperimentalStructures / get3DBeaconsModels /
# getIsoformAlphaFoldStructures). Everything here is network + pure Python (no PyMOL): safe to run in
# worker threads. Small per-entry JSON is memoised for the session on top of the disk cache.
# ----------------------------------------------------------------------------------------------
_JSON_MEMO = {}       # url -> parsed JSON (or None: PDBe answers 404 for "no mapping", never disk-cached)


def _memo_json(url):
    """_get_json + a per-session memo that also remembers definitive misses (HTTP 404: PDBe's answer for
    "no SIFTS mapping"). A network error is NOT memoised, so a transient failure is retried next time.
    The structure list, loading and the CLI ask for the same per-entry metadata many times."""
    try:
        return _JSON_MEMO[url]
    except KeyError:
        pass
    raw = _fetch_bytes(url, True, False, validate=lambda b: json.loads(b.decode("utf-8")), cache_404=b"null")
    if raw is None:
        return None
    try:
        data = json.loads(raw.decode("utf-8"))
    except Exception:
        return None
    _JSON_MEMO[url] = data
    return data


_RCSB_ENTRY_FIELDS = (
    "polymer_entities{rcsb_polymer_entity_container_identifiers{auth_asym_ids uniprot_ids}"
    "rcsb_entity_source_organism{scientific_name}"
    "polymer_entity_instances{rcsb_polymer_entity_instance_container_identifiers{auth_asym_id}"
    "rcsb_polymer_instance_info{modeled_residue_count}}}")


def _rcsb_entry_url(pid):
    # per-entry cache key (same shape as api.js RCSB_ENTRY_INSTANCES, plus entity accessions/organisms)
    return RCSB_GRAPHQL + urllib.parse.quote('{entry(entry_id:"%s"){%s}}' % (pid.upper(), _RCSB_ENTRY_FIELDS))


def _rcsb_entries(pids, chunk=100):
    """Per-entry RCSB data for the structure list, {PDBID: {"modelled": {auth_asym_id: modeled_residue_count},
    "entities": [{"chains", "uniprot", "orgs"}]}}. Entries not cached yet go out in ONE batched GraphQL
    `entries(entry_ids:[...])` query per 100 entries (the extension sends one request per chain, plus a SIFTS
    and a molecules request per entry). The batch response is disk-cached as one file (a new session reads
    one file, not one per entry) and memoised per entry; an entry the batch misses is fetched on its own."""
    pids = list(dict.fromkeys(p.upper() for p in pids))
    todo = [p for p in pids if _rcsb_entry_url(p) not in _JSON_MEMO]
    for i in range(0, len(todo), chunk):
        part = todo[i:i + chunk]
        if len(part) == 1:
            continue                                   # single entry: its own (cached) URL below
        q = "{entries(entry_ids:[%s]){rcsb_id %s}}" % (",".join('"%s"' % p for p in part), _RCSB_ENTRY_FIELDS)
        data = _get_json(RCSB_GRAPHQL + urllib.parse.quote(q))
        for e in (((data or {}).get("data") or {}).get("entries") or []):
            if e and e.get("rcsb_id"):
                _JSON_MEMO[_rcsb_entry_url(e["rcsb_id"])] = {
                    "data": {"entry": {"polymer_entities": e.get("polymer_entities") or []}}}
    out = {}
    for p in pids:
        d = _memo_json(_rcsb_entry_url(p))          # batch miss -> single-entry GET (memoised)
        modelled, ents = {}, []
        for ent in ((((d or {}).get("data") or {}).get("entry") or {}).get("polymer_entities") or []):
            ids = ent.get("rcsb_polymer_entity_container_identifiers") or {}
            ents.append({"chains": list(ids.get("auth_asym_ids") or []),
                         "uniprot": list(dict.fromkeys(ids.get("uniprot_ids") or [])),
                         "orgs": list(dict.fromkeys(o.get("scientific_name") for o in
                                                    (ent.get("rcsb_entity_source_organism") or [])
                                                    if o.get("scientific_name")))})
            for inst in ent.get("polymer_entity_instances") or []:
                ch = (inst.get("rcsb_polymer_entity_instance_container_identifiers") or {}).get("auth_asym_id")
                if ch is not None and ch not in modelled:
                    modelled[ch] = (inst.get("rcsb_polymer_instance_info") or {}).get("modeled_residue_count")
        out[p] = {"modelled": modelled, "entities": ents}
    return out


_S2A_MEMO = {}   # (pdb lower, chain) -> ({seqres: author}, {author: seqres}) or None


def _seqres_to_author(pdb_id, chain):
    """({SEQRES residue_number: author_residue_number}, inverse) for one chain, from PDBe residue_listing
    (api.js buildSeqresToAuthorMap), or None. Residues with an author insertion code are left out: PyMOL
    addresses them as e.g. '100A', which a plain author number can't select."""
    pid = pdb_id.lower()
    key = (pid, chain)
    if key in _S2A_MEMO:
        return _S2A_MEMO[key]
    data = _memo_json(PDBe_RESIDUES.format(pid))
    if data is None and PDBe_RESIDUES.format(pid) not in _JSON_MEMO:
        return None                      # network error: don't remember the miss
    found = {}
    for mol in (((data or {}).get(pid) or {}).get("molecules") or []):
        for ch in mol.get("chains") or []:
            cid = ch.get("chain_id")
            if cid is None or cid in found:      # first molecule wins (as fetchChainResidues)
                continue
            s2a = {}
            for r in ch.get("residues") or []:
                a, s = r.get("author_residue_number"), r.get("residue_number")
                if a is None or s is None or str(r.get("author_insertion_code") or "").strip():
                    continue
                s2a[int(s)] = int(a)
            found[cid] = (s2a, {a: s for s, a in s2a.items()}) if s2a else None
    for cid, v in found.items():
        _S2A_MEMO[(pid, cid)] = v
    _S2A_MEMO.setdefault(key, None)
    return _S2A_MEMO[key]


def _seg_u2r(seg, pos):
    """UniProt position -> author resi through one segment (the caller checked u_start <= pos <= u_end).
    Segments carrying a SEQRES->author map resolve per residue (gaps / jumps in author numbering)."""
    s2a = seg.get("s2a")
    if s2a is not None:
        return s2a.get(seg["seqres_start"] + (pos - seg["u_start"]))
    return seg["pdb_start"] + (pos - seg["u_start"])


def _seg_r2u(seg, resi):
    """Author resi -> UniProt position through one segment, or None if outside it."""
    a2s = seg.get("a2s")
    if a2s is not None:
        s = a2s.get(resi)
        if s is None:
            return None
        u = seg["u_start"] + (s - seg["seqres_start"])
        return u if seg["u_start"] <= u <= seg["u_end"] else None
    lo, hi = sorted((seg["pdb_start"], seg["pdb_end"]))
    if lo <= resi <= hi:
        return seg["u_start"] + (resi - seg["pdb_start"])
    return None


def _attach_seqres(seg, maps, seqres_start):
    seg["seqres_start"] = int(seqres_start)
    seg["s2a"], seg["a2s"] = maps
    return seg


def _sifts_uniprot(pdb_id):
    """PDBe SIFTS {accession: {mappings: [...]}} for an entry (memoised; {} when PDBe has none)."""
    pid = pdb_id.lower()
    return (((_memo_json(PDBe_SIFTS.format(pid)) or {}).get(pid)) or {}).get("UniProt") or {}


def _sifts_acc_entry(uni, uid):
    """Our accession's SIFTS entry: exact, base accession (isoform suffix dropped) or case-insensitive."""
    entry = uni.get(uid) or uni.get(uid.split("-")[0])
    if not entry:
        for k in uni:
            if k.upper() == uid.upper():
                return uni[k]
    return entry


def _norm_best(it):
    """(unp_start, unp_end, seqres_start, seqres_end) of a best_structures row (normalizeBestStructure)."""
    def num(*keys):
        for k in keys:
            v = it.get(k)
            if v not in (None, "", 0):
                try:
                    return int(v)
                except (TypeError, ValueError):
                    pass
        return 0
    us = num("uniprot_start", "unp_start", "start")
    ue = num("uniprot_end", "unp_end", "end")
    return us, ue, num("pdb_start", "start") or us, num("pdb_end", "end") or ue


def _best_ranges(uid, pdb_id):
    """{chain: [{u_start, u_end, seqres_start, seqres_end}]} for one entry from PDBe best_structures
    (SEQRES numbering) — the mapping the extension starts from before SIFTS refines it."""
    acc = uid.upper()
    best = _memo_json(PDBe_BEST.format(acc)) or {}
    out = {}
    for it in best.get(acc) or best.get(acc.split("-")[0]) or []:
        if (it.get("pdb_id") or "").lower() != pdb_id.lower():
            continue
        us, ue, ps, pe = _norm_best(it)
        if not (us and ue):
            continue
        r = {"u_start": us, "u_end": ue, "seqres_start": ps, "seqres_end": pe}
        rl = out.setdefault(it.get("chain_id") or "A", [])
        if r not in rl:
            rl.append(r)
    return out


def _sifts_segments(pdb_id, uid, best=None, listing=True):
    """{chain: [{u_start, u_end, pdb_start, pdb_end, ...}]} mapping `uid` onto `pdb_id`'s AUTHOR residue
    numbers (PyMOL's resi). Port of api.js decoratePdbStructures (L489-537):
      * SIFTS segment with author numbers -> linear segment;
      * SIFTS gives only the SEQRES residue_number -> per-residue SEQRES->author map from PDBe
        residue_listing (the segment carries s2a/a2s); without a listing it is flagged "label" and
        re-anchored from the model's own sequence when loaded;
      * SIFTS has no entry for this accession (chimeric / engineered entries; also when PDBe answers 404,
        where the extension keeps raw SEQRES numbers) -> the best_structures SEQRES ranges are translated
        through residue_listing and flagged "unverified" (checked against the sequence on load).
    `best` = _best_ranges() result (fetched when None and needed). listing=False skips residue_listing
    (a large download per entry): the structure LIST only needs the mapping status; the per-residue map is
    built when the entry is loaded."""
    pid = pdb_id.lower()
    uid = uid.upper()
    entry = _sifts_acc_entry(_sifts_uniprot(pid), uid)
    segs = {}
    if entry and entry.get("mappings"):
        for m in entry["mappings"]:
            ch = m.get("chain_id")
            start, end = m.get("start") or {}, m.get("end") or {}
            seg = {"u_start": m.get("unp_start"), "u_end": m.get("unp_end"),
                   "pdb_start": start.get("author_residue_number"), "pdb_end": end.get("author_residue_number")}
            if seg["pdb_start"] is None or seg["pdb_end"] is None:
                maps = _seqres_to_author(pid, ch) if (listing and start.get("residue_number") is not None) else None
                if maps:
                    s2a = maps[0]
                    if seg["pdb_start"] is None:
                        seg["pdb_start"] = s2a.get(start.get("residue_number"), m.get("unp_start"))
                    if seg["pdb_end"] is None:
                        seg["pdb_end"] = s2a.get(end.get("residue_number"), m.get("unp_end"))
                    _attach_seqres(seg, maps, start["residue_number"])
                else:
                    seg["label"] = True
                    if seg["pdb_start"] is None:
                        seg["pdb_start"] = start.get("residue_number", m.get("unp_start"))
                    if seg["pdb_end"] is None:
                        seg["pdb_end"] = end.get("residue_number", m.get("unp_end"))
            segs.setdefault(ch, []).append(seg)
        for sl in segs.values():
            sl.sort(key=lambda s: s["u_start"])
        return segs
    if best is None:
        best = _best_ranges(uid, pid)
    for ch, rl in (best or {}).items():
        maps = _seqres_to_author(pid, ch) if listing else ({}, {})
        if not maps:
            continue
        s2a = maps[0]
        for r in rl:
            seg = {"u_start": r["u_start"], "u_end": r["u_end"],
                   "pdb_start": s2a.get(r["seqres_start"], r["seqres_start"]),
                   "pdb_end": s2a.get(r["seqres_end"], r["seqres_end"]),
                   "unverified": True, "via": "residue_listing"}
            segs.setdefault(ch, []).append(_attach_seqres(seg, maps, r["seqres_start"]))
    return segs


def _segments_status(segs):
    """Mapping status line (api.js mappingStatus); segment count = the most segmented chain."""
    if not segs:
        return "Residue mapping unavailable."
    allsegs = [s for v in segs.values() for s in v]
    n = max(len(v) for v in segs.values())
    if all(s.get("via") == "sequence" for s in allsegs):
        return "Mapped by sequence (PDBe numbering re-anchored to the model's author numbering)."
    if any(s.get("via") == "residue_listing" for s in allsegs):
        return "Mapped via residue_listing (no SIFTS entry for this protein)."
    if any(s.get("s2a") is not None for s in allsegs):
        return "Mapped with SIFTS + residue_listing (%d segment%s)" % (n, "" if n == 1 else "s")
    return "Mapped with SIFTS (%d segment%s)" % (n, "" if n == 1 else "s")


def _pdb_entry_meta(pdb_id, uid, chains):
    """Chimeric flag, organisms, chain names and partner chains for OUR chains of an entry — api.js
    L477-486 + extractPartnerMappings L283. chimeric ONLY when one of OUR chains is the chimera (its entity
    spans >1 source organism, or SIFTS maps >1 accession to it); a complex that merely contains
    other-species chains (e.g. a mouse Fab on a human target) is not flagged."""
    pid = pdb_id.lower()
    uni = _sifts_uniprot(pid)
    mols = (_memo_json(PDBe_MOLECULES.format(pid)) or {}).get(pid) or []
    chains = list(chains or [])
    orgs, chim, names = [], [], {}
    for ch in chains:
        mol = next((m for m in mols if ch in (m.get("in_chains") or [])), None)
        o = list(dict.fromkeys(s.get("organism_scientific_name") for s in ((mol or {}).get("source") or [])
                               if s.get("organism_scientific_name")))
        n_acc = sum(1 for e in uni.values() if any(m.get("chain_id") == ch for m in (e.get("mappings") or [])))
        if len(o) > 1 or n_acc > 1:
            chim.append(ch)
            orgs.extend(x for x in o if x not in orgs)
        if mol:
            names[ch] = ((mol.get("molecule_name") or [""])[0]) or mol.get("molecule_type") or ""
    if not chim:   # organisms of our chains for display
        for ch in chains:
            mol = next((m for m in mols if ch in (m.get("in_chains") or [])), None)
            for s in ((mol or {}).get("source") or []):
                o = s.get("organism_scientific_name")
                if o and o not in orgs:
                    orgs.append(o)
    our_id = uid.upper()
    base = our_id.split("-")[0]
    ours = set(chains)
    partners = []
    for acc, e in uni.items():
        if acc.upper() in (our_id, base):
            continue
        for m in e.get("mappings") or []:
            if m.get("chain_id") in ours:
                continue
            partners.append({"accession": acc, "chainId": m.get("chain_id"),
                             "uniprotStart": m.get("unp_start"), "uniprotEnd": m.get("unp_end"),
                             "startAuthor": (m.get("start") or {}).get("author_residue_number"),
                             "startSeqres": (m.get("start") or {}).get("residue_number"),
                             "name": (e.get("name") or e.get("identifier") or "")})
    # RCSB entity accessions / organisms (the batched query the structure list already made; one memoised
    # request otherwise): catches chimeric constructs and partner chains PDBe has no SIFTS for (6HUO: PDBe
    # 404 + molecules lists only Homo sapiens, RCSB: α1 entity = P14867 + bovine P08219).
    try:
        ents = (_rcsb_entries([pdb_id]).get(pdb_id.upper()) or {}).get("entities") or []
    except Exception:
        ents = []
    have = {(p["accession"].upper(), p["chainId"]) for p in partners}
    for e in ents:
        mine = sorted(ours & set(e["chains"]))
        if mine:
            if len(e["orgs"]) > 1 or len(e["uniprot"]) > 1:
                chim += [c for c in mine if c not in chim]
                orgs += [o for o in e["orgs"] if o not in orgs]
            continue
        for acc in e["uniprot"]:
            if acc.upper() in (our_id, base):
                continue
            for c in e["chains"]:
                if (acc.upper(), c) not in have:
                    have.add((acc.upper(), c))
                    partners.append({"accession": acc, "chainId": c, "uniprotStart": None, "uniprotEnd": None,
                                     "startAuthor": None, "startSeqres": None, "name": ""})
    return {"chimeric": bool(chim), "chimericChains": chim, "organisms": orgs, "chainNames": names,
            "partners": partners}


def _partner_uni_to_author(pdb_id, m):
    """UniProt position -> author resi resolver for one partner-chain mapping (api.js partnerUniToAuthor)."""
    if m.get("startAuthor") is not None:
        return lambda pos: m["startAuthor"] + (pos - m["uniprotStart"])
    maps = _seqres_to_author(pdb_id, m["chainId"]) if m.get("startSeqres") is not None else None
    if maps:
        s2a = maps[0]
        return lambda pos: s2a.get(m["startSeqres"] + (pos - m["uniprotStart"]))
    return lambda pos: None


def _pdb_chain_hints(pdb_id, uid, protein_name=None):
    """Author chain ids of `pdb_id` that carry `uid`, for entries where SIFTS has no residue-level
    mapping: PDBe best_structures (per-accession chain list), else PDBe molecules (polypeptide
    entities whose name matches the UniProt protein name). [] if neither identifies a chain."""
    pid = pdb_id.lower()
    acc = uid.split("-")[0].upper()
    chains = []
    best = _memo_json(PDBe_BEST.format(acc)) or {}
    for it in best.get(acc) or []:
        ch = it.get("chain_id")
        if (it.get("pdb_id") or "").lower() == pid and ch and ch not in chains:
            chains.append(ch)
    if chains:
        return chains
    want = (protein_name or "").strip().lower()
    if not want:
        return []
    mol = _memo_json(PDBe_MOLECULES.format(pid)) or {}
    for ent in mol.get(pid) or []:
        if not str(ent.get("molecule_type", "")).lower().startswith("polypeptide"):
            continue
        names = [str(n).lower() for n in (ent.get("molecule_name") or [])]
        names.append(str(ent.get("synonym") or "").lower())
        if any(want in n for n in names if n):
            chains += [c for c in (ent.get("in_chains") or []) if c not in chains]
    return chains


def _js_round1(frac):
    """Math.round(frac*1000)/10 (JS rounds .5 up; Python's round() is banker's)."""
    return _math.floor(frac * 1000 + 0.5) / 10.0


def _js_num(x):
    """Number formatted like a JS template literal: 100 -> '100', 93.9 -> '93.9', 3.26 -> '3.26'."""
    try:
        return "%g" % float(x)
    except (TypeError, ValueError):
        return str(x)


def _short_method(m):
    """modal.js shortMethod."""
    if not m:
        return ""
    lc = str(m).lower()
    if "electron" in lc or "cryo" in lc:
        return "EM"
    if "x-ray" in lc or "diffraction" in lc or "crystallography" in lc:
        return "X-ray"
    if "nmr" in lc:
        return "NMR"
    if "neutron" in lc:
        return "Neutron"
    return "".join(w[:1].upper() for w in str(m).split(" ")) or "Exp"


def _structure_label(s):
    """Structure-selector label (modal.js stLabel)."""
    src = s.get("source")
    if src == "AlphaFold":
        return ("AlphaFold %s" % s["isoform"]) if s.get("isoform") else "AlphaFold"
    if src == "Computed":
        prov = s.get("provider") or "Computed model"
        return "%s (model, %s%%)" % (prov, _js_num(s["cov"])) if s.get("cov") else "%s (model)" % prov
    if src == "Custom":
        return "Custom: %s" % re.sub(r"\.(pdb|ent)$", "", s.get("name") or "Custom PDB", flags=re.I)
    ids = s.get("chainIds") or []
    chains = ",".join(ids) if len(ids) > 1 else (s.get("chainId") or "?")
    parts = []
    if s.get("method"):
        parts.append(_short_method(s["method"]))
    if s.get("resolution"):
        parts.append("%sÅ" % _js_num(s["resolution"]))
    if s.get("cov"):
        parts.append("%s%%" % _js_num(s["cov"]))
    return "%s-%s%s%s" % (s.get("pdbId"), chains, (" (%s)" % ", ".join(parts)) if parts else "",
                          " ⚛" if s.get("chimeric") else "")


def _chimeric_tip(s):
    if not s or not s.get("chimeric"):
        return ""
    orgs = s.get("organisms") or []
    return ("⚛ Chimeric chain — %s" % " / ".join(orgs)) if len(orgs) > 1 else \
        "⚛ Chimeric chain — our protein is fused with another sequence"


def _isoform_vsps(uniprot, vsp_ids):
    """VAR_SEQ edits defining one isoform as [{start, end, newLen}] in canonical coordinates
    (api.js isoformVsps). An empty replacement is a deletion (newLen 0)."""
    by_id = {f.get("featureId"): f for f in ((uniprot or {}).get("features") or [])
             if f.get("type") == "Alternative sequence" and f.get("featureId")}
    out = []
    for vid in vsp_ids or []:
        f = by_id.get(vid)
        loc = (f or {}).get("location") or {}
        s, e = (loc.get("start") or {}).get("value"), (loc.get("end") or {}).get("value")
        if s is None or e is None:
            continue
        alt = ((f.get("alternativeSequence") or {}).get("alternativeSequences"))
        out.append({"start": int(s), "end": int(e),
                    "newLen": len(str(alt[0])) if isinstance(alt, list) and alt else 0})
    return out


def _isoform_ranges(canon_len, vsps):
    """Canonical position <-> isoform residue segments from VSP edits (api.js buildIsoformMappedRanges):
    unedited stretches map 1:1 with a running offset; canonical positions inside an edited region have no
    isoform counterpart (annotations there are not drawn)."""
    ranges, canon, iso = [], 1, 1
    for v in sorted(vsps, key=lambda v: v["start"]):
        if v["start"] < canon:          # overlapping / duplicate edit
            continue
        if v["start"] > canon:
            n = v["start"] - canon
            ranges.append({"u_start": canon, "u_end": v["start"] - 1, "pdb_start": iso, "pdb_end": iso + n - 1})
            iso += n
        iso += v["newLen"]
        canon = v["end"] + 1
    if canon <= canon_len:
        n = canon_len - canon + 1
        ranges.append({"u_start": canon, "u_end": canon_len, "pdb_start": iso, "pdb_end": iso + n - 1})
    return ranges


def _isoform_structures(uid, seqlen=0):
    """AlphaFold models of the NON-canonical isoforms (e.g. AF-P35498-2), resolved through each isoform's
    3D-Beacons summary (api.js getIsoformAlphaFoldStructures). Numbering follows the isoform sequence; the
    canonical<->isoform mapping comes from the isoform's VSP edits (numbering "ranges"), so canonical
    annotations land on the right residues after splice edits."""
    base = uid.split("-")[0]
    # same URL as fetch_annotations: one download (per-URL lock) and disk-cached
    u = _memo_json(UNIPROT_URL.format(base)) or {}
    seqlen = seqlen or len(((u.get("sequence") or {}).get("value")) or "")
    iso_map = {}
    for c in (u.get("comments") or []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in (c.get("isoforms") or []):
                for iid in (iso.get("isoformIds") or []):
                    iso_map[iid] = iso.get("sequenceIds") or []
    targets = [i for i in iso_map if i and i not in (base, base + "-1")]
    if not targets:
        return []
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=min(DECORATE_CONCURRENCY, len(targets))) as ex:
        summaries = list(ex.map(lambda i: _memo_json(BEACONS_SUMMARY.format(i)) or {}, targets))
    out = []
    for iid, d in zip(targets, summaries):
        af = next((sm for e in (d.get("structures") or []) for sm in [(e or {}).get("summary") or e or {}]
                   if sm.get("model_url") and re.search("alphafold", sm.get("provider") or "", re.I)), None)
        if not af:
            continue
        fmt = str(af.get("model_format") or "").upper()
        vsps = _isoform_vsps(u, iso_map.get(iid))
        if seqlen:
            ranges = _isoform_ranges(seqlen, vsps) if vsps else \
                [{"u_start": 1, "u_end": seqlen, "pdb_start": 1, "pdb_end": seqlen}]
        else:
            ranges = []
        status = ("AlphaFold isoform model — canonical annotations mapped across %d segment%s; residues "
                  "inside spliced regions are not shown." % (len(ranges), "" if len(ranges) == 1 else "s")) \
            if vsps else "AlphaFold isoform model — no alternative-sequence data found; positions assumed " \
                         "identical to canonical."
        s = {"key": "AF-%s" % iid, "source": "AlphaFold", "isoform": iid, "url": af["model_url"],
             "fmt": "pdb" if fmt == "PDB" else "mmcif", "pdbId": None, "chainId": None,
             "numbering": "ranges" if ranges else "identity", "ranges": ranges, "cov": None,
             "method": "Predicted model (isoform)", "rangeText": "isoform %s" % iid, "mappingStatus": status}
        s["label"] = _structure_label(s)
        out.append(s)
    return out


def _beacons_models(uid, seqlen=0):
    """Computed / predicted models from 3D-Beacons (SWISS-MODEL, ModelArchive, PED, AlphaFill, non-F1
    AlphaFold...) — api.js get3DBeaconsModels: experimental entries and the canonical AF-<id>-F1 skipped,
    PDB/mmCIF only, host allow-list, de-duplicated by model_identifier, capped at 25. Numbering is assumed
    to follow UniProt positions."""
    url = BEACONS_SUMMARY.format(uid)
    data = _memo_json(url)
    if data is None and url not in _JSON_MEMO:      # network failure (a 404 "no models" is memoised)
        _LIST_NOTICES[uid.upper()] = ("3D-Beacons (%s) unreachable: computed models (SWISS-MODEL, "
                                      "ModelArchive, AlphaFill...) are not listed; fetch again later."
                                      % urllib.parse.urlparse(url).hostname)
    data = data or {}
    seen, out = set(), []
    for e in data.get("structures") or []:
        sm = (e or {}).get("summary") or e or {}
        murl = sm.get("model_url")
        if not murl:
            continue
        prov = sm.get("provider") or "Computed model"
        if "EXPERIMENTAL" in str(sm.get("model_category") or "").upper():
            continue
        if re.search("alphafold", prov, re.I) and \
                str(sm.get("model_identifier") or "").upper() == "AF-%s-F1" % uid.upper():
            continue
        fmt = str(sm.get("model_format") or "").upper()
        if fmt not in ("PDB", "MMCIF"):
            continue
        try:
            host = urllib.parse.urlparse(murl).hostname
        except Exception:
            continue
        if host not in BEACON_ALLOWED_HOSTS:
            continue
        mid = sm.get("model_identifier") or murl
        if mid in seen:
            continue
        seen.add(mid)
        try:
            us = int(sm.get("uniprot_start") or 1)
            ue = int(sm.get("uniprot_end") or seqlen or 0)
        except (TypeError, ValueError):
            us, ue = 1, seqlen
        if sm.get("coverage") is not None:
            cov = _js_round1(float(sm["coverage"]))
        else:
            cov = _js_round1((ue - us + 1) / float(seqlen)) if (seqlen and ue) else None
        s = {"key": "BCN-%s" % mid, "source": "Computed", "provider": prov, "url": murl,
             "fmt": "mmcif" if fmt == "MMCIF" else "pdb", "pdbId": None, "chainId": None,
             "numbering": "identity", "cov": cov, "method": "%s model" % prov,
             "resolution": sm.get("resolution"), "modelCategory": sm.get("model_category") or "",
             "rangeText": ("%d-%d" % (us, ue)) if (us and ue) else (("1-%d" % seqlen) if seqlen else ""),
             "modelPageUrl": sm.get("model_page_url") or "", "confidence": sm.get("confidence_avg_local_score"),
             "mappingStatus": "%s computed model — residue numbering assumed to follow UniProt positions." % prov}
        s["label"] = _structure_label(s)
        out.append(s)
    return out[:BEACONS_CAP]


def _experimental_structures(uid, seqlen=0):
    """Experimental PDB entries (api.js getExperimentalStructures + decoratePdbStructures): best_structures
    rows de-duplicated per (entry, chain), decorated with RCSB modelled-residue coverage and the chimeric
    flag, then all chains of one entry that carry this protein MERGED into one row (key `6HUO_A,D`,
    highest-coverage chain primary); the entry's other proteins are attached as partners.
    Listing cost is constant: best_structures + ONE batched RCSB GraphQL query (coverage, entity
    accessions and source organisms for every entry). The per-entry SIFTS / residue_listing / molecules
    requests the extension makes here (3-4 per entry: ~45 s cold for P05067's 270 entries) are deferred to
    load time (_prepare_structure), where the author-numbering map, the SIFTS-based chimeric check and the
    partner residue ranges are computed for the one entry being loaded."""
    best = _memo_json(PDBe_BEST.format(uid)) or {}
    items = best.get(uid) or best.get(uid.upper()) or []
    rows = {}
    for it in items:
        pid = (it.get("pdb_id") or it.get("pdbId") or "").upper()
        if not pid:
            continue
        ch = it.get("chain_id") or it.get("chain") or "A"
        us, ue, ps, pe = _norm_best(it)
        if not (us and ue):
            continue
        cov = _js_round1((ue - us + 1) / float(seqlen)) if seqlen else 0
        prev = rows.get((pid, ch))
        if prev is None or cov > prev["cov"]:
            rows[(pid, ch)] = {"pdbId": pid, "chainId": ch, "cov": cov, "u": (us, ue),
                               "method": it.get("experimental_method") or it.get("method") or "Experimental",
                               "resolution": it.get("resolution")}
    by_pdb = {}
    for r in rows.values():
        by_pdb.setdefault(r["pdbId"], []).append(r)
    pids = list(by_pdb)
    if len(pids) > STRUCTURE_LIST_CAP:
        print("[UFV] %s maps to %d PDB entries; listing the first %d." % (uid, len(pids), STRUCTURE_LIST_CAP))
        pids = pids[:STRUCTURE_LIST_CAP]
    if not pids:
        return []
    rcsb = _rcsb_entries(pids)
    our_id, base = uid.upper(), uid.upper().split("-")[0]
    out = []
    for pid in pids:
        info = rcsb.get(pid) or {}
        ents = info.get("entities") or []
        merged = None
        for r in by_pdb[pid]:
            ch = r["chainId"]
            mc = (info.get("modelled") or {}).get(ch)
            cov = _js_round1(mc / float(seqlen)) if (seqlen and mc) else r["cov"]
            rng = "%d-%d" % r["u"]
            if merged is None:
                merged = {"key": pid, "source": "PDB", "url": PDBe_PDB.format(pid.lower()),
                          "cifUrl": PDBe_CIF.format(pid.lower()), "fmt": "pdb", "pdbId": pid,
                          "chainId": ch, "chainIds": [ch], "numbering": "sifts", "cov": cov,
                          "method": r["method"], "resolution": r["resolution"], "rangeText": rng,
                          "mappingStatus": "Mapped with PDBe/SIFTS best-structure data; author numbering "
                                           "is resolved through SIFTS / residue_listing on load."}
            elif ch not in merged["chainIds"]:
                merged["chainIds"].append(ch)
                if (cov or 0) > (merged["cov"] or 0):        # promote the higher-coverage chain
                    merged.update(chainId=ch, cov=cov, rangeText=rng)
        merged["chainIds"].sort()
        ours = set(merged["chainIds"])
        # ⚛ only when OUR chain is the chimera: its entity spans >1 organism or >1 UniProt entry
        chim, orgs = [], []
        for e in ents:
            mine = sorted(ours & set(e["chains"]))
            if not mine:
                continue
            if len(e["orgs"]) > 1 or len(e["uniprot"]) > 1:
                chim += mine
            orgs += [o for o in e["orgs"] if o not in orgs]
        partners = []
        for e in ents:
            for acc in e["uniprot"]:
                if acc.upper() in (our_id, base):
                    continue
                for c in e["chains"]:
                    if c not in ours:
                        partners.append({"accession": acc, "chainId": c})
        merged.update(chimeric=bool(chim), chimericChains=chim, organisms=orgs, partners=partners)
        merged["key"] = "%s_%s" % (pid, ",".join(merged["chainIds"]))
        merged["label"] = _structure_label(merged)
        out.append(merged)
    return out


_STRUCT_RANK = {"AlphaFold": 0, "PDB": 1, "Computed": 2}


def _structure_rank(s):
    if s.get("source") == "AlphaFold":
        return 0.5 if s.get("isoform") else 0
    return _STRUCT_RANK.get(s.get("source"), 1)


_STRUCT_CACHE = {}   # uid -> list_structures() result (this session)
_LIST_NOTICES = {}   # uid -> one-line notice about the last listing (a source that could not be reached)


def structure_notice(uid):
    """One-line notice about the last list_structures(uid) ("" when every source answered), e.g. that
    3D-Beacons was unreachable so no computed models are listed. Both plugins show it."""
    return _LIST_NOTICES.get((uid or "").upper(), "")


def list_structures(uid, seqlen=0, force=False):
    """Available structures for an accession (for the structure selector) — api.js getStructures: the
    canonical AlphaFold model, isoform AlphaFold models (VSP-mapped), experimental PDB entries (chains of
    one entry merged; SIFTS / residue_listing numbering, chimeric flag, RCSB modelled coverage) and
    3D-Beacons computed models. Order: AlphaFold, isoforms, PDB, computed; coverage desc within each.
    Cached per accession for the session (force=True re-lists and re-reads the metadata)."""
    uid = uid.upper()
    if not force and uid in _STRUCT_CACHE:
        return list(_STRUCT_CACHE[uid])
    _LIST_NOTICES.pop(uid, None)
    if force:
        _JSON_MEMO.clear()
        _S2A_MEMO.clear()
    if not seqlen:
        seqlen = len((_CACHE.get(uid) or {}).get("sequence") or "") or \
            len((((_memo_json(UNIPROT_URL.format(uid)) or {}).get("sequence")) or {}).get("value") or "")
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as ex:
        f_af = ex.submit(_alphafold_url, uid)
        f_iso = ex.submit(_isoform_structures, uid, seqlen)
        f_pdb = ex.submit(_experimental_structures, uid, seqlen)
        f_bcn = ex.submit(_beacons_models, uid, seqlen)
        af = f_af.result()
        isoforms = f_iso.result() or []
        pdbs = f_pdb.result() or []
        beacons = f_bcn.result() or []
    out = []
    if af:
        out.append({"key": "AF-%s" % uid, "label": "AlphaFold", "source": "AlphaFold", "url": af,
                    "fmt": "mmcif" if af.lower().endswith(".cif") else "pdb", "pdbId": None, "chainId": None,
                    "numbering": "identity", "cov": 100.0, "method": "Predicted model",
                    "rangeText": ("1-%d" % seqlen) if seqlen else "full sequence",
                    "mappingStatus": "AlphaFold numbering follows UniProt sequence positions."})
    out.extend(isoforms)
    out.extend(pdbs)
    out.extend(beacons)
    out.sort(key=lambda s: (_structure_rank(s), -(s.get("cov") or 0)))   # stable, like Array.sort
    if uid in _LIST_NOTICES:
        print("[UFV] " + _LIST_NOTICES[uid])          # visible in the console (stderr under --serve)
        _STRUCT_CACHE.pop(uid, None)
    else:
        _STRUCT_CACHE[uid] = out                     # an incomplete list is not cached: retried next time
    return list(out)


# ----------------------------------------------------------------------------------------------
# Residue-number mapping (UniProt position -> object chain+resi)
# ----------------------------------------------------------------------------------------------
def _set_identity_map(obj):
    UFV_MAPS[obj] = {"mode": "identity", "chains": {}}


def _seq_offset(res, seq, min_identity=0.9, kmer=5):
    """Best constant offset (author resi - UniProt position) aligning a chain's residues `res`
    [(resv, aa)] to the UniProt sequence: candidate offsets come from shared 5-mers, then each is
    scored by identity. Returns (offset, identity, matches) if identity >= min_identity with enough
    matched residues, else None. This is how a chain is proven to BE this accession — and how its
    author numbering is related to UniProt — when SIFTS has no residue-level mapping."""
    if not res or not seq:
        return None
    from collections import Counter
    starts = {}
    for i in range(len(seq) - kmer + 1):
        starts.setdefault(seq[i:i + kmer], []).append(i + 1)
    by = dict(res)
    votes = Counter()
    for r, _aa in res:
        word = "".join(by.get(r + t, "?") for t in range(kmer))
        if "?" in word or "X" in word:
            continue
        for p in starts.get(word, ()):
            votes[r - p] += 1
    best = None
    for off in dict.fromkeys([0] + [o for o, _ in votes.most_common(3)]):
        tot = match = 0
        for r, aa in res:
            u = r - off
            if 1 <= u <= len(seq):
                tot += 1
                match += seq[u - 1] == aa
        if tot and (best is None or match > best[2]):
            best = (off, match / float(tot), match)
    need = min(20, max(5, len(res) // 2))
    if best and best[1] >= min_identity and best[2] >= need:
        return best
    return None


def _offset_segment(res, seq, off):
    us = [r - off for r, _aa in res if 1 <= r - off <= len(seq)]
    lo, hi = min(us), max(us)
    return {"u_start": lo, "u_end": hi, "pdb_start": lo + off, "pdb_end": hi + off}


def _segs_identity(res, seq, sl):
    """(mapped, identical) residue counts of a chain's residues [(resv, aa)] under segment list `sl`."""
    tot = match = 0
    n = len(seq)
    for r, aa in res:
        for s in sl:
            u = _seg_r2u(s, r)
            if u is not None:
                if 1 <= u <= n:
                    tot += 1
                    match += seq[u - 1] == aa
                break
    return tot, match


def _ensure_manual_map(obj):
    mp = UFV_MAPS.get(obj)
    if not mp or mp.get("mode") != "manual":
        mp = {"mode": "manual", "chains": {}}
        UFV_MAPS[obj] = mp
    return mp


def _uni_to_resi(obj, chain, pos):
    """Map a UniProt position to the object's author residue number on `chain` (or None)."""
    mp = UFV_MAPS.get(obj)
    if not mp:
        return pos  # default identity
    mode = mp["mode"]
    if mode == "identity":
        return None if (mp.get("exclude") and chain in mp["exclude"]) else pos
    if mode in ("sifts", "custom"):
        for seg in mp["segments"].get(chain, ()):
            if seg["u_start"] <= pos <= seg["u_end"]:
                return _seg_u2r(seg, pos)
        return None
    if mode == "ranges":            # isoform model: canonical -> isoform residue (same on every chain)
        for seg in mp["ranges"]:
            if seg["u_start"] <= pos <= seg["u_end"]:
                return _seg_u2r(seg, pos)
        return None
    if mode == "manual":
        c = mp["chains"].get(chain)
        if not c:
            return None
        resi = c["resi_start"] + (pos - c["u_start"])
        if resi < c["resi_start"]:
            return None
        if c["resi_end"] is not None and resi > c["resi_end"]:
            return None
        return resi
    return pos


def _resi_int(resi, strict=True):
    """Author residue (int or PyMOL resi string) -> int. strict: residues with an insertion code
    ('100A') return None — they are insertions relative to the reference numbering and have no
    UniProt position of their own; non-strict keeps the leading number (for sorting)."""
    if isinstance(resi, int):
        return resi
    m = re.match(r"^\s*(-?\d+)\s*$" if strict else r"^\s*(-?\d+)", str(resi))
    return int(m.group(1)) if m else None


def _resi_to_uni(obj, chain, resi):
    """Inverse of _uni_to_resi: object author residue -> UniProt position (or None). Residues with an
    insertion code are deliberately unmapped (see _resi_int)."""
    resi = _resi_int(resi)
    if resi is None:
        return None
    mp = UFV_MAPS.get(obj)
    if not mp or mp["mode"] == "identity":
        return None if (mp and mp.get("exclude") and chain in mp["exclude"]) else resi
    mode = mp["mode"]
    if mode in ("sifts", "custom", "ranges"):
        for seg in (mp["ranges"] if mode == "ranges" else mp["segments"].get(chain, ())):
            u = _seg_r2u(seg, resi)
            if u is not None:
                return u
        return None
    if mode == "manual":
        c = mp["chains"].get(chain)
        if not c:
            return None
        if resi < c["resi_start"]:
            return None
        if c["resi_end"] is not None and resi > c["resi_end"]:
            return None
        return c["u_start"] + (resi - c["resi_start"])
    return resi


# ----------------------------------------------------------------------------------------------
# Custom-structure mapping (modal.js Custom PDB, L3311-3600): per-chain auto | offset | adopt | ignore,
# per-chain accession. Works on any loaded object (your own PDB, a trajectory, a listed structure).
# ----------------------------------------------------------------------------------------------


def _same_acc(a, b):
    a, b = (a or "").upper(), (b or "").upper()
    return a == b or re.sub(r"-1$", "", a) == re.sub(r"-1$", "", b)


def _map_chain_segments(mp, chain, seqlen=0):
    """The segments a mapping descriptor (any mode) gives `chain`, as a fresh list, or None."""
    if not mp:
        return None
    mode = mp.get("mode")
    if mode in ("sifts", "custom"):
        sl = (mp.get("segments") or {}).get(chain)
        return [dict(s) for s in sl] if sl else None
    if mode == "ranges":
        return [dict(s) for s in mp.get("ranges") or []] or None
    if mode == "manual":
        c = (mp.get("chains") or {}).get(chain)
        if not c:
            return None
        u_end = c["u_start"] + (c["resi_end"] - c["resi_start"]) if c["resi_end"] is not None \
            else max(seqlen, c["u_start"])
        return [{"u_start": c["u_start"], "u_end": u_end, "pdb_start": c["resi_start"],
                 "pdb_end": c["resi_start"] + (u_end - c["u_start"])}]
    if mode == "identity" and seqlen:
        return [{"u_start": 1, "u_end": seqlen, "pdb_start": 1, "pdb_end": seqlen}]
    return None


def _map_primary_chain(mp):
    if not mp:
        return None
    if mp.get("chainId"):
        return mp["chainId"]
    segs = mp.get("segments") or mp.get("chains") or {}
    return sorted(segs)[0] if segs else None


def _seg_text(sl):
    return "; ".join("UniProt %d-%d <-> resi %s-%s" % (s["u_start"], s["u_end"], s["pdb_start"], s["pdb_end"])
                     for s in sl)


# ----------------------------------------------------------------------------------------------
# Projection helpers
# ----------------------------------------------------------------------------------------------


def _resolve_uid(uid, obj=None):
    """Accession to use: the explicit one, else the accession `obj` was loaded/mapped for, else the
    most recently fetched one."""
    if uid:
        return uid
    if obj:
        u = (_STATE.get("obj_uid") or {}).get(obj)
        if u:
            return u
    return _STATE.get("uid")


# Each point (sphere) layer is drawn on its OWN object (a copy of just the annotated Cα atoms, no
# cartoon) named ufv_<obj>_<tag>. This keeps the main structure's cartoon out of every sphere
# rebuild — large variant sets no longer crawl — and keeps the object panel to a few clean names.
_SPHERE_STATE = {}    # obj -> { tag -> {color_hex: [uniprot_positions]} }  (active sphere layers)


def _cons_token(cons):
    c = (cons or "").lower()
    if "pathogenic" in c:
        return "pathogenic"
    if "benign" in c:
        return "benign"
    if "deleterious" in c:
        return "deleterious"
    return "uncertain"
_LAYER_BUILT = {}     # obj -> {tag: {"pos": frozenset(built positions), "groups": last groups shown}}


# ----------------------------------------------------------------------------------------------
# Projection commands  (each is independently toggleable; pair with ufv_hide to remove)
# ----------------------------------------------------------------------------------------------


_SEVERITY = ["Likely pathogenic or pathogenic", "Predicted deleterious", "Uncertain significance",
             "Likely benign or benign"]
_PROV_TOKENS = {"reviewed": "UniProt reviewed", "uniprot": "UniProt reviewed", "clinvar": "ClinVar",
                "large": "Large scale studies", "largescale": "Large scale studies", "large_scale": "Large scale studies"}


def _parse_provenances(spec):
    """'reviewed clinvar' / 'UniProt reviewed;ClinVar' -> set of provenance categories (None = all)."""
    if spec in (None, "", "all"):
        return None
    if isinstance(spec, (set, list, tuple, frozenset)):
        return set(spec)
    out = set()
    for tok in re.split(r"[;,|]+|\s+", str(spec).strip()):
        t = tok.strip().lower()
        if t in _PROV_TOKENS:
            out.add(_PROV_TOKENS[t])
        elif t:
            for cat in PROVENANCE_COLORS:
                if cat.lower().startswith(t):
                    out.add(cat)
    return out


def _filter_variants(variants, tokens=None, provenances=None, diseases=None, reviewed_only=False):
    """data.js filterVariants: a variant is dropped when its (known) consequence is not active, its
    (known) provenance is not active, or — with a disease set — none of its diseases is in the set
    (an empty set shows nothing). `tokens` = consequence tokens (pathogenic/deleterious/benign/uncertain)."""
    if reviewed_only:
        provenances = {"UniProt reviewed"} if provenances is None else (set(provenances) & {"UniProt reviewed"})
    out = []
    for v in variants:
        if tokens and v.get("consequence") in CONSEQUENCE_COLORS and _cons_token(v["consequence"]) not in tokens:
            continue
        if provenances is not None and v.get("provenance") in PROVENANCE_COLORS and v["provenance"] not in provenances:
            continue
        if diseases is not None and not any(d in diseases for d in (v.get("diseases") or [])):
            continue
        out.append(v)
    return out


def _variant_colorer(ann, color="consequence", diseases=None):
    color = (color or "consequence").lower()
    if color.startswith("dis"):
        cmap = ann.get("diseaseColors") or {}

        def f(v):
            d = next((d for d in (v.get("diseases") or []) if not diseases or d in diseases), None)
            return cmap.get(d) or UNCLASSIFIED_COLOR if d else UNCLASSIFIED_COLOR
        return f
    if color.startswith("prov"):
        return lambda v: PROVENANCE_COLORS.get(v.get("provenance"), UNCLASSIFIED_COLOR)
    return lambda v: v.get("consequenceColor") or UNCLASSIFIED_COLOR


def _variant_groups(ann, tokens, reviewed_only=False, provenances=None, diseases=None, color="consequence"):
    """{colour: [positions]} for the variant sphere layer: filtered like the extension's variant window,
    ONE sphere per position coloured by its most severe variant (pathogenic > predicted deleterious >
    uncertain > benign) under the active colour profile (consequence | disease | provenance)."""
    shown = _filter_variants(ann["variants"], tokens, provenances, diseases, reviewed_only)
    col = _variant_colorer(ann, color, diseases)
    pos_col, pos_sev = {}, {}
    for v in shown:
        p = v["position"]
        sev = _SEVERITY.index(v["consequence"]) if v.get("consequence") in _SEVERITY else -1
        if p not in pos_col:
            pos_col[p], pos_sev[p] = col(v), sev
        elif sev >= 0 and (pos_sev[p] < 0 or sev < pos_sev[p]):
            pos_col[p], pos_sev[p] = col(v), sev
    groups = {}
    for p, c in pos_col.items():
        groups.setdefault(c, []).append(p)
    return groups


def _variant_positions(ann):
    """Every variant position (cached on the annotation dict): the variant layer is built over all of
    them once, so consequence / reviewed-only filters are an in-place re-colour + show/hide."""
    vp = ann.get("_varPositions")
    if vp is None:
        vp = ann["_varPositions"] = sorted({v["position"] for v in ann["variants"]})
    return vp


def _parse_diseases(ann, spec):
    """'AD1;CAA' (labels, ';' or '|' separated, case-insensitive) -> set of disease labels; None = all."""
    if spec in (None, "", "all"):
        return None
    if isinstance(spec, (set, list, tuple, frozenset)):
        return set(spec)
    if str(spec).strip().lower() == "none":
        return set()
    known = {d.lower(): d for d in (ann.get("diseaseColors") or {})}
    out = set()
    for tok in re.split(r"[;|]+", str(spec)):
        t = tok.strip()
        if t:
            out.add(known.get(t.lower(), t))
    return out


PARTNER_LAYERS = ("variants", "ptms", "sites")


def _partner_annotations(obj):
    """api.js loadPartnerAnnotations: for each partner protein of a complex (first PARTNER_ACCESSION_CAP
    accessions), its variants (disease colours), PTMs (category colours) and sites mapped onto ITS chains
    -> [{chainId, accession, byResi: {author resi: {uniPos, variants, ptms, sites}}}]. Network (call off
    the UI thread); cached on the object's map."""
    mp = UFV_MAPS.get(obj) or {}
    if "_partnerAnn" in mp:
        return mp["_partnerAnn"]
    groups = _partner_maps(obj)
    out = []
    if groups:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=min(DECORATE_CONCURRENCY, len(groups))) as ex:
            datas = list(ex.map(lambda g: _partner_protein(g[0], True), groups))
        for (acc, maps), d in zip(groups, datas):
            by_uni = {}

            def ensure(pos):
                return by_uni.setdefault(pos, {"variants": [], "ptms": [], "sites": []})
            for v in d.get("variants") or []:
                ensure(v["position"])["variants"].append(v)
            for p in d.get("ptms") or []:
                ensure(p["position"])["ptms"].append(p)
                if p.get("endPosition") and p["endPosition"] != p["position"]:
                    ensure(p["endPosition"])["ptms"].append(dict(p, position=p["endPosition"]))
            for s in d.get("sites") or []:
                for pos in range(s["position"], (s.get("endPosition") or s["position"]) + 1):
                    ensure(pos)["sites"].append(s)
            if not by_uni:
                continue
            for m in maps:
                try:
                    resolve = _partner_uni_to_author(mp.get("pdbId") or "", m)
                except Exception:
                    continue
                by_resi = {}
                for uni, data in by_uni.items():
                    if uni < m["uniprotStart"] or uni > m["uniprotEnd"]:
                        continue
                    a = resolve(uni)
                    if a is not None:
                        by_resi[int(a)] = dict(data, uniPos=uni)
                if by_resi:
                    out.append({"chainId": m["chainId"], "accession": acc, "byResi": by_resi})
    mp["_partnerAnn"] = out
    return out


def _partner_sphere_colors(ann_list, layers):
    """modal.js renderPartnerSpheres: one sphere per partner residue, coloured by the highest-priority
    enabled layer — pathogenic variants (disease colour of one with a named disease, else the most
    severe), then PTMs (category colour), then sites. -> {(chain, resi): (colour, tooltip)}"""
    out = {}
    for p in ann_list:
        for resi, info in p["byResi"].items():
            color = body = None
            if "variants" in layers:
                path = [v for v in info.get("variants") or [] if re.search(r"pathogenic|deleterious", v.get("consequence") or "", re.I)]
                if path:
                    top = next((v for v in path if v.get("diseaseColor") and v["diseaseColor"] != UNCLASSIFIED_COLOR), None) \
                        or next((v for sev in _SEVERITY for v in path if v.get("consequence") == sev), None) or path[0]
                    color = top.get("diseaseColor") or top.get("consequenceColor") or "#e53935"
                    body = "%d pathogenic variant%s" % (len(path), "" if len(path) == 1 else "s")
            if color is None and "ptms" in layers and info.get("ptms"):
                color, body = info["ptms"][0].get("color") or "#ab47bc", info["ptms"][0].get("category") or "PTM"
            if color is None and "sites" in layers and info.get("sites"):
                color, body = SITE_COLOR, info["sites"][0].get("description") or "Functional site"
            if color is not None:
                out[(p["chainId"], resi)] = (color, "%s · %s: %s (partner protein)" % (p["accession"], info["uniPos"], body))
    return out


def _parse_partner_layers(spec):
    s = (spec or "all").strip().lower()
    if s in ("all", "on", "1", ""):
        return set(PARTNER_LAYERS)
    if s in ("off", "none", "0"):
        return set()
    return {l for l in PARTNER_LAYERS if any(l.startswith(t[:3]) for t in re.split(r"[ ,;|]+", s) if t)}


_LIG_PRESENT = {}   # obj -> bool, cached so residue clicks don't re-scan for ligands each time
_AFF_META = {}                                 # obj -> {asym/auth chain: transplant meta} (AlphaFill models)


LIGAND_SIM_CAP = 20


# ----------------------------------------------------------------------------------------------
# Geometry + structure-dependent analyses (hotspots, contact hubs) + residue report / focus
# ----------------------------------------------------------------------------------------------
_GEOM_CACHE = {}      # obj -> geom list
_ANALYSIS_CACHE = {}  # (obj, uid, kind) -> tiers dict (hotspots / hubs / pockets)


def _invalidate_obj(obj):
    """Forget everything derived from `obj`'s coordinates or numbering: Cα geometry, ligand presence,
    cached analyses, and which marker layers were built (they are re-carved on next show). Call
    wherever the structure is (re)loaded, removed, or its residue mapping changes."""
    _GEOM_CACHE.pop(obj, None)
    _LIG_PRESENT.pop(obj, None)
    _LAYER_BUILT.pop(obj, None)
    _TRANSPLANT.pop(obj, None)
    _AFF_META.pop(obj, None)
    for k in [k for k in list(_ANALYSIS_CACHE) if k[0] == obj]:
        _ANALYSIS_CACHE.pop(k, None)
    _traj_forget(obj)                       # per-state stashes, coordinate arrays, trajectory results


def _invalidate_uid(uid):
    """Annotations for `uid` were re-fetched: drop analyses computed from the old ones."""
    for k in [k for k in list(_ANALYSIS_CACHE) if k[1] == uid]:
        _ANALYSIS_CACHE.pop(k, None)
    for stash in _TRAJ_STASH.values():      # ... including other trajectory states' cached analyses
        for _g, an in stash.values():
            for k in [k for k in an if k[1] == uid]:
                an.pop(k, None)
_PARTNER_DATA = {}            # (pdbId, accession) -> {"variants": [...], "ptms": [...], "sites": [...]}


def _partner_protein(acc, with_features=False):
    """Variants (and optionally PTMs + sites) of a partner accession, parsed like the extension's
    loadPartnerClassified / loadPartnerAnnotations (no AlphaMissense map). Network; memoised."""
    key = acc.upper()
    d = _PARTNER_DATA.get(key)
    if d is None or (with_features and "ptms" not in d):
        d = dict(d or {})
        if "variants" not in d:
            d["variants"] = _extract_variants(_get_json(VARIATION_URL.format(acc)), None)
            _compute_disease_colors(d["variants"])
        if with_features:
            feats = _get_json(FEATURES_URL.format(acc))
            d["ptms"] = _extract_ptms(feats)
            _group_ptm_colors(d["ptms"])
            d["sites"] = _extract_sites(feats)
        _PARTNER_DATA[key] = d
    return d


def _partner_maps(obj):
    """Partner chain mappings of `obj` grouped by accession (first PARTNER_ACCESSION_CAP accessions, in
    SIFTS order) -> [(accession, [mapping, ...])]. Only mappings with a UniProt range are usable."""
    mp = UFV_MAPS.get(obj) or {}
    by = {}
    for p in mp.get("partners") or []:
        if p.get("uniprotStart") is None or p.get("uniprotEnd") is None or not p.get("chainId"):
            continue
        by.setdefault(p["accession"], []).append(p)
    return list(by.items())[:PARTNER_ACCESSION_CAP]


def _partner_points(obj):
    """api.js loadPartnerClassified: pathogenic / benign-only residues of the partner proteins of a
    complex, as [{chainId, pdbResi, path}] in accession order. Network (call off the UI thread);
    cached on the map. [] for single-protein structures."""
    mp = UFV_MAPS.get(obj) or {}
    if "_partnerPoints" in mp:
        return mp["_partnerPoints"]
    groups = _partner_maps(obj)
    pts = []
    if groups:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=min(DECORATE_CONCURRENCY, len(groups))) as ex:
            datas = list(ex.map(lambda g: _partner_protein(g[0]), groups))
        for (acc, maps), d in zip(groups, datas):
            path, benign = _hotspot_classes(d.get("variants") or [])
            if not path and not benign:
                continue
            for m in maps:
                try:
                    resolve = _partner_uni_to_author(mp.get("pdbId") or "", m)
                except Exception:
                    continue
                for poss, is_path in ((path, True), (benign, False)):
                    for pos in poss:
                        if pos < m["uniprotStart"] or pos > m["uniprotEnd"]:
                            continue
                        a = resolve(pos)
                        if a is not None:
                            pts.append({"chainId": m["chainId"], "pdbResi": int(a), "path": is_path})
    mp["_partnerPoints"] = pts
    return pts


def _resolve_partner_points(geom, points):
    """[(is_pathogenic, geom entry)] for partner points that have a Cα in the model."""
    look = {}
    for g in geom:
        look.setdefault((g["chain"], g["resi"]), g)
    out = []
    for pp in points or []:
        g = look.get((pp["chainId"], pp["pdbResi"]))
        if g is not None:
            out.append((bool(pp["path"]), g))
    return out


def _tiers_for(kind, res):
    """Colourable {uniPos: tier} for an analysis result (pockets: q <= the sensitivity threshold)."""
    return _pocket_classes(res) if kind == "pockets" else (res or {})


def _set_pocket_q(q):
    try:
        q = float(q)
    except (TypeError, ValueError):
        return POCKET_SETTINGS["q"]
    POCKET_SETTINGS["q"] = round(min(0.40, max(0.01, q)), 4)
    return POCKET_SETTINGS["q"]


_PROTVAR_CACHE = {}


def _summ(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    return {"min": min(vals), "max": max(vals), "mean": sum(vals) / len(vals), "n": len(vals)}


AA_ALPHA = "ACDEFGHIKLMNPQRSTVWY"   # order ProtVar's per-substitution arrays follow (verified vs AM CSV)


def _protvar_by_mut(raw, wt):
    """ProtVar /score's EVE/ESM/AM arrays are the 19 non-wild-type substitutions in alphabetical
    order (verified against the labelled AlphaMissense CSV). Zip a raw array back to {mutantAA: value};
    bail (empty) on any length mismatch so we never mislabel."""
    if not raw or not wt:
        return {}
    muts = [a for a in AA_ALPHA if a != wt]
    if len(raw) != len(muts):
        return {}
    return dict(zip(muts, raw))


def _parse_protvar_score(arr):
    """ProtVar /score returns a flat array of typed objects: one CONSERV (per-position), 19 EVE and
    19 ESM and 19 AM (per-substitution, alphabetical mutant order — see _protvar_by_mut), and one M3D
    structural effect. Keep the raw per-substitution arrays so callers can label them, plus summaries."""
    conserv, eve, esm, eve_cls, m3d = None, [], [], [], None
    for o in arr or []:
        t = o.get("type")
        if t == "CONSERV":
            conserv = o.get("score")
        elif t == "EVE":
            eve.append(o.get("score"))
            eve_cls.append(o.get("eveClass"))
        elif t == "ESM":
            esm.append(o.get("score"))
        elif t == "M3D" and m3d is None:
            m3d = {"prediction": o.get("prediction"), "feature": o.get("damagingFeature")}
    return {"conservation": conserv, "m3d": m3d,
            "eveRaw": eve, "eveClsRaw": eve_cls, "esmRaw": esm,
            "eve": _summ(eve), "esm": _summ(esm),
            "evePath": sum(1 for c in eve_cls if str(c).upper() == "PATHOGENIC"),
            "eveN": sum(1 for c in eve_cls if c)}


def _parse_protvar_foldx(arr):
    """ProtVar /prediction/foldx returns per-substitution FoldX ΔΔG (kcal/mol) WITH wt/mut labels, so
    we can show the value for a specific variant and the per-position range. ΔΔG > ~2 destabilising."""
    by_mut, vals = {}, []
    for o in arr or []:
        m, d = o.get("mutatedType"), o.get("foldxDdg")
        if m and d is not None:
            by_mut[m] = d
            vals.append(d)
    return {"byMut": by_mut, "summary": _summ(vals)}


def _parse_protvar_pocket(arr):
    """ProtVar /prediction/pocket returns the predicted (potentially ligand-binding) pocket(s) that
    a residue lines — buriedness, energy, a druggability-like score, and the lining residues. Returns
    a summary {count, buriedness, score, plddt} for the best-scoring pocket, or None."""
    pockets = arr if isinstance(arr, list) else ([arr] if arr else [])
    pockets = [p for p in pockets if isinstance(p, dict)]
    if not pockets:
        return None
    best = max(pockets, key=lambda p: p.get("score") or 0)
    lining = sorted({int(r) for p in pockets for r in (p.get("resid") or []) if str(r).lstrip("-").isdigit()})
    return {"count": len(pockets), "buriedness": best.get("buriedness"),
            "score": best.get("score"), "plddt": best.get("meanPlddt"), "lining": lining,
            "pockets": _normalize_pockets(pockets)}


def _finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and x == x and abs(x) != float("inf")


def _normalize_pockets(arr):
    """api.js normalizePockets: every ProtVar (AutoSite) pocket a residue lines, with its id, buriedness,
    score, mean pLDDT, radius of gyration, energy per volume and lining residues; best score first."""
    raw = arr if isinstance(arr, list) else ([arr] if arr else [])
    out = [{"pocketId": o.get("pocketId"), "buriedness": o.get("buriedness"), "score": o.get("score"),
            "meanPlddt": o.get("meanPlddt"), "radGyration": o.get("radGyration"),
            "energyPerVol": o.get("energyPerVol"),
            "resid": [r for r in (o.get("resid") or []) if _finite(r)] if isinstance(o.get("resid"), list) else []}
           for o in raw if isinstance(o, dict)]
    out.sort(key=lambda p: -(p["score"] or 0))        # stable, like Array.prototype.sort
    return out


def _canon_acc(acc):
    return re.sub(r"-\d+$", "", str(acc or "").strip())


def _cadd_by_mut(acc, pos, wt=None):
    """api.js fetchCaddByMut: CADD (phred) per mutant amino acid for one position from ProtVar /mapping
    (ONE call returns every SNV-reachable substitution of the codon; the most deleterious codon path wins).
    -> {mutant 1-letter: caddScore}."""
    ref = wt if (wt and re.match(r"^[A-Z]$", wt)) else "A"     # ProtVar corrects a wrong ref itself
    j = _get_json(PROTVAR_MAPPING.format(_uri("%s %d %s *" % (acc, int(pos), ref))), True)
    inputs = (((j or {}).get("content") or {}).get("inputs") or []) if isinstance(j, dict) else []
    gvs = ((inputs[0] or {}).get("derivedGenomicVariants") or []) if inputs else []
    out = {}
    for gv in gvs:
        genes = (gv or {}).get("genes") or []
        g = genes[0] if genes else None
        if not g or g.get("caddScore") is None:
            continue
        isos = g.get("isoforms") or []
        iso = next((i for i in isos if i.get("canonical")), None) or (isos[0] if isos else None)
        aa1 = AA3_TO_1.get((iso or {}).get("variantAA"))
        if not aa1:
            continue
        try:
            c = float(g["caddScore"])
        except (TypeError, ValueError):
            continue
        if out.get(aa1) is None or c > out[aa1]:
            out[aa1] = c
    return out


def fetch_protvar(uid, pos, wt=None):
    """ProtVar per-residue predictors (conservation / EVE / ESM1b / M3D / FoldX / binding pockets / CADD)
    for a position — four parallel requests, as api.js fetchProtVar. Network — call off the UI thread.
    Cached per (canonical accession, position)."""
    acc = _canon_acc(uid)
    key = (acc, int(pos))
    if key in _PROTVAR_CACHE:
        return _PROTVAR_CACHE[key]
    if wt is None:
        seq = (_CACHE.get(uid) or _CACHE.get(acc) or {}).get("sequence") or ""
        wt = seq[int(pos) - 1] if 0 < int(pos) <= len(seq) else None
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as ex:
        # disk-cached like every other per-residue response (a residue revisited in a later session is
        # read back in ~1 ms instead of 3 network round trips)
        fs = ex.submit(_get_json, PROTVAR_SCORE.format(acc, pos), True)
        ff = ex.submit(_get_json, PROTVAR_FOLDX.format(acc, pos), True)
        fp = ex.submit(_get_json, PROTVAR_POCKET.format(acc, pos), True)   # shared with the pocket crawl
        fc = ex.submit(_cadd_by_mut, acc, pos, wt)
        try:
            cadd = fc.result() or {}
        except Exception:
            cadd = {}
        score = _parse_protvar_score(fs.result())
        score["caddByMut"] = cadd
        out = {"score": score,
               "foldx": _parse_protvar_foldx(ff.result()),
               "pocket": _parse_protvar_pocket(fp.result())}
    _PROTVAR_CACHE[key] = out
    return out


def _protvar_csv_summary(arr):
    """Per-position ProtVar /score summary for the CSV: conservation, mean EVE, most-deleterious
    ESM1b, and an M3D-damaging flag. None if the position has no scores."""
    if not arr:
        return None
    conserv, eve, esm, m3d = None, [], [], ""
    for o in arr:
        t = o.get("type")
        if t == "CONSERV":
            conserv = o.get("score")
        elif t == "EVE":
            eve.append(o.get("score"))
        elif t == "ESM":
            esm.append(o.get("score"))
        elif t == "M3D" and m3d == "":
            m3d = 1 if str(o.get("prediction", "")).lower().startswith("damag") else 0
    eve = [x for x in eve if isinstance(x, (int, float))]
    esm = [x for x in esm if isinstance(x, (int, float))]
    if conserv is None and not eve and not esm and m3d == "":
        return None
    return {"conservation": conserv,
            "eve_mean": sum(eve) / len(eve) if eve else None,
            "esm_min": min(esm) if esm else None,
            "m3d": m3d}


def _fetch_protvar_all(uid, length, progress=None):
    """Bulk ProtVar /score for every residue 1..length (one request per position, bounded
    concurrency) -> {pos: summary}. Slow — no bulk endpoint exists — so it's opt-in. `progress`
    (done, total) for status updates."""
    from concurrent.futures import ThreadPoolExecutor

    def one(pos):
        try:
            return pos, _protvar_csv_summary(_get_json(PROTVAR_SCORE.format(uid, pos), True))
        except Exception:
            return pos, None
    out, done = {}, 0
    with ThreadPoolExecutor(max_workers=8) as ex:
        for pos, summ in ex.map(one, range(1, length + 1)):
            if summ:
                out[pos] = summ
            done += 1
            if progress and (done % 25 == 0 or done == length):
                progress(done, length)
    return out


# ----------------------------------------------------------------------------------------------
# Services (WP5-WP7): ProtVar pockets crawl, PDBe-KB binding context, Open Targets tractability,
# AlphaFill transplants, RCSB strucmotif search, pocket evidence / ligand clash, Nearby, flags,
# B-factor export. Network helpers are pure (safe in worker threads); cmd.* only in the commands.
# ----------------------------------------------------------------------------------------------
POCKET_PALETTE = ["#26c6da", "#ff7043", "#66bb6a", "#ab47bc", "#ffa726", "#ec407a", "#29b6f6", "#9ccc65"]
POCKET_CRAWL_CONCURRENCY = 8
_ALL_POCKETS = {}     # (canonical accession, length) -> {"count", "pockets"} or None
_ALLPK_KEY = "ufv-aggregate:protvar-all-pockets/{}/{}"   # disk-cache key of a finished crawl (never fetched)


def fetch_all_pockets(uid, length, progress=None):
    """api.js fetchAllProtVarPockets: ProtVar has no protein-level pocket endpoint, so crawl
    /pocket/{acc}/{pos} over 1..length (8 in flight; each response is disk-cached, 404 = no pocket),
    de-duplicate by pocketId and sort best score first. -> {"count", "pockets"} or None. Network."""
    acc = _canon_acc(uid)
    key = (acc, int(length or 0))
    if key in _ALL_POCKETS:
        return _ALL_POCKETS[key]
    # A finished crawl is also stored as ONE disk-cache entry, so a new session reads 1 file instead of
    # `length` (2000+ small gzip files: 0.5-2.5 s on Windows). Only used with the real fetcher, and only
    # written when every residue answered (a network gap must be retried, never baked in for 7 days).
    agg = _ALLPK_KEY.format(acc, key[1]) if _fetch_bytes is _FETCH_BYTES_REAL else None
    if agg:
        raw = _cache_read(agg)
        if raw is not None:
            try:
                res = json.loads(raw.decode("utf-8"))
                if res is None or (isinstance(res, dict) and isinstance(res.get("pockets"), list)):
                    _ALL_POCKETS[key] = res
                    if progress:
                        progress(key[1], key[1])
                    return res
            except Exception:
                pass
    failed = []

    def one(pos):
        try:
            raw, data = _get_json_raw(PROTVAR_POCKET.format(acc, pos), True)
            if raw is None:
                failed.append(pos)
            return _normalize_pockets(data)
        except Exception:
            failed.append(pos)
            return []
    by_id, done = {}, 0
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=POCKET_CRAWL_CONCURRENCY) as ex:
        for pk_list in ex.map(one, range(1, key[1] + 1)):
            for pk in pk_list:
                if pk["pocketId"] is not None and pk["pocketId"] not in by_id:
                    by_id[pk["pocketId"]] = pk
            done += 1
            if progress and (done % 25 == 0 or done == key[1]):
                progress(done, key[1])
    pockets = sorted(by_id.values(), key=lambda p: -(p["score"] or 0))
    res = {"count": len(pockets), "pockets": pockets} if pockets else None
    _ALL_POCKETS[key] = res
    if agg and not failed and key[1] > 0:
        try:
            _cache_write(agg, json.dumps(res).encode("utf-8"))
        except Exception:
            pass
    return res


def _find_pocket(uid, pocket_id):
    """A ProtVar pocket already fetched (protein-wide crawl or any residue report) by its pocketId."""
    acc = _canon_acc(uid)
    sid = str(pocket_id)
    pools = [v for k, v in _ALL_POCKETS.items() if k[0] == acc and v]
    pools += [{"pockets": (v.get("pocket") or {}).get("pockets") or []}
              for k, v in _PROTVAR_CACHE.items() if k[0] == acc and v]
    for pool in pools:
        for p in pool.get("pockets") or []:
            if str(p.get("pocketId")) == sid:
                return p
    return None


def all_pockets_html(res, obj=None, uid=None):
    """Panel listing of every ProtVar pocket (modal.js buildAllPocketsPanel): palette swatch, stats,
    show / motif-search links (ufv:apk anchors)."""
    if not res or not res.get("pockets"):
        return "<i>No pockets detected across this protein.</i>"
    pl = res["pockets"]
    motif_ok = bool(obj) and _structure_motif_ok(obj)
    T = ['<table cellspacing="0" cellpadding="3" style="font-family:sans-serif; font-size:12px;">',
         '<tr><td colspan="2" style="font-size:14px;"><b style="color:#00897b;">%d pocket%s detected</b> '
         '<a href="ufv:apk:overlay" style="color:#1565c0; text-decoration:none; font-size:11px;">overlay all</a></td></tr>'
         % (len(pl), "" if len(pl) == 1 else "s")]
    for i, p in enumerate(pl):
        c = POCKET_PALETTE[i % len(POCKET_PALETTE)]
        stats = ["%d res" % len(p["resid"])]
        if _finite(p.get("score")):
            stats.append("score %d" % round(p["score"]))
        if _finite(p.get("meanPlddt")):
            stats.append("pLDDT %d" % round(p["meanPlddt"]))
        if _finite(p.get("radGyration")):
            stats.append("Rg %.1f" % p["radGyration"])
        if _finite(p.get("energyPerVol")):
            stats.append("E/V %.2f" % p["energyPerVol"])
        links = ' <a href="ufv:apk:%d" style="color:#fb8c00; text-decoration:none;">show ▸</a>' % i
        if motif_ok and len(p["resid"]) >= 2:
            links += ' <a href="ufv:apk:motif:%d" style="color:#1565c0; text-decoration:none;">motifs ↗</a>' % i
        T.append('<tr><td style="white-space:nowrap;"><span style="color:%s; font-size:14px;">■</span> '
                 '<b>Pocket %s</b></td><td>%s%s</td></tr>' % (c, _h(p.get("pocketId")), " · ".join(stats), links))
    T.append("</table>")
    return "".join(T)


# ---- PDBe-KB experimental binding context ------------------------------------------------------
_PDBEKB_CACHE = {}    # canonical accession -> {"ligandByPos", "interfaceByPos"} or None


def fetch_pdbekb(uid):
    """api.js fetchPdbeKbContext: PDBe-KB ligand_sites + interface_residues (two parallel requests,
    disk-cached), inverted to per-UniProt-position maps: ligandByPos {pos: [{code, name, drugbank, chembl,
    pdb, pdbCount}]} and interfaceByPos {pos: [{accession, name, pdb, pdbCount}]}. None when empty."""
    canon = _canon_acc(uid)
    if canon in _PDBEKB_CACHE:
        return _PDBEKB_CACHE[canon]
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=2) as ex:
        fl = ex.submit(_get_json, PDBEKB_URL.format("ligand_sites", canon))
        fi = ex.submit(_get_json, PDBEKB_URL.format("interface_residues", canon))
        lig, iface = fl.result(), fi.result()

    def add_range(m, item, residues):
        for r in residues or []:
            a = r.get("startIndex")
            b = r.get("endIndex") if r.get("endIndex") is not None else a
            if not isinstance(a, int) or not isinstance(b, int):
                continue
            entries = r.get("interactingPDBEntries") or r.get("allPDBEntries") or []
            e0 = entries[0] if entries else None
            pdb = (e0.get("pdbId") if isinstance(e0, dict) else e0) or None
            for pos in range(a, b + 1):
                m.setdefault(pos, []).append(dict(item, pdb=pdb, pdbCount=len(entries)))
    lbp, ibp = {}, {}
    for d in ((((lig or {}).get(canon) or {}).get("data")) or []) if isinstance(lig, dict) else []:
        ad = d.get("additionalData") or {}
        add_range(lbp, {"code": d.get("accession"), "name": d.get("name"), "drugbank": ad.get("drugBankId") or None,
                        "chembl": ad.get("chemblId") or None}, d.get("residues"))
    for d in ((((iface or {}).get(canon) or {}).get("data")) or []) if isinstance(iface, dict) else []:
        add_range(ibp, {"accession": d.get("accession"), "name": d.get("name")}, d.get("residues"))
    res = {"ligandByPos": lbp, "interfaceByPos": ibp} if (lbp or ibp) else None
    if lig is not None or iface is not None:          # don't memoise a network failure
        _PDBEKB_CACHE[canon] = res
    return res


def _kb_top(items, key):
    """modal.js renderStructuralContext: de-duplicate by `key` keeping the most-observed entry, most
    observed first (stable)."""
    seen = {}
    for it in items or []:
        e = seen.get(it.get(key))
        if e is None or (it.get("pdbCount") or 0) > (e.get("pdbCount") or 0):
            seen[it.get(key)] = it
    return sorted(seen.values(), key=lambda x: -(x.get("pdbCount") or 0))


def _kb_residue(kb, pos):
    """(ligands top-8 list, total, interfaces top-4 list, total) for one residue."""
    if not kb:
        return [], 0, [], 0
    lt = _kb_top(kb["ligandByPos"].get(int(pos)), "code")
    it = _kb_top(kb["interfaceByPos"].get(int(pos)), "accession")
    return lt[:8], len(lt), it[:4], len(it)


# ---- Open Targets tractability -----------------------------------------------------------------
_OT_CACHE = {}
OT_MODALITY = {"SM": "Small molecule", "AB": "Antibody", "PR": "PROTAC", "OC": "Other"}


def _post_json(url, payload, cache=True):
    """POST JSON -> parsed JSON or None; disk-cached under url + body hash (same TTL as GETs)."""
    body = json.dumps(payload, separators=(",", ":")).encode("utf-8")   # = JSON.stringify
    key = url + "#POST#" + _sha1hex(body)
    raw = _cache_read(key) if cache else None
    if raw is None:
        try:
            req = _urlreq().Request(url, data=body, headers={
                "User-Agent": _USER_AGENT, "Content-Type": "application/json", "Accept": "application/json"})
            with _urlreq().urlopen(req, timeout=60, context=_ssl_default()) as r:
                raw = r.read()
            json.loads(raw.decode("utf-8"))
            if cache:
                _cache_write(key, raw)
        except Exception:
            return None
    try:
        return json.loads(raw.decode("utf-8"))
    except Exception:
        return None


def _ensembl_gene(uid):
    """api.js ensemblGeneForUniProt: Ensembl gene id (version stripped) from UniProt's own xref."""
    j = _get_json(UNIPROT_URL.format(_canon_acc(uid)) + "?fields=xref_ensembl")
    for r in (j or {}).get("uniProtKBCrossReferences") or []:
        if r.get("database") != "Ensembl":
            continue
        g = next((x for x in (r.get("properties") or []) if x.get("key") == "GeneId"), None)
        if g and g.get("value"):
            return g["value"].split(".")[0]
    return None


def fetch_open_targets(uid):
    """api.js fetchOpenTargets: target-level tractability (small molecule / antibody / PROTAC / other),
    POSITIVE evidence only, grouped by modality. -> {symbol, ensemblId, groups:[{modality, labels}], url}
    or None. Network (UniProt xref + one GraphQL POST), cached."""
    canon = _canon_acc(uid)
    if canon in _OT_CACHE:
        return _OT_CACHE[canon]
    ensg = _ensembl_gene(canon)
    res = None
    if ensg:
        q = '{ target(ensemblId:"%s"){ approvedSymbol tractability{ label modality value } } }' % ensg
        j = _post_json(OT_GQL, {"query": q})
        t = (((j or {}).get("data") or {}).get("target")) if isinstance(j, dict) else None
        if t:
            by_mod = {}
            for x in t.get("tractability") or []:
                if x.get("value"):
                    by_mod.setdefault(x.get("modality"), []).append(x.get("label"))
            res = {"symbol": t.get("approvedSymbol"), "ensemblId": ensg,
                   "groups": [{"modality": OT_MODALITY.get(m, m), "labels": lb} for m, lb in by_mod.items()],
                   "url": "https://platform.opentargets.org/target/%s" % ensg}
    _OT_CACHE[canon] = res
    return res


def _h(x):
    return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def tractability_html(ot):
    """Rich-text line for the panel: modality buckets + link (modal.js updateTractabilityNav)."""
    if not ot or not ot.get("groups"):
        return ""
    parts = ["<b>%s:</b> %s" % (_h(g["modality"]), _h(", ".join(g["labels"]))) for g in ot["groups"]]
    return ('<b style="color:#00695c;">Tractability</b> — %s · <a href="%s">%s on Open Targets ↗</a>'
            % (" · ".join(parts), _h(ot["url"]), _h(ot.get("symbol") or "target")))


# ---- AlphaFill transplants ---------------------------------------------------------------------
_AFF_CACHE = {}


def alphafill_transplants(aff_id):
    """api.js getAlphaFillTransplants: {transplant asym_id (= the loaded model's auth chain): {identity,
    clash, compound, donorPdb, rmsd}} from alphafill.eu, or None. Network, cached."""
    if not aff_id:
        return None
    key = str(aff_id).upper()
    if key in _AFF_CACHE:
        return _AFF_CACHE[key]
    j = _get_json(ALPHAFILL_JSON.format(_uri(key)))
    out = None
    if isinstance(j, dict) and isinstance(j.get("hits"), list):
        out = {}
        for h in j["hits"]:
            al = h.get("alignment") or {}
            ident = _num_or_none(al.get("identity"))
            for t in h.get("transplants") or []:
                if not t.get("asym_id"):
                    continue
                cl = t.get("clash") or {}
                out[str(t["asym_id"])] = {"identity": ident, "clash": _num_or_none(cl.get("score")),
                                          "compound": t.get("compound_id"), "donorPdb": h.get("pdb_id"),
                                          "rmsd": _num_or_none(t.get("local_rmsd"))}
    if j is not None:
        _AFF_CACHE[key] = out
    return out


def _obj_struct(obj):
    return (_STATE.get("obj_struct") or {}).get(obj) or {}


def _is_alphafill(obj):
    st = _obj_struct(obj)
    return bool(re.search("alphafill", "%s %s %s" % (st.get("provider") or "", st.get("source") or "",
                                                     st.get("label") or ""), re.I))


def _is_af_based(obj):
    """modal.js isAlphaFoldBased: AlphaFold or AlphaFill model (no PDB id)."""
    st = _obj_struct(obj)
    return not st.get("pdbId") and (st.get("source") == "AlphaFold" or _is_alphafill(obj))


_TRANSPLANT = {}      # obj -> frozenset of AlphaFill donor/transplant chains (AlphaFold-based models only)


# ---- RCSB Structure Motif search ---------------------------------------------------------------
_LABEL_IDS = {}


def _pocket_label_ids(pdb_id):
    """api.js fetchPocketLabelIds: {"authChain|authResi": (label_asym_id, label_seq_id)} from PDBe
    residue_listing (struct_asym_id / residue_number). The model file's label_seq_id is NOT used: PDBe's
    CIF numbering differs from RCSB's deposited mmCIF."""
    pid = str(pdb_id or "").lower()
    if not pid:
        return {}
    if pid in _LABEL_IDS:
        return _LABEL_IDS[pid]
    j = _memo_json(PDBe_RESIDUES.format(pid))
    out = {}
    for mol in (((j or {}).get(pid) or {}).get("molecules") or []):
        for ch in mol.get("chains") or []:
            la, auth = ch.get("struct_asym_id"), ch.get("chain_id")
            if not la:
                continue
            for r in ch.get("residues") or []:
                if r.get("author_residue_number") is None or r.get("residue_number") is None:
                    continue
                out["%s|%s" % (auth, r["author_residue_number"])] = (la, r["residue_number"])
    if j is not None:
        _LABEL_IDS[pid] = out
    return out


def _structure_motif_ok(obj):
    """A strucmotif search is possible on PDB entries and AlphaFold / AlphaFill models."""
    return bool(_obj_struct(obj).get("pdbId") or (UFV_MAPS.get(obj) or {}).get("pdbId") or _is_af_based(obj))


def _open_url(url):
    try:
        import webbrowser
        webbrowser.open(url)
        return True
    except Exception:
        return False


# ---- pocket evidence / ligand clash (analysis.js pocketEvidence, ligandClash) -------------------
_AA_CLASS = {}
for _a in "AVLIMFWPGC":
    _AA_CLASS[_a] = "hydrophobic"
for _a in "STNQY":
    _AA_CLASS[_a] = "polar"
for _a in "DE":
    _AA_CLASS[_a] = "acidic"
for _a in "KRH":
    _AA_CLASS[_a] = "basic"


def pocket_composition(positions, sequence):
    if not sequence:
        return None
    comp = {"hydrophobic": 0, "polar": 0, "acidic": 0, "basic": 0, "aromatic": 0, "total": 0}
    for p in dict.fromkeys(int(x) for x in positions or []):
        aa = sequence[p - 1] if 0 < p <= len(sequence) else None
        if not aa:
            continue
        comp["total"] += 1
        c = _AA_CLASS.get(aa)
        if c:
            comp[c] += 1
        if aa in "FWYH":
            comp["aromatic"] += 1
    return comp if comp["total"] else None


def pocket_evidence(positions, geom, sites, sequence):
    """analysis.js pocketEvidence: mean Cα B-factor over the pocket residues (= mean pLDDT on AlphaFold /
    AlphaFill models), UniProt functional sites overlapping the pocket, and residue-class composition."""
    pset = {int(p) for p in positions or []}
    vals = [g["b"] for g in geom or [] if g.get("uniPos") is not None and g["uniPos"] in pset
            and isinstance(g.get("b"), (int, float))]
    hits = []
    for s in sites or []:
        end = s.get("endPosition") or s["position"]
        if any(p in pset for p in range(s["position"], end + 1)):
            hits.append(s)
    return {"meanB": (sum(vals) / len(vals)) if vals else None, "count": len(vals), "siteHits": hits,
            "composition": pocket_composition(pset, sequence)}


def composition_text(comp):
    if not comp:
        return ""
    t = comp["total"] or 1
    return " · ".join("%d%% %s" % (round(100.0 * comp[k] / t), k)
                      for k in ("hydrophobic", "aromatic", "acidic", "basic", "polar"))


SPIDER_LABELS = ("Hydrophobic", "Aromatic", "Acidic", "Basic", "Polar")


def composition_axes(comp):
    """Pocket-character radar values: the five residue-class fractions (independent, not a partition)."""
    if not comp:
        return []
    t = comp["total"] or 1
    return [comp[k] / t for k in ("hydrophobic", "aromatic", "acidic", "basic", "polar")]


def spider_img_html(comp, size=168):
    """<img> for the pocket-character radar. PyMOL's report browser resolves the ufvspider: URL to a drawn image."""
    ax = composition_axes(comp)
    if not ax:
        return ""
    return ('<img src="ufvspider:%d:%s" width="%d" height="%d" title="%s">'
            % (size, ",".join("%.4f" % a for a in ax), size, size, composition_text(comp)))


def plddt_band(b):
    return "very high" if b >= 90 else "confident" if b >= 70 else "low" if b >= 50 else "very low"


PLDDT_BAND_COLORS = {"very high": "#0053d6", "confident": "#4b9be0", "low": "#d4a017", "very low": "#ef7c3c"}


def aff_clash_band(c):
    return "low" if c <= 0.5 else "moderate" if c <= 0.9 else "high"


def residue_report_colors(rep):
    """{pos: variant consequence colour} for the report's Nearby chips (last variant at a position wins)."""
    want = set((rep.get("nearby") or {}).get("near") or [])
    out = {}
    for v in (_CACHE.get(rep.get("uid")) or {}).get("variants") or []:
        if v["position"] in want:
            out[v["position"]] = v.get("consequenceColor")
    return out


def _displayed_positions(obj, layer=None):
    out = set()
    for tag, groups in (_SPHERE_STATE.get(obj) or {}).items():
        if layer and tag != _HIDE_POINT_ALIASES.get(layer, layer):
            continue
        for positions in (groups or {}).values():
            out.update(int(p) for p in positions)
    return sorted(out)


def _site_cat(desc):
    d = (desc or "").lower()
    if "active" in d:
        return "active"
    if "metal" in d:
        return "metal"
    if "binding" in d:
        return "binding"
    return "other"


def _js_fixed(x, d):
    """Number.prototype.toFixed(d): round half up on the exact binary value (Python's % rounds half to
    even, which differs on exact ties)."""
    s = "%.*f" % (d, x)
    r = abs(x) * 10 ** d
    if abs(r - _math.floor(r) - 0.5) > 1e-6:      # not near a rounding tie: % already agrees with JS
        return s
    from decimal import Decimal, ROUND_HALF_UP
    try:
        return str(Decimal(x).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP))
    except Exception:
        return s


def _js_exp3(x):
    """Number.prototype.toExponential(3): 1.234e-5 / 1.234e+2 (no zero-padded exponent)."""
    if x == 0:
        return "0.000e+0"
    t = "%.10e" % abs(x)
    if not (t[5:11].startswith(("500000", "499999"))):        # not a rounding tie: % agrees with JS
        mant, ex = ("%.3e" % x).split("e")
        return "%se%s%d" % (mant, "+" if int(ex) >= 0 else "-", abs(int(ex)))
    from decimal import Decimal, ROUND_HALF_UP
    try:
        dx = Decimal(x)
        e = dx.copy_abs().adjusted()
        m = (dx.scaleb(-e)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
        if abs(m) >= 10:
            e += 1
            m = (dx.scaleb(-e)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
        return "%se%s%d" % (m, "+" if e >= 0 else "-", abs(e))
    except Exception:
        s = "%.3e" % x
        mant, ex = s.split("e")
        return "%se%s%d" % (mant, "+" if int(ex) >= 0 else "-", abs(int(ex)))


def _residue_csv_text(ann, obj=None, protvar_by_pos=None, chains=(), hot=None, hub=None, pocket=None,
                      prox=None, hot_ch=None, hub_ch=None, lig_contacts=None):
    """The per-residue CSV text (see _residue_csv) from structure results already computed: `chains`
    (mapped chains of `obj`), hotspot / hub tiers (+ per chain), constraint-pocket byPos, PTM-variant
    proximity and per-residue ligand contacts. obj=None -> sequence-only matrix (the VMD --csv backend)."""
    hot, hub, pocket, prox = hot or {}, hub or {}, pocket or {}, prox or {}
    hot_ch, hub_ch = hot_ch or {}, hub_ch or {}
    chains = list(chains or [])
    seq = ann.get("sequence", "")
    if not seq:
        return ""
    uid = ann.get("uid")
    ptm_cats = sorted({p["category"] for p in ann["ptms"]})
    diseases = sorted({p["label"] for v in ann["variants"] for p in (v.get("diseasePairs") or []) if p.get("label")})
    am_map = ann.get("amMap") or {}

    ptm_by_pos, disease_by_pos, site_by_pos = {}, {}, {}
    nvar, npath, gnomad = {}, {}, {}
    for p in ann["ptms"]:
        for pos in (p["position"], p.get("endPosition")):
            if pos:
                ptm_by_pos.setdefault(pos, set()).add(p["category"])
    for v in ann["variants"]:
        pos = v["position"]
        nvar[pos] = nvar.get(pos, 0) + 1
        if "pathogenic" in ("%s %s" % (v.get("consequence") or "", v.get("clinVar") or "")).lower():
            npath[pos] = npath.get(pos, 0) + 1
        af = v.get("gnomad")
        if isinstance(af, (int, float)) and af == af and abs(af) != float("inf"):
            gnomad[pos] = max(gnomad.get(pos, 0.0), af)
        ds = disease_by_pos.setdefault(pos, set())
        for p in v.get("diseasePairs") or []:
            if p.get("label"):
                ds.add(p["label"])
    for s in ann["sites"]:
        for pos in range(s["position"], (s.get("endPosition") or s["position"]) + 1):
            site_by_pos.setdefault(pos, set()).add(_site_cat(s["description"]))
    mut_pos = set()
    for mtg in ann.get("mutagenesis", []):
        mut_pos.update(range(mtg["position"], (mtg.get("endPosition") or mtg["position"]) + 1))
    burden = ann.get("burden") or set()

    lig_contacts = lig_contacts or {}
    multi = len(chains) > 1
    hot_num, hub_num = {"strong": 3, "moderate": 2, "weak": 1}, {"strong": 2, "moderate": 1}
    safe = lambda x: re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_]", "_", x))
    pv_headers = ["protvar_conservation", "protvar_eve_mean", "protvar_esm_min", "protvar_m3d_damaging"] if protvar_by_pos else []

    def pv_cols(pos):
        if not protvar_by_pos:
            return []
        d = protvar_by_pos.get(pos)
        if not d:
            return ["", "", "", ""]
        return [_js_fixed(d["conservation"], 3) if d.get("conservation") is not None else "",
                _js_fixed(d["eve_mean"], 3) if d.get("eve_mean") is not None else "",
                _js_fixed(d["esm_min"], 2) if d.get("esm_min") is not None else "",
                d["m3d"] if d.get("m3d") not in ("", None) else ""]

    def prox_cols(pos):
        p = prox.get(pos)
        if not p:
            return ["", "", "", "", ""]
        return [p["tier"] if p.get("tier") is not None else "", p.get("nearestVariant") or "",
                _js_fixed(p["nearestDist"], 1) if p.get("nearestDist") is not None else "",
                p["nearbyCount8A"], p["pathCount8A"]]

    def resi_for(pos, ch):
        if not obj:
            return pos
        r = _uni_to_resi(obj, ch, pos)
        return "" if r is None else r
    prox_headers = ["ptm_variant_tier", "nearest_variant_to_ptm", "ptm_variant_distance_A",
                    "nearby_variant_count_8A", "nearby_pathogenic_count_8A"]
    feat_headers = ["site", "site_active", "site_binding", "site_metal", "mutagenesis",
                    "n_variants", "n_pathogenic", "gnomad_max_af"]
    head = (["ptm_" + safe(c) for c in ptm_cats] + ["disease_" + safe(d) for d in diseases]
            + feat_headers + pv_headers + ["am_avg_score", "am_max_score", "am_n_subs"])
    if multi:
        headers = (["position", "aa"] + head + ["residue_burden", "constraint_pocket_q", "constraint_pocket_class"]
                   + prox_headers + ["nearby_ligands"]
                   + [h for c in chains for h in ("pdb_residue_%s" % c, "hotspot_tier_%s" % c, "contact_hub_tier_%s" % c)])
    else:
        headers = (["position", "aa", "pdb_residue"] + head
                   + ["hotspot_tier", "contact_hub_tier", "residue_burden", "constraint_pocket_q",
                      "constraint_pocket_class"] + prox_headers + ["nearby_ligands"])
    rows = [",".join(headers)]
    am_csv = ann.get("_csvAm")              # per-position AM columns, formatted once per annotation set
    if am_csv is None:
        am_csv = {}
        if am_map:
            abp = ann.get("amByPos") or _am_by_pos(am_map)
            for i, aa in enumerate(seq):
                ent = abp.get(i + 1)
                if not ent:
                    continue
                d_ = {m: s for (w, m, s) in ent if w == aa and m != aa and s == s and abs(s) != float("inf")}
                scores = [d_[m] for m in AM_AAS if m in d_]       # summed in AM_AAS order, as the JS does
                if scores:
                    am_csv[i + 1] = "%s,%s,%d" % (_js_fixed(sum(scores) / len(scores), 4), _js_fixed(max(scores), 4),
                                                  len(scores))
        ann["_csvAm"] = am_csv
    # one-hot blocks: most residues carry no PTM / disease -> one shared all-zero string
    zero_ptm = ",".join("0" * len(ptm_cats))
    zero_dis = ",".join("0" * len(diseases))

    def onehot(cats, have, zero):
        if not have:
            return zero
        return ",".join("1" if c in have else "0" for c in cats)
    for i, aa in enumerate(seq):
        pos = i + 1
        sc = site_by_pos.get(pos)
        af = gnomad.get(pos)
        pk = pocket.get(pos)
        parts = [str(pos), aa]
        if not multi:
            parts.append(str(resi_for(pos, chains[0]) if (obj and chains) else ("" if obj else pos)))
        if ptm_cats:
            parts.append(onehot(ptm_cats, ptm_by_pos.get(pos), zero_ptm))
        if diseases:
            parts.append(onehot(diseases, disease_by_pos.get(pos), zero_dis))
        parts.append("%d,%d,%d,%d,%d,%d,%d,%s" % (
            1 if sc else 0, 1 if sc and "active" in sc else 0, 1 if sc and "binding" in sc else 0,
            1 if sc and "metal" in sc else 0, 1 if pos in mut_pos else 0,
            nvar.get(pos, 0), npath.get(pos, 0), _js_exp3(af) if af is not None else ""))
        if protvar_by_pos:
            parts.append(",".join(str(x) for x in pv_cols(pos)))
        parts.append(am_csv.get(pos) or ",,")
        tail = ["%d" % (1 if pos in burden else 0),
                (_js_fixed(pk["q"], 4) + "," + pk["cat"]) if pk else ",",
                ",".join(str(x) for x in prox_cols(pos)) if pos in prox else ",,,,",
                ";".join(lig_contacts[pos]) if pos in lig_contacts else ""]
        if multi:
            parts.extend(tail)
            for c in chains:
                parts.append("%s,%d,%d" % (resi_for(pos, c), hot_num.get((hot_ch.get(c) or {}).get(pos), 0),
                                           hub_num.get((hub_ch.get(c) or {}).get(pos), 0)))
        else:
            parts.append("%d,%d" % (hot_num.get(hot.get(pos), 0), hub_num.get(hub.get(pos), 0)))
            parts.extend(tail)
        rows.append(",".join(parts))
    return "\n".join(rows)


def format_report(rep):
    """Render a residue_report dict as a text block for the GUI / console."""
    if not rep:
        return ""
    L = ["%s  position %d  (%s)" % (rep["uid"], rep["pos"], rep["aa"])]
    for p in rep["ptms"]:
        L.append("  PTM: %s" % p["description"])
    for s in rep["sites"]:
        L.append("  Site: %s" % s["description"])
    for d in rep["domains"]:
        L.append("  Domain: %s" % d["description"])
    for v in rep["variants"]:
        extra = []
        if v.get("clinVar"):
            extra.append("ClinVar: %s" % v["clinVar"])
        if v.get("alphaMissense") is not None:
            extra.append("AM %.2f" % v["alphaMissense"])
        if v.get("diseases"):
            extra.append("; ".join(v["diseases"]))
        L.append("  Variant %s%d%s — %s%s" % (v.get("wildType", ""), rep["pos"], v.get("mutant", ""),
                                              v["consequence"], (" [" + " | ".join(extra) + "]") if extra else ""))
    if rep["amMean"] is not None:
        L.append("  AlphaMissense (mean): %.3f" % rep["amMean"])
    if rep["amProfile"]:
        L.append("  AM top: " + ", ".join("%s%s %.2f" % (wt, mt, sc) for wt, mt, sc in rep["amProfile"]))
    if rep.get("flags"):
        L.append("  Flags: " + ", ".join("%s %s" % ("●" if on else "○", lbl) for _k, lbl, _c, on in rep["flags"]))
    nb = rep.get("nearby") or {}
    if nb.get("near"):
        L.append("  Nearby (%g Å): %s" % (nb["dist"], ", ".join(str(p) for p in nb["near"])))
    if nb.get("partners"):
        L.append("  other chains: " + ", ".join("%s:%s" % (ch, r) for ch, r, _i in nb["partners"][:40]))
    if rep.get("nearbyLigands"):
        L.append("  Nearby ligands (≤5Å): " + ", ".join(rep["nearbyLigands"]))
    pr, vr = rep.get("proxRadii") or (PROX_SETTINGS["ptm"], PROX_SETTINGS["var"])
    if rep.get("proxPtms"):
        L.append("  PTMs ≤%gÅ: " % pr + ", ".join("%d %s (%.1fÅ)" % (n["pos"], n.get("category") or "", n["dist"])
                                                   for n in rep["proxPtms"][:12]))
    if rep.get("proxVariants"):
        L.append("  Pathogenic variants ≤%gÅ: " % vr + ", ".join("%s (%.1fÅ)" % (n["label"], n["dist"])
                                                                  for n in rep["proxVariants"][:12]))
    return "\n".join(L)


def _am_color(score):
    """AlphaMissense colour: red likely-pathogenic (>=0.564), green likely-benign (<=0.34), amber
    ambiguous between — matching the extension's substitution grid."""
    if score is None:
        return "#888"
    return "#d32f2f" if score >= 0.564 else "#4caf50" if score <= 0.34 else "#f5a623"


def format_report_html(rep, expanded=False, protvar="off"):
    """Clean label/value layout for the GUI detail panel, mirroring the extension's residue card:
    a header, then colour-coded section headers with aligned 'label : value' rows (no indentation).
    `expanded` controls whether each variant's evidence (ClinVar / review / dbSNP / disease) is shown
    or folded behind a toggle link. Nearby residues are clickable (ufv:res: anchors). `protvar` is
    'off' (no section), 'loading', or a fetch_protvar() dict (renders the ProtVar predictors)."""
    if not rep:
        return ""

    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    T = ['<table cellspacing="0" cellpadding="3" '
         'style="font-family:sans-serif; font-size:12px; border-collapse:collapse;">']
    bulbs = "".join('&nbsp;<span title="%s — %s" style="color:%s; font-size:13px;">%s</span>'
                    % (esc(lbl), "flagged at this residue" if on else "computed; not flagged here",
                       col if on else "#c8c8c8", "●" if on else "○")
                    for _k, lbl, col, on in rep.get("flags") or [])
    T.append('<tr><td colspan="2" style="font-size:15px; padding-bottom:2px;">'
             '<b>%s%d</b>%s</td></tr>' % (esc(rep["aa"]), rep["pos"], bulbs))

    def hdr(title, color):
        # title is a controlled string (may contain a <span title> for a tooltip) — do NOT escape it.
        T.append('<tr><td colspan="2" style="border-top:1px solid #ccc; padding-top:5px;">'
                 '<b style="color:%s;">%s</b></td></tr>' % (color, title))

    def row(label, value, vcolor=None, bold=False):
        # label is a controlled string / pre-built HTML (variant mutations); value is caller-escaped
        v = ('<span style="color:%s;">%s</span>' % (vcolor, value)) if vcolor else value
        lab = ('<b>%s</b>' % label) if bold else ('<span style="color:#888;">%s</span>' % label)
        T.append('<tr><td style="vertical-align:top; white-space:nowrap;">%s</td>'
                 '<td style="vertical-align:top;">%s</td></tr>' % (lab, v))

    if rep["variants"]:
        hdr("Variants", "#c62828")
        has_evidence = any(v.get("clinVar") or v.get("clinVarReview") or v.get("rsIds")
                           or v.get("genomic") or v.get("diseases") for v in rep["variants"])
        for v in rep["variants"]:
            mut = "%s%d%s" % (esc(v.get("wildType", "")), rep["pos"], esc(v.get("mutant", "")))
            c = v.get("consequenceColor", "#9e9e9e")
            cons = esc(v["consequence"])
            if v.get("diseases"):       # name the disease in parentheses after the consequence
                cons += " <span style='color:#888;'>(%s)</span>" % esc("; ".join(v["diseases"]))
            row('<span style="color:%s;">%s</span>' % (c, mut), cons, c, bold=True)
            if v.get("alphaMissense") is not None:
                row("AlphaMissense", "%.2f" % v["alphaMissense"], _am_color(v["alphaMissense"]))
            if v.get("gnomad") is not None:
                af = v["gnomad"]
                rarity = "common" if af >= 0.01 else "rare" if af > 0 else "absent"
                row("gnomAD AF", "%s <span style='color:#888;'>(%s)</span>" % (_fmt_af(af), rarity),
                    "#388e3c" if af >= 0.01 else "#777")
            else:
                row("gnomAD AF", '<span title="Not observed in gnomAD — consistent with a rare variant">not in gnomAD</span>', "#999")
            if expanded:  # evidence rows, folded by default
                if v.get("clinVar"):
                    row("ClinVar", esc(v["clinVar"]))
                if v.get("clinVarReview"):
                    row("Review", esc(v["clinVarReview"]))
                if v.get("rsIds"):
                    row("dbSNP", esc(", ".join(v["rsIds"])))
                if v.get("genomic"):
                    row("Genomic", '<span style="font-family:monospace; font-size:11px;">%s</span>' % esc(v["genomic"]))
                if v.get("diseases"):
                    row("Disease", esc("; ".join(v["diseases"])))
        if has_evidence:
            link = "▾ hide evidence" if expanded else "▸ show evidence"
            T.append('<tr><td colspan="2" style="padding-top:2px;">'
                     '<a href="ufv:evi" style="color:#1565c0; text-decoration:none; font-size:11px;">%s</a>'
                     '</td></tr>' % link)

    if rep["ptms"] or rep["sites"] or rep["domains"] or rep.get("mutagenesis"):
        hdr("Features", "#1565c0")
        for p in rep["ptms"]:
            row("PTM", esc(p["description"]), p.get("color"))
        for s in rep["sites"]:
            row("Site", esc(s["description"]), SITE_COLOR)
        for d in rep["domains"]:
            row("Domain", esc(d["description"]), d.get("color"))
        for mtg in rep.get("mutagenesis", []):
            sub = "%s→%s" % (esc(mtg.get("wildType", "")), esc(", ".join(mtg.get("mutants", []))))
            val = "%s <span style='color:#888;'>%s</span>" % (sub, esc(mtg.get("effect", "")))
            row('<span title="Experimentally mutated residue (UniProt)">Mutagenesis</span>', val, MUTAG_COLOR)

    if rep.get("proximity"):
        px = rep["proximity"]
        hdr('<span title="Spatial proximity of this PTM to disease variants">PTM–variant proximity</span>', "#b87800")
        d8, d12 = px.get("d8", 8.0), px.get("d12", 12.0)
        tlabel = {1: "variant at this residue", 2: "pathogenic variant ≤ %g Å" % d8,
                  3: "pathogenic variant ≤ %g Å" % d12}
        row("Tier", tlabel.get(px.get("tier"), str(px.get("tier"))))
        if px.get("nearestVariant"):
            nd = px.get("nearestDist")
            dist = "same residue" if nd == 0 else ("%.1f Å" % nd if nd is not None else "")
            row("Nearest variant", "%s <span style='color:#888;'>(%s)</span>" % (esc(px["nearestVariant"]), dist))
        row("Within %g Å" % d8, "%d variant(s), %d pathogenic" % (px.get("nearbyCount8A", 0), px.get("pathCount8A", 0)))

    # ── Prediction: one boxed per-substitution table combining AlphaMissense (always) with ProtVar's
    # EVE / ESM1b / FoldX (when fetched). Conservation + M3D are per-position rows above it.
    pv_dict = protvar if (protvar not in ("off", "loading") and protvar) else None
    if rep.get("amProfile") or rep.get("amMean") is not None or protvar != "off":
        hdr('<span title="Per-substitution variant-effect predictors">Prediction</span>', "#6a1b9a")
        sc = (pv_dict or {}).get("score") or {}
        fx = (pv_dict or {}).get("foldx") or {}
        wt = rep["aa"]
        if sc.get("conservation") is not None:
            row('<span title="Evolutionary conservation across homologues (0 variable – 1 conserved)">Conservation</span>',
                "%.2f" % sc["conservation"])
        if sc.get("m3d"):
            md = sc["m3d"]
            feat = (" — %s" % esc(md["feature"])) if md.get("feature") else ""
            mcol = "#d32f2f" if str(md.get("prediction", "")).lower().startswith("damag") else "#388e3c"
            row('<span title="Missense3D: predicted structural consequence">M3D</span>',
                "<span style='color:%s;'>%s</span>%s" % (mcol, esc(md.get("prediction", "?")), feat))
        if protvar == "loading":
            row("EVE / CADD / ESM1b / FoldX", "<span style='color:#888;'>loading…</span>")
        am_by_mut = {m: s for (_w, m, s) in rep.get("amProfile", [])}
        show_pv = bool(pv_dict)
        eve = _protvar_by_mut(sc.get("eveRaw"), wt) if show_pv else {}
        evc = _protvar_by_mut(sc.get("eveClsRaw"), wt) if show_pv else {}
        esm = _protvar_by_mut(sc.get("esmRaw"), wt) if show_pv else {}
        ddg = fx.get("byMut") or {}
        cadd = sc.get("caddByMut") or {}
        obs = {v.get("mutant") for v in rep.get("variants", []) if v.get("mutant")}
        data = []
        for mt in [a for a in AA_ALPHA if a != wt]:
            am, e, es, d, c = am_by_mut.get(mt), eve.get(mt), esm.get(mt), ddg.get(mt), cadd.get(mt)
            if am is None and e is None and es is None and d is None and c is None:
                continue
            data.append((mt, am, e, evc.get(mt), es, d, c))
        data.sort(key=lambda r: (r[1] is None, -(r[1] if r[1] is not None else 0), -(r[2] if r[2] is not None else 0)))
        if data:
            box = "background:#f3f3f3; padding:2px 6px; text-align:center;"
            # Only show predictor columns that actually have data for this protein/position — EVE in
            # particular only covers ~3000 genes, so a column of dashes would just look like a failure.
            has_am = any(r[1] is not None for r in data)
            has_eve = any(r[2] is not None for r in data)
            has_esm = any(r[4] is not None for r in data)
            has_fx = any(r[5] is not None for r in data)
            has_cadd = any(r[6] is not None for r in data)

            def _eve(e, ec):
                return ('<span style="color:%s;">%.2f</span>' %
                        ("#d32f2f" if str(ec or "").upper() == "PATHOGENIC" else "#666", e)) if e is not None else "–"

            def _col(v, color):
                return ('<span style="color:%s;">%s</span>' % (color, v)) if color else v

            def _fx(d):          # modal.js: >=2 red, >=1 orange, <=-1 green, ~neutral uncoloured
                if d is None:
                    return "–"
                v = "%.1f" % d
                v = ("+" + v) if (d > 0 and v != "0.0") else v
                return _col(v, "#e53935" if d >= 2 else "#fb8c00" if d >= 1 else "#43a047" if d <= -1 else None)

            def _esm(x):         # modal.js: <= -10 red, <= -7.5 orange
                return _col("%.1f" % x, "#e53935" if x <= -10 else "#fb8c00" if x <= -7.5 else None)                     if x is not None else "–"

            def _cadd(c):        # modal.js: >= 25 red, >= 20 orange (phred; > 20 = top 1 %)
                return _col("%.1f" % c, "#e53935" if c >= 25 else "#fb8c00" if c >= 20 else None)                     if c is not None else "–"
            cols = []  # (header, tooltip, cell_fn(row))
            if has_am:
                cols.append(("AM", "AlphaMissense (0–1 pathogenicity)",
                             lambda r: ('<span style="color:%s;">%.2f</span>' % (_am_color(r[1]), r[1])) if r[1] is not None else "–"))
            if has_eve:
                cols.append(("EVE", "EVE evolutionary model (0–1)", lambda r: _eve(r[2], r[3])))
            if has_cadd:
                cols.append(("CADD", "CADD deleteriousness, phred-scaled (&gt;20 ≈ top 1% most deleterious, "
                                     "&gt;30 ≈ top 0.1%). Only substitutions reachable by a single nucleotide change.",
                             lambda r: _cadd(r[6])))
            if has_esm:
                cols.append(("ESM1b", "ESM-1b pathogenicity (lower / more negative = more pathogenic)",
                             lambda r: _esm(r[4])))
            if has_fx:
                cols.append(("FoldX", "FoldX ΔΔG kcal/mol; &gt;2 destabilising", lambda r: _fx(r[5])))
            head = '<td style="color:#888;">sub</td>' + "".join(
                '<td style="color:#888;" title="%s">%s</td>' % (tip, h) for h, tip, _ in cols)
            cells = ['<tr><td colspan="2"><table cellspacing="3" style="font-size:11px;"><tr>%s</tr>' % head]
            for r in data:
                mt = r[0]
                link = "https://www.ebi.ac.uk/ProtVar/query?search=%s+%s%d%s" % (rep["uid"], wt, rep["pos"], mt)
                sub = '<a href="%s" style="color:#6a1b9a; text-decoration:none;">%s%d%s</a>' % (link, esc(wt), rep["pos"], esc(mt))
                if mt in obs:
                    sub = "<b>%s •</b>" % sub
                tds = '<td style="white-space:nowrap;">%s</td>' % sub
                tds += "".join('<td style="%s">%s</td>' % (box, fn(r)) for _h, _t, fn in cols)
                cells.append("<tr>%s</tr>" % tds)
            cells.append("</table></td></tr>")
            T.append("".join(cells))

    # ── AutoPocket: one block per ProtVar (AutoSite) pocket this residue lines (modal.js renderBindingPocket)
    pk = (pv_dict or {}).get("pocket")
    if pk:
        ann_ = _CACHE.get(rep["uid"]) or {}
        geom_ = _GEOM_CACHE.get(rep.get("obj")) or []
        hdr('<span title="Predicted ligand-binding cavities this residue lines (AutoSite via ProtVar). Score ranks '
            'how pronounced / ligand-favourable the cavity is (not a druggability probability).">AutoPocket</span>'
            ' <a href="https://www.ebi.ac.uk/ProtVar/%s" style="color:#888; text-decoration:none; font-weight:normal;">'
            '⇗ ProtVar</a>' % esc(_canon_acc(rep["uid"])), "#00897b")
        pockets = pk.get("pockets") or [{"pocketId": None, "buriedness": pk.get("buriedness"), "score": pk.get("score"),
                                         "resid": pk.get("lining") or []}]
        for i, p in enumerate(pockets):
            resid = p.get("resid") or []
            stats = ["%d res" % len(resid)]
            if _finite(p.get("meanPlddt")):
                stats.append('<span title="Mean AlphaFold confidence over the pocket residues">pLDDT %d</span>'
                             % round(p["meanPlddt"]))
            if _finite(p.get("radGyration")):
                stats.append('<span title="Radius of gyration (Å): smaller = more compact cavity">Rg %.1f</span>'
                             % p["radGyration"])
            if _finite(p.get("energyPerVol")):
                stats.append('<span title="Interaction-affinity energy per pocket volume">E/V %.2f</span>' % p["energyPerVol"])
            if _finite(p.get("score")):
                stats.append('<span title="AutoSite pocket score">score %d</span>' % round(p["score"]))
            if _finite(p.get("buriedness")):
                stats.append('<span title="0 = open surface, 1 = fully buried">buried %.2f</span>' % p["buriedness"])
            pid = p.get("pocketId") if p.get("pocketId") is not None else i + 1
            links = ""
            if resid:
                links = (' <a href="ufv:pocket:%d" style="color:#fb8c00; text-decoration:none;">show ▸</a>' % i)
                if len(resid) >= 2 and rep.get("obj") and _structure_motif_ok(rep["obj"]):
                    links += (' <a href="ufv:motif:%d" style="color:#1565c0; text-decoration:none;" title="RCSB Structure '
                              'Motif Search on this pocket">Find Similar Motifs ↗</a>' % i)
            row("Pocket %s" % esc(pid), "<span style='color:#00897b;'>%s</span>%s" % (" · ".join(stats), links))
            if resid:
                ev = pocket_evidence(resid, geom_, ann_.get("sites") or [], ann_.get("sequence") or "")
                bits = []
                if ev["composition"]:
                    bits.append(spider_img_html(ev["composition"], 120))
                if ev["siteHits"]:
                    bits.append("<span style='color:%s;'>sites: %s</span>" % (SITE_COLOR, esc("; ".join(
                        "%s (%d)" % (h.get("description") or "site", h["position"]) for h in ev["siteHits"][:4]))))
                if bits:
                    row("", "<span style='color:#666; font-size:11px;'>%s</span>" % "<br>".join(bits))
        T.append('<tr><td colspan="2"><a href="ufv:allpockets" style="color:#00897b; text-decoration:none; '
                 'font-size:11px;">All pockets of this protein ▸</a></td></tr>')

    # ── Experimental binding (PDBe-KB ligand sites + PPI interfaces across every PDB entry)
    kb = rep.get("pdbekb")
    if kb and kb != "pending":
        lt, nl, it, ni = _kb_residue(kb, rep["pos"])
        if lt or it:
            hdr('<span title="Aggregated from every experimental PDB structure of this protein (PDBe-KB)">'
                'Experimental binding</span> <a href="https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/%s" '
                'style="color:#888; text-decoration:none; font-weight:normal;">⇗ PDBe-KB</a>'
                % esc(_canon_acc(rep["uid"])), "#2e7d32")
            if lt:
                row("Ligands", ", ".join(
                    '<a href="%s" title="%s" style="color:#2e7d32; text-decoration:none;">%s</a>' % (
                        ("https://www.rcsb.org/structure/%s" % str(l["pdb"]).upper()) if l.get("pdb")
                        else ("https://www.rcsb.org/ligand/%s" % esc(l["code"])),
                        esc("%s — bound in %d PDB structure%s%s" % (l.get("name") or l["code"], l["pdbCount"],
                                                                   "" if l["pdbCount"] == 1 else "s",
                                                                   (" (e.g. %s)" % str(l["pdb"]).upper()) if l.get("pdb") else "")),
                        esc(l["code"])) for l in lt) + ((" +%d" % (nl - 8)) if nl > 8 else ""))
            if it:
                row("PPI interface", "; ".join(
                    '<a href="%s" title="%s" style="color:#1565c0; text-decoration:none;">%s</a>' % (
                        ("https://www.rcsb.org/structure/%s" % str(p["pdb"]).upper()) if p.get("pdb")
                        else "https://www.uniprot.org/uniprotkb/%s/entry" % esc(p.get("accession") or ""),
                        esc("Interface seen in %d PDB complex%s" % (p["pdbCount"], "" if p["pdbCount"] == 1 else "es")),
                        esc(p.get("name") or p.get("accession") or "?")) for p in it) + ((" +%d" % (ni - 4)) if ni > 4 else ""))

    # ── Nearby (any-atom distance, slider 3-15 Å: ufv_nearby) + other chains of the complex
    nb = rep.get("nearby")
    if nb and (nb.get("near") or nb.get("partners")):
        seq = (_CACHE.get(rep["uid"]) or {}).get("sequence") or ""
        vc = rep.get("varColor") or {}
        chips = ", ".join('<a href="ufv:res:%d" style="color:%s; text-decoration:none;">%s%d</a>'
                          % (p, vc.get(p) or "#37474f", esc(seq[p - 1] if 0 < p <= len(seq) else ""), p)
                          for p in nb["near"]) or "—"
        hdr('<span title="Residues with any atom within the Nearby distance (ufv_nearby, 3-15 Å); coloured by '
            'variant consequence">Nearby (%g Å)</span>' % nb["dist"], "#455a64")
        row("Residues", chips)
        if nb.get("partners"):
            pc = []
            for ch, resi, info in nb["partners"][:40]:
                col, tip = "#607d8b", "chain %s %s" % (ch, resi)
                if info:
                    vs = info.get("variants") or []
                    if vs:
                        col = vs[0].get("diseaseColor") or vs[0].get("consequenceColor") or col
                    tip = "%s %s — %s" % (info.get("accession") or "", info.get("uniPos"),
                                          "; ".join(sorted({d for v in vs for d in (v.get("diseases") or [])})) or
                                          ("%d variant(s)" % len(vs)))
                pc.append('<span title="%s" style="color:%s;">%s:%s</span>' % (esc(tip), col, esc(ch), esc(resi)))
            row("other chains", ", ".join(pc) + ((" +%d" % (len(nb["partners"]) - 40)) if len(nb["partners"]) > 40 else ""))

    if rep.get("nearbyLigands"):
        hdr("Nearby ligands (≤5 Å)", "#6d4c41")
        row("Ligands", esc(", ".join(rep["nearbyLigands"])))

    if rep.get("proxPtms") or rep.get("proxVariants"):
        # residue-centric PTM–variant proximity (analysis.js computeResidueProximity); radii: ufv_prox
        pr, vr = rep.get("proxRadii") or (PROX_SETTINGS["ptm"], PROX_SETTINGS["var"])

        def chip(n, text, col, tip):
            # clickable: jumps the report to that residue (ufv:res anchor)
            return ('<a href="ufv:res:%d" title="%s" style="color:%s; text-decoration:none;">'
                    '%s&nbsp;(%.1f)</a>' % (n["pos"], esc(tip), col or "#37474f", esc(text), n["dist"]))
        hdr('<span title="PTMs and pathogenic variants near this residue (Cα distance) — radii: ufv_prox">'
            'PTM–variant proximity</span>', "#37474f")
        if rep.get("proxPtms"):
            row("PTMs ≤ %g Å" % pr, ", ".join(chip(n, str(n["pos"]), n.get("color"),
                                                    "%s — %s" % (n.get("category") or "", n.get("description") or ""))
                                               for n in rep["proxPtms"][:24]))
        if rep.get("proxVariants"):
            row("Pathogenic ≤ %g Å" % vr, ", ".join(chip(n, n["label"], n.get("color"),
                                                         "%s%s" % (n.get("consequence") or "",
                                                                   (" — " + "; ".join(n["diseases"])) if n.get("diseases") else ""))
                                                    for n in rep["proxVariants"][:24]))

    T.append("</table>")
    return "".join(T)


def format_ligand_html(resn, info, instances=1, copies=None, cur=0, aff=None, evidence=None, clash=None,
                       motif=False, experimental=None, near=None, nearby_dist=None, sequence=""):
    """Ligand panel (modal.js renderLigandPanel): copy stepper, AlphaFill transplant rows (donor identity,
    clash band), chemistry (name, formula, mol. weight, H-bond D/A, SMILES, InChIKey), pocket evidence
    (closest contact band, mean pocket pLDDT / B-factor, functional sites, residue-class composition,
    strucmotif link), external links (PubChem · DrugBank · 2D / 3D similarity · ChEMBL) and the in-structure
    Tanimoto grid. When the component occurs multiple times, `copies` = [(chain, resi), ...] and `cur` is
    the focused copy; a ◂ i/N ▸ stepper (ufv:lig anchors) circulates them."""
    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    info = info or {}
    n = instances if not copies else len(copies)
    T = ['<table cellspacing="0" cellpadding="3" style="font-family:sans-serif; font-size:12px;">']
    title = "%s%s" % (esc(resn), (" ×%d" % n) if n > 1 else "")
    T.append('<tr><td colspan="2" style="font-size:14px;"><b style="color:#2e7d32;">%s</b> '
             '<span style="color:#888;">(CCD)</span></td></tr>' % title)
    if copies and len(copies) > 1:
        ch, resi = copies[cur % len(copies)]
        T.append('<tr><td colspan="2" style="padding:3px 0;">'
                 '<span style="background:#e8f0fe; padding:2px 6px; border-radius:3px;">'
                 '<a href="ufv:lig:prev" style="color:#1565c0; text-decoration:none; font-size:14px;">◂</a>'
                 '&nbsp;<b>showing copy %d of %d</b> &mdash; <b style="color:#2e7d32;">chain %s · %s</b>&nbsp;'
                 '<a href="ufv:lig:next" style="color:#1565c0; text-decoration:none; font-size:14px;">▸</a>'
                 '</span> <span style="color:#888; font-size:10px;">(highlighted in 3D)</span>'
                 '</td></tr>' % (cur % len(copies) + 1, len(copies), esc(ch or "?"), esc(resi)))

    def row(label, value):
        T.append('<tr><td style="color:#888; vertical-align:top; white-space:nowrap;">%s</td>'
                 '<td style="vertical-align:top;">%s</td></tr>' % (esc(label), value))

    def sec(title, color="#6d4c41"):
        T.append('<tr><td colspan="2" style="border-top:1px solid #ccc; padding-top:5px;">'
                 '<b style="color:%s;">%s</b></td></tr>' % (color, title))
    if near is not None:
        seq = sequence or ""
        row("Nearby%s" % ((" (%g Å)" % nearby_dist) if nearby_dist else ""),
            ", ".join('<a href="ufv:res:%d" style="color:#37474f; text-decoration:none;">%s%d</a>'
                      % (p, esc(seq[p - 1] if 0 < p <= len(seq) else ""), p) for p in near) or "—")
    if aff:                                  # AlphaFill transplant evidence (information only, nothing filtered)
        if aff.get("identity") is not None:
            row("Donor identity", "%d%%%s" % (round(aff["identity"] * 100),
                                              (" (from %s)" % esc(aff["donorPdb"])) if aff.get("donorPdb") else ""))
        if aff.get("clash") is not None:
            band = aff_clash_band(aff["clash"])
            row("Clash", '<span style="color:%s;">%.2f (%s)</span>' % (
                {"low": "#388e3c", "moderate": "#f57c00", "high": "#d32f2f"}[band], aff["clash"], band))
    if info.get("name"):
        row("Name", esc(info["name"]))
    elif info:
        row("Name", "%s (no chemical record found)" % esc(resn))
    if info.get("formula"):
        row("Formula", esc(info["formula"]))
    if info.get("weight") is not None:
        row("Mol. weight", "%.1f Da" % info["weight"])
    if info.get("hbondDonors") is not None or info.get("hbondAcceptors") is not None:
        row("H-bond D/A", "%s / %s" % (info["hbondDonors"] if info.get("hbondDonors") is not None else "?",
                                       info["hbondAcceptors"] if info.get("hbondAcceptors") is not None else "?"))
    if info.get("smiles"):
        row("SMILES", '<span style="font-family:monospace; font-size:11px;">%s</span>' % esc(info["smiles"]))
    if info.get("inchikey"):
        row("InChIKey", '<span style="font-family:monospace; font-size:11px;">%s</span>' % esc(info["inchikey"]))
    ev = evidence or {}
    if clash or ev.get("meanB") is not None or ev.get("siteHits") or ev.get("composition"):
        sec("Pocket", "#00897b")
        if clash:
            col = {"severe": "#d32f2f", "suspicious": "#f57c00", "normal": "#388e3c"}[clash["band"]]
            row("Closest contact", '<span style="color:%s;">%.2f Å (%s)</span>' % (col, clash["minDist"], clash["band"]))
        if ev.get("meanB") is not None:
            tip = "mean over %d pocket residue%s" % (ev["count"], "" if ev["count"] == 1 else "s")
            if experimental:
                row("Mean pocket B-factor", '<span title="%s">%.1f</span>' % (tip, ev["meanB"]))
            else:
                band = plddt_band(ev["meanB"])
                row("Mean pocket pLDDT", '<span title="%s" style="color:%s;">%.1f (%s)</span>' % (
                    tip, PLDDT_BAND_COLORS[band], ev["meanB"], band))
        if ev.get("siteHits"):
            row("Functional sites", '<span style="color:%s;">%s</span>' % (SITE_COLOR, esc("; ".join(
                "%s (%d)" % (h.get("description") or "site", h["position"]) for h in ev["siteHits"]))))
        if ev.get("composition"):
            row("Pocket character", spider_img_html(ev["composition"]))
        if motif:
            row("Pocket similarity", '<a href="ufv:motif:lig" style="color:#1565c0; text-decoration:none;" '
                                     'title="RCSB Structure Motif Search on this pocket">Find Similar Motifs ↗</a>')
    links = info.get("links") or (ligand_links(info, info.get("cid")) if info else [])
    if links:
        sec("Chemical info")
        T.append('<tr><td colspan="2">%s</td></tr>' % " · ".join(
            '<a href="%s" style="color:#1565c0; text-decoration:none;">%s</a>' % (esc(u), esc(lb)) for lb, u in links))
    sims = info.get("similar") or []
    if sims:
        sec("Ligand similarity within structure (Tanimoto)")
        cells = ['<tr><td colspan="2"><table cellspacing="2"><tr>']
        for i, (d, s) in enumerate(sims[:LIGAND_SIM_CAP]):
            if i and i % 3 == 0:
                cells.append('</tr><tr>')
            cells.append('<td style="background:#f3f3f3; padding:2px 6px;">'
                         '<a href="ufv:lig:show:%s" title="CACTVS substructure-fingerprint Tanimoto %.2f (click to focus)" '
                         'style="color:#2e7d32; text-decoration:none;"><b>%s</b></a>&nbsp;%.2f</td>' % (esc(d), s, esc(d), s))
        cells.append('</tr></table></td></tr>')
        T.append("".join(cells))
    T.append("</table>")
    return "".join(T)


_CONS_RANK = {"pathogenic": 4, "deleterious": 3, "uncertain": 2, "benign": 1}


def _annotation_color_map(ann):
    """uniprot pos -> annotation colour (variant > PTM > site). For variants, the MOST SEVERE
    consequence at a position wins (so a pathogenic variant isn't overwritten by a co-located
    uncertain one — that's what made disease colours look absent on the focus sticks)."""
    m = {}
    for s in ann.get("sites", []):
        m[s["position"]] = s["color"]
    for p in ann.get("ptms", []):
        m[p["position"]] = p["color"]
    best = {}
    for v in ann.get("variants", []):
        pos = v["position"]
        r = _CONS_RANK.get(_cons_token(v["consequence"]), 0)
        if best.get(pos, -1) < r:
            best[pos] = r
            m[pos] = v["consequenceColor"]
    return m


def _fetch_model_text(struct):
    """Download a listed structure -> (text, ext) or (None, None). PDB entries too large for the legacy
    PDB format (most big cryo-EM assemblies) have no pdbXXXX.ent; fall back to the entry's mmCIF."""
    data = _get_structure_text(struct["url"])
    if data:
        return data, ("cif" if struct.get("fmt") == "mmcif" else "pdb")
    if struct.get("pdbId"):
        data = _get_structure_text(PDBe_CIF.format(struct["pdbId"].lower()))
        if data:
            return data, "cif"
    return None, None


def _write_temp(data, ext):
    import tempfile
    fd, path = tempfile.mkstemp(prefix="ufv_", suffix="." + ext, dir=_ufv_tmpdir())   # unique: parallel downloads
    with os.fdopen(fd, "w") as fh:
        fh.write(data)
    return path


_HIDE_POINT_ALIASES = {"ptms": "ptm", "ptm": "ptm", "variants": "var", "variant": "var", "var": "var",
                       "sites": "site", "site": "site", "mutagenesis": "mutag", "mutag": "mutag",
                       "filtered": "filt", "filt": "filt", "partners": "prt", "partner": "prt", "prt": "prt",
                       "binding": "kb", "pdbekb": "kb", "kb": "kb"}


def _structure_is_predicted(obj):
    """True for AlphaFold/computed models (B-factor column is pLDDT), False for experimental."""
    return (_STATE.get("sources") or {}).get(obj, "PDB") != "PDB"


# ----------------------------------------------------------------------------------------------
# Trajectory: multi-state objects (MD trajectories, NMR ensembles)
# ----------------------------------------------------------------------------------------------
_TRAJ_GEOM_STATE = {}   # obj -> state the cached Cα geometry / analyses belong to
_TRAJ_STASH = {}        # obj -> {state: (geom, {analysis key: value})}: other states' caches
_TRAJ_COORDS = {}       # obj -> {"sig": (nstates, natoms), "sets": {(sel, atoms): (rows, X(nstates, n, 3))}}
_TRAJ_RESULTS = {}      # obj -> {kind: result}: last run of each analysis (CSV, B-factor, panel)


def _traj_forget(obj):
    """Drop every trajectory cache of `obj` (called from _invalidate_obj)."""
    for d in (_TRAJ_GEOM_STATE, _TRAJ_STASH, _TRAJ_COORDS, _TRAJ_RESULTS):
        d.pop(obj, None)


def _traj_take(X, frames, cols=None, f32=False):
    """X[frames-1][:, cols] as float64 (float32 with f32: big contact sets) — a copy, safe to hand to a
    worker thread."""
    np = _numpy()
    fi = np.asarray(frames, dtype=int) - 1
    out = X[fi] if cols is None else X[np.ix_(fi, np.asarray(cols, dtype=int))]
    return out.astype(np.float32 if f32 else np.float64, copy=not f32 or out.base is X)


def _kabsch_batch(P, Q):
    """Optimal rotations of every frame onto a reference: P (F, N, 3), Q (N, 3) -> (R (F,3,3), pc (F,1,3),
    qc (1,3)) with (P - pc) @ R + qc ~= Q (one batched 3x3 SVD per frame)."""
    np = _numpy()
    pc = P.mean(1, keepdims=True)
    qc = Q.mean(0, keepdims=True)
    H = np.einsum("fni,nj->fij", P - pc, Q - qc)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    d[d == 0] = 1.0
    U[:, :, 2] *= d[:, None]
    return U @ Vt, pc, qc


def _fit(P, Q):
    R, pc, qc = _kabsch_batch(P, Q)
    return (P - pc) @ R + qc, (R, pc, qc)


def _apply_fit(T, X):
    R, pc, qc = T
    # `+ qc` in place: same values as `(X - pc) @ R + qc`, one (F, N, 3) temporary fewer (~90 MB for a
    # 500-frame, 15 k-atom contact set)
    if X.dtype == "float32":                    # big atom sets: stay in float32 (half the memory traffic)
        f = X.dtype.type
        Y = (X - pc.astype(f)) @ R.astype(f)
        Y += qc.astype(f)
        return Y
    Y = (X - pc) @ R
    Y += qc
    return Y


def _traj_events(x, min_scale=0.25):
    """Opening / closing runs of a pocket size series: frames above / below median ± 2 x robust sd
    (1.4826 MAD, at least `min_scale` Å) -> (state per frame +1/0/-1, [(type, i0, i1, peak)])."""
    np = _numpy()
    m = float(np.median(x))
    s = max(min_scale, 1.4826 * float(np.median(np.abs(x - m))))
    st = np.where(x > m + 2 * s, 1, np.where(x < m - 2 * s, -1, 0))
    ev, i = [], 0
    while i < len(st):
        if st[i]:
            j = i
            while j + 1 < len(st) and st[j + 1] == st[i]:
                j += 1
            seg = x[i:j + 1]
            ev.append(("open" if st[i] > 0 else "closed", i, j, float(seg.max() if st[i] > 0 else seg.min())))
            i = j + 1
        else:
            i += 1
    return st, ev


def _max_disp(Y, r, chunk=64):
    """Per atom: the largest distance from frame `r` over all frames of Y (F, N, 3) — the same values as
    np.sqrt(((Y - Y[r]) ** 2).sum(-1)).max(0), computed `chunk` frames at a time so the (F, N, 3)
    temporary (~90 MB for 500 frames x 15 k atoms) is never built."""
    np = _numpy()
    out = None
    for f0 in range(0, len(Y), chunk):
        d = np.sqrt(((Y[f0:f0 + chunk] - Y[r]) ** 2).sum(-1)).max(0)
        out = d if out is None else np.maximum(out, d)
    return out


def _contact_pairs(YA, YB, ra, rb, r, cutoff, chA=None, chB=None, upper=False):
    """Candidate atom pairs (ia, ib) that CAN come within `cutoff` in some frame — exact: with dA / dB the
    largest displacement of each atom from the reference frame (after fitting), a pair whose reference
    distance exceeds cutoff + dA + dB never gets closer. Residues are screened first as bounding spheres
    (centre, radius + displacement), then the atoms of the surviving residue pairs. ra / rb: residue id
    per atom; chA / chB: inter-chain pairs only; upper: A is B, each residue pair once."""
    np = _numpy()
    dA = _max_disp(YA, r).astype(np.float64) + 1e-4     # + float32 rounding slack
    dB = dA if YB is YA else _max_disp(YB, r).astype(np.float64) + 1e-4
    A0, B0 = YA[r].astype(np.float64), YB[r].astype(np.float64)

    def spheres(X0, d, res):
        _u, inv = np.unique(res, return_inverse=True)
        cnt = np.bincount(inv)
        cen = np.zeros((len(cnt), 3))
        np.add.at(cen, inv, X0)
        cen /= cnt[:, None]
        rad = np.zeros(len(cnt))
        np.maximum.at(rad, inv, np.sqrt(((X0 - cen[inv]) ** 2).sum(1)) + d)
        return cen, rad, np.argsort(inv, kind="stable"), np.concatenate([[0], np.cumsum(cnt)[:-1]]), cnt
    cA, radA, oA, sA, nA = spheres(A0, dA, ra)
    cB, radB, oB, sB, nB = (cA, radA, oA, sA, nA) if (YB is YA and ra is rb) else spheres(B0, dB, rb)
    ok = np.sqrt(((cA[:, None, :] - cB[None, :, :]) ** 2).sum(-1)) - radA[:, None] - radB[None, :] <= cutoff
    if chA is not None:
        ok &= chA[oA[sA]][:, None] != chB[oB[sB]][None, :]
    if upper:
        ok &= np.arange(len(cA))[:, None] < np.arange(len(cB))[None, :]
    pa, pb = np.nonzero(ok)
    if not len(pa):
        return np.zeros(0, int), np.zeros(0, int)
    m = nA[pa] * nB[pb]                      # expand residue pairs to their atom pairs (no Python loop)
    pid = np.repeat(np.arange(len(pa)), m)
    off = np.arange(int(m.sum())) - np.repeat(np.cumsum(m) - m, m)
    nbp = nB[pb][pid]
    ia = oA[sA[pa][pid] + off // nbp]
    ib = oB[sB[pb][pid] + off % nbp]
    keep = np.sqrt(((A0[ia] - B0[ib]) ** 2).sum(1)) - dA[ia] - dB[ib] <= cutoff
    return ia[keep], ib[keep]


def _contact_series(YA, YB, ia, ib, cutoff, K=10):
    """(F, P) bool: pair p within `cutoff` in frame f. Exact and windowed: in each window of K frames only
    the pairs whose distance at the window's first frame, minus both atoms' largest displacement within
    the window, is <= cutoff are evaluated frame by frame (a fraction of the candidates)."""
    np = _numpy()
    F, P = YA.shape[0], len(ia)
    C = np.zeros((F, P), bool)
    if not P:
        return C
    c2 = cutoff * cutoff
    # Only the atoms that occur in a candidate pair matter: per window their coordinates are gathered and
    # the displacement bound is computed for them alone (not for every atom of the selection). Same
    # element-wise arithmetic on the same values -> identical booleans.
    same = YB is YA
    if same:
        uA, inv = np.unique(np.concatenate([ia, ib]), return_inverse=True)
        ja, jb = inv[:P], inv[P:]
    else:
        uA, ja = np.unique(ia, return_inverse=True)
        uB, jb = np.unique(ib, return_inverse=True)
    for w0 in range(0, F, K):
        w1 = min(F, w0 + K)
        ZA = YA[w0:w1][:, uA]
        ZB = ZA if same else YB[w0:w1][:, uB]
        dA = np.sqrt(((ZA - ZA[0]) ** 2).sum(-1)).max(0)
        dB = dA if same else np.sqrt(((ZB - ZB[0]) ** 2).sum(-1)).max(0)
        dg = np.sqrt(((ZA[0, ja] - ZB[0, jb]) ** 2).sum(-1))
        live = np.nonzero(dg - dA[ja] - dB[jb] <= cutoff + 1e-3)[0]
        if len(live):
            C[w0:w1, live] = ((ZA[:, ja[live]] - ZB[:, jb[live]]) ** 2).sum(-1) <= c2
    return C


def _reduce_any(C, keys):
    """(F, P) bool grouped by integer keys (P,) -> (unique keys, (F, U) bool 'any pair of the group')."""
    np = _numpy()
    if not len(keys):
        return np.zeros(0, int), np.zeros((C.shape[0], 0), bool)
    order = np.argsort(keys, kind="stable")
    ks = keys[order]
    uk, start = np.unique(ks, return_index=True)
    R = np.maximum.reduceat(C[:, order].view(np.uint8), start, axis=1).astype(bool)
    return uk, R


def _switches(B):
    """Number of on/off transitions per column of a (F, n) bool matrix."""
    np = _numpy()
    return (np.diff(B.astype(np.int8), axis=0) != 0).sum(0) if B.shape[0] > 1 else np.zeros(B.shape[1], int)


def _traj_clusters(xyz, cols, link=8.0, min_size=4):
    """Single-linkage clusters (Cα <= `link` Å in the reference frame) of the residue columns `cols`,
    largest first, >= min_size members: the constraint-pocket class is spread over the protein, a pocket
    is one spatial cluster of it."""
    np = _numpy()
    cols = list(dict.fromkeys(cols))
    if not cols:
        return []
    P = xyz[cols]
    adj = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)) <= link
    lab = -np.ones(len(cols), int)
    for s in range(len(cols)):
        if lab[s] >= 0:
            continue
        lab[s] = s
        stack = [s]
        while stack:
            i = stack.pop()
            for j in np.nonzero(adj[i] & (lab < 0))[0]:
                lab[j] = s
                stack.append(j)
    groups = {}
    for i, l in enumerate(lab):
        groups.setdefault(l, []).append(cols[i])
    return sorted([g for g in groups.values() if len(g) >= min_size], key=lambda g: -len(g))


def _traj_calc(p):
    """Pure-numpy half: returns the result dict (series, table, per-residue values). Thread-safe."""
    t = time.time()
    k = p["kind"]
    if k == "annot":
        r = _calc_annot(p)
    else:
        r = {"rmsf": _calc_rmsf, "dist": _calc_dist, "pockets": _calc_pockets,
             "ligand": _calc_contacts, "partner": _calc_contacts, "nearby": _calc_contacts}[k](p)
    r.update({"kind": k, "obj": p["obj"], "uid": p["uid"], "frames": p["frames"], "ref": p["ref"],
              "prep_ms": p.get("prep_ms", 0.0), "calc_ms": (time.time() - t) * 1000})
    return r


def _ref_index(p):
    fr = p["frames"]
    return fr.index(p["ref"]) if p["ref"] in fr else None


def _calc_rmsf(p):
    np = _numpy()
    X, Q = p["CA"], p["CAref"]
    Y, _T = _fit(X, Q)
    rmsd = np.sqrt(((Y - Q) ** 2).sum(-1).mean(1))
    for _i in range(2):                          # refit to the average structure (standard RMSF)
        Y, _T = _fit(Y, Y.mean(0))
    M = Y.mean(0)
    rmsf = np.sqrt(((Y - M) ** 2).sum(-1).mean(0))
    rg = np.sqrt(((X - X.mean(1, keepdims=True)) ** 2).sum(-1).mean(1))
    uni = [r[4] for r in p["res"]]
    rows = []
    allv = {i: float(v) for i, v in enumerate(rmsf)}
    for g, ps in sorted(p.get("groups", {}).items()):
        inn = [allv[i] for i, u in enumerate(uni) if u is not None and u in ps]
        out = [allv[i] for i, u in enumerate(uni) if u is not None and u not in ps]
        if inn and out:
            a, b = sum(inn) / len(inn), sum(out) / len(out)
            rows.append([g, len(inn), "%.2f" % a, "%.2f" % b, "%.2f" % (a / b if b else 0.0)])
    order = np.argsort(-rmsf)[:15]
    top = ",".join("%s" % (uni[i] if uni[i] is not None else "%s/%s" % (p["res"][i][0], p["res"][i][1])) for i in order)
    return {"series": [("RMSD (Å)", rmsd), ("Rg (Å)", rg)], "rmsd": rmsd, "rmsf": rmsf, "res": p["res"],
            "table": (["Group", "n", "RMSF", "rest", "ratio"], rows), "rowmeta": [{"group": r[0]} for r in rows],
            "summary": "RMSD %.2f ± %.2f Å (max %.2f); mean RMSF %.2f Å; most mobile: %s"
                       % (rmsd.mean(), rmsd.std(), rmsd.max(), rmsf.mean(), top)}


def _calc_dist(p):
    np = _numpy()
    series, rows, meta = [], [], []
    thr = p["threshold"]
    for i, pr in enumerate(p["pairs"]):
        A, B = pr["A"], pr["B"]
        disp = None
        if pr["mode"] == "ca":
            d = np.sqrt(((A[:, 0] - B[:, 0]) ** 2).sum(-1))
            disp = (pr["idxA"][0], pr["idxB"][0])
        elif pr["mode"] == "cen":
            d = np.sqrt(((A.mean(1) - B.mean(1)) ** 2).sum(-1))
        else:
            D = np.sqrt(((A[:, :, None, :] - B[:, None, :, :]) ** 2).sum(-1))       # (F, na, nb)
            d = D.reshape(D.shape[0], -1).min(1)
            a, b = np.unravel_index(int(D.mean(0).argmin()), D.shape[1:])
            disp = (pr["idxA"][a], pr["idxB"][b])
        pr["disp"] = disp
        rng = float(d.max() - d.min())
        flag = rng > thr
        series.append((pr["label"], d))
        rows.append([pr["label"], pr["mode"], "%.2f" % d.mean(), "%.2f" % d.std(), "%.2f" % d.min(),
                     "%.2f" % d.max(), "%+.2f" % (d[-1] - d[0]), "⚑" if flag else ""])
        meta.append({"pair": i, "flag": flag})
    nflag = sum(m["flag"] for m in meta)
    return {"series": series, "pairs": [{k: v for k, v in pr.items() if k not in ("A", "B")} for pr in p["pairs"]],
            "dist": np.array([s[1] for s in series]), "threshold": thr,
            "table": (["Pair", "mode", "mean", "sd", "min", "max", "Δ", "⚑"], rows), "rowmeta": meta,
            "summary": "%d pairs; %d change by > %.1f Å" % (len(rows), nflag, thr)}


def _calc_pockets(p):
    np = _numpy()
    X, r = p["CA"], _ref_index(p)
    Y, T = _fit(X, p["CAref"])
    Yr = Y[r] if r is not None else p["CAref"]
    disp = np.sqrt(((Y - Yr) ** 2).sum(-1)).max(0)
    series, rows, meta, persist = [], [], [], {}
    for pk in p["pockets"]:
        L = np.array(pk["cols"])
        Ls = set(L.tolist())
        P = Y[:, L]
        cen = P.mean(1)
        rg = np.sqrt(((P - cen[:, None]) ** 2).sum(-1).mean(1))
        D = np.sqrt(((P[:, :, None] - P[:, None, :]) ** 2).sum(-1))
        nL = len(L)
        mpd = D.sum((1, 2)) / (nL * (nL - 1))
        # buriedness: Cα within 10 Å of each lining Cα (mean over the lining); exact candidate prefilter
        dref = np.sqrt(((Yr[L][:, None] - Yr[None]) ** 2).sum(-1))
        cand = np.nonzero(((dref - disp[L][:, None]) - disp[None]).min(0) <= 10.0)[0]
        Dn = np.sqrt(((P[:, :, None] - Y[:, cand][:, None]) ** 2).sum(-1))
        bur = ((Dn <= 10.0).sum(2) - 1).mean(1)
        # lining persistence: residues within R (max lining Cα–centroid distance in the reference) of the centroid
        cref = Yr[L].mean(0)
        R = float(np.sqrt(((Yr[L] - cref) ** 2).sum(-1)).max())
        near = np.nonzero(np.sqrt(((Yr - cref) ** 2).sum(-1)) - disp <= R)[0]
        near = np.union1d(near, L)
        inside = np.sqrt(((Y[:, near] - cen[:, None]) ** 2).sum(-1)) <= R
        frac = inside.mean(0)
        st, ev = _traj_events(mpd)
        lab = pk["label"]
        series += [("%s Rg" % lab, rg), ("%s size" % lab, mpd), ("%s buried" % lab, bur)]
        pers = {int(c): float(f) for c, f in zip(near, frac)}
        persist[lab] = pers
        lost = sorted((f, c) for c, f in pers.items() if c in Ls and f < 0.5)
        gained = sorted((-f, c) for c, f in pers.items() if c not in Ls and f >= 0.3)

        def nm(c):
            u = p["res"][c][4]
            return str(u) if u is not None else "%s/%s" % (p["res"][c][0], p["res"][c][1])
        rows.append([lab, nL, "%.2f±%.2f" % (rg.mean(), rg.std()), "%.1f–%.1f" % (mpd.min(), mpd.max()),
                     "%.1f" % bur.mean(), "%d/%d" % (sum(e[0] == "open" for e in ev), sum(e[0] == "closed" for e in ev)),
                     "%.0f%%" % (100 * (st > 0).mean()),
                     " ".join(nm(c) for _f, c in lost[:4]) or "-", " ".join(nm(c) for _f, c in gained[:4]) or "-"])
        meta.append({"pocket": lab, "cols": L.tolist(), "events": ev})
    return {"series": series, "persist": persist, "res": p["res"], "pockets": meta,
            "table": (["Pocket", "n", "Rg", "size", "buried", "open/closed", "open%", "lost", "gained"], rows),
            "rowmeta": meta, "summary": "%d pockets; size = mean pairwise lining Cα distance; buried = Cα within 10 Å"
                                       % len(meta)}


def _calc_contacts(p):
    """Contact occupancy (heavy atoms within cutoff) per residue / residue pair across the window."""
    np = _numpy()
    kind, c = p["kind"], p["cutoff"]
    _Y, T = _fit(p["CA"], p["CAref"])
    YA = _apply_fit(T, p["A"])
    YB = YA if p["same"] else _apply_fit(T, p["B"])
    r = _ref_index(p)
    if r is None:                                    # reference outside the window: prepend it
        YA = np.concatenate([p["Aref"][None], YA])
        YB = YA if p["same"] else np.concatenate([p["Bref"][None], YB])
        r, off = 0, 1
    else:
        off = 0
    rowsA, rowsB = p["rowsA"], p["rowsB"]
    resA_key = [(rw[0], rw[2]) for rw in rowsA]
    resB_key = [(rw[0], rw[2]) for rw in rowsB]
    allkeys = list(dict.fromkeys(resA_key + resB_key))
    kid = {k: i for i, k in enumerate(allkeys)}
    ra = np.array([kid[k] for k in resA_key])
    rb = ra if p["same"] else np.array([kid[k] for k in resB_key])
    if kind == "partner":
        chA = np.array([rw[0] for rw in rowsA])
        chB = chA if p["same"] else np.array([rw[0] for rw in rowsB])
        ia, ib = _contact_pairs(YA, YB, ra, rb, r, c, chA, chB, upper=p["same"])
    else:
        ia, ib = _contact_pairs(YA, YB, ra, rb, r, c)
    ncand = len(ia)
    nk = len(allkeys)
    if kind == "partner":
        pk = np.minimum(ra[ia], rb[ib]) * nk + np.maximum(ra[ia], rb[ib])
    else:
        pk = ra[ia] * nk + rb[ib]
    o = np.argsort(pk, kind="stable")                # group pairs by residue pair: reduceat needs no copy
    ia, ib, pk = ia[o], ib[o], pk[o]
    C = _contact_series(YA, YB, ia, ib, c)[off:]
    ever = C.any(0)                                  # keep only pairs that ever touch
    C, pk = C[:, ever], pk[ever]
    F = C.shape[0]
    if len(pk):                                      # residue-pair occupancy
        upk, start = np.unique(pk, return_index=True)
        Cp = np.maximum.reduceat(C.view(np.uint8), start, axis=1).astype(bool)
    else:
        upk, Cp = np.zeros(0, int), np.zeros((F, 0), bool)
    occ_pair = Cp.mean(0) if F else np.zeros(len(upk))
    pa, pb = upk // nk, upk % nk
    if kind == "partner":                            # per-residue occupancy (both sides for inter-chain)
        ures, Cr = _reduce_any(np.concatenate([Cp, Cp], axis=1), np.concatenate([pa, pb]))
    else:
        ures, Cr = _reduce_any(Cp, pa)
    occ = Cr.mean(0) if F else np.zeros(len(ures))
    sw = _switches(Cr)
    info = {}
    for rw, u in list(zip(rowsA, p["uniA"])) + list(zip(rowsB, p["uniB"])):
        info.setdefault((rw[0], rw[2]), (rw[0], rw[1], rw[3], u))
    partners = {}
    for a, b, oc in zip(pa.tolist(), pb.tolist(), occ_pair.tolist()):
        partners.setdefault(a, []).append((oc, b))
        if kind == "partner":
            partners.setdefault(b, []).append((oc, a))
    tags = p.get("tags") or {}
    ours = p.get("ours")
    res_occ, res_frames, rows, meta = {}, {}, [], []

    def lbl(i):
        ch_, resi_, resn_, _u = info[allkeys[i]]
        return "%s %s/%s" % (resn_, ch_, resi_)
    for j, (kk, o) in enumerate(zip(ures, occ)):
        if o <= 0:
            continue
        ch, resi, resn, u = info[allkeys[int(kk)]]
        if kind == "partner" and ours is not None and ch not in ours:
            continue
        key = (ch, allkeys[int(kk)][1])
        res_occ[key] = float(o)
        res_frames[key] = Cr[:, j]
        best = sorted(partners.get(int(kk), []), reverse=True)[:3]
        rows.append([u if u is not None else "", "%s/%s" % (ch, resi), resn, "%.0f%%" % (100 * o), int(sw[j]),
                     "; ".join("%s %.0f%%" % (lbl(b), 100 * oo) for oo, b in best), tags.get(u, "") if u else ""])
        meta.append({"pos": u, "key": key, "chain": ch, "resi": resi})
    order = sorted(range(len(rows)), key=lambda i: (-res_occ[meta[i]["key"]], str(rows[i][1])))
    rows = [rows[i] for i in order]
    meta = [meta[i] for i in order]
    series = []
    if kind == "ligand":                             # per ligand: residues in contact per frame
        cnt = {}
        for b in pb.tolist():
            cnt[b] = cnt.get(b, 0) + 1
        for b, _n in sorted(cnt.items(), key=lambda x: -x[1])[:6]:
            series.append(("%s contacts" % lbl(b), Cp[:, pb == b].sum(1).astype(float)))
    elif kind == "partner":
        series.append(("inter-chain residue pairs", Cp.sum(1).astype(float)))
    else:
        series.append(("residues within %.1f Å" % c, Cp.sum(1).astype(float)))
    head = {"ligand": "ligand", "partner": "other chain", "nearby": "neighbour"}[kind]
    return {"series": series, "occupancy": res_occ, "contact_frames": res_frames, "cutoff": c,
            "table": (["UniProt", "res", "aa", "occ", "switch", "top %s" % head, "tags"], rows[:400]),
            "rowmeta": meta[:400],
            "summary": "%d residues in contact (heavy atoms ≤ %.1f Å)%s; %d candidate atom pairs"
                       % (len(rows), c, (" of residue %s" % p.get("center")) if kind == "nearby" else "", ncand)}


def _calc_annot(p):
    """Annotation dynamics: per annotated residue (PTM / pathogenic variant / site / hotspot / pocket),
    the fraction of frames it touches a ligand / another chain and how often that switches."""
    np = _numpy()
    tags = p.get("tags") or {}
    out = {}
    for k in ("ligand", "partner"):
        q = p.get(k) or {}
        if "reuse" in q:
            out[k] = q["reuse"]
            continue
        if q.get("error") or "A" not in q:
            out[k] = None
            continue
        out[k] = _calc_contacts(dict(q, kind=k))
    rows, meta, frames_ser = [], [], []
    for (ch, resi, resv, resn, u) in p["res"]:
        if u is None or u not in tags:
            continue
        vals = []
        for k in ("ligand", "partner"):
            r = out.get(k)
            occ = r["occupancy"].get((ch, resv), 0.0) if r else None
            fr = r["contact_frames"].get((ch, resv)) if r else None
            vals.append((occ, fr))
        if not any(v[0] for v in vals):
            continue
        sw = sum(int(_switches(v[1][:, None])[0]) for v in vals if v[1] is not None)
        rows.append([u, "%s/%s" % (ch, resi), resn, tags[u],
                     "-" if vals[0][0] is None else "%.0f%%" % (100 * vals[0][0]),
                     "-" if vals[1][0] is None else "%.0f%%" % (100 * vals[1][0]), sw])
        meta.append({"pos": u, "key": (ch, resv)})
        both = [v[1] for v in vals if v[1] is not None]
        if both:
            frames_ser.append(np.any(np.array(both), axis=0))
    series = []
    if frames_ser:
        series.append(("annotated residues in contact", np.array(frames_ser).sum(0).astype(float)))
    pct = lambda s: float(s.rstrip("%")) if s != "-" else 0.0
    order = sorted(range(len(rows)), key=lambda i: -(pct(rows[i][4]) + pct(rows[i][5])))
    return {"series": series, "sub": out,
            "table": (["UniProt", "res", "aa", "tags", "ligand", "partner", "switch"], [rows[i] for i in order]),
            "rowmeta": [meta[i] for i in order],
            "summary": "%d annotated residues contact a ligand / another chain%s" % (
                len(rows), "" if any(out.values()) else " (no ligand and no partner chain in this object)")}


def _traj_print(r):
    hd, rows = r["table"]
    print("[UFV] %s %s: %d frames (%d-%d), ref %d — %s  [read %.0f ms, compute %.0f ms]"
          % (r["obj"], r["kind"], len(r["frames"]), r["frames"][0], r["frames"][-1], r["ref"], r["summary"],
             r["prep_ms"], r["calc_ms"]))
    if rows:
        print("      " + " | ".join(str(h) for h in hd))
        for row in rows[:20]:
            print("      " + " | ".join(str(x) for x in row))
        if len(rows) > 20:
            print("      ... %d more rows (ufv_traj_csv)" % (len(rows) - 20))


def _traj_bfactor_values(obj, mode):
    """{uniPos: value} for ufv_save_bfactor modes rmsf (Å) / ligocc / partnerocc (fraction 0-1)."""
    res = _TRAJ_RESULTS.get(obj) or {}
    if mode == "rmsf":
        r = res.get("rmsf")
        if not r:
            raise RuntimeError("run ufv_traj_rmsf first")
        return {u: float(v) for (_c, _ri, _rv, _rn, u), v in zip(r["res"], r["rmsf"]) if u is not None}
    r = res.get("ligand" if mode == "ligocc" else "partner")
    if not r:
        raise RuntimeError("run ufv_traj_contacts %s first" % ("ligand" if mode == "ligocc" else "partner"))
    out = {}
    for (ch, resv), o in r["occupancy"].items():
        u = _resi_to_uni(obj, ch, resv)
        if u is not None:
            out[u] = max(out.get(u, 0.0), float(o))
    return out


# ----------------------------------------------------------------------------------------------
# Structure download for the loader (network only; the PyMOL side loads and maps the file)
# ----------------------------------------------------------------------------------------------
def _prepare_structure(struct, uid):
    """Network half of loading a listed structure (safe off the UI thread — no PyMOL calls): download
    the model to a unique temp file and, for PDB entries, fetch the SIFTS segments (and, if SIFTS has
    none, the PDBe chain hints used by the sequence-checked fallback). -> {path, segments, hints}."""
    uid = (uid or "").upper()

    def mapping():
        segments = hints = meta = None
        if struct.get("numbering") == "sifts" and struct.get("pdbId"):
            segments = _sifts_segments(struct["pdbId"], uid)
            if not segments:
                name = ((_CACHE.get(uid) or {}).get("protnlm") or {}).get("name")
                hints = _pdb_chain_hints(struct["pdbId"], uid, name)
            try:
                meta = _pdb_entry_meta(struct["pdbId"], uid,
                                       sorted(segments or {}) or struct.get("chainIds") or [])
            except Exception:
                meta = None
        return segments, hints, meta
    # model download and numbering metadata (SIFTS, residue_listing, molecules) in parallel
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=2) as ex:
        f_map = ex.submit(mapping)
        data, ext = _fetch_model_text(struct)
        segments, hints, meta = f_map.result()
    if not data:
        return None
    path = _write_temp(data, ext)
    return {"path": path, "segments": segments, "hints": hints, "meta": meta}


# ----------------------------------------------------------------------------------------------
# Command-line backend for the VMD plugin: emit annotations / SIFTS maps as Tcl, or download a
# model.  All progress output is suppressed so stdout carries only the requested payload, which
# the Tcl side captures with `exec`.
# ----------------------------------------------------------------------------------------------
def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _emit_tcl_annotations(uid):
    ann = _quiet(fetch_annotations, uid)
    out = ["set ::ufv::a_seqlen %d" % len(ann["sequence"])]

    ptm_items = []
    ptm_rich = []
    for p in ann["ptms"]:
        ptm_items.append("{%d %s}" % (p["position"], p["color"]))
        if p["endPosition"] != p["position"]:
            ptm_items.append("{%d %s}" % (p["endPosition"], p["color"]))
        desc = p.get("description") or p.get("category") or ""
        ptm_rich.append("{%d %s %s}" % (p["position"], p["color"], _tcl_braced(desc)))
    out.append("set ::ufv::a_ptms { %s }" % " ".join(ptm_items))
    out.append("set ::ufv::a_ptms_rich { %s }" % " ".join(ptm_rich))

    def _cons_token(c):
        c = c.lower()
        if "pathogenic" in c:
            return "pathogenic"
        if "benign" in c:
            return "benign"
        if "deleterious" in c:
            return "deleterious"
        return "uncertain"
    var_items = ["{%d %s %s}" % (v["position"], v["consequenceColor"], _cons_token(v["consequence"]))
                 for v in ann["variants"]]
    out.append("set ::ufv::a_variants { %s }" % " ".join(var_items))
    var_rich = []
    for v in ann["variants"]:
        wt = v.get("wildType") or "?"
        mut = v.get("mutant") or "?"
        tok = _cons_token(v["consequence"])
        dis = "; ".join(v.get("diseases") or [])
        label = "%s%d%s  %s%s" % (wt, v["position"], mut, v["consequence"],
                                  (" — " + dis) if dis else "")
        var_rich.append("{%d %s %s %s %s}" % (
            v["position"], v["consequenceColor"], tok,
            _tcl_braced("%s%s" % (wt, mut)), _tcl_braced(label)))
    out.append("set ::ufv::a_variants_rich { %s }" % " ".join(var_rich))

    site_items = []
    site_rich = []
    for s in ann["sites"]:
        site_items.append("{%d %s}" % (s["position"], s["color"]))
        if s["endPosition"] != s["position"]:
            site_items.append("{%d %s}" % (s["endPosition"], s["color"]))
        desc = s.get("description") or "Site"
        site_rich.append("{%d %s %s}" % (s["position"], s["color"], _tcl_braced(desc)))
    out.append("set ::ufv::a_sites { %s }" % " ".join(site_items))
    out.append("set ::ufv::a_sites_rich { %s }" % " ".join(site_rich))

    dom_items = ["{%d %d %s %d}" % (d["position"], d["endPosition"], d["color"], 1 if d["isRange"] else 0)
                 for d in ann["domains"]]
    out.append("set ::ufv::a_domains { %s }" % " ".join(dom_items))

    topo_items = ["{%d %d %s}" % (t["start"], t["end"], t["color"]) for t in ann["topology"]]
    out.append("set ::ufv::a_topology { %s }" % " ".join(topo_items))

    am_items = ["{%d %.4f}" % (pos, avg) for pos, avg in sorted(ann["amMean"].items())]
    out.append("set ::ufv::a_am { %s }" % " ".join(am_items))

    mut_items = ["{%d %s}" % (m["position"], m.get("color", MUTAG_COLOR)) for m in ann.get("mutagenesis", [])]
    out.append("set ::ufv::a_mutagenesis { %s }" % " ".join(mut_items))
    out.append("set ::ufv::a_burden { %s }" % " ".join(str(p) for p in sorted(ann.get("burden") or [])))
    out.append("set ::ufv::a_seq %s" % _tcl_braced(ann["sequence"]))

    fn = ann.get("function") or {}
    fbits = []
    if fn.get("summary"):
        fbits.append("Function: " + fn["summary"])
    if fn.get("locations"):
        fbits.append("Location: " + ", ".join(fn["locations"][:4]))
    if fn.get("catalytic"):
        fbits.append("Catalytic: " + "; ".join(fn["catalytic"][:2]))
    out.append("set ::ufv::a_function %s" % _tcl_braced("\n".join(fbits)))
    pn = ann.get("protnlm") or {}
    pn_txt = ("ProtNLM (AI): " + pn["name"]) if (pn.get("isAI") and pn.get("name")) else (("Name: " + pn["name"]) if pn.get("name") else "")
    out.append("set ::ufv::a_protnlm %s" % _tcl_braced(pn_txt))

    sys.stdout.write("\n".join(out) + "\n")


def _tcl_braced(s):
    """Wrap text as a safe single Tcl word (strip characters that would unbalance brace quoting)."""
    return "{%s}" % (s or "").replace("\\", "/").replace("{", "(").replace("}", ")").replace("[", "(").replace("]", ")")


def _emit_tcl_structures(uid):
    """Emit the available structures (AlphaFold / experimental PDB / computed / isoforms) as a Tcl
    list of {key label source url fmt numbering pdbId ranges} so the VMD GUI can list and load any of them.
    `ranges` (8th field) is the canonical->isoform segment list {{u_start u_end r_start r_end} ...} for
    numbering "ranges" (AlphaFold isoform models), else "-"."""
    structs = _quiet(list_structures, uid)
    items = []
    for s in structs:
        rng = s.get("ranges") if s.get("numbering") == "ranges" else None
        rtxt = ("{%s}" % " ".join("{%d %d %d %d}" % (r["u_start"], r["u_end"], r["pdb_start"], r["pdb_end"])
                                  for r in rng)) if rng else "-"
        items.append("{%s %s %s %s %s %s %s %s}" % (
            s.get("key", ""), _tcl_braced(s.get("label", "")), s.get("source", ""),
            _tcl_braced(s.get("url", "")), s.get("fmt", "pdb"),
            s.get("numbering", "identity"), s.get("pdbId") or "-", rtxt))
    sys.stdout.write("set ::ufv::a_structures { %s }\n" % " ".join(items))


def _download_url(url, fmt):
    """Download a structure file from a URL to a temp path; print the path (for the VMD loader)."""
    data = _quiet(_get_structure_text, url)
    if not data:
        sys.exit("download failed: %s" % url)
    ext = "cif" if fmt == "mmcif" else "pdb"
    path = os.path.join(_ufv_tmpdir(), "ufv_struct.%s" % ext)
    with open(path, "w") as fh:
        fh.write(data)
    sys.stdout.write(path + "\n")


def _emit_tcl_residue(uid, pos):
    """Emit a plain-text residue report (annotations + ProtVar predictors) as a Tcl variable for the
    VMD GUI's readout panel."""
    sys.stdout.write("set ::ufv::a_residue %s\n" % _tcl_braced(_residue_text(uid, pos)))


def _residue_text(uid, pos):
    """Plain-text residue report (annotations + ProtVar predictors + PDBe-KB) for VMD: the CLI's
    --residue and the --serve 'residue' op."""
    pos = int(pos)
    ann = _quiet(fetch_annotations, uid)
    seq = ann["sequence"]
    aa = seq[pos - 1] if 0 < pos <= len(seq) else "?"
    L = ["%s%d  (%s)" % (aa, pos, uid)]
    for v in ann["varByPos"].get(pos, []):
        dis = (" — " + "; ".join(v["diseases"])) if v.get("diseases") else ""
        af = v.get("gnomad")
        afs = (" · gnomAD %s" % _fmt_af(af)) if af is not None else " · not in gnomAD"
        L.append("  Variant %s%d%s — %s%s%s" % (v.get("wildType", ""), pos, v.get("mutant", ""),
                                                v.get("consequence", ""), dis, afs))
    for s in ann["sites"]:
        if s["position"] <= pos <= s["endPosition"]:
            L.append("  Site: %s" % s["description"])
    for p in ann["ptms"]:
        if p["position"] <= pos <= p["endPosition"]:
            L.append("  PTM: %s" % p["description"])
    for m in ann.get("mutagenesis", []):
        if m["position"] <= pos <= m["endPosition"]:
            L.append("  Mutagenesis: %s->%s %s" % (m.get("wildType", ""), ",".join(m.get("mutants", [])), m.get("effect", "")))
    for d in ann["domains"]:
        if d["position"] <= pos <= d["endPosition"]:
            L.append("  Domain: %s" % d["description"])
    if ann["amMean"].get(pos) is not None:
        L.append("  AlphaMissense mean: %.3f" % ann["amMean"][pos])
    try:
        pv = _quiet(fetch_protvar, uid, pos)
        sc = pv.get("score") or {}
        bits = []
        if sc.get("conservation") is not None:
            bits.append("conservation %.2f" % sc["conservation"])
        if sc.get("m3d"):
            bits.append("M3D %s" % sc["m3d"].get("prediction", "?"))
        pk = pv.get("pocket")
        if pk:
            bits.append("binding pocket (buriedness %.2f)" % (pk.get("buriedness") or 0))
        if bits:
            L.append("  ProtVar: " + " · ".join(bits))
        for p in (pk or {}).get("pockets") or []:
            L.append("  Pocket %s: %d residues%s%s" % (
                p.get("pocketId"), len(p.get("resid") or []),
                (" · score %d" % round(p["score"])) if _finite(p.get("score")) else "",
                (" · pLDDT %d" % round(p["meanPlddt"])) if _finite(p.get("meanPlddt")) else ""))
        eve = _protvar_by_mut(sc.get("eveRaw"), aa)
        esm = _protvar_by_mut(sc.get("esmRaw"), aa)
        ddg = (pv.get("foldx") or {}).get("byMut") or {}
        cadd = sc.get("caddByMut") or {}
        for v in ann["varByPos"].get(pos, []):
            mt = v.get("mutant")
            if not mt:
                continue
            seg = []
            if mt in eve and eve[mt] is not None:
                seg.append("EVE %.2f" % eve[mt])
            if mt in esm and esm[mt] is not None:
                seg.append("ESM1b %.1f" % esm[mt])
            if mt in cadd:
                seg.append("CADD %.1f" % cadd[mt])
            if mt in ddg:
                seg.append("FoldX %+.1f" % ddg[mt])
            if seg:
                L.append("    %s%d%s: %s" % (aa, pos, mt, " · ".join(seg)))
    except Exception:
        pass
    try:
        lt, nl, it, ni = _kb_residue(_quiet(fetch_pdbekb, uid), pos)
        if lt:
            L.append("  PDBe-KB ligands: " + ", ".join("%s (%d PDB)" % (x["code"], x["pdbCount"]) for x in lt))
        if it:
            L.append("  PDBe-KB interfaces: " + "; ".join(x.get("name") or x.get("accession") or "?" for x in it))
    except Exception:
        pass
    return "\n".join(L)


def _linear_runs(seg):
    """Split a per-residue (SEQRES->author) segment into linear runs {u_start, u_end, pdb_start, pdb_end}
    so a consumer with plain linear segments (the VMD plugin) maps every residue exactly."""
    if seg.get("s2a") is None:
        return [seg]
    runs, cur = [], None
    for u in range(seg["u_start"], seg["u_end"] + 1):
        a = seg["s2a"].get(seg["seqres_start"] + (u - seg["u_start"]))
        if a is None:
            cur = None
            continue
        if cur is not None and a == cur["pdb_end"] + 1:
            cur["u_end"], cur["pdb_end"] = u, a
        else:
            cur = {"u_start": u, "u_end": u, "pdb_start": a, "pdb_end": a}
            runs.append(cur)
    return runs


def _emit_tcl_sifts(pdb_id, uid):
    """SIFTS / residue_listing segments as Tcl lists for the VMD plugin. Per-residue segments are split
    into exact linear runs; "unverified" ones (no SIFTS entry: best_structures ranges translated through
    residue_listing) are left out because there is no model here to check them against the sequence."""
    segs = _quiet(_sifts_segments, pdb_id, uid.upper())
    segs = {ch: [r for s in sl if not s.get("unverified") for r in _linear_runs(s)] for ch, sl in segs.items()}
    segs = {ch: sl for ch, sl in segs.items() if sl}
    out = ["set ::ufv::s_chains { %s }" % " ".join(segs.keys())]
    for ch, segl in segs.items():
        items = ["{%s %s %s %s}" % (s["u_start"], s["u_end"], s["pdb_start"], s["pdb_end"]) for s in segl]
        out.append("set ::ufv::s_%s { %s }" % (ch, " ".join(items)))
    sys.stdout.write("\n".join(out) + "\n")


def _emit_tcl_tiers(uid, coords_file):
    """Compute hotspot / betweenness-hub / constraint-pocket tiers from exported Cα coordinates
    (written as JSON by the VMD plugin) and emit the results as Tcl variable assignments so the
    VMD plugin can colour the cartoon without PyMOL.  The three analyses run together in one call
    so the expensive coordinate scan only happens once per fetch."""
    import json
    ann = _quiet(fetch_annotations, uid)
    with open(coords_file) as fh:
        geom = json.load(fh)
    # hotspots — need variants
    try:
        hot = _per_chain_merge(geom, lambda gs: _compute_hotspots(gs, ann["variants"]))
    except Exception:
        hot = {}
    # betweenness contact hubs
    try:
        hub = _per_chain_merge(geom, _betweenness_hubs, _HUB_RANK)
    except Exception:
        hub = {}
    # constraint pockets (AlphaMissense required); classes at the default sensitivity q <= 0.10
    pocket = {}
    if ann.get("amMap"):
        try:
            pocket = _pocket_classes(_compute_pockets(geom, ann["amMap"], ann.get("sequence") or ""),
                                     POCKET_Q_DEFAULT)
        except Exception:
            pass

    def _emit(var, d):
        items = ["{%s %s}" % (pos, tier) for pos, tier in d.items()]
        sys.stdout.write("set ::ufv::%s { %s }\n" % (var, " ".join(items)))

    _emit("a_hotspots", hot)
    _emit("a_hubs", hub)
    _emit("a_pockets", pocket)


def _emit_csv_file(uid, outfile):
    """Write the annotation CSV (without structure tiers) to outfile and emit the path as a Tcl
    variable so the VMD plugin can open a save-file dialog or display a status message."""
    ann = _quiet(fetch_annotations, uid)
    text = _residue_csv_text(ann, None, None)
    with open(outfile, "w") as fh:
        fh.write(text)
    sys.stdout.write("set ::ufv::a_csv_path %s\n" % _tcl_braced(outfile))


def _download_model(uid):
    url = _quiet(_alphafold_url, uid)
    if not url:
        sys.exit("no AlphaFold model for %s" % uid)
    data = _quiet(_get_structure_text, url)
    if not data:
        sys.exit("download failed for %s" % uid)
    path = os.path.join(_ufv_tmpdir(), "ufv_%s.pdb" % uid.upper())
    with open(path, "w") as fh:
        fh.write(data)
    sys.stdout.write(path + "\n")


# ----------------------------------------------------------------------------------------------
# Persistent request server for the VMD plugin (ufv_core.py --serve).
#
# Protocol: one JSON request per line on stdin  {"id": <any>, "op": "<name>", "args": {...}, "fmt": "json"|"tcl"}
#           one reply per line on stdout          {"id": .., "ok": true, "result": ..} | {"id": .., "ok": false, "error": ".."}
# With "fmt": "tcl" the reply line is a Tcl dict instead (id <id> ok 1 result <value> | ok 0 error <msg>), every
# value encoded as a proper Tcl list/dict word (dicts -> key/value lists, None -> {}, bools -> 1/0), pure ASCII
# (\uXXXX escapes). The Tcl side reads it with `dict get` — it is data, never `eval`ed.
# Requests run on a small thread pool, so a slow network op does not hold up a quick one; replies may come
# back out of order (match them by id). stdout carries replies only: anything else printed goes to stderr.
# Ops are registered with @_serve_op; each returns plain JSON-able data (dict / list / str / number / None).
# ----------------------------------------------------------------------------------------------
SERVE_VERSION = 1
_SERVE_OPS = {}


def _serve_op(name):
    def deco(fn):
        _SERVE_OPS[name] = fn
        return fn
    return deco


_TCL_PLAIN = re.compile(r'^[A-Za-z0-9_.,:/@%+=*!?<>^~`|\'#&()-]+$')   # needs no quoting at all
_TCL_BRACEABLE = re.compile(r'^[\x20-\x5b\x5d-\x7a|~]+$')   # printable ASCII except { } \ (no newline): {..}


def _tcl_word(s):
    """One Tcl list element for the string s: bare, {braced} (printable ASCII without braces/backslashes) or
    backslash-quoted with \\n / \\uXXXX escapes. The result is always pure ASCII on one line."""
    if s == "":
        return "{}"
    if _TCL_PLAIN.match(s):
        return s
    if _TCL_BRACEABLE.match(s):
        return "{" + s + "}"
    out = []
    for ch in s:
        o = ord(ch)
        if ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        elif ch == "\r":
            out.append("\\r")
        elif ch in ' ;"$[]{}\\':
            out.append("\\" + ch)
        elif o < 32 or o == 127:
            out.append("\\u%04x" % o)
        elif o > 126:
            out.append("\\u%04x" % o if o <= 0xFFFF else "?")
        else:
            out.append(ch)
    return "".join(out)


def _tcl_scalar(v):
    if v is None:
        return ""
    if v is True:
        return "1"
    if v is False:
        return "0"
    if isinstance(v, float):
        if v != v or v in (float("inf"), float("-inf")):
            return ""
        return repr(v)
    return str(v)


def _to_tcl(v):
    """Python data -> a Tcl list/dict string (the value itself, not wrapped as an element)."""
    if isinstance(v, dict):
        return " ".join(_tcl_elem(k) + " " + _tcl_elem(x) for k, x in v.items())
    if isinstance(v, (list, tuple, set, frozenset)):
        return " ".join(_tcl_elem(x) for x in v)
    return _tcl_scalar(v)


def _tcl_elem(v):
    """Python data -> one Tcl list element. Containers are brace-wrapped (their leaves are backslash-quoted,
    so the only unescaped braces inside are the balanced ones added here). Memoised big results (see
    _encoded) are encoded once."""
    if isinstance(v, (dict, list, tuple, set, frozenset)):
        if id(v) in _MEMO_OBJS and _MEMO_OBJS[id(v)] is v:
            return _encoded(v, "tcl", lambda o: "{" + _to_tcl(o) + "}")
        return "{" + _to_tcl(v) + "}"
    return _tcl_word(_tcl_scalar(v))


def _arg_uid(a, key="uid"):
    uid = str(a.get(key) or "").strip().upper()
    if not uid:
        raise ValueError("missing '%s'" % key)
    return uid


_VMD_ANN_MEMO = {}    # uid -> (id(ann dict), structured result): rebuilt when fetch_annotations returns a new dict
_ENC_MEMO = {}        # (id(result), fmt) -> (result, encoded text): big replies are encoded once per session
_MEMO_OBJS = {}       # id -> result objects whose encodings may be memoised (immutable once built)


def _encoded(obj, fmt, enc):
    """Encoded reply text of a memoised result object (identity-keyed; the object is kept alive in the memo)."""
    k = (id(obj), fmt)
    hit = _ENC_MEMO.get(k)
    if hit is not None and hit[0] is obj:
        return hit[1]
    txt = enc(obj)
    if len(_ENC_MEMO) > 64:
        _ENC_MEMO.clear()
    _ENC_MEMO[k] = (obj, txt)
    return txt


def _vmd_annotations(uid, force=False):
    """Structured annotation set for the VMD plugin (the --emit-tcl data, plus the fields the layers need).
    Memoised per annotation dict, so the reply encoding can be memoised too (warm calls are cheap)."""
    ann = fetch_annotations(uid, force=force)
    hit = _VMD_ANN_MEMO.get(uid)
    if hit is not None and hit[0] == id(ann) and hit[2] is ann:
        return hit[1]
    res = _vmd_annotations_build(uid, ann)
    old = _VMD_ANN_MEMO.get(uid)
    if old is not None:
        _MEMO_OBJS.pop(id(old[1]), None)
    _VMD_ANN_MEMO[uid] = (id(ann), res, ann)
    _MEMO_OBJS[id(res)] = res
    return res


def _vmd_annotations_build(uid, ann):

    def tok(c):
        c = (c or "").lower()
        return ("pathogenic" if "pathogenic" in c else "benign" if "benign" in c
                else "deleterious" if "deleterious" in c else "uncertain")
    fn = ann.get("function") or {}
    fbits = []
    if fn.get("summary"):
        fbits.append("Function: " + fn["summary"])
    if fn.get("locations"):
        fbits.append("Location: " + ", ".join(fn["locations"][:4]))
    if fn.get("catalytic"):
        fbits.append("Catalytic: " + "; ".join(fn["catalytic"][:2]))
    pn = ann.get("protnlm") or {}
    pn_txt = ("ProtNLM (AI): " + pn["name"]) if (pn.get("isAI") and pn.get("name")) else \
        (("Name: " + pn["name"]) if pn.get("name") else "")
    variants = []
    for v in ann["variants"]:
        wt, mt = v.get("wildType") or "?", v.get("mutant") or "?"
        dis = v.get("diseases") or []
        variants.append({
            "pos": v["position"], "wt": wt, "mut": mt, "consequence": v.get("consequence") or "",
            "tok": tok(v.get("consequence")), "color": v.get("consequenceColor") or UNCLASSIFIED_COLOR,
            "provenance": v.get("provenance") or "", "pcolor": v.get("provenanceColor") or UNCLASSIFIED_COLOR,
            "reviewed": bool(v.get("reviewed")), "diseases": list(dis),
            "dcolor": v.get("diseaseColor") or UNCLASSIFIED_COLOR, "am": v.get("alphaMissense"),
            "gnomad": v.get("gnomad"), "clinvar": v.get("clinVar") or "",
            "label": "%s%d%s  %s%s" % (wt, v["position"], mt, v.get("consequence") or "",
                                       (" — " + "; ".join(dis)) if dis else "")})
    return {
        "uid": ann.get("uid") or uid, "seqlen": len(ann["sequence"]), "sequence": ann["sequence"],
        "ptms": [{"pos": p["position"], "end": p["endPosition"], "color": p["color"], "type": p.get("type") or "",
                  "category": p.get("category") or "", "desc": p.get("description") or p.get("category") or "",
                  "source": p.get("source") or ""} for p in ann["ptms"]],
        "variants": variants,
        "sites": [{"pos": s["position"], "end": s["endPosition"], "color": s["color"],
                   "desc": s.get("description") or "Site"} for s in ann["sites"]],
        "domains": [{"start": d["position"], "end": d["endPosition"], "color": d["color"],
                     "isRange": bool(d["isRange"]), "type": d.get("type") or "", "desc": d.get("description") or ""}
                    for d in ann["domains"]],
        "topology": [{"start": t["start"], "end": t["end"], "color": t["color"], "label": t.get("label") or ""}
                     for t in ann["topology"]],
        "mutagenesis": [{"pos": m["position"], "end": m["endPosition"], "wt": m.get("wildType") or "",
                         "mutants": list(m.get("mutants") or []), "effect": m.get("effect") or "",
                         "color": m.get("color") or MUTAG_COLOR} for m in ann.get("mutagenesis", [])],
        "am": {p: round(a, 4) for p, a in sorted(ann["amMean"].items())},
        "burden": sorted(ann.get("burden") or []),
        "diseaseColors": dict(ann.get("diseaseColors") or {}),
        "function": "\n".join(fbits), "protnlm": pn_txt,
    }


def _vmd_structure(s):
    rng = s.get("ranges") if s.get("numbering") == "ranges" else None
    return {"key": s.get("key", ""), "label": s.get("label", ""), "source": s.get("source", ""),
            "url": s.get("url", ""), "fmt": s.get("fmt", "pdb"), "numbering": s.get("numbering", "identity"),
            "pdbId": s.get("pdbId") or "", "chainIds": list(s.get("chainIds") or []),
            "chimeric": bool(s.get("chimeric")), "coverage": s.get("coverage"),
            "ranges": [[r["u_start"], r["u_end"], r["pdb_start"], r["pdb_end"]] for r in rng] if rng else [],
            "isoform": s.get("isoform") or "", "provider": s.get("provider") or ""}


def _seg_rows(sl):
    """Segment dicts (incl. per-residue SEQRES->author ones) -> exact linear runs [[us, ue, rs, re], ...]."""
    return [[r["u_start"], r["u_end"], r["pdb_start"], r["pdb_end"]] for s in sl for r in _linear_runs(s)]


def _res_rows(residues):
    """{chain: [[resid, one-letter], ...]} (from VMD) -> {chain: [(int, str)]}."""
    out = {}
    for ch, rows in (residues or {}).items():
        out[ch] = [(int(r[0]), str(r[1])) for r in rows]
    return out


def _verify_segments_res(segs, seq, res_by, notes):
    """_verify_segments (ufv_pymol) for a model given as residue rows instead of a PyMOL object."""
    for ch in list(segs):
        sl = segs[ch]
        if not any(s.get("label") or s.get("unverified") for s in sl):
            continue
        res = res_by.get(ch) or []
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
            notes.append("chain %s: PDBe numbering did not match the model; re-anchored from the sequence "
                         "(author resi = UniProt %+d, %.0f%% identity)." % (ch, off, hit[1] * 100))
            continue
        if any(s.get("unverified") for s in sl) and not (tot and match >= 0.5 * tot):
            segs.pop(ch)
            notes.append("WARNING: chain %s: residue_listing mapping does not match the UniProt sequence "
                         "(%d/%d identical) - chain left unmapped." % (ch, match, tot))
        else:
            for s in sl:
                s.pop("unverified", None)
            notes.append("WARNING: PDBe gives only label numbering for chain %s and it could not be verified "
                         "against the sequence; positions on that chain may be offset." % ch)
    return segs


def _vmd_sifts_map(pdb_id, uid, residues=None, segs=None, hints=None, meta=None):
    """_set_sifts_map (ufv_pymol) for VMD: SIFTS segments checked against the model's residues
    {chain: [[resid, aa]]} (sequence fallback when PDBe has no residue mapping). Without residues, unverified
    segments are dropped (as --emit-sifts). -> {mode, segments {chain: [[us ue rs re]..]}, status, notes, meta}."""
    uid = uid.upper()
    notes = []
    if segs is None:
        segs = _sifts_segments(pdb_id, uid)
    segs = {ch: [dict(s) for s in sl] for ch, sl in (segs or {}).items()}
    ann = _CACHE.get(uid) or fetch_annotations(uid)
    seq = ann.get("sequence", "")
    res_by = _res_rows(residues) if residues is not None else None
    if segs and res_by is not None:
        segs = _verify_segments_res(segs, seq, res_by, notes)
    elif segs:
        segs = {ch: [s for s in sl if not s.get("unverified")] for ch, sl in segs.items()}
        segs = {ch: sl for ch, sl in segs.items() if sl}
    if segs:
        status = _segments_status(segs)
    else:
        found = {}
        if res_by:
            if hints is None:
                hints = _pdb_chain_hints(pdb_id, uid, (ann.get("protnlm") or {}).get("name"))
            cands = [c for c in hints if c in res_by] if hints else list(res_by)
            for ch in cands:
                hit = _seq_offset(res_by[ch], seq)
                if hit:
                    found[ch] = [_offset_segment(res_by[ch], seq, hit[0])]
        segs = found
        if found:
            status = "Mapped by sequence (no PDBe residue mapping)."
            notes.append("no SIFTS residue mapping for %s / %s - mapped chain(s) %s by sequence." %
                         (pdb_id.upper(), uid, ", ".join(sorted(found))))
        else:
            status = "Residue mapping unavailable."
            notes.append("WARNING: no SIFTS mapping for %s / %s and no chain matched the UniProt sequence - "
                         "set the numbering by hand (ufv_chain / ufv_chainmap)." % (pdb_id.upper(), uid))
    if meta is None:
        try:
            meta = _pdb_entry_meta(pdb_id, uid, sorted(segs))
        except Exception:
            meta = None
    m = {k: (meta or {}).get(k) for k in ("chimeric", "chimericChains", "organisms", "partners", "chainNames")}
    return {"mode": "sifts", "pdbId": pdb_id.upper(), "uid": uid,
            "segments": {ch: _seg_rows(sl) for ch, sl in segs.items()}, "status": status, "notes": notes,
            "meta": m}


def _vmd_custom_map(uid, residues, configs=None, base=None, base_mode=None):
    """_apply_custom_map / _resolve_chain_config (ufv_pymol) for VMD. residues {chain: [[resid, aa]]};
    configs {chain: {mode auto|offset|adopt|ignore, start, accession, segments (adopt: the source chain's
    [[us ue rs re]..], resolved by the caller)}}; base {chain: [[us ue rs re]..]} = the previous (non-custom)
    map, kept by 'auto' when base_mode is sifts / ranges / manual."""
    uid = uid.upper()
    ann = _CACHE.get(uid) or fetch_annotations(uid)
    seq = ann.get("sequence", "")
    res_by = _res_rows(residues)
    configs = configs or {}
    segs, notes, partners, ignored, acc, cfg_all = {}, {}, {}, [], {}, {}
    for ch in sorted(res_by):
        rows = res_by[ch]
        cfg = dict(configs.get(ch) or {})
        cfg_all[ch] = cfg
        mode = (cfg.get("mode") or "auto").lower()
        a = (cfg.get("accession") or uid).upper()
        acc[ch] = a
        sl, note = None, ""
        if not rows:
            note = "no residues"
        elif mode == "ignore":
            note = "ignored"
            ignored.append(ch)
        elif not _same_acc(a, uid):
            note = "partner chain (%s)" % a
            partners[ch] = a
        elif mode == "offset":
            start = int(cfg.get("start") or 1)
            first, last = min(r for r, _ in rows), max(r for r, _ in rows)
            sl = [[start, start + (last - first), first, last]]
            note = "offset: resi %d == UniProt %d" % (first, start)
        elif mode == "adopt":
            sl = [list(map(int, s)) for s in (cfg.get("segments") or [])] or None
            note = ("adopted from %s" % (cfg.get("src") or "source")) if sl else "adopt: source has no mapping"
        else:
            if base_mode not in (None, "", "identity", "custom") and (base or {}).get(ch):
                sl = [list(map(int, s)) for s in base[ch]]
                note = "auto: kept %s mapping" % base_mode
            else:
                hit = _seq_offset(rows, seq)
                if hit:
                    s = _offset_segment(rows, seq, hit[0])
                    sl = [[s["u_start"], s["u_end"], s["pdb_start"], s["pdb_end"]]]
                    note = "auto: sequence match (resi = UniProt %+d, %.0f%% identity)" % (hit[0], hit[1] * 100)
                else:
                    note = "auto: no sequence match - set an offset or adopt a mapping"
        notes[ch] = note
        if sl:
            segs[ch] = sl
    return {"mode": "custom", "uid": uid, "segments": segs, "notes": notes, "accessions": acc,
            "partners": partners, "ignored": ignored, "config": cfg_all,
            "status": "Custom mapping: %d of %d chain(s) mapped." % (len(segs), len(res_by))}


def _vmd_layer_groups(uid, layer, only=None, provenance=None, disease=None, color=None, reviewed=False):
    """{colour hex: [UniProt positions]} for one VMD layer, computed exactly as the PyMOL commands do.
    Cartoon layers (domains, topology, alphamissense, burden) give the per-colour cartoon groups in 'groups'
    (+ 'points' marker groups for non-range domain features and 'base' = colour of unlisted residues)."""
    ann = fetch_annotations(uid)
    layer = (layer or "").lower()
    groups, points, base = {}, {}, None
    if layer in ("ptm", "ptms"):
        for p in ann["ptms"]:
            groups.setdefault(p["color"], []).append(p["position"])
            if p["endPosition"] != p["position"]:
                groups[p["color"]].append(p["endPosition"])
    elif layer in ("var", "variants"):
        o = (only if only is not None else "pathogenic") or ""
        o = o.lower() if isinstance(o, str) else " ".join(o).lower()
        tokens = None if o.strip() in ("", "all") else set(re.split(r"[ ,]+", o.strip()))
        groups = _variant_groups(ann, tokens, bool(reviewed), _parse_provenances(provenance),
                                 _parse_diseases(ann, disease), color or "consequence")
    elif layer in ("site", "sites"):
        for s in ann["sites"]:
            groups.setdefault(s["color"], []).append(s["position"])
            if s["endPosition"] != s["position"]:
                groups[s["color"]].append(s["endPosition"])
    elif layer in ("mutag", "mutagenesis"):
        for m in ann.get("mutagenesis", []):
            groups.setdefault(m.get("color") or MUTAG_COLOR, []).extend(range(m["position"], m["endPosition"] + 1))
    elif layer in ("domain", "domains"):
        pos_color = {}
        for d in ann["domains"]:
            if d["isRange"]:
                for p in range(d["position"], d["endPosition"] + 1):
                    pos_color[p] = d["color"]
            else:
                points.setdefault(d["color"], []).append(d["position"])
        for p, c in pos_color.items():
            groups.setdefault(c, []).append(p)
    elif layer == "topology":
        pos_color = {}
        for t in ann["topology"]:
            for p in range(t["start"], t["end"] + 1):
                pos_color[p] = t["color"]
        for p, c in pos_color.items():
            groups.setdefault(c, []).append(p)
    elif layer in ("am", "alphamissense"):
        base = "#b9c2cf"
        for pos, avg in ann["amMean"].items():
            c = ("#b71c1c" if avg >= 0.78 else "#e06666" if avg >= 0.564 else "#b9c2cf" if avg >= 0.34 else "#3d85c8")
            groups.setdefault(c, []).append(pos)
    elif layer == "burden":
        if ann.get("burden"):
            groups["#e65100"] = list(ann["burden"])
    else:
        raise ValueError("unknown layer '%s'" % layer)
    srt = lambda g: {c: sorted(set(v)) for c, v in g.items()}
    return {"layer": layer, "groups": srt(groups), "points": srt(points), "base": base,
            "count": sum(len(set(v)) for v in groups.values())}


@_serve_op("ping")
def _op_ping(a):
    return {"version": SERVE_VERSION, "python": sys.version.split()[0], "executable": sys.executable,
            "pid": os.getpid(), "core": os.path.abspath(__file__), "ops": sorted(_SERVE_OPS),
            "cache": _DISK_CACHE_DIR, "numpy": _numpy() is not None}


@_serve_op("annotations")
def _op_annotations(a):
    return _vmd_annotations(_arg_uid(a), bool(a.get("force")))


@_serve_op("structures")
def _op_structures(a):
    return [_vmd_structure(s) for s in list_structures(_arg_uid(a), force=bool(a.get("force")))]


@_serve_op("fetch")
def _op_fetch(a):
    """Annotations + (unless structures=false) the structure list, in ONE request (the Fetch button)."""
    uid = _arg_uid(a)
    out = {"annotations": _vmd_annotations(uid, bool(a.get("force")))}   # (only 'annotations' replies are memo-encoded)
    if a.get("structures", True) not in (False, 0, "0", "false"):
        try:
            out["structures"] = [_vmd_structure(s) for s in list_structures(uid, len(out["annotations"]["sequence"]))]
        except Exception as e:      # annotations are still useful without the list
            out["structures"] = []
            out["structures_error"] = str(e)
        out["structures_notice"] = structure_notice(uid)
    return out


@_serve_op("structure_notice")
def _op_structure_notice(a):
    """One-line notice about the last structure listing of uid ("" = complete), e.g. 3D-Beacons unreachable."""
    return structure_notice(_arg_uid(a))


@_serve_op("download")
def _op_download(a):
    """AlphaFold model of uid -> {path} (a unique temp file)."""
    uid = _arg_uid(a)
    url = _alphafold_url(uid)
    if not url:
        raise RuntimeError("no AlphaFold model for %s" % uid)
    data = _get_structure_text(url)
    if not data:
        raise RuntimeError("download failed for %s" % uid)
    return {"path": _write_temp(data, "pdb"), "url": url}


@_serve_op("get")
def _op_get(a):
    url = str(a.get("url") or "")
    data = _get_structure_text(url)
    if not data:
        raise RuntimeError("download failed: %s" % url)
    return {"path": _write_temp(data, "cif" if a.get("fmt") == "mmcif" else "pdb")}


@_serve_op("prepare")
def _op_prepare(a):
    """Download a listed structure (by key) + its SIFTS segments in parallel -> {path, struct, segments, hints,
    meta}; then call map_sifts with the loaded model's residues (it reuses the memoised SIFTS data)."""
    uid = _arg_uid(a)
    key = str(a.get("key") or "")
    st = next((s for s in list_structures(uid) if s.get("key") == key), None)
    if st is None:
        raise ValueError("structure %s not listed for %s" % (key, uid))
    r = _prepare_structure(st, uid)
    if not r:
        raise RuntimeError("download failed for %s" % key)
    return {"path": r["path"], "struct": _vmd_structure(st),
            "segments": {ch: _seg_rows(sl) for ch, sl in (r.get("segments") or {}).items()},
            "hints": r.get("hints")}


@_serve_op("map_sifts")
def _op_map_sifts(a):
    return _vmd_sifts_map(str(a.get("pdb") or ""), _arg_uid(a), a.get("residues"))


@_serve_op("map_custom")
def _op_map_custom(a):
    return _vmd_custom_map(_arg_uid(a), a.get("residues") or {}, a.get("configs"), a.get("base"), a.get("base_mode"))


@_serve_op("seq_offset")
def _op_seq_offset(a):
    """{offset, identity, matches} of one chain's residues [[resid, aa]..] against uid's sequence, or {} ."""
    uid = _arg_uid(a)
    seq = (_CACHE.get(uid) or fetch_annotations(uid)).get("sequence", "")
    hit = _seq_offset([(int(r[0]), str(r[1])) for r in a.get("residues") or []], seq)
    return {"offset": hit[0], "identity": hit[1], "matches": hit[2]} if hit else {}


@_serve_op("layer_groups")
def _op_layer_groups(a):
    return _vmd_layer_groups(_arg_uid(a), a.get("layer"), a.get("only"), a.get("provenance"), a.get("disease"),
                             a.get("color"), a.get("reviewed"))


@_serve_op("residue")
def _op_residue(a):
    return {"pos": int(a.get("pos")), "text": _residue_text(_arg_uid(a), int(a.get("pos")))}


@_serve_op("tiers")
def _op_tiers(a):
    """Hotspot / hub / constraint-pocket tiers from Cα geometry: `coords` = path of the JSON the plugin writes
    ([{uniPos, chain, resi, x, y, z}]) or `geom` inline. q = pocket sensitivity (default 0.10)."""
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    if a.get("map") is not None or a.get("struct") is not None:
        # molecule context given (::ufv::feat::ctx): exactly the PyMOL plugin's inputs (AlphaFold PAE gate incl.
        # isoforms, partner-protein residues for hotspots) through the 'analyses' path
        r = _op_analyses(dict(a, kinds="hotspots hubs pockets"))
        q = float(a.get("q") or POCKET_Q_DEFAULT)
        return {"hotspots": r.get("hotspots") or {}, "hubs": r.get("hubs") or {},
                "pockets": _pocket_classes(r.get("pockets") or {}, q), "errors": r.get("errors") or {}}
    if a.get("geom") is not None:
        geom = a["geom"]
    else:
        with open(str(a.get("coords"))) as fh:
            geom = json.load(fh)
    out = {"hotspots": {}, "hubs": {}, "pockets": {}, "errors": {}}
    try:
        out["hotspots"] = _per_chain_merge(geom, lambda gs: _compute_hotspots(gs, ann["variants"]))
    except Exception as e:
        out["errors"]["hotspots"] = str(e)
    try:
        out["hubs"] = _per_chain_merge(geom, _betweenness_hubs, _HUB_RANK)
    except Exception as e:
        out["errors"]["hubs"] = str(e)
    if ann.get("amMap"):
        try:
            q = float(a.get("q") or POCKET_Q_DEFAULT)
            out["pockets"] = _pocket_classes(_compute_pockets(geom, ann["amMap"], ann.get("sequence") or ""), q)
        except Exception as e:
            out["errors"]["pockets"] = str(e)
    return out


@_serve_op("csv")
def _op_csv(a):
    uid = _arg_uid(a)
    out = str(a.get("out") or "")
    if not out:
        raise ValueError("missing 'out'")
    text = _residue_csv_text(fetch_annotations(uid), None, None)
    with open(out, "w") as fh:
        fh.write(text)
    return {"path": out}


# ---- VMD features (vmd/ufv_features.tcl, vmd/ufv_commands.tcl) ---------------------------------------------------
# The PyMOL glue's object-based helpers (numbering, partner chains, AlphaFold model info, analyses, motif search,
# CSV) run here on a "virtual object" that carries the VMD molecule's numbering / structure / accession from the
# request; Cα geometry (current frame) comes from the plugin's export ([{uniPos, chain, resi, x, y, z, b}]).
_VMD_AN = {}            # (vobj, geometry hash, kind) -> analysis result (memo, bounded)
_VMD_PROGRESS = {}      # progress key -> [done, total] for long crawls (pockets_all, protvar CSV)
VMD_FLAG_DEFS = (("hotspot", "Pathogenic variant hotspot", "#e53935"),
                 ("burden", "Recurrent phenotype residue", "#e65100"),
                 ("hub", "Contact-network centrality", "#6a1b9a"),
                 ("pocket", "Burial-adjusted constraint cluster", "#00897b"))


def _vmd_vobj(a, uid):
    """Register (once) the virtual object of a VMD molecule: map {mode, segments {ch: [[us ue rs re]..]},
    exclude [chains]} + struct {key source url isoform pdbId provider label} + uid. SIFTS entries get the entry's
    partner / chimeric metadata (memoised PDBe calls). -> object name."""
    m = a.get("map") or {}
    st = {k: v for k, v in (a.get("struct") or {}).items() if v not in (None, "", "-")}
    mode = str(m.get("mode") or "identity")
    segs = {str(ch): [[int(x) for x in s] for s in sl] for ch, sl in (m.get("segments") or {}).items()}
    excl = sorted(str(c) for c in (m.get("exclude") or []))
    pdb = str(m.get("pdbId") or st.get("pdbId") or "").upper()
    key = json.dumps([uid, mode, segs, excl, pdb, sorted(st.items())], sort_keys=True, default=str)
    name = "vmd:" + _sha1hex(key.encode("utf-8"))[:20]
    if name in UFV_MAPS:
        return name
    if mode == "identity":
        mp = {"mode": "identity", "chains": {}}
        if excl:
            mp["exclude"] = frozenset(excl)
    else:
        mp = {"mode": "custom", "segments": {ch: [{"u_start": s[0], "u_end": s[1], "pdb_start": s[2], "pdb_end": s[3]}
                                                  for s in sl] for ch, sl in segs.items()}}
    mp["uid"] = uid
    if pdb:
        mp["pdbId"] = pdb
        try:
            meta = _pdb_entry_meta(pdb, uid, sorted(segs))
        except Exception:
            meta = None
        for k in ("chimeric", "chimericChains", "organisms", "partners", "chainNames"):
            if meta and meta.get(k) is not None:
                mp[k] = meta.get(k)
    if st.get("isoform"):
        mp["isoform"] = st["isoform"]
    if pdb and not st.get("pdbId"):
        st["pdbId"] = pdb
    _STATE.setdefault("obj_struct", {})[name] = st
    _STATE.setdefault("obj_uid", {})[name] = uid
    UFV_MAPS[name] = mp
    return name


def _vmd_geom(a):
    """(geometry, hash) from `coords` (path of the plugin's Cα JSON) or inline `geom`; ([], '') if neither."""
    if a.get("geom") is not None:
        raw = json.dumps(a["geom"], sort_keys=True).encode("utf-8")
        geom = a["geom"]
    elif a.get("coords"):
        with open(str(a["coords"]), "rb") as fh:
            raw = fh.read()
        geom = json.loads(raw.decode("utf-8"))
    else:
        return [], ""
    for g in geom:
        g.setdefault("b", None)
        if g.get("uniPos") is not None:
            g["uniPos"] = int(g["uniPos"])
        g["resi"] = int(g["resi"])
    return geom, _sha1hex(raw)


def _vmd_radii(a):
    return (float(a.get("ptm_radius") or PROX_SETTINGS["ptm"]), float(a.get("var_radius") or PROX_SETTINGS["var"]))


def _vmd_analysis(vobj, uid, geom, ghash, kind, ann=None):
    """The PyMOL glue's _analysis for a VMD molecule (memoised per virtual object + geometry):
    hotspots / hubs -> {tiers, byChain}; pockets -> the full computePockets result (every candidate's q)."""
    key = (vobj, ghash, kind)
    hit = _VMD_AN.get(key)
    if hit is not None:
        return hit
    if kind == "hubs":
        bych = {}
        res = {"tiers": _per_chain_merge(geom, _betweenness_hubs, _HUB_RANK, bych), "byChain": bych}
    else:
        ann = ann or _CACHE.get(uid) or fetch_annotations(uid)
        if kind == "hotspots":
            variants = ann.get("variants", [])
            classes = _hotspot_classes(variants)
            partners = []
            if len(classes[0]) >= HOTSPOT_MIN_PATH and _partner_maps(vobj):
                try:
                    partners = _resolve_partner_points(geom, _partner_points(vobj))
                except Exception:
                    partners = []
            bych = {}
            res = {"tiers": _per_chain_merge(geom, lambda gs: _compute_hotspots(gs, variants, partners=partners,
                                                                                  classes=classes), _TIER_RANK, bych),
                   "byChain": bych}
        elif kind == "pockets":
            pae = None
            info = _af_model_info(vobj)
            if info and ann.get("amMap"):
                try:
                    pae = _pae_matrix(*info)
                except Exception:
                    pae = None
            res = _compute_pockets(geom, ann.get("amMap") or {}, ann.get("sequence") or "", pae)
        else:
            raise ValueError("unknown analysis '%s'" % kind)
    if len(_VMD_AN) > 96:
        _VMD_AN.clear()
    _VMD_AN[key] = res
    return res


def _vmd_pockets_reply(res):
    return {"byPos": {u: {"q": v.get("q"), "cat": v.get("cat")} for u, v in (res.get("byPos") or {}).items()},
            "hasPae": bool(res.get("hasPae")), "reason": res.get("reason"), "n": res.get("n") or 0}


@_serve_op("analyses")
def _op_analyses(a):
    """Structure analyses of a VMD molecule exactly as the PyMOL commands compute them: hotspots (partner-protein
    disease residues of a complex folded in), contact hubs, constraint pockets (PAE-gated on AlphaFold models;
    every candidate's q, so the sensitivity slider recolours without recomputing). kinds = hotspots hubs pockets."""
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    vobj = _vmd_vobj(a, uid)
    geom, gh = _vmd_geom(a)
    kinds = a.get("kinds") or "hotspots hubs pockets"
    kinds = re.split(r"[ ,]+", kinds.strip()) if isinstance(kinds, str) else list(kinds)
    out = {"errors": {}}
    for k in kinds:
        try:
            if k == "pockets" and not ann.get("amMap"):
                out["pockets"] = {"byPos": {}, "hasPae": False, "reason": "No AlphaMissense data for this protein.",
                                  "n": 0}
                continue
            r = _vmd_analysis(vobj, uid, geom, gh, k, ann)
            if k == "pockets":
                out[k] = _vmd_pockets_reply(r)
            else:
                out[k] = r["tiers"]
                out[k + "ByChain"] = r["byChain"]
        except Exception as e:
            out["errors"][k] = "%s: %s" % (type(e).__name__, e)
    return out


@_serve_op("diseases")
def _op_diseases(a):
    """ufv_diseases: [{name, count, color}] in first-appearance order (data.js getDiseaseSummary)."""
    ann = fetch_annotations(_arg_uid(a))
    return [{"name": n, "count": d["count"], "color": d["color"]}
            for n, d in _disease_summary(ann["variants"], ann.get("diseaseColors")).items()]


@_serve_op("partners")
def _op_partners(a):
    """Partner chains (other UniProt entries) of a SIFTS-mapped entry: [{accession chainId uniprotStart uniprotEnd
    name}] + chimeric / organisms of our chains."""
    uid = _arg_uid(a)
    vobj = _vmd_vobj(a, uid)
    mp = UFV_MAPS.get(vobj) or {}
    return {"partners": [{k: p.get(k) for k in ("accession", "chainId", "uniprotStart", "uniprotEnd", "name")}
                         for p in mp.get("partners") or []],
            "chimeric": bool(mp.get("chimeric")), "chimericChains": list(mp.get("chimericChains") or []),
            "organisms": list(mp.get("organisms") or []), "pdbId": mp.get("pdbId") or ""}


@_serve_op("partner_layers")
def _op_partner_layers(a):
    """ufv_partner_layers: one sphere per partner-protein residue coloured by the top enabled layer (pathogenic
    variants > PTMs > sites) -> {spheres [[chain resi hex tooltip]..], chains n, accessions n}."""
    uid = _arg_uid(a)
    vobj = _vmd_vobj(a, uid)
    want = _parse_partner_layers(a.get("layers") or "all")
    if not want or not _partner_maps(vobj):
        return {"spheres": [], "chains": 0, "accessions": len(_partner_maps(vobj)), "layers": sorted(want)}
    cols = _partner_sphere_colors(_partner_annotations(vobj), want)
    return {"spheres": [[ch, int(r), c, t] for (ch, r), (c, t) in sorted(cols.items(), key=lambda kv: (kv[0][0], kv[0][1]))],
            "chains": len({k[0] for k in cols}), "accessions": len(_partner_maps(vobj)),
            "layers": [l for l in PARTNER_LAYERS if l in want]}


@_serve_op("progress")
def _op_progress(a):
    p = _VMD_PROGRESS.get(str(a.get("key") or ""))
    return {"done": p[0], "total": p[1]} if p else {}


@_serve_op("pockets_all")
def _op_pockets_all(a):
    """ufv_pockets_all: every ProtVar (AutoSite) pocket of the protein (crawl, cached), best first, with the
    8-colour palette colour of its rank. Progress: op 'progress' with key 'pockets_all:<UID>'."""
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    key = "pockets_all:" + uid
    _VMD_PROGRESS[key] = [0, len(ann.get("sequence") or "")]

    def prog(d, t):
        _VMD_PROGRESS[key] = [d, t]
    res = fetch_all_pockets(uid, len(ann.get("sequence") or ""), prog)
    if not res:
        return {"count": 0, "pockets": []}
    pk = []
    for i, p in enumerate(res["pockets"]):
        q = {k: (p.get(k) if _finite(p.get(k)) else None) for k in
             ("score", "meanPlddt", "radGyration", "energyPerVol", "buriedness")}
        q.update(pocketId=p.get("pocketId"), resid=[int(r) for r in p.get("resid") or []],
                 color=POCKET_PALETTE[i % len(POCKET_PALETTE)])
        pk.append(q)
    return {"count": res["count"], "pockets": pk}


@_serve_op("pocket_find")
def _op_pocket_find(a):
    """A ProtVar pocket already fetched (crawl or a residue report) by id -> {pocketId resid ...} or {}."""
    p = _find_pocket(_arg_uid(a), a.get("id"))
    return {"pocketId": p.get("pocketId"), "resid": [int(r) for r in p.get("resid") or []],
            "score": p.get("score") if _finite(p.get("score")) else None} if p else {}


def _vmd_motif_url(vobj, positions, uid, geom):
    """motif_search_url (PyMOL glue) for a virtual object: up to 8 residues closest to the pocket's Cα centroid;
    PDB entries -> mmCIF label ids from residue_listing; AlphaFold / AlphaFill -> AF_AF<acc>F1, chain A."""
    st = _obj_struct(vobj)
    pdb_id = st.get("pdbId") or (UFV_MAPS.get(vobj) or {}).get("pdbId")
    is_af = not pdb_id and _is_af_based(vobj)
    if not pdb_id and not is_af:
        return None
    geo = {}
    for g in geom or []:
        if g["uniPos"] is not None:
            geo[g["uniPos"]] = g
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
        acc = _canon_acc(uid or "")
        if not acc:
            return None
        entry = "AF_AF%sF1" % acc
        rid = [{"label_asym_id": "A", "label_seq_id": g["uniPos"]} for g in picked]
    if len(rid) < 2:
        return None
    q = {"query": {"type": "terminal", "service": "strucmotif",
                   "parameters": {"value": {"entry_id": entry, "residue_ids": rid}}}, "return_type": "assembly"}
    return "https://www.rcsb.org/search?request=" + _uri(json.dumps(q, separators=(",", ":")))


@_serve_op("motif_url")
def _op_motif_url(a):
    """ufv_motif_search: RCSB strucmotif URL for positions (or pocket=<ProtVar id>) -> {url} ('' if unavailable)."""
    uid = _arg_uid(a)
    vobj = _vmd_vobj(a, uid)
    geom, _gh = _vmd_geom(a)
    pos = a.get("positions")
    if a.get("pocket") not in (None, ""):
        pk = _find_pocket(uid, a["pocket"])
        pos = pk["resid"] if pk else []
    if isinstance(pos, str):
        pos = [p for p in re.split(r"[\s,+]+", pos) if p]
    return {"url": _vmd_motif_url(vobj, pos or [], uid, geom) or "", "ok": _structure_motif_ok(vobj)}


@_serve_op("pdbekb")
def _op_pdbekb(a):
    """PDBe-KB layer groups: {ligand [pos..], interface [pos..], url} (positions seen binding a ligand / in a PPI
    interface in any PDB entry of the protein)."""
    uid = _arg_uid(a)
    kb = fetch_pdbekb(uid)
    return {"ligand": sorted(int(p) for p in (kb or {}).get("ligandByPos") or {}),
            "interface": sorted(int(p) for p in (kb or {}).get("interfaceByPos") or {}),
            "url": "https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/%s" % _canon_acc(uid), "found": bool(kb)}


@_serve_op("binding")
def _op_binding(a):
    """ufv_binding: PDBe-KB ligands (top 8) and PPI interfaces (top 4) at one residue."""
    uid = _arg_uid(a)
    lt, nl, it, ni = _kb_residue(fetch_pdbekb(uid), int(a.get("pos")))
    return {"ligands": lt, "nLigands": nl, "interfaces": it, "nInterfaces": ni,
            "url": "https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/%s" % _canon_acc(uid)}


@_serve_op("tractability")
def _op_tractability(a):
    """Open Targets tractability {symbol ensemblId groups [{modality labels}] url} or {}."""
    return fetch_open_targets(_arg_uid(a)) or {}


@_serve_op("alphafill")
def _op_alphafill(a):
    """AlphaFill transplant metadata {asym/auth chain: {identity clash compound donorPdb rmsd}} or {}."""
    return alphafill_transplants(_arg_uid(a)) or {}


@_serve_op("ligand_info")
def _op_ligand_info(a):
    """ligand_chemistry for CCD codes in structure order: {code: {name formula weight hbondDonors hbondAcceptors
    smiles inchikey drugbank cid links [[label url]..] similar [[code score]..]}}."""
    codes = a.get("codes") or []
    codes = [c for c in re.split(r"[ ,]+", codes.strip()) if c] if isinstance(codes, str) else [str(c) for c in codes]
    codes = list(dict.fromkeys(c.upper() for c in codes))
    if not codes:
        return {}

    def one(c):
        inf = ligand_info(c)
        ik = inf.get("inchikey")
        return inf, ligand_fingerprint(ik), (pubchem_cid(ik) if ik else None)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=8) as ex:
        got = dict(zip(codes, ex.map(one, codes)))
    out = {}
    for c in codes:
        inf, fp, cid = got[c]
        sims = []
        if fp:
            for d in [x for x in codes if x != c][:LIGAND_SIM_CAP]:
                t = tanimoto(fp, got[d][1])
                if t is not None and t > 0:
                    sims.append([d, round(t, 3)])
            sims.sort(key=lambda x: -x[1])
        out[c] = dict(inf, similar=sims, cid=cid, links=[list(l) for l in ligand_links(inf, cid)])
    return out


@_serve_op("pocket_evidence")
def _op_pocket_evidence(a):
    """pocket_evidence for UniProt positions (ligand pocket / ProtVar pocket): mean Cα B (pLDDT on AlphaFold
    models) + band, functional sites hit, residue-class composition (+ text)."""
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    geom, _gh = _vmd_geom(a)
    pos = a.get("positions") or []
    if isinstance(pos, str):
        pos = [p for p in re.split(r"[\s,+]+", pos) if p]
    ev = pocket_evidence([int(p) for p in pos], geom, ann.get("sites") or [], ann.get("sequence") or "")
    mb = ev.get("meanB")
    return {"meanB": mb, "count": ev["count"], "band": plddt_band(mb) if mb is not None else "",
            "sites": [{"pos": s["position"], "end": s.get("endPosition") or s["position"],
                       "desc": s.get("description") or "site"} for s in ev["siteHits"]],
            "composition": ev["composition"] or {}, "compositionText": composition_text(ev["composition"])}


def _vmd_report(a):
    """residue_report + format_report (PyMOL glue) for a VMD molecule. The plugin passes the Nearby search it ran
    in VMD (current frame), the ligands within 5 A, the Cα geometry and the pocket sensitivity q; with services=1
    ProtVar (predictors, CADD, AutoPocket) and PDBe-KB are fetched too. -> structured sections for the panel +
    the console text."""
    uid = _arg_uid(a)
    pos = int(a.get("pos"))
    ann = fetch_annotations(uid)
    vobj = _vmd_vobj(a, uid) if (a.get("map") or a.get("struct")) else None
    geom, gh = _vmd_geom(a)
    seq = ann["sequence"]
    aa = seq[pos - 1] if 0 < pos <= len(seq) else "?"
    pr, vr = _vmd_radii(a)
    rp = _residue_proximity(pos, ann["ptms"], ann["variants"], geom, pr, vr) if geom else {"ptms": [], "variants": []}
    prox_ptms = [{"pos": p["position"], "dist": d, "category": p.get("category"), "color": p.get("color"),
                  "description": p.get("description")} for p, d in rp["ptms"]]
    prox_vars = [{"pos": v["position"], "dist": d, "label": _vlabel(v), "consequence": v.get("consequence"),
                  "color": v.get("consequenceColor"), "diseases": v.get("diseases") or []} for v, d in rp["variants"]]
    proximity = None
    ptms_here = [p for p in ann["ptms"] if p["position"] == pos]
    if ptms_here and geom:
        proximity = _ptm_variant_proximity(ptms_here, ann["variants"], geom, pr, vr).get(pos)
    nbin = a.get("nearby") or {}
    near = sorted({int(p) for p in nbin.get("near") or []} - {pos})
    pann = (UFV_MAPS.get(vobj) or {}).get("_partnerAnn") if vobj else None
    plist = []
    for row_ in nbin.get("partners") or []:
        ch, resi = str(row_[0]), str(row_[1])
        info = None
        if pann:
            r = _resi_int(resi)
            p = next((x for x in pann if x["chainId"] == ch and r in x["byResi"]), None)
            if p is None:
                continue
            info = dict(p["byResi"][r], accession=p["accession"])
        plist.append((ch, resi, info))
    nb = {"dist": float(nbin.get("dist") or 5.0), "near": near, "partners": plist}
    var_color = {}
    for v in ann["variants"]:
        if v["position"] in set(near):
            var_color[v["position"]] = v.get("consequenceColor")
    # flags: only analyses already computed for this geometry (memo), pocket class at the plugin's q
    q = float(a.get("q") or POCKET_SETTINGS["q"])
    ran = a.get("ran")            # analyses the user ran on this molecule (None = any memoised one)
    ran = None if ran is None else set(re.split(r"[ ,]+", str(ran).strip())) - {""}
    flags = []
    for k, lbl, col in VMD_FLAG_DEFS:
        on = None
        kind = {"hotspot": "hotspots", "hub": "hubs", "pocket": "pockets"}.get(k)
        if k == "burden":
            on = None if ann.get("burden") is None else pos in ann["burden"]
        elif vobj and gh and (ran is None or kind in ran):
            r = _VMD_AN.get((vobj, gh, {"hotspot": "hotspots", "hub": "hubs", "pocket": "pockets"}[k]))
            if r is not None:
                on = pos in (_pocket_classes(r, q) if k == "pocket" else r["tiers"])
        if on is not None:
            flags.append((k, lbl, col, bool(on)))
    rep = {"uid": ann["uid"], "pos": pos, "aa": aa, "proximity": proximity, "proxPtms": prox_ptms,
           "proxVariants": prox_vars, "proxRadii": (pr, vr),
           "ptms": [p for p in ann["ptms"] if p["position"] <= pos <= p["endPosition"]],
           "variants": ann["varByPos"].get(pos, []),
           "sites": [s for s in ann["sites"] if s["position"] <= pos <= s["endPosition"]],
           "domains": [d for d in ann["domains"] if d["position"] <= pos <= d["endPosition"]],
           "mutagenesis": [m for m in ann.get("mutagenesis", []) if m["position"] <= pos <= m["endPosition"]],
           "amMean": ann["amMean"].get(pos), "amProfile": ann.get("amByPos", {}).get(pos, []),
           "nearbyLigands": sorted(set(str(x) for x in a.get("ligands") or [])), "obj": vobj, "nearby": nb,
           "flags": flags, "varColor": var_color}
    services = str(a.get("services", 1)) not in ("0", "false", "False", "")
    pv = kb = None
    if services:
        try:
            pv = fetch_protvar(uid, pos)
        except Exception:
            pv = None
        try:
            kb = fetch_pdbekb(uid)
        except Exception:
            kb = None
    # console text: the VMD report heading ("N1378  (P35498)", as the CLI --residue / 'residue' op) + format_report body
    text = ["%s%d  (%s)" % (aa, pos, rep["uid"])] + format_report(rep).split("\n")[1:]
    S = []

    def sec(title, color, rows):
        if rows:
            S.append({"title": title, "color": color, "rows": rows})

    def row(label, value, color="", **kw):
        d = {"label": label, "value": value, "color": color or ""}
        d.update(kw)
        return d
    rows = []
    for v in rep["variants"]:
        lab = "%s%d%s" % (v.get("wildType", ""), pos, v.get("mutant", ""))
        cons = v["consequence"] + ((" (%s)" % "; ".join(v["diseases"])) if v.get("diseases") else "")
        rows.append(row(lab, cons, v.get("consequenceColor") or UNCLASSIFIED_COLOR, bold=1))
        if v.get("alphaMissense") is not None:
            rows.append(row("AlphaMissense", "%.2f" % v["alphaMissense"], _am_color(v["alphaMissense"])))
        af = v.get("gnomad")
        rows.append(row("gnomAD AF", ("%s (%s)" % (_fmt_af(af), "common" if af >= 0.01 else "rare" if af > 0 else "absent"))
                        if af is not None else "not in gnomAD", "#388e3c" if (af or 0) >= 0.01 else "#777"))
        for k_, lab_ in (("clinVar", "ClinVar"), ("clinVarReview", "Review")):
            if v.get(k_):
                rows.append(row(lab_, str(v[k_]), evidence=1))
        if v.get("rsIds"):
            rows.append(row("dbSNP", ", ".join(v["rsIds"]), evidence=1))
    sec("Variants", "#c62828", rows)
    rows = [row("PTM", p["description"], p.get("color")) for p in rep["ptms"]]
    rows += [row("Site", s["description"], SITE_COLOR) for s in rep["sites"]]
    rows += [row("Domain", d["description"], d.get("color")) for d in rep["domains"]]
    rows += [row("Mutagenesis", "%s->%s %s" % (m.get("wildType", ""), ", ".join(m.get("mutants", [])), m.get("effect", "")),
                 MUTAG_COLOR) for m in rep["mutagenesis"]]
    sec("Features", "#1565c0", rows)
    if proximity:
        px = proximity
        tl = {1: "variant at this residue", 2: "pathogenic variant <= %g A" % pr, 3: "pathogenic variant <= %g A" % vr}
        rows = [row("Tier", tl.get(px.get("tier"), str(px.get("tier"))))]
        if px.get("nearestVariant"):
            nd = px.get("nearestDist")
            rows.append(row("Nearest variant", "%s (%s)" % (px["nearestVariant"], "same residue" if nd == 0 else
                                                            ("%.1f A" % nd if nd is not None else ""))))
        rows.append(row("Within %g A" % pr, "%d variant(s), %d pathogenic" % (px.get("nearbyCount8A", 0),
                                                                            px.get("pathCount8A", 0))))
        sec("PTM-variant proximity", "#b87800", rows)
    # prediction table (AM always; ProtVar EVE / ESM1b / CADD / FoldX when fetched)
    sc = (pv or {}).get("score") or {}
    fx = (pv or {}).get("foldx") or {}
    rows = []
    if rep["amMean"] is not None:
        rows.append(row("AlphaMissense mean", "%.3f" % rep["amMean"], _am_color(rep["amMean"])))
    if sc.get("conservation") is not None:
        rows.append(row("Conservation", "%.2f" % sc["conservation"]))
    if sc.get("m3d"):
        md = sc["m3d"]
        rows.append(row("M3D", "%s%s" % (md.get("prediction", "?"), (" - %s" % md["feature"]) if md.get("feature") else ""),
                        "#d32f2f" if str(md.get("prediction", "")).lower().startswith("damag") else "#388e3c"))
    am_by = {m: s for (_w, m, s) in rep["amProfile"]}
    eve = _protvar_by_mut(sc.get("eveRaw"), aa) if pv else {}
    esm = _protvar_by_mut(sc.get("esmRaw"), aa) if pv else {}
    ddg = fx.get("byMut") or {}
    cadd = sc.get("caddByMut") or {}
    obs = {v.get("mutant") for v in rep["variants"] if v.get("mutant")}
    table = []
    for mt in [x for x in AA_ALPHA if x != aa]:
        r_ = [mt, am_by.get(mt), eve.get(mt), esm.get(mt), cadd.get(mt), ddg.get(mt)]
        if all(x is None for x in r_[1:]):
            continue
        table.append(r_)
    table.sort(key=lambda r_: (r_[1] is None, -(r_[1] if r_[1] is not None else 0), -(r_[2] if r_[2] is not None else 0)))
    cols = [c for i, c in enumerate(("AM", "EVE", "ESM1b", "CADD", "FoldX"), 1) if any(r_[i] is not None for r_ in table)]
    idx = {"AM": 1, "EVE": 2, "ESM1b": 3, "CADD": 4, "FoldX": 5}
    sub = {"columns": ["sub"] + cols,
           "rows": [["%s%d%s" % (aa, pos, r_[0])] + [r_[idx[c]] for c in cols] for r_ in table],
           "observed": ["%s%d%s" % (aa, pos, m) for m in sorted(x for x in obs if x)],
           "link": "https://www.ebi.ac.uk/ProtVar/query?search=%s+%s%d" % (uid, aa, pos)}
    if rows or table:
        S.append({"title": "Prediction", "color": "#6a1b9a", "rows": rows, "table": sub})
    pockets = []
    pk = (pv or {}).get("pocket")
    if pk:
        rows = []
        for i, p in enumerate(pk.get("pockets") or []):
            resid = [int(r) for r in p.get("resid") or []]
            stats = ["%d res" % len(resid)]
            for k_, fmt in (("meanPlddt", "pLDDT %d"), ("radGyration", "Rg %.1f"), ("energyPerVol", "E/V %.2f"),
                            ("score", "score %d"), ("buriedness", "buried %.2f")):
                if _finite(p.get(k_)):
                    stats.append(fmt % (round(p[k_]) if "%d" in fmt else p[k_]))
            ev = pocket_evidence(resid, geom, ann.get("sites") or [], seq) if resid else None
            pid = p.get("pocketId") if p.get("pocketId") is not None else i + 1
            motif = bool(vobj) and len(resid) >= 2 and _structure_motif_ok(vobj)
            pockets.append({"pocketId": pid, "resid": resid, "stats": " - ".join(stats), "motif": motif,
                            "composition": composition_text(ev["composition"]) if ev and ev["composition"] else "",
                            "axes": composition_axes(ev["composition"]) if ev and ev["composition"] else [],
                            "sites": ["%s (%d)" % (h.get("description") or "site", h["position"]) for h in
                                      (ev["siteHits"] if ev else [])[:4]]})
            rows.append(row("Pocket %s" % pid, " - ".join(stats), "#00897b", pocket=str(pid)))
        sec("AutoPocket", "#00897b", rows)
    lt, nl, it, ni = _kb_residue(kb, pos)
    rows = []
    if lt:
        rows.append(row("Ligands", ", ".join("%s (%d PDB)" % (l["code"], l["pdbCount"]) for l in lt)
                        + ((" +%d" % (nl - 8)) if nl > 8 else ""), "#2e7d32"))
    if it:
        rows.append(row("PPI interface", "; ".join("%s (%d)" % (p.get("name") or p.get("accession") or "?", p["pdbCount"])
                                                   for p in it) + ((" +%d" % (ni - 4)) if ni > 4 else ""), "#1565c0"))
    sec("Experimental binding (PDBe-KB)", "#2e7d32", rows)
    rows = []
    if near:
        rows.append(row("Residues", ", ".join("%s%d" % (seq[p - 1] if 0 < p <= len(seq) else "", p) for p in near),
                        positions=near))
    if plist:
        rows.append(row("other chains", ", ".join("%s:%s" % (ch, r) for ch, r, _i in plist[:40])
                        + ((" +%d" % (len(plist) - 40)) if len(plist) > 40 else "")))
    sec("Nearby (%g A)" % nb["dist"], "#455a64", rows)
    if rep["nearbyLigands"]:
        sec("Nearby ligands (<= 5 A)", "#6d4c41", [row("Ligands", ", ".join(rep["nearbyLigands"]))])
    rows = []
    if prox_ptms:
        rows.append(row("PTMs <= %g A" % pr, ", ".join("%d (%.1f)" % (n["pos"], n["dist"]) for n in prox_ptms[:24]),
                        positions=[n["pos"] for n in prox_ptms[:24]]))
    if prox_vars:
        rows.append(row("Pathogenic <= %g A" % vr, ", ".join("%s (%.1f)" % (n["label"], n["dist"]) for n in prox_vars[:24]),
                        positions=[n["pos"] for n in prox_vars[:24]]))
    sec("PTM-variant proximity (residue)", "#37474f", rows)
    # console text: format_report + the services lines of the CLI report
    if pv:
        bits = []
        if sc.get("conservation") is not None:
            bits.append("conservation %.2f" % sc["conservation"])
        if sc.get("m3d"):
            bits.append("M3D %s" % sc["m3d"].get("prediction", "?"))
        if bits:
            text.append("  ProtVar: " + " · ".join(bits))
        for p in pockets:
            text.append("  Pocket %s: %s" % (p["pocketId"], p["stats"]))
        for r_ in sub["rows"]:
            text.append("    %s: %s" % (r_[0], " · ".join("%s %s" % (c, ("%.2f" % v) if isinstance(v, float) else v)
                                                       for c, v in zip(cols, r_[1:]) if v is not None)))
    if lt:
        text.append("  PDBe-KB ligands: " + ", ".join("%s (%d PDB)" % (x["code"], x["pdbCount"]) for x in lt))
    if it:
        text.append("  PDBe-KB interfaces: " + "; ".join(x.get("name") or x.get("accession") or "?" for x in it))
    return {"uid": rep["uid"], "pos": pos, "aa": aa, "title": "%s%d" % (aa, pos),
            "flags": [{"key": k, "label": l, "color": c, "on": on} for k, l, c, on in flags],
            "sections": S, "pockets": pockets, "near": near,
            "partners": [{"chain": ch, "resi": r, "accession": (i or {}).get("accession") or "",
                          "uniPos": (i or {}).get("uniPos")} for ch, r, i in plist],
            "nearbyLigands": rep["nearbyLigands"], "services": services and pv is not None,
            "proxPtms": [[n["pos"], round(n["dist"], 2)] for n in prox_ptms],
            "proxVariants": [[n["label"], n["pos"], round(n["dist"], 2)] for n in prox_vars],
            "proximity": {k: v for k, v in (proximity or {}).items() if not isinstance(v, (list, dict, tuple))},
            "text": "\n".join(t for t in text if t)}


@_serve_op("report")
def _op_report(a):
    return _vmd_report(a)


@_serve_op("cif_auth")
def _op_cif_auth(a):
    """Author chain / residue numbers of an mmCIF model's atoms in file order (first model): VMD's pdbx reader
    takes label_seq_id ('.' -> 0 for ligands; SWISS-MODEL label numbering starts at 1) and decorates multi-letter
    asym ids, while UniProt numbering and the PyMOL plugin use auth_asym_id / auth_seq_id. -> {n, chain [..],
    resid [..], ins [..], name [..]} (the plugin checks n and the atom names before applying it)."""
    path = str(a.get("path") or "")
    cols, rows, inloop, model0 = [], [], False, None
    tok = re.compile(r"'[^']*'(?=\s|$)|\"[^\"]*\"(?=\s|$)|\S+")
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split()[0][11:])
                inloop = True
                continue
            if not inloop or not cols:
                continue
            if not s or s.startswith(("#", "loop_", "_", "data_")):
                if rows:
                    break
                continue
            if not s.startswith(("ATOM", "HETATM")):
                continue
            t = [x[1:-1] if x[:1] in "'\"" else x for x in tok.findall(s)]
            if len(t) != len(cols):
                continue
            rows.append(dict(zip(cols, t)))
    if not rows:
        return {"n": 0}
    model0 = rows[0].get("pdbx_PDB_model_num")
    rows = [r for r in rows if r.get("pdbx_PDB_model_num") == model0]

    def num(r):
        for k in ("auth_seq_id", "label_seq_id"):
            v = r.get(k)
            if v not in (None, "", ".", "?"):
                try:
                    return int(v)
                except ValueError:
                    pass
        return 0
    return {"n": len(rows),
            "chain": [r.get("auth_asym_id") if r.get("auth_asym_id") not in (None, "", ".", "?") else r.get("label_asym_id", "")
                      for r in rows],
            "resid": [num(r) for r in rows],
            "ins": [("" if r.get("pdbx_PDB_ins_code") in (None, "", ".", "?") else r["pdbx_PDB_ins_code"]) for r in rows],
            "name": [r.get("label_atom_id") or r.get("auth_atom_id") or "" for r in rows]}


@_serve_op("bfactor_values")
def _op_bfactor_values(a):
    """ufv_save_bfactor values (export.js rewritePdbBeta / the PyMOL plugin's _bfactor_values), rounded to 2
    decimals: {keep 1} for plddt / bfactor (file values kept), else {keep 0, values {pos: v}}. Modes:
    alphamissense, hotspots (1 / 0.67 / 0.33), hubs (1 / 0.5), burden (1), pocket (1 buried / 0.66 exposed at q),
    selected (1 = `positions`, the residues of the plugin's active marker layers)."""
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    mode = str(a.get("mode") or "selected").strip().lower()
    mode = {"am": "alphamissense", "contacthubs": "hubs", "pockets": "pocket", "prism": "pocket"}.get(mode, mode)
    if mode in ("plddt", "bfactor"):
        return {"keep": 1, "values": {}, "mode": mode}
    vals = {}
    if mode == "alphamissense":
        vals = {int(p): float(v) for p, v in (ann.get("amMean") or {}).items()}
    elif mode == "burden":
        vals = {int(p): 1.0 for p in ann.get("burden") or ()}
    elif mode in ("hotspots", "hubs", "pocket"):
        vobj = _vmd_vobj(a, uid)
        geom, gh = _vmd_geom(a)
        if mode == "hotspots":
            vals = {p: {"strong": 1.0, "moderate": 0.67, "weak": 0.33}.get(t, 0.0)
                    for p, t in _vmd_analysis(vobj, uid, geom, gh, "hotspots", ann)["tiers"].items()}
        elif mode == "hubs":
            vals = {p: {"strong": 1.0, "moderate": 0.5}.get(t, 0.0)
                    for p, t in _vmd_analysis(vobj, uid, geom, gh, "hubs", ann)["tiers"].items()}
        else:
            q = float(a.get("q") or POCKET_SETTINGS["q"])
            vals = {p: {"pocket": 1.0, "exposed": 0.66}.get(c, 1.0) for p, c in
                    _pocket_classes(_vmd_analysis(vobj, uid, geom, gh, "pockets", ann), q).items()}
    elif mode == "selected":
        pos = a.get("positions") or []
        if isinstance(pos, str):
            pos = [p for p in re.split(r"[\s,+]+", pos) if p]
        vals = {int(p): 1.0 for p in pos}
    else:
        raise ValueError("unknown B-factor mode '%s'" % mode)
    return {"keep": 0, "values": {int(p): round(float(v), 2) for p, v in vals.items()}, "mode": mode}


@_serve_op("csv_full")
def _op_csv_full(a):
    """ufv_csv for a VMD molecule: the PyMOL plugin's per-residue matrix with the structure columns (pdb_residue
    per mapped chain, hotspot / hub tiers, constraint-pocket q / class, PTM-variant proximity at the radii,
    nearby_ligands from the plugin's VMD contact search) and optional ProtVar columns (protvar=1; progress key
    'protvar:<UID>')."""
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    out = str(a.get("out") or "")
    if not out:
        raise ValueError("missing 'out'")
    vobj = _vmd_vobj(a, uid)
    geom, gh = _vmd_geom(a)
    chains = [str(c) for c in (a.get("chains") or [])]
    hot = hub = pocket = prox = hot_ch = hub_ch = None
    if geom:
        try:
            r = _vmd_analysis(vobj, uid, geom, gh, "hotspots", ann)
            hot, hot_ch = r["tiers"], r["byChain"]
        except Exception:
            pass
        try:
            r = _vmd_analysis(vobj, uid, geom, gh, "hubs", ann)
            hub, hub_ch = r["tiers"], r["byChain"]
        except Exception:
            pass
        if ann.get("amMap"):
            try:
                pocket = (_vmd_analysis(vobj, uid, geom, gh, "pockets", ann) or {}).get("byPos") or {}
            except Exception:
                pass
        if ann["ptms"] and ann["variants"]:
            try:
                prox = _ptm_variant_proximity(ann["ptms"], ann["variants"], geom, *_vmd_radii(a))
            except Exception:
                pass
    lig = {}
    for k, v in (a.get("lig_contacts") or {}).items():
        names = [x for x in (re.split(r"[;, ]+", v) if isinstance(v, str) else v) if x]
        if names:
            lig[int(k)] = sorted(set(names))
    pv = None
    if str(a.get("protvar", 0)) not in ("0", "", "false", "False"):
        key = "protvar:" + uid
        _VMD_PROGRESS[key] = [0, len(ann["sequence"])]
        pv = _fetch_protvar_all(ann["uid"], len(ann["sequence"]), lambda d, t: _VMD_PROGRESS.__setitem__(key, [d, t]))
    text = _residue_csv_text(ann, vobj if geom else None, pv, chains, hot, hub, pocket, prox, hot_ch, hub_ch, lig)
    with open(out, "w") as fh:
        fh.write(text)
    return {"path": out, "rows": max(0, text.count("\n")), "chains": chains}


# ---- trajectory analyses for the VMD plugin (vmd/ufv_traj.tcl) -----------------------------------------------
# VMD writes the atom sets of the analysis window as DCD files (`animate write dcd ... sel`), the plugin sends
# the paths + atom bookkeeping, and the SAME numpy maths as the PyMOL plugin (_traj_calc: _calc_rmsf / _dist /
# _pockets / _contacts / _annot) runs here, so both plugins give identical numbers on the same trajectory.
_VMD_TRAJ_LAST = {}     # (key, kind) -> (signature, result): last contact runs (annotation dynamics reuses them)


def _read_dcd(path):
    """CHARMM / NAMD / VMD DCD (little endian, optional unit-cell record, no fixed atoms) -> float32
    (frames, atoms, 3). One read + strided numpy views (the plugin writes a fresh file per request)."""
    np = _numpy()
    path = str(path)
    with open(path, "rb") as fh:
        b = fh.read()

    def i4(o):
        return int.from_bytes(b[o:o + 4], "little", signed=True)
    n0 = i4(0)
    if n0 != 84 or b[4:8] != b"CORD":
        raise ValueError("%s: not a little-endian CHARMM/NAMD DCD file" % path)
    icntrl = [i4(8 + 4 * k) for k in range(20)]
    if icntrl[8]:
        raise ValueError("%s: DCD files with fixed atoms are not supported" % path)
    cell = bool(icntrl[19] and icntrl[10])
    four = bool(icntrl[19] and icntrl[11])
    p = 4 + n0 + 4
    p += 4 + i4(p) + 4                                  # title record
    natoms = i4(p + 4)
    p += 4 + i4(p) + 4
    rec = 4 * natoms + 8
    fsize = (56 if cell else 0) + 3 * rec + (rec if four else 0)
    nf = (len(b) - p) // fsize if fsize else 0
    A = np.frombuffer(b, np.uint8, nf * fsize, p).reshape(nf, fsize)
    off = 56 if cell else 0
    X = np.empty((nf, natoms, 3), np.float32)
    for k in range(3):
        s = off + k * rec + 4
        X[:, :, k] = np.ascontiguousarray(A[:, s:s + 4 * natoms]).view("<f4")
    return X


def _vmd_cols(X, cols, f32=False):
    np = _numpy()
    out = X if cols is None else X[:, np.asarray(cols, dtype=int)]
    return out.astype(np.float32 if f32 else np.float64, copy=True)


def _vmd_set(a, dcd_key, ref_key, cols, frames, ref, f32=False):
    """(window array, reference frame) of one atom set: the window DCD (one frame per window frame) and the
    reference frame from it, or from a one-frame DCD when the reference lies outside the window."""
    X = _read_dcd(a[dcd_key])
    if X.shape[0] != len(frames):
        raise ValueError("%s has %d frames, the window has %d" % (a[dcd_key], X.shape[0], len(frames)))
    W = _vmd_cols(X, cols, f32)
    if a.get(ref_key):
        R = _vmd_cols(_read_dcd(a[ref_key])[:1], cols, f32)[0]
    else:
        R = W[frames.index(ref)].copy()
    return W, R


def _vmd_traj_groups(uid, a, with_lig=True):
    """_traj_groups (PyMOL glue) without cmd: hotspots / constraint-pocket classes and the ligand pocket come
    from the plugin (its per-frame tier cache / a VMD `within` query)."""
    ann = (_CACHE.get(uid) or {}) if uid else {}
    if uid and not ann:
        try:
            ann = _quiet(fetch_annotations, uid)
        except Exception:
            ann = {}
    g = {}
    g["PTM sites"] = {int(p["position"]) for p in ann.get("ptms") or [] if p.get("position")}
    g["Pathogenic variants"] = {int(v["position"]) for v in ann.get("variants") or []
                                if v.get("position") and _is_pathogenic(v)}
    s = set()
    for f in ann.get("sites") or []:
        x, y = int(f.get("position") or 0), int(f.get("endPosition") or f.get("position") or 0)
        if x and y - x <= 20:
            s.update(range(x, y + 1))
    g["Functional sites"] = s
    if uid:
        if a.get("hot"):
            g["Hotspots"] = {int(p) for p in a["hot"]}
        if a.get("cp"):
            g["Constraint pocket"] = {int(p) for p, _c in a["cp"]}
        acc = _canon_acc(uid)
        pv = [v for k, v in _ALL_POCKETS.items() if k[0] == acc and v]
        if pv:
            g["ProtVar pocket"] = {int(r) for pk in pv[0]["pockets"][:3] for r in pk.get("resid") or []}
    if with_lig and a.get("ligpocket"):
        g["Ligand pocket"] = {int(p) for p in a["ligpocket"]}
    return {k: v for k, v in g.items() if v}


def _vmd_traj_tags(uid, a):
    short = {"PTM sites": "PTM", "Pathogenic variants": "Path", "Functional sites": "Site", "Hotspots": "Hot",
             "Constraint pocket": "CP", "ProtVar pocket": "PV"}
    out = {}
    for g, ps in _vmd_traj_groups(uid, a, with_lig=False).items():
        for p in ps:
            out.setdefault(p, []).append(short.get(g, g))
    return {p: ",".join(v) for p, v in out.items()}


def _vmd_traj_base(a):
    """The common half of the prep dict (_traj_prep): window, reference, deduplicated Cα, annotation tags."""
    if _numpy() is None:
        raise RuntimeError("trajectory analyses need numpy")
    uid = str(a.get("uid") or "").strip().upper()
    frames = [int(f) for f in a.get("frames") or []]
    if not frames:
        raise ValueError("empty frame window")
    ref = int(a.get("ref", frames[0]))
    p = {"kind": a.get("kind"), "obj": str(a.get("key") or ""), "uid": uid, "frames": frames, "ref": ref,
         "t0": time.time(), "cutoff": float(a.get("cutoff") or 4.5), "threshold": float(a.get("threshold") or 2.0),
         "spec": a.get("spec")}
    p["res"] = [(str(r[0]), str(r[1]), int(r[2]), str(r[3]), None if r[4] in (None, "") else int(r[4]))
                for r in a.get("res") or []]
    if not p["res"]:
        raise RuntimeError("no polymer CA atoms")
    p["CA"], p["CAref"] = _vmd_set(a, "ca_dcd", "ca_ref", a.get("ca_cols"), frames, ref)
    p["tags"] = _vmd_traj_tags(uid, a)
    return p


def _vmd_contact_prep(p, q, kind):
    """Contact atom sets from the plugin (rows = [chain, resi, resv, resn, name], the window DCD of each set)."""
    q = q or {}
    if q.get("error"):
        raise RuntimeError(q["error"])
    same = bool(q.get("same"))
    p["same"] = same
    p["A"], p["Aref"] = _vmd_set(q, "dcdA", "refA", None, p["frames"], p["ref"], f32=True)
    if same:
        p["B"], p["Bref"] = p["A"], p["Aref"]
    else:
        p["B"], p["Bref"] = _vmd_set(q, "dcdB", "refB", None, p["frames"], p["ref"], f32=True)

    def rows(key):
        return [(str(r[0]), str(r[1]), int(r[2]), str(r[3]), str(r[4])) for r in q.get(key) or []]
    p["rowsA"] = rows("rowsA")
    p["rowsB"] = p["rowsA"] if same else rows("rowsB")
    uni = lambda l: [None if u in (None, "") else int(u) for u in l or []]
    p["uniA"] = uni(q.get("uniA"))
    p["uniB"] = p["uniA"] if same else uni(q.get("uniB"))
    p["ours"] = set(str(c) for c in q["ours"]) if q.get("ours") is not None else None
    if q.get("center") not in (None, ""):
        p["center"] = str(q["center"])
    if q.get("cutoff") not in (None, ""):
        p["cutoff"] = float(q["cutoff"])
    p["kind"] = kind
    return p


def _vmd_prep_pockets(p, a):
    """_prep_pockets without cmd: cp = constraint-pocket classes [[pos, class]..] from the plugin's tier cache,
    pv = fetched ProtVar pockets, lig = ligpockets [{label, keys [[chain, resv]..]}] (VMD query in the
    reference frame), 120+124+156 = UniProt positions."""
    col_of = {}
    for i, (ch, resi, resv, resn, u) in enumerate(p["res"]):
        if u is not None:
            col_of.setdefault(u, i)
    key_col = {(r[0], r[2]): i for i, r in enumerate(p["res"])}
    uid = p["uid"]
    cls = {int(u): str(c) for u, c in a.get("cp") or []}
    pockets = []
    for part in re.split(r"[;\n]+", (p.get("spec") or "auto").strip()):
        part = part.strip()
        pl = part.lower()
        if not part:
            continue
        if pl in ("auto", "cp"):
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
        if pl in ("auto", "lig"):
            for lp in a.get("ligpockets") or []:
                pockets.append((str(lp["label"]), [("k", str(k[0]), int(k[1])) for k in lp.get("keys") or []]))
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


def _jsonable(o):
    np = _numpy()
    if isinstance(o, dict):
        return {(k if isinstance(k, str) else str(k)): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_jsonable(v) for v in o]
    if np is not None and isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if np is not None and isinstance(o, np.generic):
        return o.item()
    return o


def _vmd_traj_reply(r):
    """Result dict -> JSON-able reply (series as float lists; per-kind extras for colouring / display / CSV)."""
    hd, rows = r["table"]
    out = {"kind": r["kind"], "frames": r["frames"], "ref": r["ref"], "summary": r["summary"],
           "series": [[n, [float(v) for v in s]] for n, s in r.get("series") or []],
           "header": list(hd), "rows": _jsonable(rows), "rowmeta": _jsonable(r.get("rowmeta") or []),
           "prep_ms": round(r.get("prep_ms", 0.0), 1), "calc_ms": round(r.get("calc_ms", 0.0), 1)}
    k = r["kind"]
    if k == "rmsf":
        out["rmsd"] = [float(v) for v in r["rmsd"]]
        out["rmsf"] = [float(v) for v in r["rmsf"]]
    elif k == "dist":
        out["pairs"] = _jsonable(r["pairs"])
        out["threshold"] = r["threshold"]
    elif k == "pockets":
        out["persist"] = {lab: [[int(c), float(f)] for c, f in d.items()] for lab, d in r["persist"].items()}
    elif k in ("ligand", "partner", "nearby"):
        out["occupancy"] = [[ch, int(rv), float(o)] for (ch, rv), o in r["occupancy"].items()]
        out["cutoff"] = r["cutoff"]
    return out


def _vmd_traj(a):
    p = _vmd_traj_base(a)
    kind = p["kind"]
    uid = p["uid"]
    if kind == "rmsf":
        p["groups"] = _vmd_traj_groups(uid, a)
    elif kind == "dist":
        Xu = _vmd_cols(_read_dcd(a["dcd"]), None)        # pair atoms: window frames only (no reference needed)
        if Xu.shape[0] != len(p["frames"]):
            raise ValueError("%s has %d frames, the window has %d" % (a["dcd"], Xu.shape[0], len(p["frames"])))
        pairs = []
        for pr in a.get("pairs") or []:
            pairs.append({"label": str(pr["label"]), "tok": str(pr.get("tok") or pr["label"]), "mode": str(pr["mode"]),
                          "A": Xu[:, [int(c) for c in pr["A"]]], "B": Xu[:, [int(c) for c in pr["B"]]],
                          "idxA": [int(x) for x in pr["idxA"]], "idxB": [int(x) for x in pr["idxB"]]})
        if not pairs:
            raise RuntimeError("no distance pairs resolved (spec: %s)" % p.get("spec"))
        p["pairs"] = pairs
    elif kind == "pockets":
        _vmd_prep_pockets(p, a)
    elif kind in ("ligand", "partner", "nearby"):
        _vmd_contact_prep(p, a, kind)
    elif kind == "annot":
        for k in ("ligand", "partner"):
            sig = (tuple(p["frames"]), p["ref"], p["cutoff"])
            last = _VMD_TRAJ_LAST.get((p["obj"], k))
            if last and last[0] == sig and not (a.get(k) or {}).get("dcdA"):
                p[k] = {"reuse": last[1]}
                continue
            try:
                p[k] = _vmd_contact_prep(dict(p), a.get(k), k)
            except Exception as e:
                p[k] = {"error": str(e)}
    else:
        raise ValueError("unknown trajectory analysis '%s'" % kind)
    p["prep_ms"] = (time.time() - p["t0"]) * 1000
    r = _traj_calc(p)
    if kind in ("ligand", "partner"):
        _VMD_TRAJ_LAST[(p["obj"], kind)] = ((tuple(p["frames"]), p["ref"], p["cutoff"]), r)
    elif kind == "annot":
        for k in ("ligand", "partner"):
            if r["sub"].get(k) is not None:
                _VMD_TRAJ_LAST[(p["obj"], k)] = ((tuple(p["frames"]), p["ref"], p["cutoff"]), r["sub"][k])
    return _vmd_traj_reply(r)


@_serve_op("traj")
def _op_traj(a):
    """One trajectory analysis (kind rmsf | dist | pockets | ligand | partner | nearby | annot) on DCD files
    written by the VMD plugin; same maths and output as the PyMOL plugin's ufv_traj_* commands."""
    return _vmd_traj(a)


@_serve_op("traj_pad")
def _op_traj_pad(a):
    """Contact candidate pad (_prep_contacts): 2 x the 95th-percentile Cα displacement after fitting to the
    reference + 6 Å (<= 25 Å) -> {pad}."""
    np = _numpy()
    p = _vmd_traj_base(dict(a, kind="pad"))
    Y, _T = _fit(p["CA"], p["CAref"])
    ri = _ref_index(p)
    Yr = Y[ri] if ri is not None else p["CAref"]
    dmax = float(np.percentile(np.sqrt(((Y - Yr) ** 2).sum(-1)).max(0), 95)) if len(Y) > 1 else 0.0
    return {"pad": min(25.0, 2.0 * dmax + 6.0)}


@_serve_op("traj_ctx")
def _op_traj_ctx(a):
    """Annotation context for trajectory summaries: {groups {label [pos..]}, tags {pos tag}}."""
    uid = str(a.get("uid") or "").strip().upper()
    g = _vmd_traj_groups(uid, a)
    return {"groups": {k: sorted(v) for k, v in g.items()}, "tags": {str(k): v for k, v in _vmd_traj_tags(uid, a).items()}}


@_serve_op("echo")
def _op_echo(a):
    """Test helper: the arguments back (round-trip check of the request / reply encodings)."""
    return a


@_serve_op("sleep")
def _op_sleep(a):
    """Test helper: wait `ms` milliseconds (shows that other requests are not blocked)."""
    time.sleep(float(a.get("ms") or 0) / 1000.0)
    return {"slept": float(a.get("ms") or 0)}


# ---- annotation search (v3) ----
# ONE annotation list for both panels (PyMOL Qt / VMD Tk): rows per category with stable ids, and ONE
# matcher, so a query selects exactly the same rows in both plugins. Pure functions (no PyMOL / Tk, no
# network): ligand / partner / PDBe-KB rows are built from data the caller passes in `extra`.
#   rows = ann_rows(ann, cat, extra)           cat in ANN_CATS, or "all" (every category, ANN_CATS order)
#   ids  = ann_search(rows, query)             matching row ids, in row order
#   sugg = ann_suggest(rows, partial)          completions for the last (partial) term of a query
#   q    = ann_filter_query(tokens, provenances, diseases, reviewed_only)   old variant filters -> query
# Query syntax: space-separated terms are ANDed. A term is free text (case-insensitive substring of the row
# text: category, consequence, provenance, disease names / ids, description, PTM type, ligand code ...), a
# position (1378) or range (100-200), a mutation (N1378H, R2K, N1378), a "quoted phrase", or a facet:
#   cons:pathogenic|deleterious|benign|uncertain   src:uniprot|clinvar|large (src:reviewed = reviewed only)
#   dis:"Label"|ID      cat:ptm|var|...|<PTM category>      pos:100-200
# '|' separates alternatives inside one term (OR). A facet value 'none' matches nothing (the old filters'
# empty sets). Serve ops: ann_rows, search, filter_query.
ANN_CATS = ("ptm", "var", "site", "mutag", "dom", "lig", "prt", "kb")
ANN_CAT_LABELS = {"ptm": "PTMs", "var": "Variants", "site": "Sites", "mutag": "Mutagenesis", "dom": "Domains",
                  "lig": "Ligands", "prt": "Partners", "kb": "PDBe-KB"}
ANN_CAT_SHORT = {"ptm": "PTM", "var": "Var", "site": "Site", "mutag": "Mut", "dom": "Dom", "lig": "Lig",
                 "prt": "Prt", "kb": "KB"}
_ANN_CAT_ALIAS = {"ptm": "ptm", "ptms": "ptm", "var": "var", "variant": "var", "variants": "var", "site": "site",
                  "sites": "site", "mutag": "mutag", "mutagenesis": "mutag", "dom": "dom", "domain": "dom",
                  "domains": "dom", "lig": "lig", "ligand": "lig", "ligands": "lig", "prt": "prt", "partner": "prt",
                  "partners": "prt", "kb": "kb", "pdbekb": "kb", "pdbe-kb": "kb", "binding": "kb"}
ANN_CONS_TOKENS = ("pathogenic", "deleterious", "benign", "uncertain")
ANN_SRC_TOKENS = ("uniprot", "clinvar", "large")
_ANN_SRC_OF = {"UniProt reviewed": "uniprot", "ClinVar": "clinvar", "Large scale studies": "large"}
_ANN_SRC_CAT = {v: k for k, v in _ANN_SRC_OF.items()}
_ANN_SRC_ALIAS = {"uniprot": "uniprot", "reviewed": "uniprot", "clinvar": "clinvar", "large": "large",
                  "large-scale": "large", "largescale": "large", "large_scale": "large"}
_ANN_FACETS = {"cons": "cons", "consequence": "cons", "src": "src", "source": "src", "prov": "src",
               "provenance": "src", "dis": "dis", "disease": "dis", "cat": "cat", "category": "cat",
               "type": "cat", "pos": "pos", "position": "pos"}
_ANN_ROWS_MEMO = {}          # (id(ann), cat) -> (ann, rows)   (ann-derived categories only)
_RE_ANN_NUM = re.compile(r"^\d+$")
_RE_ANN_RANGE = re.compile(r"^(\d+)\s*-\s*(\d+)$")
_RE_ANN_MUT = re.compile(r"^([A-Za-z])(\d+)([A-Za-z*]?)$")


def ann_category(name):
    """'variants' / 'Var' / 'PDBe-KB' ... -> category key (ANN_CATS) or None."""
    s = str(name or "").strip().lower()
    if s in _ANN_CAT_ALIAS:
        return _ANN_CAT_ALIAS[s]
    for k, lbl in ANN_CAT_LABELS.items():
        if lbl.lower() == s:
            return k
    return None


def _ann_row(rid, cat, pos, end, label, color, hay, pts, poslabel=None, **facets):
    return {"id": rid, "cat": cat, "pos": pos, "end": end if end is not None else pos, "label": label,
            "color": color, "poslabel": poslabel if poslabel is not None else ("" if pos is None else str(pos)),
            "pts": pts, "hay": hay.lower(), "f": facets}


def _ann_build(ann, cat, extra=None):
    rows = []
    extra = extra or {}
    if cat == "ptm":
        for i, p in enumerate(ann.get("ptms") or []):
            pos, end = p["position"], p.get("endPosition") or p["position"]
            pts = [pos] + ([end] if end != pos else [])
            rows.append(_ann_row("ptm:%d" % i, cat, pos, end, "%s — %s" % (p.get("category"), p.get("description")),
                                 p.get("color"), "ptm ptms %s %s %s" % (p.get("category") or "", p.get("description") or "",
                                                                         p.get("type") or ""),
                                 pts, ("%d-%d" % (pos, end)) if end != pos else None,
                                 ptmcat=(p.get("category") or "").lower(), ptype=(p.get("type") or "").lower()))
    elif cat == "var":
        for i, v in enumerate(ann.get("variants") or []):
            pos = v["position"]
            wt, mt = v.get("wildType") or "", v.get("mutant") or ""
            dis = v.get("diseases") or []
            cons = v.get("consequence") or ""
            src = _ANN_SRC_OF.get(v.get("provenance"))
            label = "%s%d%s · %s%s" % (wt, pos, mt, cons, (" · " + "; ".join(dis)) if dis else "")
            hay = " ".join(["variant variants", "%s%d%s" % (wt, pos, mt), cons, v.get("provenance") or "",
                            "large-scale" if src == "large" else "", " ".join(dis),
                            " ".join(v.get("diseaseIds") or []), v.get("description") or "", v.get("clinVar") or "",
                            " ".join(v.get("rsIds") or [])])
            rows.append(_ann_row("var:%d" % i, cat, pos, pos, label, v.get("consequenceColor"), hay, [pos],
                                 cons=_cons_token(cons) if cons in CONSEQUENCE_COLORS else None, src=src,
                                 dis=frozenset(x.lower() for x in list(dis) + list(v.get("diseaseIds") or [])),
                                 wt=wt.upper(), muts=frozenset([mt.upper()])))
    elif cat == "site":
        for i, s in enumerate(ann.get("sites") or []):
            pos, end = s["position"], s.get("endPosition") or s["position"]
            rows.append(_ann_row("site:%d" % i, cat, pos, end, s.get("description") or "Site", s.get("color") or SITE_COLOR,
                                 "site sites " + (s.get("description") or ""), [pos] + ([end] if end != pos else []),
                                 ("%d-%d" % (pos, end)) if end != pos else None))
    elif cat == "mutag":
        for i, m in enumerate(ann.get("mutagenesis") or []):
            pos, end = m["position"], m.get("endPosition") or m["position"]
            sub = "%s→%s" % (m.get("wildType", ""), ", ".join(m.get("mutants", [])))
            rows.append(_ann_row("mutag:%d" % i, cat, pos, end, "%s · %s" % (sub, m.get("effect", "")),
                                 m.get("color") or MUTAG_COLOR,
                                 "mutagenesis %s %s %s" % (sub, " ".join(m.get("mutants", [])), m.get("effect", "")),
                                 list(range(pos, end + 1)), ("%d-%d" % (pos, end)) if end != pos else None,
                                 wt=(m.get("wildType") or "").upper(),
                                 muts=frozenset(x.upper() for x in m.get("mutants", []))))
    elif cat == "dom":
        for i, d in enumerate(ann.get("domains") or []):
            pos, end = d["position"], d.get("endPosition") or d["position"]
            rng = d.get("isRange") and end != pos
            rows.append(_ann_row("dom:%d" % i, cat, pos, end, d.get("description") or d.get("type") or "Domain",
                                 d.get("color"), "domain domains %s %s" % (d.get("type") or "", d.get("description") or ""),
                                 list(range(pos, end + 1)) if d.get("isRange") else [pos],
                                 ("%d-%d" % (pos, end)) if rng else None))
    elif cat == "lig":
        for L in extra.get("lig") or []:
            resn, cnt = str(L.get("resn")), int(L.get("count") or 1)
            nm = L.get("name") or ""
            rows.append(_ann_row("lig:%s" % resn, cat, None, None, "%s%s%s" % (
                resn, (" ×%d" % cnt) if cnt > 1 else "", (" — " + nm) if nm else ""), L.get("color") or "#2e7d32",
                "ligand ligands %s %s" % (resn, nm), [], resn, resn=resn.upper()))
    elif cat == "prt":
        for P in extra.get("prt") or []:
            ch, resi = str(P.get("chain")), str(P.get("resi"))
            rows.append(_ann_row("prt:%s:%s" % (ch, resi), cat, None, None, P.get("text") or "", P.get("color"),
                                 "partner partners " + (P.get("text") or ""), [], "%s:%s" % (ch, resi),
                                 chain=ch, resi=resi))
    elif cat == "kb":
        kb = extra.get("kb") or {}
        lbp = {int(k): v for k, v in (kb.get("ligandByPos") or {}).items()}
        ibp = {int(k): v for k, v in (kb.get("interfaceByPos") or {}).items()}
        cols = {"ligand": "#2e7d32", "interface": "#1565c0", "both": "#00838f"}      # == ufv_pymol KB_COLORS
        for p in sorted(set(lbp) | set(ibp)):
            ls = sorted(lbp.get(p) or [], key=lambda x: -(x.get("pdbCount") or 0))
            its = sorted(ibp.get(p) or [], key=lambda x: -(x.get("pdbCount") or 0))
            bits = []
            if ls:
                bits.append("ligand: " + ", ".join(str(x.get("code")) for x in ls[:6]) + (" +%d" % (len(ls) - 6) if len(ls) > 6 else ""))
            if its:
                bits.append("interface: " + ", ".join(str(x.get("name") or x.get("accession")) for x in its[:3])
                            + (" +%d" % (len(its) - 3) if len(its) > 3 else ""))
            kind = "both" if (ls and its) else "ligand" if ls else "interface"
            rows.append(_ann_row("kb:%d" % p, cat, p, p, " · ".join(bits), cols[kind],
                                 "pdbe-kb binding %s %s %s" % (" ".join(bits), " ".join(str(x.get("name") or "") for x in ls),
                                                               " ".join(str(x.get("accession") or "") for x in its)),
                                 [p], kb=kind))
    return rows


def ann_rows(ann, cat, extra=None):
    """Rows of one category (or 'all'). Each row: {id, cat, pos, end, poslabel, label, color, pts (UniProt
    positions drawn in 3D), hay (lower-case search text), f (facets)}. Ids are stable for one annotation set:
    '<cat>:<index in the annotation list>' (ptm / var / site / mutag / dom), 'lig:<RESN>', 'prt:<chain>:<resi>',
    'kb:<pos>'. ann-derived categories are memoised per annotation dict."""
    if cat in (None, "", "all"):
        out = []
        for c in ANN_CATS:
            out.extend(ann_rows(ann, c, extra))
        return out
    cat = ann_category(cat) or cat
    if cat in ("lig", "prt", "kb"):
        return _ann_build(ann or {}, cat, extra)
    key = (id(ann), cat)
    hit = _ANN_ROWS_MEMO.get(key)
    if hit is not None and hit[0] is ann:
        return hit[1]
    rows = _ann_build(ann or {}, cat)
    if len(_ANN_ROWS_MEMO) > 64:
        _ANN_ROWS_MEMO.clear()
    _ANN_ROWS_MEMO[key] = (ann, rows)
    return rows


def _ann_lex(q):
    """Query -> [(raw term text, [(alt text, quote_at)])]: whitespace splits terms and '|' alternatives,
    both only outside double quotes; quote_at = index of the first quote in the alt (-1: unquoted)."""
    terms, alts, buf, qat, raw, inq = [], [], [], -1, [], False

    def flush():
        if buf or alts or qat >= 0:
            alts.append(("".join(buf), qat))
            terms.append(("".join(raw), alts))
    for ch in (q or ""):
        if ch == '"':
            inq = not inq
            if qat < 0:
                qat = len(buf)
            raw.append(ch)
            continue
        if not inq and ch.isspace():
            flush()
            alts, buf, qat, raw = [], [], -1, []
            continue
        raw.append(ch)
        if not inq and ch == "|":
            alts.append(("".join(buf), qat)); buf, qat = [], -1
            continue
        buf.append(ch)
    flush()                    # last term (also an unterminated "quote)
    return terms


def _ann_alt(text, qat, facet):
    """One alternative -> predicate(row) (or None for an empty alternative)."""
    i = text.find(":")
    if i > 0 and (qat < 0 or i < qat) and text[:i].lower() in _ANN_FACETS:
        facet, text = _ANN_FACETS[text[:i].lower()], text[i + 1:]
    v = text.strip()
    if not v and facet is None:
        return None, facet
    lv = v.lower()
    if facet == "cons":
        tok = next((t for t in ANN_CONS_TOKENS if t.startswith(lv)), None) if lv and lv != "none" else None
        if lv in ("vus",):
            tok = "uncertain"
        return (lambda r: tok is not None and r["f"].get("cons") == tok), facet
    if facet == "src":
        tok = _ANN_SRC_ALIAS.get(lv) or next((t for t in ANN_SRC_TOKENS if lv and t.startswith(lv)), None)
        return (lambda r: tok is not None and r["f"].get("src") == tok), facet
    if facet == "dis":
        return (lambda r: lv != "none" and lv in (r["f"].get("dis") or ())), facet
    if facet == "cat":
        c = ann_category(lv)
        if c:
            return (lambda r: r["cat"] == c), facet
        return (lambda r: bool(lv) and (r["f"].get("ptmcat") or "").startswith(lv)
                or bool(lv) and r["f"].get("ptype") == lv), facet
    if facet == "pos" or (qat < 0 and (_RE_ANN_NUM.match(v) or _RE_ANN_RANGE.match(v))):
        m = _RE_ANN_RANGE.match(v)
        if m:
            lo, hi = sorted((int(m.group(1)), int(m.group(2))))
        elif _RE_ANN_NUM.match(v):
            lo = hi = int(v)
        else:
            return (lambda r: False), facet
        return (lambda r: r["pos"] is not None and r["pos"] <= hi and r["end"] >= lo), facet
    m = _RE_ANN_MUT.match(v) if qat < 0 else None
    if m:
        wt, p, mt = m.group(1).upper(), int(m.group(2)), m.group(3).upper()

        def mut(r, wt=wt, p=p, mt=mt, lv=lv):
            f = r["f"]
            if "wt" in f and r["pos"] is not None and r["pos"] <= p <= r["end"] and f["wt"][:1] == wt \
                    and (not mt or mt in f.get("muts", ())):
                return True
            return lv in r["hay"]
        return mut, facet
    return (lambda r: lv in r["hay"]), facet


def _ann_compile(q):
    terms = []
    for _raw, alts in _ann_lex(q):
        preds, facet = [], None
        for text, qat in alts:
            p, facet = _ann_alt(text, qat, facet)
            if p is not None:
                preds.append(p)
        if preds:
            terms.append(preds)
    return terms


def ann_search(rows, query):
    """Ids of the rows matching `query` (see the section comment), in row order. Empty query = all rows."""
    terms = _ann_compile(query)
    if not terms:
        return [r["id"] for r in rows]
    if all(len(t) == 1 for t in terms):
        preds = [t[0] for t in terms]
        return [r["id"] for r in rows if all(p(r) for p in preds)]
    return [r["id"] for r in rows if all(any(p(r) for p in t) for t in terms)]


def ann_last_term(query):
    """(text before the last term, the last term) — the part a completion replaces ('' after a space)."""
    q = query or ""
    if not q or (q[-1].isspace() and q.count('"') % 2 == 0):
        return q, ""
    terms = _ann_lex(q)
    if not terms:
        return q, ""
    last = terms[-1][0]                 # the raw text of the last term is a suffix of q
    return (q[:len(q) - len(last)], last) if q.endswith(last) else (q, "")


def _ann_quote(s):
    s = str(s).replace('"', "")
    return '"%s"' % s if re.search(r'[\s|:]', s) else s


def _ann_vocab(rows):
    """Completion vocabulary of a row set: [(insert text, match key, weight)]."""
    cnt = {}

    def add(txt, w=1):
        cnt[txt] = cnt.get(txt, 0) + w
    cats = {r["cat"] for r in rows}
    if "var" in cats:
        for r in rows:
            if r["cat"] == "var":
                if r["f"].get("cons"):
                    add("cons:" + r["f"]["cons"])
                if r["f"].get("src"):
                    add("src:" + r["f"]["src"])
        for w in ("pathogenic", "likely", "benign", "uncertain", "deleterious", "uniprot", "clinvar", "large-scale",
                  "reviewed"):
            add(w, 0)
        add("src:reviewed", 0)
    seen_dis = {}
    for r in rows:
        if r["cat"] == "var":
            lab = r["label"].split(" · ", 2)
            for d in (lab[2].split("; ") if len(lab) > 2 else []):
                seen_dis[d] = seen_dis.get(d, 0) + 1
        elif r["cat"] == "ptm":
            if r["f"].get("ptmcat"):
                add(_ann_quote(r["f"]["ptmcat"]))
            if r["f"].get("ptype"):
                add(r["f"]["ptype"])
        elif r["cat"] == "lig":
            add(r["f"]["resn"].lower())
    for d, n in seen_dis.items():
        add("dis:" + _ann_quote(d), n)
    for c in cats:
        add("cat:" + c, 0)
    out = []
    for txt, n in cnt.items():
        out.append((txt, txt.lower().replace('"', ""), n))
    return out


_ANN_VOCAB_MEMO = {}


def ann_suggest(rows, partial, limit=12):
    """Completions for a partial term (consequence / source words and facets, disease names, PTM categories /
    types, ligand codes, cat: facets): prefix matches first (then word-start matches inside multi-word
    entries), most frequent first. Each is a whole term to put in place of `partial`."""
    p = (partial or "").strip().lower().replace('"', "")
    if not p:
        return []
    key = (id(rows), len(rows))
    hit = _ANN_VOCAB_MEMO.get(key)
    if hit is None or hit[0] is not rows:
        if len(_ANN_VOCAB_MEMO) > 16:
            _ANN_VOCAB_MEMO.clear()
        hit = _ANN_VOCAB_MEMO[key] = (rows, _ann_vocab(rows))
    a, b = [], []
    for txt, k, n in hit[1]:
        if k == p:
            continue
        bare = k.split(":", 1)[1] if ":" in k and k.split(":", 1)[0] in _ANN_FACETS else k
        if k.startswith(p) or bare.startswith(p):
            a.append((-n, len(k), txt))
        elif len(p) >= 2 and re.search(r"(^|[\s\-(/])" + re.escape(p), bare):
            b.append((-n, len(k), txt))
    return [t for _n, _l, t in sorted(a)[:limit]] + [t for _n, _l, t in sorted(b)[:max(0, limit - len(a))]]


def ann_filter_query(tokens=None, provenances=None, diseases=None, reviewed_only=False):
    """The old variant filters (_filter_variants arguments) as search terms selecting the same rows:
    consequence tokens (None / empty / all four = no term), provenance categories (None = all; empty set =
    src:none), reviewed-only (src:reviewed, ANDed), diseases (None = all; empty set = dis:none)."""
    parts = []
    toks = [t for t in ANN_CONS_TOKENS if t in (tokens or ())]
    if toks and len(toks) < len(ANN_CONS_TOKENS):
        parts.append("cons:" + "|".join(toks))
    if provenances is not None:
        ps = [_ANN_SRC_OF[c] for c in PROVENANCE_COLORS if c in provenances]
        if len(ps) < len(PROVENANCE_COLORS):
            parts.append("src:" + ("|".join(ps) if ps else "none"))
    if reviewed_only:
        parts.append("src:reviewed")
    if diseases is not None:
        parts.append("dis:" + ("|".join(_ann_quote(d) for d in sorted(diseases)) if diseases else "none"))
    return " ".join(parts)


def ann_query_facets(query):
    """Facet state of a query for the Filters menu: {cons, src, dis: set of values (None = no term),
    reviewed: bool}. Only single-facet terms count (cons:a|b ...); 'src:reviewed' alone = reviewed-only."""
    out = {"cons": None, "src": None, "dis": None, "reviewed": False}
    for raw, alts in _ann_lex(query):
        text, qat = alts[0]
        i = text.find(":")
        if not (i > 0 and (qat < 0 or i < qat) and text[:i].lower() in _ANN_FACETS):
            continue
        fac = _ANN_FACETS[text[:i].lower()]
        vals = [text[i + 1:].strip()] + [t.strip() for t, _q in alts[1:]]
        vals = [v for v in vals if v and v.lower() != "none"]
        if fac == "src" and [v.lower() for v in vals] == ["reviewed"]:
            out["reviewed"] = True
        elif fac == "cons":
            out["cons"] = {next((t for t in ANN_CONS_TOKENS if t.startswith(v.lower())), v.lower()) for v in vals}
        elif fac == "src":
            out["src"] = {_ANN_SRC_ALIAS.get(v.lower(), v.lower()) for v in vals}
        elif fac == "dis":
            out["dis"] = set(vals)
    return out


def ann_query_set(query, key, values):
    """`query` with the term for facet `key` (cons | src | dis | reviewed) replaced by `values` (None or an
    'everything' set removes it; reviewed: values = bool). Other terms are kept as typed."""
    keep = []
    for raw, alts in _ann_lex(query):
        text = alts[0][0] if alts else ""
        i = text.find(":")
        fac = _ANN_FACETS.get(text[:i].lower()) if i > 0 and (alts[0][1] < 0 or i < alts[0][1]) else None
        is_rev = fac == "src" and raw.lower().strip() == "src:reviewed"
        if (key == "reviewed" and is_rev) or (key != "reviewed" and fac == key and not is_rev):
            continue
        keep.append(raw)
    if key == "reviewed":
        if values:
            keep.append("src:reviewed")
    elif values is not None:
        vals = list(values)
        allv = ANN_CONS_TOKENS if key == "cons" else ANN_SRC_TOKENS if key == "src" else None
        if allv is not None:
            vals = [t for t in allv if t in vals]
            if len(vals) < len(allv):
                keep.append("%s:%s" % (key, "|".join(vals) if vals else "none"))
        else:
            keep.append("dis:" + ("|".join(_ann_quote(v) for v in vals) if vals else "none"))
    return " ".join(keep)


def ann_variant_groups(ann, variants, color="consequence", diseases=None):
    """{colour: [positions]} for a list of variants: one sphere per position, its most severe variant's
    colour under the profile (== _variant_groups on an already filtered list)."""
    col = _variant_colorer(ann, color, diseases)
    pos_col, pos_sev = {}, {}
    for v in variants:
        p = v["position"]
        sev = _SEVERITY.index(v["consequence"]) if v.get("consequence") in _SEVERITY else -1
        if p not in pos_col:
            pos_col[p], pos_sev[p] = col(v), sev
        elif sev >= 0 and (pos_sev[p] < 0 or sev < pos_sev[p]):
            pos_col[p], pos_sev[p] = col(v), sev
    groups = {}
    for p, c in pos_col.items():
        groups.setdefault(c, []).append(p)
    return groups


def _ann_public(r):
    return {k: r[k] for k in ("id", "cat", "pos", "end", "poslabel", "label", "color", "pts")}


def _ann_op_rows(a):
    uid = _arg_uid(a)
    ann = fetch_annotations(uid)
    extra = dict(a.get("extra") or {})
    if extra.get("kb") in (True, 1, "1", "true"):
        extra["kb"] = fetch_pdbekb(uid) or {}
    return ann, ann_rows(ann, a.get("cat") or "all", extra)


@_serve_op("ann_rows")
def _op_ann_rows(a):
    """Annotation-list rows of a category ('all' = every category): [{id, cat, pos, end, poslabel, label,
    color, pts}]. extra = {lig: [{resn, count, name}], prt: [{chain, resi, color, text}], kb: true|<pdbekb dict>}."""
    return [_ann_public(r) for r in _ann_op_rows(a)[1]]


@_serve_op("search")
def _op_search(a):
    """The shared annotation matcher: {uid, cat, query, extra?, suggest?} -> {ids, count, total, suggest}.
    suggest = completions for the query's last term (unless suggest=false)."""
    _ann, rows = _ann_op_rows(a)
    q = a.get("query") or ""
    ids = ann_search(rows, q)
    out = {"ids": ids, "count": len(ids), "total": len(rows)}
    if a.get("suggest", True) not in (False, 0, "0", "false"):
        out["suggest"] = ann_suggest(rows, ann_last_term(q)[1])
    return out


@_serve_op("filter_query")
def _op_filter_query(a):
    """Old variant filters -> the equivalent search terms (see ann_filter_query)."""
    def st(x):
        return None if x is None else set(x if isinstance(x, (list, tuple, set)) else str(x).split())
    return ann_filter_query(st(a.get("tokens")), st(a.get("provenances")), st(a.get("diseases")),
                            bool(a.get("reviewed")))
# ---- end annotation search (v3) ----


def _serve(stdin=None, stdout=None, workers=4):
    """JSON-lines request loop (see the section comment). Ends at stdin EOF or a 'shutdown' op."""
    import traceback
    from concurrent.futures import ThreadPoolExecutor
    if stdin is None and stdout is None:
        # Own private handles for the request pipes and point fds 0/1 at NUL / stderr. On Windows a synchronous
        # read pending on the stdin pipe blocks any other thread that touches that handle (e.g. numpy / platform
        # probing stdin while an op runs) - a deadlock; and C-level writes to fd 1 would corrupt the replies.
        try:
            sys.stdout.flush()
            fin, fout = os.dup(0), os.dup(1)
            nul = os.open(os.devnull, os.O_RDONLY)
            os.dup2(nul, 0)
            os.close(nul)
            os.dup2(2, 1)
            stdin = io.open(fin, "r", encoding="utf-8", errors="replace", newline=None)
            stdout = io.open(fout, "w", encoding="ascii", errors="replace", newline="\n")
            sys.stdin = open(os.devnull)
        except Exception:
            stdin = stdout = None
    stdin = stdin or sys.stdin
    out = stdout or sys.stdout
    sys.stdout = sys.stderr            # stray prints from library code must not corrupt the reply stream
    _numpy()                           # import numpy now (main thread), not inside the first op
    lock = threading.Lock()

    def reply(fmt, rid, ok, payload):
        memo = ok and _MEMO_OBJS.get(id(payload)) is payload
        if fmt == "tcl":
            body = _tcl_elem(payload)
            line = "id %s ok %d %s %s" % (_tcl_elem(rid), 1 if ok else 0, "result" if ok else "error", body)
        else:
            body = _encoded(payload, "json", lambda o: json.dumps(o, ensure_ascii=True, default=str)) if memo                 else json.dumps(payload, ensure_ascii=True, default=str)
            line = '{"id": %s, "ok": %s, "%s": %s}' % (json.dumps(rid, default=str), "true" if ok else "false",
                                                      "result" if ok else "error", body)
        with lock:
            out.write(line + "\n")
            out.flush()

    def run(req):
        rid, fmt = req.get("id"), req.get("fmt") or "json"
        fn = _SERVE_OPS.get(req.get("op"))
        try:
            if fn is None:
                raise ValueError("unknown op '%s'" % req.get("op"))
            res = fn(req.get("args") or {})
        except Exception as e:
            sys.stderr.write("[ufv serve] op %s failed:\n%s" % (req.get("op"), traceback.format_exc()))
            reply(fmt, rid, False, "%s: %s" % (type(e).__name__, e))
            return
        try:
            reply(fmt, rid, True, res)
        except Exception as e:
            reply(fmt, rid, False, "reply encoding failed: %s" % e)

    ex = ThreadPoolExecutor(max_workers=workers)
    try:
        for line in stdin:
            line = line.strip()
            if not line:
                continue
            try:
                req = json.loads(line)
                if not isinstance(req, dict):
                    raise ValueError("request must be a JSON object")
            except Exception as e:
                reply("json", None, False, "bad request: %s" % e)
                continue
            if req.get("op") == "shutdown":
                reply(req.get("fmt") or "json", req.get("id"), True, "bye")
                break
            ex.submit(run, req)
    finally:
        ex.shutdown(wait=True)      # finish requests already accepted (stdin EOF = client done sending)
    out.flush()


if __name__ == "__main__":
    _args = sys.argv[1:]
    if len(_args) >= 2 and _args[0] == "--emit-tcl":
        _emit_tcl_annotations(_args[1])
    elif len(_args) >= 3 and _args[0] == "--emit-sifts":
        _emit_tcl_sifts(_args[1], _args[2])
    elif len(_args) >= 2 and _args[0] == "--download":
        _download_model(_args[1])
    elif len(_args) >= 2 and _args[0] == "--structures":
        _emit_tcl_structures(_args[1])
    elif len(_args) >= 3 and _args[0] == "--get":
        _download_url(_args[1], _args[2] if len(_args) > 2 else "pdb")
    elif len(_args) >= 3 and _args[0] == "--residue":
        _emit_tcl_residue(_args[1], _args[2])
    elif len(_args) >= 3 and _args[0] == "--compute-tiers":
        _emit_tcl_tiers(_args[1], _args[2])
    elif len(_args) >= 3 and _args[0] == "--csv":
        _emit_csv_file(_args[1], _args[2])
    elif len(_args) >= 1 and _args[0] == "--serve":
        _serve()
        sys.stdout.flush()
        os._exit(0)                 # do not wait for worker threads still busy with abandoned requests
    else:
        sys.stderr.write(
            "3D Feature Viewer — backend for the VMD plugin\n"
            "usage:\n"
            "  python ufv_core.py --emit-tcl       <uniprot_id>\n"
            "  python ufv_core.py --emit-sifts     <pdb_id> <uniprot_id>\n"
            "  python ufv_core.py --structures     <uniprot_id>\n"
            "  python ufv_core.py --download        <uniprot_id>\n"
            "  python ufv_core.py --get            <url> <fmt>\n"
            "  python ufv_core.py --residue        <uniprot_id> <pos>\n"
            "  python ufv_core.py --compute-tiers  <uniprot_id> <coords.json>\n"
            "  python ufv_core.py --csv            <uniprot_id> <outfile.csv>\n")
        sys.exit(2)
