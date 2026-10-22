# =============================================================================================
# ufv_commands.tcl - the user-level ufv_* commands on the backend / palette / layer foundation.
# OWNER: features agent (was foundation). Trajectory / GUI code only calls these.
# Network work is implemented once as an async operation (::ufv::*_async ... cb) and the console
# commands wait for it with ::ufv::await, so the GUI can call the async form and never block.
# New feature commands (diseases, partners, PDBe-KB, ProtVar pockets, tractability, AlphaFill, motif search,
# exports) are in ufv_features.tcl. The GUI-facing API of both files is documented in API.md ("Features API").
# =============================================================================================

namespace eval ::ufv {
    if {![info exists pocket_q]} { variable pocket_q 0.10 }             ;# constraint-pocket sensitivity (FDR q)
    if {![info exists prox]}     { variable prox [dict create ptm 8.0 var 12.0] }  ;# PTM / variant radii (A)
    if {![info exists aff_min]}  { variable aff_min 0.25 }              ;# AlphaFill donor-identity cutoff
}
namespace eval ::ufv::feat {
    variable seq 0
    if {![array exists acmap]} { array set acmap {} }
}

# ION_CODES of the extension (viewer-molstar.js): drawn as small spheres in the ligand layer
set ::ufv::ion_codes {NA K LI RB CS MG CA SR BA ZN FE FE2 MN MN3 CU CU1 CO NI CD HG CL BR IOD FLO F AL PT AU AG PB SO4 PO4}
# Polymer as PyMOL's `polymer`: VMD's protein flag needs a complete, bonded backbone (incomplete residues in
# crystal structures, or residues next to overlapping AlphaFill transplants, drop out), so standard residue names
# count too. Ligands = everything else except water.
set ::ufv::poly_base "(protein or nucleic or resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE SEC PYL HID HIE HIP HSD HSE HSP CYX DA DC DG DT DU)"
set ::ufv::lig_base "not $::ufv::poly_base and not water"

# ---------------------------------------------------------------------------------------------
# Fetch (annotations + structure list in ONE backend request)
# ---------------------------------------------------------------------------------------------
proc ::ufv::fetch_async {uid force cb} {
    set uid [string toupper [string trim $uid]]
    if {$uid eq ""} { error "no accession given" }
    ::ufv::status "fetching $uid ..."
    ::ufv::backend::async fetch [dict create uid $uid force [expr {$force ? 1 : 0}]] \
        [list ::ufv::_fetched $uid $cb]
}

proc ::ufv::_fetched {uid cb ok res} {
    if {!$ok} {
        ::ufv::status "fetch $uid failed: $res"
        {*}$cb 0 $res
        return
    }
    set a [dict get $res annotations]
    set ::ufv::ann($uid) $a
    set ::ufv::structs($uid) [dict get $res structures]
    unset -nocomplain ::ufv::feat::acmap($uid)
    array unset ::ufv::annot::rowcache $uid,*
    set ::ufv::uid $uid
    ::ufv::status [format "%s: %d PTMs, %d variants, %d sites, %d domains, %d topology, %d mutagenesis, %d structures" \
        $uid [llength [dict get $a ptms]] [llength [dict get $a variants]] [llength [dict get $a sites]] \
        [llength [dict get $a domains]] [llength [dict get $a topology]] [llength [dict get $a mutagenesis]] \
        [llength [dict get $res structures]]]
    # a structure source that could not be reached (3D-Beacons down): one visible line, not a silently short list
    if {[dict exists $res structures_notice] && [dict get $res structures_notice] ne ""} {
        ::ufv::status [dict get $res structures_notice]
    }
    ::ufv::emit fetched $uid
    {*}$cb 1 $uid
}

# ufv_fetch <uid> ?refresh?
proc ufv_fetch {uid {refresh ""}} {
    set force [expr {$refresh in {refresh -refresh force 1}}]
    return [::ufv::await [list ::ufv::fetch_async $uid $force]]
}

# Make sure uid's annotations are present (fetches when missing). "" = current accession.
proc ::ufv::need {{uid ""}} {
    if {$uid eq ""} { set uid $::ufv::uid }
    if {$uid eq ""} { error "no accession fetched yet - run: ufv_fetch <uniprot_id>" }
    set uid [string toupper $uid]
    if {![info exists ::ufv::ann($uid)]} { ufv_fetch $uid }
    return $uid
}

# Accession of a molecule (its own, else the current one).
proc ::ufv::mol_uid {mol} {
    set u [::ufv::st::get $mol uid]
    if {$u eq ""} { set u [::ufv::need] ; ::ufv::st::put $mol uid $u }
    return $u
}

# ufv_structures ?uid?  - list the structures of an accession (AlphaFold, isoforms, PDB entries with their chains
# merged, 3D-Beacons computed models incl. SWISS-MODEL / ModelArchive / AlphaFill). Returns the list of dicts.
proc ufv_structures {{uid ""}} {
    set uid [::ufv::need $uid]
    if {![info exists ::ufv::structs($uid)]} {
        set ::ufv::structs($uid) [::ufv::backend::call structures [dict create uid $uid]]
        set note [::ufv::backend::call structure_notice [dict create uid $uid]]
        if {$note ne ""} { ::ufv::status $note }
    }
    puts "\[UFV\] [llength $::ufv::structs($uid)] structures for $uid:"
    foreach s $::ufv::structs($uid) {
        puts [format "    %-26s %s" [dict get $s key] [dict get $s label]]
    }
    return $::ufv::structs($uid)
}

# Structure-list entry by key (case-insensitive, like the PyMOL plugin), "" if not listed.
proc ::ufv::struct_entry {uid key} {
    if {![info exists ::ufv::structs($uid)]} { return "" }
    foreach s $::ufv::structs($uid) {
        if {[string equal -nocase [dict get $s key] $key]} { return $s }
    }
    return ""
}

# ---------------------------------------------------------------------------------------------
# Loading: AlphaFold (identity), isoform models (ranges), PDB entries (chains merged, SIFTS verified against
# the model, chimeric / partner metadata), 3D-Beacons computed models (identity; AlphaFill transplant chains
# excluded from the numbering).
# ---------------------------------------------------------------------------------------------
proc ::ufv::_newmol {path name} {
    # VMD re-fits every molecule's view when a further molecule is loaded: keep the view (first load: VMD's fit)
    set m [::ufv::view::keep { mol new $path waitfor all }]
    catch { mol rename $m $name }
    ::ufv::st::init $m
    set ::ufv::molid $m
    ::ufv::cartoon::ensure $m 1
    ::ufv::ss::start $m
    return $m
}

# ---------------------------------------------------------------------------------------------
# Secondary structure without freezing the UI. VMD runs STRIDE synchronously (an external process + file I/O,
# ~0.35 s for 2000 residues) the first time a cartoon of a new molecule is drawn. For the molecules the plugin
# loads, the same STRIDE binary ($env(STRIDE_BIN), the one VMD's ssrecalc uses) runs on the protein atoms as a
# pipe read from the event loop; its assignment is set per residue in order (identical to VMD's own run: checked
# atom by atom). Until it arrives the protein carries coil (that stops VMD's own run) and the base cartoon is
# hidden. Any mismatch / failure / no binary -> VMD's own `mol ssrecalc` (the old behaviour).
#   ::ufv::ss::when mol script   runs script once the assignment is in (now if nothing is pending).
# ---------------------------------------------------------------------------------------------
namespace eval ::ufv::ss {
    variable pend
    variable seq 0
    array set pend {}
}
proc ::ufv::ss::bin {} {
    if {![info exists ::env(STRIDE_BIN)]} { return "" }
    set b [string trim $::env(STRIDE_BIN) "\" "]
    if {$b eq "" || ![file isfile $b]} { return "" }
    return $b
}
proc ::ufv::ss::pending {mol} { variable pend; return [info exists pend($mol)] }
# Synchronous part (cheap): coil placeholder + hidden base cartoon. The protein is written and STRIDE started in
# the next event-loop slice (keeps `mol new` and the PDB write in separate slices).
proc ::ufv::ss::start {mol} {
    variable pend
    if {[info exists pend($mol)]} return
    if {[bin] eq ""} return
    set sel [atomselect $mol protein]
    set n [$sel num]
    if {$n} { $sel set structure C }
    $sel delete
    if {!$n} return
    set base [lindex [::ufv::layers::reps $mol cartoon] 0]
    if {$base ne ""} { mol showrep $mol $base 0 }
    set pend($mol) [dict create fh "" tmp "" buf "" waiters {} \
        timer [after 30000 [list ::ufv::ss::_finish $mol 0]]]
    after 0 [list ::ufv::ss::_run $mol]
}
proc ::ufv::ss::_run {mol} {
    variable pend; variable seq
    if {![info exists pend($mol)]} return
    if {$mol ni [molinfo list]} { _finish $mol 0; return }
    set tmp [file join [::ufv::backend::tmpdir] "ufv_ss_[pid]_${mol}_[incr seq].pdb"]
    dict set pend($mol) tmp $tmp
    if {[catch {
        ::ufv::with_sel sel $mol protein { $sel writepdb $tmp }
        set fh [open |[list [bin] $tmp 2>@1] r]
    } err]} { _finish $mol 0; return }
    dict set pend($mol) fh $fh
    fconfigure $fh -blocking 0
    fileevent $fh readable [list ::ufv::ss::_read $mol]
}
proc ::ufv::ss::_read {mol} {
    variable pend
    if {![info exists pend($mol)]} return
    set fh [dict get $pend($mol) fh]
    if {[catch { read $fh } chunk]} { _finish $mol 0; return }
    dict append pend($mol) buf $chunk
    if {[eof $fh]} { _finish $mol 1 }
}
proc ::ufv::ss::_finish {mol ok} {
    variable pend
    if {![info exists pend($mol)]} return
    set p $pend($mol)
    unset pend($mol)
    catch { after cancel [dict get $p timer] }
    set fh [dict get $p fh]
    if {$fh ne ""} {
        catch { fileevent $fh readable {} }
        catch { close $fh }
    }
    if {[dict get $p tmp] ne ""} { catch { file delete [dict get $p tmp] } }
    if {$mol in [molinfo list]} {
        set codes {}
        if {$ok} {
            foreach l [split [dict get $p buf] \n] { if {[string range $l 0 2] eq "ASG"} { lappend codes [string index $l 24] } }
        }
        if {![apply_codes $mol $codes]} { catch { mol ssrecalc $mol } }
        if {[::ufv::st::get $mol {layers cartoon visible} 1]} {
            set base [lindex [::ufv::layers::reps $mol cartoon] 0]
            if {$base ne ""} { mol showrep $mol $base 1 }
        }
    }
    foreach w [dict get $p waiters] {
        if {[catch { uplevel #0 $w } err]} { puts "\[UFV\] $err" }
    }
}
# STRIDE's ASG records come one per residue in file order (= atom order of the protein selection).
proc ::ufv::ss::apply_codes {mol codes} {
    if {![llength $codes]} { return 0 }
    set sel [atomselect $mol protein]
    set res [$sel get residue]
    set order [dict create]
    foreach r $res { if {![dict exists $order $r]} { dict set order $r [dict size $order] } }
    set ok [expr {[dict size $order] == [llength $codes]}]
    if {$ok} {
        foreach c $codes { if {$c ni {H G I E B T C}} { set ok 0; break } }
    }
    if {$ok} { $sel set structure [lmap r $res { lindex $codes [dict get $order $r] }] }
    $sel delete
    return $ok
}
proc ::ufv::ss::when {mol script} {
    variable pend
    if {[info exists pend($mol)]} { dict lappend pend($mol) waiters $script; return }
    uplevel #0 $script
}

# Trailing callback + leading options (-name N): async procs keep the callback LAST (for ::ufv::await).
proc ::ufv::_split_cb {argl {defaults {-name ""}}} {
    set cb [lindex $argl end]
    set o [dict merge $defaults [lrange $argl 0 end-1]]
    return [list $o $cb]
}

# AlphaFold model of uid: download + load + identity numbering (+ annotations fetched if missing).
#   ::ufv::load_af_async uid ?-name N? cb
proc ::ufv::load_af_async {uid args} {
    lassign [::ufv::_split_cb $args] o cb
    set uid [string toupper [string trim $uid]]
    ::ufv::status "downloading AlphaFold model of $uid ..."
    ::ufv::backend::async download [dict create uid $uid] [list ::ufv::_af_downloaded $uid [dict get $o -name] $cb]
}
proc ::ufv::_af_downloaded {uid name cb ok res} {
    if {!$ok} { ::ufv::status "download failed for $uid: $res"; {*}$cb 0 $res; return }
    if {[catch {
        set m [::ufv::_newmol [dict get $res path] [expr {$name ne "" ? $name : "AF-$uid"}]]
        catch { file delete [dict get $res path] }
        ::ufv::st::put $m uid $uid
        ::ufv::st::put $m source AlphaFold
        ::ufv::st::put $m struct [dict create key AF-$uid label "AlphaFold $uid" source AlphaFold numbering identity \
            url [dict get $res url] pdbId "" isoform "" provider "" fmt pdb]
        ::ufv::map::install $m [::ufv::feat::identity_map $m]
    } err]} { {*}$cb 0 $err; return }
    set next [list ::ufv::_af_loaded $m $uid $cb]
    if {[info exists ::ufv::ann($uid)]} { {*}$next 1 $uid } else { ::ufv::fetch_async $uid 0 $next }
}
proc ::ufv::_af_loaded {m uid cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    set ::ufv::uid $uid
    ::ufv::status "loaded AlphaFold $uid as mol $m - add layers: ufv_ptms / ufv_variants / ufv_sites / ufv_domains ..."
    ::ufv::emit loaded $m AlphaFold
    ::ufv::ss::when $m [list {*}$cb 1 $m]
}

# Looks like a structure-list key (AF-<iso>, BCN-<id>, <pdb>_<chains>) rather than an accession?
proc ::ufv::_is_key {s} {
    return [expr {[regexp {^(AF-|BCN-)} $s] || [regexp {_|,} $s]}]
}

# ufv_load <uid|key> ?name?   AlphaFold model of an accession, or any structure listed by ufv_structures.
proc ufv_load {uid {name ""}} {
    if {[::ufv::_is_key $uid]} { return [ufv_load_structure $uid "" $name] }
    return [::ufv::await [list ::ufv::load_af_async $uid -name $name]]
}

# A listed structure (key from ufv_structures): download (+ SIFTS in parallel), load, map.
#   ::ufv::load_structure_async key uid ?-name N? cb
proc ::ufv::load_structure_async {key uid args} {
    lassign [::ufv::_split_cb $args] o cb
    if {$uid eq ""} { set uid $::ufv::uid }
    if {$uid eq ""} { {*}$cb 0 "fetch an accession first"; return }
    set uid [string toupper $uid]
    set s [::ufv::struct_entry $uid $key]
    if {$s ne ""} { set key [dict get $s key] }
    ::ufv::status "downloading $key ..."
    ::ufv::backend::async prepare [dict create uid $uid key $key] [list ::ufv::_prepared $uid $key [dict get $o -name] $cb]
}
proc ::ufv::_prepared {uid key name cb ok res} {
    if {!$ok} { ::ufv::status "loading $key failed: $res"; {*}$cb 0 $res; return }
    # mmCIF: VMD's reader takes label numbering; fetch the author chain / residue ids first (UniProt numbering)
    if {[string tolower [file extension [dict get $res path]]] eq ".cif"} {
        ::ufv::backend::async cif_auth [dict create path [dict get $res path]] [list ::ufv::_prepared2 $uid $key $name $cb $res]
        return
    }
    ::ufv::_prepared2 $uid $key $name $cb $res 0 ""
}
# Author chain / resid / insertion from the backend's mmCIF parse (checked: same atom count and names).
proc ::ufv::feat::apply_auth {mol auth} {
    if {![dict exists $auth n] || [dict get $auth n] != [molinfo $mol get numatoms]} { return 0 }
    set all [atomselect $mol all]
    set ok [expr {[$all get name] eq [dict get $auth name]}]
    if {$ok} {
        $all set chain [dict get $auth chain]
        $all set resid [dict get $auth resid]
        catch { $all set insertion [lmap i [dict get $auth ins] { expr {$i eq "" ? " " : $i} }] }
        ::ufv::st::put $mol chains ""
    }
    $all delete
    return $ok
}
proc ::ufv::_prepared2 {uid key name cb res aok auth} {
    set st [dict get $res struct]
    if {[catch {
        set m [::ufv::_newmol [dict get $res path] [expr {$name ne "" ? $name : [string map {, _} $key]}]]
        if {$aok && ![::ufv::feat::apply_auth $m $auth]} {
            ::ufv::status "mol $m: mmCIF author numbering not applied (atom list mismatch) - label numbering kept."
        }
        catch { file delete [dict get $res path] }
        ::ufv::st::put $m uid $uid
        ::ufv::st::put $m source [dict get $st source]
        ::ufv::st::put $m struct $st
    } err]} { {*}$cb 0 $err; return }
    set ::ufv::uid $uid
    switch -- [dict get $st numbering] {
        sifts {
            ::ufv::map_sifts_async $m [dict get $st pdbId] $uid [list ::ufv::_struct_mapped $m $key $cb]
            return
        }
        ranges {
            set segs [dict create]
            foreach ch [::ufv::feat::our_chains $m] { dict set segs $ch [dict get $st ranges] }
            ::ufv::map::install $m [dict create mode ranges segments $segs status "isoform ranges"]
        }
        default { ::ufv::map::install $m [::ufv::feat::identity_map $m] }
    }
    ::ufv::_struct_mapped $m $key $cb 1 ""
}
proc ::ufv::_struct_mapped {m key cb ok res} {
    if {!$ok} { ::ufv::status "loaded $key as mol $m, but numbering failed: $res" }
    set mp [::ufv::map::get $m]
    if {[dict exists $mp meta chimeric] && [dict get $mp meta chimeric] eq "1"} {
        ::ufv::status "\[chimeric\] mol $m: chain(s) [join [dict get $mp meta chimericChains] ,] span\
            [join [dict get $mp meta organisms] { / }] - residues outside the mapping stay unannotated (ufv_chimera highlights them)."
    }
    set tx [::ufv::feat::transplant_chains $m]
    if {[llength $tx]} { ::ufv::status "mol $m: AlphaFill transplant chain(s) [join $tx {, }] are not numbered as this protein." }
    ::ufv::status "loaded $key as mol $m ([::ufv::map::describe $m])."
    ::ufv::emit loaded $m $key
    ::ufv::ss::when $m [list {*}$cb 1 $m]
}
proc ufv_load_structure {key {uid ""} {name ""}} {
    return [::ufv::await [list ::ufv::load_structure_async $key $uid -name $name]]
}

# ufv_use ?molid? ?uid? ?key?   adopt an already loaded molecule (your own PDB, a trajectory) as the current one;
#                                 with a ufv_structures key (e.g. AF-P35498) the molecule is that structure (a local
#                                 copy / MD topology of it): its entry (source, AlphaFold PAE gate, motif search) applies.
# ufv_use <key> ?uid?            (PyMOL form) download + load a structure listed by ufv_structures.
proc ufv_use {{molid ""} {uid ""} {key ""}} {
    if {$molid ne "" && $molid ne "top" && ![string is integer -strict $molid]} {
        return [ufv_load_structure $molid $uid]
    }
    set m [::ufv::resolve $molid]
    ::ufv::st::init $m
    set ::ufv::molid $m
    if {$uid ne ""} { set uid [::ufv::need $uid]; ::ufv::st::put $m uid $uid }
    if {$key ne ""} {
        set u [::ufv::mol_uid $m]
        if {$u eq ""} { error "ufv_use: give the accession with the structure key" }
        if {![info exists ::ufv::structs($u)]} { ufv_structures $u }
        set s [::ufv::struct_entry $u $key]
        if {$s eq ""} { error "ufv_use: '$key' is not listed by ufv_structures $u" }
        ::ufv::st::put $m struct $s
        ::ufv::st::put $m source [dict get $s source]
    }
    ::ufv::cartoon::ensure $m
    ::ufv::status "current molecule: $m ([::ufv::map::describe $m])"
    return $m
}

# ---------------------------------------------------------------------------------------------
# Numbering
# ---------------------------------------------------------------------------------------------
proc ::ufv::map_sifts_async {mol pdb uid cb} {
    if {$uid eq ""} { set uid [::ufv::mol_uid $mol] }
    set args [dict create uid $uid pdb $pdb residues [::ufv::json::raw [::ufv::json::residues [::ufv::map::residues $mol]]]]
    ::ufv::backend::async map_sifts $args [list ::ufv::_sifts_mapped $mol $cb]
}
proc ::ufv::_sifts_mapped {mol cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; return }
    foreach n [dict get $res notes] { puts "\[UFV\] $n" }
    ::ufv::map::install $mol $res
    {*}$cb 1 [dict get $res status]
}

# ufv_map identity | sifts <pdb> | ranges {{us ue rs re} ...} | auto     (current molecule; -mol N to choose)
proc ufv_map {args} {
    set mol ""
    if {[set i [lsearch -exact $args -mol]] >= 0} {
        set mol [lindex $args $i+1]; set args [lreplace $args $i $i+1]
    }
    set mol [::ufv::resolve $mol]
    set mode [lindex $args 0]
    switch -- $mode {
        sifts {
            set pdb [lindex $args 1]
            if {$pdb eq ""} { error "usage: ufv_map sifts <pdb_id>" }
            ::ufv::await [list ::ufv::map_sifts_async $mol $pdb ""]
        }
        ranges {
            set segs [dict create]
            foreach ch [::ufv::map::protein_chains $mol] { dict set segs $ch [lindex $args 1] }
            ::ufv::map::install $mol [dict create mode ranges segments $segs]
        }
        auto - custom { return [ufv_custom $mol] }
        identity - "" { ::ufv::map::install $mol [::ufv::feat::identity_map $mol] }
        default { error "usage: ufv_map identity | sifts <pdb> | ranges <segments> | auto ?-mol N?" }
    }
    ::ufv::status "mol $mol numbering: [::ufv::map::describe $mol]"
    return [::ufv::map::get $mol]
}

# ufv_chain <chain> <uniprot_start> ?resid_start? ?resid_end?  (manual anchor; other manual chains kept)
proc ufv_chain {chain ustart {rstart ""} {rend ""} {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {$rstart eq ""} { set rstart $ustart }
    set m [::ufv::map::get $mol]
    set segs [expr {[dict get $m mode] eq "manual" ? [dict get $m segments] : [dict create]}]
    dict set segs $chain [::ufv::map::manual_segments $mol $chain $ustart $rstart $rend]
    ::ufv::map::install $mol [dict create mode manual segments $segs]
    ::ufv::status "chain $chain: resid $rstart == UniProt $ustart[expr {$rend ne {} ? " (..$rend)" : {}}]"
}

# Custom per-chain mapping (the extension's Custom PDB): every chain auto | offset | adopt | ignore.
proc ::ufv::map_custom_async {mol configs cb} {
    set uid [::ufv::mol_uid $mol]
    set m [::ufv::map::get $mol]
    # 'auto' keeps a known (non-custom) map: the base. The entry metadata (pdbId, partners, chimeric) is kept.
    if {[dict get $m mode] eq "custom"} {
        set base [dict get $m base]; set basemode [dict get $m basemode]
        set old [dict get $m config]
        set keep [expr {[dict exists $m basemeta] ? [dict get $m basemeta] : {}}]
    } else {
        set base [dict get $m segments]; set basemode [dict get $m mode]; set old [dict create]
        set keep [dict create]
        foreach k {meta pdbId} { if {[dict exists $m $k]} { dict set keep $k [dict get $m $k] } }
    }
    set cfg $old
    dict for {ch c} $configs {
        dict set cfg $ch $c
        # adopt: resolve the source molecule's segments here ({mol N chain X})
        if {[dict exists $c mode] && [dict get $c mode] eq "adopt"} {
            set src [dict get $c src]; set sch [expr {[dict exists $c src_chain] ? [dict get $c src_chain] : $ch}]
            set sm [::ufv::map::get $src]
            if {[dict get $sm mode] eq "identity"} {
                set sl [list [list 1 100000 1 100000]]
            } else {
                set sl [expr {[dict exists $sm segments $sch] ? [dict get $sm segments $sch] : {}}]
            }
            dict set cfg $ch segments $sl
            dict set cfg $ch src "mol $src/$sch"
        }
    }
    set args [dict create uid $uid basemode $basemode \
        residues [::ufv::json::raw [::ufv::json::residues [::ufv::map::residues $mol]]] \
        base [::ufv::json::raw [::ufv::json::segments $base]] configs [::ufv::json::raw [::ufv::json::configs $cfg]]]
    ::ufv::backend::async map_custom [dict replace $args base_mode $basemode] \
        [list ::ufv::_custom_mapped $mol $base $basemode $keep $cb]
}
proc ::ufv::_custom_mapped {mol base basemode keep cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    dict set res base $base
    dict set res basemode $basemode
    dict set res basemeta $keep
    dict for {k v} $keep { dict set res $k $v }
    ::ufv::map::install $mol $res
    puts "\[UFV\] mol $mol - per-chain mapping for [dict get $res uid]:"
    dict for {ch note} [dict get $res notes] {
        puts [format "    %-5s %-11s %-7s %s" $ch [dict get $res accessions $ch] \
            [expr {[dict exists $res config $ch mode] ? [dict get $res config $ch mode] : "auto"}] $note]
    }
    {*}$cb 1 [dict get $res status]
}

proc ufv_custom {{molid ""} {uid ""}} {
    set mol [::ufv::resolve $molid]
    if {$uid ne ""} { ::ufv::st::put $mol uid [::ufv::need $uid] }
    set cfg [dict create]
    foreach ch [::ufv::map::protein_chains $mol] { dict set cfg $ch [dict create mode auto] }
    return [::ufv::await [list ::ufv::map_custom_async $mol $cfg]]
}

# ufv_chainmap <chain|*> auto|offset|adopt|ignore ?value? ?accession? ?molid?
#   offset: value = UniProt position of the chain's first residue; adopt: value = molid[/chain]
#   accession: assign the chain to another UniProt entry (a partner chain, not annotated with this protein)
proc ufv_chainmap {chain mode {value ""} {accession ""} {molid ""}} {
    set mol [::ufv::resolve $molid]
    set chs [expr {$chain eq "*" ? [::ufv::map::protein_chains $mol] : [list $chain]}]
    set cfg [dict create]
    foreach ch $chs {
        set c [dict create mode $mode]
        switch -- $mode {
            offset { dict set c start $value }
            adopt  { lassign [split $value /] src sch; dict set c src $src; if {$sch ne ""} { dict set c src_chain $sch } }
            auto - ignore {}
            default { error "mode must be auto | offset | adopt | ignore" }
        }
        if {$accession ne ""} { dict set c accession [string toupper $accession] }
        dict set cfg $ch $c
    }
    return [::ufv::await [list ::ufv::map_custom_async $mol $cfg]]
}

# ---------------------------------------------------------------------------------------------
# Molecule context for the backend's feature ops (virtual object: numbering + structure + Ca geometry)
# ---------------------------------------------------------------------------------------------
proc ::ufv::feat::jstrs {l} {
    set p {}
    foreach x $l { lappend p [::ufv::json::str $x] }
    return "\[[join $p ,]\]"
}
proc ::ufv::feat::jints {l} {
    set p {}
    foreach x $l { if {[string is integer -strict $x]} { lappend p $x } }
    return "\[[join $p ,]\]"
}
# Every value as a JSON string (structure fields such as a pdbId "1E12" must not become numbers).
proc ::ufv::feat::jstrobj {d} {
    set p {}
    dict for {k v} $d { lappend p "[::ufv::json::str $k]:[::ufv::json::str $v]" }
    return "{[join $p ,]}"
}

proc ::ufv::feat::struct {mol} {
    set st [::ufv::st::get $mol struct]
    set out [dict create]
    foreach k {key label source url fmt numbering pdbId isoform provider} {
        if {[dict exists $st $k]} { dict set out $k [dict get $st $k] }
    }
    return $out
}

# Is the molecule an AlphaFold-based model (AlphaFold or AlphaFill, no PDB id)? (modal.js isAlphaFoldBased)
proc ::ufv::feat::is_alphafill {mol} {
    set st [::ufv::st::get $mol struct]
    set txt ""
    foreach k {provider source label} { if {[dict exists $st $k]} { append txt " " [dict get $st $k] } }
    return [regexp -nocase alphafill $txt]
}
proc ::ufv::feat::af_based {mol} {
    set st [::ufv::st::get $mol struct]
    if {[dict exists $st pdbId] && [dict get $st pdbId] ni {"" -}} { return 0 }
    set src [expr {[dict exists $st source] ? [dict get $st source] : [::ufv::st::get $mol source]}]
    return [expr {$src eq "AlphaFold" || [is_alphafill $mol]}]
}

# AlphaFill transplant chains (viewer-molstar.js _focusNeighbourhood): chains carrying a non-water HETATM that
# are not the largest chain - AlphaFold-based models ONLY (in an experimental complex every glycosylated subunit
# carries HETATMs). Cached per molecule.
proc ::ufv::feat::transplant_chains {mol} {
    set c [::ufv::st::get $mol transplant "-"]
    if {$c ne "-"} { return $c }
    set out {}
    if {[af_based $mol]} {
        set het [::ufv::with_sel s $mol "hetero and not water" { lsort -unique [$s get chain] }]
        if {[llength $het]} {
            set n [dict create]
            ::ufv::with_sel s $mol "not water" {
                set seen [dict create]
                foreach {ch ri} [join [$s get {chain residue}]] {
                    if {![dict exists $seen $ri]} { dict set seen $ri 1; dict incr n $ch }
                }
            }
            set main ""; set best -1
            dict for {ch k} $n { if {$k > $best} { set best $k; set main $ch } }
            foreach ch $het { if {$ch ne $main} { lappend out $ch } }
        }
    }
    ::ufv::st::put $mol transplant $out
    return $out
}

# Chains that are this protein: protein chains minus AlphaFill transplant chains.
proc ::ufv::feat::our_chains {mol} {
    set tx [transplant_chains $mol]
    set out {}
    foreach ch [::ufv::map::protein_chains $mol] { if {$ch ni $tx} { lappend out $ch } }
    if {![llength $out]} { set out [::ufv::map::protein_chains $mol] }
    return $out
}

proc ::ufv::feat::identity_map {mol} {
    return [dict create mode identity chains [our_chains $mol]]
}

# " and not chain X Y" for the transplant chains ("" when none)
proc ::ufv::feat::notx {mol} {
    set tx [transplant_chains $mol]
    if {![llength $tx]} { return "" }
    set t {}
    foreach ch $tx { lappend t [::ufv::map::chain_token $ch] }
    return " and not chain [join $t { }]"
}

proc ::ufv::feat::mapjson {mol} {
    set m [::ufv::map::get $mol]
    set mode [dict get $m mode]
    set p [list "\"mode\":[::ufv::json::str $mode]"]
    if {$mode eq "identity"} {
        set ex {}
        foreach ch [::ufv::map::protein_chains $mol] { if {$ch ni [dict get $m chains]} { lappend ex $ch } }
        lappend p "\"exclude\":[jstrs $ex]"
    } else {
        lappend p "\"segments\":[::ufv::json::segments [dict get $m segments]]"
    }
    set pdb ""
    if {[dict exists $m pdbId]} { set pdb [dict get $m pdbId] }
    if {$pdb ne ""} { lappend p "\"pdbId\":[::ufv::json::str $pdb]" }
    return "{[join $p ,]}"
}

# Ca geometry of the CURRENT frame as the backend's JSON ([{uniPos chain resi x y z b}], model order, first
# (chain, resid) wins, insertion-code residues unmapped - as the PyMOL plugin's _ca_geometry). Written once per
# numbering version + frame (+ frame count) and reused; returns the file path.
proc ::ufv::feat::coords {mol} {
    set key [list [::ufv::st::get $mol mapver 0] [molinfo $mol get frame] [molinfo $mol get numframes]]
    set c [::ufv::st::get $mol {cache coords} ""]
    if {$c ne "" && [lindex $c 0] eq $key && [file exists [lindex $c 1]]} { return [lindex $c 1] }
    variable seq
    set path [file join [::ufv::backend::tmpdir] "ufv_ca_[pid]_${mol}_[incr seq].json"]
    ::ufv::export_ca $mol $path
    set old [::ufv::st::get $mol {cache coordfiles} {}]
    lappend old $path
    # keep the last 3 files (a request may still be reading the previous one)
    while {[llength $old] > 3} { catch { file delete [lindex $old 0] }; set old [lrange $old 1 end] }
    ::ufv::st::put $mol {cache coordfiles} $old
    ::ufv::st::put $mol {cache coords} [list $key $path]
    return $path
}

# Request arguments describing a molecule for the feature ops: uid, map, struct (+ coords).
proc ::ufv::feat::ctx {mol {with_coords 1}} {
    set a [dict create uid [::ufv::mol_uid $mol] map [::ufv::json::raw [mapjson $mol]] \
               struct [::ufv::json::raw [jstrobj [struct $mol]]]]
    if {$with_coords} { dict set a coords [coords $mol] }
    dict set a ptm_radius [dict get $::ufv::prox ptm]
    dict set a var_radius [dict get $::ufv::prox var]
    dict set a q $::ufv::pocket_q
    return $a
}

# uniprot pos -> annotation colour (variant > PTM > site; most severe variant wins): _annotation_color_map.
proc ::ufv::feat::anno_colors {uid} {
    variable acmap
    if {[info exists acmap($uid)]} { return $acmap($uid) }
    set a [::ufv::annotations $uid]
    set m [dict create]
    foreach s [dict get $a sites] { dict set m [dict get $s pos] [dict get $s color] }
    foreach p [dict get $a ptms] { dict set m [dict get $p pos] [dict get $p color] }
    set rank [dict create pathogenic 4 deleterious 3 uncertain 2 benign 1]
    set best [dict create]
    foreach v [dict get $a variants] {
        set p [dict get $v pos]
        set r [expr {[dict exists $rank [dict get $v tok]] ? [dict get $rank [dict get $v tok]] : 0}]
        if {![dict exists $best $p] || [dict get $best $p] < $r} {
            dict set best $p $r
            dict set m $p [dict get $v color]
        }
    }
    set acmap($uid) $m
    return $m
}

# ---------------------------------------------------------------------------------------------
# Marker layers (CA spheres): groups computed by the backend exactly like the PyMOL commands
# ---------------------------------------------------------------------------------------------
proc ::ufv::groups {uid layer {opts {}}} {
    return [::ufv::backend::call layer_groups [dict merge [dict create uid $uid layer $layer] $opts]]
}

# Build (or rebuild) a backend-grouped marker layer. Returns positions shown.
# The annotation categories (ptm var site mutag) go through the checked-row state (::ufv::annot below): the rows
# matching the options are checked and the layer is drawn from them. Other layers (dom) keep the backend groups.
proc ::ufv::marker_layer {mol name layer {opts {}} {style {VDW 0.6 12}}} {
    set cat [::ufv::annot::cat_of $name]
    if {$cat in {ptm var site mutag}} { return [::ufv::annot::apply_opts $mol $cat $opts] }
    set uid [::ufv::mol_uid $mol]
    set g [::ufv::groups $uid $layer $opts]
    return [::ufv::marker_apply $mol $name $layer $opts $style $g]
}
proc ::ufv::marker_apply {mol name layer opts style g} {
    set rebuild [list ::ufv::marker_layer $mol $name $layer $opts $style]
    ::ufv::st::put $mol [list shown $name] [dict get $g groups]
    return [::ufv::layers::markers $mol $name [dict get $g groups] -style $style -rebuild $rebuild]
}
# Non-blocking form for the GUI: ::ufv::marker_async mol name layer opts cb   (cb ok npositions)
proc ::ufv::marker_async {mol name layer opts cb {style {VDW 0.6 12}}} {
    set cat [::ufv::annot::cat_of $name]
    if {$cat in {ptm var site mutag}} {
        after 0 [list ::ufv::_marker_local $mol $cat $opts $cb]
        return
    }
    set uid [::ufv::mol_uid $mol]
    ::ufv::backend::async layer_groups [dict merge [dict create uid $uid layer $layer] $opts] \
        [list ::ufv::_marker_got $mol $name $layer $opts $style $cb]
}
proc ::ufv::_marker_got {mol name layer opts style cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; return }
    if {[catch { ::ufv::marker_apply $mol $name $layer $opts $style $res } n]} { {*}$cb 0 $n; return }
    {*}$cb 1 $n
}
proc ::ufv::_marker_local {mol cat opts cb} {
    if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; return }
    if {[catch { ::ufv::annot::apply_opts $mol $cat $opts } n]} { {*}$cb 0 $n; return }
    {*}$cb 1 $n
}

proc ufv_ptms {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::marker_layer $mol ptm ptms]
    ::ufv::status "mol $mol: PTM layer on ($n positions)."
    return $n
}

# Variant options from the command line: consequence words, -provenance, -disease, -color, -reviewed, -only.
proc ::ufv::variant_opts {argl} {
    set molid ""; set only {}; set opts [dict create]
    for {set i 0} {$i < [llength $argl]} {incr i} {
        set a [lindex $argl $i]
        switch -- $a {
            -provenance - -disease - -color { dict set opts [string range $a 1 end] [lindex $argl [incr i]] }
            -only { lappend only {*}[string map {, " "} [lindex $argl [incr i]]] }
            -reviewed { dict set opts reviewed 1 }
            pathogenic - benign - uncertain - deleterious - all { lappend only $a }
            default { set molid $a }
        }
    }
    if {[llength $only]} { dict set opts only [join $only " "] }
    return [list $molid $opts]
}

# ufv_variants ?pathogenic|deleterious|uncertain|benign|all ...? ?-provenance reviewed|clinvar|large?
#              ?-disease "A;B"|none? ?-color consequence|disease|provenance? ?-reviewed? ?molid?
# = check exactly the matching variant rows (the panel's list follows).
proc ufv_variants {args} {
    lassign [::ufv::variant_opts $args] molid opts
    set mol [::ufv::resolve $molid]
    set n [::ufv::marker_layer $mol var variants $opts]
    set filt [expr {[dict exists $opts only] ? [string map {" " ,} [dict get $opts only]] : "pathogenic"}]
    set extra ""
    foreach k {provenance disease color} { if {[dict exists $opts $k]} { append extra "; $k [dict get $opts $k]" } }
    if {[dict exists $opts reviewed]} { append extra "; reviewed only" }
    ::ufv::status "mol $mol: variant layer on ($n positions, $filt$extra)."
    return $n
}

proc ufv_sites {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::marker_layer $mol site sites]
    ::ufv::status "mol $mol: site layer on ($n positions)."
    return $n
}

proc ufv_mutagenesis {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::marker_layer $mol mutag mutagenesis]
    ::ufv::status "mol $mol: mutagenesis layer on ($n positions)."
    return $n
}

# ---------------------------------------------------------------------------------------------
# Annotations v3: the checked rows ARE the 3D. One checked-row set per (molecule, accession, category); every
# category is ONE layer (one rep per colour), rebuilt from its checked set in one batched call. Console commands
# check / uncheck rows (ufv_ptms = all PTM rows, ufv_variants = exactly the matching variant rows, ufv_hide = none,
# ufv_clear = everything off); the panel's list shows the same state (event `annot molid category`).
#   state: st {annsel UID cat} = dict rowid -> 1   {annprev UID cat} = ids before ufv_hide (ufv_show restores)
#          {annvc UID} = variant colour context {color consequence|disease|provenance diseases *|{names}}
#          {annkb UID} = PDBe-KB kind of the last ufv_pdbekb (colours)
#          partner_rows = partner spheres of every layer {chain resi hex text} (the Partners rows)
#   ids (ufv_core ann_rows, same as the PyMOL panel): <cat>:<index into the annotation list> (ptm var site mutag
#          dom), lig:<CODE>, prt:<chain>:<resi>, kb:<pos>
# ---------------------------------------------------------------------------------------------
namespace eval ::ufv::annot {
    # categories, ids, row texts and the search matcher are ufv_core's (ann_rows / search ops), shared with the
    # PyMOL panel: ptm var site mutag dom lig prt kb
    variable cats {ptm var site mutag dom lig prt kb}
    variable info [dict create \
        ptm   {label PTMs short PTM layer ptm} \
        var   {label Variants short Var layer var} \
        site  {label Sites short Site layer site} \
        mutag {label Mutagenesis short Mut layer mutag} \
        dom   {label Domains short Dom layer ann:dom} \
        lig   {label Ligands short Lig layer lig} \
        prt   {label Partners short Prt layer feat:partner} \
        kb    {label PDBe-KB short KB layer feat:kb}]
    # the annotation list each ann-derived category comes from
    variable annkey {ptm ptms var variants site sites mutag mutagenesis dom domains}
    variable sev [dict create "Likely pathogenic or pathogenic" 0 "Predicted deleterious" 1 "Uncertain significance" 2 \
        "Likely benign or benign" 3]
    variable provcol [dict create "UniProt reviewed" #42a5f5 ClinVar #ab47bc "Large scale studies" #26c6da]
    variable provtok [dict create reviewed "UniProt reviewed" uniprot "UniProt reviewed" clinvar ClinVar \
        large "Large scale studies" largescale "Large scale studies" large_scale "Large scale studies"]
    variable unclassified #9e9e9e
    if {![array exists rowcache]} { array set rowcache {} }
}
proc ::ufv::annot::label {cat} { variable info; return [dict get $info $cat label] }
proc ::ufv::annot::short {cat} { variable info; return [dict get $info $cat short] }
proc ::ufv::annot::layer {cat} { variable info; return [dict get $info $cat layer] }
# category of a layer name / category key / label / alias ("" when it is not an annotation category)
proc ::ufv::annot::cat_of {name} {
    variable info
    set n [string tolower [string trim $name]]
    dict for {c d} $info {
        if {$n eq $c || $n eq [string tolower [dict get $d label]] || $n eq [dict get $d layer]} { return $c }
    }
    # ("dom" alone is the cartoon's domain-point layer; the category is "domains")
    set al {ptms ptm variant var variants var sites site mutagenesis mutag domain dom domains dom ligand lig
        ligands lig partner prt partners prt pdbekb kb pdbe-kb kb binding kb}
    if {[dict exists $al $n]} { return [dict get $al $n] }
    return ""
}
proc ::ufv::annot::_uid {mol} {
    set u [expr {$mol ne "" && [::ufv::st::has $mol] ? [::ufv::st::get $mol uid] : ""}]
    if {$u eq ""} { set u $::ufv::uid }
    return [expr {$u eq "" ? "-" : $u}]
}
proc ::ufv::annot::_key {mol cat} { return [list annsel [_uid $mol] $cat] }
proc ::ufv::annot::checked {mol cat} { return [::ufv::st::get $mol [_key $mol $cat] {}] }
proc ::ufv::annot::ids {mol cat} { return [dict keys [checked $mol $cat]] }
proc ::ufv::annot::count {mol cat} { return [dict size [checked $mol $cat]] }
proc ::ufv::annot::counts {mol} {
    variable cats
    set out [dict create]
    foreach c $cats { dict set out $c [count $mol $c] }
    return $out
}
proc ::ufv::annot::total {mol} {
    set n 0
    dict for {c k} [counts $mol] { incr n $k }
    return $n
}

# ---- row ids (no backend: console commands) and row data the backend cannot know (ligands / partners / PDBe-KB)
proc ::ufv::annot::all_ids {mol cat} {
    variable annkey
    set uid [_uid $mol]
    if {[dict exists $annkey $cat]} {
        if {![info exists ::ufv::ann($uid)]} { return {} }
        set n [llength [dict get $::ufv::ann($uid) [dict get $annkey $cat]]]
        set out {}
        for {set i 0} {$i < $n} {incr i} { lappend out $cat:$i }
        return $out
    }
    switch -- $cat {
        lig { return [lmap c [dict keys [lig_codes $mol]] { format lig:%s $c }] }
        prt { return [lmap s [::ufv::st::get $mol partner_rows {}] { format prt:%s:%s [lindex $s 0] [lindex $s 1] }] }
        kb {
            if {![info exists ::ufv::feat::kb($uid)]} { return {} }
            set r $::ufv::feat::kb($uid)
            return [lmap p [lsort -integer -unique [concat [dict get $r ligand] [dict get $r interface]]] { format kb:%s $p }]
        }
    }
    return {}
}
# ligand codes of a molecule -> copies (sorted by code, as the PyMOL plugin's enumerate_ligands)
proc ::ufv::annot::lig_codes {mol} {
    if {$mol eq "" || $mol ni [molinfo list]} { return {} }
    set d [dict create]
    foreach c [::ufv::ligand_copies $mol] { dict incr d [string toupper [lindex $c 0]] }
    set out [dict create]
    foreach k [lsort [dict keys $d]] { dict set out $k [dict get $d $k] }
    return $out
}
# Extra rows for the ann_rows / search ops: {lig [{resn count name}] prt [{chain resi color text}] kb true}
# (prt / kb only once their data is here: partner_rows / ::ufv::feat::kb).
proc ::ufv::annot::extra_json {mol} {
    set parts {}
    set L {}
    dict for {code n} [lig_codes $mol] {
        set nm ""
        if {[info exists ::ufv::feat::liginfo($code)]} { catch { set nm [dict get $::ufv::feat::liginfo($code) name] } }
        lappend L [::ufv::json::obj [dict create resn [::ufv::json::raw [::ufv::json::str $code]] count $n name [::ufv::json::raw [::ufv::json::str $nm]]]]
    }
    lappend parts "\"lig\":\[[join $L ,]\]"
    if {$mol ne "" && [::ufv::st::has $mol] && [::ufv::st::get $mol partner_rows "-"] ne "-"} {
        set P {}
        foreach s [::ufv::st::get $mol partner_rows {}] {
            lassign $s ch ri hex txt
            lappend P "{\"chain\":[::ufv::json::str $ch],\"resi\":[::ufv::json::str $ri],\"color\":[::ufv::json::str $hex],\"text\":[::ufv::json::str $txt]}"
        }
        lappend parts "\"prt\":\[[join $P ,]\]"
    }
    if {[info exists ::ufv::feat::kb([_uid $mol])]} { lappend parts "\"kb\":true" }
    return "{[join $parts ,]}"
}
# Signature of the extra data (the panel rebuilds its list when it changes).
proc ::ufv::annot::extra_sig {mol} {
    set kb [info exists ::ufv::feat::kb([_uid $mol])]
    set prt [expr {$mol ne "" && [::ufv::st::has $mol] ? [llength [::ufv::st::get $mol partner_rows "-"]] : "-"}]
    return [list [_uid $mol] $mol [lig_codes $mol] $prt $kb]
}
# Rows of the panel's list (backend ann_rows, so the texts are the PyMOL panel's): cb ok {row dicts}
#   row: id cat pos end poslabel label color pts
proc ::ufv::annot::rows_async {mol cb} {
    set a [dict create uid [_uid $mol] cat all extra [::ufv::json::raw [extra_json $mol]]]
    ::ufv::backend::async ann_rows $a $cb
}
# The shared matcher: cb ok {ids count total suggest}
proc ::ufv::annot::search_async {mol cat query cb {suggest 1}} {
    set a [dict create uid [_uid $mol] cat $cat query [::ufv::json::raw [::ufv::json::str $query]] \
        extra [::ufv::json::raw [extra_json $mol]] suggest [::ufv::json::raw [expr {$suggest ? "true" : "false"}]]]
    ::ufv::backend::async search $a $cb
}
proc ::ufv::annot::search {mol cat query} {
    set a [dict create uid [_uid $mol] cat $cat query [::ufv::json::raw [::ufv::json::str $query]] \
        extra [::ufv::json::raw [extra_json $mol]] suggest [::ufv::json::raw false]]
    return [dict get [::ufv::backend::call search $a] ids]
}

# Data a category needs before it has rows (network): PDBe-KB (per accession), partner residues (per molecule).
proc ::ufv::annot::needs_load {mol cat} {
    switch -- $cat {
        kb  { return [expr {![info exists ::ufv::feat::kb([_uid $mol])]}] }
        prt { return [expr {$mol ne "" && [::ufv::st::has $mol] && [::ufv::st::get $mol partner_rows "-"] eq "-"}] }
    }
    return 0
}
# cb ok nrows (fetches only; nothing is drawn)
proc ::ufv::annot::load_async {mol cat cb} {
    switch -- $cat {
        kb {
            set uid [_uid $mol]
            if {![info exists ::ufv::feat::kb($uid)]} {
                ::ufv::backend::async pdbekb [dict create uid $uid] [list ::ufv::annot::_kb_loaded $mol $uid $cb]
                return
            }
        }
        prt {
            if {[needs_load $mol prt]} {
                set a [::ufv::feat::ctx $mol 0]
                dict set a layers all
                ::ufv::backend::async partner_layers $a [list ::ufv::annot::_prt_loaded $mol $cb]
                return
            }
        }
    }
    after 0 [list {*}$cb 1 [llength [all_ids $mol $cat]]]
}
proc ::ufv::annot::_kb_loaded {mol uid cb ok res} {
    if {$ok} { set ::ufv::feat::kb($uid) $res }
    {*}$cb $ok [expr {$ok ? [llength [all_ids $mol kb]] : $res}]
}
proc ::ufv::annot::_prt_loaded {mol cb ok res} {
    if {$ok && $mol in [molinfo list] && [::ufv::st::has $mol]} {
        ::ufv::st::put $mol partner_rows [dict get $res spheres]
        if {[::ufv::st::get $mol partner_spheres "-"] eq "-"} { ::ufv::st::put $mol partner_spheres [dict get $res spheres] }
    }
    {*}$cb $ok [expr {$ok ? [llength [all_ids $mol prt]] : $res}]
}

# ---- state changes (each = ONE batched layer rebuild per category + event `annot molid cat`)
proc ::ufv::annot::put {mol cat ids {emit 1}} {
    if {[::ufv::st::has $mol]} { ::ufv::st::remove $mol [list annhid [_uid $mol] $cat] }
    set d [dict create]
    foreach id $ids { dict set d $id 1 }
    ::ufv::st::put $mol [_key $mol $cat] $d
    set n [render $mol $cat]
    if {$emit} { ::ufv::emit annot $mol $cat }
    return $n
}
# check (on 1) / uncheck (on 0) some rows; the others keep their state. changes = {cat {ids on} ...}
proc ::ufv::annot::change {mol changes} {
    set cats {}
    ::ufv::layers::batch {
        dict for {cat spec} $changes {
            lassign $spec ids on
            if {[::ufv::st::has $mol]} { ::ufv::st::remove $mol [list annhid [_uid $mol] $cat] }
            set d [checked $mol $cat]
            if {$on} { foreach id $ids { dict set d $id 1 } } else { foreach id $ids { dict unset d $id } }
            ::ufv::st::put $mol [_key $mol $cat] $d
            render $mol $cat
            lappend cats $cat
        }
    }
    foreach cat $cats { ::ufv::emit annot $mol $cat }
}
# ufv_hide: uncheck every row of the category. The layer's reps are only hidden (state annhid = the rows + their
# groups), so ufv_show brings the same rows back without a rebuild when nothing changed meanwhile (colours,
# numbering, frame-independent marker selections); any other change of the category rebuilds / deletes it.
proc ::ufv::annot::hide {mol cat} {
    set ids [ids $mol $cat]
    if {[llength $ids]} { ::ufv::st::put $mol [list annprev [_uid $mol] $cat] $ids }
    set lname [layer $cat]
    if {[llength $ids] && $cat ni {lig prt} && [::ufv::layers::exists $mol $lname]} {
        ::ufv::st::put $mol [_key $mol $cat] {}
        ::ufv::st::put $mol [list annhid [_uid $mol] $cat] [list $ids [::ufv::st::get $mol [list shown $lname] {}]]
        ::ufv::st::remove $mol [list shown $lname]
        ::ufv::layers::show $mol $lname 0
        ::ufv::emit annot $mol $cat
        return 0
    }
    return [put $mol $cat {}]
}
proc ::ufv::annot::unhide {mol cat} {
    set prev [::ufv::st::get $mol [list annprev [_uid $mol] $cat] {}]
    if {![llength $prev]} { return 0 }
    set lname [layer $cat]
    set hid [::ufv::st::get $mol [list annhid [_uid $mol] $cat] {}]
    if {$hid ne "" && [lindex $hid 0] eq $prev && [::ufv::layers::exists $mol $lname] && ![::ufv::layers::visible $mol $lname]             && [lindex $hid 1] ne "" && [groups $mol $cat $prev] eq [lindex $hid 1]} {
        ::ufv::st::remove $mol [list annhid [_uid $mol] $cat]
        set d [dict create]
        foreach id $prev { dict set d $id 1 }
        ::ufv::st::put $mol [_key $mol $cat] $d
        ::ufv::st::put $mol [list shown $lname] [lindex $hid 1]
        ::ufv::layers::show $mol $lname 1
        ::ufv::emit annot $mol $cat
        return 1
    }
    put $mol $cat $prev
    return 1
}
# Uncheck every row of every category (one batch; cartoon colouring and camera untouched).
proc ::ufv::annot::reset {mol} {
    variable cats
    ::ufv::layers::batch {
        foreach c $cats {
            if {[count $mol $c]} {
                put $mol $c {} 0
            } elseif {[::ufv::layers::exists $mol [layer $c]]} {
                ::ufv::layers::delete $mol [layer $c]
            }
        }
    }
    ::ufv::emit annot $mol *
}
# Console options -> checked rows (ptm / site / mutag: every row; var: the filter of ufv_variants).
proc ::ufv::annot::apply_opts {mol cat opts} {
    set uid [::ufv::mol_uid $mol]
    set a [::ufv::annotations $uid]
    if {$cat eq "var"} {
        set col [string tolower [expr {[dict exists $opts color] ? [dict get $opts color] : ""}]]
        if {![regexp {^(con|dis|pro)} $col]} { set col consequence }
        set vc [dict create color $col diseases [expr {[dict exists $opts disease] ? [parse_dis $a [dict get $opts disease]] : "*"}]]
        ::ufv::st::put $mol [list annvc $uid] $vc
        return [put $mol var [var_match $a $opts]]
    }
    return [put $mol $cat [all_ids $mol $cat]]
}
# Variant colour profile of the panel (consequence | disease | provenance), redraws the checked variants.
proc ::ufv::annot::var_color {mol profile} {
    set uid [_uid $mol]
    set vc [::ufv::st::get $mol [list annvc $uid] [dict create diseases *]]
    dict set vc color $profile
    ::ufv::st::put $mol [list annvc $uid] $vc
    if {[count $mol var]} { render $mol var }
}

# ---- the old variant filters (data.js filterVariants / ufv_core._filter_variants), exactly
proc ::ufv::annot::parse_prov {spec} {
    variable provcol; variable provtok
    if {$spec in {"" all}} { return "*" }
    set out {}
    foreach tok [regexp -all -inline {[^;,|\s]+} $spec] {
        set t [string tolower $tok]
        if {[dict exists $provtok $t]} {
            set c [dict get $provtok $t]
            if {$c ni $out} { lappend out $c }
            continue
        }
        foreach c [dict keys $provcol] { if {[string first $t [string tolower $c]] == 0 && $c ni $out} { lappend out $c } }
    }
    return $out
}
proc ::ufv::annot::parse_dis {a spec} {
    if {$spec in {"" all}} { return "*" }
    if {[string tolower [string trim $spec]] eq "none"} { return {} }
    set known [dict create]
    foreach d [dict keys [dict get $a diseaseColors]] { dict set known [string tolower $d] $d }
    set out {}
    foreach tok [split $spec ";|"] {
        set t [string trim $tok]
        if {$t eq ""} continue
        set d [expr {[dict exists $known [string tolower $t]] ? [dict get $known [string tolower $t]] : $t}]
        if {$d ni $out} { lappend out $d }
    }
    return $out
}
proc ::ufv::annot::var_match {a opts} {
    variable sev; variable provcol
    set o [string tolower [string trim [expr {[dict exists $opts only] ? [dict get $opts only] : "pathogenic"}]]]
    set toks {}
    if {$o ni {"" all}} { set toks [regexp -all -inline {[^ ,]+} $o] }
    set provs [expr {[dict exists $opts provenance] ? [parse_prov [dict get $opts provenance]] : "*"}]
    if {[dict exists $opts reviewed] && [dict get $opts reviewed] ni {0 "" false no}} {
        set provs [expr {$provs eq "*" || "UniProt reviewed" in $provs ? [list "UniProt reviewed"] : {}}]
    }
    set dis [expr {[dict exists $opts disease] ? [parse_dis $a [dict get $opts disease]] : "*"}]
    set out {}; set i -1
    foreach v [dict get $a variants] {
        incr i
        set c [dict get $v consequence]
        if {[llength $toks] && [dict exists $sev $c] && [dict get $v tok] ni $toks} continue
        set pv [dict get $v provenance]
        if {$provs ne "*" && [dict exists $provcol $pv] && $pv ni $provs} continue
        if {$dis ne "*"} {
            set hit 0
            foreach d [dict get $v diseases] { if {$d in $dis} { set hit 1; break } }
            if {!$hit} continue
        }
        lappend out var:$i
    }
    return $out
}

# ---- drawing (one layer per category, one rep per colour)
proc ::ufv::annot::_idx {ids} {
    set out {}
    foreach id $ids { lappend out [string range $id [string first : $id]+1 end] }
    return [lsort -integer $out]
}
proc ::ufv::annot::render {mol cat} {
    if {$mol ni [molinfo list] || ![::ufv::st::has $mol]} { return 0 }
    set d [checked $mol $cat]
    set lname [layer $cat]
    if {$cat eq "lig"} {
        if {![dict size $d]} { if {[::ufv::layers::exists $mol lig]} { ::ufv::layers::delete $mol lig }; return 0 }
        return [::ufv::ligand_layer $mol [lmap k [dict keys $d] { string range $k 4 end }]]
    }
    if {$cat eq "prt"} { return [_render_partners $mol [dict keys $d]] }
    if {![dict size $d]} {
        ::ufv::st::remove $mol [list shown $lname]
        if {[::ufv::layers::exists $mol $lname]} { ::ufv::layers::delete $mol $lname }
        return 0
    }
    set g [groups $mol $cat [dict keys $d]]
    ::ufv::st::put $mol [list shown $lname] $g
    return [::ufv::layers::markers $mol $lname $g -style {VDW 0.6 12} -rebuild [list ::ufv::annot::render $mol $cat]]
}
# {hex {positions}} of checked rows (the rows' pts / colours, as ufv_core.ann_rows and the PyMOL panel)
proc ::ufv::annot::groups {mol cat ids} {
    variable unclassified; variable annkey
    set uid [_uid $mol]
    set g [dict create]
    if {$cat eq "kb"} { return [_kb_groups $mol $uid $ids] }
    if {![info exists ::ufv::ann($uid)]} { return $g }
    set a $::ufv::ann($uid)
    if {$cat eq "var"} { return [_var_groups $mol $a [_idx $ids]] }
    set list [dict get $a [dict get $annkey $cat]]
    foreach i [_idx $ids] {
        set x [lindex $list $i]
        if {$x eq ""} continue
        switch -- $cat {
            ptm - site {
                set hex [dict get $x color]; if {$hex eq "" && $cat eq "site"} { set hex #ffb300 }
                set p [dict get $x pos]; set e [dict get $x end]
                dict lappend g $hex $p
                if {$e != $p} { dict lappend g $hex $e }
            }
            mutag {
                set hex [dict get $x color]; if {$hex eq ""} { set hex $unclassified }
                for {set p [dict get $x pos]} {$p <= [dict get $x end]} {incr p} { dict lappend g $hex $p }
            }
            dom {
                set s [dict get $x start]
                set e [expr {[dict get $x isRange] ? [dict get $x end] : $s}]
                for {set p $s} {$p <= $e} {incr p} { dict lappend g [dict get $x color] $p }
            }
        }
    }
    return $g
}
# one sphere per position, coloured by its most severe checked variant under the colour profile
proc ::ufv::annot::_var_groups {mol a idx} {
    variable sev; variable provcol; variable unclassified
    set vc [::ufv::st::get $mol [list annvc [_uid $mol]] {}]
    set color [string tolower [expr {[dict exists $vc color] ? [dict get $vc color] : "consequence"}]]
    set dset [expr {[dict exists $vc diseases] ? [dict get $vc diseases] : "*"}]
    set anyd [expr {$dset eq "*" || ![llength $dset]}]
    set dcol [dict get $a diseaseColors]
    set vars [dict get $a variants]
    set pc [dict create]; set ps [dict create]
    foreach i $idx {
        set v [lindex $vars $i]
        if {$v eq ""} continue
        set p [dict get $v pos]
        set c [dict get $v consequence]
        set s [expr {[dict exists $sev $c] ? [dict get $sev $c] : -1}]
        if {[dict exists $pc $p]} {
            if {$s < 0} continue
            set s0 [dict get $ps $p]
            if {$s0 >= 0 && $s >= $s0} continue
        }
        if {[string match dis* $color]} {
            set hex $unclassified
            foreach dn [dict get $v diseases] {
                if {$anyd || $dn in $dset} {
                    if {[dict exists $dcol $dn] && [dict get $dcol $dn] ne ""} { set hex [dict get $dcol $dn] }
                    break
                }
            }
        } elseif {[string match prov* $color]} {
            set pv [dict get $v provenance]
            set hex [expr {[dict exists $provcol $pv] ? [dict get $provcol $pv] : $unclassified}]
        } else {
            set hex [dict get $v color]
            if {$hex eq ""} { set hex $unclassified }
        }
        dict set pc $p $hex; dict set ps $p $s
    }
    set g [dict create]
    dict for {p hex} $pc { dict lappend g $hex $p }
    return $g
}
proc ::ufv::annot::_kb_groups {mol uid ids} {
    if {![info exists ::ufv::feat::kb($uid)]} { return {} }
    set r $::ufv::feat::kb($uid)
    set kind [::ufv::st::get $mol [list annkb $uid] both]
    set L [dict create]; foreach p [dict get $r ligand] { dict set L $p 1 }
    set I [dict create]; foreach p [dict get $r interface] { dict set I $p 1 }
    set fixed [expr {$kind in {ligand ligands} ? "ligand" : ($kind in {interface interfaces ppi} ? "interface" : "")}]
    set g [dict create]
    foreach p [_idx $ids] {
        if {$fixed ne ""} {
            set c $fixed
        } else {
            set c [expr {[dict exists $L $p] && [dict exists $I $p] ? "both" : ([dict exists $L $p] ? "ligand" : "interface")}]
        }
        dict lappend g [dict get $::ufv::kb_colors $c] $p
    }
    return $g
}
# partner rows prt:<chain>:<resi> -> spheres (colour of the last ufv_partner_layers set, else the row's colour)
proc ::ufv::annot::_render_partners {mol keys} {
    set want [dict create]; foreach k $keys { dict set want $k 1 }
    set cols [dict create]
    foreach s [::ufv::st::get $mol partner_rows {}] { dict set cols "prt:[lindex $s 0]:[lindex $s 1]" [lindex $s 2] }
    foreach s [::ufv::st::get $mol partner_spheres {}] { dict set cols "prt:[lindex $s 0]:[lindex $s 1]" [lindex $s 2] }
    set by [dict create]; set n 0
    dict for {k hex} $cols {
        if {![dict exists $want $k]} continue
        lassign [split $k :] - ch ri
        dict lappend by [::ufv::palette::id $hex] $ch $ri
        incr n
    }
    if {![dict size $by]} {
        if {[::ufv::layers::exists $mol feat:partner]} { ::ufv::layers::delete $mol feat:partner }
        return 0
    }
    set groups {}
    dict for {cid ks} $by { lappend groups [list $cid "name CA and [::ufv::feat::keysel $ks]"] }
    ::ufv::layers::define $mol feat:partner $groups -style {VDW 0.6 12} -kind marker -rebuild [list ::ufv::annot::render $mol prt]
    return $n
}

# ufv_reset_annotations ?molid?  - uncheck every annotation row of every category (cartoon colouring + view kept)
proc ufv_reset_annotations {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::annot::total $mol]
    ::ufv::annot::reset $mol
    ::ufv::status "mol $mol: $n annotation row(s) unchecked in every category."
    return $n
}
# ufv_show_annotations ?molid? category ?ids? ?mode? ?-search query?
#   Script the annotation list: check exactly (mode set), add or remove rows of one category (ptms variants sites
#   mutagenesis domains ligands partners pdbekb | all). ids = all | none | row ids / positions separated by space,
#   '+', ';' or ',' (var:12 or 12 = row 12 of that category; ligand codes for ligands); -search <query> = the
#   panel's search syntax (e.g. "cons:pathogenic src:clinvar 1300-1500"). Returns {cat checked ..}.
proc ufv_show_annotations {args} {
    set search ""
    if {[set i [lsearch -exact $args -search]] >= 0} {
        set search [lindex $args $i+1]
        set args [lreplace $args $i $i+1]
    }
    set molid ""
    if {[llength $args] && ([string is integer -strict [lindex $args 0]] || [lindex $args 0] eq "top")} {
        set molid [lindex $args 0]; set args [lrange $args 1 end]
    }
    lassign $args category ids mode
    if {$category eq ""} { error "usage: ufv_show_annotations ?molid? category ?ids? ?mode? ?-search query?" }
    if {$ids eq ""} { set ids all }
    if {$mode eq ""} { set mode set }
    set mol [::ufv::resolve $molid]
    set cat [expr {[string tolower $category] eq "all" ? "all" : [::ufv::annot::cat_of $category]}]
    if {$cat eq ""} { error "unknown annotation category '$category' (ptms variants sites mutagenesis domains ligands partners pdbekb all)" }
    set cats [expr {$cat eq "all" ? $::ufv::annot::cats : [list $cat]}]
    set known [dict create]
    foreach c $cats { foreach id [::ufv::annot::all_ids $mol $c] { dict set known $id 1 } }
    set spec [string trim $ids]
    if {$search ne ""} {
        set want [::ufv::annot::search $mol $cat $search]
    } elseif {[string tolower $spec] in {"" all}} {
        set want [dict keys $known]
    } elseif {[string tolower $spec] eq "none"} {
        set want {}
    } else {
        set want {}
        foreach t [regexp -all -inline {[^\s+;,]+} $spec] {
            if {[string first : $t] < 0 && $cat ne "all"} { set t "$cat:[expr {$cat eq "lig" ? [string toupper $t] : $t}]" }
            if {[dict exists $known $t]} { lappend want $t } else { puts "\[UFV\] no row '$t' in $cat." }
        }
    }
    set mode [string tolower $mode]
    set out [dict create]
    ::ufv::layers::batch {
        foreach c $cats {
            set sub [lsearch -all -inline -glob $want $c:*]
            if {$mode eq "set"} {
                ::ufv::annot::put $mol $c $sub
            } else {
                ::ufv::annot::change $mol [dict create $c [list $sub [expr {![string match rem* $mode]}]]]
            }
            dict set out $c [::ufv::annot::count $mol $c]
        }
    }
    ::ufv::status "mol $mol: [join [lmap c $cats { format "%s %d" [::ufv::annot::short $c] [dict get $out $c] }] {, }]."
    return $out
}

# ---------------------------------------------------------------------------------------------
# Cartoon colourings (exclusive)
# ---------------------------------------------------------------------------------------------
# Categorical cartoon mode from backend groups (domains / topology / alphamissense / burden).
proc ::ufv::cartoon_layer {mol mode layer} {
    if {[::ufv::cartoon::mode $mol] eq $mode && [::ufv::layers::exists $mol "cart:$mode"]} { return -1 }
    if {[::ufv::cartoon::reshow $mol $mode]} { return -1 }
    set uid [::ufv::mol_uid $mol]
    set g [::ufv::groups $uid $layer]
    set rebuild [list ::ufv::cartoon_rebuild $mol $mode $layer]
    set n [::ufv::cartoon::categorical $mol $mode [dict get $g groups] [dict get $g base] $rebuild]
    if {$layer in {domains domain}} {
        ::ufv::layers::markers $mol dom [dict get $g points] -rebuild [list ::ufv::marker_layer $mol dom domains {} {VDW 0.7 12}] \
            -style {VDW 0.7 12}
    }
    return [dict get $g count]
}
# rebuild after renumbering: rebuild the data, but keep the mode shown only if it is the current one
proc ::ufv::cartoon_rebuild {mol mode layer} {
    set cur [::ufv::cartoon::mode $mol]
    set uid [::ufv::mol_uid $mol]
    set g [::ufv::groups $uid $layer]
    ::ufv::cartoon::categorical $mol $mode [dict get $g groups] [dict get $g base] [list ::ufv::cartoon_rebuild $mol $mode $layer]
    if {$cur ne $mode} {
        if {$cur eq "none"} { ::ufv::cartoon::plain $mol } elseif {![::ufv::cartoon::reshow $mol $cur]} { ::ufv::cartoon::plain $mol }
    }
}
proc ufv_domains {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::cartoon_layer $mol domains domains]
    ::ufv::status "mol $mol: cartoon coloured by domains / regions."
    return $n
}
proc ufv_topology {{molid ""}} {
    set mol [::ufv::resolve $molid]
    ::ufv::cartoon_layer $mol topology topology
    ::ufv::status "mol $mol: cartoon coloured by membrane topology."
}
proc ufv_alphamissense {{molid ""}} {
    set mol [::ufv::resolve $molid]
    if {![dict size [dict get [::ufv::annotations [::ufv::mol_uid $mol]] am]]} {
        ::ufv::status "no AlphaMissense scores."; return 0
    }
    ::ufv::cartoon_layer $mol alphamissense alphamissense
    ::ufv::status "mol $mol: cartoon coloured by mean AlphaMissense."
}
proc ufv_burden {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set b [dict get [::ufv::annotations [::ufv::mol_uid $mol]] burden]
    if {![llength $b]} { ::ufv::status "no burden-positive residues."; return 0 }
    ::ufv::cartoon_layer $mol burden burden
    ::ufv::status "mol $mol: mutation burden ([llength $b] residues)."
}

# pLDDT: the AlphaFold DB bands (categorical, exact AFDB colours) on the B-factor column.
proc ufv_plddt {{molid ""}} {
    set mol [::ufv::resolve $molid]
    if {![::ufv::cartoon::reshow $mol plddt]} {
        ::ufv::layers::batch {
            set b [::ufv::cartoon::base_idx $mol]
            ::ufv::cartoon::_hide_modes $mol plddt
            set groups [list \
                [list [::ufv::palette::id #0053d6] "protein and beta >= 90"] \
                [list [::ufv::palette::id #65cbf3] "protein and beta >= 70 and beta < 90"] \
                [list [::ufv::palette::id #ffdb13] "protein and beta >= 50 and beta < 70"] \
                [list [::ufv::palette::id #ff7d45] "protein and beta < 50"]]
            ::ufv::layers::define $mol cart:plddt $groups -style $::ufv::cartoon::style -kind cartoon
            mol modselect $b $mol "none"
        }
        ::ufv::st::put $mol cartoon [dict create mode plddt kind categorical basesel none basehex ""]
        ::ufv::st::put $mol {modes plddt} [dict create basesel none basehex ""]
        ::ufv::emit cartoon $mol plddt
    }
    ::ufv::status "mol $mol: coloured by pLDDT (AlphaFold DB bands)."
}

# B-factor (experimental): continuous colour scale on Beta.
proc ufv_bfactor {{molid ""}} {
    set mol [::ufv::resolve $molid]
    ::ufv::cartoon::continuous $mol bfactor Beta
    ::ufv::status "mol $mol: coloured by B-factor (colour scale, auto range)."
}

# ufv_cartoon plain|<mode>   (plain = remove the cartoon colouring)
proc ufv_cartoon {{mode plain} {molid ""}} {
    set mol [::ufv::resolve $molid]
    switch -- $mode {
        plain - none - off { ::ufv::cartoon::plain $mol }
        domains { ufv_domains $mol }  topology { ufv_topology $mol }  alphamissense - am { ufv_alphamissense $mol }
        burden { ufv_burden $mol }  plddt { ufv_plddt $mol }  bfactor { ufv_bfactor $mol }
        hotspots { ufv_hotspots $mol }  hubs - contacthubs { ufv_contacthubs $mol }  pockets { ufv_pockets $mol }
        chimera { ufv_chimera $mol }  allpockets { ufv_pockets_all $mol }
        default { if {![::ufv::cartoon::reshow $mol $mode]} { error "unknown / unbuilt cartoon mode '$mode'" } }
    }
    return [::ufv::cartoon::mode $mol]
}

# ---------------------------------------------------------------------------------------------
# Structure analyses: hotspots / contact hubs / constraint pockets (async, cached per uid + numbering + frame).
# One backend request computes the three (as the PyMOL plugin's _analysis, per virtual object + geometry):
# hotspots fold in partner-protein disease residues of a complex, pockets are PAE-gated on AlphaFold models
# and keep every candidate's q, so ufv_pocket_q recolours without recomputing.
# ---------------------------------------------------------------------------------------------
proc ::ufv::export_ca {mol path} {
    set items {}
    ::ufv::with_sel s $mol "$::ufv::poly_base and name CA" {
        set seen [dict create]
        foreach {ch ri ins x y z b} [join [$s get {chain resid insertion x y z beta}]] {
            if {[dict exists $seen $ch|$ri]} continue
            dict set seen $ch|$ri 1
            set u [expr {[string trim $ins] eq "" ? [::ufv::map::r2u $mol $ch $ri] : ""}]
            lappend items [format {{"uniPos":%s,"chain":%s,"resi":%d,"x":%.3f,"y":%.3f,"z":%.3f,"b":%.2f}} \
                [expr {$u eq "" ? "null" : $u}] [::ufv::json::str $ch] $ri $x $y $z $b]
        }
    }
    set fh [open $path w]
    puts $fh "\[[join $items ,\n]\]"
    close $fh
    return [llength $items]
}

# {pos cat} of the constraint-pocket candidates with q <= the sensitivity threshold (_pocket_classes).
proc ::ufv::pocket_classes {res {q ""}} {
    if {$q eq ""} { set q $::ufv::pocket_q }
    set out [dict create]
    if {![dict exists $res byPos]} { return $out }
    dict for {p v} [dict get $res byPos] {
        set pq [expr {[dict exists $v q] && [dict get $v q] ne "" ? [dict get $v q] : 1.0}]
        if {$pq <= $q} { dict set out $p [dict get $v cat] }
    }
    return $out
}

# One analysis kind (hotspots | hubs | pockets) of the CURRENT frame -> cb ok result
#   hotspots / hubs: {tiers {pos tier} byChain {chain {pos tier}}}; pockets: {byPos {pos {q cat}} hasPae reason n}
proc ::ufv::analysis_async {mol kind cb} {
    if {$kind ni {hotspots hubs pockets}} { after 0 [list {*}$cb 0 "unknown analysis '$kind'"]; return }
    set uid [::ufv::mol_uid $mol]
    set key [list $uid [::ufv::st::get $mol mapver 0] [molinfo $mol get frame]]
    set c [::ufv::st::get $mol [list cache an:$kind] ""]
    if {$c ne "" && [lindex $c 0] eq $key} { after 0 [list {*}$cb 1 [lindex $c 1]]; return }
    # one request in flight per key: later callers queue on it
    set qk [list $mol $key]
    if {[info exists ::ufv::feat::anq($qk)]} { lappend ::ufv::feat::anq($qk) [list $kind $cb]; return }
    set ::ufv::feat::anq($qk) [list [list $kind $cb]]
    if {[catch { set a [::ufv::feat::ctx $mol 1] } err]} { ::ufv::_analysis_done $mol $key $qk 0 $err; return }
    dict set a kinds "hotspots hubs pockets"
    ::ufv::status "computing hotspots / hubs / pockets ..."
    ::ufv::backend::async analyses $a [list ::ufv::_analysis_done $mol $key $qk]
}
proc ::ufv::_analysis_done {mol key qk ok res} {
    set waiters [expr {[info exists ::ufv::feat::anq($qk)] ? $::ufv::feat::anq($qk) : {}}]
    unset -nocomplain ::ufv::feat::anq($qk)
    if {$ok && $mol in [molinfo list]} {
        foreach kind {hotspots hubs pockets} {
            if {[dict exists $res errors $kind] || ![dict exists $res $kind]} continue
            set r [dict get $res $kind]
            if {$kind ne "pockets"} { set r [dict create tiers $r byChain [dict get $res ${kind}ByChain]] }
            ::ufv::st::put $mol [list cache an:$kind] [list $key $r]
        }
        # legacy key (trajectory module / tests): {uid mapver frame q} -> {hotspots hubs pockets errors}
        set leg [dict create hotspots {} hubs {} pockets {} errors [dict get $res errors]]
        foreach kind {hotspots hubs} { if {[dict exists $res $kind]} { dict set leg $kind [dict get $res $kind] } }
        if {[dict exists $res pockets]} { dict set leg pockets [::ufv::pocket_classes [dict get $res pockets]] }
        ::ufv::st::put $mol {cache tiers} [list [list {*}$key $::ufv::pocket_q] $leg]
        ::ufv::emit tiers $mol
    }
    foreach w $waiters {
        lassign $w kind cb
        if {!$ok} { {*}$cb 0 $res; continue }
        if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; continue }
        if {[dict exists $res errors $kind]} { {*}$cb 0 [dict get $res errors $kind]; continue }
        {*}$cb 1 [lindex [::ufv::st::get $mol [list cache an:$kind]] 1]
    }
}

# Legacy combined form: cb ok {hotspots {pos tier} hubs {pos tier} pockets {pos class (at q)} errors {}}
proc ::ufv::tiers_async {mol q cb} {
    ::ufv::analysis_async $mol hotspots [list ::ufv::_tiers_join $mol $q $cb]
}
proc ::ufv::_tiers_join {mol q cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    set out [dict create errors {}]
    foreach kind {hotspots hubs pockets} {
        set c [lindex [::ufv::st::get $mol [list cache an:$kind] ""] 1]
        if {$kind eq "pockets"} { dict set out pockets [::ufv::pocket_classes $c $q] } \
        else { dict set out $kind [expr {$c eq "" ? {} : [dict get $c tiers]}] }
    }
    {*}$cb 1 $out
}

set ::ufv::tier_colors [dict create \
    hotspots {strong #b71c1c moderate #e64a19 weak #ffa726} \
    hubs     {strong #6a1b9a moderate #ab47bc} \
    pockets  {pocket #00897b exposed #8e24aa}]

# Colour the cartoon from an analysis result (no backend call). Returns the number of residues coloured.
proc ::ufv::tier_apply {mol kind r {q ""}} {
    set tiers [expr {$kind eq "pockets" ? [::ufv::pocket_classes $r $q] : [dict get $r tiers]}]
    set hexgroups [dict create]
    dict for {tier hex} [dict get $::ufv::tier_colors $kind] {
        dict for {pos t} $tiers { if {$t eq $tier} { dict lappend hexgroups $hex $pos } }
    }
    ::ufv::cartoon::categorical $mol $kind $hexgroups "" [list ::ufv::tier_rebuild $mol $kind]
    set ran [::ufv::st::get $mol ran {}]
    if {$kind ni $ran} { ::ufv::st::put $mol ran [lappend ran $kind] }
    return [dict size $tiers]
}
# after renumbering: recompute in the background; re-colour only if the mode is still the one shown
proc ::ufv::tier_rebuild {mol kind} {
    ::ufv::analysis_async $mol $kind [list ::ufv::_tier_rebuilt $mol $kind]
}
proc ::ufv::_tier_rebuilt {mol kind ok res} {
    if {!$ok || $mol ni [molinfo list]} return
    if {[::ufv::cartoon::mode $mol] eq $kind} {
        ::ufv::tier_apply $mol $kind $res
    } elseif {[::ufv::layers::exists $mol cart:$kind]} {
        ::ufv::layers::delete $mol cart:$kind
        ::ufv::st::remove $mol [list modes $kind]
    }
}

# Async colouring for the GUI: ::ufv::tier_async mol hotspots|hubs|pockets cb  (cb ok nresidues)
proc ::ufv::tier_async {mol kind cb} {
    ::ufv::analysis_async $mol $kind [list ::ufv::_tier_got $mol $kind $cb]
}
proc ::ufv::_tier_got {mol kind cb ok res} {
    if {!$ok} { ::ufv::status "$kind failed: $res"; {*}$cb 0 $res; return }
    if {[catch { ::ufv::tier_apply $mol $kind $res } n]} { {*}$cb 0 $n; return }
    {*}$cb 1 $n
}

proc ::ufv::tier_mode {mol kind {q ""}} {
    if {$q ne ""} { ::ufv::pocket_q_set $q }
    set r [::ufv::await [list ::ufv::analysis_async $mol $kind]]
    return [::ufv::tier_apply $mol $kind $r]
}

proc ufv_hotspots {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::tier_mode $mol hotspots]
    ::ufv::status "mol $mol: $n hotspot residues."
    return $n
}
proc ufv_contacthubs {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::tier_mode $mol hubs]
    ::ufv::status "mol $mol: $n contact-hub residues."
    return $n
}
# ufv_pockets ?molid? ?q?  - constraint pockets at sensitivity q (default the current ufv_pocket_q)
proc ufv_pockets {{molid ""} {q ""}} {
    set mol [::ufv::resolve $molid]
    if {![dict size [dict get [::ufv::annotations [::ufv::mol_uid $mol]] am]]} {
        ::ufv::status "no AlphaMissense data - constraint-pocket needs it."; return 0
    }
    set n [::ufv::tier_mode $mol pockets $q]
    set r [lindex [::ufv::st::get $mol {cache an:pockets}] 1]
    ::ufv::status [format "mol %s: %d constraint-pocket residues at q <= %.2f (%d candidates%s)%s." $mol $n \
        $::ufv::pocket_q [dict size [dict get $r byPos]] [expr {[dict get $r hasPae] ? ", PAE-gated" : ""}] \
        [expr {[dict get $r reason] ne "" ? " - [dict get $r reason]" : ""}]]
    return $n
}

# Sensitivity (FDR q, 0.01-0.40). Molecules showing the constraint-pocket colouring are recoloured from the
# cached result (no recompute). Returns the clamped value.
proc ::ufv::pocket_q_set {q {molid ""}} {
    set q [expr {double($q)}]
    if {$q < 0.01} { set q 0.01 } elseif {$q > 0.40} { set q 0.40 }
    set ::ufv::pocket_q [expr {round($q * 10000) / 10000.0}]
    set mols [expr {$molid ne "" ? [list [::ufv::resolve $molid]] : [molinfo list]}]
    foreach m $mols {
        if {![::ufv::st::has $m] || [::ufv::cartoon::mode $m] ne "pockets"} continue
        set c [::ufv::st::get $m {cache an:pockets} ""]
        if {$c eq ""} continue
        set n [::ufv::tier_apply $m pockets [lindex $c 1]]
        ::ufv::status [format "mol %s: %d constraint-pocket residues at q <= %.2f." $m $n $::ufv::pocket_q]
    }
    return $::ufv::pocket_q
}
proc ufv_pocket_q {{q ""} {molid ""}} {
    if {$q eq ""} {
        ::ufv::status [format "constraint-pocket sensitivity q <= %.2f" $::ufv::pocket_q]
        return $::ufv::pocket_q
    }
    return [::ufv::pocket_q_set $q $molid]
}

# ufv_prox ?ptm_radius? ?var_radius?  - PTM-variant proximity radii (A, 2-30; defaults 8 / 12): report lists and
# the PTM-centric CSV tiers.
proc ufv_prox {{ptm ""} {var ""}} {
    foreach {k v} [list ptm $ptm var $var] {
        if {$v ne "" && [string is double -strict $v]} {
            dict set ::ufv::prox $k [expr {min(30.0, max(2.0, double($v)))}]
        }
    }
    ::ufv::status [format "proximity radii: PTMs <= %.1f A, pathogenic variants <= %.1f A" \
        [dict get $::ufv::prox ptm] [dict get $::ufv::prox var]]
    return $::ufv::prox
}

# ---------------------------------------------------------------------------------------------
# Layer visibility
# ---------------------------------------------------------------------------------------------
set ::ufv::layer_alias [dict create ptms ptm ptm ptm variants var variant var var var sites site site site \
    mutagenesis mutag mutag mutag ligands lig ligand lig lig lig focus focus pocket pocket dom dom \
    domains cart:domains topology cart:topology alphamissense cart:alphamissense am cart:alphamissense \
    burden cart:burden plddt cart:plddt hotspots cart:hotspots hubs cart:hubs contacthubs cart:hubs \
    pockets cart:pockets chimera cart:chimera allpockets cart:allpockets bfactor cart:bfactor \
    partner feat:partner partners feat:partner prt feat:partner binding feat:kb kb feat:kb pdbekb feat:kb]

proc ::ufv::layer_name {name} {
    if {[dict exists $::ufv::layer_alias $name]} { return [dict get $::ufv::layer_alias $name] }
    return $name
}

# ufv_hide <layer|cartoon|all> ?molid?  /  ufv_show <layer> ?molid?
# Annotation layers (ptm var site mutag lig feat:partner feat:kb ann:dom): ufv_hide unchecks every row of the
# category (remembered: ufv_show checks them again).
proc ufv_hide {name {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {$name in {cartoon colouring coloring}} { ::ufv::cartoon::plain $mol; return }
    if {$name eq "all"} {
        ::ufv::layers::batch {
            foreach c $::ufv::annot::cats { if {[::ufv::annot::count $mol $c]} { ::ufv::annot::hide $mol $c } }
            foreach n [::ufv::layers::names $mol] { if {$n ne "cartoon"} { ::ufv::layers::show $mol $n 0 } }
        }
        return
    }
    set n [::ufv::layer_name $name]
    set cat [expr {$n eq "dom" ? "" : [::ufv::annot::cat_of $n]}]
    if {$cat ne ""} { ::ufv::annot::hide $mol $cat; return }
    if {[string match cart:* $n]} {
        set mode [string range $n 5 end]
        if {[::ufv::cartoon::mode $mol] eq $mode} { ::ufv::cartoon::plain $mol }
        return
    }
    if {$n eq "focus"} { ::ufv::feat::clear_focus $mol; return }
    if {![::ufv::layers::show $mol $n 0]} { ::ufv::status "no layer '$name' on mol $mol." }
}
proc ufv_show {name {molid ""}} {
    set mol [::ufv::resolve $molid]
    set n [::ufv::layer_name $name]
    if {[string match cart:* $n]} { return [ufv_cartoon [string range $n 5 end] $mol] }
    set cat [expr {$n eq "dom" ? "" : [::ufv::annot::cat_of $n]}]
    if {$cat ne "" && ([::ufv::annot::count $mol $cat] || [::ufv::annot::unhide $mol $cat])} { return }
    if {![::ufv::layers::show $mol $n 1]} { ::ufv::status "no layer '$name' on mol $mol (build it first)." }
}
proc ufv_layers {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set out {}
    foreach n [::ufv::layers::names $mol] {
        lappend out $n [list visible [::ufv::layers::visible $mol $n] reps [::ufv::layers::reps $mol $n]]
    }
    lappend out cartoon_mode [::ufv::cartoon::mode $mol]
    return $out
}

# Remove every UFV layer (only reps the plugin created; your own reps stay).
proc ufv_clear {{molid ""}} {
    set mol [::ufv::resolve $molid]
    ::ufv::st::remove $mol focus
    ::ufv::st::remove $mol shown
    ::ufv::st::remove $mol ran
    foreach k {annsel annprev annhid annvc annkb} { ::ufv::st::remove $mol $k }
    ::ufv::layers::clear $mol
    ::ufv::emit annot $mol *
    ::ufv::status "cleared UFV layers on mol $mol (your own representations are kept)."
}

# ---------------------------------------------------------------------------------------------
# Ligands (organic copies = one licorice layer, ions = one sphere layer; AlphaFill identity filter)
# ---------------------------------------------------------------------------------------------
# {resname chain resid} copies of non-polymer, non-water residues in model order (one `get`).
proc ::ufv::ligand_copies {mol} {
    set out {}
    ::ufv::with_sel s $mol $::ufv::lig_base {
        set seen [dict create]
        foreach {rn ch ri} [join [$s get {resname chain resid}]] {
            if {![dict exists $seen "$rn $ch $ri"]} { dict set seen "$rn $ch $ri" 1; lappend out [list $rn $ch $ri] }
        }
    }
    return $out
}

# AlphaFill donor-identity filter: chains of transplants below ::ufv::aff_min (-> " and not chain ...")
proc ::ufv::feat::aff_low {mol} {
    set meta [::ufv::st::get $mol aff ""]
    set low {}; set total 0
    if {$meta eq ""} { return [list {} 0] }
    dict for {ch m} $meta {
        if {![dict exists $m identity] || [dict get $m identity] eq ""} continue
        incr total
        if {[dict get $m identity] < $::ufv::aff_min} { lappend low $ch }
    }
    return [list [lsort $low] $total]
}

# codes: "*" = every copy, else the ligand codes of the checked ligand rows (lig:<code>). Returns copies drawn.
proc ::ufv::ligand_layer {mol {codes *}} {
    set copies [::ufv::ligand_copies $mol]
    set sub ""
    if {$codes ne "*"} {
        set want [dict create]; foreach k $codes { dict set want [string toupper $k] 1 }
        set rns {}; set n 0
        foreach c $copies {
            set rn [lindex $c 0]
            if {![dict exists $want [string toupper $rn]]} continue
            if {$rn ni $rns} { lappend rns $rn }
            incr n
        }
        if {!$n} { set copies {} } else { set sub " and resname [join [lmap r $rns { list $r }] { }]"; set ncopies $n }
    }
    if {![llength $copies]} {
        if {[::ufv::layers::exists $mol lig]} { ::ufv::layers::delete $mol lig }
        return 0
    }
    set ions {}
    foreach c $copies { if {[lindex $c 0] in $::ufv::ion_codes && [lindex $c 0] ni $ions} { lappend ions [lindex $c 0] } }
    lassign [::ufv::feat::aff_low $mol] low total
    set base "($::ufv::lig_base)"
    if {[llength $low]} {
        set t {}; foreach ch $low { lappend t [::ufv::map::chain_token $ch] }
        append base " and not chain [join $t { }]"
    }
    append base $sub
    set org $base
    if {[llength $ions]} { append org " and not resname [join $ions { }]" }
    set groups [list [list [::ufv::palette::id #4caf50] "$org and element C"] [list Element "$org and not element C"]]
    if {[llength $ions]} { lappend groups [list Element "$base and resname [join $ions { }]" {VDW 0.3 12}] }
    ::ufv::layers::define $mol lig $groups -style {Licorice 0.3 12 12} -kind other -rebuild [list ::ufv::annot::render $mol lig]
    ::ufv::st::put $mol lig [dict create copies [llength $copies] ions $ions low $low total $total]
    return [expr {[info exists ncopies] ? $ncopies : [llength $copies]}]
}

# ::ufv::ligands_async mol cb : AlphaFill transplant metadata fetched first (once) on AlphaFill models
proc ::ufv::ligands_async {mol cb} {
    if {[::ufv::feat::is_alphafill $mol] && [::ufv::st::get $mol aff "-"] eq "-"} {
        ::ufv::backend::async alphafill [dict create uid [::ufv::mol_uid $mol]] [list ::ufv::_aff_got $mol $cb]
        return
    }
    after 0 [list ::ufv::_ligands_build $mol $cb]
}
proc ::ufv::_aff_got {mol cb ok res} {
    ::ufv::st::put $mol aff [expr {$ok ? $res : {}}]
    ::ufv::_ligands_build $mol $cb
}
proc ::ufv::_ligands_build {mol cb} {
    if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; return }
    if {[catch { ::ufv::annot::put $mol lig [::ufv::annot::all_ids $mol lig] } n]} { {*}$cb 0 $n; return }
    if {!$n} { ::ufv::status "no ligands in mol $mol."; {*}$cb 1 0; return }
    set L [::ufv::st::get $mol lig]
    set org 0; foreach c [::ufv::ligand_copies $mol] { if {[lindex $c 0] ni $::ufv::ion_codes} { incr org } }
    set msg "mol $mol: ligands shown ($n copies, $org organic, [llength [dict get $L ions]] ion type(s))"
    if {[dict get $L total]} {
        append msg [format ", AlphaFill transplants %d/%d at identity >= %d%%" [expr {[dict get $L total] - [llength [dict get $L low]]}] \
            [dict get $L total] [expr {round($::ufv::aff_min * 100)}]]
    }
    ::ufv::status "$msg."
    {*}$cb 1 $n
}

# ufv_ligands ?molid? ?uid? ?min_identity?
proc ufv_ligands {{molid ""} {uid ""} {min_identity ""}} {
    set mol [::ufv::resolve $molid]
    if {$uid ne ""} { ::ufv::st::put $mol uid [::ufv::need $uid] }
    if {$min_identity ne ""} { set v [expr {double($min_identity)}]; set ::ufv::aff_min [expr {$v > 1 ? $v / 100.0 : $v}] }
    return [::ufv::await [list ::ufv::ligands_async $mol]]
}
proc ufv_ligands_hide {{molid ""}} {
    set mol [::ufv::resolve $molid]
    ::ufv::annot::hide $mol lig
}

# ufv_ligand_identity pct ?molid?  - AlphaFill: show only transplants with donor identity >= pct % (25 = all)
proc ufv_ligand_identity {{pct 25} {molid ""}} {
    set v [expr {double($pct)}]
    set ::ufv::aff_min [expr {$v > 1 ? $v / 100.0 : $v}]
    set mol [::ufv::resolve $molid]
    if {[::ufv::st::get $mol aff ""] in {"" "-"}} {
        ::ufv::status "mol $mol: no AlphaFill transplant metadata (load the AlphaFill model and run ufv_ligands)."
        return ""
    }
    if {[::ufv::annot::count $mol lig]} { ::ufv::annot::render $mol lig }
    lassign [::ufv::feat::aff_low $mol] low total
    set shown [expr {$total - [llength $low]}]
    ::ufv::status [format "mol %s: %d/%d transplants with donor identity >= %d%%." $mol $shown $total [expr {round($::ufv::aff_min * 100)}]]
    return [list $shown $total]
}

# ---------------------------------------------------------------------------------------------
# Focus / Nearby (current frame: atomselect evaluates `within` at the frame shown when it is created)
# ---------------------------------------------------------------------------------------------
proc ufv_nearby {{dist ""} {molid ""}} {
    if {$dist eq ""} { return $::ufv::nearby }
    set d [expr {double($dist)}]
    if {$d < 3} { set d 3.0 } elseif {$d > 15} { set d 15.0 }
    set ::ufv::nearby $d
    # re-focus the current residue / ligand focus with the new radius
    set mols [expr {$molid ne "" ? [list [::ufv::resolve $molid]] : [molinfo list]}]
    foreach m $mols {
        if {![::ufv::st::has $m]} continue
        set f [::ufv::st::get $m focus ""]
        if {$f eq "" || [dict get $f kind] ni {residue ligand}} continue
        catch {
            if {[dict get $f kind] eq "ligand"} { ::ufv::feat::ligand_focus $m {*}[dict get $f ligand] 0 } \
            else { ::ufv::feat::focus $m [dict get $f pos] 0 }
        }
    }
    ::ufv::status "Nearby distance $d A."
    return $d
}

# Explicit zoom (focus / range / ligand / pocket): deterministic; every molecule keeps sharing one view.
proc ::ufv::zoom_to {mol selstr {buffer 3}} {
    return [::ufv::view::zoom $mol $selstr $buffer]
}

# Residue keys {chain resid ...} of an atom selection text (current frame) -> selection "(chain A and resid ..)"
proc ::ufv::feat::keysel {keys} {
    set by [dict create]
    foreach {ch ri} $keys { dict lappend by $ch $ri }
    set cl {}
    dict for {ch rs} $by { lappend cl "(chain [::ufv::map::chain_token $ch] and resid [::ufv::map::compress $rs])" }
    return [expr {[llength $cl] ? "([join $cl { or }])" : "none"}]
}

# Nearby query for one UniProt position at the current frame (viewer-molstar.js _focusNeighbourhood):
#   residues with any polymer atom within d (transplant chains excluded) -> dict
#   sel (the residue), keys {ch ri ..} (neighbour residues, incl. itself), near (UniProt positions, itself
#   excluded), partners {{chain resid}..} (residues of OTHER-protein chains, mapped-structure modes only)
proc ::ufv::feat::nearby_query {mol pos {d ""}} {
    if {$d eq ""} { set d $::ufv::nearby }
    set s [::ufv::map::sel $mol [list $pos]]
    if {$s eq ""} { return "" }
    set keys {}; set near {}; set partners {}
    set mode [::ufv::map::mode $mol]
    set ours [expr {$mode in {sifts custom manual} ? [::ufv::map::chains $mol] : ""}]
    ::ufv::with_sel n $mol "same residue as ($::ufv::poly_base[notx $mol] and within $d of ($s))" {
        set seen [dict create]
        foreach {ch ri nm} [join [$n get {chain resid name}]] {
            if {![dict exists $seen $ch|$ri]} { dict set seen $ch|$ri 0; lappend keys $ch $ri }
            if {$nm ne "CA" || [dict get $seen $ch|$ri]} continue
            dict set seen $ch|$ri 1
            set u [::ufv::map::r2u $mol $ch $ri]
            if {$u ne ""} { if {$u != $pos && $u ni $near} { lappend near $u } } \
            elseif {$ours ne "" && $ch ni $ours} { lappend partners [list $ch $ri] }
        }
    }
    set partners [lsort -index 1 -integer $partners]
    set partners [lsort -index 0 $partners]
    return [dict create sel $s keys $keys near [lsort -integer $near] partners $partners dist $d]
}

# Hide the visible marker layers while zoomed in (declutter), without touching their visible flag.
proc ::ufv::feat::hide_markers {mol} {
    set hidden {}
    foreach n {ptm var site dom mutag ann:dom} {
        if {[::ufv::layers::visible $mol $n]} {
            foreach i [::ufv::layers::reps $mol $n] { mol showrep $mol $i 0 }
            lappend hidden $n
        }
    }
    return $hidden
}
proc ::ufv::feat::restore_markers {mol names} {
    foreach n $names {
        if {[::ufv::layers::visible $mol $n]} { foreach i [::ufv::layers::reps $mol $n] { mol showrep $mol $i 1 } }
    }
}
proc ::ufv::feat::clear_focus {mol} {
    if {![::ufv::st::has $mol]} return
    set f [::ufv::st::get $mol focus ""]
    ::ufv::layers::batch {
        if {[::ufv::layers::exists $mol focus]} { ::ufv::layers::delete $mol focus }
        if {[::ufv::layers::exists $mol pocket]} { ::ufv::layers::delete $mol pocket }
        if {$f ne ""} {
            restore_markers $mol [dict get $f hidden]
            if {[dict exists $f lighidden] && [dict get $f lighidden] && [::ufv::layers::visible $mol lig]} {
                foreach i [::ufv::layers::reps $mol lig] { mol showrep $mol $i 1 }
            }
        }
    }
    ::ufv::st::remove $mol focus
}

# Licorice groups for neighbour residues coloured by annotation (variant / PTM / site) or grey.
proc ::ufv::feat::anno_groups {mol keys uid {style {Licorice 0.2 12 12}}} {
    set cmap [anno_colors $uid]
    set by [dict create]
    foreach {ch ri} $keys {
        set u [::ufv::map::r2u $mol $ch $ri]
        set hex [expr {$u ne "" && [dict exists $cmap $u] ? [dict get $cmap $u] : "#b3b3b3"}]
        dict lappend by $hex $ch $ri
    }
    set groups {}
    dict for {hex ks} $by { lappend groups [list [::ufv::palette::id $hex] "[keysel $ks] and not hydrogen" $style] }
    return $groups
}

# Residue focus: the residue (ball-and-stick, carbons in its annotation colour) + its Nearby neighbourhood
# (licorice, annotation colours) + pocket ligands within d+3 A when the ligand layer is off. Fixed residue sets
# from the current frame (the PyMOL plugin's named selections). Returns the focus dict ("" if not mapped).
# declutter 0: the marker layers stay shown (the panel's single click: a transient highlight, camera untouched).
proc ::ufv::feat::focus {mol pos {zoom 1} {declutter 1}} {
    set uid [::ufv::mol_uid $mol]
    set q [nearby_query $mol $pos]
    if {$q eq ""} { ::ufv::status "position $pos not mapped on mol $mol."; return "" }
    set s [dict get $q sel]
    set d [dict get $q dist]
    set prev [::ufv::st::get $mol focus ""]
    set hidden [expr {$prev ne "" ? [dict get $prev hidden] : {}}]
    ::ufv::layers::batch {
        if {[::ufv::layers::exists $mol pocket]} { ::ufv::layers::delete $mol pocket }
        if {$declutter} { foreach n [hide_markers $mol] { if {$n ni $hidden} { lappend hidden $n } } }
        set cmap [anno_colors $uid]
        set hex [expr {[dict exists $cmap $pos] ? [dict get $cmap $pos] : "#d9d9d9"}]
        set cid [::ufv::palette::id $hex]
        # neighbours (without the focused residue)
        set own {}
        ::ufv::with_sel r $mol $s { foreach {ch ri} [join [$r get {chain resid}]] { lappend own "$ch|$ri" } }
        set nk {}
        foreach {ch ri} [dict get $q keys] { if {"$ch|$ri" ni $own} { lappend nk $ch $ri } }
        set groups [anno_groups $mol $nk $uid]
        lappend groups [list $cid "($s) and element C and not name CA" {Licorice 0.3 12 12}]
        lappend groups [list Element "($s) and not element C and not hydrogen" {Licorice 0.3 12 12}]
        lappend groups [list $cid "($s) and name CA" {VDW 0.5 12}]
        set ligkeys {}
        if {![::ufv::layers::visible $mol lig]} {
            ::ufv::with_sel l $mol "same residue as ($::ufv::lig_base and within [expr {$d + 3}] of ($s))" {
                foreach {ch ri} [join [$l get {chain resid}]] { lappend ligkeys $ch $ri }
            }
            if {[llength $ligkeys]} { lappend groups [list Element "($::ufv::lig_base) and [keysel $ligkeys]" {Licorice 0.25 12 12}] }
        }
        ::ufv::layers::define $mol focus $groups -style {Licorice 0.2 12 12} -kind other -rebuild [list ::ufv::feat::refocus $mol]
    }
    set f [dict create kind residue pos $pos dist $d near [dict get $q near] partners [dict get $q partners] \
               keys [dict get $q keys] hidden $hidden sel $s mapver [::ufv::st::get $mol mapver 0] frame [molinfo $mol get frame]                declutter $declutter]
    ::ufv::st::put $mol focus $f
    if {$zoom} { ::ufv::zoom_to $mol [keysel [dict get $q keys]] 3 }
    ::ufv::emit focus $mol $pos
    return $f
}
proc ::ufv::feat::refocus {mol} {
    set f [::ufv::st::get $mol focus ""]
    if {$f eq ""} { if {[::ufv::layers::exists $mol focus]} { ::ufv::layers::delete $mol focus }; return }
    set dc [expr {[dict exists $f declutter] ? [dict get $f declutter] : 1}]
    if {[dict get $f kind] eq "ligand"} { ligand_focus $mol {*}[dict get $f ligand] 0 $dc } else { focus $mol [dict get $f pos] 0 $dc }
}

# ufv_focus <pos> ?molid? ?zoom?   - returns the selection text of the residue ("" if not mapped)
proc ufv_focus {pos {molid ""} {zoom 1}} {
    set mol [::ufv::resolve $molid]
    set f [::ufv::feat::focus $mol $pos $zoom]
    if {$f eq ""} { return "" }
    ::ufv::status "mol $mol: focused on UniProt $pos ([llength [dict get $f near]] residues within [dict get $f dist] A)."
    return [dict get $f sel]
}

# Ligand focus: one copy (chain + resid) or every copy of a component, + its Nearby pocket (licorice by
# annotation). The ligand layer is hidden meanwhile. Returns the focus dict.
proc ::ufv::feat::ligand_focus {mol resn {chain ""} {resid ""} {zoom 1} {declutter 1}} {
    set uid [::ufv::mol_uid $mol]
    set lig "resname [list $resn]"
    if {$chain ne "" && $resid ne ""} {
        append lig " and chain [::ufv::map::chain_token $chain] and resid [::ufv::map::compress [list $resid]]"
    }
    set n [::ufv::with_sel s $mol "$::ufv::lig_base and $lig" { $s num }]
    if {!$n} { ::ufv::status "no ligand $resn $chain $resid in mol $mol."; return "" }
    set d $::ufv::nearby
    set prev [::ufv::st::get $mol focus ""]
    set hidden [expr {$prev ne "" ? [dict get $prev hidden] : {}}]
    set keys {}; set near {}
    ::ufv::with_sel p $mol "same residue as ($::ufv::poly_base[notx $mol] and within $d of ($::ufv::lig_base and $lig))" {
        set seen [dict create]
        foreach {ch ri nm} [join [$p get {chain resid name}]] {
            if {![dict exists $seen $ch|$ri]} { dict set seen $ch|$ri 1; lappend keys $ch $ri }
            if {$nm eq "CA"} {
                set u [::ufv::map::r2u $mol $ch $ri]
                if {$u ne "" && $u ni $near} { lappend near $u }
            }
        }
    }
    ::ufv::layers::batch {
        if {[::ufv::layers::exists $mol pocket]} { ::ufv::layers::delete $mol pocket }
        set lighidden 0
        if {$declutter && [::ufv::layers::visible $mol lig]} {
            foreach i [::ufv::layers::reps $mol lig] { mol showrep $mol $i 0 }
            set lighidden 1
        }
        set groups [anno_groups $mol $keys $uid]
        lappend groups [list [::ufv::palette::id #4caf50] "($::ufv::lig_base and $lig) and element C" {Licorice 0.35 12 12}]
        lappend groups [list Element "($::ufv::lig_base and $lig) and not element C" {Licorice 0.35 12 12}]
        ::ufv::layers::define $mol focus $groups -style {Licorice 0.2 12 12} -kind other -rebuild [list ::ufv::feat::refocus $mol]
    }
    set f [dict create kind ligand ligand [list $resn $chain $resid] dist $d near [lsort -integer $near] keys $keys \
               hidden $hidden lighidden $lighidden sel "$::ufv::lig_base and $lig" partners {} declutter $declutter]
    ::ufv::st::put $mol focus $f
    if {$zoom} { ::ufv::zoom_to $mol "$::ufv::lig_base and $lig" 6 }
    ::ufv::emit focus $mol $resn
    return $f
}

# ufv_ligand_focus <resn> ?chain? ?resid? ?molid?  - returns the pocket's UniProt positions
proc ufv_ligand_focus {resn {chain ""} {resid ""} {molid ""}} {
    set mol [::ufv::resolve $molid]
    set f [::ufv::feat::ligand_focus $mol $resn $chain $resid]
    if {$f eq ""} { return "" }
    ::ufv::status "mol $mol: ligand $resn $chain $resid focused ([llength [dict get $f near]] pocket residues within [dict get $f dist] A)."
    return [dict get $f near]
}

proc ufv_resetview {{molid ""}} {
    set mol [::ufv::resolve $molid]
    ::ufv::feat::clear_focus $mol
    ::ufv::view::reset $mol
    ::ufv::status "view reset."
}

# ufv_show_pocket <positions | pocket=<id>> ?molid?  - lining residues as orange licorice + a translucent
# surface, zoomed. pocket=<id> takes a ProtVar pocket already fetched (residue report or ufv_pockets_all).
proc ufv_show_pocket {positions {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {[regexp {^pocket=(.+)$} $positions -> pid]} {
        set pk [::ufv::feat::find_pocket [::ufv::mol_uid $mol] $pid]
        if {$pk eq ""} { ::ufv::status "ProtVar pocket $pid not fetched yet (ufv_pockets_all or a residue report)."; return 0 }
        set positions [dict get $pk resid]
    }
    set positions [string map {+ " " , " "} $positions]
    set s [::ufv::map::sel $mol $positions]
    if {$s eq ""} { ::ufv::status "no pocket residues mapped."; return 0 }
    ::ufv::layers::batch {
        ::ufv::feat::clear_focus $mol
        ::ufv::layers::define $mol pocket [list \
            [list [::ufv::palette::id #fb8c00] "same residue as ($s) and not hydrogen" {Licorice 0.25 12 12}] \
            [list [::ufv::palette::id #ffd180] "same residue as ($s)" {QuickSurf 1.0 0.5 1.0 1.0}]] \
            -style {Licorice 0.25 12 12} -kind other -rebuild [list ufv_show_pocket $positions $mol]
        catch { mol modmaterial [lindex [::ufv::layers::reps $mol pocket] 1] $mol Transparent }
    }
    ::ufv::st::put $mol focus [dict create kind pocket positions $positions hidden {} near $positions partners {}]
    ::ufv::zoom_to $mol "same residue as ($s)" 2
    ::ufv::status "pocket lining shown ([llength $positions] residues)."
    return [llength $positions]
}

# Superpose every other molecule onto refmol. Molecules that both carry a UniProt numbering are fitted on the
# Ca atoms of their COMMON UniProt positions (sequence-aware, like PyMOL's align); otherwise a backbone fit when
# the atom counts match, else a centroid shift. Returns {molid {rmsd n} ...}.
proc ::ufv::feat::ca_by_uni {mol} {
    set out [dict create]
    ::ufv::with_sel s $mol "$::ufv::poly_base and name CA" {
        foreach {i ch ri ins} [join [$s get {index chain resid insertion}]] {
            if {[string trim $ins] ne ""} continue
            set u [::ufv::map::r2u $mol $ch $ri]
            if {$u ne "" && ![dict exists $out $u]} { dict set out $u $i }
        }
    }
    return $out
}
proc ufv_align {{refmol ""}} {
    set ref [::ufv::resolve $refmol]
    set out [dict create]
    set rmap [expr {[::ufv::st::has $ref] ? [::ufv::feat::ca_by_uni $ref] : {}}]
    foreach m [molinfo list] {
        if {$m == $ref} continue
        set done 0
        if {[dict size $rmap] && [::ufv::st::has $m]} {
            set mmap [::ufv::feat::ca_by_uni $m]
            set pairs {}
            dict for {u i} $mmap { if {[dict exists $rmap $u]} { lappend pairs [list $i [dict get $rmap $u]] } }
            if {[llength $pairs] >= 3} {
                set pairs [lsort -integer -index 0 $pairs]
                set mi {}; set ri {}
                foreach p $pairs { lappend mi [lindex $p 0]; lappend ri [lindex $p 1] }
                set rs [lsort -integer $ri]
                set ord {}
                foreach r $ri { lappend ord [lsearch -sorted -integer $rs $r] }
                set sm [atomselect $m "index $mi"]; set sr [atomselect $ref "index $rs"]; set all [atomselect $m all]
                if {![catch {
                    $all move [measure fit $sm $sr order $ord]
                    set mx [$sm get {x y z}]; set rx [$sr get {x y z}]
                    set sum 0.0
                    foreach a $mx o $ord { set sum [expr {$sum + pow([veclength [vecsub $a [lindex $rx $o]]], 2)}] }
                    set rmsd [expr {sqrt($sum / [llength $mx])}]
                } err]} {
                    dict set out $m [list $rmsd [llength $pairs]]
                    puts [format "\[UFV\] aligned mol %s onto mol %s: RMSD %.2f A over %d Ca (common UniProt positions)." \
                        $m $ref $rmsd [llength $pairs]]
                    set done 1
                }
                $sm delete; $sr delete; $all delete
            }
        }
        if {$done} continue
        set rs [atomselect $ref "protein and backbone"]
        set ms [atomselect $m "protein and backbone"]
        set all [atomselect $m all]
        catch {
            if {[$ms num] > 0 && [$ms num] == [$rs num]} {
                $all move [measure fit $ms $rs]
                dict set out $m [list [measure rmsd $ms $rs] [$ms num]]
            } elseif {[$ms num] > 0 && [$rs num] > 0} {
                $all moveby [vecsub [measure center $rs] [measure center $ms]]
                dict set out $m [list "" 0]
            }
        }
        $rs delete; $ms delete; $all delete
    }
    ::ufv::status "aligned [dict size $out] mol(s) to mol $ref."
    return $out
}

# ---------------------------------------------------------------------------------------------
# Residue report, CSV, links
# ---------------------------------------------------------------------------------------------
# Accession-only text report (no molecule needed): the backend's plain residue text.
proc ::ufv::residue_async {pos uid cb} {
    if {$uid eq ""} { set uid [::ufv::need] }
    ::ufv::backend::async residue [dict create uid $uid pos $pos] [list ::ufv::_residue_done $pos $cb]
}
proc ::ufv::_residue_done {pos cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    ::ufv::emit report $pos [dict get $res text]
    {*}$cb 1 [dict get $res text]
}

# Rich residue report of a molecule (structured for the panel; see API.md):
#   ::ufv::report_async mol pos cb ?services?
# Nearby (current frame) and ligands within 5 A are computed here in VMD; the backend adds proximity, flags of the
# analyses already run, and (services=1) ProtVar predictors / CADD / AutoPocket + PDBe-KB. With services the quick
# part is emitted first (event report_data), then the full one; cb gets the final dict (its 'text' = console report).
proc ::ufv::report_async {mol pos cb {services 1}} {
    set uid [::ufv::mol_uid $mol]
    set f [::ufv::st::get $mol focus ""]
    if {$f ne "" && [dict get $f kind] eq "residue" && [dict get $f pos] == $pos && [dict get $f dist] == $::ufv::nearby
        && [dict get $f mapver] == [::ufv::st::get $mol mapver 0] && [dict get $f frame] == [molinfo $mol get frame]} {
        set nb [dict create near [dict get $f near] partners [dict get $f partners] dist [dict get $f dist]]
    } else {
        set nb [::ufv::feat::nearby_query $mol $pos]
        if {$nb eq ""} { set nb [dict create near {} partners {} dist $::ufv::nearby] }
    }
    set ligs {}
    set s [::ufv::map::sel $mol [list $pos]]
    if {$s ne ""} {
        ::ufv::with_sel l $mol "($::ufv::lig_base) and within 5 of ($s)" { set ligs [lsort -unique [$l get resname]] }
        # PyMOL `around 5` excludes the residue's own atoms: ligands only (never the residue itself)
    }
    set pj {}
    foreach p [dict get $nb partners] { lappend pj "\[[::ufv::json::str [lindex $p 0]],[::ufv::json::str [lindex $p 1]]\]" }
    set nbj "{\"dist\":[dict get $nb dist],\"near\":[::ufv::feat::jints [dict get $nb near]],\"partners\":\[[join $pj ,]\]}"
    set a [::ufv::feat::ctx $mol 1]
    dict set a pos $pos
    dict set a nearby [::ufv::json::raw $nbj]
    dict set a ligands [::ufv::json::raw [::ufv::feat::jstrs $ligs]]
    dict set a ran [join [::ufv::st::get $mol ran {}] " "]
    set tok [incr ::ufv::feat::seq]
    ::ufv::st::put $mol reportseq $tok
    if {$services} {
        dict set a services 0
        ::ufv::backend::async report $a [list ::ufv::_report_part $mol $pos $tok]
    }
    dict set a services $services
    ::ufv::backend::async report $a [list ::ufv::_report_done $mol $pos $tok $cb]
}
proc ::ufv::_report_part {mol pos tok ok res} {
    if {!$ok || ![::ufv::st::has $mol] || [::ufv::st::get $mol reportseq] != $tok} return
    if {[::ufv::st::get $mol reportfull] == $tok} return
    ::ufv::emit report_data $mol $pos $res
}
proc ::ufv::_report_done {mol pos tok cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    # a newer report of this molecule was requested meanwhile (replies come back in any order): only this
    # request's caller gets it; the panel / hooks keep showing the newest report
    if {[::ufv::st::has $mol] && [::ufv::st::get $mol reportseq] != $tok} { {*}$cb 1 $res; return }
    if {[::ufv::st::has $mol]} { ::ufv::st::put $mol reportfull $tok }
    ::ufv::emit report_data $mol $pos $res
    ::ufv::emit report $pos [dict get $res text]
    {*}$cb 1 $res
}

# ufv_report <pos> ?molid? ?services?  - focus the residue, print the report, return the report dict
proc ufv_report {pos {molid ""} {services 1}} {
    if {![llength [molinfo list]]} {
        set txt [::ufv::await [list ::ufv::residue_async $pos [::ufv::need]]]
        puts $txt
        return [dict create pos $pos text $txt]
    }
    set mol [::ufv::resolve $molid]
    catch { ::ufv::feat::focus $mol $pos }
    set r [::ufv::await [list ::ufv::report_async $mol $pos] ]
    if {![llength [::ufv::hook list report]]} { puts [dict get $r text] }
    return $r
}
# ufv_residue <pos> ?molid?  - report text (returned + printed) and the residue focused
proc ufv_residue {pos {molid ""}} {
    return [dict get [ufv_report $pos $molid] text]
}

# Per-residue ligand contacts (CSV nearby_ligands): {pos {CCD ..}} for ligands within 5 A of a residue's atoms
# (current frame; transplant chains excluded on the polymer side). One `measure contacts`.
proc ::ufv::feat::lig_contacts {mol {cutoff 5.0}} {
    set out [dict create]
    set nl [::ufv::with_sel l $mol $::ufv::lig_base { $l num }]
    if {!$nl} { return $out }
    set L [atomselect $mol $::ufv::lig_base]
    set P [atomselect $mol "$::ufv::poly_base[notx $mol] and within $cutoff of ($::ufv::lig_base)"]
    set res [dict create]
    if {[$P num]} {
        lassign [measure contacts $cutoff $L $P] li pi
        if {[llength $li]} {
            set ln [dict create]; foreach i [$L list] n [$L get resname] { dict set ln $i $n }
            set pk [dict create]; foreach i [$P list] k [$P get {chain resid}] { dict set pk $i $k }
            foreach a $li b $pi {
                set k [dict get $pk $b]
                set n [dict get $ln $a]
                if {![dict exists $res $k] || $n ni [dict get $res $k]} { dict lappend res $k $n }
            }
        }
    }
    $L delete; $P delete
    dict for {k names} $res {
        lassign $k ch ri
        set u [::ufv::map::r2u $mol $ch $ri]
        if {$u eq ""} continue
        foreach n $names { if {![dict exists $out $u] || $n ni [dict get $out $u]} { dict lappend out $u $n } }
    }
    return $out
}

# ::ufv::csv_async mol out protvar cb  (mol "" = sequence-only matrix)   -> cb ok path
proc ::ufv::csv_async {mol out protvar cb} {
    if {$mol eq ""} {
        set uid [::ufv::need]
        ::ufv::backend::async csv [dict create uid $uid out [file normalize $out]] [list ::ufv::_csv_done $cb]
        return
    }
    set a [::ufv::feat::ctx $mol 1]
    dict set a out [file normalize $out]
    dict set a protvar [expr {$protvar ? 1 : 0}]
    dict set a chains [::ufv::json::raw [::ufv::feat::jstrs [::ufv::map::chains $mol]]]
    set lc [::ufv::feat::lig_contacts $mol]
    set p {}
    dict for {u names} $lc { lappend p "[::ufv::json::str $u]:[::ufv::feat::jstrs [lsort $names]]" }
    dict set a lig_contacts [::ufv::json::raw "{[join $p ,]}"]
    if {$protvar} { ::ufv::feat::poll protvar:[dict get $a uid] "ProtVar" }
    ::ufv::backend::async csv_full $a [list ::ufv::_csv_done $cb]
}
proc ::ufv::_csv_done {cb ok res} {
    if {!$ok} { ::ufv::status "CSV failed: $res"; {*}$cb 0 $res; return }
    ::ufv::status "CSV written: [dict get $res path]"
    {*}$cb 1 [dict get $res path]
}

# ufv_csv ?filename? ?molid? ?uid? ?protvar?  - the per-residue matrix (structure tiers when a molecule is loaded;
# protvar=1 adds ProtVar conservation / EVE / ESM1b / M3D - one request per residue, slow)
proc ufv_csv {{outfile ""} {molid ""} {uid ""} {protvar 0}} {
    if {$uid ne ""} { set uid [::ufv::need $uid] }
    set mol ""
    if {[llength [molinfo list]]} {
        set mol [::ufv::resolve $molid]
        if {$uid ne ""} { ::ufv::st::put $mol uid $uid }
        set uid [::ufv::mol_uid $mol]
    } else {
        set uid [::ufv::need $uid]
    }
    if {$outfile eq ""} { set outfile [file normalize "${uid}_residue_annotations.csv"] }
    return [::ufv::await [list ::ufv::csv_async $mol $outfile $protvar]]
}

# Poll a backend progress key while a long request runs: status line + event `progress key done total`.
proc ::ufv::feat::poll {key label} {
    set ::ufv::feat::polling($key) 1
    after 800 [list ::ufv::feat::_poll $key $label]
}
proc ::ufv::feat::_poll {key label} {
    if {![info exists ::ufv::feat::polling($key)]} return
    ::ufv::backend::async progress [dict create key $key] [list ::ufv::feat::_polled $key $label]
}
proc ::ufv::feat::_polled {key label ok res} {
    if {![info exists ::ufv::feat::polling($key)]} return
    if {$ok && [dict exists $res done]} {
        ::ufv::emit progress $key [dict get $res done] [dict get $res total]
        ::ufv::status "$label: [dict get $res done]/[dict get $res total] residues"
        if {[dict get $res done] >= [dict get $res total]} { unset -nocomplain ::ufv::feat::polling($key); return }
    }
    after 1000 [list ::ufv::feat::_poll $key $label]
}
proc ::ufv::feat::poll_stop {key} { unset -nocomplain ::ufv::feat::polling($key) }

# Open a URL in the system browser (Windows: rundll32 - safe with & and ? in the URL).
proc ::ufv::open_url {url} {
    switch -- $::tcl_platform(platform) {
        windows { exec {*}[auto_execok rundll32] url.dll,FileProtocolHandler $url & }
        default {
            if {$::tcl_platform(os) eq "Darwin"} { exec open $url & } else { exec xdg-open $url & }
        }
    }
}
proc ufv_open_uniprot {{uid ""}} {
    if {$uid eq ""} { set uid $::ufv::uid }
    if {$uid eq ""} { error "no accession" }
    ::ufv::open_url "https://www.uniprot.org/uniprotkb/$uid/entry"
}
