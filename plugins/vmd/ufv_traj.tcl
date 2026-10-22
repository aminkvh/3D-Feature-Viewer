# =============================================================================================
# ufv_traj.tcl - OWNER: the TRAJECTORY agent. Trajectory analyses for multi-frame molecules (MD trajectories,
# NMR ensembles): RMSD / Rg / RMSF, distance time series, pocket dynamics, ligand / partner / Nearby contact
# occupancy, annotation dynamics, per-frame structure tiers, CSV + B-factor export, and the embeddable Tk panel
# (::ufv::traj::panel). Parity with the PyMOL plugin's ufv_traj_* commands (same definitions, same numbers).
#
# Sourced by ufv_vmd.tcl after ufv_backend / ufv_palette / ufv_layers / ufv_commands / ufv_features.
# Built on the documented foundation API only (API.md). Backend ops used: traj, traj_pad, traj_ctx; structure tiers go
# through ufv_hotspots / ufv_pockets (op `analyses` with the molecule context ::ufv::feat::ctx, PAE-gated like PyMOL).
#
# Engines (::ufv::traj::cfg(engine) auto | native | backend):
#   native  - VMD `measure`: RMSD / Rg / RMSF (fit + move on a hidden scratch molecule holding the window's
#             Calpha frames, never the user's coordinates), Calpha / centroid distances, contact occupancy
#             (`measure contacts` per frame). Chunked: the analysis runs as a coroutine that yields to the
#             event loop every ~40 ms (VMD / Tk stay responsive, ::ufv::traj::cancel stops it).
#   backend - the window's atom sets are written with `animate write dcd ... sel` and ufv_core.py runs the
#             PyMOL plugin's numpy maths on them (bit-identical to the PyMOL plugin).
#   auto    - native where VMD is fast (rmsf, ca/centroid distances, contacts); backend for pocket dynamics,
#             closest-heavy-atom distances and annotation dynamics context.
# Display: RMSF / occupancy go into the `user` field of every frame + ONE continuous cartoon mode; distances
# are VMD Bonds labels (atom pairs, follow frames natively) or ONE graphics set redrawn by a debounced
# vmd_frame trace (centroid pairs). No per-frame objects.
# Frames are VMD frames (0-based). Window: start / stop (end) / stride / ref (first = first window frame).
# =============================================================================================

namespace eval ::ufv::traj {
    variable cfg
    if {![array exists cfg]} {
        array set cfg {start 0 stop end stride 1 ref first contact 4.5 dthr 2.0 engine auto show_pairs 3}
    }
    variable R              ;# mol -> dict kind -> result (dict order = recency)
    variable RF             ;# "mol,kind" -> dict "chain resv" -> touched window-frame indices (contacts)
    variable jobs           ;# coroutine -> dict (mol kind cancel)
    variable seq 0
    variable tslice 0
    variable nfile 0
    variable ctx            ;# mol -> {sig result} annotation context
    variable tiercache      ;# "mol,uid,mapver,frame,q" -> tiers result (per-frame tier cache)
    variable gfx            ;# mol -> graphics ids of the centroid-distance set
    variable cen            ;# mol -> list of {idxA idxB} centroid pairs drawn
    variable labels         ;# mol -> atom pairs of the Bonds labels this module added
    variable trace_after ""
    variable pw ""          ;# the panel widget
    variable P              ;# panel state
    if {![array exists R]} { array set R {} }
    if {![array exists RF]} { array set RF {} }
    if {![array exists jobs]} { array set jobs {} }
    if {![array exists ctx]} { array set ctx {} }
    if {![array exists tiercache]} { array set tiercache {} }
    if {![array exists gfx]} { array set gfx {} }
    if {![array exists cen]} { array set cen {} }
    if {![array exists labels]} { array set labels {} }
    if {![array exists P]} { array set P {} }
    variable POLY "(protein or nucleic)"
    variable LIG "(not (protein or nucleic) and not water)"
    variable KINDS {rmsf dist pockets ligand partner annot nearby}
    variable LABELS {
        {RMSD / RMSF} rmsf {}
        {Distances} dist {auto | 120~245; A/45~lig:ATP; 10+12~200}
        {Pocket dynamics} pockets {auto | cp; pv; lig; 120+124+156}
        {Ligand contacts} ligand {}
        {Partner contacts} partner {}
        {Nearby occupancy} nearby {UniProt position}
        {Annotation dynamics} annot {}
    }
}

# ---------------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_mean {l} {
    if {![llength $l]} { return 0.0 }
    set s 0.0; foreach v $l { set s [expr {$s + $v}] }
    return [expr {$s / [llength $l]}]
}
proc ::ufv::traj::_std {l} {
    if {![llength $l]} { return 0.0 }
    set m [_mean $l]; set s 0.0
    foreach v $l { set s [expr {$s + ($v - $m) * ($v - $m)}] }
    return [expr {sqrt($s / [llength $l])}]
}
proc ::ufv::traj::_min {l} { return [tcl::mathfunc::min {*}$l] }
proc ::ufv::traj::_max {l} { return [tcl::mathfunc::max {*}$l] }
# numpy.percentile (linear interpolation)
proc ::ufv::traj::_percentile {vals q} {
    set s [lsort -real $vals]; set n [llength $s]
    if {!$n} { return 0.0 }
    set pos [expr {$q / 100.0 * ($n - 1)}]
    set i [expr {int(floor($pos))}]; set fr [expr {$pos - $i}]
    set a [lindex $s $i]; set b [lindex $s [expr {min($i + 1, $n - 1)}]]
    return [expr {$a + ($b - $a) * $fr}]
}
# on/off transitions of a sorted list of touched window-frame indices over F frames
proc ::ufv::traj::_switches {idx F} {
    if {$F < 2 || ![llength $idx]} { return 0 }
    set n [expr {[lindex $idx 0] > 0 ? 1 : 0}]
    set prev [lindex $idx 0]
    foreach t [lrange $idx 1 end] { if {$t - $prev > 1} { incr n 2 }; set prev $t }
    if {$prev < $F - 1} { incr n }
    return $n
}
proc ::ufv::traj::_jstr {s} { return [::ufv::json::str $s] }
proc ::ufv::traj::_jints {l} { return "\[[join $l ,]\]" }
proc ::ufv::traj::_jstrs {l} { set o {}; foreach s $l { lappend o [_jstr $s] }; return "\[[join $o ,]\]" }
proc ::ufv::traj::_jnull {v} { return [expr {$v eq "" ? "null" : $v}] }
proc ::ufv::traj::_tmp {ext} {
    variable nfile
    return [file join [::ufv::backend::tmpdir] "ufv_traj_[pid]_[incr nfile].$ext"]
}
proc ::ufv::traj::_pct {o} { return [format "%.0f%%" [expr {100.0 * $o}]] }

# ---------------------------------------------------------------------------------------------
# Window, Calpha table, atom sets
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::nframes {mol} { return [molinfo $mol get numframes] }

# Analysis window: list of VMD frames (start..stop every stride; stop "end" = last frame).
proc ::ufv::traj::frames {mol {start ""} {stop ""} {stride ""}} {
    variable cfg
    set n [nframes $mol]
    if {$n < 1} { error "mol $mol has no coordinates" }
    foreach v {start stop stride} { if {[set $v] eq ""} { set $v $cfg($v) } }
    set s [expr {$start in {"" first} ? 0 : int($start)}]
    set e [expr {$stop in {"" end last} ? $n - 1 : int($stop)}]
    set k [expr {max(1, int($stride))}]
    set s [expr {min(max(0, $s), $n - 1)}]
    set e [expr {min(max($s, $e), $n - 1)}]
    set out {}
    for {set f $s} {$f <= $e} {incr f $k} { lappend out $f }
    return $out
}

# Reference frame ("first" = first frame of the window).
proc ::ufv::traj::ref {mol frames {ref ""}} {
    variable cfg
    if {$ref eq ""} { set ref $cfg(ref) }
    if {$ref in {"" first}} { return [lindex $frames 0] }
    return [expr {min(max(0, int($ref)), [nframes $mol] - 1)}]
}

# Deduplicated polymer Calpha (first per chain + resid, like the PyMOL plugin):
#   dict res {{chain resi resv resn uni} ..}  idx {atom index ..}  cols {column in the all-CA selection ..}
proc ::ufv::traj::ca {mol} {
    return [::ufv::st::cached $mol traj_ca {
        set res {}; set idx {}; set cols {}; set seen [dict create]; set i 0
        ::ufv::with_sel s $mol "$::ufv::traj::POLY and name CA" {
            foreach {ix ch ri rn ins} [join [$s get {index chain resid resname insertion}]] {
                if {![dict exists $seen "$ch $ri"]} {
                    dict set seen "$ch $ri" 1
                    set ins [string trim $ins]
                    lappend res [list $ch "$ri$ins" $ri $rn [::ufv::map::r2u $mol $ch $ri]]
                    lappend idx $ix; lappend cols $i
                }
                incr i
            }
        }
        dict create res $res idx $idx cols $cols
    }]
}

proc ::ufv::traj::_res_json {res} {
    set o {}
    foreach r $res {
        lassign $r ch ri rv rn u
        lappend o "\[[_jstr $ch],[_jstr $ri],$rv,[_jstr $rn],[_jnull $u]\]"
    }
    return "\[[join $o ,]\]"
}

# Write the window frames (and the reference frame when it lies outside the window) of a selection as DCD.
# -> {window.dcd ref.dcd|""}
proc ::ufv::traj::_write_dcd {mol seltext frames ref} {
    set d [_tmp dcd]
    set s [lindex $frames 0]; set e [lindex $frames end]
    set k [expr {[llength $frames] > 1 ? [lindex $frames 1] - $s : 1}]
    set rp ""
    ::ufv::with_sel sel $mol $seltext {
        if {![$sel num]} { error "empty atom set: $seltext" }
        animate write dcd $d beg $s end $e skip $k sel $sel waitfor all $mol
        if {$ref ne "" && $ref ni $frames} {
            set rp [_tmp dcd]
            animate write dcd $rp beg $ref end $ref skip 1 sel $sel waitfor all $mol
        }
    }
    return [list $d $rp]
}

proc ::ufv::traj::_idxsel {idx} { return "index [::ufv::map::compress $idx]" }

# ---------------------------------------------------------------------------------------------
# Jobs: every analysis runs as a coroutine; _slice yields to the event loop every ~40 ms.
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_slice {{force 0}} {
    variable tslice
    variable jobs
    set co [info coroutine]
    if {$co eq ""} return
    if {!$force && [clock milliseconds] - $tslice < 40} return
    after 0 [list $co ""]
    yield
    set tslice [clock milliseconds]
    if {[info exists jobs($co)] && [dict get $jobs($co) cancel]} { error "cancelled" }
}

# Backend request from inside a job (coroutine: waits without blocking; plain call otherwise).
proc ::ufv::traj::_backend {op argdict} {
    set co [info coroutine]
    if {$co eq ""} { return [::ufv::backend::call $op $argdict 600000] }
    ::ufv::backend::async $op $argdict [list ::ufv::traj::_resume $co]
    lassign [yield] ok res
    if {!$ok} { error "backend $op: $res" }
    return $res
}
proc ::ufv::traj::_resume {co ok res} {
    if {[llength [info commands $co]]} { $co [list $ok $res] }
}

# Run one analysis asynchronously: {*}$cb ok result. kind = rmsf dist pockets ligand partner nearby annot.
#   opts: engine color show cutoff threshold frames ref start stop stride
proc ::ufv::traj::run_async {mol kind spec opts cb} {
    variable seq
    variable jobs
    variable KINDS
    if {$kind ni $KINDS} { error "unknown trajectory analysis '$kind' ($KINDS)" }
    set co ::ufv::traj::job[incr seq]
    set jobs($co) [dict create mol $mol kind $kind cancel 0]
    coroutine $co ::ufv::traj::_job $co $mol $kind $spec $opts $cb
    return $co
}
proc ::ufv::traj::_job {co mol kind spec opts cb} {
    variable jobs
    variable tslice
    after 0 [list $co ""]
    yield
    set tslice [clock milliseconds]
    set t0 [clock milliseconds]
    set code [catch { analyse $mol $kind $spec $opts } res]
    unset -nocomplain jobs($co)
    if {$code} {
        if {$res ne "cancelled"} { ::ufv::status "trajectory $kind: $res" }
        after 0 [list {*}$cb 0 $res]
    } else {
        dict set res wall_ms [expr {[clock milliseconds] - $t0}]
        ::ufv::emit traj $mol $kind
        after 0 [list {*}$cb 1 $res]
    }
}
proc ::ufv::traj::cancel {{co ""}} {
    variable jobs
    foreach c [expr {$co eq "" ? [array names jobs] : [list $co]}] {
        if {[info exists jobs($c)]} { dict set jobs($c) cancel 1 }
    }
}
proc ::ufv::traj::running {} { variable jobs; return [array size jobs] }

# Synchronous form (console commands): the event loop keeps running while it waits.
proc ::ufv::traj::run {mol kind {spec ""} {opts {}}} {
    return [::ufv::await [list ::ufv::traj::run_async $mol $kind $spec $opts] 900000]
}

# ---------------------------------------------------------------------------------------------
# Annotation context (tags / groups): hotspots + constraint-pocket classes of the CURRENT frame (if computed)
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_tiers_now {mol} {
    variable tiercache
    set uid [::ufv::mol_uid $mol]
    set ver [::ufv::st::get $mol mapver 0]
    set fr [molinfo $mol get frame]
    set c [::ufv::st::get $mol {cache tiers} ""]
    if {$c ne ""} {
        lassign [lindex $c 0] cu cv cf cq
        if {$cu eq $uid && $cv == $ver && $cf == $fr} { return [lindex $c 1] }
    }
    foreach k [array names tiercache "$mol,$uid,$ver,$fr,*"] { return $tiercache($k) }
    return {}
}

proc ::ufv::traj::_lig_near {mol frame cutoff ligsel {noh 1}} {
    # UniProt positions of polymer residues within cutoff of a ligand selection (one frame)
    set h [expr {$noh ? "and noh" : ""}]
    set out [dict create]
    ::ufv::with_sel s $mol "name CA and same residue as ($::ufv::traj::POLY $h and within $cutoff of ($ligsel))" {
        $s frame $frame; $s update
        foreach {ch ri} [join [$s get {chain resid}]] {
            set u [::ufv::map::r2u $mol $ch $ri]
            if {$u ne ""} { dict set out $u 1 }
        }
    }
    return [lsort -integer [dict keys $out]]
}

proc ::ufv::traj::_has_lig {mol} {
    return [::ufv::with_sel s $mol $::ufv::traj::LIG { expr {[$s num] > 0} }]
}

# {groups {label {pos..}} tags {pos tag}} (backend traj_ctx; cached per mol + inputs)
proc ::ufv::traj::context {mol {with_lig 1}} {
    variable ctx
    set uid [::ufv::mol_uid $mol]
    set t [_tiers_now $mol]
    set a [dict create uid $uid]
    set hot {}; set cp {}
    if {$t ne ""} {
        set hot [lsort -integer [dict keys [dict get $t hotspots]]]
        dict for {p c} [dict get $t pockets] { lappend cp "\[$p,[_jstr $c]\]" }
    }
    dict set a hot [::ufv::json::raw [_jints $hot]]
    dict set a cp [::ufv::json::raw "\[[join $cp ,]\]"]
    set lp {}
    if {$with_lig && [_has_lig $mol]} {
        set lp [_lig_near $mol [molinfo $mol get frame] 5 "$::ufv::traj::LIG and noh"]
    }
    dict set a ligpocket [::ufv::json::raw [_jints $lp]]
    set sig [list $uid [::ufv::st::get $mol mapver 0] $hot $cp $lp]
    if {[info exists ctx($mol)] && [lindex $ctx($mol) 0] eq $sig} { return [lindex $ctx($mol) 1] }
    if {$uid eq ""} {
        set r [dict create groups {} tags {}]
    } else {
        set r [_backend traj_ctx $a]
    }
    dict set r args $a
    set ctx($mol) [list $sig $r]
    return $r
}

# ---------------------------------------------------------------------------------------------
# The analysis driver (runs inside a job)
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::analyse {mol kind spec opts} {
    variable cfg
    set n [nframes $mol]
    if {$n < 2 && $kind ne "annot"} { ::ufv::status "mol $mol has a single frame: analysing it as a 1-frame window." }
    set fr [expr {[dict exists $opts frames] ? [dict get $opts frames] : [frames $mol]}]
    set rf [ref $mol $fr [expr {[dict exists $opts ref] ? [dict get $opts ref] : ""}]]
    set engine [expr {[dict exists $opts engine] ? [dict get $opts engine] : $cfg(engine)}]
    if {$engine eq "auto"} {
        # per benchmark (500 frames): native rmsf / distances / Nearby are 2-3x faster than the round trip;
        # contact occupancy on large complexes is faster in numpy (exact pruning) than `measure contacts`
        set engine [dict get {rmsf native dist native nearby native pockets backend ligand backend partner backend annot backend} $kind]
    }
    set cut [expr {[dict exists $opts cutoff] && [dict get $opts cutoff] ne "" ? double([dict get $opts cutoff]) : double($cfg(contact))}]
    set thr [expr {[dict exists $opts threshold] && [dict get $opts threshold] ne "" ? double([dict get $opts threshold]) : double($cfg(dthr))}]
    set t0 [clock milliseconds]
    set ca [ca $mol]
    if {![llength [dict get $ca res]]} { error "mol $mol has no polymer CA atoms" }
    set P [dict create mol $mol kind $kind frames $fr ref $rf engine $engine cutoff $cut threshold $thr \
        spec $spec ca $ca uid [::ufv::mol_uid $mol] t0 $t0]
    switch -- $kind {
        rmsf    { set r [_rmsf $P] }
        dist    { set r [_dist $P] }
        pockets { set r [_pockets $P] }
        ligand - partner - nearby { set r [_contacts $P] }
        annot   { set r [_annot $P] }
    }
    dict set r kind $kind
    dict set r mol $mol
    dict set r frames $fr
    dict set r ref $rf
    dict set r res [dict get $ca res]
    if {![dict exists $r engine]} { dict set r engine $engine }
    _store $mol $kind $r
    set color [expr {![dict exists $opts color] || [dict get $opts color]}]
    set show [expr {![dict exists $opts show] || [dict get $opts show]}]
    apply $r $color $show
    return $r
}

proc ::ufv::traj::_store {mol kind r} {
    variable R
    if {![info exists R($mol)]} { set R($mol) [dict create] }
    dict unset R($mol) $kind
    dict set R($mol) $kind $r
}
proc ::ufv::traj::result {mol kind} {
    variable R
    if {[info exists R($mol)] && [dict exists $R($mol) $kind]} { return [dict get $R($mol) $kind] }
    return ""
}
proc ::ufv::traj::results {mol} {
    variable R
    return [expr {[info exists R($mol)] ? $R($mol) : {}}]
}

# Arguments common to every backend traj request: window, reference, Calpha DCD + table, context.
proc ::ufv::traj::_base_args {P files} {
    upvar 1 $files tmp
    set mol [dict get $P mol]
    set ca [dict get $P ca]
    lassign [_write_dcd $mol [_idxsel [dict get $ca idx]] [dict get $P frames] [dict get $P ref]] d rp
    lappend tmp $d $rp
    set a [dict create kind [dict get $P kind] key "mol$mol" uid [dict get $P uid] \
        frames [::ufv::json::raw [_jints [dict get $P frames]]] ref [dict get $P ref] \
        cutoff [dict get $P cutoff] threshold [dict get $P threshold] ca_dcd $d \
        res [::ufv::json::raw [_res_json [dict get $ca res]]]]
    if {$rp ne ""} { dict set a ca_ref $rp }
    if {[dict get $P spec] ne ""} { dict set a spec [::ufv::json::raw [_jstr [dict get $P spec]]] }
    return $a
}
proc ::ufv::traj::_cleanup {files} { foreach f $files { if {$f ne ""} { catch { file delete $f } } } }

proc ::ufv::traj::_backend_run {a files} {
    set code [catch { _backend traj $a } r]
    _cleanup $files
    if {$code} { error $r }
    dict set r engine backend
    return $r
}

# ---------------------------------------------------------------------------------------------
# RMSD / Rg / RMSF
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_rmsf {P} {
    set mol [dict get $P mol]
    set c [context $mol]
    if {[dict get $P engine] eq "backend"} {
        set files {}
        set a [_base_args $P files]
        foreach k {hot cp ligpocket} { dict set a $k [dict get $c args $k] }
        return [_backend_run $a $files]
    }
    set fr [dict get $P frames]; set ref [dict get $P ref]
    lassign [_scratch_fit $mol [dict get $P ca] $fr $ref] rmsd rg rmsf
    set res [dict get $P ca res]
    # annotation-aware table (same groups as the PyMOL plugin)
    set rows {}; set meta {}
    foreach g [lsort [dict keys [dict get $c groups]]] {
        set ps [dict create]; foreach p [dict get $c groups $g] { dict set ps $p 1 }
        set inn {}; set out {}
        foreach r $res v $rmsf {
            set u [lindex $r 4]
            if {$u eq ""} continue
            if {[dict exists $ps $u]} { lappend inn $v } else { lappend out $v }
        }
        if {[llength $inn] && [llength $out]} {
            set x [_mean $inn]; set y [_mean $out]
            lappend rows [list $g [llength $inn] [format %.2f $x] [format %.2f $y] [format %.2f [expr {$y != 0 ? $x / $y : 0.0}]]]
            lappend meta [dict create group $g]
        }
    }
    set top {}
    foreach i [lrange [lsort -real -decreasing -indices $rmsf] 0 14] {
        set r [lindex $res $i]
        lappend top [expr {[lindex $r 4] ne "" ? [lindex $r 4] : "[lindex $r 0]/[lindex $r 1]"}]
    }
    set sum [format "RMSD %.2f \u00b1 %.2f \u00c5 (max %.2f); mean RMSF %.2f \u00c5; most mobile: %s" \
        [_mean $rmsd] [_std $rmsd] [_max $rmsd] [_mean $rmsf] [join $top ,]]
    return [dict create series [list [list "RMSD (\u00c5)" $rmsd] [list "Rg (\u00c5)" $rg]] rmsd $rmsd rmsf $rmsf \
        header {Group n RMSF rest ratio} rows $rows rowmeta $meta summary $sum engine native \
        calc_ms [expr {[clock milliseconds] - [dict get $P t0]}]]
}

# Fit on a hidden scratch molecule: frame 0 = reference, frames 1..F = the window's Calpha.
# RMSD to the reference after fitting, Rg, then RMSF after two refits to the average (PyMOL plugin's _calc_rmsf).
proc ::ufv::traj::_scratch_fit {mol ca frames ref} {
    set sel [_idxsel [dict get $ca idx]]
    set pdb [_tmp pdb]
    lassign [_write_dcd $mol $sel $frames ""] dcd _x
    # the DCD write, the PDB write and the scratch load each get their own event-loop slice
    if {[catch { _slice } err]} { _cleanup [list $dcd]; error $err }
    ::ufv::with_sel s $mol $sel { $s frame $ref; $s writepdb $pdb }
    if {[catch { _slice } err]} { _cleanup [list $pdb $dcd]; error $err }
    set top [molinfo top]
    # `mol new` re-fits every molecule's view (even for a hidden scratch molecule): keep it exactly
    set sm [::ufv::view::keep { mol new $pdb type pdb autobonds off waitfor all }]
    set code [catch {
        mol delrep 0 $sm
        mol off $sm
        mol addfile $dcd type dcd waitfor all molid $sm
        if {$top >= 0 && $top in [molinfo list]} { mol top $top }
        set F [llength $frames]
        set all [atomselect $sm all]
        set r0 [atomselect $sm all frame 0]
        set rmsd {}; set rg {}
        for {set i 1} {$i <= $F} {incr i} {
            $all frame $i
            lappend rg [measure rgyr $all]
            $all move [measure fit $all $r0]
            lappend rmsd [measure rmsd $all $r0]
            if {$i % 50 == 0} { _slice }
        }
        for {set it 0} {$it < 2} {incr it} {
            $r0 set {x y z} [measure avpos $all first 1 last $F]
            for {set i 1} {$i <= $F} {incr i} {
                $all frame $i
                $all move [measure fit $all $r0]
                if {$i % 50 == 0} { _slice }
            }
        }
        set rmsf [measure rmsf $all first 1 last $F]
        $all delete; $r0 delete
    } err]
    catch { mol delete $sm }
    if {$top >= 0 && $top in [molinfo list]} { catch { mol top $top } }
    _cleanup [list $pdb $dcd]
    if {$code} { error $err }
    return [list $rmsd $rg $rmsf]
}

# ---------------------------------------------------------------------------------------------
# Distances
# ---------------------------------------------------------------------------------------------
# 'A/120' / '120' (UniProt position; first mapped chain that has it) -> "chain X and resid R" or ""
proc ::ufv::traj::_res_tok {mol tok} {
    if {![regexp {^\s*(?:([A-Za-z0-9]+)[/:])?(-?\d+)\s*$} $tok -> ch p]} { return "" }
    foreach c [expr {$ch ne "" ? [list $ch] : [::ufv::map::chains $mol]}] {
        set r [::ufv::map::u2r $mol $c $p]
        if {$r ne ""} { return "chain [::ufv::map::chain_token $c] and resid [::ufv::map::compress [list $r]]" }
    }
    return ""
}
proc ::ufv::traj::_lig_tok {tok} {
    set parts [split [string range $tok 4 end] /]
    set s "$::ufv::traj::LIG and noh and resname [lindex $parts 0]"
    if {[llength $parts] > 1 && [lindex $parts 1] ne ""} { append s " and chain [::ufv::map::chain_token [lindex $parts 1]]" }
    if {[llength $parts] > 2 && [lindex $parts 2] ne ""} { append s " and resid [lindex $parts 2]" }
    return $s
}

# Ligand residues {resn chain resid nheavy} with >= min heavy atoms, largest first (stable).
proc ::ufv::traj::ligand_ids {mol {minheavy 1}} {
    set c [dict create]
    ::ufv::with_sel s $mol "$::ufv::traj::LIG and noh" {
        foreach {rn ch ri} [join [$s get {resname chain resid}]] { dict incr c [list $rn $ch $ri] }
    }
    set out {}
    dict for {k n} $c { if {$n >= $minheavy} { lappend out [list {*}$k $n] } }
    return [lsort -integer -decreasing -index 3 $out]
}

# Suggested pairs (PyMOL _traj_auto_pairs): PTM <-> nearest pathogenic variant (Calpha <= 15 A, closest 6),
# the 2 largest ligands (>= 5 heavy atoms) <-> up to 4 lining residues (annotated first). Current frame.
proc ::ufv::traj::auto_pairs {mol {maxn 10}} {
    set c [context $mol 0]
    set fr [molinfo $mol get frame]
    set ca [ca $mol]
    set xyz [dict create]
    ::ufv::with_sel s $mol [_idxsel [dict get $ca idx]] {
        $s frame $fr
        foreach r [dict get $ca res] p [$s get {x y z}] {
            set u [lindex $r 4]
            if {$u ne "" && ![dict exists $xyz $u]} { dict set xyz $u $p }
        }
    }
    set g [dict get $c groups]
    set ptms {}; set path {}
    if {[dict exists $g "PTM sites"]} { foreach p [dict get $g "PTM sites"] { if {[dict exists $xyz $p]} { lappend ptms $p } } }
    if {[dict exists $g "Pathogenic variants"]} { foreach p [dict get $g "Pathogenic variants"] { if {[dict exists $xyz $p]} { lappend path $p } } }
    set ptms [lsort -integer -unique $ptms]; set path [lsort -integer -unique $path]
    set out {}
    if {[llength $ptms] && [llength $path]} {
        set cand {}
        foreach t $ptms {
            set best ""; set bp ""
            set xt [dict get $xyz $t]
            foreach q $path {
                if {$q == $t} continue
                set d [veclength [vecsub [dict get $xyz $q] $xt]]
                if {$best eq "" || $d < $best} { set best $d; set bp $q }
            }
            if {$best ne "" && $best <= 15.0} { lappend cand [list $best $t $bp] }
        }
        set seen [dict create]
        foreach x [lsort -real -index 0 [lsort -integer -index 2 [lsort -integer -index 1 $cand]]] {
            lassign $x d a b
            if {[dict exists $seen "$a $b"]} continue
            dict set seen "$a $b" 1
            lappend out [list "$a~$b" "PTM $a \u2194 pathogenic $b"]
            if {[llength $out] >= 6} break
        }
    }
    set tags [dict get $c tags]
    foreach l [lrange [ligand_ids $mol 5] 0 1] {
        lassign $l rn ch ri
        set lt "lig:$rn/$ch/$ri"
        set near [_lig_near $mol $fr 4.5 [_lig_tok $lt]]
        set near [lsort -command [list ::ufv::traj::_tagfirst $tags] $near]
        foreach q [lrange $near 0 3] { lappend out [list "$q~$lt" "$rn \u2194 lining $q"] }
    }
    return [lrange $out 0 [expr {$maxn - 1}]]
}
proc ::ufv::traj::_tagfirst {tags a b} {
    set ta [expr {![dict exists $tags $a]}]; set tb [expr {![dict exists $tags $b]}]
    if {$ta != $tb} { return [expr {$ta - $tb}] }
    return [expr {$a - $b}]
}

# Parse a pair spec into {label tok mode selA selB idxA idxB} (atom indices; mode ca | cen | min).
proc ::ufv::traj::parse_pairs {mol spec} {
    set toks {}
    foreach part [split [string map {; \n , \n} $spec] \n] {
        set part [string trim $part]
        if {$part eq ""} continue
        if {[string tolower $part] eq "auto"} { lappend toks {*}[auto_pairs $mol] } else { lappend toks [list $part $part] }
    }
    set pairs {}
    foreach t $toks {
        lassign $t tok label
        if {[string first ~ $tok] >= 0} {
            set i [string first ~ $tok]
            set a [string range $tok 0 [expr {$i - 1}]]; set b [string range $tok [expr {$i + 1}] end]
        } elseif {![regexp {^(.+?)\s*-\s*(.+)$} $tok -> a b]} {
            ::ufv::status "distance pair '$tok': use A~B (e.g. 120~245, A/120~lig:ATP, 10+12+15~200)."
            continue
        }
        set eps {}; set bad 0
        foreach e [list [string trim $a] [string trim $b]] {
            set el [string tolower $e]
            if {[string match lig:* $el]} {
                lappend eps [list lig [_lig_tok $e]]
            } elseif {[string match sel:* $el]} {
                lappend eps [list cen [string range $e 4 end]]
            } elseif {[string first + $e] >= 0 || [string index $e 0] eq "\{"} {
                set ss {}
                foreach x [regexp -all -inline {(?:[A-Za-z0-9]+/)?\d+} $e] {
                    set s [_res_tok $mol $x]
                    if {$s ne ""} { lappend ss "($s and name CA)" }
                }
                lappend eps [list cen [expr {[llength $ss] ? [join $ss { or }] : ""}]]
            } else {
                set s [_res_tok $mol $e]
                lappend eps [list res [expr {$s ne "" ? "($s) and $::ufv::traj::POLY and noh" : ""}]]
            }
            if {[lindex $eps end 1] eq ""} { set bad 1 }
        }
        if {$bad} { ::ufv::status "distance pair '$tok': endpoint not found in mol $mol."; continue }
        set atoms {}
        foreach ep $eps {
            lassign $ep k s
            set ix {}; set nm {}
            ::ufv::with_sel x $mol $s { set ix [$x get index]; set nm [$x get name] }
            set caidx {}
            if {$k eq "res"} { foreach i $ix n $nm { if {$n eq "CA"} { lappend caidx $i } } }
            lappend atoms [list $k $ix $caidx]
        }
        lassign $atoms ea eb
        lassign $ea ka ia caa; lassign $eb kb ib cab
        if {![llength $ia] || ![llength $ib]} { ::ufv::status "distance pair '$tok': no atoms."; continue }
        if {$ka eq "res" && $kb eq "res" && [llength $caa] && [llength $cab]} {
            set mode ca; set A [lrange $caa 0 0]; set B [lrange $cab 0 0]
        } elseif {"cen" in [list $ka $kb]} {
            set mode cen
            set A [expr {$ka eq "res" && [llength $caa] ? [lrange $caa 0 0] : $ia}]
            set B [expr {$kb eq "res" && [llength $cab] ? [lrange $cab 0 0] : $ib}]
        } else {
            set mode min; set A $ia; set B $ib
        }
        lappend pairs [dict create label $label tok $tok mode $mode idxA $A idxB $B]
    }
    if {![llength $pairs]} { error "no distance pairs resolved (spec: $spec)" }
    return $pairs
}

# mean position of atoms (sequential sum / n, like numpy's mean over the atom axis)
proc ::ufv::traj::_centroid {X col idx} {
    if {[llength $idx] == 1} { return [lindex $X [dict get $col [lindex $idx 0]]] }
    set sx 0.0; set sy 0.0; set sz 0.0
    foreach i $idx {
        lassign [lindex $X [dict get $col $i]] x y z
        set sx [expr {$sx + $x}]; set sy [expr {$sy + $y}]; set sz [expr {$sz + $z}]
    }
    set n [llength $idx]
    return [list [expr {$sx / $n}] [expr {$sy / $n}] [expr {$sz / $n}]]
}

proc ::ufv::traj::_dist {P} {
    set mol [dict get $P mol]
    set spec [dict get $P spec]; if {$spec eq ""} { set spec auto }
    set pairs [parse_pairs $mol $spec]
    set fr [dict get $P frames]
    set engine [dict get $P engine]
    set thr [dict get $P threshold]
    # which pairs go to the backend: closest-heavy-atom pairs always (numpy); all with engine backend
    set bk {}; set nat {}
    set i 0
    foreach pr $pairs {
        if {$engine eq "backend" || [dict get $pr mode] eq "min"} { lappend bk $i } else { lappend nat $i }
        incr i
    }
    set out [dict create]
    if {[llength $bk]} {
        set union [dict create]
        foreach i $bk { foreach x [concat [dict get [lindex $pairs $i] idxA] [dict get [lindex $pairs $i] idxB]] { dict set union $x 1 } }
        set uidx [lsort -integer [dict keys $union]]
        set col [dict create]; set j 0; foreach x $uidx { dict set col $x $j; incr j }
        set files {}
        set a [_base_args $P files]
        lassign [_write_dcd $mol [_idxsel $uidx] $fr [lindex $fr 0]] d _x
        lappend files $d
        dict set a dcd $d
        set pj {}
        foreach i $bk {
            set pr [lindex $pairs $i]
            set ca {}; foreach x [dict get $pr idxA] { lappend ca [dict get $col $x] }
            set cb {}; foreach x [dict get $pr idxB] { lappend cb [dict get $col $x] }
            lappend pj "{\"label\":[_jstr [dict get $pr label]],\"tok\":[_jstr [dict get $pr tok]],\"mode\":\"[dict get $pr mode]\",\"A\":[_jints $ca],\"B\":[_jints $cb],\"idxA\":[_jints [dict get $pr idxA]],\"idxB\":[_jints [dict get $pr idxB]]}"
        }
        dict set a pairs [::ufv::json::raw "\[[join $pj ,]\]"]
        set r [_backend_run $a $files]
        foreach i $bk s [dict get $r series] row [dict get $r rows] p [dict get $r pairs] {
            dict set out $i [list $s $row [dict get $p disp]]
        }
    }
    if {[llength $nat]} {
        set union [dict create]
        foreach i $nat { foreach x [concat [dict get [lindex $pairs $i] idxA] [dict get [lindex $pairs $i] idxB]] { dict set union $x 1 } }
        set uidx [lsort -integer [dict keys $union]]
        set col [dict create]; set j 0; foreach x $uidx { dict set col $x $j; incr j }
        set series [dict create]
        ::ufv::with_sel s $mol [_idxsel $uidx] {
            set k 0
            foreach f $fr {
                $s frame $f
                set X [$s get {x y z}]
                foreach i $nat {
                    set pr [lindex $pairs $i]
                    lassign [_centroid $X $col [dict get $pr idxA]] ax ay az
                    lassign [_centroid $X $col [dict get $pr idxB]] bx by bz
                    dict lappend series $i [expr {sqrt(($ax-$bx)*($ax-$bx) + ($ay-$by)*($ay-$by) + ($az-$bz)*($az-$bz))}]
                }
                if {[incr k] % 50 == 0} { _slice }
            }
        }
        foreach i $nat {
            set pr [lindex $pairs $i]
            set d [dict get $series $i]
            set fl [expr {[_max $d] - [_min $d] > $thr}]
            set row [list [dict get $pr label] [dict get $pr mode] [format %.2f [_mean $d]] [format %.2f [_std $d]] \
                [format %.2f [_min $d]] [format %.2f [_max $d]] [format %+.2f [expr {[lindex $d end] - [lindex $d 0]}]] \
                [expr {$fl ? "\u2691" : ""}]]
            set disp [expr {[dict get $pr mode] eq "ca" ? [list [lindex [dict get $pr idxA] 0] [lindex [dict get $pr idxB] 0]] : ""}]
            dict set out $i [list [list [dict get $pr label] $d] $row $disp]
        }
    }
    set series {}; set rows {}; set meta {}; set pl {}; set nflag 0
    set i 0
    foreach pr $pairs {
        lassign [dict get $out $i] s row disp
        lappend series $s; lappend rows $row
        set fl [expr {[lindex $row 7] ne ""}]
        incr nflag $fl
        lappend meta [dict create pair $i flag $fl]
        dict set pr disp $disp
        lappend pl $pr
        incr i
    }
    set sum [format "%d pairs; %d change by > %.1f \u00c5" [llength $rows] $nflag $thr]
    return [dict create series $series pairs $pl threshold $thr header [list Pair mode mean sd min max "\u0394" "\u2691"] \
        rows $rows rowmeta $meta summary $sum engine [expr {[llength $bk] && ![llength $nat] ? "backend" : ([llength $bk] ? "native+backend" : "native")}] \
        calc_ms [expr {[clock milliseconds] - [dict get $P t0]}]]
}

# ---------------------------------------------------------------------------------------------
# Pocket dynamics (backend: same maths as the PyMOL plugin)
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_pockets {P} {
    set mol [dict get $P mol]
    set c [context $mol 0]
    set files {}
    set a [_base_args $P files]
    foreach k {hot cp} { dict set a $k [dict get $c args $k] }
    set spec [dict get $P spec]; if {$spec eq ""} { set spec auto }
    dict set a spec [::ufv::json::raw [_jstr $spec]]
    # ligand pockets: residues within 5 A of each of the 4 largest ligands in the reference frame
    set lp {}
    if {[regexp -nocase {(^|;)\s*(auto|lig)\s*($|;)} $spec]} {
        foreach l [lrange [ligand_ids $mol 5] 0 3] {
            lassign $l rn ch ri
            set keys {}
            ::ufv::with_sel s $mol "name CA and same residue as ($::ufv::traj::POLY and within 5 of ([_lig_tok "lig:$rn/$ch/$ri"]))" {
                $s frame [dict get $P ref]; $s update
                foreach {kc kr} [join [$s get {chain resid}]] { lappend keys "\[[_jstr $kc],$kr\]" }
            }
            lappend lp "{\"label\":[_jstr "$rn $ch/$ri pocket"],\"keys\":\[[join $keys ,]\]}"
        }
    }
    dict set a ligpockets [::ufv::json::raw "\[[join $lp ,]\]"]
    return [_backend_run $a $files]
}

# ---------------------------------------------------------------------------------------------
# Contacts: ligand / partner / nearby
# ---------------------------------------------------------------------------------------------
# Candidate pad (PyMOL _prep_contacts): 2 x p95 Calpha displacement after fitting + 6 A (<= 25 A).
proc ::ufv::traj::_pad {P} {
    set files {}
    set a [_base_args $P files]
    set code [catch { _backend traj_pad $a } r]
    _cleanup $files
    if {$code} { error $r }
    return [dict get $r pad]
}

# Atom sets of a contact analysis (selected in the reference frame, like the PyMOL plugin):
#   -> dict selA selB same idxA idxB rowsA rowsB uniA uniB ours center cutoff
proc ::ufv::traj::contact_sets {P kind {pad ""}} {
    set mol [dict get $P mol]
    set c [dict get $P cutoff]
    set st [dict get $P ref]
    set POLY $::ufv::traj::POLY
    if {$pad eq ""} { set pad [_pad $P] }
    set S [dict create cutoff $c pad $pad]
    if {$kind eq "ligand"} {
        if {![_has_lig $mol]} { error "mol $mol has no ligands" }
        set ligs "$::ufv::traj::LIG and noh"
        set aset "same residue as ($POLY and noh and within [expr {$c + $pad}] of ($ligs))"
        set idx [::ufv::with_sel s $mol $aset { $s frame $st; $s update; $s get index }]
        set selA [expr {[llength $idx] ? "[_idxsel $idx] and noh" : "none"}]
        set selB $ligs
    } elseif {$kind eq "partner"} {
        set chains [::ufv::with_sel s $mol $POLY { lsort -unique [$s get chain] }]
        if {[llength $chains] < 2} { error "mol $mol has a single polymer chain (no partner / inter-chain contacts)" }
        set ours {}
        foreach x [::ufv::map::chains $mol] { if {$x in $chains} { lappend ours $x } }
        set sub [expr {[llength $ours] && [llength $ours] < [llength $chains]}]
        set parts {}
        foreach ch [expr {$sub ? $ours : $chains}] {
            set t [::ufv::map::chain_token $ch]
            set own "$POLY and chain $t"; set other "$POLY and not chain $t"
            lappend parts "(($own) and within [expr {$c + $pad}] of ($other))"
            if {$sub} { lappend parts "(($other) and within [expr {$c + $pad}] of ($own))" }
        }
        set idx [::ufv::with_sel s $mol "same residue as ([join $parts { or }])" { $s frame $st; $s update; $s get index }]
        set selB [expr {[llength $idx] ? "[_idxsel $idx] and noh" : "none"}]
        set selA [expr {$sub ? "$selB and chain [join [lmap x $ours {::ufv::map::chain_token $x}] { }]" : $selB}]
        dict set S ours $ours
    } else {
        set pos [dict get $P spec]
        if {$pos eq ""} { error "nearby: give a residue (UniProt position)" }
        set s0 [_res_tok $mol $pos]
        if {$s0 eq ""} { error "nearby: position $pos is not mapped on mol $mol" }
        set d [expr {double($::ufv::nearby)}]
        set c $d
        dict set S cutoff $d
        dict set S center $pos
        set idx [::ufv::with_sel s $mol "same residue as (noh and not water and not ($s0) and within [expr {$d + $pad}] of (($s0) and noh))" {
            $s frame $st; $s update; $s get index
        }]
        set selA "($s0) and noh"
        set selB [expr {[llength $idx] ? "[_idxsel $idx] and noh" : "none"}]
    }
    dict set S selA $selA; dict set S selB $selB
    dict set S same [expr {$selA eq $selB}]
    foreach k {A B} {
        set rows {}; set uni {}; set ix {}
        ::ufv::with_sel s $mol [dict get $S sel$k] {
            foreach {i ch ri ins rn nm} [join [$s get {index chain resid insertion resname name}]] {
                lappend ix $i
                lappend rows [list $ch "$ri[string trim $ins]" $ri $rn $nm]
                lappend uni [::ufv::map::r2u $mol $ch $ri]
            }
        }
        dict set S idx$k $ix; dict set S rows$k $rows; dict set S uni$k $uni
    }
    if {![llength [dict get $S idxA]] || ![llength [dict get $S idxB]]} { error "no $kind contacts possible in mol $mol" }
    if {![dict exists $S ours]} { dict set S ours [::ufv::map::chains $mol] }
    return $S
}

proc ::ufv::traj::_rows_json {rows} {
    set o {}
    foreach r $rows { lassign $r ch ri rv rn nm; lappend o "\[[_jstr $ch],[_jstr $ri],$rv,[_jstr $rn],[_jstr $nm]\]" }
    return "\[[join $o ,]\]"
}
proc ::ufv::traj::_uni_json {l} { return "\[[join [lmap u $l {expr {$u eq "" ? "null" : $u}}] ,]\]" }

# JSON object of a contact set for the backend (writes its DCDs; file paths appended to files)
# Backend request arguments of a contact set (writes its window DCDs; paths appended to files)
proc ::ufv::traj::_contact_args {P S filesvar} {
    upvar 1 $filesvar files
    set mol [dict get $P mol]
    lassign [_write_dcd $mol [dict get $S selA] [dict get $P frames] [dict get $P ref]] dA rA
    lappend files $dA $rA
    set a [dict create dcdA $dA rowsA [::ufv::json::raw [_rows_json [dict get $S rowsA]]] \
        uniA [::ufv::json::raw [_uni_json [dict get $S uniA]]]]
    if {$rA ne ""} { dict set a refA $rA }
    if {[dict get $S same]} {
        dict set a same [::ufv::json::raw true]
    } else {
        lassign [_write_dcd $mol [dict get $S selB] [dict get $P frames] [dict get $P ref]] dB rB
        lappend files $dB $rB
        dict set a dcdB $dB
        dict set a rowsB [::ufv::json::raw [_rows_json [dict get $S rowsB]]]
        dict set a uniB [::ufv::json::raw [_uni_json [dict get $S uniB]]]
        if {$rB ne ""} { dict set a refB $rB }
    }
    dict set a ours [::ufv::json::raw [_jstrs [dict get $S ours]]]
    dict set a cutoff [dict get $S cutoff]
    if {[dict exists $S center]} { dict set a center [::ufv::json::raw [_jstr [dict get $S center]]] }
    return $a
}

proc ::ufv::traj::_contacts {P} {
    set mol [dict get $P mol]
    set kind [dict get $P kind]
    set c [context $mol 0]
    set S [contact_sets $P $kind]
    if {[dict get $P engine] eq "backend"} {
        set files {}
        set a [_base_args $P files]
        set a [dict merge $a [_contact_args $P $S files]]
        foreach k {hot cp} { dict set a $k [dict get $c args $k] }
        return [_backend_run $a $files]
    }
    return [_contacts_native $P $S [dict get $c tags] $kind]
}

# Native contact occupancy: `measure contacts` per window frame on the same atom sets as the PyMOL plugin,
# reduced to residue pairs; table / series / summary built exactly like _calc_contacts.
proc ::ufv::traj::_contacts_native {P S tags kind} {
    variable RF
    set mol [dict get $P mol]
    set fr [dict get $P frames]
    set F [llength $fr]
    set cut [dict get $S cutoff]
    # residue ids: allkeys in order of first appearance in rowsA then rowsB
    set kid [dict create]; set allkeys {}; set info {}
    foreach side {A B} {
        foreach row [dict get $S rows$side] u [dict get $S uni$side] {
            set k "[lindex $row 0] [lindex $row 2]"
            if {![dict exists $kid $k]} {
                dict set kid $k [llength $allkeys]
                lappend allkeys $k
                lappend info [list [lindex $row 0] [lindex $row 1] [lindex $row 3] $u]
            }
        }
    }
    set nk [llength $allkeys]
    # atom index -> residue id of the A side / B side (lists indexed by atom index), chain per atom
    set natoms [molinfo $mol get numatoms]
    set kA [lrepeat $natoms -1]
    set kB [lrepeat $natoms -1]
    set a2ch [dict create]
    foreach side {A B} {
        set kv k$side
        foreach i [dict get $S idx$side] row [dict get $S rows$side] {
            lset $kv $i [dict get $kid "[lindex $row 0] [lindex $row 2]"]
            dict set a2ch $i [lindex $row 0]
        }
    }
    # atom-set pairs for `measure contacts`: partner = each of our chains vs the other chains (inter-chain only,
    # each chain pair once); ligand / nearby = A vs B
    set pairsels {}
    set code [catch {
        if {$kind eq "partner"} {
            set chs [expr {[dict get $S same] ? [lsort -unique [lmap r [dict get $S rowsB] { lindex $r 0 }]] : [dict get $S ours]}]
            set done {}
            foreach ch $chs {
                set t [::ufv::map::chain_token $ch]
                lappend done $t
                set s1 [atomselect $mol "([dict get $S selA]) and chain $t"]
                set s2 [atomselect $mol "([dict get $S selB]) and not chain [join $done { }]"]
                lappend pairsels [list $s1 $s2]
            }
        } else {
            lappend pairsels [list [atomselect $mol [dict get $S selA]] [atomselect $mol [dict get $S selB]]]
        }
        # bonded A-B pairs (`measure contacts` skips bonded atoms; the PyMOL plugin counts them): checked per frame
        set bp {}
        foreach ps $pairsels {
            lassign $ps s1 s2
            set in2 [dict create]; foreach i [$s2 get index] { dict set in2 $i 1 }
            foreach i [$s1 get index] bl [$s1 getbonds] {
                foreach j $bl { if {[dict exists $in2 $j]} { lappend bp [list $i $j] } }
            }
        }
        set bsel ""
        if {[llength $bp]} {
            set bidx [lsort -integer -unique [concat {*}$bp]]
            set bsel [atomselect $mol "index [::ufv::map::compress $bidx]"]
            set bcol [dict create]; set j 0; foreach i $bidx { dict set bcol $i $j; incr j }
        }
    } err]
    if {$code} { foreach ps $pairsels { foreach s $ps { catch { $s delete } } }; error $err }
    set perframe {}
    set cnt [dict create]
    set natc 0
    set c2 [expr {$cut * $cut}]
    set code [catch {
        set t 0
        foreach f $fr {
            set la {}; set lb {}
            foreach ps $pairsels {
                lassign $ps s1 s2
                if {![$s1 num] || ![$s2 num]} continue
                $s1 frame $f; $s2 frame $f
                lassign [measure contacts $cut $s1 $s2] x y
                lappend la {*}$x; lappend lb {*}$y
            }
            if {$bsel ne ""} {
                $bsel frame $f
                set X [$bsel get {x y z}]
                foreach p $bp {
                    lassign $p a b
                    set d [vecsub [lindex $X [dict get $bcol $a]] [lindex $X [dict get $bcol $b]]]
                    if {[vecdot $d $d] <= $c2} { lappend la $a; lappend lb $b }
                }
            }
            incr natc [llength $la]
            if {$kind eq "partner"} {
                set keys [lmap a $la b $lb {
                    set x [lindex $kB $a]; set y [lindex $kB $b]
                    expr {$x < $y ? $x * $nk + $y : $y * $nk + $x}
                }]
            } else {
                set keys [lmap a $la b $lb {
                    set x [lindex $kA $a]
                    expr {$x >= 0 ? $x * $nk + [lindex $kB $b] : [lindex $kA $b] * $nk + [lindex $kB $a]}
                }]
            }
            set keys [lsort -integer -unique $keys]
            lappend perframe $keys
            foreach k $keys { dict incr cnt $k }
            if {[incr t] % 25 == 0} { _slice }
        }
    } err]
    foreach ps $pairsels { foreach s $ps { catch { $s delete } } }
    if {$bsel ne ""} { $bsel delete }
    if {$code} { error $err }
    # residue-pair occupancy (pairs sorted by key, like np.unique)
    set upk [lsort -integer [dict keys $cnt]]
    set partners [dict create]
    foreach k $upk {
        set a [expr {$k / $nk}]; set b [expr {$k % $nk}]
        set oc [expr {double([dict get $cnt $k]) / $F}]
        dict lappend partners $a [list $oc $b]
        if {$kind eq "partner"} { dict lappend partners $b [list $oc $a] }
    }
    # per-residue touched frames
    set touched [dict create]
    set t 0
    foreach keys $perframe {
        set rs [dict create]
        foreach k $keys {
            dict set rs [expr {$k / $nk}] 1
            if {$kind eq "partner"} { dict set rs [expr {$k % $nk}] 1 }
        }
        foreach r [dict keys $rs] { dict lappend touched $r $t }
        incr t
    }
    set ours [dict get $S ours]
    set rows {}; set meta {}; set occd [dict create]; set rframes [dict create]
    foreach r [lsort -integer [dict keys $touched]] {
        set idx [dict get $touched $r]
        set o [expr {double([llength $idx]) / $F}]
        if {$o <= 0} continue
        lassign [lindex $info $r] ch resi resn u
        if {$kind eq "partner" && $ch ni $ours} continue
        set key [lindex $allkeys $r]
        dict set occd $key $o
        dict set rframes $key $idx
        set best [lrange [lsort -real -decreasing -index 0 [lsort -integer -decreasing -index 1 [dict get $partners $r]]] 0 2]
        set top [join [lmap x $best { format "%s %s" [_lbl $info [lindex $x 1]] [_pct [lindex $x 0]] }] "; "]
        lappend rows [list $u "$ch/$resi" $resn [_pct $o] [_switches $idx $F] $top \
            [expr {$u ne "" && [dict exists $tags $u] ? [dict get $tags $u] : ""}]]
        lappend meta [dict create pos $u key $key chain $ch resi $resi occ $o]
    }
    # order: occupancy desc, then "chain/resi" (stable sorts)
    set o1 [lsort -ascii -indices [lmap x $rows { lindex $x 1 }]]
    set o2 [lsort -real -decreasing -indices [lmap i $o1 { dict get [lindex $meta $i] occ }]]
    set ord [lmap j $o2 { lindex $o1 $j }]
    set rows [lmap i $ord { lindex $rows $i }]
    set meta [lmap i $ord { lindex $meta $i }]
    # series
    set series {}
    if {$kind eq "ligand"} {
        set bc [dict create]
        foreach k $upk { dict incr bc [expr {$k % $nk}] }
        set bl {}; dict for {b n} $bc { lappend bl [list $b $n] }
        foreach x [lrange [lsort -integer -decreasing -index 1 $bl] 0 5] {
            set b [lindex $x 0]
            set s {}
            foreach keys $perframe { set n 0; foreach k $keys { if {$k % $nk == $b} { incr n } }; lappend s [expr {double($n)}] }
            lappend series [list "[_lbl $info $b] contacts" $s]
        }
    } else {
        set s {}; foreach keys $perframe { lappend s [expr {double([llength $keys])}] }
        lappend series [list [expr {$kind eq "partner" ? "inter-chain residue pairs" : [format "residues within %.1f \u00c5" $cut]}] $s]
    }
    set head [dict get {ligand ligand partner {other chain} nearby neighbour} $kind]
    set sum [format "%d residues in contact (heavy atoms \u2264 %.1f \u00c5)%s; %d atom contacts (measure contacts)" \
        [llength $rows] $cut [expr {$kind eq "nearby" ? " of residue [dict get $S center]" : ""}] $natc]
    set RF($mol,$kind) [list [list $fr [dict get $P ref] [dict get $P cutoff]] $rframes]
    set occl {}; dict for {k o} $occd { lappend occl [list {*}$k $o] }
    return [dict create series $series occupancy $occl cutoff $cut header [list UniProt res aa occ switch "top $head" tags] \
        rows [lrange $rows 0 399] rowmeta [lrange $meta 0 399] summary $sum engine native \
        calc_ms [expr {[clock milliseconds] - [dict get $P t0]}]]
}
proc ::ufv::traj::_lbl {info i} {
    lassign [lindex $info $i] ch resi resn
    return "$resn $ch/$resi"
}

# ---------------------------------------------------------------------------------------------
# Annotation dynamics (ligand / partner occupancy of annotated residues; reuses the last contact runs)
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_annot {P} {
    if {[dict get $P engine] eq "backend"} { return [_annot_backend $P] }
    return [_annot_native $P]
}

# Backend annotation dynamics: the backend reuses its last ligand / partner runs of the same window.
proc ::ufv::traj::_annot_backend {P} {
    set mol [dict get $P mol]
    set c [context $mol 0]
    set files {}
    set a [_base_args $P files]
    foreach k {hot cp} { dict set a $k [dict get $c args $k] }
    set pad ""
    foreach k {ligand partner} {
        set last [result $mol $k]
        if {$last ne "" && [dict get $last engine] eq "backend" && [dict get $last frames] eq [dict get $P frames] \
                && [dict get $last ref] == [dict get $P ref] && [dict get $last cutoff] == [dict get $P cutoff]} {
            continue
        }
        if {[catch {
            if {$pad eq ""} { set pad [_pad $P] }
            set Q [dict replace $P kind $k]
            set sub [::ufv::json::obj [_contact_args $Q [contact_sets $Q $k $pad] files]]
        } err]} {
            set sub [::ufv::json::obj [dict create error $err]]
        }
        dict set a $k [::ufv::json::raw $sub]
    }
    return [_backend_run $a $files]
}

proc ::ufv::traj::_annot_native {P} {
    variable RF
    set mol [dict get $P mol]
    set fr [dict get $P frames]; set F [llength $fr]
    set c [context $mol 0]
    set tags [dict get $c tags]
    set sig [list $fr [dict get $P ref] [dict get $P cutoff]]
    set sub [dict create]
    set pad ""
    foreach k {ligand partner} {
        if {[info exists RF($mol,$k)] && [lindex $RF($mol,$k) 0] eq $sig} {
            dict set sub $k [lindex $RF($mol,$k) 1]; continue
        }
        if {[catch {
            if {$pad eq ""} { set pad [_pad $P] }
            set Q [dict replace $P kind $k]
            set S [contact_sets $Q $k $pad]
            set r [_contacts_native $Q $S $tags $k]
            _store $mol $k [dict merge $r [dict create kind $k mol $mol frames $fr ref [dict get $P ref] res [dict get $P ca res]]]
            dict set sub $k [lindex $RF($mol,$k) 1]
        } err]} {
            dict set sub $k ""
        }
    }
    set rows {}; set meta {}; set ser {}
    foreach r [dict get $P ca res] {
        lassign $r ch resi resv resn u
        if {$u eq "" || ![dict exists $tags $u]} continue
        set vals {}; set any 0; set sw 0; set union {}
        foreach k {ligand partner} {
            set d [dict get $sub $k]
            if {$d eq ""} { lappend vals ""; continue }
            set idx [expr {[dict exists $d "$ch $resv"] ? [dict get $d "$ch $resv"] : {}}]
            set o [expr {double([llength $idx]) / $F}]
            lappend vals $o
            if {$o > 0} { set any 1 }
            incr sw [_switches $idx $F]
            lappend union {*}$idx
        }
        if {!$any} continue
        lappend rows [list $u "$ch/$resi" $resn [dict get $tags $u] \
            [expr {[lindex $vals 0] eq "" ? "-" : [_pct [lindex $vals 0]]}] \
            [expr {[lindex $vals 1] eq "" ? "-" : [_pct [lindex $vals 1]]}] $sw]
        lappend meta [dict create pos $u key "$ch $resv"]
        lappend ser [lsort -integer -unique $union]
    }
    set series {}
    if {[llength $ser]} {
        set cnt [lrepeat $F 0]
        foreach l $ser { foreach t $l { lset cnt $t [expr {[lindex $cnt $t] + 1}] } }
        lappend series [list "annotated residues in contact" [lmap x $cnt { expr {double($x)} }]]
    }
    set key {}
    foreach row $rows {
        set a [expr {[lindex $row 4] eq "-" ? 0.0 : double([string trimright [lindex $row 4] %])}]
        set b [expr {[lindex $row 5] eq "-" ? 0.0 : double([string trimright [lindex $row 5] %])}]
        lappend key [expr {$a + $b}]
    }
    set order [lsort -real -decreasing -indices $key]
    set none [expr {[dict get $sub ligand] eq "" && [dict get $sub partner] eq ""}]
    return [dict create series $series header {UniProt res aa tags ligand partner switch} \
        rows [lmap i $order { lindex $rows $i }] rowmeta [lmap i $order { lindex $meta $i }] \
        summary [format "%d annotated residues contact a ligand / another chain%s" [llength $rows] \
            [expr {$none ? " (no ligand and no partner chain in this object)" : ""}]] engine native \
        calc_ms [expr {[clock milliseconds] - [dict get $P t0]}]]
}

# ---------------------------------------------------------------------------------------------
# Apply: colours (user field of every frame + one continuous mode), distance display
# ---------------------------------------------------------------------------------------------
# Write per-residue values {"chain resv" value} into `user` of every polymer atom in EVERY frame (user is
# per-frame in VMD), others 0. The value list is built once; each frame costs one `$sel set` (~0.2 ms).
proc ::ufv::traj::write_user {mol keyval {default 0.0}} {
    ::ufv::with_sel s $mol $::ufv::traj::POLY {
        set vals {}
        foreach {ch ri} [join [$s get {chain resid}]] {
            lappend vals [expr {[dict exists $keyval "$ch $ri"] ? [dict get $keyval "$ch $ri"] : $default}]
        }
        set n [nframes $mol]
        for {set f 0} {$f < $n} {incr f} {
            $s frame $f
            $s set user $vals
            if {$f % 100 == 99} { _slice }
        }
    }
}

proc ::ufv::traj::apply {r {color 1} {show 1}} {
    set mol [dict get $r mol]
    set kind [dict get $r kind]
    if {$kind eq "rmsf" && $color} {
        set kv [dict create]
        foreach rr [dict get $r res] v [dict get $r rmsf] { dict set kv "[lindex $rr 0] [lindex $rr 2]" $v }
        write_user $mol $kv
        set hi [_percentile [dict get $r rmsf] 98]
        if {$hi <= 0} { set hi 1.0 }
        ::ufv::cartoon::continuous $mol traj_rmsf User 0.0 $hi BWR    ;# BWR borrowed only while shown
        ::ufv::status "mol $mol: cartoon coloured by RMSF (blue rigid -> red mobile, 0-[format %.2f $hi] \u00c5)."
    } elseif {$kind in {ligand partner} && $color} {
        set kv [dict create]
        foreach o [dict get $r occupancy] { dict set kv "[lindex $o 0] [lindex $o 1]" [lindex $o 2] }
        write_user $mol $kv
        ::ufv::cartoon::continuous $mol "traj_$kind" User -1.0 1.0 BWR
        ::ufv::status "mol $mol: cartoon coloured by $kind contact occupancy (white 0 -> red 100%)."
    } elseif {$kind eq "dist" && $show} {
        set rng {}
        foreach s [dict get $r series] { set v [lindex $s 1]; lappend rng [expr {[_max $v] - [_min $v]}] }
        set order [lsort -real -decreasing -indices $rng]
        set fl {}; foreach i $order { if {[dict get [lindex [dict get $r rowmeta] $i] flag]} { lappend fl $i } }
        set which [lrange [expr {[llength $fl] ? $fl : $order}] 0 [expr {$::ufv::traj::cfg(show_pairs) - 1}]]
        show_pairs $mol $which $r
    }
    return $r
}

# Distance display for the chosen pairs (indices into the last dist result): atom pairs as VMD Bonds labels
# (they follow the frame natively), centroid pairs as ONE graphics set redrawn by the debounced frame trace.
proc ::ufv::traj::show_pairs {mol which {r ""}} {
    variable cen
    if {$r eq ""} { set r [result $mol dist] }
    if {$r eq ""} { error "run ufv_traj_dist first" }
    variable labels
    clear_pairs $mol
    set cl {}
    set labels($mol) {}
    foreach i [lrange $which 0 11] {
        set pr [lindex [dict get $r pairs] $i]
        set d [dict get $pr disp]
        if {[llength $d] == 2} {
            label add Bonds $mol/[lindex $d 0] $mol/[lindex $d 1]
            lappend labels($mol) [lsort -integer $d]
        } else {
            lappend cl [list [dict get $pr idxA] [dict get $pr idxB]]
        }
    }
    set cen($mol) $cl
    if {[llength $cl]} { _trace_on; _draw_cen $mol }
    return [llength $which]
}

# Remove this plugin's distance display (only labels between pairs of the last dist result; other labels stay).
proc ::ufv::traj::clear_pairs {mol} {
    variable cen
    variable gfx
    variable labels
    if {[info exists labels($mol)] && [llength $labels($mol)]} {
        set mine [dict create]
        foreach d $labels($mol) { dict set mine $d 1 }
        set labels($mol) {}
        set i 0; set del {}
        foreach l [label list Bonds] {
            lassign $l a b
            if {[lindex $a 0] == $mol && [lindex $b 0] == $mol && [dict exists $mine [lsort -integer [list [lindex $a 1] [lindex $b 1]]]]} { lappend del $i }
            incr i
        }
        foreach i [lsort -integer -decreasing $del] { label delete Bonds $i }
    }
    if {[info exists gfx($mol)]} { foreach g $gfx($mol) { catch { graphics $mol delete $g } } }
    set gfx($mol) {}
    set cen($mol) {}
}

proc ::ufv::traj::_draw_cen {mol} {
    variable cen
    variable gfx
    if {$mol ni [molinfo list]} return
    if {[info exists gfx($mol)]} { foreach g $gfx($mol) { catch { graphics $mol delete $g } } }
    set gfx($mol) {}
    if {![info exists cen($mol)] || ![llength $cen($mol)]} return
    set f [molinfo $mol get frame]
    lappend gfx($mol) [graphics $mol color orange]
    foreach p $cen($mol) {
        lassign $p A B
        set pa [::ufv::with_sel s $mol [_idxsel $A] { $s frame $f; measure center $s }]
        set pb [::ufv::with_sel s $mol [_idxsel $B] { $s frame $f; measure center $s }]
        lappend gfx($mol) [graphics $mol line $pa $pb width 3 style dashed]
        lappend gfx($mol) [graphics $mol text [vecscale 0.5 [vecadd $pa $pb]] [format "%.2f" [veclength [vecsub $pa $pb]]] size 1.0]
    }
}

# ---------------------------------------------------------------------------------------------
# Frame trace (debounced): centroid-distance graphics + panel cursor. Installed only when needed.
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::_trace_on {} {
    if {![llength [info commands ::molinfo]]} return
    if {[lsearch -exact [trace info variable ::vmd_frame] {write ::ufv::traj::_ontrace}] < 0} {
        trace add variable ::vmd_frame write ::ufv::traj::_ontrace
    }
}
proc ::ufv::traj::_ontrace {n1 n2 op} {
    variable trace_after
    if {$trace_after ne ""} { after cancel $trace_after }
    set trace_after [after 30 ::ufv::traj::_onframe]
}
proc ::ufv::traj::_onframe {} {
    variable trace_after ""
    variable cen
    foreach mol [array names cen] { if {[llength $cen($mol)]} { catch { _draw_cen $mol } } }
    catch { _panel_cursor }
}

# Layers that depend on the frame (focus / Nearby `within` selections) are re-evaluated on every frame.
proc ::ufv::traj::follow {mol} {
    if {[nframes $mol] < 2} { return 0 }
    set n 0
    foreach name [::ufv::layers::names $mol] {
        foreach i [::ufv::layers::reps $mol $name] {
            if {[regexp {within|same residue} [lindex [molinfo $mol get "{selection $i}"] 0]]} {
                mol selupdate $i $mol 1; incr n
            }
        }
    }
    return $n
}
proc ::ufv::traj::_onlayers {mol args} { if {$mol in [molinfo list]} { catch { follow $mol } } }
::ufv::hook add layers ::ufv::traj::_onlayers
proc ::ufv::traj::_ontraj_read {n1 n2 op} {
    catch { follow $n2 }
    catch { _panel_refresh }
}
if {[llength [info commands ::molinfo]] && ![info exists ::ufv::traj::_rt_on]} {
    set ::ufv::traj::_rt_on 1
    trace add variable ::vmd_trajectory_read write ::ufv::traj::_ontraj_read
}

# ---------------------------------------------------------------------------------------------
# Structure tiers on one frame (cached per frame) / CSV / B-factor PDB
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::tiers {mol kind {frame ""}} {
    variable tiercache
    set fn [dict get {hotspots ufv_hotspots hubs ufv_contacthubs contacthubs ufv_contacthubs pockets ufv_pockets} $kind]
    if {$frame ne "" && $frame != [molinfo $mol get frame]} { animate goto $frame }
    set fr [molinfo $mol get frame]
    set uid [::ufv::mol_uid $mol]
    set ver [::ufv::st::get $mol mapver 0]
    set k "$mol,$uid,$ver,$fr,"
    if {[info exists tiercache($k)]} { ::ufv::st::put $mol {cache tiers} [list [list $uid $ver $fr ""] $tiercache($k)] }
    set n [$fn $mol]
    set c [::ufv::st::get $mol {cache tiers} ""]
    if {$c ne ""} { set tiercache($k) [lindex $c 1] }
    ::ufv::status "(computed on frame $fr of mol $mol)"
    return $fr
}

proc ::ufv::traj::csv {mol prefix} {
    variable KINDS
    set res [results $mol]
    if {![dict size $res]} { error "no trajectory results for mol $mol yet" }
    set cols {}; set data [dict create]; set frs [dict create]
    foreach kind $KINDS {
        if {![dict exists $res $kind]} continue
        set r [dict get $res $kind]
        foreach s [dict get $r series] {
            lassign $s name arr
            set c "$kind:$name"
            lappend cols $c
            foreach f [dict get $r frames] v $arr { if {$v ne ""} { dict set data "$f\t$c" $v; dict set frs $f 1 } }
        }
    }
    set f1 "${prefix}_frames.csv"
    set fh [open $f1 w]; fconfigure $fh -encoding utf-8
    puts $fh [_csvrow [concat frame $cols]]
    foreach f [lsort -integer [dict keys $frs]] {
        set row [list $f]
        foreach c $cols { lappend row [expr {[dict exists $data "$f\t$c"] ? [format %.4f [dict get $data "$f\t$c"]] : ""}] }
        puts $fh [_csvrow $row]
    }
    close $fh
    set info [dict create]; set rmsf [dict create]; set pers [dict create]
    if {[dict exists $res rmsf]} {
        set r [dict get $res rmsf]
        foreach rr [dict get $r res] v [dict get $r rmsf] {
            lassign $rr c ri rv rn u
            dict set rmsf "$c $rv" $v
            dict set info "$c $rv" [list $c $ri $rn $u]
        }
    }
    set occ [dict create ligand {} partner {}]
    foreach k {ligand partner} {
        if {[dict exists $res $k]} {
            set d [dict create]
            foreach o [dict get $res $k occupancy] { dict set d "[lindex $o 0] [lindex $o 1]" [lindex $o 2] }
            dict set occ $k $d
        }
    }
    if {[dict exists $res pockets]} {
        set r [dict get $res pockets]
        dict for {lab l} [dict get $r persist] {
            foreach x $l {
                lassign $x col f
                lassign [lindex [dict get $r res] $col] c ri rv rn u
                if {![dict exists $info "$c $rv"]} { dict set info "$c $rv" [list $c $ri $rn $u] }
                dict lappend pers "$c $rv" [format "%s=%.2f" $lab $f]
            }
        }
    }
    set tags [dict get [context $mol 0] tags]
    foreach k {ligand partner} {
        dict for {key o} [dict get $occ $k] {
            if {![dict exists $info $key]} { dict set info $key [list [lindex $key 0] [lindex $key 1] "" [::ufv::map::r2u $mol {*}$key]] }
        }
    }
    set keys [lsort -ascii -index 0 [lsort -integer -index 1 [dict keys $info]]]
    set f2 "${prefix}_residues.csv"
    set fh [open $f2 w]; fconfigure $fh -encoding utf-8
    puts $fh [_csvrow {chain resi resn uniprot_pos rmsf ligand_occupancy partner_occupancy pocket_persistence tags}]
    foreach k $keys {
        lassign [dict get $info $k] c ri rn u
        set lo [dict get $occ ligand]; set po [dict get $occ partner]
        puts $fh [_csvrow [list $c $ri $rn $u \
            [expr {[dict exists $rmsf $k] ? [format %.3f [dict get $rmsf $k]] : ""}] \
            [expr {[dict exists $lo $k] ? [format %.3f [dict get $lo $k]] : ""}] \
            [expr {[dict exists $po $k] ? [format %.3f [dict get $po $k]] : ""}] \
            [expr {[dict exists $pers $k] ? [join [dict get $pers $k] ";"] : ""}] \
            [expr {$u ne "" && [dict exists $tags $u] ? [dict get $tags $u] : ""}]]]
    }
    close $fh
    ::ufv::status "wrote $f1 and $f2."
    return [list $f1 $f2]
}
proc ::ufv::traj::_csvrow {l} {
    set o {}
    foreach v $l {
        if {[regexp {[,"\n]} $v]} { set v "\"[string map {\" \"\"} $v]\"" }
        lappend o $v
    }
    return [join $o ,]
}

# B-factor PDB of the current frame: beta = RMSF (A) | ligand / partner occupancy (0-1) per residue; the
# molecule's own beta values are restored afterwards.
proc ::ufv::traj::bfactor {mol path {mode rmsf}} {
    set kv [dict create]
    if {$mode eq "rmsf"} {
        set r [result $mol rmsf]
        if {$r eq ""} { error "run ufv_traj_rmsf first" }
        foreach rr [dict get $r res] v [dict get $r rmsf] { dict set kv "[lindex $rr 0] [lindex $rr 2]" $v }
    } else {
        set k [dict get {ligocc ligand partnerocc partner ligand ligand partner partner} $mode]
        set r [result $mol $k]
        if {$r eq ""} { error "run ufv_traj_contacts $k first" }
        foreach o [dict get $r occupancy] { dict set kv "[lindex $o 0] [lindex $o 1]" [lindex $o 2] }
    }
    ::ufv::with_sel s $mol all {
        set old [$s get beta]
        set vals {}
        foreach {ch ri} [join [$s get {chain resid}]] {
            lappend vals [expr {[dict exists $kv "$ch $ri"] ? [dict get $kv "$ch $ri"] : 0.0}]
        }
        $s set beta $vals
        set code [catch { $s writepdb $path } err]
        $s set beta $old
        if {$code} { error $err }
    }
    ::ufv::status "wrote $path (beta = $mode)."
    return $path
}

# ---------------------------------------------------------------------------------------------
# Console commands (PyMOL ufv_traj_* parity). molid is optional and last everywhere.
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::print {r} {
    set fr [dict get $r frames]
    puts [format "\[UFV\] mol %s %s: %d frames (%d-%d), ref %d - %s  \[%s, %d ms\]" [dict get $r mol] [dict get $r kind] \
        [llength $fr] [lindex $fr 0] [lindex $fr end] [dict get $r ref] [dict get $r summary] [dict get $r engine] \
        [expr {[dict exists $r wall_ms] ? [dict get $r wall_ms] : 0}]]
    set rows [dict get $r rows]
    if {[llength $rows]} {
        puts "      [join [dict get $r header] { | }]"
        foreach row [lrange $rows 0 19] { puts "      [join $row { | }]" }
        if {[llength $rows] > 20} { puts "      ... [expr {[llength $rows] - 20}] more rows (ufv_traj_csv)" }
    }
}
proc ::ufv::traj::_cmd {kind spec opts molid} {
    set mol [::ufv::resolve $molid]
    set r [run $mol $kind $spec $opts]
    print $r
    return $r
}

proc ufv_traj {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set fr [::ufv::traj::frames $mol]
    set msg [format "mol %s: %d frame(s), showing %d; window %d-%d step %s (%d frames), ref %d; engine %s; results: %s" \
        $mol [::ufv::traj::nframes $mol] [molinfo $mol get frame] [lindex $fr 0] [lindex $fr end] \
        $::ufv::traj::cfg(stride) [llength $fr] [::ufv::traj::ref $mol $fr] $::ufv::traj::cfg(engine) \
        [expr {[dict size [::ufv::traj::results $mol]] ? [join [lsort [dict keys [::ufv::traj::results $mol]]] {, }] : "none"}]]
    ::ufv::status $msg
    return [::ufv::traj::nframes $mol]
}
# ufv_traj_frames ?start? ?stop|end? ?stride? ?ref|first?
proc ufv_traj_frames {{start ""} {stop ""} {stride ""} {ref ""}} {
    foreach k {start stop stride ref} {
        set v [set $k]
        if {$v eq ""} continue
        if {$k eq "stop" && $v in {end last 0 -1}} { set v end }
        if {$k eq "ref" && $v in {first -1}} { set v first }
        if {$v ni {end first} && ![string is integer -strict $v]} { error "ufv_traj_frames: $k must be an integer (or end / first)" }
        set ::ufv::traj::cfg($k) $v
    }
    set ::ufv::traj::cfg(stride) [expr {max(1, $::ufv::traj::cfg(stride))}]
    set ::ufv::traj::cfg(start) [expr {max(0, $::ufv::traj::cfg(start))}]
    ::ufv::status "trajectory window: start $::ufv::traj::cfg(start), stop $::ufv::traj::cfg(stop), stride $::ufv::traj::cfg(stride), ref $::ufv::traj::cfg(ref)"
    return [array get ::ufv::traj::cfg]
}
proc ufv_traj_engine {{engine ""}} {
    if {$engine eq ""} { return $::ufv::traj::cfg(engine) }
    if {$engine ni {auto native backend}} { error "engine: auto | native | backend" }
    set ::ufv::traj::cfg(engine) $engine
}
proc ufv_traj_rmsf {{molid ""} {color 1}} { ::ufv::traj::_cmd rmsf "" [dict create color $color] $molid }
proc ufv_traj_dist {{pairs auto} {threshold ""} {show 1} {molid ""}} {
    ::ufv::traj::_cmd dist $pairs [dict create threshold $threshold show $show] $molid
}
proc ufv_traj_pockets {{pockets auto} {molid ""}} { ::ufv::traj::_cmd pockets $pockets {} $molid }
proc ufv_traj_contacts {{kind ligand} {cutoff ""} {color 1} {molid ""}} {
    if {$kind ni {ligand partner}} { error "ufv_traj_contacts ligand|partner ?cutoff? ?color? ?molid?" }
    ::ufv::traj::_cmd $kind "" [dict create cutoff $cutoff color $color] $molid
}
proc ufv_traj_annot {{molid ""}} { ::ufv::traj::_cmd annot "" {} $molid }
proc ufv_traj_nearby {pos {molid ""}} { ::ufv::traj::_cmd nearby $pos {} $molid }
proc ufv_traj_tiers {{kind hotspots} {frame ""} {molid ""}} {
    if {$kind ni {hotspots hubs contacthubs pockets}} { error "ufv_traj_tiers hotspots|hubs|pockets ?frame? ?molid?" }
    return [::ufv::traj::tiers [::ufv::resolve $molid] $kind $frame]
}
proc ufv_traj_csv {{prefix ""} {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {$prefix eq ""} { set prefix [file join [pwd] "mol${mol}_traj"] }
    return [::ufv::traj::csv $mol $prefix]
}
proc ufv_traj_bfactor {{path ""} {mode rmsf} {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {$path eq ""} { set path [file join [pwd] "mol${mol}_$mode.pdb"] }
    return [::ufv::traj::bfactor $mol $path $mode]
}
proc ufv_traj_show {{pairs ""} {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {$pairs in {"" none off}} { ::ufv::traj::clear_pairs $mol; return 0 }
    return [::ufv::traj::show_pairs $mol $pairs]
}
proc ufv_traj_cancel {} { ::ufv::traj::cancel }

# ---------------------------------------------------------------------------------------------
# Tk panel: ::ufv::traj::panel <parent>  -> the panel frame (embedded by ufv_gui.tcl)
# ---------------------------------------------------------------------------------------------
proc ::ufv::traj::panel {parent} {
    variable pw
    variable P
    variable LABELS
    if {![llength [info commands ::winfo]]} { error "the trajectory panel needs Tk" }
    set w [expr {$parent eq "." ? ".ufvtraj" : "$parent.ufvtraj"}]
    if {[winfo exists $w]} { destroy $w }
    set pw $w
    array unset P
    set P(kind) rmsf; set P(param) ""; set P(series) ""; set P(data) ""; set P(mol) ""; set P(job) ""
    set P(nf) ""; set P(hint) ""
    # flat white look (matches the ufv_gui panel it is embedded in)
    ::ttk::style configure UfvTraj.TFrame -background white
    ::ttk::style configure UfvTraj.TLabel -background white
    ::ttk::frame $w -padding 2 -style UfvTraj.TFrame
    # window
    set f [::ttk::frame $w.win -style UfvTraj.TFrame]
    ::ttk::label $f.l -text "Frames" -style UfvTraj.TLabel
    ::ttk::entry $f.s -width 5 -textvariable ::ufv::traj::cfg(start)
    ::ttk::label $f.t -text "to" -style UfvTraj.TLabel
    ::ttk::entry $f.e -width 5 -textvariable ::ufv::traj::cfg(stop)
    ::ttk::label $f.k -text "step" -style UfvTraj.TLabel
    ::ttk::spinbox $f.st -from 1 -to 100000 -width 4 -textvariable ::ufv::traj::cfg(stride)
    ::ttk::label $f.r -text "ref" -style UfvTraj.TLabel
    ::ttk::entry $f.re -width 5 -textvariable ::ufv::traj::cfg(ref)
    ::ttk::label $f.n -textvariable ::ufv::traj::P(nf) -foreground #666 -style UfvTraj.TLabel
    pack $f.l $f.s $f.t $f.e $f.k $f.st $f.r $f.re $f.n -side left -padx 1
    # analysis + actions
    set a [::ttk::frame $w.run -style UfvTraj.TFrame]
    set names {}; foreach {lab k hint} $LABELS { lappend names $lab }
    ::ttk::combobox $a.kind -state readonly -width 18 -values $names
    $a.kind set [lindex $names 0]
    bind $a.kind <<ComboboxSelected>> ::ufv::traj::_panel_kind
    ::ttk::entry $a.param -width 12 -textvariable ::ufv::traj::P(param)
    ::ttk::button $a.go -text Run -width 6 -command ::ufv::traj::_panel_run
    ::ttk::button $a.csv -text CSV -width 4 -command ::ufv::traj::_panel_csv
    ::ttk::button $a.b -text B-PDB -width 6 -command ::ufv::traj::_panel_bpdb
    pack $a.kind $a.param -side left -padx 1
    pack $a.b $a.csv $a.go -side right -padx 1
    ::ttk::label $w.hint -textvariable ::ufv::traj::P(hint) -foreground #666 -style UfvTraj.TLabel
    ::ttk::label $w.sum -text "" -wraplength 340 -justify left -style UfvTraj.TLabel
    # plot
    set s [::ttk::frame $w.ser -style UfvTraj.TFrame]
    ::ttk::label $s.l -text "Plot" -style UfvTraj.TLabel
    ::ttk::combobox $s.c -state readonly -width 30
    bind $s.c <<ComboboxSelected>> ::ufv::traj::_panel_plot
    pack $s.l $s.c -side left -padx 1
    canvas $w.plot -height 120 -width 340 -background white -highlightthickness 1 -highlightbackground #e0e5e8
    bind $w.plot <Button-1> [list ::ufv::traj::_panel_click %x]
    bind $w.plot <B1-Motion> [list ::ufv::traj::_panel_click %x]
    bind $w.plot <Configure> ::ufv::traj::_panel_plot
    # table
    set t [::ttk::frame $w.tab -style UfvTraj.TFrame]
    ::ttk::treeview $t.tv -show headings -height 5 -selectmode browse -yscrollcommand [list $t.sb set]
    ::ttk::scrollbar $t.sb -orient vertical -command [list $t.tv yview]
    bind $t.tv <<TreeviewSelect>> ::ufv::traj::_panel_row
    pack $t.sb -side right -fill y
    pack $t.tv -side left -fill both -expand 1
    grid $w.win -sticky w
    grid $w.run -sticky ew
    grid $w.hint -sticky w
    grid $w.sum -sticky w
    grid $w.ser -sticky w
    grid $w.plot -sticky ew
    grid $w.tab -sticky nsew
    grid columnconfigure $w 0 -weight 1
    grid rowconfigure $w 6 -weight 1
    _trace_on
    ::ufv::hook add traj ::ufv::traj::_panel_result
    ::ufv::hook add loaded ::ufv::traj::_panel_loaded
    _panel_kind
    _panel_refresh
    return $w
}
proc ufv_traj_gui {} {
    if {![llength [info commands ::toplevel]]} { error "ufv_traj_gui needs Tk" }
    if {[winfo exists .ufvtrajtop]} { destroy .ufvtrajtop }
    toplevel .ufvtrajtop -bg white
    wm title .ufvtrajtop "UFV trajectory"
    pack [::ufv::traj::panel .ufvtrajtop] -fill both -expand 1
    # wheel over the panel's combos / spins never changes their value (ufv_gui helper, when loaded)
    catch { ::ufv::gui::wheelfix .ufvtrajtop }
    return .ufvtrajtop
}

proc ::ufv::traj::_panel_ok {} { variable pw; return [expr {$pw ne "" && [llength [info commands ::winfo]] && [winfo exists $pw]}] }
proc ::ufv::traj::_panel_mol {} {
    if {[catch { ::ufv::resolve "" } m]} { return "" }
    return $m
}
proc ::ufv::traj::_panel_loaded {args} { _panel_refresh }
proc ::ufv::traj::_panel_refresh {} {
    variable pw
    variable P
    if {![_panel_ok]} return
    set m [_panel_mol]
    if {$m eq ""} { set P(nf) "no molecule"; $pw.run.go state disabled; return }
    set n [nframes $m]
    set P(nf) "of $n (mol $m)"
    if {$n < 2} { $pw.run.go state disabled; set P(nf) "mol $m: 1 frame - load a trajectory (mol addfile)" } \
    else { $pw.run.go state !disabled }
}
proc ::ufv::traj::_panel_kind {} {
    variable pw
    variable P
    variable LABELS
    if {![_panel_ok]} return
    set lab [$pw.run.kind get]
    foreach {l k hint} $LABELS { if {$l eq $lab} { set P(kind) $k; set P(hint) $hint } }
    if {$P(kind) in {dist pockets} && $P(param) eq ""} { set P(param) auto }
    if {$P(kind) ni {dist pockets nearby}} { set P(param) ""; $pw.run.param state disabled } else { $pw.run.param state !disabled }
}
proc ::ufv::traj::_panel_run {} {
    variable pw
    variable P
    if {![_panel_ok]} return
    if {$P(job) ne ""} { cancel $P(job); return }
    set m [_panel_mol]
    if {$m eq ""} return
    set P(mol) $m
    $pw.sum configure -text "running $P(kind) ..."
    $pw.run.go configure -text Cancel
    set P(job) [run_async $m $P(kind) $P(param) {} ::ufv::traj::_panel_done]
}
proc ::ufv::traj::_panel_done {ok r} {
    variable pw
    variable P
    set P(job) ""
    if {![_panel_ok]} return
    $pw.run.go configure -text Run
    if {!$ok} { $pw.sum configure -text "error: $r"; return }
    _panel_show $r
}
proc ::ufv::traj::_panel_result {mol kind} {
    variable P
    if {![_panel_ok] || $P(job) ne ""} return
    set r [result $mol $kind]
    if {$r ne ""} { _panel_show $r }
}
proc ::ufv::traj::_panel_show {r} {
    variable pw
    variable P
    set P(data) $r
    set P(mol) [dict get $r mol]
    $pw.sum configure -text "[dict get $r summary]  \[[dict get $r engine], [expr {[dict exists $r wall_ms] ? [dict get $r wall_ms] : [dict get $r calc_ms]}] ms\]"
    set names [lmap s [dict get $r series] { lindex $s 0 }]
    $pw.ser.c configure -values $names
    $pw.ser.c set [lindex $names 0]
    set tv $pw.tab.tv
    $tv delete [$tv children {}]
    set hd [dict get $r header]
    set ids {}; set i 0; foreach h $hd { lappend ids c$i; incr i }
    $tv configure -columns $ids
    set i 0
    foreach h $hd {
        $tv heading c$i -text $h
        $tv column c$i -width [expr {$i == 0 ? 110 : ($h in {top {top ligand} {top {other chain}} {top neighbour}} ? 160 : 52)}] -stretch [expr {$i == 0}]
        incr i
    }
    set i 0
    foreach row [dict get $r rows] { $tv insert {} end -id r$i -values $row; incr i }
    _panel_plot
}
proc ::ufv::traj::_panel_plot {} {
    variable pw
    variable P
    if {![_panel_ok] || $P(data) eq ""} return
    set c $pw.plot
    $c delete all
    set name [$pw.ser.c get]
    set ys {}
    foreach s [dict get $P(data) series] { if {[lindex $s 0] eq $name} { set ys [lindex $s 1] } }
    set xs [dict get $P(data) frames]
    if {[llength $ys] < 1} return
    set W [winfo width $c]; set H [winfo height $c]
    if {$W < 50} { set W [$c cget -width] }
    if {$H < 50} { set H [$c cget -height] }
    set l 38; set r 6; set t 6; set b 16
    set x0 [lindex $xs 0]; set x1 [lindex $xs end]; if {$x1 == $x0} { set x1 [expr {$x0 + 1}] }
    set y0 [_min $ys]; set y1 [_max $ys]; if {$y1 == $y0} { set y1 [expr {$y0 + 1}] }
    set P(geo) [list $l $r $t $b $W $H $x0 $x1]
    set pts {}
    foreach x $xs y $ys {
        lappend pts [expr {$l + ($x - $x0) * ($W - $l - $r) / double($x1 - $x0)}] [expr {$H - $b - ($y - $y0) * ($H - $t - $b) / double($y1 - $y0)}]
    }
    $c create rectangle $l $t [expr {$W - $r}] [expr {$H - $b}] -outline #ccc
    if {[llength $pts] >= 4} { $c create line {*}$pts -fill #1f5fa8 -width 1 } else { $c create oval {*}$pts {*}$pts -outline #1f5fa8 }
    $c create text 2 $t -anchor nw -text [format %.2f $y1] -font {TkDefaultFont 7} -fill #555
    $c create text 2 [expr {$H - $b}] -anchor sw -text [format %.2f $y0] -font {TkDefaultFont 7} -fill #555
    $c create text $l [expr {$H - 2}] -anchor sw -text $x0 -font {TkDefaultFont 7} -fill #555
    $c create text [expr {$W - $r}] [expr {$H - 2}] -anchor se -text $x1 -font {TkDefaultFont 7} -fill #555
    $c create line 0 0 0 0 -fill red -tags cursor
    _panel_cursor
}
proc ::ufv::traj::_panel_cursor {} {
    variable pw
    variable P
    if {![_panel_ok] || $P(data) eq "" || ![info exists P(geo)]} return
    set m $P(mol)
    if {$m ni [molinfo list]} return
    lassign $P(geo) l r t b W H x0 x1
    set f [molinfo $m get frame]
    set x [expr {$l + ($f - $x0) * ($W - $l - $r) / double($x1 - $x0)}]
    $pw.plot coords cursor $x $t $x [expr {$H - $b}]
}
proc ::ufv::traj::_panel_click {px} {
    variable P
    if {$P(data) eq "" || ![info exists P(geo)]} return
    lassign $P(geo) l r t b W H x0 x1
    set f [expr {$x0 + ($px - $l) * double($x1 - $x0) / ($W - $l - $r)}]
    # nearest analysed frame
    set best ""; set bd 1e9
    foreach x [dict get $P(data) frames] { if {abs($x - $f) < $bd} { set bd [expr {abs($x - $f)}]; set best $x } }
    if {$best ne "" && $P(mol) in [molinfo list]} { animate goto $best; _panel_cursor }
}
proc ::ufv::traj::_panel_row {} {
    variable pw
    variable P
    set sel [$pw.tab.tv selection]
    if {$sel eq "" || $P(data) eq ""} return
    set i [string range [lindex $sel 0] 1 end]
    row_action $P(data) $i
    # plot the matching series
    set kind [dict get $P(data) kind]
    set names [lmap s [dict get $P(data) series] { lindex $s 0 }]
    set want ""
    if {$kind eq "dist"} { set want [lindex $names $i] }
    if {$kind eq "pockets"} { set want "[dict get [lindex [dict get $P(data) rowmeta] $i] pocket] size" }
    if {$want ne "" && $want in $names} { $pw.ser.c set $want; _panel_plot }
}
# Row action (also usable without Tk): dist -> draw that pair; pockets / contacts / groups -> highlight residues.
proc ::ufv::traj::row_action {r i} {
    set mol [dict get $r mol]
    set kind [dict get $r kind]
    set meta [lindex [dict get $r rowmeta] $i]
    switch -- $kind {
        dist { return [show_pairs $mol [list $i] $r] }
        pockets {
            set keys [lmap c [dict get $meta cols] { lindex [dict get $r res] $c }]
            return [_highlight $mol [lmap k $keys { list [lindex $k 0] [lindex $k 2] }]]
        }
        rmsf {
            set ctx [context $mol 0]
            set g [dict get $meta group]
            if {![dict exists $ctx groups $g]} { return 0 }
            set s [::ufv::map::sel $mol [dict get $ctx groups $g]]
            if {$s eq ""} { return 0 }
            ::ufv::layers::define $mol traj:hl [list [list [::ufv::palette::id #e65100] "same residue as ($s)"]] -style {Licorice 0.3 12 12} -kind other
            return 1
        }
        default {
            set k [dict get $meta key]
            set n [_highlight $mol [list $k]]
            ::ufv::zoom_to $mol "chain [::ufv::map::chain_token [lindex $k 0]] and resid [::ufv::map::compress [list [lindex $k 1]]]"
            return $n
        }
    }
}
proc ::ufv::traj::_highlight {mol keys} {
    set cl {}
    foreach k $keys { lappend cl "(chain [::ufv::map::chain_token [lindex $k 0]] and resid [::ufv::map::compress [list [lindex $k 1]]])" }
    if {![llength $cl]} { return 0 }
    ::ufv::layers::define $mol traj:hl [list [list [::ufv::palette::id #e65100] "$::ufv::traj::POLY and ([join $cl { or }])"]] -style {Licorice 0.3 12 12} -kind other
    return [llength $keys]
}
proc ::ufv::traj::_panel_csv {} {
    variable P
    if {$P(mol) eq ""} return
    set f [tk_getSaveFile -title "Trajectory CSV prefix" -initialfile "mol$P(mol)_traj"]
    if {$f eq ""} return
    if {[catch { csv $P(mol) [file rootname $f] } e]} { tk_messageBox -message $e }
}
proc ::ufv::traj::_panel_bpdb {} {
    variable P
    if {$P(mol) eq ""} return
    set mode rmsf
    if {$P(data) ne "" && [dict get $P(data) kind] in {ligand partner}} { set mode [dict get $P(data) kind] }
    set f [tk_getSaveFile -title "B-factor PDB ($mode)" -defaultextension .pdb -initialfile "mol$P(mol)_$mode.pdb"]
    if {$f eq ""} return
    if {[catch { bfactor $P(mol) $f $mode } e]} { tk_messageBox -message $e }
}
