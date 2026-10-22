# =============================================================================================
# ufv_layers.tcl - per-molecule state, UniProt <-> resid numbering, and the layer (representation) model.
# OWNER: foundation. Feature / trajectory / GUI code builds on the APIs below (see API.md).
#
#   ::ufv::st::*      per-molecule state dict (uid, map, layers, cartoon mode, caches)
#   ::ufv::map::*     numbering: identity | sifts | ranges | manual | custom  ->  per-chain linear segments,
#                     O(1) dict lookups UniProt->resid and resid->UniProt, selection strings
#   ::ufv::layers::*  a layer = a named set of reps (one per colour id) found by their VMD repname
#                     (stable when other reps are deleted); show/hide = `mol showrep`; rebuilt in place
#   ::ufv::cartoon::* ONE base cartoon rep per molecule; exclusive cartoon colouring modes:
#                     categorical (per-colour cartoon reps + base excludes them) or continuous
#                     (base rep `mol modcolor` Beta/Occupancy/User.. + `mol scaleminmax`)
# =============================================================================================

namespace eval ::ufv {
    if {![info exists molid]} { variable molid "" }       ;# current molecule ("" = top)
    if {![info exists uid]}   { variable uid "" }         ;# current accession
    variable ann                                          ;# uid -> annotation dict (backend 'annotations')
    variable structs                                      ;# uid -> structure list (backend 'structures')
    variable M                                            ;# molid -> state dict
    variable U2R                                          ;# "molid,chain" -> dict uniprot -> resid
    variable R2U                                          ;# "molid,chain" -> dict resid -> uniprot
    variable hooks                                        ;# event -> list of command prefixes
    if {![array exists ann]}     { array set ann {} }
    if {![array exists structs]} { array set structs {} }
    if {![array exists M]}       { array set M {} }
    if {![array exists U2R]}     { array set U2R {} }
    if {![array exists R2U]}     { array set R2U {} }
    if {![array exists hooks]}   { array set hooks {} }
    if {![info exists nearby]}   { variable nearby 5.0 }
}

# ---------------------------------------------------------------------------------------------
# Status / hooks (Tk-free: the GUI subscribes with ::ufv::hook add status <cmd>)
# ---------------------------------------------------------------------------------------------
# Events: status {msg}, fetched {uid}, loaded {molid key}, mapped {molid}, layers {molid}, cartoon {molid mode},
#         report {pos text}, tiers {molid}. Handlers are called as {*}$cmd {*}$args; errors are reported, not raised.
proc ::ufv::hook {op event {cmd ""}} {
    variable hooks
    switch -- $op {
        add    { if {![info exists hooks($event)] || $cmd ni $hooks($event)} { lappend hooks($event) $cmd } }
        remove { if {[info exists hooks($event)]} { set hooks($event) [lsearch -all -inline -not -exact $hooks($event) $cmd] } }
        list   { return [expr {[info exists hooks($event)] ? $hooks($event) : {}}] }
        default { error "::ufv::hook add|remove|list event ?cmd?" }
    }
}

proc ::ufv::emit {event args} {
    variable hooks
    if {![info exists hooks($event)]} return
    foreach cmd $hooks($event) {
        if {[catch { uplevel #0 [list {*}$cmd {*}$args] } err]} { puts "\[UFV\] hook $event: $err" }
    }
}

proc ::ufv::status {msg} {
    puts "\[UFV\] $msg"
    ::ufv::emit status $msg
}

# Molecule argument: "" / top -> the current UFV molecule if it still exists, else VMD's top molecule.
proc ::ufv::resolve {{molid ""}} {
    if {$molid eq "" || $molid eq "top"} {
        if {$::ufv::molid ne "" && $::ufv::molid in [molinfo list]} { return $::ufv::molid }
        set t [molinfo top]
        if {$t < 0} { error "no molecule loaded" }
        return $t
    }
    if {$molid ni [molinfo list]} { error "no molecule $molid" }
    return $molid
}

# atomselect with guaranteed cleanup: ::ufv::with_sel s $mol "protein" { $s get resid }
proc ::ufv::with_sel {var mol text body} {
    upvar 1 $var s
    set s [atomselect $mol $text]
    set code [catch { uplevel 1 $body } res]
    catch { $s delete }
    if {$code == 1} { return -code error -errorinfo $::errorInfo $res }
    return $res
}

# Annotation dict of uid (fetched once; sync). ::ufv::fetch_async fills it without blocking.
proc ::ufv::annotations {{uid ""}} {
    variable ann
    if {$uid eq ""} { set uid $::ufv::uid }
    if {$uid eq ""} { error "no accession fetched yet - run: ufv_fetch <uniprot_id>" }
    set uid [string toupper $uid]
    if {![info exists ann($uid)]} {
        set ann($uid) [::ufv::backend::call annotations [dict create uid $uid]]
    }
    return $ann($uid)
}

# ---------------------------------------------------------------------------------------------
# Per-molecule state
# ---------------------------------------------------------------------------------------------
namespace eval ::ufv::st {}

proc ::ufv::st::init {mol} {
    if {[info exists ::ufv::M($mol)]} return
    set ::ufv::M($mol) [dict create uid "" map [dict create mode identity segments {} chains {}] mapver 0 \
                            layers [dict create] cartoon [dict create mode none] chains {} cache [dict create]]
}
proc ::ufv::st::has {mol} { return [info exists ::ufv::M($mol)] }
proc ::ufv::st::get {mol key {default ""}} {
    init $mol
    if {[dict exists $::ufv::M($mol) {*}$key]} { return [dict get $::ufv::M($mol) {*}$key] }
    return $default
}
proc ::ufv::st::put {mol key value} {
    init $mol
    dict set ::ufv::M($mol) {*}$key $value
}
proc ::ufv::st::remove {mol key} {
    init $mol
    if {[dict exists $::ufv::M($mol) {*}$key]} { dict unset ::ufv::M($mol) {*}$key }
}
proc ::ufv::st::drop {mol} {
    unset -nocomplain ::ufv::M($mol)
    array unset ::ufv::U2R "$mol,*"
    array unset ::ufv::R2U "$mol,*"
    if {$::ufv::molid eq $mol} { set ::ufv::molid "" }
}
# Cache slot invalidated by every numbering change: ::ufv::st::cached $mol key { script }
proc ::ufv::st::cached {mol key script} {
    set ver [get $mol mapver 0]
    set c [get $mol [list cache $key] ""]
    if {$c ne "" && [lindex $c 0] == $ver} { return [lindex $c 1] }
    set v [uplevel 1 $script]
    put $mol [list cache $key] [list $ver $v]
    return $v
}

# Forget the state of deleted molecules (VMD sets vmd_initialize_structure(molid) 0 on deletion).
proc ::ufv::st::_moltrace {name1 name2 op} {
    if {$name2 ne "" && [info exists ::vmd_initialize_structure($name2)] && !$::vmd_initialize_structure($name2)} {
        ::ufv::palette::scale_release $name2    ;# a deleted molecule's continuous colouring gives the scale back
        catch { foreach f [::ufv::st::get $name2 {cache coordfiles} {}] { file delete $f } }   ;# Calpha JSON temp files
        ::ufv::st::drop $name2
    }
}
if {[llength [info commands molinfo]] && ![info exists ::ufv::_moltrace_on]} {
    set ::ufv::_moltrace_on 1
    trace add variable ::vmd_initialize_structure write ::ufv::st::_moltrace
}

# ---------------------------------------------------------------------------------------------
# Numbering
# ---------------------------------------------------------------------------------------------
# A map dict:  mode identity|sifts|ranges|manual|custom
#              segments {chain {{u_start u_end resid_start resid_end} ...}}   (all modes except identity)
#              chains   {chain ...}   identity: the chains numbered by UniProt position (default all protein)
#              + status, notes, meta, pdbId, config/accessions/partners (custom), base/basemode (custom 'auto')
namespace eval ::ufv::map {}

# Protein chains of a molecule (cached per molecule; structure changes do not happen after load).
proc ::ufv::map::protein_chains {mol} {
    set c [::ufv::st::get $mol chains]
    if {$c ne ""} { return $c }
    set c [::ufv::with_sel s $mol "protein and name CA" { lsort -unique [$s get chain] }]
    if {![llength $c]} { set c [::ufv::with_sel s $mol "protein" { lsort -unique [$s get chain] }] }
    ::ufv::st::put $mol chains $c
    return $c
}

# {chain {{resid one-letter} ...}} of the protein CA atoms (insertion-code residues skipped): for the backend.
proc ::ufv::map::residues {mol} {
    variable aa3
    set out [dict create]
    ::ufv::with_sel s $mol "protein and name CA" {
        set seen [dict create]
        foreach {ch ri rn ins} [join [$s get {chain resid resname insertion}]] {
            if {[string trim $ins] ne ""} continue
            if {[dict exists $seen $ch,$ri]} continue
            dict set seen $ch,$ri 1
            set a [expr {[dict exists $aa3 $rn] ? [dict get $aa3 $rn] : "X"}]
            dict lappend out $ch [list $ri $a]
        }
    }
    return $out
}
set ::ufv::map::aa3 [dict create ALA A ARG R ASN N ASP D CYS C GLN Q GLU E GLY G HIS H ILE I LEU L LYS K \
    MET M PHE F PRO P SER S THR T TRP W TYR Y VAL V MSE M SEC U PYL O HID H HIE H HIP H HSD H HSE H HSP H CYX C]

# Install a map on a molecule: clears the lookup caches, bumps mapver and rebuilds the visible layers.
proc ::ufv::map::install {mol m} {
    if {![dict exists $m mode]} { error "map without mode" }
    if {![dict exists $m segments]} { dict set m segments {} }
    if {[dict get $m mode] eq "identity" && (![dict exists $m chains] || [dict get $m chains] eq "")} {
        dict set m chains [protein_chains $mol]
    }
    ::ufv::st::put $mol map $m
    ::ufv::st::put $mol mapver [expr {[::ufv::st::get $mol mapver 0] + 1}]
    array unset ::ufv::U2R "$mol,*"
    array unset ::ufv::R2U "$mol,*"
    ::ufv::layers::remapped $mol
    ::ufv::emit mapped $mol
}

proc ::ufv::map::get {mol} {
    set m [::ufv::st::get $mol map]
    if {[dict get $m mode] eq "identity" && (![dict exists $m chains] || ![llength [dict get $m chains]])} {
        dict set m chains [protein_chains $mol]       ;# default map of a fresh molecule: identity, all chains
        ::ufv::st::put $mol map $m
    }
    return $m
}
proc ::ufv::map::mode {mol} { return [::ufv::st::get $mol {map mode} identity] }

# Chains that carry UniProt positions.
proc ::ufv::map::chains {mol} {
    set m [get $mol]
    if {[dict get $m mode] eq "identity"} { return [dict get $m chains] }
    return [dict keys [dict get $m segments]]
}

proc ::ufv::map::_build {mol ch} {
    set m [get $mol]
    set u2r [dict create]; set r2u [dict create]
    if {[dict exists $m segments $ch]} {
        foreach seg [dict get $m segments $ch] {
            lassign $seg us ue rs re
            set off [expr {$rs - $us}]
            for {set u $us} {$u <= $ue} {incr u} {
                set r [expr {$u + $off}]
                if {![dict exists $u2r $u]} { dict set u2r $u $r }
                if {![dict exists $r2u $r]} { dict set r2u $r $u }
            }
        }
    }
    set ::ufv::U2R($mol,$ch) $u2r
    set ::ufv::R2U($mol,$ch) $r2u
}

proc ::ufv::map::_identity {mol ch} {
    set m [get $mol]
    return [expr {[dict get $m mode] eq "identity" && $ch in [dict get $m chains]}]
}

# UniProt position -> resid on chain ("" when unmapped).
proc ::ufv::map::u2r {mol ch pos} {
    if {[_identity $mol $ch]} { return $pos }
    if {![info exists ::ufv::U2R($mol,$ch)]} { _build $mol $ch }
    if {[dict exists $::ufv::U2R($mol,$ch) $pos]} { return [dict get $::ufv::U2R($mol,$ch) $pos] }
    return ""
}

# resid -> UniProt position ("" when unmapped).
proc ::ufv::map::r2u {mol ch resid} {
    if {[_identity $mol $ch]} { return $resid }
    if {![info exists ::ufv::R2U($mol,$ch)]} { _build $mol $ch }
    if {[dict exists $::ufv::R2U($mol,$ch) $resid]} { return [dict get $::ufv::R2U($mol,$ch) $resid] }
    return ""
}

# Bulk: positions -> resids on one chain (unmapped dropped).
proc ::ufv::map::resids {mol ch positions} {
    if {[_identity $mol $ch]} { return $positions }
    if {![info exists ::ufv::U2R($mol,$ch)]} { _build $mol $ch }
    set d $::ufv::U2R($mol,$ch)
    set out {}
    foreach p $positions { if {[dict exists $d $p]} { lappend out [dict get $d $p] } }
    return $out
}

# Sorted integer list -> VMD resid tokens with ranges ("1 to 5 9 12 to 14"); negatives quoted.
proc ::ufv::map::compress {ints} {
    set ints [lsort -integer -unique $ints]
    set out {}
    set n [llength $ints]
    set i 0
    while {$i < $n} {
        set a [lindex $ints $i]; set j $i
        while {$j + 1 < $n && [lindex $ints [expr {$j+1}]] == [lindex $ints $j] + 1} { incr j }
        set b [lindex $ints $j]
        set qa [expr {$a < 0 ? "\"$a\"" : $a}]
        set qb [expr {$b < 0 ? "\"$b\"" : $b}]
        if {$b - $a >= 2} { lappend out "$qa to $qb" } elseif {$b > $a} { lappend out $qa $qb } else { lappend out $qa }
        set i [expr {$j + 1}]
    }
    return [join $out " "]
}

proc ::ufv::map::chain_token {ch} {
    if {[regexp {^[A-Za-z0-9]+$} $ch]} { return $ch }
    return "'$ch'"
}

# Selection for UniProt positions on every mapped chain ("" if nothing maps).
#   -ca 1        restrict to CA (markers)          -within <text>  extra restriction (default "protein")
proc ::ufv::map::sel {mol positions args} {
    array set o {-ca 0 -within protein}
    array set o $args
    if {![llength $positions]} { return "" }
    set m [get $mol]
    set clauses {}
    set chs [chains $mol]
    if {[dict get $m mode] eq "identity" && [lsort $chs] eq [lsort [protein_chains $mol]]} {
        lappend clauses "resid [compress $positions]"
    } else {
        foreach ch $chs {
            set r [resids $mol $ch $positions]
            if {[llength $r]} { lappend clauses "(chain [chain_token $ch] and resid [compress $r])" }
        }
    }
    if {![llength $clauses]} { return "" }
    set s "([join $clauses { or }])"
    if {$o(-within) ne ""} { set s "($o(-within)) and $s" }
    if {$o(-ca)} { append s " and name CA" }
    return $s
}

# Per-atom UniProt position list for `atomselect $mol $text` (cached per mapver): for field writes.
proc ::ufv::map::atom_uni {mol {text protein}} {
    return [::ufv::st::cached $mol "atom_uni:$text" {
        ::ufv::with_sel s $mol $text {
            set out {}
            foreach {ch ri} [join [$s get {chain resid}]] { lappend out [::ufv::map::r2u $mol $ch $ri] }
            set out
        }
    }]
}

# Manual anchor -> segments: resid rstart == UniProt ustart, up to resid rend; open-ended (no rend) runs to the
# end of the UniProt sequence (or at least the chain's last residue), like the PyMOL plugin's resi_end None.
proc ::ufv::map::manual_segments {mol ch ustart rstart {rend ""}} {
    if {$rend eq ""} {
        set rend [::ufv::with_sel s $mol "protein and chain [chain_token $ch]" {
            set r [$s get resid]
            expr {[llength $r] ? [tcl::mathfunc::max {*}$r] : $rstart}
        }]
        set u [::ufv::st::get $mol uid]
        if {$u eq ""} { set u $::ufv::uid }
        if {$u ne "" && [info exists ::ufv::ann($u)]} {
            set rend [expr {max($rend, $rstart + [dict get $::ufv::ann($u) seqlen] - $ustart)}]
        }
    }
    if {$rend < $rstart} { return {} }
    return [list [list $ustart [expr {$ustart + $rend - $rstart}] $rstart $rend]]
}

# Human-readable one-line summary of a molecule's numbering.
proc ::ufv::map::describe {mol} {
    set m [get $mol]
    set mode [dict get $m mode]
    if {$mode eq "identity"} { return "identity (resid == UniProt) on chain(s) [join [dict get $m chains] {, }]" }
    set parts {}
    dict for {ch segs} [dict get $m segments] {
        set t {}
        foreach s $segs { lassign $s us ue rs re; lappend t "$us-$ue<->$rs-$re" }
        lappend parts "$ch: [join $t {; }]"
    }
    return "$mode: [join $parts { | }]"
}

# ---------------------------------------------------------------------------------------------
# Layers
# ---------------------------------------------------------------------------------------------
# Layer dict (st key {layers <name>}): reps {repname ...} visible 0|1 kind marker|cartoon|base|other
#   rebuild <script or ""> (re-run on numbering change)  groups <the last group spec>
# Group spec: list of {colour selection ?style?} where colour is a colour id (integer) or a VMD colouring
#   method string (e.g. Name, Beta, "User2"), selection is VMD selection text, style overrides -style.
namespace eval ::ufv::layers {
    variable batchdepth 0
}

# Run script with the display frozen (one redraw at the end; nested calls are fine).
proc ::ufv::layers::batch {script} {
    variable batchdepth
    set has [llength [info commands display]]
    if {$has && $batchdepth == 0} { catch { display update off } }
    incr batchdepth
    set code [catch { uplevel 1 $script } res]
    incr batchdepth -1
    if {$has && $batchdepth == 0} { catch { display update on } }
    if {$code == 1} { return -code error -errorinfo $::errorInfo $res }
    return $res
}

# Current rep indices of a layer (reps deleted by the user are dropped from the layer).
proc ::ufv::layers::reps {mol name} {
    set names [::ufv::st::get $mol [list layers $name reps] {}]
    set idx {}; set keep {}
    foreach rn $names {
        set i [mol repindex $mol $rn]
        if {$i >= 0} { lappend idx $i; lappend keep $rn }
    }
    if {[llength $keep] != [llength $names]} { ::ufv::st::put $mol [list layers $name reps] $keep }
    return $idx
}

proc ::ufv::layers::exists {mol name} { return [expr {[llength [reps $mol $name]] > 0}] }
proc ::ufv::layers::names {mol} { return [dict keys [::ufv::st::get $mol layers [dict create]]] }
proc ::ufv::layers::visible {mol name} { return [::ufv::st::get $mol [list layers $name visible] 0] }

proc ::ufv::layers::_color {mol idx colour} {
    if {[string is integer -strict $colour]} { mol modcolor $idx $mol ColorID $colour } else { mol modcolor $idx $mol {*}$colour }
}

# Create or update a layer in place. groups = list of {colour selection}; options:
#   -style {VDW 0.6 12}  -material Opaque  -kind marker  -visible 1  -rebuild <script>
# Existing reps are reused (modselect / modcolor / modstyle) - no rep churn when a layer is rebuilt.
proc ::ufv::layers::define {mol name groups args} {
    array set o {-style {VDW 0.6 12} -material Opaque -kind marker -visible 1 -rebuild ""}
    array set o $args
    set groups [lsearch -all -inline -not -index 1 $groups ""]
    batch {
        set old [reps $mol $name]
        set names [::ufv::st::get $mol [list layers $name reps] {}]
        set keep {}
        set i 0
        foreach g $groups {
            lassign $g colour seltext gstyle
            set style [expr {$gstyle ne "" ? $gstyle : $o(-style)}]
            if {$i < [llength $old]} {
                set idx [lindex $old $i]
                mol modstyle $idx $mol {*}$style
                mol modselect $idx $mol $seltext
                _color $mol $idx $colour
                mol modmaterial $idx $mol $o(-material)
                lappend keep [lindex $names $i]
            } else {
                mol representation {*}$style
                mol selection $seltext
                if {[string is integer -strict $colour]} { mol color ColorID $colour } else { mol color {*}$colour }
                mol material $o(-material)
                mol addrep $mol
                set idx [expr {[molinfo $mol get numreps] - 1}]
                lappend keep [mol repname $mol $idx]
            }
            mol showrep $mol $idx $o(-visible)
            incr i
        }
        # surplus reps (fewer groups than before): delete highest index first
        foreach idx [lsort -integer -decreasing [lrange $old $i end]] { mol delrep $idx $mol }
        ::ufv::st::put $mol [list layers $name] [dict create reps $keep visible $o(-visible) kind $o(-kind) \
            rebuild $o(-rebuild) groups $groups style $o(-style)]
    }
    ::ufv::emit layers $mol
    return [llength $keep]
}

proc ::ufv::layers::show {mol name {on 1}} {
    set idx [reps $mol $name]
    if {![llength $idx]} {
        if {[dict exists [::ufv::st::get $mol layers] $name]} { ::ufv::st::put $mol [list layers $name visible] [expr {$on ? 1 : 0}] }
        return 0
    }
    batch { foreach i $idx { mol showrep $mol $i $on } }
    ::ufv::st::put $mol [list layers $name visible] [expr {$on ? 1 : 0}]
    ::ufv::emit layers $mol
    return 1
}

proc ::ufv::layers::toggle {mol name} {
    set v [visible $mol $name]
    show $mol $name [expr {!$v}]
    return [expr {!$v}]
}

proc ::ufv::layers::delete {mol name} {
    batch {
        foreach i [lsort -integer -decreasing [reps $mol $name]] { mol delrep $i $mol }
    }
    ::ufv::st::remove $mol [list layers $name]
    ::ufv::emit layers $mol
}

# Delete every UFV layer of a molecule (user reps are never touched); the base cartoon goes back to plain.
proc ::ufv::layers::clear {mol} {
    batch {
        foreach n [names $mol] { if {$n ne "cartoon"} { delete $mol $n } }
        if {[exists $mol cartoon]} { ::ufv::cartoon::plain $mol }
    }
}

# Numbering changed: rebuild layers that know how (their -rebuild script), drop the others.
proc ::ufv::layers::remapped {mol} {
    batch {
        foreach n [names $mol] {
            if {$n eq "cartoon"} continue
            set script [::ufv::st::get $mol [list layers $n rebuild] ""]
            set vis [visible $mol $n]
            if {$script eq ""} {
                if {[string match cart:* $n]} { delete $mol $n }
                continue
            }
            if {[catch { uplevel #0 $script } err]} {
                puts "\[UFV\] layer $n not rebuilt after renumbering: $err"
                delete $mol $n
            } elseif {!$vis} {
                show $mol $n 0
            }
        }
        set mode [::ufv::cartoon::mode $mol]
        if {$mode ne "none" && ![exists $mol "cart:$mode"] && [::ufv::st::get $mol {cartoon kind}] eq "categorical"} {
            ::ufv::cartoon::plain $mol
        }
    }
}

# Marker layer from {hex {positions}} groups (UniProt positions; CA markers).
# Returns the number of distinct positions that mapped onto the structure.
proc ::ufv::layers::markers {mol name hexgroups args} {
    array set o [list -style {VDW 0.6 12} -rebuild "" -visible 1]
    array set o $args
    set groups {}
    set shown 0
    dict for {cid pos} [::ufv::palette::group $hexgroups] {
        set s [::ufv::map::sel $mol $pos -ca 1]
        if {$s eq ""} continue
        lappend groups [list $cid $s]
    }
    if {![llength $groups]} {
        if {[exists $mol $name]} { delete $mol $name }
        return 0
    }
    set all {}
    dict for {hex pos} $hexgroups { lappend all {*}$pos }
    foreach ch [::ufv::map::chains $mol] { incr shown [llength [::ufv::map::resids $mol $ch [lsort -integer -unique $all]]] }
    define $mol $name $groups -style $o(-style) -kind marker -rebuild $o(-rebuild) -visible $o(-visible)
    return $shown
}

# Write per-residue values (dict UniProt position -> value) into an atom field (user, user2-4, beta,
# occupancy) of every protein atom; unmapped / missing residues get $default. One `$sel set` call.
proc ::ufv::layers::write_field {mol field values {default 0.0} {text protein}} {
    set uni [::ufv::map::atom_uni $mol $text]
    set out {}
    foreach u $uni {
        if {$u ne "" && [dict exists $values $u]} { lappend out [dict get $values $u] } else { lappend out $default }
    }
    ::ufv::with_sel s $mol $text { $s set $field $out }
    return [llength $out]
}

# ---------------------------------------------------------------------------------------------
# Cartoon: one base rep ("cartoon" layer) + exclusive colouring modes
# ---------------------------------------------------------------------------------------------
namespace eval ::ufv::cartoon {
    variable fields {Beta Occupancy User User2 User3 User4}
    # Every UFV cartoon rep. VMD's default NewCartoon is 0.3 12 4.5 0; resolution 8 looks the same at normal
    # zoom with 53% less geometry (Wavefront export 7.8 -> 3.7 MB for 300 residues). Static scenes are
    # vsync-bound either way; on a trajectory the per-frame rebuild drops from ~14-16 to ~9-10 ms/frame
    # (P35498 AF, 500 frames, 28 reps).
    variable style {NewCartoon 0.3 8 4.5 0}
}

# Make sure the base cartoon exists. replace_default=1 deletes VMD's default rep 0 (fresh loads only).
proc ::ufv::cartoon::ensure {mol {replace_default 0}} {
    if {[::ufv::layers::exists $mol cartoon]} return
    ::ufv::layers::batch {
        if {$replace_default && [molinfo $mol get numreps] == 1} { mol delrep 0 $mol }
        ::ufv::layers::define $mol cartoon [list [list [::ufv::palette::base] protein]] \
            -style $::ufv::cartoon::style -kind base -material Opaque
    }
    ::ufv::st::put $mol cartoon [dict create mode none kind none]
}

proc ::ufv::cartoon::base_idx {mol} {
    ensure $mol
    return [lindex [::ufv::layers::reps $mol cartoon] 0]
}

proc ::ufv::cartoon::mode {mol} { return [::ufv::st::get $mol {cartoon mode} none] }

proc ::ufv::cartoon::_hide_modes {mol {except ""}} {
    foreach n [::ufv::layers::names $mol] {
        if {[string match cart:* $n] && $n ne "cart:$except"} { ::ufv::layers::show $mol $n 0 }
    }
}

# Plain cartoon (no colouring mode).
proc ::ufv::cartoon::plain {mol} {
    ::ufv::palette::scale_release $mol
    ::ufv::layers::batch {
        set b [base_idx $mol]
        _hide_modes $mol
        mol modselect $b $mol protein
        mol modcolor $b $mol ColorID [::ufv::palette::base]
        mol showrep $mol $b 1
    }
    ::ufv::st::put $mol cartoon [dict create mode none kind none]
    ::ufv::emit cartoon $mol none
}

# Categorical mode: hexgroups {hex {UniProt positions}} -> per-colour cartoon reps (layer cart:<mode>);
# the base rep shows the remaining residues in basehex (or the neutral colour).
# Re-selecting an already built mode with the same data only shows it (fast path).
proc ::ufv::cartoon::categorical {mol mode hexgroups {basehex ""} {rebuild ""}} {
    ::ufv::palette::scale_release $mol           ;# categorical modes use colour ids, not the scale
    ::ufv::layers::batch {
        set b [base_idx $mol]
        _hide_modes $mol $mode
        set groups {}
        set covered {}
        dict for {cid pos} [::ufv::palette::group $hexgroups] {
            set s [::ufv::map::sel $mol $pos]
            if {$s eq ""} continue
            lappend groups [list $cid $s]
            lappend covered "($s)"
        }
        ::ufv::layers::define $mol "cart:$mode" $groups -style $::ufv::cartoon::style -kind cartoon -rebuild $rebuild
        set bsel [expr {[llength $covered] ? "protein and not ([join $covered { or }])" : "protein"}]
        mol modselect $b $mol $bsel
        mol modcolor $b $mol ColorID [expr {$basehex ne "" ? [::ufv::palette::id $basehex] : [::ufv::palette::base]}]
        mol showrep $mol $b 1
    }
    ::ufv::st::put $mol cartoon [dict create mode $mode kind categorical basesel $bsel basehex $basehex]
    ::ufv::st::put $mol [list modes $mode] [dict create basesel $bsel basehex $basehex]
    ::ufv::emit cartoon $mol $mode
    return [llength $groups]
}

# Re-show a categorical mode that was built before (no rebuild). Returns 0 if it has to be built.
proc ::ufv::cartoon::reshow {mol mode} {
    if {![::ufv::layers::exists $mol "cart:$mode"]} { return 0 }
    set last [::ufv::st::get $mol [list modes $mode] ""]
    if {$last eq ""} { return 0 }
    ::ufv::palette::scale_release $mol
    ::ufv::layers::batch {
        set b [base_idx $mol]
        _hide_modes $mol $mode
        ::ufv::layers::show $mol "cart:$mode" 1
        mol modselect $b $mol [dict get $last basesel]
        mol modcolor $b $mol ColorID [expr {[dict get $last basehex] ne "" ? [::ufv::palette::id [dict get $last basehex]] : [::ufv::palette::base]}]
    }
    ::ufv::st::put $mol cartoon [dict merge $last [dict create mode $mode kind categorical]]
    ::ufv::emit cartoon $mol $mode
    return 1
}

# Continuous mode: colour the base cartoon by an atom field with the colour scale.
#   field: Beta | Occupancy | User | User2 | User3 | User4     min/max "" = auto (data range)
#   method: colour-scale method to borrow while this colouring is active (e.g. BWR; the scale is GLOBAL in VMD,
#   see ::ufv::palette::scale_acquire); "" = keep the user's method.
proc ::ufv::cartoon::continuous {mol mode field {min ""} {max ""} {method ""}} {
    variable fields
    if {$field ni $fields} { error "field must be one of $fields" }
    if {$method eq ""} { ::ufv::palette::scale_release $mol } else { ::ufv::palette::scale_acquire $mol $method }
    ::ufv::layers::batch {
        set b [base_idx $mol]
        _hide_modes $mol
        mol modselect $b $mol protein
        mol modcolor $b $mol $field
        if {$min ne "" && $max ne ""} { mol scaleminmax $mol $b $min $max } else { mol scaleminmax $mol $b auto }
        mol showrep $mol $b 1
    }
    ::ufv::st::put $mol cartoon [dict create mode $mode kind continuous field $field min $min max $max]
    ::ufv::emit cartoon $mol $mode
}

# ---------------------------------------------------------------------------------------------
# View (camera) stability. The view only moves on EXPLICIT view actions: focus / zoom to a residue, range, ligand
# or pocket, Reset view, and VMD's own fit when the very first molecule is loaded. Everything else (layers,
# colourings, numbering, loads of further molecules, trajectory analyses, ...) leaves every molecule's
# center / rotate / scale / global matrix exactly as it was.
#   Measured (VMD 1.9.4): of the operations the plugin uses, only `mol new <file>` changes the view - it re-fits
#   EVERY molecule when a molecule after the first is loaded (also for hidden scratch molecules). Reps, colours,
#   mol top / on / off / delete / addfile, animate, ssrecalc, labels and graphics never touch the matrices.
#   So `mol new` runs inside ::ufv::view::keep, which puts every matrix back exactly and gives the new molecule
#   the view the others share.
#   Focus zoom: deterministic (selection extent + buffer, the current orientation kept, like PyMOL's zoom); the
#   view before the first focus of a chain is remembered and Reset view returns to it exactly when the camera was
#   not moved since (PyMOL plugin _focus_zoom / ufv_resetview).
# ---------------------------------------------------------------------------------------------
namespace eval ::ufv::view {
    variable depth 0
    variable saved ""
    variable savedtop -1
    variable pre ""
    variable post ""
    variable keys {center_matrix rotate_matrix scale_matrix global_matrix}
}
proc ::ufv::view::ok {} { return [expr {[llength [info commands ::molinfo]] > 0}] }

# Every molecule's view: dict molid -> {center rotate scale global}
proc ::ufv::view::get {} {
    variable keys
    set v [dict create]
    if {![ok]} { return $v }
    foreach m [molinfo list] { dict set v $m [molinfo $m get $keys] }
    return $v
}

# Put a saved view back (exactly). Molecules that are not in it (loaded since) get the reference molecule's view.
proc ::ufv::view::put {v {ref ""}} {
    variable keys
    if {![ok] || ![dict size $v]} return
    if {$ref eq "" || ![dict exists $v $ref]} { set ref [lindex [dict keys $v] 0] }
    foreach m [molinfo list] {
        set want [expr {[dict exists $v $m] ? [dict get $v $m] : [dict get $v $ref]}]
        if {[molinfo $m get $keys] ne $want} { molinfo $m set $keys $want }
    }
}

# Run script; afterwards every view matrix is what it was before (nesting-safe; synchronous scripts only).
proc ::ufv::view::keep {script} {
    variable depth; variable saved; variable savedtop
    if {$depth > 0 || ![ok]} { return [uplevel 1 $script] }
    set saved [get]
    set savedtop [molinfo top]
    incr depth
    set code [catch { uplevel 1 $script } res]
    set einfo $::errorInfo
    incr depth -1
    set v $saved; set saved ""
    catch { put $v $savedtop }
    if {$code == 1} { return -code error -errorinfo $einfo $res }
    return $res
}

# An explicit view action inside keep: its result is the view to end on.
proc ::ufv::view::commit {} {
    variable depth; variable saved
    if {$depth > 0} { set saved [get] }
}

proc ::ufv::view::_ident {} { return {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} }

# Set the view of every (non-fixed) molecule: centre on point c (molecule coordinates; molecules share one frame,
# as after VMD's own resetview), rotation R (4x4), uniform scale s, no global translation.
proc ::ufv::view::_apply {mol c R s} {
    variable keys
    set C [transoffset [vecscale -1.0 $c]]
    set S [list [list $s 0 0 0] [list 0 $s 0 0] [list 0 0 $s 0] {0 0 0 1}]
    set want [list $C $R $S [_ident]]
    foreach m [molinfo list] {
        if {$m != $mol && ![catch { molinfo $m get fixed } fx] && $fx} continue
        molinfo $m set $keys $want
    }
}

# Fit a selection: bounding-box centre, scale 1.5 / (largest extent + 2 buffer) (VMD's resetview uses 1.5 / extent),
# the target molecule's current orientation. Returns 0 when the selection is empty.
proc ::ufv::view::fit {mol selstr {buffer 3}} {
    if {![ok] || $mol ni [molinfo list]} { return 0 }
    set sel [atomselect $mol $selstr]
    set n [$sel num]
    if {$n} { lassign [measure minmax $sel] mn mx }
    $sel delete
    if {!$n} { return 0 }
    set c [vecscale 0.5 [vecadd $mn $mx]]
    set span 0.0
    foreach d [vecsub $mx $mn] { if {$d > $span} { set span $d } }
    set span [expr {$span + 2.0 * $buffer}]
    if {$span < 1.0} { set span 1.0 }
    _apply $mol $c [lindex [molinfo $mol get rotate_matrix] 0] [expr {1.5 / $span}]
    return 1
}

# Explicit focus zoom (residue / range / ligand / pocket).
proc ::ufv::view::zoom {mol selstr {buffer 3}} {
    variable pre; variable post
    if {![ok]} { return 0 }
    set cur [get]
    if {$post eq "" || $cur ne $post} { set pre $cur }
    if {![fit $mol $selstr $buffer]} { return 0 }
    set post [get]
    commit
    return 1
}

# Explicit Reset view: back to the view before the focus chain when the camera was not moved since, else the
# whole molecule (orientation kept).
proc ::ufv::view::reset {mol} {
    variable pre; variable post
    if {![ok]} return
    set p $pre; set q $post
    set pre ""; set post ""
    if {$p ne "" && $q ne "" && [get] eq $q} {
        put $p
    } else {
        fit $mol all 0
    }
    commit
}
