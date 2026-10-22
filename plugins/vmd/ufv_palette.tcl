# =============================================================================================
# ufv_palette.tcl - fixed colour strategy for the VMD plugin.
# OWNER: foundation. Everyone maps a hex colour to a VMD colour id with ::ufv::palette::id <hex>.
#
# VMD has 33 named colour ids (0-32); ids >= 33 are the colour scale. The plugin:
#   * leaves ids 0-16 untouched (they are still used as nearest-colour candidates, with their CURRENT rgb);
#   * defines ids 17-32 ONCE at load (::ufv::palette::init) with the extension's main categorical colours and
#     never re-assigns them - so a layer built earlier can never be silently recoloured by a later one
#     (the old next_color recycling bug);
#   * any other hex (PTM categories, domains, diseases, topology...) is mapped to the NEAREST id among the
#     16 slots + ids 0-15 (white 8 / black 16 excluded), memoised;
#   * continuous colourings (B-factor, RMSF, occupancy, AM score...) do not use ids: they go through VMD's
#     colour scale (Beta / Occupancy / User..User4 colouring + `mol scaleminmax`), see ufv_layers.tcl.
# Note: the colour-scale METHOD is global in VMD (all molecules and every scale-coloured rep). The plugin only
# borrows it while a UFV continuous colouring that asks for a method is active (scale_acquire / scale_release):
# the user's method is saved on the first borrow and restored when the last such colouring goes away (cleared,
# replaced, ufv_clear / ufv_hide, molecule deleted) - unless the user changed the method meanwhile, which wins.
# =============================================================================================

namespace eval ::ufv::palette {
    # id hex role  (role names are stable API: ::ufv::palette::role <name>)
    variable slots {
        17 #d55e00 pathogenic
        18 #e6a817 deleterious
        19 #0072b2 benign
        20 #9e9e9e uncertain
        21 #fbc02d site
        22 #6d4c41 mutagenesis
        23 #8b5cf6 phosphorylation
        24 #06b6d4 acetylation
        25 #16a34a methylation
        26 #e879f9 ubiquitination
        27 #b71c1c am_pathogenic
        28 #e06666 am_likely_pathogenic
        29 #b9c2cf am_ambiguous
        30 #3d85c8 am_benign
        31 #0053d6 plddt_very_high
        32 #ff7d45 plddt_very_low
    }
    variable base_id 6           ;# neutral cartoon colour (VMD silver, untouched)
    variable exclude {8 16}      ;# never chosen as nearest (white, black)
    variable cache               ;# hex -> id
    variable rgb                 ;# id -> {r g b} 0-255 (candidates)
    variable roles               ;# role -> id
    variable initialised 0
    array set cache {}
    array set rgb {}
    array set roles {}
}

proc ::ufv::palette::hex_rgb {hex} {
    set h [string trimleft [string tolower $hex] "#"]
    if {[string length $h] == 3} {
        set h [string index $h 0][string index $h 0][string index $h 1][string index $h 1][string index $h 2][string index $h 2]
    }
    if {![regexp {^[0-9a-f]{6}$} $h]} { error "bad colour '$hex'" }
    scan $h %2x%2x%2x r g b
    return [list $r $g $b]
}

# Define slots 17-32 (idempotent: call once at load; re-running restores them if a user changed them).
proc ::ufv::palette::init {} {
    variable slots
    variable rgb
    variable roles
    variable cache
    variable exclude
    variable initialised
    array unset cache
    array unset rgb
    set have_vmd [llength [info commands color]]
    foreach {id hex role} $slots {
        lassign [hex_rgb $hex] r g b
        if {$have_vmd} { color change rgb $id [expr {$r/255.0}] [expr {$g/255.0}] [expr {$b/255.0}] }
        set rgb($id) [list $r $g $b]
        set roles($role) $id
        set cache([string tolower $hex]) $id
    }
    # ids 0-15 as candidates, with their current rgb (users may have changed them)
    for {set id 0} {$id <= 15} {incr id} {
        if {$id in $exclude} continue
        if {$have_vmd && [llength [info commands colorinfo]] && ![catch {colorinfo rgb $id} c]} {
            lassign $c r g b
            set rgb($id) [list [expr {round($r*255)}] [expr {round($g*255)}] [expr {round($b*255)}]]
        }
    }
    set initialised 1
}

# hex -> colour id (exact slot, else nearest candidate by the "redmean" perceptual distance). Memoised.
proc ::ufv::palette::id {hex} {
    variable cache
    variable rgb
    variable initialised
    if {!$initialised} { init }
    set k [string tolower $hex]
    if {[info exists cache($k)]} { return $cache($k) }
    lassign [hex_rgb $hex] r g b
    set best -1; set bd 1e30
    foreach cid [array names rgb] {
        lassign $rgb($cid) r2 g2 b2
        set rm [expr {($r + $r2) / 2.0}]
        set dr [expr {$r - $r2}]; set dg [expr {$g - $g2}]; set db [expr {$b - $b2}]
        set d [expr {(2 + $rm/256.0)*$dr*$dr + 4*$dg*$dg + (2 + (255 - $rm)/256.0)*$db*$db}]
        # prefer our slots on ties, then lower ids (deterministic)
        if {$d < $bd || ($d == $bd && $cid > 16 && $best <= 16)} { set bd $d; set best $cid }
    }
    set cache($k) $best
    return $best
}

proc ::ufv::palette::role {name} {
    variable roles
    variable initialised
    if {!$initialised} { init }
    if {![info exists roles($name)]} { error "unknown palette role '$name'" }
    return $roles($name)
}

proc ::ufv::palette::base {} { variable base_id; return $base_id }

# Group {hex {positions}} pairs by colour id: -> dict id -> positions (merged, unsorted).
proc ::ufv::palette::group {hexgroups} {
    set out [dict create]
    dict for {hex pos} $hexgroups {
        set id [id $hex]
        dict lappend out $id {*}$pos
    }
    return $out
}

# Borrowing the GLOBAL colour-scale method for a molecule's continuous colouring.
#   SC(users)  molecules whose active UFV colouring borrowed the method
#   SC(saved)  the method before the first borrow (restored when the last user releases)
#   SC(ours)   the method we set ("" = the user changed it while borrowed: never override / restore it)
namespace eval ::ufv::palette { variable SC; array set SC {users {} saved "" ours ""} }
proc ::ufv::palette::_scale_now {} { return [lindex [colorinfo scale method] 0] }
proc ::ufv::palette::scale_acquire {mol method} {
    variable SC
    if {![llength [info commands color]] || $method eq ""} return
    set cur [_scale_now]
    if {![llength $SC(users)]} {
        set SC(saved) $cur
        set SC(ours) $method
        if {$cur ne $method} { color scale method $method }
    } elseif {$SC(ours) ne "" && $cur eq $SC(ours)} {
        if {$method ne $cur} { color scale method $method; set SC(ours) $method }
    } else {
        set SC(ours) ""                      ;# the user picked a method while ours was active: theirs wins
    }
    if {$mol ni $SC(users)} { lappend SC(users) $mol }
}
proc ::ufv::palette::scale_release {mol} {
    variable SC
    if {$mol ni $SC(users)} return
    set SC(users) [lsearch -all -inline -not -exact $SC(users) $mol]
    if {[llength $SC(users)]} return
    if {[llength [info commands color]] && $SC(ours) ne "" && $SC(saved) ne ""} {
        set cur [_scale_now]
        if {$cur eq $SC(ours) && $cur ne $SC(saved)} { color scale method $SC(saved) }
    }
    set SC(saved) ""; set SC(ours) ""
}
proc ::ufv::palette::scale_users {} { variable SC; return $SC(users) }

# VMD colour-scale method (GLOBAL in VMD), set unconditionally. Returns the previous method. Colourings use
# scale_acquire / scale_release instead, so the user's method comes back.
proc ::ufv::palette::scale {{method ""}} {
    if {![llength [info commands color]]} { return "" }
    set old [lindex [colorinfo scale method] 0]
    if {$method ne ""} { color scale method $method }
    return $old
}
