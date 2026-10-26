# =============================================================================================
# ufv_gui.tcl - OWNER: the GUI agent. Scope: the Tk panel (ufv_gui), using the async ::ufv::*_async operations and ::ufv::hook events.
#
# Sourced by ufv_vmd.tcl after ufv_backend / ufv_palette / ufv_layers / ufv_commands / ufv_features / ufv_traj.
# Build ONLY on the documented foundation API (see API.md): ::ufv::backend::async|call, ::ufv::await,
# ::ufv::palette::id, ::ufv::map::*, ::ufv::layers::*, ::ufv::cartoon::*, ::ufv::st::*, ::ufv::hook / emit.
# Do not edit the foundation files from here; new backend ops go into ufv_core.py with @_serve_op (additive).
#
# Layout (UI v2, compact): scrolling body of collapsible sections, each with a one-line summary that stays visible
# when collapsed (Protein | Structures | Colouring | Trajectory | Annotations | Residue report) | optional settings
# drawer (report radii, apply-to-all) | fixed footer (status line, settings, Export menu, Clear).
# Annotations v3: ONE list (category drop-down + search + check boxes) - a checked row is shown in 3D.
# Primary actions are buttons; secondary ones live in small menus (the "more" buttons, Filters, Export).
# Every network action goes through ::ufv::backend::async (or a ::ufv::*_async op); the panel never vwaits.
# The panel follows console commands through ::ufv::hook (status fetched loaded mapped layers cartoon report).
# Test entry points: ::ufv::gui::* procs (fetch, load_selected, ann_set_cat, ann_toggle, set_cartoon, show_residue, ...).
# Source is ASCII-only: non-ASCII glyphs are written as \u escapes inside double quotes.
# =============================================================================================

namespace eval ::ufv::gui {
    variable W .ufv
    variable v          ;# widget variables (entries, combos, spins, chips)
    variable S          ;# misc state
    variable full       ;# elided label -> full text
    variable tipt       ;# widget -> tooltip text
    variable linkcmd    ;# report link tag -> command
    variable flowh      ;# toggle grid -> list of hidden children
    variable sec        ;# section key -> frame / open flag
    array set v {}; array set S {}; array set full {}; array set tipt {}; array set linkcmd {}
    array set flowh {}; array set sec {}

    # one accent colour (the PyMOL panel's teal), neutral greys, flat borders
    variable C
    array set C {
        teal #00796b  tealdk #005b50  tealbg #e0f2ef  hdrbg #f3f6f7  fg #263238  grey #6b777e  line #e0e5e8
        border #c9d2d7  btn #edf1f3  hover #e0e7ea  sel #d5ece8  drawer #f6f8f9  on #d3ebe6  led #9aa7ae
        status #a33333  ok #37474f  link #1565c0  bg #ffffff  chipbg #ffffff
        pathogenic #d55e00  deleterious #e6a817  benign #0072b2  uncertain #9e9e9e
        reviewed #42a5f5  clinvar #ab47bc  large #26c6da
    }
    variable U
    array set U [list down "\u25be" right "\u25b8" ext "\u2197" x "\u2715" no "\u2116" on "\u2611" off "\u2610" \
        dot "\u00b7" ell "\u2026" ang "\u00c5" ge "\u2265" le "\u2264" dash "\u2014" bullet "\u25cf" more "\u22ef" \
        gear "\u2699" reset "\u21ba" pick "\u2316" part "\u25a3"]

    # cartoon combo: label -> mode (ufv_cartoon names); "All pockets" only when the features module provides it
    variable cartmodes {None plain Domains domains Topology topology pLDDT plddt B-factor bfactor AlphaMissense alphamissense
        Burden burden Hotspots hotspots {Contact hubs} hubs {Constraint pocket} pockets {All pockets} allpockets}
}

# ---------------------------------------------------------------------------------------------
# Look: fonts, ttk theme + styles, flat widget helpers
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::fonts {} {
    if {"ufvSymSmallB" in [font names]} return
    foreach n [font names] { if {[string match ufv* $n]} { catch { font delete $n } } }
    array set f [font actual TkDefaultFont]
    set fam $f(-family); set sz $f(-size)
    if {$::tcl_platform(platform) eq "windows" && "Segoe UI" in [font families]} { set fam "Segoe UI"; set sz 9 }
    if {$sz < 0} { set sz 9 }
    font create ufvNorm  -family $fam -size $sz
    font create ufvBold  -family $fam -size $sz -weight bold
    font create ufvSmall -family $fam -size [expr {$sz - 1}]
    font create ufvSmallB -family $fam -size [expr {$sz - 1}] -weight bold
    font create ufvH1    -family $fam -size [expr {$sz + 3}] -weight bold
    font create ufvTitle -family $fam -size [expr {$sz + 1}] -weight bold
    font create ufvLink  -family $fam -size $sz -underline 1
    # Glyphs Segoe UI lacks (arrows, gear, ballot boxes, spinner, ...) would make Tk search every installed font
    # for a fallback the first time one is drawn (~110 ms, measured inside the panel build): widgets that show
    # them use Segoe UI Symbol directly (it has them all, and Segoe's Latin letters).
    set sym $fam
    if {$::tcl_platform(platform) eq "windows" && "Segoe UI Symbol" in [font families]} { set sym "Segoe UI Symbol" }
    font create ufvIcon  -family $sym -size [expr {$sz + 1}]
    font create ufvSym   -family $sym -size $sz
    font create ufvSymB  -family $sym -size $sz -weight bold
    font create ufvSymSmall  -family $sym -size [expr {$sz - 1}]
    font create ufvSymSmallB -family $sym -size [expr {$sz - 1}] -weight bold
}

# ttk: the platform theme (vista on Windows; clam instead of the bare X11 'default') + a few Ufv.* styles.
proc ::ufv::gui::theme {} {
    variable C
    if {[catch { ttk::style theme use } cur]} { set cur $::ttk::currentTheme }
    set names [ttk::style theme names]
    if {$::tcl_platform(platform) eq "windows"} {
        if {"vista" in $names && $cur ne "vista"} { catch { ttk::style theme use vista } }
    } elseif {$cur in {default classic} && "clam" in $names} {
        catch { ttk::style theme use clam }
    }
    ttk::style configure Ufv.Treeview -rowheight 20 -font ufvSmall -background white -fieldbackground white -borderwidth 0
    ttk::style map Ufv.Treeview -background [list selected $C(sel)] -foreground [list selected $C(fg)]
    ttk::style configure Ufv.Treeview.Heading -font ufvSmallB -padding {4 1}
    # annotation list: check boxes are symbol glyphs (Segoe UI Symbol has them + Latin); heading = tri-state box
    ttk::style configure Ufv.Ann.Treeview -rowheight 19 -font ufvSymSmall -background white -fieldbackground white -borderwidth 0
    ttk::style map Ufv.Ann.Treeview -background [list selected $C(sel)] -foreground [list selected $C(fg)]
    ttk::style configure Ufv.Ann.Treeview.Heading -font ufvSymSmallB -padding {3 1}
    # loaded-molecule list: tick boxes / delete cross are symbol glyphs
    ttk::style configure Ufv.Struct.Treeview -rowheight 20 -font ufvSymSmall -background white -fieldbackground white -borderwidth 0
    ttk::style map Ufv.Struct.Treeview -background [list selected $C(sel)] -foreground [list selected $C(fg)]
    # uniform 24 px control height (the flat buttons are 24 px)
    ttk::style configure Ufv.TCombobox -padding {5 3 4 4}
    ttk::style configure Ufv.TSpinbox -padding {4 2 2 2}
    ttk::style configure Ufv.TEntry -padding {5 2 4 3}
}

# hover tint for the flat classic buttons (Tk on Windows has no hover state for them)
proc ::ufv::gui::_hover {w on} {
    variable C
    if {![winfo exists $w] || [$w cget -state] eq "disabled"} return
    set key [expr {[info exists ::ufv::gui::S(hoverbg,$w)] ? $::ufv::gui::S(hoverbg,$w) : [list $C(btn) $C(hover)]}]
    $w configure -bg [lindex $key [expr {$on ? 1 : 0}]]
}
proc ::ufv::gui::_flat {w} {
    bind $w <Enter> +[list ::ufv::gui::_hover $w 1]
    bind $w <Leave> +[list ::ufv::gui::_hover $w 0]
}

# Flat filled button (no 3D border): light grey, darker on hover.  -primary 1 = accent fill, white bold text.
proc ::ufv::gui::btn {w text cmd {tip ""} args} {
    variable C
    set primary 0
    if {[dict exists $args -primary]} { set primary [dict get $args -primary]; dict unset args -primary }
    if {$primary} {
        button $w -text $text -command $cmd -relief flat -bd 0 -overrelief flat -padx 12 -pady 3 -highlightthickness 0 \
            -bg $C(teal) -fg white -activebackground $C(tealdk) -activeforeground white -font ufvBold -cursor hand2 {*}$args
        set ::ufv::gui::S(hoverbg,$w) [list $C(teal) $C(tealdk)]
    } else {
        button $w -text $text -command $cmd -relief flat -bd 0 -overrelief flat -padx 9 -pady 3 -highlightthickness 0 \
            -bg $C(btn) -activebackground $C(sel) -fg $C(fg) -activeforeground $C(fg) -font ufvNorm -cursor hand2 {*}$args
    }
    _flat $w
    if {$tip ne ""} { tip $w $tip }
    return $w
}

# Small "more" / dropdown button with its menu ($w.m, returned).
proc ::ufv::gui::mbtn {w text {tip ""} args} {
    variable C
    menubutton $w -text $text -menu $w.m -relief flat -bd 0 -padx 7 -pady 3 -highlightthickness 0 \
        -bg $C(btn) -activebackground $C(hover) -fg $C(fg) -activeforeground $C(fg) -font ufvNorm -cursor hand2 \
        -direction below -indicatoron 0 {*}$args
    menu $w.m -tearoff 0
    _flat $w
    if {$tip ne ""} { tip $w $tip }
    return $w.m
}

# Indicator images for toggles: off = grey outline, on = filled accent (rounded 9 px square).
proc ::ufv::gui::leds {} {
    variable C
    if {"ufvledon" in [image names]} return
    foreach {img col fill} [list ufvledoff $C(led) 0 ufvledon $C(teal) 1] {
        image create photo $img -width 9 -height 9
        if {$fill} {
            $img put $col -to 1 0 8 9
            $img put $col -to 0 1 9 8
        } else {
            $img put $col -to 1 0 8 1; $img put $col -to 1 8 8 9
            $img put $col -to 0 1 1 8; $img put $col -to 8 1 9 8
        }
    }
}

# Toggle (checkbutton without the native indicator): grey fill + outline LED when off, accent tint + filled LED when on.
#   -led 0 : no LED (icon toggles).
proc ::ufv::gui::chip {w text var cmd {tip ""} {colour ""} args} {
    variable C
    if {$colour eq ""} { set colour $C(fg) }
    set led 1
    if {[dict exists $args -led]} { set led [dict get $args -led]; dict unset args -led }
    leds
    checkbutton $w -text $text -variable $var -command $cmd -indicatoron 0 -relief flat -offrelief flat -overrelief flat \
        -bd 0 -padx 8 -pady 3 -highlightthickness 0 -font ufvNorm -bg $C(btn) -selectcolor $C(on) \
        -activebackground $C(hover) -fg $colour -activeforeground $colour -cursor hand2 {*}$args
    if {$led} { $w configure -image ufvledoff -selectimage ufvledon -compound left -padx 5 }
    set ::ufv::gui::S(chipfg,$w) $colour
    _flat $w
    trace add variable $var write [list ::ufv::gui::_chipvar $w]
    _chipstate $w
    if {$tip ne ""} { tip $w $tip }
    return $w
}
proc ::ufv::gui::_chipvar {w args} { after idle [list ::ufv::gui::_chipstate $w] }
proc ::ufv::gui::_chipstate {w} {
    variable C
    if {![winfo exists $w]} return
    set on [expr {[catch { set [$w cget -variable] } x] ? 0 : $x}]
    if {$on eq "" || ![string is integer -strict $on]} { set on 0 }
    set fg $::ufv::gui::S(chipfg,$w)
    set want [expr {$on && $fg eq $C(fg) ? $C(tealdk) : $fg}]
    if {[$w cget -fg] ne $want} { $w configure -fg $want -activeforeground $want }
}

# Label that elides its text to the available width (full text in the tooltip).
proc ::ufv::gui::elabel {w args} {
    label $w -anchor w -justify left -width 1 {*}$args
    bind $w <Configure> [list ::ufv::gui::_elide $w]
    return $w
}
proc ::ufv::gui::settext {w text} {
    variable full
    if {![winfo exists $w]} return
    set full($w) $text
    _elide $w
}
proc ::ufv::gui::_elide {w} {
    variable full; variable tipt
    if {![winfo exists $w] || ![info exists full($w)]} return
    set t $full($w)
    set wd [expr {[winfo width $w] - 6}]
    set f [$w cget -font]
    set one [lindex [split $t \n] 0]
    if {$wd < 20 || [font measure $f $one] <= $wd} {
        set show $one
    } else {
        set lo 0; set hi [string length $one]
        while {$lo < $hi} {
            set mid [expr {($lo + $hi + 1) / 2}]
            if {[font measure $f "[string range $one 0 $mid-1]\u2026"] <= $wd} { set lo $mid } else { set hi [expr {$mid - 1}] }
        }
        set show "[string range $one 0 $lo-1]\u2026"
    }
    if {[$w cget -text] ne $show} { $w configure -text $show }
    if {[info exists ::ufv::gui::S(fixtip,$w)]} {
        set tipt($w) $::ufv::gui::S(fixtip,$w)
    } elseif {$show ne $t} { set tipt($w) $t } else { unset -nocomplain tipt($w) }
    if {![info exists ::ufv::gui::S(tipbound,$w)]} { set ::ufv::gui::S(tipbound,$w) 1; _tipbind $w }
}

# Tooltips
proc ::ufv::gui::tip {w text} {
    variable tipt
    set tipt($w) $text
    _tipbind $w
}
proc ::ufv::gui::_tipbind {w} {
    bind $w <Enter> +[list ::ufv::gui::_tipsched $w]
    bind $w <Leave> +[list ::ufv::gui::_tiphide]
    bind $w <ButtonPress> +[list ::ufv::gui::_tiphide]
}
proc ::ufv::gui::_tipsched {w} {
    variable S
    catch { after cancel $S(tipafter) }
    set S(tipafter) [after 600 [list ::ufv::gui::_tipshow $w]]
}
proc ::ufv::gui::_tipshow {w} {
    variable tipt
    if {![winfo exists $w] || ![info exists tipt($w)] || $tipt($w) eq ""} return
    set t .ufvtip
    catch { destroy $t }
    toplevel $t -bg #90a4ae
    wm overrideredirect $t 1
    catch { wm attributes $t -topmost 1 }
    label $t.l -text $tipt($w) -bg #fbfbf7 -fg #263238 -font ufvSmall -justify left -wraplength 360 -padx 6 -pady 3
    pack $t.l -padx 1 -pady 1
    set x [expr {[winfo pointerx $w] + 12}]; set y [expr {[winfo pointery $w] + 16}]
    wm geometry $t +$x+$y
}
proc ::ufv::gui::_tiphide {} {
    variable S
    catch { after cancel $S(tipafter) }
    catch { destroy .ufvtip }
}

# Toggle grid: equal-width cells in as many columns as fit the width (re-gridded on resize); hidden cells skipped.
proc ::ufv::gui::tgrid {f} {
    variable C
    frame $f -bg $C(bg)
    bind $f <Configure> [list ::ufv::gui::_regrid $f]
    return $f
}
proc ::ufv::gui::flowshow {f child on} {
    variable flowh
    set h [expr {[info exists flowh($f)] ? $flowh($f) : {}}]
    set h [lsearch -all -inline -not -exact $h $child]
    if {!$on} { lappend h $child }
    if {[info exists flowh($f)] && $flowh($f) eq $h} return
    set flowh($f) $h
    _regrid $f 1
}
proc ::ufv::gui::_regrid {f {force 0}} {
    variable flowh; variable S
    if {![winfo exists $f]} return
    set W [winfo width $f]
    if {$W <= 1} { set W [expr {max(300, [winfo width [winfo parent $f]] - 16)}] }
    set hidden [expr {[info exists flowh($f)] ? $flowh($f) : {}}]
    set kids {}; set mw 0
    foreach c [winfo children $f] {
        if {$c in $hidden} { grid forget $c; continue }
        lappend kids $c
        set mw [expr {max($mw, [winfo reqwidth $c])}]
    }
    set gap 4
    set n [expr {max(1, min([llength $kids], ($W + $gap) / max(1, $mw + $gap)))}]
    set sig [list $n $kids]
    if {!$force && [info exists S(grid,$f)] && $S(grid,$f) eq $sig} return
    set S(grid,$f) $sig
    set i 0
    foreach c $kids {
        set col [expr {$i % $n}]
        grid $c -row [expr {$i / $n}] -column $col -sticky ew -padx [list 0 [expr {$col < $n - 1 ? $gap : 0}]] -pady {0 4}
        incr i
    }
    for {set c 0} {$c < 8} {incr c} {
        if {$c < $n} { grid columnconfigure $f $c -weight 1 -uniform t } else { grid columnconfigure $f $c -weight 0 -uniform {} }
    }
    fit_later
}

# ---------------------------------------------------------------------------------------------
# Collapsible sections (flat header: chevron, title, grey one-line summary; 1 px separator)
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::section {key title tipText {expand 0}} {
    variable C; variable U; variable sec; variable S
    set f $S(body).s_$key
    frame $f -bg $C(bg)
    set h [frame $f.h -bg $C(hdrbg) -cursor hand2 -height 24]
    label $h.a -text $U(down) -bg $C(hdrbg) -fg $C(teal) -font ufvSymB -width 2
    label $h.t -text $title -bg $C(hdrbg) -fg $C(fg) -font ufvBold
    elabel $h.s -bg $C(hdrbg) -fg $C(grey) -font ufvSmall
    pack $h.a -side left -padx {4 0}
    pack $h.t -side left
    pack $h.s -side left -fill x -expand 1 -padx {8 6}
    frame $f.sep -height 1 -bg $C(line)
    pack $f.sep -fill x -side top
    pack $h -fill x -side top -ipady 2
    set b [frame $f.b -bg $C(bg) -padx 8 -pady 5]
    pack $b -fill both -expand 1 -side top
    foreach w [list $h $h.a $h.t $h.s] { bind $w <Button-1> [list ::ufv::gui::toggle_section $key] }
    tip $h.t $tipText
    set sec($key) $f
    set sec($key,open) 1
    set sec($key,expand) $expand
    pack $f -fill [expr {$expand ? "both" : "x"}] -expand $expand -side top
    return $b
}
# Small flat tool button inside a section header (right side).
proc ::ufv::gui::hdrtool {key w} {
    variable sec
    pack $w -in $sec($key).h -side right -padx {0 4} -pady 1 -before $sec($key).h.s
}
proc ::ufv::gui::toggle_section {key {open ""}} {
    variable sec; variable U
    set f $sec($key)
    if {$open eq ""} { set open [expr {!$sec($key,open)}] }
    if {$open == $sec($key,open)} return
    set sec($key,open) $open
    if {$open} {
        pack $f.b -fill both -expand 1 -side top
        $f.h.a configure -text $U(down)
        if {$sec($key,expand)} { pack configure $f -fill both -expand 1 }
    } else {
        pack forget $f.b
        $f.h.a configure -text $U(right)
        if {$sec($key,expand)} { pack configure $f -fill x -expand 0 }
    }
    fit_later
}
proc ::ufv::gui::summary {key text} {
    variable sec
    if {[info exists sec($key)]} { settext $sec($key).h.s $text }
}

# Scroll body: inner frame as wide as the canvas, at least as tall as the canvas (so the report fills the rest).
proc ::ufv::gui::fit_later {} {
    variable S
    if {[info exists S(fitpending)]} return
    set S(fitpending) 1
    # a timer, not 'after idle': nested geometry managers propagate their requests over several idle rounds
    after 15 ::ufv::gui::fit
}
proc ::ufv::gui::fit {} {
    variable S
    unset -nocomplain S(fitpending)
    if {![info exists S(canvas)]} return
    set c $S(canvas)
    if {![winfo exists $c]} return
    set cw [winfo width $c]; set ch [winfo height $c]
    if {$cw <= 1} return
    set rh [winfo reqheight $S(body)]
    set h [expr {max($ch, $rh)}]
    $c itemconfigure inner -width $cw -height $h
    $c configure -scrollregion [list 0 0 $cw $h]
    if {$rh <= $ch} {
        $c yview moveto 0
        if {[winfo ismapped $S(vsb)]} { grid remove $S(vsb) }
    } elseif {![winfo ismapped $S(vsb)]} {
        grid $S(vsb)
    }
    # late geometry propagation (nested grids / the trajectory panel): re-fit once more if the body grew / shrank
    catch { after cancel $S(fitcheck) }
    set S(fitcheck) [after 150 [list ::ufv::gui::_fit_check $rh]]
}
proc ::ufv::gui::_fit_check {rh} {
    variable S
    unset -nocomplain S(fitcheck)
    if {[alive] && [winfo reqheight $S(body)] != $rh} { fit_later }
}
# Mouse wheel: anywhere over the panel it scrolls the panel body.  Value widgets (combobox, spinbox, scale) never
# change on the wheel - their wheel goes to the panel; lists / the report scroll themselves while they can.
proc ::ufv::gui::_wheel {w delta} {
    variable S
    if {![info exists S(canvas)] || ![winfo exists $S(canvas)]} return
    if {![winfo exists $w] || [winfo toplevel $w] ne $::ufv::gui::W} return
    if {[winfo class $w] in {Text Treeview Listbox}} return
    if {[lindex [$S(canvas) yview] 0] == 0 && [lindex [$S(canvas) yview] 1] == 1} return
    set n [expr {-$delta / 120}]
    if {$n == 0 && $delta != 0} { set n [expr {$delta > 0 ? -1 : 1}] }
    $S(canvas) yview scroll $n units
}
proc ::ufv::gui::_wheel_tag {w delta} {
    # bound on the UfvWheel tag (before the class tag); return break unless the widget should handle it itself
    if {[winfo class $w] in {Text Treeview Listbox}} {
        lassign [$w yview] a b
        if {!($a == 0 && $b == 1)} { return }
    }
    _wheel_any $w $delta
    return -code break
}
proc ::ufv::gui::_wheel_any {w delta} {
    variable S
    if {![info exists S(canvas)] || ![winfo exists $S(canvas)] || ![winfo exists $w] || [winfo toplevel $w] ne $::ufv::gui::W} return
    if {[lindex [$S(canvas) yview] 0] == 0 && [lindex [$S(canvas) yview] 1] == 1} return
    set n [expr {-$delta / 120}]
    if {$n == 0 && $delta != 0} { set n [expr {$delta > 0 ? -1 : 1}] }
    $S(canvas) yview scroll $n units
}
# Put the UfvWheel tag in front of the class tag of every value / list widget under $top (idempotent).
proc ::ufv::gui::wheelfix {top} {
    if {![winfo exists $top]} return
    if {![llength [bind UfvWheel <MouseWheel>]]} {
        bind UfvWheel <MouseWheel> {::ufv::gui::_wheel_tag %W %D}
        bind UfvWheel <Shift-MouseWheel> {break}
        if {[tk windowingsystem] eq "x11"} {
            bind UfvWheel <Button-4> {::ufv::gui::_wheel_tag %W 120}
            bind UfvWheel <Button-5> {::ufv::gui::_wheel_tag %W -120}
        }
    }
    set todo [list $top]
    while {[llength $todo]} {
        set w [lindex $todo 0]; set todo [lrange $todo 1 end]
        lappend todo {*}[winfo children $w]
        if {[winfo class $w] ni {TCombobox TSpinbox Spinbox TScale Scale Menubutton TMenubutton Text Treeview Listbox}} continue
        set tags [bindtags $w]
        if {"UfvWheel" in $tags} continue
        bindtags $w [linsert $tags [lsearch -exact $tags [winfo class $w]] UfvWheel]
    }
}

# ---------------------------------------------------------------------------------------------
# State helpers
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::alive {} { variable W; return [expr {[llength [info commands winfo]] && [winfo exists $W]}] }

# Active molecule ("" when none).
proc ::ufv::gui::amol {} {
    if {[catch { ::ufv::resolve "" } m]} { return "" }
    return $m
}
# Molecules the layer / colouring actions apply to ("all" = every UFV molecule).
proc ::ufv::gui::targets {} {
    variable v
    set m [amol]
    if {$m eq ""} { return {} }
    if {$v(all)} {
        set out {}
        foreach x [molinfo list] { if {[::ufv::st::has $x]} { lappend out $x } }
        if {$m ni $out} { lappend out $m }
        return $out
    }
    return [list $m]
}
proc ::ufv::gui::uid_of {m} {
    set u ""
    if {$m ne "" && [::ufv::st::has $m]} { set u [::ufv::st::get $m uid] }
    if {$u eq ""} { set u $::ufv::uid }
    return $u
}
proc ::ufv::gui::ann {{uid ""}} {
    if {$uid eq ""} { set uid [uid_of [amol]] }
    if {$uid ne "" && [info exists ::ufv::ann($uid)]} { return $::ufv::ann($uid) }
    return ""
}
proc ::ufv::gui::dg {d key {def ""}} {
    if {[dict exists $d $key]} { return [dict get $d $key] }
    return $def
}
proc ::ufv::gui::status {msg {err 0}} {
    variable S; variable C
    if {![alive]} return
    set w $S(status)
    if {!$err} { set err [regexp -nocase {fail|error|cannot|unreachable|not mapped|no molecule} $msg] }
    $w configure -fg [expr {$err ? $C(status) : $C(ok)}]
    settext $w $msg
}
proc ::ufv::gui::fail {what res} {
    status "$what: $res" 1
}
# Run a console-level command from a button; errors go to the status line, never to bgerror.
proc ::ufv::gui::run {what script} {
    if {[catch { uplevel #0 $script } res]} { fail $what $res; return -code ok "" }
    return $res
}

# Busy spinner (polls the backend's pending count only while something is pending).
proc ::ufv::gui::busy {} {
    variable S
    if {[info exists S(spin)]} return
    set S(spin) 0
    _spin
}
proc ::ufv::gui::_spin {} {
    variable S
    if {![alive]} { unset -nocomplain S(spin); return }
    set n [::ufv::backend::busy]
    if {$n == 0 && ![info exists S(work)]} { $S(spinner) configure -text ""; unset -nocomplain S(spin); return }
    set frames [list "\u25d0" "\u25d3" "\u25d1" "\u25d2"]
    $S(spinner) configure -text "[lindex $frames [expr {[incr S(spin)] % 4}]] $n"
    after 150 ::ufv::gui::_spin
}

# ---------------------------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------------------------
proc ufv_gui {} {
    if {![llength [info commands toplevel]]} { error "ufv_gui needs Tk (VMD started without graphics?)" }
    set W $::ufv::gui::W
    if {[winfo exists $W]} { wm deiconify $W; raise $W; return $W }
    ::ufv::gui::build
    return $W
}

proc ::ufv::gui::build {} {
    variable W; variable C; variable U; variable v; variable S; variable cartmodes; variable markers
    set t0 [clock microseconds]
    fonts
    array set v [list uid $::ufv::uid struct "" cart None q 0.10 all 0 vcolor consequence kind "" \
        filter "" ptmA 8 varA 12 nearby [expr {int($::ufv::nearby)}] pick 0 protnlm 0 settings 0]
    set S(anncat) all

    toplevel $W -bg $C(bg)
    wm title $W "3D Feature Viewer"
    # Fixed window size: set once (remembered size, else 420 x screen), never driven by the content afterwards.
    # Everything that changes (sections, report, table, trajectory panel) lives in the scrolling body or elides.
    wm geometry $W [geom_initial]
    wm minsize $W 340 360
    bind $W <Configure> {if {"%W" eq $::ufv::gui::W} { ::ufv::gui::_geom_seen %w %h }}
    wm protocol $W WM_DELETE_WINDOW [list wm withdraw $W]
    theme

    # ---- footer (fixed): status | spinner | settings | Export | Clear
    set ft [frame $W.ft -bg $C(hdrbg)]
    elabel $ft.st -bg $C(hdrbg) -fg $C(ok) -font ufvSmall
    label $ft.sp -bg $C(hdrbg) -fg $C(teal) -font ufvSymSmall -width 3 -anchor e
    chip $ft.set $U(gear) ::ufv::gui::v(settings) ::ufv::gui::toggle_settings \
        "Settings: report radii (PTM / variant / Nearby), apply layers to all molecules" "" -led 0 -font ufvIcon -padx 5 -pady 0
    set m [mbtn $ft.exp "Export $U(down)" "Export the annotations or the current colouring" -font ufvSym]
    $m add command -label "CSV - per-residue annotations" -command [list ::ufv::gui::export_csv 0]
    $m add command -label "CSV + ProtVar columns (slow)" -command [list ::ufv::gui::export_csv 1]
    $m add separator
    $m add command -label "B-factor PDB (current colouring)" -command [list ::ufv::gui::save_bfactor]
    set S(exportmenu) $m
    btn $ft.clr "Clear" [list ::ufv::gui::clear] "Remove every UFV layer / colouring (your own representations are kept)"
    pack $ft.clr -side right -padx {4 6} -pady 4 -fill y
    pack $ft.exp $ft.set -side right -padx {4 0} -pady 4 -fill y
    pack $ft.sp -side right
    pack $ft.st -side left -fill x -expand 1 -padx {8 0}
    set S(status) $ft.st; set S(spinner) $ft.sp

    # ---- settings drawer (hidden until the gear is on)
    set sd [frame $W.set -bg $C(drawer) -padx 8 -pady 4]
    set r [frame $sd.r -bg $C(drawer)]
    label $r.h -text "Report radii ($U(ang))" -bg $C(drawer) -fg $C(grey) -font ufvSmall
    pack $r.h -side left -padx {0 6}
    foreach {k lbl lo hi tipText} [list ptmA "PTM" 2 30 "Report: PTMs within this CA distance" \
            varA "Var" 2 30 "Report: pathogenic variants within this CA distance" \
            nearby "Nearby" 3 15 "Residue / ligand focus neighbourhood (any-atom distance) and the report's Nearby list - ufv_nearby"] {
        label $r.l$k -text $lbl -bg $C(drawer) -font ufvSmall
        ttk::spinbox $r.$k -from $lo -to $hi -increment 1 -width 3 -textvariable ::ufv::gui::v($k) \
            -command [list ::ufv::gui::on_radius $k] -style Ufv.TSpinbox
        bind $r.$k <Return> [list ::ufv::gui::on_radius $k]
        pack $r.l$k $r.$k -side left -padx {4 0}
        tip $r.l$k $tipText; tip $r.$k $tipText
    }
    pack $r -fill x
    checkbutton $sd.all -text "Layers & colouring act on all loaded molecules" -variable ::ufv::gui::v(all) \
        -bg $C(drawer) -activebackground $C(drawer) -font ufvSmall -highlightthickness 0 -anchor w -command ::ufv::gui::sync_later
    tip $sd.all "Apply to all - layers / colouring act on every loaded UFV molecule, not only the active one"
    pack $sd.all -fill x -pady {3 0}
    set S(settings) $sd

    # ---- scrolling body
    set mid [frame $W.mid -bg $C(bg)]
    set c [canvas $mid.c -highlightthickness 0 -bd 0 -bg $C(bg) -yscrollcommand [list $mid.sb set]]
    ttk::scrollbar $mid.sb -orient vertical -command [list $c yview]
    grid $c -row 0 -column 0 -sticky nsew
    grid $mid.sb -row 0 -column 1 -sticky ns
    grid columnconfigure $mid 0 -weight 1
    grid rowconfigure $mid 0 -weight 1
    set S(canvas) $c; set S(vsb) $mid.sb
    set S(body) [frame $c.f -bg $C(bg)]
    $c create window 0 0 -window $S(body) -anchor nw -tags inner
    bind $c <Configure> ::ufv::gui::fit_later
    bind $W <MouseWheel> {::ufv::gui::_wheel %W %D}

    pack $ft -side bottom -fill x
    frame $W.fsep -height 1 -bg $C(line); pack $W.fsep -side bottom -fill x
    pack $mid -side top -fill both -expand 1

    # ---- Protein: accession + Fetch + more (UniProt page, AI name, re-fetch)
    set b [section protein "Protein" "UniProt accession, entry summary, ProtNLM name, tractability"]
    set r [frame $b.r -bg $C(bg)]
    ttk::entry $r.e -textvariable ::ufv::gui::v(uid) -width 12 -style Ufv.TEntry
    bind $r.e <Return> ::ufv::gui::fetch
    tip $r.e "UniProt accession - press Enter to fetch"
    btn $r.f "Fetch" ::ufv::gui::fetch "Fetch UniProt annotations and the list of structures (Enter)" -primary 1
    set m [mbtn $r.mo $U(more) "More: UniProt page, AI (ProtNLM) name, re-fetch" -font ufvSymB -padx 7]
    $m add command -label "Open UniProt entry $U(ext)" -command ::ufv::gui::open_uniprot
    $m add checkbutton -label "Show AI (ProtNLM) name" -variable ::ufv::gui::v(protnlm) -command ::ufv::gui::show_protnlm
    $m add separator
    $m add command -label "Re-fetch (ignore cache)" -command {::ufv::gui::fetch 1}
    set S(protmenu) $m
    grid $r.e $r.f $r.mo -sticky nsew
    grid configure $r.f $r.mo -padx {4 0}
    grid columnconfigure $r 0 -weight 1
    pack $r -fill x
    elabel $b.info -bg $C(bg) -fg $C(grey) -font ufvSmall
    label $b.pn -bg $C(bg) -fg #5b3a8e -font ufvSmall -justify left -anchor w
    elabel $b.tr -bg $C(bg) -fg $C(link) -font ufvSmall
    pack $b.info -fill x -pady {4 0}
    bind $b <Configure> [list ::ufv::gui::_wrap $b.pn %w]
    set S(info) $b.info; set S(protnlm) $b.pn; set S(tract) $b.tr

    # ---- Structures: combo + Load + more (Load all, Align, Numbering, unmapped)
    set b [section structures "Structures" "Choose / load / align structures; click a loaded molecule to make it active"]
    set r [frame $b.r -bg $C(bg)]
    ttk::combobox $r.c -textvariable ::ufv::gui::v(struct) -state readonly -width 8 -style Ufv.TCombobox
    tip $r.c "Structures of this accession (AlphaFold, PDB, AlphaFill, isoforms)"
    btn $r.l "Load" ::ufv::gui::load_selected "Download and load the selected structure"
    set m [mbtn $r.mo $U(more) "More: load all, align, residue numbering, unmapped residues" -font ufvSymB -padx 7]
    $m add command -label "Load all structures" -command ::ufv::gui::load_all
    $m add command -label "Align all onto the active molecule" -command ::ufv::gui::align
    $m add command -label "Residue numbering$U(ell)" -command ::ufv::gui::numbering_dialog
    $m add separator
    $m add command -label "Highlight unmapped residues (chimera)" -command ::ufv::gui::chimera_toggle -state disabled
    set S(structmenu) $m; set S(chimidx) [$m index end]
    grid $r.c $r.l $r.mo -sticky nsew
    grid configure $r.l $r.mo -padx {4 0}
    grid columnconfigure $r 0 -weight 1
    pack $r -fill x
    set S(structcombo) $r.c
    ttk::treeview $b.list -columns {vis name map x} -show {} -height 1 -selectmode browse -style Ufv.Struct.Treeview
    $b.list column vis -width 22 -stretch 0 -anchor center
    $b.list column name -width 150 -stretch 1
    $b.list column map -width 90 -stretch 0
    $b.list column x -width 22 -stretch 0 -anchor center
    $b.list tag configure active -font ufvSymSmallB -foreground $C(tealdk)
    bind $b.list <ButtonRelease-1> {::ufv::gui::_structclick %W %x %y}
    tip $b.list "Loaded molecules: tick = shown, click a name to make it active, $U(x) deletes"
    set S(structlist) $b.list
    set r2 [frame $b.m -bg $C(bg)]
    elabel $r2.l -bg $C(bg) -fg $C(grey) -font ufvSmall
    pack $r2.l -side left -fill x -expand 1
    set S(maplabel) $r2.l

    # ---- Colouring: cartoon colouring (one drop-down), AlphaFill identity (AlphaFill models only)
    set b [section layers "Colouring" "Cartoon colouring of the active structure (the annotation spheres are in Annotations)"]
    set f0 [frame $b.f0 -bg $C(bg)]
    label $f0.l -text "Colour" -bg $C(bg) -fg $C(grey) -font ufvSmall
    set labels {}
    foreach {lbl md} $cartmodes { lappend labels $lbl }
    ttk::combobox $f0.c -textvariable ::ufv::gui::v(cart) -values $labels -state readonly -width 15 -style Ufv.TCombobox \
        -height 12
    bind $f0.c <<ComboboxSelected>> ::ufv::gui::apply_cartoon
    tip $f0.c "Cartoon colouring of the active structure (analyses run in the background, then are cached)"
    label $f0.ql -text "q $U(le)" -bg $C(bg) -fg $C(grey) -font ufvSmall
    ttk::spinbox $f0.q -from 0.01 -to 0.40 -increment 0.01 -format %.2f -width 5 -textvariable ::ufv::gui::v(q) \
        -command ::ufv::gui::on_q -style Ufv.TSpinbox
    bind $f0.q <Return> ::ufv::gui::on_q
    tip $f0.q "Constraint-pocket sensitivity: FDR threshold (BH q)"
    grid $f0.l $f0.c $f0.ql $f0.q -sticky nsew
    grid configure $f0.c -padx {6 0}
    grid configure $f0.ql -padx {6 2}
    grid columnconfigure $f0 1 -weight 1
    grid remove $f0.ql $f0.q
    set S(cartrow) $f0; set S(qspin) $f0.q
    pack $f0 -fill x
    set fa [frame $b.aff -bg $C(bg)]
    label $fa.l -text "AlphaFill donor identity" -bg $C(bg) -fg $C(grey) -font ufvSmall
    ttk::combobox $fa.c -textvariable ::ufv::gui::v(aff) -state readonly -width 9 -style Ufv.TCombobox \
        -values [lmap p {25 30 40 50 60 70 80 90} { format "%s %d%%" $U(ge) $p }]
    set v(aff) [format "%s 25%%" $U(ge)]
    bind $fa.c <<ComboboxSelected>> ::ufv::gui::aff_select
    tip $fa.c "AlphaFill: show transplants whose donor has at least this sequence identity (25% = all) - ufv_ligand_identity"
    pack $fa.l -side left
    pack $fa.c -side left -padx {6 0}
    set S(affrow) $fa

    # ---- Trajectory (only when the active molecule has more than one frame)
    set b [section traj "Trajectory" "Trajectory / ensemble analyses of the active molecule"]
    set S(trajbody) $b
    pack forget $::ufv::gui::sec(traj)

    # ---- Annotations v3: category drop-down + search + reset + Filters; the list with check boxes (checked = in 3D)
    set b [section annotations "Annotations" "Tick a row to show it in 3D (the header box ticks every row shown). Click a row for its report, double-click (or Enter) to zoom; Space ticks the selected rows. Search: free text, positions (1378, 100-200), mutations (N1378H), \"phrases\", cons: src: dis: cat: (a|b = either)."]
    set r [frame $b.r -bg $C(bg)]
    ttk::combobox $r.k -textvariable ::ufv::gui::v(kind) -state readonly -width 16 -style Ufv.TCombobox -height 10
    bind $r.k <<ComboboxSelected>> ::ufv::gui::ann_cat_selected
    tip $r.k "Category (checked / total rows shown in 3D) - switching never changes the 3D"
    ttk::entry $r.f -textvariable ::ufv::gui::v(filter) -width 6 -style Ufv.TEntry
    tip $r.f "Search: words (AND), 1378 or 100-200, N1378H, \"a phrase\", cons:pathogenic|likely src:clinvar dis:... cat:..."
    bind $r.f <KeyRelease> {::ufv::gui::ann_search_later %K}
    bind $r.f <Return> {::ufv::gui::ann_search_now}
    bind $r.f <Down> {::ufv::gui::sugg_down; break}
    bind $r.f <Escape> {::ufv::gui::sugg_hide}
    bind $r.f <FocusOut> {after 150 {if {[focus] ne ".ufvsugg.l"} ::ufv::gui::sugg_hide}}
    btn $r.x $U(reset) ::ufv::gui::ann_reset "Uncheck all categories" -font ufvIcon -padx 5 -pady 0
    set fm [mbtn $r.mo "Filters $U(down)" "Variant filters: add / remove search terms (consequence, source, reviewed, diseases) and the sphere colours" -font ufvSym]
    $fm configure -postcommand ::ufv::gui::filters_post
    set S(filtmenu) $fm
    grid $r.k $r.f $r.x $r.mo -sticky nsew
    grid configure $r.f $r.x $r.mo -padx {4 0}
    # the drop-down and the search share the spare width (2:3); the drop-down's list is as wide as its longest entry
    grid columnconfigure $r 0 -weight 2
    grid columnconfigure $r 1 -weight 3
    pack $r -fill x
    set S(anncombo) $r.k; set S(annentry) $r.f; set S(annreset) $r.x
    if {![catch { ttk::combobox::PopdownWindow $r.k } pd]} { bind $pd <Map> [list ::ufv::gui::_popwide $r.k $pd] }
    set tf [frame $b.t -bg $C(line) -padx 1 -pady 1]
    ttk::treeview $tf.tv -columns {chk cat pos ann} -displaycolumns {chk cat pos ann} -show headings \
        -height [expr {[winfo screenheight .] < 900 ? 6 : 8}] -selectmode extended -style Ufv.Ann.Treeview \
        -yscrollcommand [list $tf.sb set]
    $tf.tv heading chk -text $U(off) -anchor center -command ::ufv::gui::ann_header_click
    $tf.tv heading cat -text "Cat" -anchor w
    $tf.tv heading pos -text "Pos" -anchor w
    $tf.tv heading ann -text "Annotation" -anchor w
    $tf.tv column chk -width 24 -minwidth 24 -stretch 0 -anchor center
    $tf.tv column cat -width 40 -minwidth 30 -stretch 0
    $tf.tv column pos -width 50 -minwidth 30 -stretch 0
    $tf.tv column ann -width 200 -stretch 1
    # the pool: every row lives under this detached item; the view = the root's child list
    $tf.tv insert {} end -id pool
    $tf.tv detach pool
    ttk::scrollbar $tf.sb -orient vertical -command [list $tf.tv yview]
    pack $tf.sb -side right -fill y
    pack $tf.tv -side left -fill both -expand 1
    pack $tf -fill x -pady {6 0}
    bind $tf.tv <Button-1> {if {[::ufv::gui::ann_click %W %x %y]} break}
    bind $tf.tv <<TreeviewSelect>> ::ufv::gui::on_row
    # a quick second click on a check box is a toggle too (Tk reports it as a double click, not a Button-1)
    bind $tf.tv <Double-Button-1> {if {[::ufv::gui::ann_click %W %x %y]} break; if {[%W identify region %x %y] eq "cell"} { ::ufv::gui::ann_activate [%W identify item %x %y] }}
    bind $tf.tv <Return> {::ufv::gui::ann_activate; break}
    bind $tf.tv <space> {::ufv::gui::ann_space; break}
    bind $tf.tv <Button-3> {::ufv::gui::_table_menu %W %x %y %X %Y}
    set S(table) $tf.tv
    menu $tf.tv.m -tearoff 0
    set S(tablemenu) $tf.tv.m
    ann_columns

    # ---- Residue report (takes the remaining height); Pick 3D + Reset view in the header
    set b [section report "Residue report" "Report of the clicked / picked residue or ligand" 1]
    set h $::ufv::gui::sec(report).h
    btn $h.reset $U(reset) ::ufv::gui::reset_view "Reset view - clear the focus and zoom back to the whole structure" -font ufvIcon -padx 5 -pady 0
    chip $h.pick "Pick 3D" ::ufv::gui::v(pick) ::ufv::gui::toggle_pick "Click any atom in the VMD window to report its residue / ligand"
    $h.pick configure -font ufvSmall -padx 6 -pady 0
    hdrtool report $h.reset
    hdrtool report $h.pick
    set S(pickbtn) $h.pick; set S(resetbtn) $h.reset
    set rf [frame $b.rf -bg $C(line) -padx 1 -pady 1]
    text $rf.t -wrap word -height 5 -width 30 -font ufvNorm -relief flat -bd 0 -highlightthickness 0 -padx 8 -pady 6 \
        -bg white -fg $C(fg) -cursor arrow -yscrollcommand [list $rf.sb set] -spacing1 1
    ttk::scrollbar $rf.sb -orient vertical -command [list $rf.t yview]
    pack $rf.sb -side right -fill y
    pack $rf.t -side left -fill both -expand 1
    pack $rf -fill both -expand 1
    set S(report) $rf.t
    report_tags $rf.t
    report_placeholder

    # hooks + traces
    foreach ev {status fetched loaded mapped layers cartoon report report_data progress tiers annot} {
        ::ufv::hook add $ev [list ::ufv::gui::_on_$ev]
    }
    if {[llength [info commands molinfo]]} {
        trace add variable ::vmd_initialize_structure write ::ufv::gui::_moltrace
        trace add variable ::vmd_trajectory_read write ::ufv::gui::_moltrace
    }
    bind $W <Destroy> {if {"%W" eq $::ufv::gui::W} ::ufv::gui::_destroyed}

    # initial content (a session may already have data)
    if {$::ufv::uid ne "" && [info exists ::ufv::ann($::ufv::uid)]} { _on_fetched $::ufv::uid }
    refresh_structs
    sync_layers
    traj_check
    ann_follow
    summary report "click a residue / row, or Pick 3D"
    status "Ready - enter a UniProt accession and Fetch."
    wheelfix $W
    if {[tk windowingsystem] eq "x11"} {
        bind $W <Button-4> {::ufv::gui::_wheel %W 120}
        bind $W <Button-5> {::ufv::gui::_wheel %W -120}
    }
    update idletasks
    fit
    set S(buildms) [expr {([clock microseconds] - $t0) / 1000.0}]
    # Warm the backend while the user types an accession: the first request spawns Python (~0.2-0.5 s, partly
    # synchronous interpreter discovery), so do it now rather than inside the first Fetch.
    if {![::ufv::backend::running]} { after 300 { catch { ::ufv::backend::async ping {} {} } } }
    return $W
}

# ---- window size preference (user resizes are kept for the session and in a small prefs file)
proc ::ufv::gui::geom_file {} {
    set d [expr {[info exists ::env(UFV_CACHE_DIR)] && $::env(UFV_CACHE_DIR) ne "" ? $::env(UFV_CACHE_DIR) : [file normalize ~]}]
    return [file join $d ufv_vmd_panel.geom]
}
proc ::ufv::gui::geom_initial {} {
    variable geom
    if {![info exists geom]} {
        catch {
            set fh [open [geom_file]]; set g [string trim [read $fh]]; close $fh
            if {[regexp {^\d+x\d+$} $g]} { set geom $g }
        }
    }
    set sh [expr {[winfo screenheight .] - 60}]; set sw [winfo screenwidth .]
    if {[info exists geom] && [regexp {^(\d+)x(\d+)$} $geom -> w h]} {
        return "[expr {max(340, min($w, $sw))}]x[expr {max(360, min($h, $sh))}]"
    }
    return 420x[expr {min(800, $sh)}]
}
proc ::ufv::gui::_geom_seen {w h} {
    variable geom; variable S
    if {$w < 100 || $h < 100} return
    if {[info exists geom] && $geom eq "${w}x$h"} return
    set geom ${w}x$h
    catch { after cancel $S(geomafter) }
    set S(geomafter) [after 1000 ::ufv::gui::_geom_write]
}
proc ::ufv::gui::_geom_write {} {
    variable geom
    if {![info exists geom]} return
    catch { set fh [open [geom_file] w]; puts $fh $geom; close $fh }
}

proc ::ufv::gui::swatch {key hex} {
    set img ufvsw_$key
    if {$img ni [image names]} {
        image create photo $img -width 10 -height 10
        catch { $img put $hex -to 0 0 10 10 }
    }
    return $img
}
# Settings drawer (gear in the footer).
proc ::ufv::gui::toggle_settings {} {
    variable v; variable S; variable W
    if {$v(settings)} {
        pack $S(settings) -side bottom -fill x -before $W.fsep
        frame $W.ssep -height 1 -bg $::ufv::gui::C(line)
        pack $W.ssep -side bottom -fill x -before $S(settings)
    } else {
        pack forget $S(settings)
        catch { destroy $W.ssep }
    }
    fit_later
}

proc ::ufv::gui::_wrap {w width} {
    if {[winfo exists $w]} { $w configure -wraplength [expr {max(100, $width - 16)}] }
}

proc ::ufv::gui::_destroyed {} {
    foreach ev {status fetched loaded mapped layers cartoon report report_data progress tiers annot} {
        ::ufv::hook remove $ev [list ::ufv::gui::_on_$ev]
    }
    catch { destroy .ufvsugg }
    catch { trace remove variable ::vmd_initialize_structure write ::ufv::gui::_moltrace }
    catch { trace remove variable ::vmd_trajectory_read write ::ufv::gui::_moltrace }
    catch { trace remove variable ::vmd_pick_atom write ::ufv::gui::_picked }
    variable S; variable v
    catch { _geom_write }
    foreach k [array names S *after*] { catch { after cancel $S($k) } }
    foreach tr [array names v] {
        foreach t [trace info variable ::ufv::gui::v($tr)] { catch { trace remove variable ::ufv::gui::v($tr) {*}$t } }
    }
    array unset S
}

# ---------------------------------------------------------------------------------------------
# Hooks (console commands keep the panel in sync)
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::_on_status {msg} { if {[alive]} { status $msg } }
proc ::ufv::gui::_on_fetched {uid} {
    variable v; variable S
    if {![alive]} return
    set v(uid) $uid
    set a [ann $uid]
    if {$a eq ""} return
    set nv [llength [dict get $a variants]]
    set np 0
    foreach x [dict get $a variants] { if {[dg $x tok] eq "pathogenic"} { incr np } }
    set rng 0; foreach d [dict get $a domains] { incr rng }
    set txt [format "%s %s %d aa %s PTM %d %s Var %d (path %d) %s Site %d %s Mutag %d %s Dom %d %s AM %d" \
        $uid "\u00b7" [dict get $a seqlen] "\u00b7" [llength [dict get $a ptms]] "\u00b7" $nv $np "\u00b7" \
        [llength [dict get $a sites]] "\u00b7" [llength [dict get $a mutagenesis]] "\u00b7" $rng "\u00b7" [dict size [dict get $a am]]]
    settext $S(info) $txt
    set fn [dg $a function]
    set S(fixtip,$S(info)) "$txt\n\n$fn"
    set ::ufv::gui::tipt($S(info)) $S(fixtip,$S(info))
    summary protein "$uid \u00b7 [dict get $a seqlen] aa"
    show_protnlm
    # structure combo
    set labels {}; set S(keys) {}
    if {[info exists ::ufv::structs($uid)]} {
        foreach s $::ufv::structs($uid) { lappend labels [dg $s label [dg $s key]]; lappend S(keys) [dg $s key] }
    }
    $S(structcombo) configure -values $labels
    if {[llength $labels] && [lsearch -exact $labels $v(struct)] < 0} { set v(struct) [lindex $labels 0] }
    if {![llength $labels]} { set v(struct) "" }
    refresh_structs
    # the disease menu + annotation table in the next event-loop slice (keeps the fetch reply's slice short)
    catch { after cancel $S(fetchrest) }
    set S(fetchrest) [after 0 [list ::ufv::gui::_fetched_rest $uid]]
}
proc ::ufv::gui::_fetched_rest {uid} {
    variable S
    unset -nocomplain S(fetchrest)
    if {![alive] || [ann $uid] eq ""} return
    ann_follow
    tractability $uid
}
proc ::ufv::gui::_on_loaded {m key} {
    if {![alive]} return
    table_follow
    refresh_structs
    sync_layers
    traj_check
}
proc ::ufv::gui::_on_mapped {m} {
    variable S
    if {![alive]} return
    unset -nocomplain S(cov,$m)
    refresh_structs
    if {[winfo exists .ufv.num] && $S(nummol) == $m} { numbering_fill $m }
}
proc ::ufv::gui::_on_layers {m} { sync_later }
proc ::ufv::gui::_on_cartoon {m mode} { sync_later }
proc ::ufv::gui::_on_tiers {m} {}
# Long backend crawls (ProtVar pockets / CSV+ProtVar): progress in the status line.
proc ::ufv::gui::_on_progress {key done total} {
    if {![alive]} return
    set what [expr {[string match pockets_all:* $key] ? "ProtVar pockets" : ([string match protvar:* $key] ? "ProtVar (CSV)" : $key)}]
    status "$what: $done / $total ..."
}
proc ::ufv::gui::_moltrace {args} {
    variable S
    if {![alive]} return
    catch { after cancel $S(molafter) }
    set S(molafter) [after 50 { ::ufv::gui::refresh_structs; ::ufv::gui::traj_check; ::ufv::gui::sync_layers; ::ufv::gui::table_follow }]
}
proc ::ufv::gui::sync_later {} {
    variable S
    if {![alive] || [info exists S(syncpending)]} return
    set S(syncpending) 1
    after idle ::ufv::gui::sync_layers
}

# Cartoon combo / q spin / AlphaFill row / summaries from the layer state of the active molecule.
proc ::ufv::gui::sync_layers {} {
    variable S; variable v; variable cartmodes
    unset -nocomplain S(syncpending)
    if {![alive]} return
    set m [amol]
    set mode [expr {$m ne "" && [::ufv::st::has $m] ? [::ufv::cartoon::mode $m] : "none"}]
    set lbl None
    foreach {l md} $cartmodes { if {$md eq $mode} { set lbl $l } }
    if {$mode ni {none plain} && $lbl eq "None"} {
        set lbl [expr {[dict exists {allpockets "All pockets" chimera "Unmapped (chimera)"} $mode] ? [dict get {allpockets "All pockets" chimera "Unmapped (chimera)"} $mode] : $mode}]
    }
    set v(cart) $lbl
    # constraint-pocket sensitivity spin only while that mode is shown
    set qon [expr {$mode eq "pockets"}]
    if {$qon && ![winfo ismapped $S(qspin)]} { grid $S(cartrow).ql $S(qspin) } elseif {!$qon} { grid remove $S(cartrow).ql $S(qspin) }
    summary layers "[expr {$lbl ne "None" ? $lbl : "plain cartoon"}][expr {$v(all) ? " \u00b7 all molecules" : ""}]"
    # "All pockets" colouring only when the features module provides it
    set apok [llength [info commands ::ufv::pockets_all_async]]
    set labels {}
    foreach {l md} $cartmodes { if {$md ne "allpockets" || $apok} { lappend labels $l } }
    if {[$S(cartrow).c cget -values] ne $labels} { $S(cartrow).c configure -values $labels }
    # AlphaFill donor identity: AlphaFill models only
    set aff [expr {$m ne "" && [::ufv::st::has $m] && [llength [info commands ::ufv::feat::is_alphafill]] && ![catch { ::ufv::feat::is_alphafill $m } x] && $x}]
    if {$aff && [winfo manager $S(affrow)] eq ""} { pack $S(affrow) -fill x -pady {4 0} -after $S(cartrow) } elseif {!$aff} { pack forget $S(affrow) }
    # unmapped (chimera) residues: structures mapped by SIFTS / custom numbering
    set chim [expr {$m ne "" && [::ufv::st::has $m] && [llength [info commands ufv_chimera]] && [::ufv::map::mode $m] ne "identity"}]
    $S(structmenu) entryconfigure $S(chimidx) -state [expr {$chim ? "normal" : "disabled"}]
    ann_follow
    fit_later
}

# ---------------------------------------------------------------------------------------------
# Protein
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::fetch {{force 0}} {
    variable v; variable S
    set uid [string toupper [string trim $v(uid)]]
    if {$uid eq ""} { status "Enter a UniProt accession." 1; return }
    set v(uid) $uid
    busy
    set S(fetch_t0) [clock milliseconds]
    if {[catch { ::ufv::fetch_async $uid $force [list ::ufv::gui::_fetched $uid] } err]} { fail "fetch" $err }
}
proc ::ufv::gui::_fetched {uid ok res} {
    variable S
    if {!$ok} { fail "fetch $uid" $res; return }
    set S(fetch_ms) [expr {[clock milliseconds] - $S(fetch_t0)}]
    # the 'fetched' hook already refreshed the panel
}
proc ::ufv::gui::open_uniprot {} {
    variable v
    set u [string toupper [string trim $v(uid)]]
    if {$u eq ""} { set u $::ufv::uid }
    if {$u eq ""} { status "No accession." 1; return }
    open_url "https://www.uniprot.org/uniprotkb/$u/entry"
}
proc ::ufv::gui::open_url {url} {
    variable S
    set S(lasturl) $url
    if {[info exists S(nobrowser)]} { status "would open $url"; return }
    if {[catch { ::ufv::open_url $url } err]} { fail "open $url" $err }
}
proc ::ufv::gui::show_protnlm {} {
    variable v; variable S
    set a [ann]
    set t [expr {$a ne "" ? [dg $a protnlm] : ""}]
    if {$v(protnlm) && $t ne ""} {
        $S(protnlm) configure -text $t
        pack $S(protnlm) -fill x -after $S(info)
    } else {
        pack forget $S(protnlm)
        if {$v(protnlm) && $a ne ""} { status "No ProtNLM name for this entry." }
    }
    fit_later
}
# Open Targets tractability (features module, optional): one line under the summary.
proc ::ufv::gui::tractability {uid} {
    variable S
    pack forget $S(tract)
    if {[llength [info commands ::ufv::tractability_async]]} {
        catch { ::ufv::tractability_async $uid [list ::ufv::gui::_tract $uid] }
    }
}
proc ::ufv::gui::_tract {uid ok res} {
    variable S; variable U; variable tipt
    if {!$ok || ![alive] || $uid ne $::ufv::gui::v(uid)} return
    if {[catch { dict size $res } n] || !$n || ![llength [dg $res groups]]} { return }
    set parts {}; set tip {}
    foreach g [dg $res groups] {
        lappend parts "[dg $g modality] ([llength [dg $g labels]])"
        lappend tip "[dg $g modality]: [join [dg $g labels] {, }]"
    }
    set S(tracturl) [dg $res url]
    settext $S(tract) "Tractability $U(dash) [join $parts " $U(dot) "] $U(dot) [dg $res symbol] $U(ext)"
    $S(tract) configure -cursor hand2
    bind $S(tract) <Button-1> { ::ufv::gui::open_url $::ufv::gui::S(tracturl) }
    pack $S(tract) -fill x -after $S(info)
    set S(fixtip,$S(tract)) "Open Targets tractability (positive evidence) - click to open\n[join $tip \n]"
    set tipt($S(tract)) $S(fixtip,$S(tract))
    fit_later
}

# ---------------------------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::_selected_key {} {
    variable v; variable S
    set labels [$S(structcombo) cget -values]
    set i [lsearch -exact $labels $v(struct)]
    if {$i < 0 || ![info exists S(keys)]} { return "" }
    return [lindex $S(keys) $i]
}
proc ::ufv::gui::load_selected {} {
    variable v
    set uid [string toupper [string trim $v(uid)]]
    if {$uid eq "" || ![info exists ::ufv::ann($uid)]} { status "Fetch an accession first." 1; return }
    set key [_selected_key]
    if {$key eq ""} { status "No structure selected." 1; return }
    load_key $key $uid ""
}
proc ::ufv::gui::load_key {key uid then} {
    variable S
    busy
    set S(load_t0) [clock milliseconds]
    if {[catch { ::ufv::load_structure_async $key $uid [list ::ufv::gui::_loaded $key $then] } err]} {
        fail "load $key" $err
    }
}
proc ::ufv::gui::_loaded {key then ok res} {
    variable S
    if {!$ok} { fail "load $key" $res } else { set S(load_ms) [expr {[clock milliseconds] - $S(load_t0)}]; set S(lastmol) $res }
    if {$then ne ""} { uplevel #0 $then }
}
proc ::ufv::gui::load_all {} {
    variable v; variable S
    set uid [string toupper [string trim $v(uid)]]
    if {$uid eq "" || ![info exists S(keys)] || ![llength $S(keys)]} { status "Fetch an accession first." 1; return }
    set S(loadall) [list $uid $S(keys) 0]
    _load_next
}
proc ::ufv::gui::_load_next {} {
    variable S
    if {![info exists S(loadall)]} return
    lassign $S(loadall) uid keys i
    if {$i >= [llength $keys]} { unset S(loadall); status "Loaded [llength $keys] structures."; return }
    set S(loadall) [list $uid $keys [expr {$i + 1}]]
    status "loading [expr {$i + 1}] / [llength $keys]: [lindex $keys $i] ..."
    load_key [lindex $keys $i] $uid ::ufv::gui::_load_next
}
proc ::ufv::gui::align {} {
    set m [amol]
    if {$m eq ""} { status "No molecule loaded." 1; return }
    run "align" [list ufv_align $m]
}

# Loaded-molecule list (tick = drawn, bold = active, x = delete).
proc ::ufv::gui::refresh_structs {} {
    variable S; variable U
    if {![alive]} return
    set tv $S(structlist)
    $tv delete [$tv children {}]
    set act [amol]
    set n 0
    foreach m [molinfo list] {
        set name [molinfo $m get name]
        set vis [expr {[molinfo $m get drawn] ? $U(on) : $U(off)}]
        set mp [expr {[::ufv::st::has $m] ? [::ufv::map::mode $m] : "-"}]
        set fr [molinfo $m get numframes]
        if {$fr > 1} { append mp " \u00b7 $fr fr" }
        $tv insert {} end -id m$m -values [list $vis "$m: $name" $mp $U(x)] -tags [expr {$m eq $act ? "active" : ""}]
        incr n
    }
    $tv configure -height [expr {min(4, max(1, $n))}]
    if {$n && ![winfo ismapped $tv]} { pack $tv -fill x -pady {4 0} -after [winfo parent $S(structcombo)] }
    if {!$n} { pack forget $tv }
    if {$act ne "" && [::ufv::st::has $act]} {
        settext $S(maplabel) "Numbering: [::ufv::map::describe $act]"
        set cov [coverage $act]
        summary structures "[molinfo $act get name][expr {$cov ne {} ? " \u00b7 $cov% modelled" : {}}] \u00b7 $n loaded"
    } elseif {$act ne ""} {
        settext $S(maplabel) "mol $act is not a UFV molecule yet - use $U(no) or ufv_use"
        summary structures "[molinfo $act get name] \u00b7 $n loaded"
    } else {
        settext $S(maplabel) ""
        summary structures "nothing loaded"
    }
    # numbering line only when it says more than the list's "identity" column
    set mf [winfo parent $S(maplabel)]
    set show [expr {$act ne "" && !([::ufv::st::has $act] && [::ufv::map::mode $act] eq "identity")}]
    if {$show && [winfo manager $mf] eq ""} { pack $mf -fill x -pady {4 0} } elseif {!$show} { pack forget $mf }
    fit_later
}
proc ::ufv::gui::coverage {m} {
    variable S
    set ver [::ufv::st::get $m mapver 0]
    if {[info exists S(cov,$m)] && [lindex $S(cov,$m) 0] == $ver} { return [lindex $S(cov,$m) 1] }
    set a [ann [uid_of $m]]
    if {$a eq ""} { return "" }
    set seen [dict create]
    if {[catch {
        ::ufv::with_sel s $m "protein and name CA" {
            foreach {ch ri} [join [$s get {chain resid}]] {
                set u [::ufv::map::r2u $m $ch $ri]
                if {$u ne "" && $u >= 1 && $u <= [dict get $a seqlen]} { dict set seen $u 1 }
            }
        }
    }]} { return "" }
    set pct [expr {round(100.0 * [dict size $seen] / max(1, [dict get $a seqlen]))}]
    set S(cov,$m) [list $ver $pct]
    return $pct
}
proc ::ufv::gui::_structclick {tv x y} {
    set item [$tv identify item $x $y]
    if {$item eq ""} return
    set col [$tv identify column $x $y]
    set m [string range $item 1 end]
    if {$m ni [molinfo list]} { refresh_structs; return }
    switch -- $col {
        "#1" { mol [expr {[molinfo $m get drawn] ? "off" : "on"}] $m; refresh_structs }
        "#4" { mol delete $m; status "deleted mol $m." }
        default { set_active $m }
    }
}
proc ::ufv::gui::set_active {m} {
    ::ufv::st::init $m
    set ::ufv::molid $m
    set u [::ufv::st::get $m uid]
    if {$u ne ""} {
        set ::ufv::uid $u
        if {$::ufv::gui::v(uid) ne $u && [info exists ::ufv::ann($u)]} { _on_fetched $u }
    }
    catch { ::ufv::cartoon::ensure $m }
    table_follow
    refresh_structs
    sync_layers
    traj_check
    status "active molecule: $m ([molinfo $m get name])"
}

# ---------------------------------------------------------------------------------------------
# Numbering dialog (per chain: auto | offset | adopt | ignore + accession; identity / SIFTS)
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::numbering_dialog {} {
    variable S; variable C; variable v
    set m [amol]
    if {$m eq ""} { status "No molecule loaded." 1; return }
    ::ufv::st::init $m
    set d .ufv.num
    catch { destroy $d }
    toplevel $d -bg $C(bg)
    wm title $d "Numbering - mol $m"
    wm transient $d .ufv
    set S(nummol) $m
    set uid [uid_of $m]
    set v(num,uid) $uid
    set st [::ufv::st::get $m struct]
    set v(num,pdb) [expr {$st ne "" && [dg $st pdbId] ne "" ? [dg $st pdbId] : [string range [molinfo $m get name] 0 3]}]
    label $d.cur -bg $C(bg) -fg #455a64 -font ufvSmall -justify left -anchor w -wraplength 520
    pack $d.cur -fill x -padx 8 -pady {6 2}
    set r [frame $d.r -bg $C(bg)]
    label $r.al -text "Accession" -bg $C(bg) -font ufvSmall
    ttk::entry $r.a -textvariable ::ufv::gui::v(num,uid) -width 10
    btn $r.id "Identity" [list ::ufv::gui::num_identity] "resid == UniProt position on every protein chain"
    label $r.pl -text "PDB" -bg $C(bg) -font ufvSmall
    ttk::entry $r.p -textvariable ::ufv::gui::v(num,pdb) -width 6
    btn $r.sf "SIFTS" [list ::ufv::gui::num_sifts] "Map with SIFTS for this PDB id (checked against the model's residues)"
    pack $r.al $r.a $r.id $r.pl $r.p $r.sf -side left -padx 2
    pack $r -fill x -padx 6 -pady 2
    frame $d.t -bg $C(bg)
    pack $d.t -fill both -expand 1 -padx 8 -pady 4
    set bb [frame $d.b -bg $C(bg)]
    btn $bb.ap "Apply per-chain mapping" [list ::ufv::gui::num_apply] "auto = keep the known map, else find the offset from the sequence; offset = UniProt position of the chain's first residue; adopt = copy molid/chain's mapping; ignore = no mapping"
    btn $bb.cl "Close" [list destroy $d]
    pack $bb.ap -side left
    pack $bb.cl -side right
    pack $bb -fill x -padx 8 -pady {2 8}
    numbering_fill $m
    return $d
}
proc ::ufv::gui::numbering_fill {m} {
    variable v; variable C
    set d .ufv.num
    if {![winfo exists $d]} return
    $d.cur configure -text "Current: [::ufv::map::describe $m]"
    set t $d.t
    foreach w [winfo children $t] { destroy $w }
    set mp [::ufv::map::get $m]
    set cfg [expr {[dict exists $mp config] ? [dict get $mp config] : {}}]
    set notes [expr {[dict exists $mp notes] && [catch {dict size [dict get $mp notes]}] == 0 ? [dict get $mp notes] : {}}]
    set accs [expr {[dict exists $mp accessions] ? [dict get $mp accessions] : {}}]
    set col 0
    foreach h {Chain Residues Mode Value Accession Note} {
        label $t.h$col -text $h -bg $C(bg) -font ufvSmallB -anchor w
        grid $t.h$col -row 0 -column $col -sticky w -padx 2
        incr col
    }
    set res [::ufv::map::residues $m]
    set row 1
    dict for {ch rl} $res {
        set k [string map {" " _} $ch]
        set first [lindex $rl 0 0]; set last [lindex $rl end 0]
        set c [expr {[dict exists $cfg $ch] ? [dict get $cfg $ch] : {}}]
        set v(num,$k,mode) [expr {[dict exists $c mode] ? [dict get $c mode] : "auto"}]
        set v(num,$k,val) [expr {[dict exists $c start] ? [dict get $c start] : ""}]
        if {$v(num,$k,mode) eq "adopt" && [dict exists $c src]} { set v(num,$k,val) [string map {"mol " ""} [dict get $c src]] }
        set v(num,$k,acc) [expr {[dict exists $c accession] ? [dict get $c accession] : ""}]
        label $t.c$row -text $ch -bg $C(bg) -font ufvBold
        label $t.r$row -text "[llength $rl] ($first-$last)" -bg $C(bg) -font ufvSmall
        ttk::combobox $t.m$row -values {auto offset adopt ignore} -state readonly -width 7 -textvariable ::ufv::gui::v(num,$k,mode)
        ttk::entry $t.v$row -width 8 -textvariable ::ufv::gui::v(num,$k,val)
        ttk::entry $t.a$row -width 10 -textvariable ::ufv::gui::v(num,$k,acc)
        set note ""
        if {[dict exists $notes $ch]} { set note [dict get $notes $ch] }
        if {[dict exists $accs $ch] && [dict get $accs $ch] ne ""} { set note "[dict get $accs $ch] $note" }
        label $t.n$row -text $note -bg $C(bg) -fg $C(grey) -font ufvSmall -anchor w -wraplength 220 -justify left
        grid $t.c$row $t.r$row $t.m$row $t.v$row $t.a$row $t.n$row -row $row -sticky w -padx 2 -pady 1
        tip $t.v$row "offset: UniProt position of the chain's first residue ($first); adopt: molid or molid/chain"
        tip $t.a$row "Other accession for this chain (partner / chimera); empty = [uid_of $m]"
        incr row
    }
    grid columnconfigure $t 5 -weight 1
    wheelfix $d
}
proc ::ufv::gui::_num_uid {m} {
    variable v
    set u [string toupper [string trim $v(num,uid)]]
    if {$u ne "" && $u ne [::ufv::st::get $m uid]} {
        if {![info exists ::ufv::ann($u)]} { status "fetch $u first (Protein section)." 1; return -code error "no annotations for $u" }
        ::ufv::st::put $m uid $u
    }
    return $u
}
proc ::ufv::gui::num_identity {} {
    variable S
    set m $S(nummol)
    if {[catch { _num_uid $m }]} return
    run "numbering" [list ::ufv::map::install $m [dict create mode identity]]
    status "mol $m numbering: identity."
}
proc ::ufv::gui::num_sifts {} {
    variable S; variable v
    set m $S(nummol)
    if {[catch { _num_uid $m } u]} return
    set pdb [string trim $v(num,pdb)]
    if {$pdb eq ""} { status "Enter a PDB id." 1; return }
    busy
    status "SIFTS mapping $pdb ..."
    if {[catch { ::ufv::map_sifts_async $m $pdb $u [list ::ufv::gui::_num_done $m "SIFTS $pdb"] } err]} { fail SIFTS $err }
}
proc ::ufv::gui::num_apply {} {
    variable S; variable v
    set m $S(nummol)
    if {[catch { _num_uid $m }]} return
    set cfg [dict create]
    foreach ch [dict keys [::ufv::map::residues $m]] {
        set k [string map {" " _} $ch]
        set mode $v(num,$k,mode); set val [string trim $v(num,$k,val)]; set acc [string toupper [string trim $v(num,$k,acc)]]
        set c [dict create mode $mode]
        switch -- $mode {
            offset {
                if {![string is integer -strict $val]} { status "chain $ch: offset needs the UniProt position of its first residue." 1; return }
                dict set c start $val
            }
            adopt {
                lassign [split $val /] src sch
                if {$src ni [molinfo list]} { status "chain $ch: adopt needs an existing molid (molid or molid/chain)." 1; return }
                dict set c src $src
                if {$sch ne ""} { dict set c src_chain $sch }
            }
        }
        if {$acc ne ""} { dict set c accession $acc }
        dict set cfg $ch $c
    }
    busy
    status "applying per-chain mapping ..."
    if {[catch { ::ufv::map_custom_async $m $cfg [list ::ufv::gui::_num_done $m "per-chain mapping"] } err]} { fail mapping $err }
}
proc ::ufv::gui::_num_done {m what ok res} {
    if {!$ok} { fail $what $res; return }
    status "mol $m: $what applied ($res)."
    refresh_structs
    numbering_fill $m
}

# ---------------------------------------------------------------------------------------------
# Layers & colouring
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::groups_async {uid layer opts cb} {
    busy
    ::ufv::backend::async layer_groups [dict merge [dict create uid $uid layer $layer] $opts] $cb
}

proc ::ufv::gui::apply_cartoon {} {
    variable v; variable cartmodes
    set mode plain
    foreach {l m} $cartmodes { if {$l eq $v(cart)} { set mode $m } }
    set_cartoon $mode
}
proc ::ufv::gui::set_cartoon {mode} {
    variable v; variable S
    set ms [targets]
    if {![llength $ms]} { status "Load a structure first." 1; sync_later; return }
    if {$mode eq "allpockets"} { all_pockets; sync_later; return }
    foreach m $ms {
        set S(cart_t0) [clock milliseconds]
        switch -- $mode {
            plain - none { ::ufv::cartoon::plain $m }
            plddt   { run pLDDT [list ufv_plddt $m] }
            bfactor { run B-factor [list ufv_bfactor $m] }
            domains - topology - alphamissense - burden {
                if {[::ufv::cartoon::reshow $m $mode]} continue
                set a [ann [uid_of $m]]
                if {$a eq ""} { status "Fetch an accession first." 1; break }
                if {$mode eq "alphamissense" && ![dict size [dict get $a am]]} { status "No AlphaMissense scores."; continue }
                if {$mode eq "burden" && ![llength [dict get $a burden]]} { status "No burden-positive residues."; continue }
                groups_async [uid_of $m] $mode {} [list ::ufv::gui::_cart_groups $m $mode]
            }
            hotspots - hubs - pockets {
                busy
                status "computing $mode ..."
                if {[llength [info commands ::ufv::tier_async]]} {
                    if {$mode eq "pockets" && [llength [info commands ::ufv::pocket_q_set]]} { catch { ::ufv::pocket_q_set $v(q) $m } }
                    if {[catch { ::ufv::tier_async $m $mode [list ::ufv::gui::_cart_tier_n $m $mode] } err]} { fail $mode $err }
                } else {
                    set q [expr {$mode eq "pockets" ? [format %.2f $v(q)] : ""}]
                    if {[catch { ::ufv::tiers_async $m $q [list ::ufv::gui::_cart_tiers $m $mode $q] } err]} { fail $mode $err }
                }
            }
            default {
                if {[catch { ufv_cartoon $mode $m } err]} { fail "cartoon $mode" $err }
            }
        }
    }
    sync_later
}
proc ::ufv::gui::_cart_groups {m mode ok res} {
    variable S
    if {!$ok} { fail "cartoon $mode" $res; sync_later; return }
    if {$m ni [molinfo list]} return
    set rebuild [list ::ufv::cartoon_rebuild $m $mode $mode]
    if {[catch {
        ::ufv::cartoon::categorical $m $mode [dict get $res groups] [dict get $res base] $rebuild
        if {$mode eq "domains"} {
            ::ufv::layers::markers $m dom [dict get $res points] -style {VDW 0.7 12} \
                -rebuild [list ::ufv::marker_layer $m dom domains {} {VDW 0.7 12}]
        }
    } err]} { fail "cartoon $mode" $err; sync_later; return }
    set S(cart_ms) [expr {[clock milliseconds] - $S(cart_t0)}]
    status "mol $m: cartoon coloured by $mode ([dict get $res count] residues)."
    sync_later
}
proc ::ufv::gui::_cart_tiers {m kind q ok res} {
    variable S
    if {!$ok} { fail $kind $res; sync_later; return }
    if {$m ni [molinfo list]} return
    set tiers [dict get $res $kind]
    set hexgroups [dict create]
    dict for {tier hex} [dict get $::ufv::tier_colors $kind] {
        dict for {pos t} $tiers { if {$t eq $tier} { dict lappend hexgroups $hex $pos } }
    }
    if {[catch { ::ufv::cartoon::categorical $m $kind $hexgroups "" [list ::ufv::tier_mode $m $kind $q] } err]} {
        fail $kind $err; sync_later; return
    }
    set S(cart_ms) [expr {[clock milliseconds] - $S(cart_t0)}]
    status "mol $m: $kind ([dict size $tiers] residues)."
    sync_later
}
proc ::ufv::gui::on_q {} {
    variable S
    catch { after cancel $S(qafter) }
    set S(qafter) [after 250 ::ufv::gui::_q_apply]
}
proc ::ufv::gui::all_pockets {} {
    set m [amol]
    if {$m eq ""} { status "Load a structure first." 1; return }
    if {![llength [info commands ::ufv::pockets_all_async]]} { status "All pockets: not available in this build." 1; return }
    busy
    set ::ufv::gui::S(work) 1
    if {[catch { ::ufv::pockets_all_async $m [list ::ufv::gui::_all_pockets $m] } err]} { unset -nocomplain ::ufv::gui::S(work); fail "All pockets" $err }
}
proc ::ufv::gui::_all_pockets {m ok res} {
    variable U
    unset -nocomplain ::ufv::gui::S(work)
    if {!$ok} { fail "All pockets" $res; return }
    sync_later
    if {![alive]} return
    set ps [dg $res pockets]
    rclear
    rput "ProtVar pockets" h1; rput "  [uid_of $m] $U(dot) [llength $ps] pocket(s)\n" sub
    if {![llength $ps]} { rput "No ProtVar pockets detected across the protein.\n" muted; rdone; return }
    foreach p $ps {
        set id [dg $p pocketId]; set col [dg $p color "#00897b"]
        rput "$U(bullet) " [ctag $col]
        rput "Pocket $id  " b
        set st {}
        foreach {k f} {score "score %.0f" meanPlddt "pLDDT %.0f" radGyration "Rg %.1f" buriedness "buried %.2f"} {
            set x [dg $p $k]
            if {$x ne "" && [string is double -strict $x]} { lappend st [format $f $x] }
        }
        rput "[llength [dg $p resid]] res $U(dot) [join $st " $U(dot) "]  " muted
        rlink "show" [list ::ufv::gui::show_pocket $m [dg $p resid]]
        if {[llength [info commands ufv_motif_search]]} {
            rput "  "; rlink "motif $U(ext)" [list ::ufv::gui::motif $m "pocket=$id"]
        }
        rput "\n"
    }
    rdone
    summary report "ProtVar pockets $U(dot) [llength $ps]"
    status "All pockets: [llength $ps] ProtVar pocket(s) coloured on mol $m."
}

# ---------------------------------------------------------------------------------------------
# Trajectory section (the trajectory module's panel, shown when the active molecule has > 1 frame)
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::traj_check {} {
    variable S; variable sec
    if {![alive]} return
    set m [amol]
    set n [expr {$m ne "" ? [molinfo $m get numframes] : 0}]
    set f $sec(traj)
    if {$n > 1} {
        if {![info exists S(trajbuilt)]} {
            if {[llength [info commands ::ufv::traj::panel]]} {
                if {[catch {
                    set tw [::ufv::traj::panel $S(trajbody)]
                    if {[winfo exists $tw] && [winfo manager $tw] eq ""} { pack $tw -fill both -expand 1 }
                    set S(trajpanel) $tw
                    wheelfix $tw
                } err]} {
                    puts "\[UFV\] trajectory panel: $err"
                    label $S(trajbody).err -text "Trajectory panel failed: $err" -bg white -fg #a33 -font ufvSmall -wraplength 360 -justify left
                    pack $S(trajbody).err -fill x
                }
            } else {
                label $S(trajbody).na -text "Trajectory tools are not available in this build (ufv_traj.tcl)." \
                    -bg $::ufv::gui::C(bg) -fg $::ufv::gui::C(grey) -font ufvSmall
                pack $S(trajbody).na -fill x
            }
            set S(trajbuilt) 1
        }
        if {![winfo ismapped $f] && [lsearch -exact [pack slaves $S(body)] $f] < 0} {
            pack $f -fill x -side top -after $sec(layers)
        }
        summary traj "[molinfo $m get name] \u00b7 $n frames"
        if {[llength [info commands ::ufv::traj::panel_mol]]} { catch { ::ufv::traj::panel_mol $m } }
    } elseif {[lsearch -exact [pack slaves $S(body)] $f] >= 0} {
        pack forget $f
    }
    fit_later
}

# ---------------------------------------------------------------------------------------------
# Annotations v3: ONE list - category drop-down + search + check boxes. A checked row is shown in 3D (state:
# ::ufv::annot, per molecule + accession + category; console commands change the same state). Rows, ids and the
# search matcher come from ufv_core (ann_rows / search ops), shared with the PyMOL panel.
# The tree holds every row of every category once (inserted in chunks under a detached pool item); the category
# drop-down and the search only change the root's child list (ONE `children` call).
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::ann_mol {} { return [amol] }
# Rebuild the rows when the accession / molecule / extra row data (ligands, partners, PDBe-KB) changed.
proc ::ufv::gui::ann_follow {} {
    variable S
    if {![alive]} return
    set m [amol]
    set uid [uid_of $m]
    if {$uid eq "" || ![info exists ::ufv::ann($uid)]} { set key [list $uid $m] } else {
        set key [expr {$m ne "" ? [::ufv::annot::extra_sig $m] : [list $uid]}]
        lappend key [string length [dict get $::ufv::ann($uid) sequence]] [llength [dict get $::ufv::ann($uid) variants]]
    }
    if {[info exists S(annkey)] && $S(annkey) eq $key} { ann_counts; return }
    set S(annkey) $key
    ann_reload
}
proc ::ufv::gui::table_follow {} { ann_follow }
proc ::ufv::gui::ann_reload {} {
    variable S
    set tv $S(table)
    set gen [incr S(anngen)]
    catch { after cancel $S(annins) }
    # only the rows inserted so far exist (the insert runs in slices; a reload can come in the middle of it)
    if {[info exists S(ins)] && [dict size $S(ins)]} { $tv delete [dict keys $S(ins)] }
    set S(allids) {}; set S(ins) [dict create]; set S(vis) {}; set S(rowsready) 0
    array unset S rc,*; array unset S rp,*; array unset S on,*; array unset S catids,*
    $tv children pool {}
    $tv children {} {}
    set m [amol]
    set uid [uid_of $m]
    if {$uid eq "" || ![info exists ::ufv::ann($uid)]} { ann_counts; ann_header; return }
    set S(annt0) [clock milliseconds]
    if {[catch { ::ufv::annot::rows_async $m [list ::ufv::gui::_ann_rows $gen $m] } err]} { fail "annotation list" $err }
}
proc ::ufv::gui::_ann_rows {gen m ok rows} {
    variable S; variable C
    if {![alive] || $gen != $S(anngen)} return
    if {!$ok} { fail "annotation list" $rows; return }
    set S(rowsms) [expr {[clock milliseconds] - $S(annt0)}]
    foreach c $::ufv::annot::cats { set S(catids,$c) {} }
    set ids {}
    foreach r $rows {
        set id [dict get $r id]
        set c [dict get $r cat]
        set S(rc,$id) $c; set S(rp,$id) $r
        lappend S(catids,$c) $id
        lappend ids $id
    }
    set S(allids) $ids
    set S(rowsready) 1
    # check-box state of every row (a dict lookup per row)
    ann_load_state $m
    ann_counts
    # rows of the current view first, so the visible list fills in the first slice
    ann_view_ids
    set first $S(vis)
    set fd [dict create]; foreach id $first { dict set fd $id 1 }
    set rest {}; foreach id $ids { if {![dict exists $fd $id]} { lappend rest $id } }
    set S(annins) [after 0 [list ::ufv::gui::_ann_insert $gen [concat $first $rest] 0]]
    # a search typed before the rows arrived
    if {[string trim $::ufv::gui::v(filter)] ne ""} { ann_search_now }
}
proc ::ufv::gui::ann_load_state {m} {
    variable S
    foreach c $::ufv::annot::cats {
        set d [expr {$m ne "" && [::ufv::st::has $m] ? [::ufv::annot::checked $m $c] : {}}]
        foreach id $S(catids,$c) { set S(on,$id) [dict exists $d $id] }
    }
}
proc ::ufv::gui::_ann_insert {gen order i} {
    variable S; variable U; variable C
    if {![alive] || $gen != $S(anngen)} return
    set tv $S(table)
    set end [expr {$i + 400}]
    foreach id [lrange $order $i $end-1] {
        set r $S(rp,$id)
        set col [dict get $r color]
        if {$col eq ""} { set col $C(fg) }
        set tag c[string trimleft $col #]
        if {![info exists S(tag,$tag)]} { $tv tag configure $tag -foreground $col; set S(tag,$tag) 1 }
        $tv insert pool end -id $id -tags $tag \
            -values [list [expr {$S(on,$id) ? $U(on) : $U(off)}] [::ufv::annot::short $S(rc,$id)] [dict get $r poslabel] [dict get $r label]]
        dict set S(ins) $id 1
    }
    set n [llength $order]
    if {$end < $n} {
        set S(annins) [after 1 [list ::ufv::gui::_ann_insert $gen $order $end]]
    } else {
        set S(insms) [expr {[clock milliseconds] - $S(annt0)}]
    }
    ann_apply_vis
}
# Ids of the current view without a search: the category's rows (All = every row).
proc ::ufv::gui::ann_view_ids {} {
    variable S
    if {![info exists S(rowsready)] || !$S(rowsready)} { set S(vis) {}; return }
    if {[string trim $::ufv::gui::v(filter)] ne "" && [info exists S(searchres)] && [lindex $S(searchres) 0] eq [ann_search_key]} {
        set S(vis) [lindex $S(searchres) 1]
        return
    }
    set S(vis) [expr {$S(anncat) eq "all" ? $S(allids) : $S(catids,$S(anncat))}]
}
proc ::ufv::gui::ann_apply_vis {} {
    variable S
    if {![alive]} return
    set ins $S(ins)
    set show {}
    foreach id $S(vis) { if {[dict exists $ins $id]} { lappend show $id } }
    $S(table) children {} $show
    set S(shown) [llength $show]
    ann_header
}
# Tri-state check-all box (visible rows) + "N shown" in the header.
proc ::ufv::gui::ann_header {} {
    variable S; variable U
    if {![alive]} return
    set n 0; set tot 0
    foreach id $S(vis) { incr tot; if {$S(on,$id)} { incr n } }
    set g [expr {$tot == 0 || $n == 0 ? $U(off) : ($n == $tot ? $U(on) : $U(part))}]
    set S(tri) [expr {$tot == 0 || $n == 0 ? 0 : ($n == $tot ? 1 : 2)}]
    $S(table) heading chk -text $g
    $S(table) heading ann -text "Annotation \u00b7 $tot row[expr {$tot == 1 ? {} : {s}}] \u00b7 $n shown"
    set S(hdrshown) $n
}
# Drop-down values (with checked / total counts), section summary, reset tooltip.
proc ::ufv::gui::ann_counts {} {
    variable S; variable v
    if {![alive]} return
    set m [amol]
    set vals {}; set order {}
    set tot 0; set totn 0; set parts {}
    foreach c $::ufv::annot::cats {
        set n [expr {[info exists S(catids,$c)] ? [llength $S(catids,$c)] : 0}]
        set k [expr {$m ne "" && [::ufv::st::has $m] ? [::ufv::annot::count $m $c] : 0}]
        if {$k} { lappend parts "[::ufv::annot::short $c] $k" }
        incr tot $k; incr totn $n
        if {$n > 0} {
            lappend vals [format "%s (%d / %d shown)" [::ufv::annot::label $c] $k $n]; lappend order $c
        } elseif {$m ne "" && [::ufv::st::has $m] && [ann_loadable $m $c]} {
            lappend vals "[::ufv::annot::label $c] (fetch\u2026)"; lappend order $c
        }
    }
    set vals [linsert $vals 0 [format "All (%d / %d shown)" $tot $totn]]
    set order [linsert $order 0 all]
    set S(catorder) $order
    if {[$S(anncombo) cget -values] ne $vals} { $S(anncombo) configure -values $vals }
    set i [lsearch -exact $order $S(anncat)]
    if {$i < 0} { set S(anncat) all; set i 0; ann_columns }
    set v(kind) [lindex $vals $i]
    set S(anntotal) $tot
    summary annotations [expr {[llength $parts] ? [join $parts " \u00b7 "] : "nothing shown in 3D"}]
    set ::ufv::gui::tipt($S(annreset)) "Uncheck all categories \u2014 $tot shown in 3D"
}
# Partners / PDBe-KB rows need a network fetch the first time (only offered where they can exist).
# Widen a combobox's pop-down list to its longest entry (ttk sizes it to the combobox), kept on the screen.
proc ::ufv::gui::_popwide {cb pd} {
    if {![winfo exists $pd]} return
    set f [ttk::style lookup TCombobox -font]
    if {$f eq ""} { set f TkDefaultFont }
    set need 0
    foreach s [$cb cget -values] { set need [expr {max($need, [font measure $f $s])}] }
    incr need 30
    set w [winfo width $pd]
    if {$need <= $w} return
    set x [winfo rootx $pd]
    set x [expr {min($x, [winfo screenwidth $pd] - $need)}]
    wm geometry $pd ${need}x[winfo height $pd]+$x+[winfo rooty $pd]
}
proc ::ufv::gui::ann_loadable {m c} {
    switch -- $c {
        kb  { return [expr {[::ufv::annot::needs_load $m kb] && [uid_of $m] ne ""}] }
        prt { return [expr {[::ufv::annot::needs_load $m prt] && [::ufv::map::mode $m] in {sifts custom manual}}] }
    }
    return 0
}
proc ::ufv::gui::ann_columns {} {
    variable S
    $S(table) configure -displaycolumns [expr {$S(anncat) eq "all" ? {chk cat pos ann} : {chk pos ann}}]
}
# Drop-down: switching category never changes the 3D.
proc ::ufv::gui::ann_cat_selected {} {
    variable S
    set i [$S(anncombo) current]
    if {$i < 0} return
    ann_set_cat [lindex $S(catorder) $i]
}
proc ::ufv::gui::ann_set_cat {cat} {
    variable S
    if {$cat ne "all" && $cat ni $::ufv::annot::cats} { set cat [::ufv::annot::cat_of $cat] }
    if {$cat eq ""} return
    set S(anncat) $cat
    ann_columns
    set m [amol]
    if {$cat ne "all" && $m ne "" && [::ufv::st::has $m] && [::ufv::annot::needs_load $m $cat]} {
        status "fetching [::ufv::annot::label $cat] \u2026"
        busy
        if {[catch { ::ufv::annot::load_async $m $cat [list ::ufv::gui::_ann_loaded $cat] } err]} { fail [::ufv::annot::label $cat] $err }
    }
    ann_counts
    if {[string trim $::ufv::gui::v(filter)] ne ""} { ann_search_now; return }
    ann_view_ids
    ann_apply_vis
}
proc ::ufv::gui::_ann_loaded {cat ok res} {
    if {![alive]} return
    if {!$ok} { fail [::ufv::annot::label $cat] $res; return }
    status "[::ufv::annot::label $cat]: $res row(s)."
    ann_follow
}

# ---- search (debounced 120 ms; the backend's matcher, so both panels select the same rows)
proc ::ufv::gui::ann_search_key {} {
    variable S
    return [list $S(anncat) [string trim $::ufv::gui::v(filter)] [expr {[info exists S(annkey)] ? $S(annkey) : ""}]]
}
proc ::ufv::gui::ann_search_later {{k ""}} {
    variable S
    if {$k in {Up Down Return Escape Tab}} return
    catch { after cancel $S(searchafter) }
    set S(searchpending) 1
    set S(searchafter) [after 120 ::ufv::gui::ann_search_now]
}
proc ::ufv::gui::ann_search_now {} {
    variable S; variable v
    catch { after cancel $S(searchafter) }
    unset -nocomplain S(searchafter)
    set q [string trim $v(filter)]
    set m [amol]
    if {$q eq "" || ![info exists S(rowsready)] || !$S(rowsready) || [uid_of $m] eq ""} {
        unset -nocomplain S(searchpending)
        sugg_hide
        ann_view_ids
        ann_apply_vis
        return
    }
    set S(searchpending) 1
    set tok [incr S(searchtok)]
    set S(search_t0) [clock microseconds]
    if {[catch { ::ufv::annot::search_async $m $S(anncat) $v(filter) [list ::ufv::gui::_ann_searched $tok [ann_search_key]] } err]} {
        unset -nocomplain S(searchpending); fail search $err
    }
}
proc ::ufv::gui::_ann_searched {tok key ok res} {
    variable S
    if {![alive] || $tok != $S(searchtok)} return
    unset -nocomplain S(searchpending)
    if {!$ok} { fail search $res; return }
    set S(searchres) [list $key [dict get $res ids]]
    set S(search_ms) [expr {([clock microseconds] - $S(search_t0)) / 1000.0}]
    if {$key ne [ann_search_key]} return
    set t0 [clock microseconds]
    ann_view_ids
    ann_apply_vis
    set S(filter_ms) [expr {([clock microseconds] - $t0) / 1000.0}]
    set sg [expr {[dict exists $res suggest] ? [dict get $res suggest] : {}}]
    set S(suggest) $sg
    if {[focus] eq $S(annentry)} { sugg_show $sg } else { sugg_hide }
}

# ---- suggestions: a small pop-down list under the search box (Down / Return / click accept, Escape closes)
proc ::ufv::gui::sugg_show {items} {
    variable S; variable C
    set w .ufvsugg
    if {![llength $items]} { sugg_hide; return }
    if {![winfo exists $w]} {
        toplevel $w -bg $C(border) -bd 0
        wm overrideredirect $w 1
        catch { wm attributes $w -topmost 1 }
        listbox $w.l -font ufvSmall -relief flat -bd 0 -highlightthickness 0 -activestyle none -exportselection 0 \
            -selectbackground $C(sel) -selectforeground $C(fg) -bg white -fg $C(fg)
        pack $w.l -padx 1 -pady 1 -fill both -expand 1
        bind $w.l <ButtonRelease-1> {::ufv::gui::sugg_accept [%W get [%W nearest %y]]}
        bind $w.l <Return> {::ufv::gui::sugg_accept [%W get active]}
        bind $w.l <Escape> {::ufv::gui::sugg_hide; focus $::ufv::gui::S(annentry)}
    }
    $w.l delete 0 end
    foreach it $items { $w.l insert end $it }
    $w.l configure -height [expr {min(8, [llength $items])}]
    set e $S(annentry)
    set x [winfo rootx $e]; set y [expr {[winfo rooty $e] + [winfo height $e]}]
    wm geometry $w [winfo width $e]x[expr {min(8, [llength $items]) * 16 + 4}]+$x+$y
    wm deiconify $w; raise $w
    set S(suggshown) $items
}
proc ::ufv::gui::sugg_hide {} {
    variable S
    unset -nocomplain S(suggshown)
    if {[winfo exists .ufvsugg]} { wm withdraw .ufvsugg }
}
proc ::ufv::gui::sugg_down {} {
    if {![winfo exists .ufvsugg] || ![winfo ismapped .ufvsugg]} return
    focus .ufvsugg.l
    .ufvsugg.l selection clear 0 end; .ufvsugg.l selection set 0; .ufvsugg.l activate 0
}
# Replace the last (partial) term of the query with the suggestion.
proc ::ufv::gui::sugg_accept {txt} {
    variable v; variable S
    if {$txt eq ""} return
    set q $v(filter)
    set pre [lindex [qlast $q] 0]
    set v(filter) "$pre$txt "
    sugg_hide
    focus $S(annentry); $S(annentry) icursor end
    ann_search_now
}
# Query lexer (ufv_core._ann_lex): terms split on whitespace and '|' outside double quotes.
#   -> list of {raw {{text qat} ..}}   (qat = index of the first quote in the alternative, -1 unquoted)
proc ::ufv::gui::qlex {q} {
    set terms {}; set alts {}; set buf ""; set qat -1; set raw ""; set inq 0
    foreach ch [split $q ""] {
        if {$ch eq "\""} {
            set inq [expr {!$inq}]
            if {$qat < 0} { set qat [string length $buf] }
            append raw $ch
            continue
        }
        if {!$inq && [string is space $ch]} {
            if {$buf ne "" || [llength $alts] || $qat >= 0} { lappend alts [list $buf $qat]; lappend terms [list $raw $alts] }
            set alts {}; set buf ""; set qat -1; set raw ""
            continue
        }
        append raw $ch
        if {!$inq && $ch eq "|"} { lappend alts [list $buf $qat]; set buf ""; set qat -1; continue }
        append buf $ch
    }
    if {$buf ne "" || [llength $alts] || $qat >= 0} { lappend alts [list $buf $qat]; lappend terms [list $raw $alts] }
    return $terms
}
# {text before the last term, last term} (ufv_core.ann_last_term)
proc ::ufv::gui::qlast {q} {
    if {$q eq "" || ([string is space [string index $q end]] && [regexp -all {"} $q] % 2 == 0)} { return [list $q ""] }
    set t [qlex $q]
    if {![llength $t]} { return [list $q ""] }
    set last [lindex $t end 0]
    if {[string range $q end-[expr {[string length $last] - 1}] end] eq $last} {
        return [list [string range $q 0 end-[string length $last]] $last]
    }
    return [list $q ""]
}
proc ::ufv::gui::_qfacet {text qat} {
    set i [string first : $text]
    if {$i > 0 && ($qat < 0 || $i < $qat)} {
        set f [string tolower [string range $text 0 $i-1]]
        set fac {cons cons consequence cons src src source src prov src provenance src dis dis disease dis cat cat category cat type cat pos pos position pos}
        if {[dict exists $fac $f]} { return [list [dict get $fac $f] [string range $text $i+1 end]] }
    }
    return ""
}
# Facet state of the query (ufv_core.ann_query_facets): cons / src / dis = value list or "*" (no term), reviewed 0/1
proc ::ufv::gui::qfacets {q} {
    set out [dict create cons * src * dis * reviewed 0]
    foreach t [qlex $q] {
        lassign $t raw alts
        lassign [lindex $alts 0] text qat
        set fv [_qfacet $text $qat]
        if {$fv eq ""} continue
        lassign $fv fac first
        set vals [list [string trim $first]]
        foreach a [lrange $alts 1 end] { lappend vals [string trim [lindex $a 0]] }
        set vals [lsearch -all -inline -not -regexp $vals {^(|[Nn][Oo][Nn][Ee])$}]
        if {$fac eq "src" && [string tolower $vals] eq "reviewed"} { dict set out reviewed 1; continue }
        switch -- $fac {
            cons {
                set l {}
                foreach x $vals {
                    set x [string tolower $x]; set k $x
                    foreach c {pathogenic deleterious benign uncertain} { if {[string first $x $c] == 0} { set k $c; break } }
                    lappend l $k
                }
                dict set out cons $l
            }
            src {
                set l {}
                foreach x $vals {
                    set x [string tolower $x]
                    lappend l [expr {[dict exists {uniprot uniprot reviewed uniprot clinvar clinvar large large large-scale large largescale large large_scale large} $x] ? [dict get {uniprot uniprot reviewed uniprot clinvar clinvar large large large-scale large largescale large large_scale large} $x] : $x}]
                }
                dict set out src $l
            }
            dis { dict set out dis $vals }
        }
    }
    return $out
}
proc ::ufv::gui::_qquote {s} {
    set s [string map {\" ""} $s]
    return [expr {[regexp {[\s|:]} $s] ? "\"$s\"" : $s}]
}
# The query with the term of facet key (cons | src | dis | reviewed) replaced (ufv_core.ann_query_set).
#   values "*" removes the term; cons / src with every value also removes it; reviewed: 0/1.
proc ::ufv::gui::qset {q key values} {
    set keep {}
    foreach t [qlex $q] {
        lassign $t raw alts
        lassign [lindex $alts 0] text qat
        set fv [_qfacet $text $qat]
        set fac [lindex $fv 0]
        set isrev [expr {$fac eq "src" && [string tolower [string trim $raw]] eq "src:reviewed"}]
        if {($key eq "reviewed" && $isrev) || ($key ne "reviewed" && $fac eq $key && !$isrev)} continue
        lappend keep $raw
    }
    if {$key eq "reviewed"} {
        if {$values} { lappend keep src:reviewed }
    } elseif {$values ne "*"} {
        set all [dict get {cons {pathogenic deleterious benign uncertain} src {uniprot clinvar large} dis {}} $key]
        if {[llength $all]} {
            set vals {}; foreach x $all { if {$x in $values} { lappend vals $x } }
            if {[llength $vals] < [llength $all]} { lappend keep "$key:[expr {[llength $vals] ? [join $vals |] : {none}}]" }
        } else {
            lappend keep "dis:[expr {[llength $values] ? [join [lmap x $values { _qquote $x }] |] : {none}}]"
        }
    }
    return [join $keep " "]
}

# ---- "Filters" menu: inserts / removes search terms (the search text IS the filter)
proc ::ufv::gui::filters_post {} {
    variable S; variable v; variable C
    set m $S(filtmenu)
    $m delete 0 end
    set f [qfacets $v(filter)]
    $m add command -label "Consequence" -state disabled
    foreach {tok lbl} {pathogenic "Pathogenic / likely pathogenic" deleterious "Predicted deleterious" benign "Benign / likely benign" uncertain "Uncertain significance"} {
        set v(fc,$tok) [expr {[dict get $f cons] eq "*" || $tok in [dict get $f cons]}]
        $m add checkbutton -label $lbl -variable ::ufv::gui::v(fc,$tok) -command [list ::ufv::gui::filter_toggle cons] \
            -image [swatch $tok $C($tok)] -compound left
    }
    $m add separator
    $m add command -label "Source" -state disabled
    foreach {tok lbl sw} {uniprot "UniProt reviewed" reviewed clinvar "ClinVar" clinvar large "Large-scale studies" large} {
        set v(fs,$tok) [expr {[dict get $f src] eq "*" || $tok in [dict get $f src]}]
        $m add checkbutton -label $lbl -variable ::ufv::gui::v(fs,$tok) -command [list ::ufv::gui::filter_toggle src] \
            -image [swatch p$sw $C($sw)] -compound left
    }
    set v(frev) [dict get $f reviewed]
    $m add checkbutton -label "Reviewed variants only" -variable ::ufv::gui::v(frev) -command [list ::ufv::gui::filter_toggle reviewed]
    $m add separator
    # diseases of this accession (cascade)
    set dm $m.dis
    if {![winfo exists $dm]} { menu $dm -tearoff 0 }
    $dm delete 0 end
    set S(fdis) {}
    set a [ann]
    if {$a ne ""} {
        set cnt [dict create]
        foreach x [dict get $a variants] { foreach d [dg $x diseases] { dict incr cnt $d } }
        set colors [dg $a diseaseColors]
        set sel [dict get $f dis]
        set lsel [lmap x [expr {$sel eq "*" ? {} : $sel}] { string tolower $x }]
        $dm add command -label "All diseases" -command [list ::ufv::gui::filter_dis_all 1]
        $dm add command -label "None" -command [list ::ufv::gui::filter_dis_all 0]
        $dm add separator
        set i 0
        foreach d [dict keys $cnt] {
            set v(fd,$i) [expr {$sel eq "*" || [string tolower $d] in $lsel}]
            set hex [expr {[dict exists $colors $d] ? [dict get $colors $d] : "#9e9e9e"}]
            $dm add checkbutton -label "$d ([dict get $cnt $d])" -variable ::ufv::gui::v(fd,$i) -command [list ::ufv::gui::filter_toggle dis] \
                -image [swatch dis[string trimleft $hex #] $hex] -compound left
            lappend S(fdis) $i $d
            incr i
        }
    }
    $m add cascade -label "Diseases" -menu $dm -state [expr {[llength $S(fdis)] ? "normal" : "disabled"}]
    $m add separator
    $m add command -label "Colour variant spheres by" -state disabled
    set m0 [amol]
    set cur consequence
    if {$m0 ne "" && [::ufv::st::has $m0]} {
        set vc [::ufv::st::get $m0 [list annvc [::ufv::annot::_uid $m0]] {}]
        if {[dict exists $vc color]} { set cur [dict get $vc color] }
    }
    foreach p {consequence disease provenance} { if {[string match [string range $p 0 2]* $cur]} { set v(vcolor) $p } }
    foreach p {consequence disease provenance} {
        $m add radiobutton -label "  $p" -value $p -variable ::ufv::gui::v(vcolor) -command ::ufv::gui::filter_color
    }
    $m add separator
    $m add command -label "Clear the search" -command ::ufv::gui::ann_clear_search
}
proc ::ufv::gui::filter_toggle {key} { filter_apply [filter_q $key] }
# The query after a Filters-menu change of facet key (cons src reviewed dis), from the menu's tick state.
proc ::ufv::gui::filter_q {key} {
    variable v; variable S
    switch -- $key {
        cons {
            set on {}; foreach t {pathogenic deleterious benign uncertain} { if {$v(fc,$t)} { lappend on $t } }
            # nothing ticked = no consequence filter (the old filters showed every consequence then)
            set q [qset $v(filter) cons [expr {[llength $on] ? $on : "*"}]]
        }
        src {
            set on {}; foreach t {uniprot clinvar large} { if {$v(fs,$t)} { lappend on $t } }
            set q [qset $v(filter) src $on]
        }
        reviewed { set q [qset $v(filter) reviewed $v(frev)] }
        dis {
            set on {}; set n 0
            foreach {i d} $S(fdis) { incr n; if {$v(fd,$i)} { lappend on $d } }
            set q [qset $v(filter) dis [expr {[llength $on] == $n ? "*" : $on}]]
        }
        default { error "unknown filter facet '$key'" }
    }
    return $q
}
proc ::ufv::gui::filter_dis_all {on} {
    variable v; variable S
    foreach {i d} $S(fdis) { set v(fd,$i) $on }
    filter_apply [qset $v(filter) dis [expr {$on ? "*" : {}}]]
}
proc ::ufv::gui::filter_apply {q} {
    variable v; variable S
    set v(filter) $q
    # variant filters: show the variant rows (unless All is shown)
    if {$S(anncat) ni {var all} && [regexp {(^|\s)(cons|src|dis)[:]} $q]} { set S(anncat) var; ann_columns; ann_counts }
    ann_search_now
}
proc ::ufv::gui::filter_color {} {
    variable v
    foreach m [targets] { if {[::ufv::st::has $m]} { ::ufv::annot::var_color $m $v(vcolor) } }
    status "variant spheres coloured by $v(vcolor)."
}
proc ::ufv::gui::ann_clear_search {} {
    variable v
    set v(filter) ""
    ann_search_now
}

# ---- check boxes (click on the box column, the header box, Space; multi-selection toggles together)
proc ::ufv::gui::ann_toggle {ids on} {
    variable S; variable v
    set m [amol]
    if {$m eq "" || ![::ufv::st::has $m]} { status "Load a structure first." 1; return }
    set ch [dict create]
    foreach id $ids { if {[info exists S(rc,$id)]} { dict lappend ch $S(rc,$id) $id } }
    set changes [dict create]
    dict for {c l} $ch { dict set changes $c [list $l $on] }
    if {![dict size $changes]} return
    set t0 [clock microseconds]
    set ms [list $m]
    if {$v(all)} {
        foreach x [targets] { if {$x != $m && [::ufv::st::has $x] && [uid_of $x] eq [uid_of $m]} { lappend ms $x } }
    }
    foreach x $ms {
        set cc $changes
        if {$x != $m} { foreach c {lig prt} { dict unset cc $c } }
        if {[dict size $cc] && [catch { ::ufv::annot::change $x $cc } err]} { fail "annotations" $err }
    }
    set S(toggle_ms) [expr {([clock microseconds] - $t0) / 1000.0}]
}
# annot event (panel or console): refresh the check boxes of that category, the header, counts and summary
proc ::ufv::gui::_on_annot {m cat} {
    variable S
    if {![alive] || $m != [amol]} return
    if {$cat in {prt kb lig *}} { after idle ::ufv::gui::ann_follow }
    ann_sync $cat
}
proc ::ufv::gui::ann_sync {{cat *}} {
    variable S; variable U
    if {![alive] || ![info exists S(rowsready)] || !$S(rowsready)} { ann_counts; return }
    set m [amol]
    set tv $S(table)
    set ins $S(ins)
    foreach c [expr {$cat eq "*" ? $::ufv::annot::cats : [list $cat]}] {
        if {![info exists S(catids,$c)]} continue
        set d [expr {$m ne "" && [::ufv::st::has $m] ? [::ufv::annot::checked $m $c] : {}}]
        foreach id $S(catids,$c) {
            set on [dict exists $d $id]
            if {$on != $S(on,$id)} {
                set S(on,$id) $on
                if {[dict exists $ins $id]} { $tv set $id chk [expr {$on ? $U(on) : $U(off)}] }
            }
        }
    }
    ann_header
    ann_counts
}
proc ::ufv::gui::ann_click {W x y} {
    variable S
    if {[$W identify region $x $y] ne "cell" || [$W identify column $x $y] ne "#1"} { return 0 }
    set it [$W identify item $x $y]
    if {$it eq "" || ![info exists S(on,$it)]} { return 1 }
    set sel [$W selection]
    set ids [expr {$it in $sel && [llength $sel] > 1 ? $sel : [list $it]}]
    ann_toggle $ids [expr {!$S(on,$it)}]
    return 1
}
proc ::ufv::gui::ann_header_click {} {
    variable S
    ann_header
    if {![llength $S(vis)]} return
    # all visible rows checked -> uncheck them, else check them all (one batched update)
    ann_toggle $S(vis) [expr {$S(tri) != 1}]
}
proc ::ufv::gui::ann_space {} {
    variable S
    set sel [$S(table) selection]
    if {![llength $sel]} { set sel [list [$S(table) focus]] }
    set sel [lsearch -all -inline -not $sel ""]
    if {![llength $sel]} return
    set any 0
    foreach id $sel { if {[info exists S(on,$id)] && !$S(on,$id)} { set any 1; break } }
    ann_toggle $sel $any
}
proc ::ufv::gui::ann_reset {} {
    foreach m [targets] { if {[::ufv::st::has $m]} { run reset [list ufv_reset_annotations $m] } }
}

# ---- clicks: single = report + transient highlight (camera never moves); double / Return = zoom
proc ::ufv::gui::on_row {} {
    variable S
    set sel [$S(table) selection]
    if {[llength $sel] != 1 || ![info exists S(rp,$sel)]} return
    ann_show $sel 0
}
proc ::ufv::gui::ann_activate {{it ""}} {
    variable S
    if {$it eq ""} { set it [$S(table) focus] }
    if {$it eq "" || ![info exists S(rp,$it)]} return
    ann_show $it 1
}
proc ::ufv::gui::ann_show {id zoom} {
    variable S
    set r $S(rp,$id)
    set m [amol]
    set c [dict get $r cat]
    switch -- $c {
        lig {
            if {$m eq ""} return
            set code [string range $id 4 end]
            set cp [lsearch -inline -index 0 [::ufv::ligand_copies $m] $code]
            if {$cp eq ""} { set cp [lsearch -inline -nocase -index 0 [::ufv::ligand_copies $m] $code] }
            if {$cp eq ""} { status "no copy of $code in mol $m." 1; return }
            lassign $cp rn ch ri
            show_ligand $m $ch $ri $rn $zoom
        }
        prt {
            lassign [split $id :] - ch ri
            show_partner $m $ch $ri [dict get $r label] $zoom
        }
        default {
            set p [dict get $r pos]; set e [dict get $r end]
            if {$p eq ""} return
            if {$e ne "" && $e != $p && $c in {dom mutag site ptm}} { show_range $p $e $zoom [dict get $r color] } else {
                if {$zoom} { zoom_residue $p $m } else { show_residue $p $m }
            }
        }
    }
}
proc ::ufv::gui::zoom_residue {pos m} {
    variable S
    if {$m eq "" || ![::ufv::st::has $m]} return
    if {![info exists S(curpos)] || $S(curpos) ne $pos} { show_residue $pos $m }
    if {[catch { ::ufv::feat::focus $m $pos 1 } err]} { status "focus: $err" 1 }
}
proc ::ufv::gui::_table_menu {tv x y X Y} {
    variable S
    set it [$tv identify item $x $y]
    set m $S(tablemenu)
    $m delete 0 end
    if {$it ne ""} {
        if {$it ni [$tv selection]} { $tv selection set $it }
        $m add command -label "Zoom to this row" -command [list ::ufv::gui::ann_activate $it]
        $m add separator
        $m add command -label "Check selected rows" -command { ::ufv::gui::ann_toggle [$::ufv::gui::S(table) selection] 1 }
        $m add command -label "Uncheck selected rows" -command { ::ufv::gui::ann_toggle [$::ufv::gui::S(table) selection] 0 }
        $m add separator
    }
    $m add command -label "Check every row shown" -command { ::ufv::gui::ann_toggle $::ufv::gui::S(vis) 1 }
    $m add command -label "Uncheck every row shown" -command { ::ufv::gui::ann_toggle $::ufv::gui::S(vis) 0 }
    $m add command -label "Clear the search" -command ::ufv::gui::ann_clear_search
    tk_popup $m $X $Y
}

# ---------------------------------------------------------------------------------------------
# Residue / ligand report
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::report_tags {t} {
    variable C
    $t tag configure h1 -font ufvH1 -foreground $C(fg) -spacing3 2
    $t tag configure sub -font ufvSmall -foreground $C(grey)
    $t tag configure h2 -font ufvBold -foreground $C(teal) -spacing1 6 -spacing3 1
    $t tag configure h2red -font ufvBold -foreground #c62828 -spacing1 6 -spacing3 1
    $t tag configure key -foreground $C(grey)
    $t tag configure b -font ufvBold
    $t tag configure small -font ufvSmall -foreground #455a64 -lmargin1 14 -lmargin2 14
    $t tag configure ind -lmargin1 8 -lmargin2 18
    $t tag configure link -foreground $C(link) -font ufvLink
    $t tag configure muted -foreground $C(grey)
    $t tag configure err -foreground $C(status)
    foreach k {pathogenic deleterious benign uncertain} { $t tag configure $k -foreground $C($k) }
    $t tag configure amhi -foreground #b71c1c -font ufvBold
    $t tag configure ammid -foreground #e06666 -font ufvBold
    $t tag configure amlo -foreground #3d85c8 -font ufvBold
    $t tag bind link <Enter> [list $t configure -cursor hand2]
    $t tag bind link <Leave> [list $t configure -cursor arrow]
}
proc ::ufv::gui::rclear {} {
    variable S; variable linkcmd
    set t $S(report)
    $t configure -state normal
    $t delete 1.0 end
    foreach tg [$t tag names] { if {[string match L* $tg] && [string is integer -strict [string range $tg 1 end]]} { $t tag delete $tg } }
    array unset linkcmd
    set S(nlink) 0
    foreach w [winfo children $t] { if {[string match *.spider* $w]} { destroy $w } }
}
# Pocket-character radar (same plot as the extension) embedded in the report. fracs = hydrophobic aromatic acidic
# basic polar as 0..1 (independent fractions, not a partition).
proc ::ufv::gui::rspider {fracs {size 150}} {
    variable S
    set t $S(report)
    set w $t.spider[incr S(nspider)]
    set bg [$t cget -background]
    canvas $w -width $size -height $size -bg $bg -highlightthickness 0 -bd 0
    set cx [expr {$size / 2.0}]; set R [expr {$size * 48.0 / 168}]; set n [llength $fracs]
    set labels {Hydrophobic Aromatic Acidic Basic Polar}
    set pi 3.14159265358979
    foreach f {0.5 1.0} {
        set pts {}
        for {set i 0} {$i < $n} {incr i} {
            set a [expr {-$pi / 2 + $i * 2 * $pi / $n}]
            lappend pts [expr {$cx + $R * $f * cos($a)}] [expr {$cx + $R * $f * sin($a)}]
        }
        $w create polygon {*}$pts -fill {} -outline #cfd6dd
    }
    set poly {}
    for {set i 0} {$i < $n} {incr i} {
        set a [expr {-$pi / 2 + $i * 2 * $pi / $n}]
        $w create line $cx $cx [expr {$cx + $R * cos($a)}] [expr {$cx + $R * sin($a)}] -fill #cfd6dd
        $w create text [expr {$cx + ($R + $size * 15.0 / 168) * cos($a)}] [expr {$cx + ($R + $size * 15.0 / 168) * sin($a)}]             -text [lindex $labels $i] -fill #5f6b76 -font ufvSmall
        set f [expr {max(0.0, min(1.0, double([lindex $fracs $i])))}]
        lappend poly [expr {$cx + $R * $f * cos($a)}] [expr {$cx + $R * $f * sin($a)}]
    }
    $w create polygon {*}$poly -fill #99d5cf -outline #00bfa5 -width 2 -stipple gray50
    $w create polygon {*}$poly -fill {} -outline #00bfa5 -width 2
    foreach {x y} $poly { $w create oval [expr {$x - 2}] [expr {$y - 2}] [expr {$x + 2}] [expr {$y + 2}] -fill #00bfa5 -outline #00bfa5 }
    $t window create end -window $w -padx 4 -pady 2
}
proc ::ufv::gui::rput {text {tags {}}} {
    variable S
    $S(report) insert end $text $tags
}
proc ::ufv::gui::rlink {text cmd {extra {}}} {
    variable S; variable linkcmd
    set tg L[incr S(nlink)]
    set linkcmd($tg) $cmd
    $S(report) insert end $text [list link $tg {*}$extra]
    $S(report) tag bind $tg <Button-1> [list ::ufv::gui::_link $tg]
}
proc ::ufv::gui::_link {tg} {
    variable linkcmd
    if {[info exists linkcmd($tg)]} { after idle [list ::ufv::gui::run link $linkcmd($tg)] }
}
proc ::ufv::gui::rdone {} {
    variable S
    $S(report) configure -state disabled
    $S(report) yview moveto 0
    toggle_section report 1
}
proc ::ufv::gui::report_placeholder {} {
    rclear
    rput "No residue selected.\n" muted
    rput "Click a row in Annotations, turn on Pick 3D and click an atom, or run ufv_residue <pos>.\n" sub
    rdone
}

# Report for a UniProt position: transient highlight (focus layer, markers kept, camera NOT moved) + backend
# report (rendered by the 'report' hook). Zooming is an explicit action (double-click / Enter: zoom_residue).
proc ::ufv::gui::show_residue {pos {m ""}} {
    variable S; variable U
    if {$m eq ""} { set m [amol] }
    set uid [uid_of $m]
    if {$uid eq "" || ![info exists ::ufv::ann($uid)]} { status "Fetch an accession first." 1; return }
    set S(curpos) $pos; set S(curmol) $m; set S(curlig) ""; set S(richpos) ""
    rclear
    rput "$pos" h1; rput "  loading report $U(ell)\n" sub
    rdone
    busy
    set S(rep_t0) [clock milliseconds]
    if {$m ne "" && [::ufv::st::has $m] && [llength [info commands ::ufv::report_async]]} {
        if {[llength [info commands ::ufv::feat::focus]]} {
            if {[catch { ::ufv::feat::focus $m $pos 0 0 } err]} { status "focus: $err" 1 }
        } elseif {[catch { ufv_focus $pos $m 0 } err]} { status "focus: $err" 1 }
        if {[catch { ::ufv::report_async $m $pos [list ::ufv::gui::_residue_cb $pos] } err]} { fail "report" $err }
        return
    }
    if {$m ne ""} { if {[catch { ufv_focus $pos $m 0 } err]} { status "focus: $err" 1 } }
    if {[catch { ::ufv::residue_async $pos $uid [list ::ufv::gui::_residue_cb $pos] } err]} { fail "report" $err }
}
proc ::ufv::gui::_residue_cb {pos ok res} {
    variable S
    if {[info exists S(rep_t0)]} { set S(rep_ms) [expr {[clock milliseconds] - $S(rep_t0)}] }
    if {!$ok} { fail "report $pos" $res }
}

# Render the backend's report text (+ structure-derived blocks computed here).
proc ::ufv::gui::render_report {pos text} {
    variable S; variable U; variable v
    if {![alive]} return
    if {![info exists S(curpos)] || $S(curpos) ne $pos} { set S(curpos) $pos; set S(curmol) [amol]; set S(curlig) "" }
    set m $S(curmol)
    if {$m ne "" && $m ni [molinfo list]} { set m [amol] }
    set uid [uid_of $m]
    rclear
    set lines [split [string trimright $text] \n]
    set head [string trim [lindex $lines 0]]
    regexp {^(\S+)\s*(.*)$} $head -> h1 rest
    if {![info exists h1]} { set h1 $pos; set rest "" }
    rput $h1 h1
    set hasPath 0
    foreach l $lines { if {[string match -nocase "*Variant*pathogenic*" $l] && ![string match -nocase "*Likely benign*" $l]} { set hasPath 1 } }
    if {$hasPath} { rput " $U(bullet)" pathogenic }
    rput "  $rest\n" sub
    set heading ""
    set heads {Variant Variants PTM PTMs Domain {Domains & regions} Region {Domains & regions} Site Sites
        Mutagenesis Mutagenesis AlphaMissense Scores ProtVar Scores Pocket Pockets PDBe-KB PDBe-KB Disease Diseases
        Tractability Tractability Ligand Ligands Partner Partners}
    foreach l [lrange $lines 1 end] {
        if {[string trim $l] eq ""} continue
        set lead [expr {[string length $l] - [string length [string trimleft $l]]}]
        set s [string trim $l]
        set w0 ""
        regexp {^([A-Za-z][A-Za-z-]*)} $s -> w0
        if {$lead <= 2 && $w0 ne "" && [dict exists $heads $w0]} {
            set hh [dict get $heads $w0]
            if {$hh ne $heading} { set heading $hh; rput "$hh\n" [expr {$hh eq "Variants" ? "h2red" : "h2"}] }
        }
        if {$lead > 2} {
            # indented detail line (per-substitution predictors etc.)
            if {[regexp {^([A-Z]\d+[A-Z*]):\s*(.*)$} $s -> mut rest2]} {
                if {$heading ne "Predictors"} { set heading Predictors; rput "Per-substitution predictors\n" h2 }
                rput $mut "small b"; rput ": $rest2\n" small
            } else { rput "$s\n" small }
        } elseif {[regexp "^Variant (\\S+) \u2014 (\[^\u2014\]+)(?: \u2014 (.*))?$" $s -> mut cons more]} {
            set cons [string trim $cons]
            set tok [cons_tok $cons]
            rput "  "
            rlink $mut [list ::ufv::gui::open_url "https://www.ebi.ac.uk/ProtVar/query?search=$uid+$mut"] [list b $tok]
            rput "  $cons" $tok
            if {$more ne ""} { rput "  $more" muted }
            rput "\n"
        } elseif {[regexp {^AlphaMissense mean:\s*([0-9.]+)(.*)$} $s -> am rest2]} {
            set tg [expr {$am >= 0.564 ? "amhi" : ($am < 0.34 ? "amlo" : "ammid")}]
            rput "  AlphaMissense mean " key; rput $am $tg; rput "$rest2\n"
        } elseif {[regexp {^([A-Za-z][A-Za-z -]*?):\s*(.*)$} $s -> k val]} {
            rput "  $k " key; rput "$val\n"
        } else {
            rput "  $s\n"
        }
    }
    if {$m ne "" && [::ufv::st::has $m]} { report_structure $m $pos }
    # links
    rput "\n"
    rlink "UniProt $U(ext)" [list ::ufv::gui::open_url "https://www.uniprot.org/uniprotkb/$uid/entry#sequences"]
    rput "   "
    rlink "AlphaFold $U(ext)" [list ::ufv::gui::open_url "https://alphafold.ebi.ac.uk/entry/$uid"]
    rput "   "
    rlink "ProtVar $U(ext)" [list ::ufv::gui::open_url "https://www.ebi.ac.uk/ProtVar/query?search=$uid+$pos"]
    rput "\n"
    rdone
    summary report "$h1 \u00b7 $uid"
}
proc ::ufv::gui::cons_tok {c} {
    set c [string tolower $c]
    if {[string match *pathogenic* $c]} { return pathogenic }
    if {[string match *benign* $c]} { return benign }
    if {[string match *deleterious* $c]} { return deleterious }
    return uncertain
}

# CA coordinates of the protein (cached per molecule / numbering / frame): "chain resid" -> xyz.
proc ::ufv::gui::ca_xyz {m} {
    variable S
    set key [list [::ufv::st::get $m mapver 0] [molinfo $m get frame]]
    if {[info exists S(ca,$m)] && [lindex $S(ca,$m) 0] eq $key} { return [lindex $S(ca,$m) 1] }
    set d [dict create]
    ::ufv::with_sel s $m "protein and name CA" {
        foreach {ch ri x y z} [join [$s get {chain resid x y z}]] { dict set d "$ch $ri" [list $x $y $z] }
    }
    set S(ca,$m) [list $key $d]
    return $d
}
# CA positions of a UniProt position on every mapped chain.
proc ::ufv::gui::uni_xyz {m ca pos} {
    set out {}
    foreach ch [::ufv::map::chains $m] {
        set r [::ufv::map::u2r $m $ch $pos]
        if {$r ne "" && [dict exists $ca "$ch $r"]} { lappend out [dict get $ca "$ch $r"] }
    }
    return $out
}
proc ::ufv::gui::mindist {xa xb} {
    set best ""
    foreach a $xa { foreach b $xb {
        set d [vecdist $a $b]
        if {$best eq "" || $d < $best} { set best $d }
    } }
    return $best
}

# Structure-derived blocks: where the residue is, Nearby residues / ligands, PTMs and pathogenic variants within N A.
proc ::ufv::gui::report_structure {m pos} {
    variable v; variable U
    set a [ann [uid_of $m]]
    set where {}
    foreach ch [::ufv::map::chains $m] {
        set r [::ufv::map::u2r $m $ch $pos]
        if {$r ne ""} { lappend where "$ch:$r" }
    }
    rput "Structure\n" h2
    if {![llength $where]} {
        rput "  not modelled in mol $m ([molinfo $m get name])\n" muted
        return
    }
    rput "  mol $m [molinfo $m get name] " key; rput "\u00b7 chain:resid [join $where {, }]\n"
    set s [::ufv::map::sel $m [list $pos]]
    set d $v(nearby)
    # Nearby (any atom within d A)
    set prot {}; set het {}
    catch {
        ::ufv::with_sel x $m "same residue as ((not water) and within $d of ($s)) and not ($s)" {
            set seen [dict create]
            foreach {ch ri rn isp} [join [$x get {chain resid resname protein}]] {
                if {[dict exists $seen "$ch $ri $rn"]} continue
                dict set seen "$ch $ri $rn" 1
                if {$isp} { lappend prot [list $ch $ri $rn] } else { lappend het [list $ch $ri $rn] }
            }
        }
    }
    rput "Nearby (\u2264 $d $U(ang))\n" h2
    if {![llength $prot] && ![llength $het]} { rput "  nothing\n" muted }
    if {[llength $prot]} {
        rput "  "
        set i 0
        foreach p $prot {
            lassign $p ch ri rn
            set u [::ufv::map::r2u $m $ch $ri]
            set lbl [expr {$u ne {} ? "[aa1 $rn]$u" : "$ch:[aa1 $rn]$ri"}]
            if {$u ne ""} { rlink $lbl [list ::ufv::gui::show_residue $u $m] } else { rput $lbl muted }
            rput "  "
            incr i
        }
        rput "\n"
    }
    foreach h $het {
        lassign $h ch ri rn
        rput "  "; rlink "$rn $ch:$ri" [list ::ufv::gui::show_ligand $m $ch $ri $rn] b; rput "  ligand / ion\n" muted
    }
    if {$a eq ""} return
    # PTMs and pathogenic variants by CA distance
    set ca [ca_xyz $m]
    set me [uni_xyz $m $ca $pos]
    if {![llength $me]} return
    foreach {kind cut} [list ptm $v(ptmA) var $v(varA)] {
        set hits {}
        if {$kind eq "ptm"} {
            set lab [dict create]
            foreach x [dict get $a ptms] { dict lappend lab [dg $x pos] [dg $x desc] }
        } else {
            set lab [dict create]
            foreach x [dict get $a variants] { if {[dg $x tok] eq "pathogenic"} { dict lappend lab [dg $x pos] "[dg $x wt][dg $x pos][dg $x mut]" } }
        }
        dict for {p labels} $lab {
            if {$p == $pos} continue
            set xyz [uni_xyz $m $ca $p]
            if {![llength $xyz]} continue
            set dd [mindist $me $xyz]
            if {$dd ne "" && $dd <= $cut} { lappend hits [list $dd $p [lsort -unique $labels]] }
        }
        set hits [lsort -real -index 0 $hits]
        rput [expr {$kind eq "ptm" ? "PTMs within $cut $U(ang)\n" : "Pathogenic variants within $cut $U(ang)\n"}] [expr {$kind eq "ptm" ? "h2" : "h2red"}]
        if {![llength $hits]} { rput "  none\n" muted; continue }
        set i 0
        foreach h $hits {
            if {$i >= 12} { rput "  + [expr {[llength $hits] - 12}] more\n" muted; break }
            lassign $h dd p labels
            rput "  "
            rlink $p [list ::ufv::gui::show_residue $p $m] [expr {$kind eq "var" ? "pathogenic" : ""}]
            rput "  [format %.1f $dd] $U(ang)  " muted
            rput "[join [lrange $labels 0 3] {, }][expr {[llength $labels] > 4 ? " $U(ell)" : {}}]\n"
            incr i
        }
    }
}
proc ::ufv::gui::aa1 {rn} {
    set m {ALA A ARG R ASN N ASP D CYS C GLN Q GLU E GLY G HIS H HSD H HSE H HSP H HIE H HID H ILE I LEU L LYS K MET M
        PHE F PRO P SER S THR T TRP W TYR Y VAL V SEC U PYL O MSE M}
    if {[dict exists $m $rn]} { return [dict get $m $rn] }
    return "$rn"
}

# Ligand report (structure ligand / ion): highlight it + pocket residues + links (zoom 1: also zoom to it).
proc ::ufv::gui::show_ligand {m ch ri rn {zoom 0}} {
    variable S; variable v; variable U
    if {$m ni [molinfo list]} return
    set S(curlig) [list $m $ch $ri $rn]; set S(curpos) ""; set S(richpos) ""
    set ls "chain [::ufv::map::chain_token $ch] and resid [::ufv::map::compress [list $ri]] and resname $rn"
    set d $v(nearby)
    if {[llength [info commands ::ufv::feat::ligand_focus]]} {
        if {[catch { ::ufv::feat::ligand_focus $m $rn $ch $ri $zoom [expr {$zoom ? 1 : 0}] } err]} { status "ligand focus: $err" 1 }
    } elseif {[catch {
        ::ufv::layers::define $m focus [list [list Name "same residue as (($ls) or ((protein or nucleic) and within $d of ($ls)))"]] \
            -style {Licorice 0.3 12 12} -kind other
        if {$zoom} { ::ufv::zoom_to $m $ls }
    } err]} { status "ligand focus: $err" 1 }
    set S(ligpend) ""
    if {[llength [info commands ::ufv::ligand_data_async]]} {
        set S(ligpend) $S(curlig)
        busy
        if {[catch { ::ufv::ligand_data_async $m $rn $ch $ri [list ::ufv::gui::_ligand_data $m $ch $ri $rn] } err]} {
            set S(ligpend) ""; fail "ligand data" $err
        }
        catch { after cancel $S(ligafter) }
        # safety net only: ligand_data_async always answers within its own 30 s join timeout (partial data)
        set S(ligafter) [after 35000 [list ::ufv::gui::_ligand_timeout $S(curlig)]]
    }
    render_ligand $m $ch $ri $rn ""
    status "mol $m: ligand $rn $ch:$ri."
}
proc ::ufv::gui::_ligand_more {m ch ri rn ok res} {
    variable S
    if {!$ok || ![alive] || $S(curlig) ne [list $m $ch $ri $rn]} return
    set txt $res
    if {[catch { dict size $res }] == 0 && [dict exists $res text]} { set txt [dict get $res text] }
    if {$txt eq ""} return
    set t $S(report)
    $t configure -state normal
    $t insert end "\nMore\n" h2 "$txt\n" small
    $t configure -state disabled
}

# Range row (domain / region / ranged feature): highlight the range (zoom 1: zoom to it).
proc ::ufv::gui::show_range {s e {zoom 1} {hex ""}} {
    set m [amol]
    if {$m eq ""} { status "Load a structure first." 1; return }
    set ps {}; for {set p $s} {$p <= $e} {incr p} { lappend ps $p }
    set sel [::ufv::map::sel $m $ps]
    if {$sel eq ""} { status "$s-$e not modelled in mol $m."; return }
    hl_sel $m "($sel) and name CA" [expr {$hex eq "" ? "#00796b" : $hex}] {Trace 0.5 12}
    if {$zoom} { ::ufv::zoom_to $m $sel; status "mol $m: zoomed to $s-$e." } else { status "mol $m: $s-$e highlighted (double-click to zoom)." }
}
# Partner-protein residue row: highlight it (zoom 1: zoom to it).
proc ::ufv::gui::show_partner {m ch ri text {zoom 0}} {
    if {$m eq "" || $m ni [molinfo list]} return
    set sel "chain [::ufv::map::chain_token $ch] and resid [::ufv::map::compress [list $ri]]"
    hl_sel $m "($sel) and not hydrogen" #7e57c2 {Licorice 0.3 12 12}
    if {$zoom} { ::ufv::zoom_to $m $sel 6 }
    status "$ch:$ri $text"
}
# Transient highlight (the focus layer): drops the previous focus; the camera is not moved.
proc ::ufv::gui::hl_sel {m sel hex style} {
    variable S
    set S(curpos) ""; set S(curlig) ""; set S(richpos) ""
    if {[llength [info commands ::ufv::feat::clear_focus]]} { catch { ::ufv::feat::clear_focus $m } }
    catch { ::ufv::layers::define $m focus [list [list [::ufv::palette::id $hex] $sel]] -style $style -kind other }
}

proc ::ufv::gui::on_radius {k} {
    variable v; variable S
    catch { after cancel $S(radafter) }
    set S(radafter) [after 250 [list ::ufv::gui::_radius $k]]
}
proc ::ufv::gui::_radius {k} {
    variable v; variable S
    if {![string is double -strict $v($k)]} return
    if {$k eq "nearby"} {
        set m [amol]
        if {[catch { ufv_nearby $v(nearby) $m }]} { catch { ufv_nearby $v(nearby) } }
        set v(nearby) [expr {int($::ufv::nearby)}]
    }
    if {$k in {ptmA varA} && [llength [info commands ufv_prox]]} { catch { ufv_prox $v(ptmA) $v(varA) } }
    if {[info exists S(curlig)] && $S(curlig) ne ""} { show_ligand {*}$S(curlig) 0 } elseif {[info exists S(curpos)] && $S(curpos) ne ""} { show_residue $S(curpos) $S(curmol) }
}

# Pick 3D: VMD pick mode + trace on vmd_pick_atom.
proc ::ufv::gui::toggle_pick {} {
    variable v
    if {$v(pick)} {
        trace add variable ::vmd_pick_atom write ::ufv::gui::_picked
        catch { mouse mode pick }
        status "Pick 3D on - click an atom in the VMD window."
    } else {
        catch { trace remove variable ::vmd_pick_atom write ::ufv::gui::_picked }
        catch { mouse mode rotate }
        status "Pick 3D off."
    }
}
proc ::ufv::gui::_picked {args} {
    variable S
    if {![alive]} return
    set m $::vmd_pick_mol; set i $::vmd_pick_atom
    after idle [list ::ufv::gui::pick_atom $m $i]
}
proc ::ufv::gui::pick_atom {m i} {
    if {$m ni [molinfo list]} return
    lassign [lindex [::ufv::with_sel x $m "index $i" { $x get {chain resid resname protein water} }] 0] ch ri rn isp isw
    if {$isw} { status "picked water $ch:$ri."; return }
    if {!$isp} { show_ligand $m $ch $ri $rn; return }
    if {[::ufv::st::has $m] && $m != [amol]} { set_active $m }
    if {![::ufv::st::has $m]} { status "mol $m is not a UFV molecule - use $::ufv::gui::U(no) / ufv_use." 1; return }
    set u [::ufv::map::r2u $m $ch $ri]
    if {$u eq ""} { status "picked $rn $ch:$ri - no UniProt mapping for this residue." 1; return }
    show_residue $u $m
}
proc ::ufv::gui::reset_view {} {
    set m [amol]
    if {$m eq ""} return
    run "reset view" [list ufv_resetview $m]
}

# ---------------------------------------------------------------------------------------------
# Footer: CSV / B-factor PDB / Clear
# ---------------------------------------------------------------------------------------------
proc ::ufv::gui::export_csv {protvar {path ""}} {
    variable v; variable S
    set m [amol]
    set uid [uid_of $m]
    if {$uid eq ""} { set uid [string toupper [string trim $v(uid)]] }
    if {$uid eq "" || ![info exists ::ufv::ann($uid)]} { status "Fetch an accession first." 1; return }
    if {$path eq ""} {
        set path [tk_getSaveFile -parent .ufv -title "Export CSV" -defaultextension .csv \
            -initialfile "${uid}_residue_annotations[expr {$protvar ? {_protvar} : {}}].csv" -filetypes {{CSV .csv} {All *}}]
        if {$path eq ""} return
    }
    status [expr {$protvar ? "writing CSV with ProtVar columns (one request per residue - slow) ..." : "writing CSV ..."}]
    busy
    set S(work) 1
    if {[llength [info commands ::ufv::csv_async]]} {
        if {[catch { ::ufv::csv_async [expr {$m ne "" && [::ufv::st::has $m] ? $m : ""}] $path $protvar [list ::ufv::gui::_csv_path] } err]} {
            unset -nocomplain S(work); fail CSV $err
        }
        return
    }
    set a [dict create uid $uid out [file normalize $path]]
    if {$protvar} { dict set a protvar 1 }
    ::ufv::backend::async csv $a [list ::ufv::gui::_csv_done $path]
}
proc ::ufv::gui::_csv_done {path ok res} {
    variable S
    unset -nocomplain S(work)
    if {!$ok} { fail "CSV" $res; return }
    set S(lastcsv) [dict get $res path]
    status "CSV written: [dict get $res path]"
}

# Current cartoon colouring -> per-atom values for the B-factor column.
proc ::ufv::gui::bfactor_values {m} {
    set mode [::ufv::cartoon::mode $m]
    set a [ann [uid_of $m]]
    set vals [dict create]
    switch -- $mode {
        none - plddt - bfactor { return "" }
        alphamissense { if {$a ne ""} { set vals [dict get $a am] } }
        burden { if {$a ne ""} { foreach p [dict get $a burden] { dict set vals $p 1 } } }
        hotspots - hubs - pockets {
            set c [::ufv::st::get $m {cache tiers} ""]
            if {$c ne ""} {
                set rank {strong 3 moderate 2 weak 1 pocket 1 exposed 2}
                dict for {p t} [dict get [lindex $c 1] $mode] { dict set vals $p [expr {[dict exists $rank $t] ? [dict get $rank $t] : 1}] }
            }
        }
        default {
            # categorical: 1-based group index of the cart:<mode> layer (per colour)
            set gi 0
            foreach g [::ufv::st::get $m [list layers cart:$mode groups] {}] {
                incr gi
                set sel [lindex $g 1]
                ::ufv::with_sel x $m "($sel) and name CA" {
                    foreach {ch ri} [join [$x get {chain resid}]] {
                        set u [::ufv::map::r2u $m $ch $ri]
                        if {$u ne ""} { dict set vals $u $gi }
                    }
                }
            }
        }
    }
    return $vals
}
proc ::ufv::gui::save_bfactor {{path ""}} {
    set m [amol]
    if {$m eq ""} { status "Load a structure first." 1; return }
    if {$path eq ""} {
        set path [tk_getSaveFile -parent .ufv -title "Save B-factor PDB" -defaultextension .pdb \
            -initialfile "[regsub -all {[^A-Za-z0-9_.-]} [molinfo $m get name] _]_ufv.pdb" -filetypes {{PDB .pdb} {All *}}]
        if {$path eq ""} return
    }
    set mode [::ufv::cartoon::mode $m]
    if {[llength [info commands ufv_save_bfactor]]} {
        # features module: the extension's export.js values (AM / hotspots / hubs / burden / pocket / pLDDT / B;
        # otherwise 1 = residues of the visible marker layers)
        set bm [expr {$mode in {alphamissense hotspots hubs burden plddt bfactor} ? $mode : ($mode eq "pockets" ? "pocket" : "selected")}]
        if {[catch { ufv_save_bfactor $path $m $bm } err]} { fail "B-factor PDB" $err; return }
        status "saved $path (B-factor = $bm)."
        return $path
    }
    set vals [bfactor_values $m]
    if {[catch {
        ::ufv::with_sel all $m all {
            if {$vals ne ""} {
                set old [$all get beta]
                ::ufv::layers::write_field $m beta $vals 0.0 all
                set code [catch { $all writepdb $path } err]
                $all set beta $old
                if {$code} { error $err }
            } else {
                $all writepdb $path
            }
        }
    } err]} { fail "B-factor PDB" $err; return }
    status "saved $path (B-factor = [expr {$vals eq {} ? {original B-factor / pLDDT} : $mode}])."
    return $path
}
proc ::ufv::gui::clear {} {
    variable v
    foreach m [targets] { run clear [list ufv_clear $m] }
    set v(cart) None
    report_placeholder
    summary report "click a residue / row, or Pick 3D"
    sync_later
}

# Generic completion of an async layer op (cb ok result): errors to the status line, panel resynced.
proc ::ufv::gui::_async_done {what ok res} {
    if {!$ok} { fail $what $res }
    sync_later
}

proc ::ufv::gui::_marker_n {m lname blayer ok n} {
    if {!$ok} { fail "$lname layer" $n; sync_later; return }
    if {$n == 0} { status "mol $m: no $blayer positions map onto this structure." } else { status "mol $m: $blayer layer on ($n positions)." }
    sync_later
}

proc ::ufv::gui::_cart_tier_n {m kind ok n} {
    variable S
    if {!$ok} { fail $kind $n; sync_later; return }
    set S(cart_ms) [expr {[clock milliseconds] - $S(cart_t0)}]
    status "mol $m: $kind ($n residues)."
    sync_later
}

proc ::ufv::gui::_q_apply {} {
    variable v
    if {![alive] || ![string is double -strict $v(q)]} return
    if {[llength [info commands ufv_pocket_q]]} {
        # recolours molecules showing the constraint pockets from the cached candidates (no recompute)
        foreach m [targets] { catch { ufv_pocket_q $v(q) $m } }
        sync_later
    } elseif {$v(cart) eq "Constraint pocket"} {
        set_cartoon pockets
    }
}

proc ::ufv::gui::show_pocket {m resid} {
    if {$m ni [molinfo list]} return
    run "pocket" [list ufv_show_pocket [join $resid +] $m]
}

proc ::ufv::gui::motif {m spec} {
    variable S
    if {[info exists S(nobrowser)]} { set S(lasturl) "motif:$spec"; status "would search motif $spec"; return }
    if {![llength [info commands ::ufv::motif_async]]} { run motif [list ufv_motif_search $spec $m]; return }
    set pocket ""; set pos $spec
    regexp {^pocket=(.+)$} $spec -> pocket
    if {$pocket ne ""} { set pos "" }
    busy
    ::ufv::motif_async $m $pos $pocket [list ::ufv::gui::_motif_url]
}

proc ::ufv::gui::_motif_url {ok url} {
    if {!$ok} { fail "motif search" $url; return }
    if {$url eq ""} { status "motif search unavailable (needs >= 2 mapped residues on a PDB entry or AlphaFold model)." 1; return }
    open_url $url
}

proc ::ufv::gui::_on_report {pos text} {
    variable S
    if {![alive]} return
    # the structured report of this position is already shown (report_data came first)
    if {[info exists S(richpos)] && $S(richpos) eq $pos} return
    render_report $pos $text
}

proc ::ufv::gui::_on_report_data {m pos res} {
    variable S
    if {![alive]} return
    if {[info exists S(rep_t0)]} { set S(rep_ms) [expr {[clock milliseconds] - $S(rep_t0)}] }
    set S(curpos) $pos; set S(curmol) $m; set S(curlig) ""; set S(richpos) $pos
    render_rich $m $pos $res
}

# colour tag (foreground) for the report text
proc ::ufv::gui::ctag {hex {bold 0}} {
    variable S
    if {$hex eq ""} { return "" }
    set t fg[string trimleft $hex #][expr {$bold ? "b" : ""}]
    if {![info exists S(ctag,$t)]} {
        $S(report) tag configure $t -foreground $hex
        if {$bold} { $S(report) tag configure $t -font ufvBold }
        set S(ctag,$t) 1
    }
    return $t
}

proc ::ufv::gui::htag {hex} {
    variable S
    set t h2[string trimleft $hex #]
    if {![info exists S(htag,$t)]} {
        $S(report) tag configure $t -font ufvBold -foreground $hex -spacing1 6 -spacing3 1
        $S(report) tag raise $t
        set S(htag,$t) 1
    }
    return $t
}

proc ::ufv::gui::fmtnum {x} {
    if {$x eq "" || ![string is double -strict $x]} { return "-" }
    if {abs($x) >= 10} { return [format %.1f $x] }
    return [format %.2f $x]
}

# Structured report (::ufv::report_async / event report_data): title + analysis flags + coloured sections, rows with
# clickable positions / pockets, the per-substitution prediction table, links.
proc ::ufv::gui::render_rich {m pos res} {
    variable S; variable U; variable C
    set t $S(report)
    set uid [dg $res uid [uid_of $m]]
    rclear
    rput [dg $res title $pos] h1
    rput "  $uid" sub
    foreach f [dg $res flags] {
        if {[dg $f on]} { rput "  $U(bullet) [dg $f label]" [ctag [dg $f color] 1] } else { rput "  [dg $f label]" muted }
    }
    rput "\n"
    if {$m ne "" && $m in [molinfo list] && [::ufv::st::has $m]} {
        set where {}
        foreach ch [::ufv::map::chains $m] { set r [::ufv::map::u2r $m $ch $pos]; if {$r ne ""} { lappend where "$ch:$r" } }
        rput "mol $m [molinfo $m get name] $U(dot) [expr {[llength $where] ? "chain:resid [join $where {, }]" : "not modelled"}]\n" sub
    }
    foreach sec [dg $res sections] {
        rput "[dg $sec title]\n" [htag [dg $sec color $C(teal)]]
        foreach row [dg $sec rows] {
            set lab [dg $row label]; set val [dg $row value]; set col [dg $row color]
            if {[dg $row evidence 0]} { rput "    $lab  " small; rput "$val\n" small; continue }
            rput "  "
            if {[dg $row bold 0]} { rput $lab [ctag [expr {$col ne "" ? $col : $C(fg)}] 1] } else { rput $lab key }
            rput "  "
            set ps [dg $row positions]
            set toks [split $val ,]
            if {[llength $ps] && [llength $toks] == [llength $ps]} {
                set i 0
                foreach tok $toks p $ps {
                    if {$i} { rput "," }
                    regexp {^(\s*)(.*)$} $tok -> sp body
                    rput $sp
                    rlink $body [list ::ufv::gui::show_residue $p $m] [ctag $col]
                    incr i
                }
            } else {
                rput $val [ctag $col]
            }
            if {[dg $row pocket] ne ""} {
                set pid [dg $row pocket]
                foreach pk [dg $res pockets] {
                    if {[dg $pk pocketId] ne $pid} continue
                    rput "  "; rlink "show" [list ::ufv::gui::show_pocket $m [dg $pk resid]]
                    if {[dg $pk motif 0] && [llength [info commands ufv_motif_search]]} {
                        rput "  "; rlink "motif $U(ext)" [list ::ufv::gui::motif $m "pocket=$pid"]
                    }
                    if {[dg $pk composition] ne ""} { rput "\n    [dg $pk composition]" small }
                    if {[llength [dg $pk sites]]} { rput "\n    sites: [join [dg $pk sites] {; }]" small }
                }
            }
            rput "\n"
        }
        set tb [dg $sec table]
        if {$tb ne "" && [llength [dg $tb rows]]} {
            set cols [dg $tb columns]
            set tabs {}; set x 64
            foreach c [lrange $cols 1 end] { lappend tabs $x; incr x 52 }
            $t tag configure tbl -tabs $tabs -font ufvSmall -lmargin1 10 -lmargin2 10
            $t tag configure tblh -tabs $tabs -font ufvSmallB -foreground $C(grey) -lmargin1 10 -lmargin2 10
            rput "[join $cols \t]\n" tblh
            set obs [dg $tb observed]
            foreach r [dg $tb rows] {
                set sub [lindex $r 0]
                rput $sub [expr {$sub in $obs ? {tbl b pathogenic} : {tbl}}]
                foreach c [lrange $cols 1 end] x [lrange $r 1 end] {
                    set tg tbl
                    if {$c eq "AM" && [string is double -strict $x]} {
                        set tg [list tbl [expr {$x >= 0.564 ? "amhi" : ($x < 0.34 ? "amlo" : "ammid")}]]
                    }
                    rput "\t" tbl; rput [fmtnum $x] $tg
                }
                rput "\n" tbl
            }
            if {[dg $tb link] ne ""} {
                rput "  "; rlink "ProtVar $U(ext)" [list ::ufv::gui::open_url [dg $tb link]]; rput "  (bold = observed variant)\n" muted
            }
        }
    }
    if {![dg $res services 0] && [info exists S(rep_t0)]} {
        rput "\nfetching ProtVar predictors / pockets / PDBe-KB $U(ell)\n" muted
    }
    rput "\n"
    rlink "UniProt $U(ext)" [list ::ufv::gui::open_url "https://www.uniprot.org/uniprotkb/$uid/entry#sequences"]
    rput "   "
    rlink "AlphaFold $U(ext)" [list ::ufv::gui::open_url "https://alphafold.ebi.ac.uk/entry/$uid"]
    rput "   "
    rlink "ProtVar $U(ext)" [list ::ufv::gui::open_url "https://www.ebi.ac.uk/ProtVar/query?search=$uid+$pos"]
    rput "\n"
    rdone
    summary report "[dg $res title $pos] $U(dot) $uid"
}

proc ::ufv::gui::_ligand_timeout {key} {
    variable S
    if {![alive] || $S(ligpend) ne $key || $S(curlig) ne $key} return
    set S(ligpend) ""
    render_ligand {*}$key ""
    status "ligand data: no reply within 35 s (chemistry / evidence not shown)." 1
}
proc ::ufv::gui::_ligand_data {m ch ri rn ok res} {
    variable S
    if {$S(ligpend) eq [list $m $ch $ri $rn]} { set S(ligpend) ""; catch { after cancel $S(ligafter) } }
    if {!$ok} { fail "ligand data" $res; return }
    if {![alive] || $S(curlig) ne [list $m $ch $ri $rn] || $m ni [molinfo list]} return
    render_ligand $m $ch $ri $rn $res
}

proc ::ufv::gui::render_ligand {m ch ri rn data} {
    variable S; variable v; variable U; variable C
    set ls "chain [::ufv::map::chain_token $ch] and resid [::ufv::map::compress [list $ri]] and resname $rn"
    set d [expr {$data ne "" ? [dg $data dist $v(nearby)] : $v(nearby)}]
    rclear
    set info [expr {$data ne "" ? [dg $data info] : ""}]
    rput $rn h1; rput "  chain $ch $U(dot) resid $ri $U(dot) mol $m\n" sub
    if {$info ne "" && [dg $info name] ne ""} { rput "[dg $info name]\n" b }
    set nat 0
    catch { ::ufv::with_sel x $m $ls { set nat [$x num] } }
    rput "Ligand\n" h2
    rput "  atoms " key; rput "$nat"
    if {$data ne ""} {
        set cps [dg $data copies]
        rput "   copies " key; rput "[llength $cps]\n"
        if {[llength $cps] > 1} {
            rput "  "
            set i 0
            foreach c $cps {
                if {$i >= 16} { rput " $U(ell)" muted; break }
                lassign $c cch cri
                if {$cch eq $ch && $cri == $ri} { rput "$cch:$cri" b } else { rlink "$cch:$cri" [list ::ufv::gui::show_ligand $m $cch $cri $rn] }
                rput "  "; incr i
            }
            rput "\n"
        }
    } else { rput "\n" }
    if {$info ne ""} {
        foreach {k lab} {formula Formula weight "Mol. weight" hbondDonors "H-bond donors" hbondAcceptors "H-bond acceptors" inchikey InChIKey} {
            if {[dg $info $k] ne ""} { rput "  $lab " key; rput "[dg $info $k]\n" }
        }
        if {[llength [dg $info similar]]} {
            rput "  Similar in structure " key
            set tx {}; foreach s [dg $info similar] { lappend tx "[lindex $s 0] [format %.2f [lindex $s 1]]" }
            rput "[join $tx {, }]\n"
        }
    }
    if {$data ne ""} {
        set aff [dg $data aff]
        if {$aff ne ""} {
            rput "AlphaFill transplant\n" h2
            rput "  donor " key; rput "[dg $aff donorPdb]"
            rput "   identity " key; rput "[dg $aff identity]"
            rput "   clash " key; rput "[dg $aff clash][expr {[dg $aff band] ne {} ? " ([dg $aff band])" : {}}]\n"
        }
        set cl [dg $data clash]
        if {$cl ne ""} {
            set band [dg $cl band]
            rput "  closest protein contact " key
            rput [format "%.2f %s (%s)\n" [dg $cl minDist] $U(ang) $band] [expr {$band eq "severe" ? "err" : ($band eq "suspicious" ? "deleterious" : "")}]
        }
        set ev [dg $data evidence]
        if {$ev ne ""} {
            rput "Pocket evidence\n" h2
            if {[dg $ev meanB] ne ""} {
                rput [expr {[dg $data experimental] ? "  mean B-factor " : "  mean pLDDT "}] key
                rput [format %.1f [dg $ev meanB]] b
                rput [expr {[dg $data experimental] || [dg $ev band] eq "" ? "\n" : " ([dg $ev band])\n"}]
            }
            if {[llength [dg $ev sites]]} {
                rput "  functional sites " key
                set i 0
                foreach s [dg $ev sites] {
                    if {$i} { rput "; " }
                    rput "[dg $s desc] "; rlink [dg $s pos] [list ::ufv::gui::show_residue [dg $s pos] $m]
                    incr i
                }
                rput "\n"
            }
            if {[dg $ev compositionText] ne ""} { rput "  residue character " key; rput "[dg $ev compositionText]\n" }
        }
    }
    # pocket residues (UniProt numbering where mapped)
    rput "Pocket residues ($U(le) $d $U(ang))\n" h2
    set pocket {}
    catch {
        ::ufv::with_sel x $m "protein and name CA and same residue as (protein and within $d of ($ls))" {
            foreach {c r rnn} [join [$x get {chain resid resname}]] { lappend pocket [list $c $r $rnn] }
        }
    }
    if {![llength $pocket]} { rput "  none\n" muted } else {
        rput "  "
        foreach p $pocket {
            lassign $p c r rnn
            set u [::ufv::map::r2u $m $c $r]
            set lbl [expr {$u ne {} ? "[aa1 $rnn]$u" : "$c:[aa1 $rnn]$r"}]
            if {$u ne ""} { rlink $lbl [list ::ufv::gui::show_residue $u $m] } else { rput $lbl muted }
            rput "  "
        }
        rput "\n"
    }
    if {$data ne "" && [dg $data motif 0] && [llength [dg $data near]] >= 2 && [llength [info commands ufv_motif_search]]} {
        rput "  "; rlink "RCSB motif search $U(ext)" [list ::ufv::gui::motif $m [join [dg $data near] +]]; rput "\n"
    }
    rput "\n"
    set links {}
    if {$info ne ""} { set links [dg $info links] }
    if {![llength $links]} {
        set links [list [list RCSB "https://www.rcsb.org/ligand/$rn"] [list PDBe "https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/$rn"]]
    }
    foreach l $links { rlink "[lindex $l 0] $U(ext)" [list ::ufv::gui::open_url [lindex $l 1]]; rput "   " }
    rput "\n"
    if {$data eq "" && [info exists S(ligpend)] && $S(ligpend) eq [list $m $ch $ri $rn]} { rput "loading ligand data $U(ell)\n" muted }
    rdone
    summary report "$rn $ch:$ri (ligand)"
}

proc ::ufv::gui::_csv_path {ok res} {
    variable S
    unset -nocomplain S(work)
    if {!$ok} { fail "CSV" $res; return }
    set S(lastcsv) $res
    status "CSV written: $res"
}

# ---- AlphaFill donor-identity filter (shown for AlphaFill models)
proc ::ufv::gui::aff_select {} {
    variable v
    if {![regexp {(\d+)} $v(aff) -> pct]} return
    foreach m [targets] { run "AlphaFill identity" [list ufv_ligand_identity $pct $m] }
    sync_later
}

# ---- unmapped (chimeric / tag) residues of our chains
proc ::ufv::gui::chimera_toggle {} {
    set m [amol]
    if {$m eq ""} return
    if {[::ufv::cartoon::mode $m] eq "chimera"} { ::ufv::cartoon::plain $m } else { run chimera [list ufv_chimera $m] }
    sync_later
}
