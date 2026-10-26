# =============================================================================================
# 3D Feature Viewer - VMD plugin (loader)
# =============================================================================================
# Projects UniProt residue annotations (PTMs, disease variants, AlphaMissense, functional sites, domains,
# membrane topology, mutagenesis, mutation burden, structure analyses) onto structures and trajectories in
# VMD, or downloads AlphaFold / experimental / computed models directly. Data come from ufv_core.py (the
# backend shared with the PyMOL plugin), run once as a persistent Python process (`ufv_core.py --serve`).
#
# Install:  source <folder>/vmd/ufv_vmd.tcl      (ufv_core.py next to this file or one folder up)
#           or put the vmd/ folder in a VMD plugin directory (pkgIndex.tcl: package require ufv)
# Python:   auto-detected (py -3, python3, python; >= 3.8). Override: ufv_python <path>.
#
# Modules (sourced in this order; each file has one owner):
#   ufv_backend.tcl   backend process + JSON requests + async/sync calls            (foundation)
#   ufv_palette.tcl   fixed colour ids 17-32 + nearest-colour mapping              (foundation)
#   ufv_layers.tcl    per-molecule state, numbering maps, layer + cartoon model    (foundation)
#   ufv_commands.tcl  the ufv_* commands ported onto the foundation                (foundation -> features)
#   ufv_features.tcl  new annotation features                                       (features agent)
#   ufv_traj.tcl      trajectory analyses                                           (trajectory agent)
#   ufv_gui.tcl       Tk panel                                                      (GUI agent)
#
# Startup cost: when this file is sourced by the user's vmd.rc (the install snippet sets ::ufv_dir), only a
# stub is set up (~1 ms instead of ~15-30 ms): the ufv_* commands and the Extensions entry exist, and the
# first call of any of them sources the modules (::ufv::load_all). `source ufv_vmd.tcl` from a script or the
# console, and `package require ufv`, load everything at once. Override: env UFV_LAZY=1 (stub) / 0 (all).
# =============================================================================================

namespace eval ::ufv {
    variable version 2.1
    # (file normalize costs ~1 ms per call on Windows: only for a relative path here, else in load_all)
    variable dir [file dirname [expr {[file pathtype [info script]] eq "absolute" ? [file join [info script]] : [file normalize [info script]]}]]
    variable modules {ufv_backend ufv_palette ufv_layers ufv_commands ufv_features ufv_traj ufv_gui}
    # global commands that exist before the modules are loaded (stubs; the first call loads everything)
    variable stubs {ufv_align ufv_alphamissense ufv_backend ufv_bfactor ufv_binding ufv_burden ufv_cartoon
        ufv_chain ufv_chainmap ufv_chimera ufv_clear ufv_contacthubs ufv_core ufv_csv ufv_custom ufv_diseases
        ufv_domains ufv_fetch ufv_focus ufv_gui ufv_hide ufv_hotspots ufv_info ufv_layers ufv_ligand_focus
        ufv_ligand_identity ufv_ligand_info ufv_ligands ufv_ligands_hide ufv_load ufv_load_annotations ufv_load_structure ufv_map
        ufv_motif_search ufv_mutagenesis ufv_nearby ufv_open_uniprot ufv_partner_layers ufv_partners ufv_pdbekb
        ufv_plddt ufv_pocket_q ufv_pockets ufv_pockets_all ufv_prox ufv_ptms ufv_python ufv_report
        ufv_reset_annotations ufv_resetview ufv_residue ufv_save_bfactor ufv_select ufv_show ufv_show_annotations ufv_show_pocket ufv_sites ufv_structures
        ufv_topology ufv_tractability ufv_traj ufv_traj_annot ufv_traj_bfactor ufv_traj_cancel
        ufv_traj_contacts ufv_traj_csv ufv_traj_dist ufv_traj_engine ufv_traj_frames ufv_traj_gui
        ufv_traj_nearby ufv_traj_pockets ufv_traj_rmsf ufv_traj_show ufv_traj_tiers ufv_unload_annotations ufv_use ufv_variants}
}

# Source every module once (idempotent). Returns 1 when this call loaded them.
proc ::ufv::load_all {} {
    variable loaded
    if {[info exists loaded]} { return 0 }
    set loaded 1
    set t0 [clock microseconds]
    set ::ufv::dir [file normalize $::ufv::dir]
    foreach m $::ufv::modules { uplevel #0 [list source [file join $::ufv::dir $m.tcl]] }
    # Colour slots 17-32 are defined once, here.
    if {[llength [info commands color]]} { ::ufv::palette::init }
    # Fallback panel until ufv_gui.tcl provides one (VMD's extension menu needs a proc returning a toplevel).
    if {![llength [info procs ::ufv_gui]] || [string match "*::ufv::_stub*" [info body ::ufv_gui]]} {
        proc ::ufv_gui {} {
            if {![llength [info commands toplevel]]} { error "ufv_gui needs Tk (VMD started without graphics?)" }
            if {[winfo exists .ufv]} { wm deiconify .ufv; raise .ufv; return .ufv }
            toplevel .ufv
            wm title .ufv "3D Feature Viewer"
            label .ufv.l -justify left -padx 12 -pady 12 -text \
                "The panel is being rebuilt.\nUse the console commands: ufv_fetch, ufv_load, ufv_ptms, ufv_variants, ..."
            pack .ufv.l
            return .ufv
        }
    }
    variable load_ms [expr {([clock microseconds] - $t0) / 1000.0}]
    return 1
}

# Stub for a global command: load the modules, then run the real command (which the modules defined).
proc ::ufv::_stub {name args} {
    ::ufv::load_all
    if {[string match "*::ufv::_stub*" [info body ::$name]]} { error "$name is not available in this plugin version" }
    return [uplevel 1 [list ::$name {*}$args]]
}

package provide ufv $::ufv::version

if {[info exists ::env(UFV_LAZY)] && $::env(UFV_LAZY) ne ""} {
    set ::ufv::lazy [string is true -strict $::env(UFV_LAZY)]
} else {
    # the vmd.rc install snippet sets ufv_dir to the folder holding vmd/ufv_vmd.tcl
    set ::ufv::lazy [expr {[info exists ::ufv_dir] && [info level] == 0 &&
        [string equal -nocase [file join $::ufv_dir vmd] $::ufv::dir]}]
}
if {$::ufv::lazy && ![info exists ::ufv::loaded]} {
    foreach ::ufv::_c $::ufv::stubs {
        if {![llength [info commands ::$::ufv::_c]]} { proc ::$::ufv::_c args "tailcall ::ufv::_stub $::ufv::_c {*}\$args" }
    }
    unset ::ufv::_c
} else {
    ::ufv::load_all
}

if {[llength [info commands vmd_install_extension]] && [llength [info commands toplevel]]} {
    catch { vmd_install_extension ufv ufv_gui "Analysis/3D Feature Viewer" }
}

puts "\[UFV\] 3D Feature Viewer $::ufv::version (VMD) [expr {$::ufv::lazy && ![info exists ::ufv::loaded] ? {ready} : {loaded}}] - type 'ufv_gui' for the panel."
puts "      ufv_fetch <uid> | ufv_load <uid> | ufv_structures | ufv_load_structure <key> | ufv_use ?molid? ?uid?"
puts "      ufv_map identity|sifts <pdb>|ranges|auto | ufv_chain | ufv_chainmap | ufv_custom"
puts "      ufv_ptms ufv_variants ufv_sites ufv_mutagenesis ufv_ligands | ufv_domains ufv_topology ufv_alphamissense"
puts "      ufv_burden ufv_plddt ufv_bfactor ufv_hotspots ufv_contacthubs ufv_pockets | ufv_cartoon plain"
puts "      ufv_hide/ufv_show <layer> | ufv_layers | ufv_clear | ufv_focus <pos> | ufv_residue <pos> | ufv_csv"
puts "      ufv_python <path> | ufv_backend ?ping|restart?"
