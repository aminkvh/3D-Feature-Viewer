# VMD / Tcl package index: lets `package require ufv` load the plugin when this folder is on auto_path
# (e.g. copied into VMD's plugins/noarch/tcl/ufv2.1).
package ifneeded ufv 2.1 [list source [file join $dir ufv_vmd.tcl]]
