# =============================================================================================
# ufv_features.tcl - OWNER: the FEATURES agent. New annotation features on the foundation: diseases, chimera /
# partners / partner layers, ProtVar pockets (all pockets, show, motif search), PDBe-KB binding, Open Targets
# tractability, ligand panel data (AlphaFill identity / clash, pocket evidence, chemistry), info, B-factor export,
# selections.
#
# Sourced by ufv_vmd.tcl after ufv_backend / ufv_palette / ufv_layers / ufv_commands.
# Every network step is an async operation (... cb, cb called as `ok result`); the console commands wait for it
# with ::ufv::await. Layers: one rep per colour (::ufv::layers), never per annotation. See API.md "Features API".
# =============================================================================================

namespace eval ::ufv::feat {
    if {![array exists kb]} { array set kb {} }
}
if {![array exists ::ufv::allpockets]} { array set ::ufv::allpockets {} }

set ::ufv::chimera_color #8fa5c0
set ::ufv::kb_colors [dict create ligand #2e7d32 interface #1565c0 both #00838f]

# Small join helper (all_of): run several async steps, then {*}$done results (dict name -> {ok result}).
# Steps still pending after timeout_ms are reported as {0 "timed out"} (a late reply is dropped).
proc ::ufv::feat::all_of {steps done {timeout_ms 30000}} {
    set id [incr ::ufv::feat::seq]
    set ::ufv::feat::J($id) [dict create left [expr {[llength $steps] / 2}] res [dict create] done $done         names [dict keys $steps] aid [after $timeout_ms [list ::ufv::feat::_join_timeout $id]]]
    foreach {name cmd} $steps {
        if {[catch { uplevel #0 [linsert $cmd end [list ::ufv::feat::_joined $id $name]] } err]} {
            ::ufv::feat::_joined $id $name 0 $err
        }
    }
    if {![llength $steps]} { after 0 [list ::ufv::feat::_joined $id "" 1 ""] }
}
proc ::ufv::feat::_join_timeout {id} {
    if {![info exists ::ufv::feat::J($id)]} return
    foreach n [dict get $::ufv::feat::J($id) names] {
        if {![dict exists $::ufv::feat::J($id) res $n]} { ::ufv::feat::_joined $id $n 0 "timed out" }
    }
}
proc ::ufv::feat::_joined {id name ok res} {
    if {![info exists ::ufv::feat::J($id)]} return
    if {$name ne ""} {
        if {[dict exists $::ufv::feat::J($id) res $name]} return
        dict set ::ufv::feat::J($id) res $name [list $ok $res]
        dict incr ::ufv::feat::J($id) left -1
    }
    if {[dict get $::ufv::feat::J($id) left] > 0} return
    set j $::ufv::feat::J($id)
    unset ::ufv::feat::J($id)
    after cancel [dict get $j aid]
    {*}[dict get $j done] [dict get $j res]
}

# Categorical cartoon mode from explicit selections (residues without a UniProt position, e.g. chimera):
# groups {{colour-id seltext} ..}; same state handling as ::ufv::cartoon::categorical.
proc ::ufv::feat::cartoon_sel {mol mode groups {rebuild ""}} {
    ::ufv::layers::batch {
        set b [::ufv::cartoon::base_idx $mol]
        ::ufv::cartoon::_hide_modes $mol $mode
        set covered {}
        foreach g $groups { lappend covered "([lindex $g 1])" }
        ::ufv::layers::define $mol "cart:$mode" $groups -style $::ufv::cartoon::style -kind cartoon -rebuild $rebuild
        set bsel [expr {[llength $covered] ? "protein and not ([join $covered { or }])" : "protein"}]
        mol modselect $b $mol $bsel
        mol modcolor $b $mol ColorID [::ufv::palette::base]
        mol showrep $mol $b 1
    }
    ::ufv::st::put $mol cartoon [dict create mode $mode kind categorical basesel $bsel basehex ""]
    ::ufv::st::put $mol [list modes $mode] [dict create basesel $bsel basehex ""]
    ::ufv::emit cartoon $mol $mode
    return [llength $groups]
}
# drop a cartoon mode (layer + remembered base); plain cartoon if it was the one shown
proc ::ufv::feat::cartoon_drop {mol mode} {
    if {[::ufv::cartoon::mode $mol] eq $mode} { ::ufv::cartoon::plain $mol }
    if {[::ufv::layers::exists $mol cart:$mode]} { ::ufv::layers::delete $mol cart:$mode }
    ::ufv::st::remove $mol [list modes $mode]
}

# ---------------------------------------------------------------------------------------------
# Info / diseases
# ---------------------------------------------------------------------------------------------
# ufv_info ?uid?  - summary of the cached annotations (returns a dict)
proc ufv_info {{uid ""}} {
    set uid [::ufv::need $uid]
    set a [::ufv::annotations $uid]
    set npath 0
    foreach v [dict get $a variants] { if {[string match -nocase *pathogenic* [dict get $v consequence]]} { incr npath } }
    set d [dict create uid $uid seqlen [dict get $a seqlen] ptms [llength [dict get $a ptms]] \
        variants [llength [dict get $a variants]] sites [llength [dict get $a sites]] domains [llength [dict get $a domains]] \
        topology [llength [dict get $a topology]] mutagenesis [llength [dict get $a mutagenesis]] pathogenic $npath]
    puts [format "\[UFV\] %s  seq=%daa  PTMs=%d  variants=%d  sites=%d  domains=%d  topology=%d" $uid \
        [dict get $a seqlen] [dict get $d ptms] [dict get $d variants] [dict get $d sites] [dict get $d domains] [dict get $d topology]]
    puts "\[UFV\] pathogenic/likely-pathogenic variants: $npath"
    if {[dict get $a protnlm] ne ""} { puts "\[UFV\] [dict get $a protnlm]" }
    return $d
}

# ::ufv::diseases_async uid cb -> {{name count color} ..} in first-appearance order (+ Unclassified)
proc ::ufv::diseases_async {uid cb} {
    ::ufv::backend::async diseases [dict create uid [string toupper $uid]] $cb
}
# ufv_diseases ?uid?  - the disease summary used by the variant colour profile / filter
proc ufv_diseases {{uid ""}} {
    set uid [::ufv::need $uid]
    set l [::ufv::await [list ::ufv::diseases_async $uid]]
    foreach d $l {
        puts [format "  %-48s %5d  %s" [string range [dict get $d name] 0 47] [dict get $d count] [dict get $d color]]
    }
    return $l
}

# ---------------------------------------------------------------------------------------------
# Chimera / partners
# ---------------------------------------------------------------------------------------------
# ufv_chimera ?molid?  - residues of OUR chains without a UniProt position (chimeric / fusion segments, tags,
# linkers) in blue-grey, as a cartoon colouring (cart:chimera). Returns the number of residues.
proc ufv_chimera {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set m [::ufv::map::get $mol]
    set uid [::ufv::st::get $mol uid]
    if {[dict get $m mode] eq "identity"} {
        ::ufv::feat::cartoon_drop $mol chimera
        ::ufv::status "mol $mol uses identity numbering: every residue maps to UniProt."
        return 0
    }
    set chains [::ufv::map::chains $mol]
    set keys {}
    ::ufv::with_sel s $mol "$::ufv::poly_base and name CA" {
        set seen [dict create]
        foreach {ch ri ins} [join [$s get {chain resid insertion}]] {
            if {$ch ni $chains || [dict exists $seen $ch|$ri]} continue
            dict set seen $ch|$ri 1
            if {[string trim $ins] ne "" || [::ufv::map::r2u $mol $ch $ri] eq ""} { lappend keys $ch $ri }
        }
    }
    if {![llength $keys]} {
        ::ufv::feat::cartoon_drop $mol chimera
        ::ufv::status "mol $mol: every modelled residue of chain(s) [join $chains ,] maps to $uid."
        return 0
    }
    ::ufv::feat::cartoon_sel $mol chimera [list [list [::ufv::palette::id $::ufv::chimera_color] "protein and [::ufv::feat::keysel $keys]"]] \
        [list ufv_chimera $mol]
    set n [expr {[llength $keys] / 2}]
    set chs {}; foreach {ch ri} $keys { if {$ch ni $chs} { lappend chs $ch } }
    set orgs [expr {[dict exists $m meta organisms] ? [dict get $m meta organisms] : {}}]
    set chim [expr {[dict exists $m meta chimeric] && [dict get $m meta chimeric] eq "1"}]
    ::ufv::status "mol $mol: $n residue(s) of chain(s) [join [lsort $chs] ,] have no $uid mapping (blue-grey)[expr {$chim && [llength $orgs] ? " - chimeric: [join $orgs { / }]" : ""}]."
    return $n
}

# ::ufv::partners_async mol cb -> {partners {{accession chainId uniprotStart uniprotEnd name}..} custom {chain acc}
#                                   chimeric 0|1 chimericChains organisms pdbId}
proc ::ufv::partners_async {mol cb} {
    set a [::ufv::feat::ctx $mol 0]
    ::ufv::backend::async partners $a [list ::ufv::feat::_partners_got $mol $cb]
}
proc ::ufv::feat::_partners_got {mol cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    set m [::ufv::map::get $mol]
    dict set res custom [expr {[dict exists $m partners] && [dict get $m mode] eq "custom" ? [dict get $m partners] : {}}]
    ::ufv::st::put $mol partnerinfo $res
    {*}$cb 1 $res
}
# ufv_partners ?molid?  - other proteins of the loaded PDB entry (accession -> chains, UniProt range) + chains
# assigned to another accession with ufv_chainmap
proc ufv_partners {{molid ""}} {
    set mol [::ufv::resolve $molid]
    set r [::ufv::await [list ::ufv::partners_async $mol]]
    set by [dict create]
    foreach p [dict get $r partners] { dict lappend by [dict get $p accession] $p }
    if {![dict size $by] && ![dict size [dict get $r custom]]} {
        ::ufv::status "mol $mol: no partner chains (other UniProt entries) known."
        return $r
    }
    puts "\[UFV\] mol $mol partner chains[expr {[dict size $by] ? " ([dict size $by] accessions)" : ""}]:"
    dict for {acc ms} $by {
        set chs {}; set rng {}
        foreach p $ms {
            if {[dict get $p chainId] ni $chs} { lappend chs [dict get $p chainId] }
            if {[dict get $p uniprotStart] ne ""} {
                set t "[dict get $p uniprotStart]-[dict get $p uniprotEnd]"
                if {$t ni $rng} { lappend rng $t }
            }
        }
        set nm [dict get [lindex $ms 0] name]
        puts [format "    %-12s chains %-10s %s%s" $acc [join [lsort $chs] ,] \
            [expr {[llength $rng] ? "UniProt [join [lsort $rng] {, }]" : ""}] [expr {$nm ne "" ? "  $nm" : ""}]]
    }
    dict for {ch acc} [dict get $r custom] { puts [format "    %-12s chain  %-10s (assigned with ufv_chainmap)" $acc $ch] }
    return $r
}

# ::ufv::partner_layers_async mol layers cb  (layers = all | off | any of variants ptms sites) -> cb ok nspheres
# One sphere per partner-protein residue (Ca), coloured by the top enabled layer (pathogenic variants > PTMs >
# sites): ONE layer feat:partner with one rep per colour.
proc ::ufv::partner_layers_async {mol layers cb} {
    set a [::ufv::feat::ctx $mol 0]
    dict set a layers $layers
    ::ufv::status "fetching partner-protein annotations ..."
    ::ufv::backend::async partner_layers $a [list ::ufv::feat::_partner_got $mol $layers $cb]
}
proc ::ufv::feat::_partner_got {mol layers cb ok res} {
    if {!$ok} { ::ufv::status "partner layers failed: $res"; {*}$cb 0 $res; return }
    if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; return }
    # = check exactly these partner rows (prt:<chain>:<resi>); the colours are this layer set's. With every layer,
    # the spheres are also the Partners rows of the panel's list.
    ::ufv::st::put $mol partner_spheres [dict get $res spheres]
    if {[llength [dict get $res layers]] == 3} { ::ufv::st::put $mol partner_rows [dict get $res spheres] }
    ::ufv::annot::put $mol prt [lmap s [dict get $res spheres] { format prt:%s:%s [lindex $s 0] [lindex $s 1] }]
    if {![llength [dict get $res spheres]]} {
        ::ufv::status [expr {[dict get $res accessions] ? "mol $mol: partner layer off." : "mol $mol has no partner proteins with a UniProt mapping."}]
        {*}$cb 1 0
        return
    }
    set n [llength [dict get $res spheres]]
    ::ufv::status "mol $mol: partner layer on ([join [dict get $res layers] {, }]): $n residue(s) on [dict get $res chains] partner chain(s)."
    {*}$cb 1 $n
}
proc ::ufv::feat::partner_rebuild {mol layers} {
    # partner chains are addressed by chain + resid (not UniProt): a renumbering of OUR chains does not move them
    ::ufv::annot::render $mol prt
}
# ufv_partner_layers ?layers? ?molid?
proc ufv_partner_layers {{layers all} {molid ""}} {
    if {[string is integer -strict $layers] && $molid eq ""} { set molid $layers; set layers all }
    set mol [::ufv::resolve $molid]
    return [::ufv::await [list ::ufv::partner_layers_async $mol $layers]]
}

# ---------------------------------------------------------------------------------------------
# ProtVar pockets: every pocket of the protein (crawl), show one, RCSB motif search
# ---------------------------------------------------------------------------------------------
# A ProtVar pocket already fetched (this session's crawl, else the backend's residue-report cache) -> dict | ""
proc ::ufv::feat::find_pocket {uid id} {
    if {[info exists ::ufv::allpockets($uid)]} {
        foreach p $::ufv::allpockets($uid) { if {[dict get $p pocketId] eq $id} { return $p } }
    }
    set r [::ufv::backend::call pocket_find [dict create uid $uid id $id]]
    return [expr {[dict size $r] ? $r : ""}]
}

# ::ufv::pockets_all_async mol cb  (mol "" = no colouring) -> cb ok {count pockets {{pocketId score meanPlddt
#   radGyration energyPerVol buriedness resid color} ..}}; progress: status + event `progress pockets_all:<UID> done total`
proc ::ufv::pockets_all_async {mol cb} {
    set uid [expr {$mol ne "" ? [::ufv::mol_uid $mol] : [::ufv::need]}]
    ::ufv::status "scanning ProtVar pockets of $uid (one request per residue the first time) ..."
    ::ufv::feat::poll pockets_all:$uid "ProtVar pockets"
    ::ufv::backend::async pockets_all [dict create uid $uid] [list ::ufv::feat::_pa_done $mol $uid $cb]
}
proc ::ufv::feat::_pa_done {mol uid cb ok res} {
    poll_stop pockets_all:$uid
    if {!$ok} { ::ufv::status "pocket scan failed: $res"; {*}$cb 0 $res; return }
    set ::ufv::allpockets($uid) [dict get $res pockets]
    if {$mol ne "" && $mol in [molinfo list] && [dict get $res count]} { allpockets_apply $mol }
    {*}$cb 1 $res
}
# Colour every pocket's lining residues (palette by pocket rank, best pocket wins overlaps): cart:allpockets.
proc ::ufv::feat::allpockets_apply {mol} {
    set uid [::ufv::mol_uid $mol]
    if {![info exists ::ufv::allpockets($uid)]} { cartoon_drop $mol allpockets; return 0 }
    set taken [dict create]; set hexgroups [dict create]
    foreach p $::ufv::allpockets($uid) {
        foreach r [dict get $p resid] {
            if {[dict exists $taken $r]} continue
            dict set taken $r 1
            dict lappend hexgroups [dict get $p color] $r
        }
    }
    ::ufv::cartoon::categorical $mol allpockets $hexgroups "" [list ::ufv::feat::allpockets_apply $mol]
    return [dict size $taken]
}
proc ::ufv::feat::fmt {v f} { return [expr {$v eq "" ? "-" : [format $f $v]}] }
# ufv_pockets_all ?molid?
proc ufv_pockets_all {{molid ""}} {
    set mol [expr {[llength [molinfo list]] ? [::ufv::resolve $molid] : ""}]
    set res [::ufv::await [list ::ufv::pockets_all_async $mol]]
    set uid [expr {$mol ne "" ? [::ufv::mol_uid $mol] : $::ufv::uid}]
    if {![dict get $res count]} { ::ufv::status "$uid: no ProtVar pockets detected across the protein."; return $res }
    foreach p [dict get $res pockets] {
        puts [format "  pocket %-6s %s  score %s  res %d  pLDDT %s  Rg %s" [dict get $p pocketId] [dict get $p color] \
            [::ufv::feat::fmt [dict get $p score] %.0f] [llength [dict get $p resid]] \
            [::ufv::feat::fmt [dict get $p meanPlddt] %.0f] [::ufv::feat::fmt [dict get $p radGyration] %.1f]]
    }
    ::ufv::status "$uid: [dict get $res count] ProtVar pocket(s)[expr {$mol ne "" ? " coloured on mol $mol" : ""}]."
    return $res
}

# ::ufv::motif_async mol positions pocket cb -> cb ok url ("" when unavailable: needs >= 2 mapped residues on a
# PDB entry or an AlphaFold / AlphaFill model)
proc ::ufv::motif_async {mol positions pocket cb} {
    set a [::ufv::feat::ctx $mol 1]
    dict set a positions [string map {, + " " +} [string trim $positions]]
    if {$pocket ne ""} { dict set a pocket $pocket }
    ::ufv::backend::async motif_url $a [list apply {{cb ok r} { {*}$cb $ok [expr {$ok ? [dict get $r url] : $r}] }} $cb]
}
# ufv_motif_search <positions | pocket=<id>> ?molid? ?browser?   - RCSB Structure Motif Search URL (opened)
proc ufv_motif_search {spec {molid ""} {browser 1}} {
    set mol [::ufv::resolve $molid]
    set pocket ""; set pos $spec
    if {[regexp {^pocket=(.+)$} $spec -> pocket]} { set pos "" }
    set url [::ufv::await [list ::ufv::motif_async $mol $pos $pocket]]
    if {$url eq ""} {
        ::ufv::status "motif search unavailable (needs >= 2 mapped residues on a PDB entry or AlphaFold/AlphaFill model)."
        return ""
    }
    if {$browser} { catch { ::ufv::open_url $url } }
    ::ufv::status "RCSB strucmotif: $url"
    return $url
}

# ---------------------------------------------------------------------------------------------
# PDBe-KB (experimental binding across every PDB entry), Open Targets tractability
# ---------------------------------------------------------------------------------------------
# ::ufv::pdbekb_async mol kind cb  (kind ligand | interface | both) -> cb ok npositions; layer feat:kb
proc ::ufv::pdbekb_async {mol kind cb} {
    set uid [::ufv::mol_uid $mol]
    if {[info exists ::ufv::feat::kb($uid)]} { after 0 [list ::ufv::feat::_kb_got $mol $kind $cb 1 $::ufv::feat::kb($uid)]; return }
    ::ufv::backend::async pdbekb [dict create uid $uid] [list ::ufv::feat::_kb_got $mol $kind $cb]
}
proc ::ufv::feat::_kb_got {mol kind cb ok res} {
    if {!$ok} { {*}$cb 0 $res; return }
    if {$mol ni [molinfo list]} { {*}$cb 0 "molecule $mol was deleted"; return }
    set ::ufv::feat::kb([::ufv::mol_uid $mol]) $res
    {*}$cb 1 [kb_apply $mol $kind]
}
proc ::ufv::feat::kb_apply {mol kind} {
    set uid [::ufv::mol_uid $mol]
    if {![info exists ::ufv::feat::kb($uid)]} { set ::ufv::feat::kb($uid) [::ufv::backend::call pdbekb [dict create uid $uid]] }
    set r $::ufv::feat::kb($uid)
    set kind [string tolower $kind]
    set lig [expr {$kind in {both ligand ligands} ? [dict get $r ligand] : {}}]
    set ifc [expr {$kind in {both interface interfaces ppi} ? [dict get $r interface] : {}}]
    # = check the PDBe-KB rows of this kind (coloured by the kind, as before); the layer is drawn from the rows
    ::ufv::st::put $mol [list annkb $uid] $kind
    set n [::ufv::annot::put $mol kb [lmap p [lsort -integer -unique [concat $lig $ifc]] { format kb:%s $p }]]
    ::ufv::status "mol $mol: PDBe-KB $kind layer - $n residue(s) ([llength $lig] ligand-site, [llength $ifc] interface)."
    return $n
}
# ufv_pdbekb ?ligand|interface|both? ?molid?
proc ufv_pdbekb {{kind both} {molid ""}} {
    if {[string is integer -strict $kind] && $molid eq ""} { set molid $kind; set kind both }
    set mol [::ufv::resolve $molid]
    return [::ufv::await [list ::ufv::pdbekb_async $mol $kind]]
}

# ::ufv::binding_async uid pos cb -> {ligands {{code name pdb pdbCount ..}..} nLigands interfaces nInterfaces url}
proc ::ufv::binding_async {uid pos cb} {
    ::ufv::backend::async binding [dict create uid [string toupper $uid] pos $pos] $cb
}
# ufv_binding <pos> ?uid?
proc ufv_binding {pos {uid ""}} {
    set uid [::ufv::need $uid]
    set r [::ufv::await [list ::ufv::binding_async $uid $pos]]
    puts "\[UFV\] $uid $pos - PDBe-KB experimental binding ([dict get $r url])"
    set lt [dict get $r ligands]; set it [dict get $r interfaces]
    if {![llength $lt] && ![llength $it]} { puts "  no ligand / interface evidence at this residue." }
    if {[llength $lt]} {
        set t {}
        foreach l $lt {
            lappend t "[dict get $l code] ([dict get $l pdbCount] PDB[expr {[dict get $l pdb] ne "" ? ", e.g. [string toupper [dict get $l pdb]]" : ""}])"
        }
        set nl [dict get $r nLigands]
        puts "  Ligands: [join $t {, }][expr {$nl > 8 ? " +[expr {$nl - 8}]" : ""}]"
    }
    if {[llength $it]} {
        set t {}
        foreach p $it {
            set nm [dict get $p name]; if {$nm eq ""} { set nm [dict get $p accession] }
            lappend t "$nm ([dict get $p pdbCount])"
        }
        set ni [dict get $r nInterfaces]
        puts "  PPI interface: [join $t {; }][expr {$ni > 4 ? " +[expr {$ni - 4}]" : ""}]"
    }
    return $r
}

# ::ufv::tractability_async uid cb -> {symbol ensemblId groups {{modality labels}..} url} or {}
proc ::ufv::tractability_async {uid cb} {
    ::ufv::backend::async tractability [dict create uid [string toupper $uid]] $cb
}
# ufv_tractability ?uid?
proc ufv_tractability {{uid ""}} {
    set uid [::ufv::need $uid]
    set ot [::ufv::await [list ::ufv::tractability_async $uid]]
    if {![dict size $ot] || ![llength [dict get $ot groups]]} {
        ::ufv::status "$uid: no Open Targets tractability evidence."
        return $ot
    }
    puts "\[UFV\] [dict get $ot symbol] ([dict get $ot ensemblId]) tractability - [dict get $ot url]"
    foreach g [dict get $ot groups] { puts [format "  %-15s %s" [dict get $g modality] [join [dict get $g labels] {, }]] }
    return $ot
}

# ---------------------------------------------------------------------------------------------
# Ligand panel data (modal.js renderLigandPanel): copies, AlphaFill transplant meta, closest contact (clash),
# pocket evidence (mean pLDDT / B, sites, composition), chemistry + links + in-structure similarity
# ---------------------------------------------------------------------------------------------
# analysis.js ligandClash: closest heavy-atom distance ligand copy (HETATM) <-> protein (non-HETATM), current frame
proc ::ufv::feat::ligand_clash {mol resn chain resid} {
    set lsel "hetero and resname [list $resn] and not hydrogen and resid [::ufv::map::compress [list $resid]]"
    if {$chain ne ""} { append lsel " and chain [::ufv::map::chain_token $chain]" }
    set L [::ufv::with_sel s $mol $lsel { $s get {x y z} }]
    if {![llength $L]} { return "" }
    set P {}
    foreach r {6 15 60} {
        set P [::ufv::with_sel s $mol "not hetero and not hydrogen and within $r of ($lsel)" { $s get {x y z} }]
        if {[llength $P]} break
    }
    if {![llength $P]} { set P [::ufv::with_sel s $mol "not hetero and not hydrogen" { $s get {x y z} }] }
    if {![llength $P]} { return "" }
    set best 1e30
    foreach a $L { foreach b $P { set d [veclength [vecsub $a $b]]; if {$d < $best} { set best $d } } }
    return [dict create minDist $best band [expr {$best < 2.0 ? "severe" : ($best < 2.6 ? "suspicious" : "normal")}]]
}

# ::ufv::ligand_data_async mol resn chain resid cb -> dict
#   resn chain resid copies {{chain resid}..} near {pos..} dist clash {minDist band} aff {identity clash band
#   compound donorPdb rmsd} evidence {meanB band count sites composition compositionText} info {name formula weight
#   hbondDonors hbondAcceptors smiles inchikey drugbank cid links similar} motif 0|1 experimental 0|1
proc ::ufv::ligand_data_async {mol resn chain resid cb} {
    set copies {}; set codes {}
    foreach c [::ufv::ligand_copies $mol] {
        lassign $c rn ch ri
        if {$rn ni $codes && $rn ne "HOH"} { lappend codes $rn }
        if {$rn eq $resn} { lappend copies [list $ch $ri] }
    }
    set f [::ufv::st::get $mol focus ""]
    if {$f ne "" && [dict get $f kind] eq "ligand" && [dict get $f ligand] eq [list $resn $chain $resid] && [dict get $f dist] == $::ufv::nearby} {
        set near [dict get $f near]
    } else {
        set lig "resname [list $resn]"
        if {$chain ne ""} { append lig " and chain [::ufv::map::chain_token $chain]" }
        if {$resid ne ""} { append lig " and resid [::ufv::map::compress [list $resid]]" }
        set near {}
        ::ufv::with_sel p $mol "name CA and same residue as ($::ufv::poly_base[::ufv::feat::notx $mol] and within $::ufv::nearby of ($::ufv::lig_base and $lig))" {
            foreach {ch ri} [join [$p get {chain resid}]] {
                set u [::ufv::map::r2u $mol $ch $ri]
                if {$u ne "" && $u ni $near} { lappend near $u }
            }
        }
        set near [lsort -integer $near]
    }
    set out [dict create resn $resn chain $chain resid $resid copies $copies near $near dist $::ufv::nearby]
    dict set out clash [expr {$resid ne "" ? [::ufv::feat::ligand_clash $mol $resn $chain $resid] : ""}]
    set aff ""
    set meta [::ufv::st::get $mol aff ""]
    if {$chain ne "" && $meta ni {"" "-"} && [dict exists $meta $chain]} {
        set aff [dict get $meta $chain]
        if {[dict get $aff clash] ne ""} {
            set c [dict get $aff clash]
            dict set aff band [expr {$c <= 0.5 ? "low" : ($c <= 0.9 ? "moderate" : "high")}]
        }
    }
    dict set out aff $aff
    set st [::ufv::st::get $mol struct]
    set src [expr {[dict exists $st source] ? [dict get $st source] : [::ufv::st::get $mol source PDB]}]
    dict set out experimental [expr {$src in {"" PDB}}]
    set pdb [expr {[dict exists $st pdbId] ? [dict get $st pdbId] : ""}]
    dict set out motif [expr {[llength $near] >= 2 && ($pdb ni {"" -} || [::ufv::feat::af_based $mol])}]
    set steps {}
    if {[info exists ::ufv::feat::liginfo($resn)]} {
        dict set out info $::ufv::feat::liginfo($resn)
    } else {
        lappend steps info [list ::ufv::feat::liginfo_async $mol [linsert [lsearch -all -inline -not -exact $codes $resn] 0 $resn]]
    }
    if {[llength $near]} {
        set a [::ufv::feat::ctx $mol 1]
        dict set a positions [join $near +]
        lappend steps evidence [list ::ufv::backend::async pocket_evidence $a]
    }
    ::ufv::feat::all_of $steps [list ::ufv::feat::_ligand_data $out $resn $cb]
}
# Chemistry of ligand codes (first = the one shown; similarity among the structure's codes). Results are cached
# per code and announced with event `ligand_info {molid code info}` even when the panel request timed out.
proc ::ufv::feat::liginfo_async {mol codes cb} {
    ::ufv::backend::async ligand_info [dict create codes [join $codes " "]] [list ::ufv::feat::_liginfo_got $mol $codes $cb]
}
proc ::ufv::feat::_liginfo_got {mol codes cb ok res} {
    if {$ok} {
        dict for {code inf} $res { set ::ufv::feat::liginfo($code) $inf }
        if {[dict exists $res [lindex $codes 0]]} { ::ufv::emit ligand_info $mol [lindex $codes 0] [dict get $res [lindex $codes 0]] }
    }
    {*}$cb $ok $res
}
proc ::ufv::feat::_ligand_data {out resn cb res} {
    dict set out evidence ""
    if {![dict exists $out info]} { dict set out info "" }
    if {[dict exists $res evidence] && [lindex [dict get $res evidence] 0]} { dict set out evidence [lindex [dict get $res evidence] 1] }
    if {[dict exists $res info] && [lindex [dict get $res info] 0]} {
        set all [lindex [dict get $res info] 1]
        if {[dict exists $all $resn]} { dict set out info [dict get $all $resn] }
    }
    {*}$cb 1 $out
}

# ufv_ligand_info <resn> ?chain? ?resid? ?molid?  - the ligand panel as console text (returns the dict)
proc ufv_ligand_info {resn {chain ""} {resid ""} {molid ""}} {
    set mol [::ufv::resolve $molid]
    if {$resid eq ""} {
        foreach c [::ufv::ligand_copies $mol] { if {[lindex $c 0] eq $resn} { lassign $c - chain resid; break } }
    }
    set d [::ufv::await [list ::ufv::ligand_data_async $mol $resn $chain $resid]]
    set inf [dict get $d info]
    puts "\[UFV\] $resn $chain $resid ([llength [dict get $d copies]] cop[expr {[llength [dict get $d copies]] == 1 ? "y" : "ies"}])"
    if {$inf ne ""} {
        foreach {k lab} {name Name formula Formula weight "Mol. weight" hbondDonors "H-bond donors" hbondAcceptors "H-bond acceptors" inchikey InChIKey} {
            if {[dict exists $inf $k] && [dict get $inf $k] ne ""} { puts [format "  %-17s %s" $lab [dict get $inf $k]] }
        }
        foreach l [dict get $inf links] { puts [format "  %-17s %s" [lindex $l 0] [lindex $l 1]] }
        if {[llength [dict get $inf similar]]} {
            set t {}; foreach s [dict get $inf similar] { lappend t "[lindex $s 0] [format %.2f [lindex $s 1]]" }
            puts "  Similar in structure: [join $t {, }]"
        }
    }
    if {[dict get $d aff] ne ""} {
        set a [dict get $d aff]
        puts [format "  AlphaFill donor %s identity %s clash %s%s" [dict get $a donorPdb] [dict get $a identity] [dict get $a clash] \
            [expr {[dict exists $a band] ? " ([dict get $a band])" : ""}]]
    }
    if {[dict get $d clash] ne ""} { puts [format "  Closest contact   %.2f A (%s)" [dict get $d clash minDist] [dict get $d clash band]] }
    set ev [dict get $d evidence]
    if {$ev ne ""} {
        if {[dict get $ev meanB] ne ""} {
            puts [format "  Pocket %s       %.1f%s" [expr {[dict get $d experimental] ? "B-factor" : "pLDDT"}] [dict get $ev meanB] \
                [expr {[dict get $d experimental] ? "" : " ([dict get $ev band])"}]]
        }
        if {[llength [dict get $ev sites]]} {
            set t {}; foreach s [dict get $ev sites] { lappend t "[dict get $s desc] ([dict get $s pos])" }
            puts "  Functional sites  [join $t {; }]"
        }
        if {[dict get $ev compositionText] ne ""} { puts "  Residue character [dict get $ev compositionText]" }
    }
    puts "  Pocket (<= [dict get $d dist] A): [join [dict get $d near] {, }]"
    return $d
}

# ---------------------------------------------------------------------------------------------
# Exports: B-factor-encoded PDB (export.js rewritePdbBeta), selections (formatSelection)
# ---------------------------------------------------------------------------------------------
# UniProt positions of the residues drawn by the visible marker layers (all, or one: ptms | variants | sites |
# mutagenesis | binding | dom)
proc ::ufv::displayed_positions {mol {layer ""}} {
    set want [expr {$layer eq "" ? "" : [::ufv::layer_name $layer]}]
    set out [dict create]
    dict for {name groups} [::ufv::st::get $mol shown [dict create]] {
        if {$want ne "" && $name ne $want} continue
        if {![::ufv::layers::visible $mol $name]} continue
        dict for {hex pos} $groups { foreach p $pos { dict set out $p 1 } }
    }
    return [lsort -integer [dict keys $out]]
}

# ufv_save_bfactor ?filename? ?molid? ?mode?  - PDB with an analysis in the B-factor column: alphamissense,
# hotspots (1/0.67/0.33), hubs (1/0.5), burden (1), pocket (1 buried / 0.66 exposed, current q), plddt / bfactor
# (file values kept), selected (1 = residues of the visible marker layers). Unmapped atoms get 0; the molecule's
# own B-factors are restored afterwards.
proc ufv_save_bfactor {{filename ""} {molid ""} {mode selected}} {
    set mol [::ufv::resolve $molid]
    set mode [string tolower [string trim $mode]]
    if {$mode eq "am"} { set mode alphamissense }
    set a [::ufv::feat::ctx $mol [expr {$mode in {hotspots hubs contacthubs pocket pockets}}]]
    dict set a mode $mode
    if {$mode eq "selected"} { dict set a positions [join [::ufv::displayed_positions $mol] +] }
    set r [::ufv::backend::call bfactor_values $a]
    set keep [dict get $r keep]
    set vals [dict get $r values]
    if {$filename eq ""} {
        set nm [regsub -all {[^A-Za-z0-9_.-]} [molinfo $mol get name] _]
        set filename [file join [pwd] "${nm}_[dict get $r mode]_bfactor.pdb"]
    }
    set ours [::ufv::map::chains $mol]
    set all [atomselect $mol all]
    set saved [$all get beta]
    set cache [dict create]
    set out {}
    foreach k [$all get {chain resid insertion}] b $saved {
        if {![dict exists $cache $k]} {
            lassign $k ch ri ins
            set u [expr {$ch in $ours && [string trim $ins] eq "" ? [::ufv::map::r2u $mol $ch $ri] : ""}]
            dict set cache $k [expr {$u eq "" ? "0" : ($keep ? "keep" : ([dict exists $vals $u] ? [dict get $vals $u] : 0.0))}]
        }
        set v [dict get $cache $k]
        lappend out [expr {$v eq "keep" ? $b : $v}]
    }
    if {[catch {
        $all set beta $out
        $all writepdb $filename
    } err]} { $all set beta $saved; $all delete; error $err }
    $all set beta $saved
    $all delete
    ::ufv::status "wrote $filename (B-factor = [dict get $r mode])."
    return $filename
}

# ufv_select ?name? ?molid? ?layer?  - a VMD selection macro of the residues the marker layers show (all visible
# layers, or one); prints and returns the selection text.
proc ufv_select {{name ufv_sel} {molid ""} {layer ""}} {
    set mol [::ufv::resolve $molid]
    set pos [::ufv::displayed_positions $mol $layer]
    set s [::ufv::map::sel $mol $pos]
    if {$s eq ""} { ::ufv::status "nothing displayed to select."; return "" }
    catch { atomselect macro $name $s }
    set rs {}
    foreach ch [::ufv::map::chains $mol] { lappend rs {*}[::ufv::map::resids $mol $ch $pos] }
    set rs [lsort -integer -unique $rs]
    ::ufv::status "$name: [llength $rs] residue(s) - resid [::ufv::map::compress $rs]"
    return $s
}

# ---------------------------------------------------------------------------------------------
# Error guard: every async entry point and backend-reply handler calls its callback with `0 message` when it
# raises (a Tcl error must never leave a GUI request without an answer). The original proc is kept as <name>__impl.
# ---------------------------------------------------------------------------------------------
proc ::ufv::feat::guard_proc {name cbidx} {
    if {![llength [info procs $name]] || [llength [info procs ${name}__impl]]} return
    rename $name ${name}__impl
    proc $name args [format {
        if {[catch { uplevel 1 [list %s__impl {*}$args] } res opts]} {
            set cb [lindex $args %s]
            if {$cb ne "" && [llength $args] > 0} {
                ::ufv::status "%s: $res"
                after 0 [list {*}$cb 0 $res]
                return
            }
            return -options $opts $res
        }
        return $res
    } $name $cbidx [namespace tail $name]]
}
foreach {p i} {
    ::ufv::fetch_async end  ::ufv::load_af_async end  ::ufv::load_structure_async end  ::ufv::map_sifts_async end
    ::ufv::map_custom_async end  ::ufv::marker_async 4  ::ufv::analysis_async end  ::ufv::tiers_async end
    ::ufv::tier_async end  ::ufv::ligands_async end  ::ufv::report_async 2  ::ufv::csv_async end
    ::ufv::residue_async end  ::ufv::pockets_all_async end  ::ufv::motif_async end  ::ufv::partner_layers_async end
    ::ufv::partners_async end  ::ufv::pdbekb_async end  ::ufv::ligand_data_async end  ::ufv::diseases_async end
    ::ufv::binding_async end  ::ufv::tractability_async end
    ::ufv::_fetched end-2  ::ufv::_af_downloaded end-2  ::ufv::_af_loaded end-2  ::ufv::_prepared end-2
    ::ufv::_prepared2 3  ::ufv::_struct_mapped end-2  ::ufv::_sifts_mapped end-2  ::ufv::_custom_mapped end-2
    ::ufv::_marker_got end-2  ::ufv::_tier_got end-2  ::ufv::_tiers_join end-2  ::ufv::_aff_got end-2
    ::ufv::_ligands_build 1  ::ufv::_report_done end-2  ::ufv::_csv_done end-2  ::ufv::_residue_done end-2
    ::ufv::feat::_partners_got end-2  ::ufv::feat::_partner_got end-2  ::ufv::feat::_pa_done end-2
    ::ufv::feat::_kb_got end-2  ::ufv::feat::_ligand_data 2  ::ufv::feat::_liginfo_got end-2
} { ::ufv::feat::guard_proc $p $i }
