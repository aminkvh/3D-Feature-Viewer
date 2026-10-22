# =============================================================================================
# ufv_backend.tcl - Python backend channel (ufv_core.py --serve) for the VMD plugin.
# OWNER: foundation. Other modules call ::ufv::backend::async / call / await; do not touch the channel.
#
#   ::ufv::backend::async op ?argdict? ?callback?   -> request id; later: {*}$callback ok result
#   ::ufv::backend::call  op ?argdict? ?timeout_ms? -> result (vwait, Tk/VMD events keep running) or error
#   ::ufv::await cmdprefix ?timeout_ms?             -> runs `{*}$cmdprefix $cb` and waits for $cb ok result
#
# The backend is ONE persistent process (python -u ufv_core.py --serve), opened with `open |... r+`, fileevent
# driven. Requests are JSON lines; replies are Tcl dict lines (id <n> ok 0|1 result|error <value>) read with
# `dict get` - backend output is data and is never eval'ed. Coordinates go through files (JSON / PDB / DCD).
# A crashed backend is restarted on the next request (pending requests fail with an error).
# =============================================================================================

namespace eval ::ufv {
    if {![info exists dir]} { variable dir [file dirname [file normalize [info script]]] }
}

namespace eval ::ufv::backend {
    variable python_user ""   ;# set by ufv_python (path, or a command list such as {py -3})
    variable python ""        ;# resolved interpreter: absolute path of a python >= 3.8 executable
    variable python_cached 0  ;# 1 = $python came from the prefs file (not probed in this session)
    variable core_user ""     ;# set by ufv_core
    variable chan ""
    variable nextid 0
    variable tick 0           ;# bumped on every reply / failure (sync waits vwait on it)
    variable pending          ;# id -> {callback t0_ms op}
    variable done             ;# id -> {ok payload}  (sync waits)
    variable waiting          ;# id -> 1             (sync waits)
    variable stats            ;# op -> last round-trip ms
    variable crashes {}       ;# clock seconds of recent unexpected exits
    variable timeout 300000   ;# default sync timeout (ms)
    variable logfile ""
    variable started 0        ;# processes started (tests)
    array set pending {}
    array set done {}
    array set waiting {}
    array set stats {}
}

# ---------------------------------------------------------------------------------------------
# JSON request encoding (requests only; replies are Tcl dicts)
# ---------------------------------------------------------------------------------------------
namespace eval ::ufv::json {}

proc ::ufv::json::str {s} {
    set m [list \\ \\\\ \" \\\" \n \\n \r \\r \t \\t \b \\b \f \\f]
    set s [string map $m $s]
    # remaining control characters
    if {[regexp {[\x00-\x1f]} $s]} {
        set out ""
        foreach ch [split $s ""] {
            scan $ch %c c
            if {$c < 32} { append out [format "\\u%04x" $c] } else { append out $ch }
        }
        set s $out
    }
    return "\"$s\""
}

# Mark a value as pre-encoded JSON (used inside obj / arr).
proc ::ufv::json::raw {j} { return "\u0001$j" }

# Auto-typed scalar: strict JSON numbers stay numbers, everything else is a string; raw passes through.
proc ::ufv::json::val {v} {
    if {[string index $v 0] eq "\u0001"} { return [string range $v 1 end] }
    if {[regexp {^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?$} $v]} { return $v }
    return [str $v]
}

proc ::ufv::json::obj {d} {
    set parts {}
    dict for {k v} $d { lappend parts "[str $k]:[val $v]" }
    return "{[join $parts ,]}"
}

proc ::ufv::json::arr {l} {
    set parts {}
    foreach v $l { lappend parts [val $v] }
    return "\[[join $parts ,]\]"
}

# {chain {{resid aa} ...} ...} -> {"A":[[1,"M"],...],...}  (residue rows for map_sifts / map_custom)
proc ::ufv::json::residues {d} {
    set parts {}
    dict for {ch rows} $d {
        set r {}
        foreach row $rows { lappend r "\[[lindex $row 0],[str [lindex $row 1]]\]" }
        lappend parts "[str $ch]:\[[join $r ,]\]"
    }
    return "{[join $parts ,]}"
}

# {chain {{us ue rs re} ...}} -> {"A":[[us,ue,rs,re],...]}
proc ::ufv::json::segments {d} {
    set parts {}
    dict for {ch segs} $d {
        set r {}
        foreach s $segs { lappend r "\[[join $s ,]\]" }
        lappend parts "[str $ch]:\[[join $r ,]\]"
    }
    return "{[join $parts ,]}"
}

# Custom-map chain configs {chain {mode m start n accession a segments {{us ue rs re}..} src s}} -> JSON
proc ::ufv::json::configs {d} {
    set parts {}
    dict for {ch c} $d {
        set kv {}
        dict for {k v} $c {
            if {$k eq "segments"} {
                set r {}
                foreach s $v { lappend r "\[[join $s ,]\]" }
                lappend kv "[str $k]:\[[join $r ,]\]"
            } else {
                lappend kv "[str $k]:[val $v]"
            }
        }
        lappend parts "[str $ch]:{[join $kv ,]}"
    }
    return "{[join $parts ,]}"
}

# ---------------------------------------------------------------------------------------------
# Interpreter / core discovery
# ---------------------------------------------------------------------------------------------
# Probe one candidate (command list). Returns the absolute executable of a python >= 3.8, or "".
proc ::ufv::backend::probe {cmdlist} {
    set exe [auto_execok [lindex $cmdlist 0]]
    if {$exe eq "" && [file executable [lindex $cmdlist 0]]} { set exe [list [lindex $cmdlist 0]] }
    if {$exe eq ""} { return "" }
    set code {import sys; print(sys.executable); print(int(sys.version_info >= (3, 8))); print(sys.version.split()[0])}
    if {[catch {exec {*}$exe {*}[lrange $cmdlist 1 end] -c $code 2>@1} out]} { return "" }
    set lines [split [string map [list \r ""] [string trim $out]] \n]
    if {[llength $lines] < 3} { return "" }
    lassign [lrange $lines end-2 end] path ok ver
    if {$ok ne "1" || ![file exists $path]} { return "" }
    return [file normalize $path]
}

proc ::ufv::backend::candidates {} {
    variable python_user
    set c {}
    if {$python_user ne ""} { lappend c $python_user }
    if {[info exists ::env(UFV_PYTHON)] && $::env(UFV_PYTHON) ne ""} { lappend c $::env(UFV_PYTHON) }
    if {$::tcl_platform(platform) eq "windows"} { lappend c {py -3} }
    lappend c python3 python
    return $c
}

# Resolved interpreter (auto-detected once; ufv_python overrides). Errors with a clear message if none.
# The auto-detected path is remembered in a small prefs file: probing runs `py -3 -c ...` synchronously
# (~150 ms) and would otherwise be paid by the first request of every session. A cached path that no longer
# exists is ignored; a backend that dies right away forgets it (died) and the next start probes again.
proc ::ufv::backend::pyfile {} {
    set d [expr {[info exists ::env(UFV_CACHE_DIR)] && $::env(UFV_CACHE_DIR) ne "" ? $::env(UFV_CACHE_DIR) : [file normalize ~]}]
    return [file join $d ufv_vmd_python.txt]
}
proc ::ufv::backend::python {} {
    variable python
    variable python_user
    variable python_cached
    if {$python ne ""} { return $python }
    set auto [expr {$python_user eq "" && !([info exists ::env(UFV_PYTHON)] && $::env(UFV_PYTHON) ne "")}]
    if {$auto && ![catch { set fh [open [pyfile]]; set p [string trim [read $fh]]; close $fh }] &&
        $p ne "" && [file isfile $p] && [file executable $p]} {
        set python_cached 1
        set python $p
        return $p
    }
    set tried {}
    foreach cand [candidates] {
        lappend tried [join $cand " "]
        set p [probe $cand]
        if {$p ne ""} {
            set python $p
            set python_cached 0
            if {$auto} { catch { set fh [open [pyfile] w]; puts $fh $p; close $fh } }
            return $p
        }
        if {$python_user ne "" && $cand eq $python_user} {
            puts "\[UFV\] WARNING: ufv_python '[join $python_user]' is not a Python >= 3.8 - trying the defaults."
        }
    }
    error "no Python >= 3.8 found (tried: [join $tried {, }]). Install Python 3 or set it with: ufv_python <path to python.exe>"
}

proc ::ufv::backend::core {} {
    variable core_user
    set c {}
    if {$core_user ne ""} { lappend c $core_user }
    if {[info exists ::env(UFV_CORE)] && $::env(UFV_CORE) ne ""} { lappend c $::env(UFV_CORE) }
    lappend c [file join $::ufv::dir ufv_core.py] [file join [file dirname $::ufv::dir] ufv_core.py]
    foreach p $c { if {[file isfile $p]} { return [file normalize $p] } }
    error "ufv_core.py not found (looked in: [join $c {, }]). Set it with: ufv_core <path to ufv_core.py>"
}

proc ::ufv::backend::tmpdir {} {
    foreach v {TEMP TMP TMPDIR} {
        if {[info exists ::env($v)] && [string length $::env($v)] > 3 && [file isdirectory $::env($v)]} {
            return [file normalize $::env($v)]
        }
    }
    return [expr {$::tcl_platform(platform) eq "windows" ? "C:/Windows/Temp" : "/tmp"}]
}

# ---------------------------------------------------------------------------------------------
# Process lifecycle
# ---------------------------------------------------------------------------------------------
proc ::ufv::backend::running {} { variable chan; return [expr {$chan ne ""}] }

proc ::ufv::backend::start {} {
    variable chan
    variable logfile
    variable started
    variable crashes
    if {$chan ne ""} { return $chan }
    set now [clock seconds]
    set recent {}
    foreach t $crashes { if {$now - $t < 60} { lappend recent $t } }
    set crashes $recent
    if {[llength $crashes] >= 5} {
        error "the Python backend crashed [llength $crashes] times in the last minute - see $logfile (ufv_backend restart to retry)"
    }
    set py [python]
    set core [core]
    set logfile [file join [tmpdir] "ufv_backend_[pid].log"]
    set chan [open |[list $py -u $core --serve 2> $logfile] r+]
    fconfigure $chan -blocking 0 -buffering line -translation {auto lf} -encoding utf-8
    fileevent $chan readable [list ::ufv::backend::readable $chan]
    incr started
    if {$started == 1} { after idle ::ufv::backend::sweep }
    return $chan
}

# Temp files of earlier sessions (VMD killed / closed without cleanup): the plugin's own names only, older than
# an hour, never this process's. Files still in use elsewhere are locked (backend log) or regenerated on
# demand (Calpha JSON: ::ufv::feat::coords checks the file exists), so a stale delete is harmless.
# At most once an hour across sessions (stamp file): stat-ing a few hundred files costs 20-40 ms.
proc ::ufv::backend::sweep {{age 3600} {force 0}} {
    set now [clock seconds]
    set stamp [file join [tmpdir] ufv_vmd_sweep.stamp]
    if {!$force && ![catch { file mtime $stamp } t] && $now - $t < $age} { return 0 }
    catch { close [open $stamp w] }
    set me [pid]
    set n 0
    foreach f [glob -nocomplain -directory [tmpdir] -types f ufv_backend_*.log ufv_ca_*.json ufv_traj_*.* ufv_ss_*.pdb] {
        if {![regexp {^ufv_(?:backend|ca|traj|ss)_(\d+)[._]} [file tail $f] -> p] || $p == $me} continue
        if {[catch { file mtime $f } t] || $now - $t < $age} continue
        if {![catch { file delete $f }]} { incr n }
    }
    return $n
}

# Unexpected exit / write failure: fail every pending request, forget the channel (restart on next request).
proc ::ufv::backend::died {reason {crash 1}} {
    variable chan
    variable pending
    variable crashes
    if {$chan ne ""} {
        catch { fileevent $chan readable {} }
        catch { fconfigure $chan -blocking 1 }
        catch { close $chan }
        set chan ""
    }
    if {$crash} {
        lappend crashes [clock seconds]
        forget_python
    }
    foreach id [array names pending] { finish $id 0 "backend exited: $reason" }
}
# A remembered interpreter is dropped after a crash (it may be gone or replaced): the next start probes again.
proc ::ufv::backend::forget_python {} {
    variable python
    variable python_cached
    if {[info exists python_cached] && $python_cached} {
        set python ""
        set python_cached 0
        catch { file delete [pyfile] }
    }
}

proc ::ufv::backend::stop {} {
    variable chan
    if {$chan eq ""} return
    catch { puts $chan {{"id":0,"op":"shutdown"}}; flush $chan }
    died "stopped" 0
}

proc ::ufv::backend::restart {} {
    variable crashes
    variable python
    set crashes {}
    stop
    forget_python
    set python ""
    start
}

# ---------------------------------------------------------------------------------------------
# Requests
# ---------------------------------------------------------------------------------------------
proc ::ufv::backend::async {op {argdict {}} {callback ""}} {
    set id [incr ::ufv::backend::nextid]
    send $id $op $argdict $callback
    return $id
}

proc ::ufv::backend::send {id op argdict callback} {
    variable chan
    variable pending
    set req "{\"id\":$id,\"op\":[::ufv::json::str $op],\"fmt\":\"tcl\",\"args\":[::ufv::json::obj $argdict]}"
    set pending($id) [list $callback [clock milliseconds] $op]
    if {[catch { start } err]} {
        finish $id 0 $err
        return
    }
    if {[catch { puts $chan $req; flush $chan }]} {
        # the process died since the last request: restart once and resend
        set keep $pending($id); unset pending($id)
        died "write failed"
        set pending($id) $keep
        if {[catch { start; puts $chan $req; flush $chan } err]} { finish $id 0 "backend unavailable: $err" }
    }
}

proc ::ufv::backend::readable {ch} {
    variable chan
    if {$ch ne $chan} { catch { close $ch }; return }
    while {[gets $ch line] >= 0} {
        if {$line eq ""} continue
        if {[catch { handle $line } err]} {
            puts "\[UFV\] backend reply ignored: $err"
        }
    }
    if {[eof $ch]} { died "process exited (see $::ufv::backend::logfile)" }
}

proc ::ufv::backend::handle {line} {
    # a reply is a Tcl dict: id <n> ok 0|1 result|error <value>   (data only - never eval'ed)
    if {[catch { dict size $line }] || ![dict exists $line id] || ![dict exists $line ok]} {
        error "malformed reply: [string range $line 0 120]"
    }
    set ok [dict get $line ok]
    if {$ok} { set payload [dict get $line result] } else { set payload [dict get $line error] }
    finish [dict get $line id] $ok $payload
}

proc ::ufv::backend::finish {id ok payload} {
    variable pending
    variable waiting
    variable done
    variable tick
    variable stats
    if {![info exists pending($id)]} return       ;# cancelled or timed out
    lassign $pending($id) cb t0 op
    unset pending($id)
    set stats($op) [expr {[clock milliseconds] - $t0}]
    if {[info exists waiting($id)]} { set done($id) [list $ok $payload] }
    if {$cb ne ""} { after 0 [list ::ufv::backend::dispatch $cb $ok $payload] }
    incr tick
}

proc ::ufv::backend::dispatch {cb ok payload} {
    if {[catch { uplevel #0 [linsert $cb end $ok $payload] } err]} {
        puts "\[UFV\] callback error: $err"
        if {[info exists ::errorInfo]} { puts $::errorInfo }
    }
}

proc ::ufv::backend::cancel {id} {
    variable pending
    unset -nocomplain pending($id)
}

proc ::ufv::backend::busy {} { variable pending; return [array size pending] }

# Synchronous call (scripting / console). The Tcl event loop keeps running while waiting.
proc ::ufv::backend::call {op {argdict {}} {timeout ""}} {
    variable waiting
    variable done
    if {$timeout eq ""} { set timeout $::ufv::backend::timeout }
    set id [incr ::ufv::backend::nextid]
    set waiting($id) 1
    send $id $op $argdict ""
    set aid [after $timeout [list ::ufv::backend::expire $id]]
    while {![info exists done($id)]} { vwait ::ufv::backend::tick }
    after cancel $aid
    lassign $done($id) ok payload
    unset done($id) waiting($id)
    if {$ok eq "timeout"} { error "backend $op timed out after $timeout ms" }
    if {!$ok} { error "backend $op: $payload" }
    return $payload
}

proc ::ufv::backend::expire {id} {
    variable pending
    variable done
    variable tick
    unset -nocomplain pending($id)
    set done($id) [list timeout ""]
    incr tick
}

# Generic "wait for an async operation": cmdprefix is called with one extra argument, a callback that
# the operation invokes as {*}$cb ok result. Returns result or raises the error.
proc ::ufv::await {cmdprefix {timeout ""}} {
    variable _aw
    variable _awn
    if {![info exists _awn]} { set _awn 0 }
    if {$timeout eq ""} { set timeout $::ufv::backend::timeout }
    set k [incr _awn]
    set aid [after $timeout [list ::ufv::_awaited $k timeout "timed out after $timeout ms"]]
    if {[catch { uplevel #0 [linsert $cmdprefix end [list ::ufv::_awaited $k]] } err]} {
        after cancel $aid
        unset -nocomplain _aw($k)
        error $err
    }
    while {![info exists _aw($k)]} { vwait ::ufv::backend::tick }
    after cancel $aid
    lassign $_aw($k) ok res
    unset _aw($k)
    if {$ok eq "timeout"} { error $res }
    if {!$ok} { error $res }
    return $res
}

proc ::ufv::_awaited {k ok res} {
    variable _aw
    if {[info exists _aw($k)]} return
    set _aw($k) [list $ok $res]
    incr ::ufv::backend::tick
}

# One-shot CLI fallback (no server): exec python ufv_core.py <args>. Only for diagnostics.
proc ::ufv::backend::oneshot {args} {
    return [exec [python] [core] {*}$args]
}

# ---------------------------------------------------------------------------------------------
# User commands
# ---------------------------------------------------------------------------------------------
proc ufv_python {{path ""}} {
    if {$path eq ""} { return [::ufv::backend::python] }
    set ::ufv::backend::python_user $path
    set ::ufv::backend::python ""
    set p [::ufv::backend::probe $path]
    if {$p eq ""} {
        puts "\[UFV\] WARNING: '$path' is not a Python >= 3.8; auto-detection will be used."
    } else {
        set ::ufv::backend::python $p
        puts "\[UFV\] python = $p"
    }
    ::ufv::backend::stop
    return $p
}

proc ufv_core {{path ""}} {
    if {$path eq ""} { return [::ufv::backend::core] }
    set ::ufv::backend::core_user $path
    ::ufv::backend::stop
    puts "\[UFV\] ufv_core = $path"
    return $path
}

# ufv_backend ?status|restart|stop|ping?
proc ufv_backend {{what status}} {
    switch -- $what {
        restart { ::ufv::backend::restart; return [ufv_backend ping] }
        stop    { ::ufv::backend::stop; return stopped }
        ping    {
            set t0 [clock milliseconds]
            set r [::ufv::backend::call ping {} 20000]
            puts "\[UFV\] backend: python [dict get $r python] ([dict get $r executable]), pid [dict get $r pid],\
                  [expr {[clock milliseconds] - $t0}] ms"
            return $r
        }
        default {
            set py [expr {$::ufv::backend::python ne "" ? $::ufv::backend::python : "(not resolved yet)"}]
            set s [list running [::ufv::backend::running] python $py pending [::ufv::backend::busy] \
                       started $::ufv::backend::started log $::ufv::backend::logfile]
            puts "\[UFV\] backend: $s"
            return $s
        }
    }
}
