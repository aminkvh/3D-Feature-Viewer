# Installs the VMD plugin for the current user (Windows). Run from the unzipped folder:
#   powershell -ExecutionPolicy Bypass -File install-vmd.ps1
# Copies the plugin to %USERPROFILE%\ufv-vmd, adds it to %USERPROFILE%\vmd.rc (backed up first) and creates the
# Start-menu shortcut "VMD (3D Feature Viewer)". The stock VMD shortcut ignores a user vmd.rc, so use the new one.
$ErrorActionPreference = "Stop"
$src = $PSScriptRoot
$dst = Join-Path $env:USERPROFILE "ufv-vmd"
$rc  = Join-Path $env:USERPROFILE "vmd.rc"

if (-not (Test-Path (Join-Path $src "vmd\ufv_vmd.tcl"))) { throw "Run this script from the unzipped vmd-plugin folder." }
if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
New-Item -ItemType Directory -Force (Join-Path $dst "vmd") | Out-Null
Copy-Item (Join-Path $src "ufv_core.py") $dst
Copy-Item (Join-Path $src "vmd\*.tcl") (Join-Path $dst "vmd")

$block = @'
# >>> 3D Feature Viewer for UniProt >>>
# A user vmd.rc replaces VMD's own: keep its defaults (lights, menus).
if {[info exists env(VMDDIR)] && [file exists [file join $env(VMDDIR) vmd.rc]]} { source [file join $env(VMDDIR) vmd.rc] }
set ufv_dir [file join $env(USERPROFILE) ufv-vmd]
if {![info exists env(UFV_NO_AUTOLOAD)] && [file exists $ufv_dir/vmd/ufv_vmd.tcl] &&[catch { source $ufv_dir/vmd/ufv_vmd.tcl } e]} { puts "UFV: $e" }
catch { vmd_install_extension ufv ufv_gui "Analysis/3D Feature Viewer for UniProt" }
# <<< 3D Feature Viewer for UniProt <<<
'@
$old = ""
if (Test-Path $rc) {
    Copy-Item $rc "$rc.bak" -Force
    $old = [regex]::Replace((Get-Content $rc -Raw), "(?s)# >>> 3D Feature Viewer for UniProt >>>.*?# <<< 3D Feature Viewer for UniProt <<<\r?\n?", "").TrimEnd()
    if ($old) { $old += "`r`n" }
}
Set-Content -Path $rc -Value ($old + $block) -Encoding ascii -NoNewline

$vmd = @("$env:ProgramFiles\VMD\vmd.exe", "${env:ProgramFiles(x86)}\VMD\vmd.exe",
         "${env:ProgramFiles(x86)}\University of Illinois\VMD\vmd.exe") | Where-Object { Test-Path $_ } | Select-Object -First 1
if ($vmd) {
    $lnk = Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs\VMD (3D Feature Viewer).lnk"
    $l = (New-Object -ComObject WScript.Shell).CreateShortcut($lnk)
    $l.TargetPath = $vmd; $l.Arguments = "-startup `"$rc`""; $l.WorkingDirectory = $env:USERPROFILE
    $l.IconLocation = "$vmd,0"; $l.Save()
    Write-Host "Done. Start VMD from the Start-menu shortcut 'VMD (3D Feature Viewer)'."
} else {
    Write-Host "Done. vmd.exe not found in Program Files: start VMD with  vmd.exe -startup `"$rc`""
}
Write-Host "Then open Extensions > Analysis > 3D Feature Viewer for UniProt."
