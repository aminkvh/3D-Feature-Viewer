# Installs the 3D Feature Viewer plugin into every VMD found on this PC (VMD 1.9.x and VMD 2).
#   powershell -ExecutionPolicy Bypass -File install-vmd.ps1              install / update
#   powershell -ExecutionPolicy Bypass -File install-vmd.ps1 -Uninstall    remove
#   powershell -ExecutionPolicy Bypass -File install-vmd.ps1 -VmdDir "D:\Apps\VMD"   also use this VMD folder
#   (add -Only to use just the folder(s) given with -VmdDir)
# The plugin is copied to %USERPROFILE%\3d-feature-viewer and registered in each VMD's own vmd.rc (backed up
# first), so your normal VMD shortcut shows Extensions > Analysis > 3D Feature Viewer. No extra shortcut is made.
# Program Files is protected, so Windows will ask once for administrator permission.
param(
    [switch]$Uninstall,
    [string[]]$VmdDir = @(),
    [switch]$Only,
    [string]$UserHome = $env:USERPROFILE,
    [switch]$Elevated
)
$ErrorActionPreference = "Stop"
$src   = $PSScriptRoot
$dst   = Join-Path $UserHome "3d-feature-viewer"
$begin = "# >>> 3D Feature Viewer >>>"
$end   = "# <<< 3D Feature Viewer <<<"

function Find-Vmd {
    if ($Only) { $roots = @() } else { $roots = @($env:ProgramFiles, ${env:ProgramFiles(x86)}, (Join-Path $env:LOCALAPPDATA "Programs")) | Where-Object { $_ -and (Test-Path $_) } }
    $dirs = @($VmdDir)
    foreach ($r in $roots) {
        $dirs += Get-ChildItem $r -Directory -Filter "VMD*" -ErrorAction SilentlyContinue | ForEach-Object FullName
        $dirs += Get-ChildItem (Join-Path $r "University of Illinois") -Directory -ErrorAction SilentlyContinue | ForEach-Object FullName
    }
    if (-not $Only) { $dirs += Get-ItemProperty "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\*", "HKLM:\SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\*" -ErrorAction SilentlyContinue |
        Where-Object { $_.DisplayName -like "VMD*" -and $_.InstallLocation } | ForEach-Object { $_.InstallLocation.TrimEnd('\') } }
    $dirs | Where-Object { $_ -and (Test-Path (Join-Path $_ "vmd.exe")) } | ForEach-Object { (Resolve-Path $_).Path } | Select-Object -Unique
}

function Strip-Block($text) {
    # drops this plugin's block, and the block written by 2.1.0 under its old name
    foreach ($n in "3D Feature Viewer", "3D Feature Viewer for UniProt") {
        $text = [regex]::Replace($text, "(?s)# >>> $n >>>.*?# <<< $n <<<\r?\n?", "")
    }
    $text.TrimEnd()
}

$vmds = @(Find-Vmd)
if (-not $Uninstall -and -not (Test-Path (Join-Path $src "vmd\ufv_vmd.tcl"))) { throw "Run this script from the unzipped vmd-plugin folder." }
if (-not $vmds) {
    Write-Host "No VMD found in Program Files. Run again with  -VmdDir `"<folder that contains vmd.exe>`"  or follow INSTALL.txt." -ForegroundColor Yellow
    exit 1
}

# Program Files is not writable for a normal user: relaunch once with administrator rights.
$isAdmin = ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
if (-not $isAdmin -and -not $Elevated) {
    Write-Host "Windows will ask for administrator permission to edit the VMD startup file..."
    $a = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", "`"$PSCommandPath`"", "-Elevated", "-UserHome", "`"$UserHome`"")
    if ($Uninstall) { $a += "-Uninstall" }
    foreach ($d in $VmdDir) { $a += @("-VmdDir", "`"$d`"") }
    if ($Only) { $a += "-Only" }
    Start-Process powershell -Verb RunAs -Wait -ArgumentList $a
    exit
}

if ($Uninstall) {
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
} else {
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
    New-Item -ItemType Directory -Force (Join-Path $dst "vmd") | Out-Null
    Copy-Item (Join-Path $src "ufv_core.py") $dst
    Copy-Item (Join-Path $src "vmd\*.tcl") (Join-Path $dst "vmd")
}

$tclDir = $dst.Replace('\', '/')
$block = @"
$begin
set ::ufv_dir {$tclDir}
if {![info exists env(UFV_NO_AUTOLOAD)] && [file exists `$::ufv_dir/vmd/ufv_vmd.tcl] && [catch { source `$::ufv_dir/vmd/ufv_vmd.tcl } e]} { puts "3D Feature Viewer: `$e" }
catch { vmd_install_extension ufv ufv_gui "Analysis/3D Feature Viewer" }
$end
"@

foreach ($v in $vmds) {
    $rc = Join-Path $v "vmd.rc"
    $old = ""
    if (Test-Path $rc) {
        if (-not (Test-Path "$rc.3dfv-backup")) { Copy-Item $rc "$rc.3dfv-backup" }
        $old = Strip-Block (Get-Content $rc -Raw)
        if ($old) { $old += "`r`n" }
    }
    $new = if ($Uninstall) { $old } else { $old + $block + "`r`n" }
    Set-Content -Path $rc -Value $new -Encoding ascii -NoNewline
    Write-Host ("{0}: {1}" -f $(if ($Uninstall) { "removed from" } else { "installed in" }), $v) -ForegroundColor Green
}

# Clean up what 2.1.0 left behind: the extra "VMD (3D Feature Viewer)" shortcut, the old folder, the old ~\vmd.rc entry.
$lnk = Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs\VMD (3D Feature Viewer).lnk"
$lnk2 = Join-Path (Join-Path $UserHome "AppData\Roaming") "Microsoft\Windows\Start Menu\Programs\VMD (3D Feature Viewer).lnk"
foreach ($l in $lnk, $lnk2) { if (Test-Path $l) { Remove-Item $l -Force } }
$oldDir = Join-Path $UserHome "ufv-vmd"
if (Test-Path $oldDir) { Remove-Item $oldDir -Recurse -Force }
$userRc = Join-Path $UserHome "vmd.rc"
if (Test-Path $userRc) {
    $t = Get-Content $userRc -Raw
    if ($t -match "3D Feature Viewer") {
        $rest = Strip-Block $t
        if ($rest) { Set-Content $userRc ($rest + "`r`n") -Encoding ascii -NoNewline } else { Remove-Item $userRc -Force }
    }
}

if ($Uninstall) { Write-Host "Done. Restart VMD." } else {
    Write-Host "Done. Restart VMD, then open Extensions > Analysis > 3D Feature Viewer."
}
if ($Elevated) { Start-Sleep -Seconds 4 }
