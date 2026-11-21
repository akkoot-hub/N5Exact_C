# Run inside Windows PowerShell 5.1/7 on a machine with VS 2026 C++ tools.
param([string]$MSBuildPath = '')
$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$project = Join-Path $root 'N5Exact.vcxproj'
if (-not (Test-Path -LiteralPath $project)) { throw "Missing project: $project" }

if ([string]::IsNullOrWhiteSpace($MSBuildPath)) {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (Test-Path -LiteralPath $vswhere) {
        $installation = & $vswhere -latest -products '*' -version '[18.0,19.0)' -requires 'Microsoft.Component.MSBuild' 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64' -property installationPath
        if ($installation) {
            $candidate = Join-Path ($installation | Select-Object -First 1) 'MSBuild\Current\Bin\MSBuild.exe'
            if (Test-Path -LiteralPath $candidate) { $MSBuildPath = $candidate }
        }
    }
}
if ([string]::IsNullOrWhiteSpace($MSBuildPath)) {
    $cmd = Get-Command 'MSBuild.exe' -ErrorAction SilentlyContinue
    if ($cmd) { $MSBuildPath = $cmd.Source }
}
if ([string]::IsNullOrWhiteSpace($MSBuildPath) -or -not (Test-Path -LiteralPath $MSBuildPath)) {
    throw 'MSBuild.exe not found. Install Visual Studio 2026 (18.x) with Desktop development with C++ / MSVC v145, or pass -MSBuildPath <full path>.'
}

Write-Host "MSBUILD = $MSBuildPath"
Write-Host 'BUILD = Release | x64 | Visual Studio 2026 v145 | ISO C17'
& $MSBuildPath $project '/m' '/p:Configuration=Release' '/p:Platform=x64' '/p:PreferredToolArchitecture=x64' '/verbosity:minimal'
$code = $LASTEXITCODE
if ($code -ne 0) { throw "MSBuild failed: exit $code" }
$exe = Join-Path $root 'x64\Release\N5Exact.exe'
if (-not (Test-Path -LiteralPath $exe)) { throw "Build exit 0 but executable missing: $exe" }
Write-Host "RELEASE_BUILD_PASS: $exe"
