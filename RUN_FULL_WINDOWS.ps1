param(
    [ValidatePattern('^(?i:auto|[1-9][0-9]*)$')][string]$Threads = 'auto',
    [switch]$NoBuild,
    [string]$MSBuildPath = ''
)
$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
if (-not $NoBuild) {
    & (Join-Path $root 'BUILD_RELEASE.ps1') -MSBuildPath $MSBuildPath
}
$exe = Join-Path $root 'x64\Release\N5Exact.exe'
if (-not (Test-Path -LiteralPath $exe)) {
    throw "Release executable not found: $exe"
}
if ($Threads -match '^(?i:auto)$') {
    Write-Host 'Running N5 exact production with automatic participating-thread selection.'
} else {
    Write-Host "Running N5 exact production with $Threads requested participating threads."
}
$threadArgument = if ($Threads -match '^(?i:auto)$') { 'auto' } else { $Threads }
& $exe --threads $threadArgument
$exitCode = $LASTEXITCODE
if ($exitCode -ne 0) {
    throw "Production process exited with code $exitCode."
}
