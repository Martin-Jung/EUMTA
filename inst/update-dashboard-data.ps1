param(
  [ValidatePattern('^\d{2}$')]
  [string]$Version = "02"
)

$ErrorActionPreference = "Stop"
$repositoryRoot = Split-Path -Parent $PSScriptRoot
$builder = Join-Path $repositoryRoot "build-dashboard-data.ps1"

& $builder -Version $Version
