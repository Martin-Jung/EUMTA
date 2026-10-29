param(
  [ValidatePattern('^\d{2}$')]
  [string]$Version = "02"
)

# Run after updating the export CSV files to refresh the file://-compatible dashboard bundle.
$ErrorActionPreference = "Stop"
$exportPath = Join-Path $PSScriptRoot "export"
function Import-TargetCsvs {
  param([string]$Pattern)
  $files = @(Get-ChildItem -Path (Join-Path $exportPath $Pattern) -ErrorAction SilentlyContinue | Sort-Object Name)
  if ($files.Count -eq 0) { throw "No export files found for version $Version using pattern $Pattern" }
  $files | ForEach-Object { Import-Csv $_.FullName }
}

$datasets = [ordered]@{
  overall = @(Import-TargetCsvs "MTA_*_MS_All_All_$Version.csv")
  art12 = @(Import-TargetCsvs "MTA_*_EU_Art12_species_$Version.csv")
  art17Species = @(Import-TargetCsvs "MTA_*_MS_Art17_species_$Version.csv")
  art17Habitats = @(Import-TargetCsvs "MTA_*_MS_Art17_habitats_$Version.csv")
  countries = @(Import-Csv (Join-Path $exportPath "Overall_statistics_MS.csv"))
  statistics = @(Import-Csv (Join-Path $exportPath "Overall_statistics.csv"))
  areas = @(Import-Csv (Join-Path $exportPath "Overall_areastatistics.csv"))
  targets = @(Import-TargetCsvs "Overall_targets*_MS_$Version.csv")
}
$json = ConvertTo-Json -InputObject $datasets -Depth 10 -Compress
$jsonVersion = ConvertTo-Json -InputObject $Version -Compress
$outputPath = Join-Path $PSScriptRoot "dashboard-data.js"
[System.IO.File]::WriteAllText($outputPath, "window.MTA_VERSION = $jsonVersion;`nwindow.MTA_DATA = $json;`n", [System.Text.UTF8Encoding]::new($false))