$ErrorActionPreference='Stop'
Set-Location $PSScriptRoot
if (-not (Test-Path -LiteralPath 'runtime/guidewire_core.dll')) {
    throw 'Build the numerical core first: ./build_native.ps1 -FetchToolchain'
}
python -m PyInstaller --noconfirm --clean --onedir --windowed --name GuidewireSimulator --paths viewer --hidden-import native_engine --hidden-import xray_viewer --hidden-import native_core --hidden-import reference_diagnostics --hidden-import scene_parameters --hidden-import PIL.ImageTk --hidden-import PIL.ImageGrab --exclude-module matplotlib --exclude-module scipy portable_launcher.py
if($LASTEXITCODE -ne 0){throw 'Executable build failed'}
$out=Join-Path $PSScriptRoot 'dist/GuidewireSimulator'
foreach($name in @('assets','config')){Copy-Item -LiteralPath $name -Destination $out -Recurse -Force}
New-Item -ItemType Directory -Path (Join-Path $out 'runtime') -Force | Out-Null
Copy-Item -LiteralPath 'runtime/guidewire_core.dll' -Destination (Join-Path $out 'runtime/guidewire_core.dll') -Force
Copy-Item -LiteralPath 'START.bat' -Destination $out -Force
Write-Host 'Copy source, third-party licenses and notices alongside the binary when distributing.'
