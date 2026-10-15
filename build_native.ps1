param([switch]$FetchToolchain)
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$compiler = Join-Path $PSScriptRoot 'runtime/toolchain/llvm-mingw-20261006-ucrt-x86_64/bin/clang++.exe'
if (-not (Test-Path -LiteralPath $compiler)) {
    if (-not $FetchToolchain) { throw 'Portable C++ compiler missing. Run ./build_native.ps1 -FetchToolchain once.' }
    New-Item -ItemType Directory -Force -Path 'runtime/toolchain' | Out-Null
    $archive = Join-Path $PSScriptRoot 'runtime/llvm-mingw.zip'
    Invoke-WebRequest 'https://github.com/mstorsjo/llvm-mingw/releases/download/20261006/llvm-mingw-20261006-ucrt-x86_64.zip' -OutFile $archive
    Expand-Archive -LiteralPath $archive -DestinationPath 'runtime/toolchain' -Force
}
# Generated functions and upstream math headers are committed, so a clone can
# compile without downloading SOFA or using the installed SOFA DLLs.
foreach ($target in @('probe', 'kernel_probe', 'contact_probe', 'geometry_probe', 'scene_probe')) {
    & $compiler -std=c++20 -O2 -static -I native/shim -I upstream/include "native/$target.cpp" -o "runtime/$target.exe"
    if ($LASTEXITCODE -ne 0) { throw "Native build failed: $target" }
}
& $compiler -std=c++20 -O2 -static -shared -I native/shim -I upstream/include native/bridge.cpp -o runtime/guidewire_core.dll
if ($LASTEXITCODE -ne 0) { throw 'Native build failed: guidewire_core.dll' }
