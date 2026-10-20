#!/bin/sh
# Gathers primerpicker-gui.exe and everything it needs into a folder
# that runs on a Windows computer without MSYS2: the Qt libraries and
# plugins (windeployqt), then every DLL of the MSYS2 environment that
# the executables and DLLs of the folder load, until none is missing.
#
# usage (MSYS2 UCRT64 shell): sh packaging/windows/deploy-gui.sh EXE FOLDER
set -eu

exe=$1
folder=$2
prefix=${MINGW_PREFIX:?run from an MSYS2 shell (UCRT64)}

mkdir -p "$folder"
cp "$exe" "$folder/"
# MSYS2 names it windeployqt6 (or windeployqt-qt6 in older packages)
deploy=$(command -v windeployqt6 || command -v windeployqt-qt6 || command -v windeployqt)
"$deploy" --release --no-translations --no-system-d3d-compiler --no-opengl-sw \
    "$folder/$(basename "$exe")"

# MSYS2 DLLs (under $prefix/bin) loaded by the files of the folder
missing_dlls() {
    find "$folder" \( -name '*.exe' -o -name '*.dll' \) -print0 | xargs -0 ldd 2> /dev/null |
        awk -v prefix="$prefix/bin/" 'index($3, prefix) == 1 { print $3 }' | sort -u |
        while read -r dll; do
            [ -e "$folder/$(basename "$dll")" ] || printf '%s\n' "$dll"
        done
}

while :; do
    dlls=$(missing_dlls)
    [ -z "$dlls" ] && break
    printf '%s\n' "$dlls" | while read -r dll; do
        echo "adding $(basename "$dll")"
        cp "$dll" "$folder/"
    done
done
