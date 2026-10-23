#!/bin/sh
# Packs primerpicker-gui and the Qt libraries and plugins it needs into
# an AppImage (one executable file), with linuxdeploy and its Qt
# plugin. The tools are downloaded once, pinned by release and
# checksum, into TOOLS (default: ./appimage-tools). Qt is found through
# its qmake: QMAKE, or the qmake of PATH.
#
# Build on the oldest supported system (the glibc floor, proposal
# §3.8): the AppImage runs on systems with that glibc or a newer one.
#
# usage: sh packaging/linux/build-appimage.sh EXE OUTPUT.AppImage
set -eu

exe=$1
output=$2
here=$(cd "$(dirname "$0")" && pwd)
tools=${TOOLS:-$PWD/appimage-tools}

linuxdeploy_url=https://github.com/linuxdeploy/linuxdeploy/releases/download/1-alpha-20251107-1/linuxdeploy-x86_64.AppImage
linuxdeploy_sha256=c20cd71e3a4e3b80c3483cef793cda3f4e990aca14014d23c544ca3ce1270b4d
plugin_qt_url=https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/1-alpha-20250213-1/linuxdeploy-plugin-qt-x86_64.AppImage
plugin_qt_sha256=15106be885c1c48a021198e7e1e9a48ce9d02a86dd0a1848f00bdbf3c1c92724
# the AppImage runtime (the start of the AppImage file); without it,
# linuxdeploy downloads the runtime of the day
runtime_url=https://github.com/AppImage/type2-runtime/releases/download/20251108/runtime-x86_64
runtime_sha256=2fca8b443c92510f1483a883f60061ad09b46b978b2631c807cd873a47ec260d

fetch() {
    url=$1
    sha256=$2
    file=$tools/$(basename "$url")
    [ -e "$file" ] || curl --fail --silent --show-error --location --output "$file" "$url"
    echo "$sha256  $file" | sha256sum --check --quiet
    chmod +x "$file"
}

mkdir -p "$tools"
fetch "$linuxdeploy_url" "$linuxdeploy_sha256"
fetch "$plugin_qt_url" "$plugin_qt_sha256"
fetch "$runtime_url" "$runtime_sha256"

# linuxdeploy finds its plugins in PATH; the tools are AppImages
# themselves, extracted at each run (no FUSE needed, as in containers)
PATH=$tools:$PATH
export PATH
export APPIMAGE_EXTRACT_AND_RUN=1
export QMAKE="${QMAKE:-$(command -v qmake6 || command -v qmake)}"
# Only the xcb platform plugin (X11; XWayland in Wayland sessions): the
# Qt plugin of linuxdeploy deploys the wayland platform plugin without
# the shell integration plugins that it needs to show a window

appdir=$(mktemp -d)/AppDir
LDAI_OUTPUT=$output
LDAI_RUNTIME_FILE=$tools/runtime-x86_64
export LDAI_OUTPUT LDAI_RUNTIME_FILE
linuxdeploy-x86_64.AppImage --appdir "$appdir" \
    --executable "$exe" \
    --desktop-file "$here/primerpicker.desktop" \
    --icon-file "$here/primerpicker.svg" \
    --plugin qt \
    --output appimage
rm -rf "$(dirname "$appdir")"
