#!/bin/sh
# Packs the primerpicker-gui.app bundle into a disk image: the Qt
# frameworks and plugins are copied into the bundle (macdeployqt), the
# bundle is signed ad hoc (macdeployqt invalidates the signatures that
# the linker adds, and Apple silicon Macs refuse unsigned code), and
# the image holds PrimerPicker.app, a link to /Applications and the
# given documents. Without an Apple Developer ID, the program is
# neither signed by a known developer nor notarized (proposal §3.8).
#
# usage: QT_BIN=... sh packaging/macos/build-dmg.sh APP OUTPUT.dmg VOLUME_NAME [FILE...]
set -eu

app=$1
output=$2
volume=$3
shift 3

staging=$(mktemp -d)
bundle=$staging/PrimerPicker.app
cp -R "$app" "$bundle"
# QT_BIN: the bin folder of Qt (macdeployqt; the software bill of
# materials is in ../sbom)
qt_bin=${QT_BIN:?the bin folder of Qt}
deploy=$qt_bin/macdeployqt
"$deploy" "$bundle" -verbose=1
# the licences (packaging/licenses/), with the software bill of
# materials of Qt, in the bundle (before it is signed)
here=$(cd "$(dirname "$0")" && pwd)
licences=$bundle/Contents/Resources/licenses
mkdir -p "$licences/qt-sbom"
cp "$here/../../LICENSE" "$here/../licenses/THIRD-PARTY.md" "$here/../licenses/LGPL-3.0.txt" "$licences/"
cp "$qt_bin"/../sbom/*.spdx "$licences/qt-sbom/"
codesign --force --deep --sign - "$bundle"
codesign --verify --deep --strict "$bundle"

ln -s /Applications "$staging/Applications"
for file in "$@"; do
    cp "$file" "$staging/"
done
# no AppleDouble (._*) files in the image
COPYFILE_DISABLE=1 hdiutil create -volname "$volume" -srcfolder "$staging" -format UDZO -ov "$output"
rm -rf "$staging"
