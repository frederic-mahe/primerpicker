#!/bin/sh
# Checks a primerpicker-gui AppImage on the system that runs this
# script: it reports the expected version, and the program starts and
# keeps running (a missing library or plugin stops it at once). Needs
# an X server (DISPLAY; xvfb-run on a server).
#
# usage: sh packaging/linux/check-appimage.sh APPIMAGE VERSION
set -eu

appimage=$1
version=$2
seconds=10

# no FUSE in containers: the AppImage is extracted at each run
export APPIMAGE_EXTRACT_AND_RUN=1

reported=$("$appimage" --version)
if [ "$reported" != "primerpicker $version" ]; then
    echo "error: $appimage reports '$reported', not 'primerpicker $version'" >&2
    exit 1
fi

"$appimage" &
pid=$!
sleep "$seconds"
if ! kill -0 "$pid" 2> /dev/null; then
    wait "$pid" || status=$?
    echo "error: $appimage stopped within $seconds s (exit status ${status:-0})" >&2
    exit 1
fi
kill "$pid"
echo "$appimage: version $version, still running after $seconds s"
