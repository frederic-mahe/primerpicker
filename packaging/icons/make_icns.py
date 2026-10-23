"""Writes a macOS icon file (.icns) from PNG images: a header, then one
entry (4-byte type, 4-byte big-endian length, PNG data) per size. The
types are those of PNG entries (macOS 10.7 and later); a "@2x" entry
holds twice as many pixels as its size in points.

usage: python3 make_icns.py OUTPUT.icns PNG_DIR
(PNG_DIR holds primerpicker-SIZE.png, SIZE in pixels)
"""
import struct
import sys
from pathlib import Path

# type -> width in pixels
ENTRIES = {
    "icp4": 16, "icp5": 32, "icp6": 64, "ic07": 128, "ic08": 256,
    "ic09": 512, "ic10": 1024, "ic11": 32, "ic12": 64, "ic13": 256,
    "ic14": 512,
}
# the 16 and 32 point sizes, also at @2x, show the small icon
SMALL_POINTS_ENTRIES = {"icp4", "icp5", "ic11", "ic12"}


def main():
    output = Path(sys.argv[1])
    pngs = Path(sys.argv[2])
    body = b""
    for kind, pixels in ENTRIES.items():
        name = f"primerpicker-small-{pixels}.png" if kind in SMALL_POINTS_ENTRIES else f"primerpicker-{pixels}.png"
        data = (pngs / name).read_bytes()
        body += kind.encode("ascii") + struct.pack(">I", 8 + len(data)) + data
    output.write_bytes(b"icns" + struct.pack(">I", 8 + len(body)) + body)


if __name__ == "__main__":
    main()
