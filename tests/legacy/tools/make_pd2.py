"""Write a Direct Entry primer file (.PD2) of Primer Designer 2.0.

Usage: python3 -I make_pd2.py PRIMER [EXTENSION] > NAME.PD2

The layout is the one of tests/legacy/README.md ("Primer files .PD2"),
with the built-in criteria. Loaded with Primer/Create -> Retrieve from
disk, such files put primers in the Active Primer List with fewer keys
than Direct Entry (the autotype command line holds about 250
characters). Check: `make_pd2.py GTAAAACGACGGCCAGT GAATTC` writes
primers/M13EXT.PD2, byte for byte.
"""
import struct
import sys

HEADER = b'\x1aSES' + bytes([0x03, 0x01, 0, 0, 0x05, 0, 0, 0, 0, 0])
TEXT, NUMBER = 8, 9
EXTENSION, DIRECT_ENTRY = 0, 1
LENGTH, LENGTH_WITHOUT_EXTENSION, EXTENSION_LENGTH, POSITION, COMPLEMENT = 0, 1, 2, 3, 4
GC, TM, HAIRPIN, DIMERS, RUNS = 0x10, 0x11, 0x12, 0x13, 0x14
END = b'\xff\xff' + bytes(4)
DIRECT_ENTRY_POSITION = 1


def tag(group, field):
    return bytes([field, group])


def text(field, value):
    data = value.encode('ascii') + b'\0'
    if len(data) % 2:
        data += b'\0'
    return tag(TEXT, field) + struct.pack('<I', len(data)) + data


def number(field, packed):
    return tag(NUMBER, field) + bytes(4) + packed


def pd2(primer, extension=''):
    records = [text(EXTENSION, extension)] if extension else []
    records.append(text(DIRECT_ENTRY, primer))
    records += [
        number(LENGTH, struct.pack('<i', len(extension) + len(primer))),
        number(LENGTH_WITHOUT_EXTENSION, struct.pack('<i', len(primer))),
        number(EXTENSION_LENGTH, struct.pack('<i', len(extension))),
        number(POSITION, struct.pack('<i', DIRECT_ENTRY_POSITION)),
        number(COMPLEMENT, struct.pack('<i', 0)),
        # built-in criteria: GC max, min; Tm max, min; hairpin cutoff;
        # adjacent bases, 3' matches; the two run thresholds
        number(GC, struct.pack('<hh', 60, 50)),
        number(TM, struct.pack('<hh', 80, 55)),
        number(HAIRPIN, struct.pack('<f', 0.0)),
        number(DIMERS, struct.pack('<hh', 7, 3)),
        number(RUNS, struct.pack('<hh', 3, 3)),
    ]
    return HEADER + b''.join(records) + END


if __name__ == '__main__':
    sys.stdout.buffer.write(pd2(sys.argv[1].upper(), sys.argv[2].upper() if len(sys.argv) > 2 else ''))
