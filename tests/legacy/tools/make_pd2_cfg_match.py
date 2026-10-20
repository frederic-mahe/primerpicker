"""Write a settings file (pd2.cfg) with other Match Active Primer settings.

Usage: python3 -I make_pd2_cfg_match.py [--gc N] [--tm N] [--start N]
           [--end N] [--criteria VALUES] [--length N]
           [--offsets-raw A,B,C,D] > PD2.CFG

The file of tools/make_pd2_cfg.py (built-in defaults, CRITERIA as there),
with the Match Active Primer settings replaced (layout: tests/legacy/
README.md, "Settings file pd2.cfg"): offset 10 "match GC to within +/-"
(%, default 5), 12 "match Tm to within +/-" (degrees C, default 10), 14
and 16 the search range "Start: Target -" and "End: Target +" (bases,
default 50 and 50), as README.md describes them. --offsets-raw writes
the four int16 at offsets 14, 16, 18 and 20 as given (to find out which
field each one holds; it replaces --start and --end). With no option,
the output is that of make_pd2_cfg.py.

Recorded 2026-10-20 (tests/legacy/match_variations/README.txt): the
program shows the value at offset 16 as "Start: Target -" and the one
at offset 14 as "End: Target +", the reverse of README.md's table;
--start and --end follow README.md's table (--start writes offset 14).
"""
import argparse
import importlib.util
import pathlib
import struct
import sys

GC_OFFSET = 10  # int16 little-endian, then Tm, start, end


def base_settings(length, criteria):
    path = pathlib.Path(__file__).with_name('make_pd2_cfg.py')
    spec = importlib.util.spec_from_file_location('make_pd2_cfg', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.settings(length=length, criteria=criteria or module.DEFAULT_CRITERIA)


def match_settings(gc, tm, start, end, length=20, criteria=None, raw=None):
    data = bytearray(base_settings(length, criteria))
    data[GC_OFFSET:GC_OFFSET + 8] = struct.pack('<4h', gc, tm, start, end)
    if raw is not None:
        data[GC_OFFSET + 4:GC_OFFSET + 12] = struct.pack('<4h', *raw)
    return bytes(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--gc', type=int, default=5)
    parser.add_argument('--tm', type=int, default=10)
    parser.add_argument('--start', type=int, default=50)
    parser.add_argument('--end', type=int, default=50)
    parser.add_argument('--length', type=int, default=20)
    parser.add_argument('--criteria', default=None)
    parser.add_argument('--offsets-raw', default=None)
    arguments = parser.parse_args()
    raw = None
    if arguments.offsets_raw is not None:
        raw = [int(value) for value in arguments.offsets_raw.split(',')]
        if len(raw) != 4:
            parser.error('--offsets-raw needs four values')
    sys.stdout.buffer.write(match_settings(arguments.gc, arguments.tm, arguments.start, arguments.end,
                                           arguments.length, arguments.criteria, raw))


if __name__ == '__main__':
    main()
