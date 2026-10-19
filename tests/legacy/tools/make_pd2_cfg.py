"""Write a settings file (pd2.cfg) of Primer Designer 2.0 for probes.

Usage: python3 -I make_pd2_cfg.py [--length N] [--exclude-tm Y|N]
           [--sort-tm Y|N] [--weights DIGITS] > PD2.CFG

The built-in defaults (tests/legacy/README.md, "Built-in defaults"),
with the layout of "Settings file pd2.cfg", and the changes asked for.
WEIGHTS: five digits in the order of the weights form (GC mismatch, Tm
mismatch, 3' dimers, any dimers, false priming); the file stores the
false-priming weight third. Checks: `--length 25` writes
config/PD2_LENGTH25.CFG and `--weights 12345` writes
config/PD2_WEIGHTS12345.CFG, byte for byte.
"""
import argparse
import struct
import sys

UNKNOWN = bytes(4)  # zero in every recorded file
TRAILER = bytes(8)  # offsets 47 to 54, zero in every recorded file


def settings(length=20, exclude_tm=True, sort_tm=True, weights='11111'):
    gc_weight, tm_weight, dimer_3prime, dimer_any, false_priming = (int(digit) for digit in weights)
    return (UNKNOWN
            + struct.pack('<6h', length, int(exclude_tm), int(sort_tm), 5, 10, 50)
            + struct.pack('<5h', 50, 100, 100, 50, 60)    # Match end; Find Pair ranges; GC min, max
            + struct.pack('<2h', 55, 80)                  # Tm min, max
            + struct.pack('<f', 0.0)                      # hairpin cutoff
            + struct.pack('<4h', 3, 7, 3, 3)              # 3' dimers, adjacent bases, runs, 3' G/C runs
            + bytes([gc_weight, tm_weight, false_priming, dimer_3prime, dimer_any])
            + TRAILER)


def yes(text):
    if text not in ('Y', 'N'):
        raise argparse.ArgumentTypeError('Y or N')
    return text == 'Y'


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--length', type=int, default=20)
    parser.add_argument('--exclude-tm', type=yes, default=True)
    parser.add_argument('--sort-tm', type=yes, default=True)
    parser.add_argument('--weights', default='11111')
    arguments = parser.parse_args()
    sys.stdout.buffer.write(settings(arguments.length, arguments.exclude_tm, arguments.sort_tm, arguments.weights))


if __name__ == '__main__':
    main()
