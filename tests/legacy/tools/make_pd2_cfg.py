"""Write a settings file (pd2.cfg) of Primer Designer 2.0 for probes.

Usage: python3 -I make_pd2_cfg.py [--length N] [--exclude-tm Y|N]
           [--sort-tm Y|N] [--weights DIGITS] [--criteria VALUES] > PD2.CFG

The built-in defaults (tests/legacy/README.md, "Built-in defaults"),
with the layout of "Settings file pd2.cfg", and the changes asked for.
WEIGHTS: five digits in the order of the weights form (GC mismatch, Tm
mismatch, 3' dimers, any dimers, false priming); the file stores the
false-priming weight third. Checks: `--length 25` writes
config/PD2_LENGTH25.CFG and `--weights 12345` writes
config/PD2_WEIGHTS12345.CFG, byte for byte.
CRITERIA: nine comma-separated values in the order of the Criteria form
(GC min, GC max, Tm min, Tm max, hairpin cutoff in kcal, 3' dimers,
adjacent bases, runs, 3' G/C runs); default 50,60,55,80,0,3,7,3,3.
`--criteria 0,100,0,99,-99,99,99,99,99` writes config/PD2_RELAXED.CFG.
"""
import argparse
import struct
import sys

UNKNOWN = bytes(4)  # zero in every recorded file
TRAILER = bytes(8)  # offsets 47 to 54, zero in every recorded file


DEFAULT_CRITERIA = '50,60,55,80,0,3,7,3,3'


def settings(length=20, exclude_tm=True, sort_tm=True, weights='11111', criteria=DEFAULT_CRITERIA):
    gc_min, gc_max, tm_min, tm_max, hairpin, dimer_3prime_matches, adjacent, runs, gc_runs = criteria.split(',')
    gc_weight, tm_weight, dimer_3prime, dimer_any, false_priming = (int(digit) for digit in weights)
    return (UNKNOWN
            + struct.pack('<6h', length, int(exclude_tm), int(sort_tm), 5, 10, 50)
            + struct.pack('<5h', 50, 100, 100, int(gc_min), int(gc_max))  # Match end; Find Pair ranges; GC
            + struct.pack('<2h', int(tm_min), int(tm_max))                 # Tm min, max
            + struct.pack('<f', float(hairpin))                            # hairpin cutoff
            + struct.pack('<4h', int(dimer_3prime_matches), int(adjacent), int(runs), int(gc_runs))
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
    parser.add_argument('--criteria', default=DEFAULT_CRITERIA)
    arguments = parser.parse_args()
    sys.stdout.buffer.write(settings(arguments.length, arguments.exclude_tm, arguments.sort_tm, arguments.weights,
                                     arguments.criteria))


if __name__ == '__main__':
    main()
