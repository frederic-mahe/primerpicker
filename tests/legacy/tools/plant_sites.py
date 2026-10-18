"""Write a random DNA molecule with planted copies of some of its parts.

The molecule is made by make_molecule.py (equal segments, each with
its own GC fraction); each COPY then overwrites part of it with a copy
of another part, for instance the 3' end of a primer, to create
false-priming sites, or with a given sequence, for instance a
restriction site.

Usage: python3 -I plant_sites.py NAME SEED GC_FRACTIONS SEGMENT_LENGTH COPY... > FILE
  COPY: FIRST:LENGTH:DESTINATION (1-based positions on the top strand)
  copies LENGTH bases from FIRST to DESTINATION; SEQUENCE@DESTINATION
  writes SEQUENCE (A, C, G, T) at DESTINATION; copies are made in the
  order given
  for instance: plant_sites.py TFP 3236 0.55,0.55,0.55,0.55,0.55,0.55 150 88:11:620 GAATTC@300
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python3 -I leaves it out
import make_molecule  # noqa: E402


def plant(sequence, copies):
    bases = list(sequence)
    for source, destination in copies:
        if isinstance(source, str):
            planted = list(source.lower())
        else:
            first, length = source
            planted = bases[first - 1:first - 1 + length]
        bases[destination - 1:destination - 1 + len(planted)] = planted
    return ''.join(bases)


def parse_copy(spec):
    if '@' in spec:
        sequence, destination = spec.split('@')
        if not sequence or set(sequence) - set('ACGT'):
            raise SystemExit(f'not a sequence of A, C, G and T: {spec}')
        return sequence, int(destination)
    first, length, destination = (int(value) for value in spec.split(':'))
    return (first, length), destination


def main(name, seed, gc_fractions, segment_length, *copy_specs):
    gcs = [float(gc) for gc in gc_fractions.split(',')]
    copies = [parse_copy(spec) for spec in copy_specs]
    sequence = plant(make_molecule.make_sequence(int(seed), gcs, int(segment_length)), copies)
    description = f'Random test molecule (seed {seed}, GC {gc_fractions}), copies {" ".join(copy_specs)}.'
    sys.stdout.buffer.write(make_molecule.genbank(name, description, sequence).encode('ascii'))


if __name__ == '__main__':
    main(*sys.argv[1:])
