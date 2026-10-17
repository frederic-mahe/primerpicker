"""Write a random DNA molecule as a GenBank file that PRIMER.exe loads.

The sequence is made of equal segments, each with its own GC fraction.

Usage: python3 -I make_molecule.py NAME SEED GC_FRACTIONS SEGMENT_LENGTH > FILE
  for instance: make_molecule.py TGRAD 1201 0.25,0.4,0.55,0.7,0.55,0.4,0.25,0.8 150

PRIMER.exe lists files in alphabetical order, and tools/pair_product.sh
selects the file listed after PRIMER.EXE: name the file accordingly
(TEST600.GB, TGRAD.GB).
"""
import random
import sys

BASES_PER_LINE = 60
BASES_PER_BLOCK = 10


def make_sequence(seed, gc_fractions, segment_length):
    rnd = random.Random(seed)
    return ''.join((rnd.choice('gc') if rnd.random() < gc else rnd.choice('at'))
                   for gc in gc_fractions for _ in range(segment_length))


def genbank(name, description, sequence):
    lines = [f'LOCUS       {name:<12s} {len(sequence)} bp    DNA     linear',
             f'DEFINITION  {description}',
             'FEATURES             Location/Qualifiers',
             'ORIGIN']
    for start in range(0, len(sequence), BASES_PER_LINE):
        end = min(start + BASES_PER_LINE, len(sequence))
        blocks = [sequence[i:i + BASES_PER_BLOCK] for i in range(start, end, BASES_PER_BLOCK)]
        lines.append(f'{start + 1:9d} ' + ' '.join(blocks))
    lines.append('//')
    return '\r\n'.join(lines) + '\r\n'  # DOS line ends


def main(name, seed, gc_fractions, segment_length):
    gcs = [float(gc) for gc in gc_fractions.split(',')]
    sequence = make_sequence(int(seed), gcs, int(segment_length))
    description = f'Random test molecule (seed {seed}, GC {gc_fractions}).'
    sys.stdout.buffer.write(genbank(name, description, sequence).encode('ascii'))


if __name__ == '__main__':
    main(*sys.argv[1:])
