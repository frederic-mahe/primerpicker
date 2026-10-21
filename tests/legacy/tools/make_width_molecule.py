"""Write picking/TWIDTH.GB, the molecule of the value-line probes.

300 bases: 100 of A and T with one G at position 50 (a 20-base primer
there has 5 % GC), 100 of G and C (100 % GC, Tm above 99 °C), then 100
with a GC fraction of 0.5; random choices from seed 2110.

Usage: python3 -I make_width_molecule.py > TWIDTH.GB
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python3 -I leaves it out
import make_molecule  # noqa: E402

SEED = 2110
STRETCH = 100
SINGLE_G = 50  # position of the only G of the first stretch


def sequence():
    rnd = random.Random(SEED)
    at_rich = ''.join(rnd.choice('at') for _ in range(STRETCH))
    at_rich = at_rich[:SINGLE_G - 1] + 'g' + at_rich[SINGLE_G:]
    gc_rich = ''.join(rnd.choice('gc') for _ in range(STRETCH))
    mixed = ''.join((rnd.choice('gc') if rnd.random() < 0.5 else rnd.choice('at')) for _ in range(STRETCH))
    return at_rich + gc_rich + mixed


if __name__ == '__main__':
    text = make_molecule.genbank('TWIDTH', 'value line widths: AT-rich (one G at 50), GC-rich, mixed', sequence())
    sys.stdout.buffer.write(text.encode('ascii'))
