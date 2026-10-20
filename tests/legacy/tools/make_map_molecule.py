"""Write a random molecule with copies of a stretch of one of its primers.

For probes of the molecule map of the Uniqueness and False Priming
reports: the sites must sit at known positions, on known strands.

Usage: python3 -I make_map_molecule.py NAME SEED LENGTH FIRST STRETCH SITE... > FILE
  LENGTH bases, GC 50 %; FIRST:LAST (1-based, inclusive) is the stretch
  of the molecule to copy (for instance 55:66, inside the primer at 51);
  each SITE is DESTINATION followed by t (a copy on the top strand,
  written from DESTINATION on) or c (its reverse complement, written
  from DESTINATION on: a site on the complement strand)
  for instance: make_map_molecule.py TMAP660 6601 660 55:66 100t 171c
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python3 -I leaves it out
import make_molecule  # noqa: E402

COMPLEMENT = str.maketrans('acgt', 'tgca')


def plant(sequence, stretch, sites):
    bases = list(sequence)
    first, last = stretch
    copied = sequence[first - 1:last]
    for destination, strand in sites:
        planted = copied if strand == 't' else copied.translate(COMPLEMENT)[::-1]
        if destination - 1 + len(planted) > len(bases):
            raise SystemExit(f'site {destination}{strand} goes beyond the molecule')
        bases[destination - 1:destination - 1 + len(planted)] = planted
    return ''.join(bases)


def parse_site(spec):
    if len(spec) < 2 or spec[-1] not in 'tc' or not spec[:-1].isdigit():
        raise SystemExit(f'not a site (DESTINATION then t or c): {spec}')
    return int(spec[:-1]), spec[-1]


def main(name, seed, length, stretch, *site_specs):
    first, last = (int(value) for value in stretch.split(':'))
    sites = [parse_site(spec) for spec in site_specs]
    sequence = plant(make_molecule.make_sequence(int(seed), [0.5], int(length)), (first, last), sites)
    description = f'Random test molecule (seed {seed}), copies of {stretch} at {" ".join(site_specs)}.'
    sys.stdout.buffer.write(make_molecule.genbank(name, description, sequence).encode('ascii'))


if __name__ == '__main__':
    main(*sys.argv[1:])
