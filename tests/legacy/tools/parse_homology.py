"""Turn recorded Uniqueness and False Priming screens into a table.

Usage: python3 -I parse_homology.py MOLECULE PRIMER POSITION STRAND KIND RUN_DIR
  MOLECULE, PRIMER, POSITION, STRAND (top or complement) describe the
  analysed primer; KIND is "uniqueness" or "false_priming"; RUN_DIR
  holds screens.txt (a recording that ends on that screen).

Prints one line per listed site: molecule, primer, position, strand,
kind, the homologous bases as shown (without the dots), the site
position (BpPos), the site strand (top, or complement for "(C)"), and
the Tm; or one line with "-" fields when the program reports no site
("Primer is unique").
"""
import re
import sys
from pathlib import Path

SITE = re.compile(r"5'\.*([ACGT]+)(?:\.*3')?\s+(\d+)( \(C\))?\s+(-?\d+)")


def sites(screens):
    for block in screens.split('==== '):
        if 'BpPos' not in block:
            continue
        found = [SITE.search(line) for line in block.split('\n')]
        found = [match.groups() for match in found if match is not None]
        if found:
            return found
    return []


def main(molecule, primer, position, strand, kind, run_dir):
    screens = (Path(run_dir) / 'screens.txt').read_text(encoding='utf-8')
    found = sites(screens)
    prefix = [molecule, primer, position, strand, kind]
    if not found:
        if 'is unique' not in screens and 'No ' not in screens:
            print(f'{run_dir}: no site and no "unique" message', file=sys.stderr)
            return 1
        print('\t'.join(prefix + ['-', '-', '-', '-']))
        return 0
    for bases, site_position, complement, tm in found:
        print('\t'.join(prefix + [bases, site_position, 'complement' if complement else 'top', tm]))
    return 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
