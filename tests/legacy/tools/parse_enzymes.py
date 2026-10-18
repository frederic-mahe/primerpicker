"""Turn a recorded True Sites or Almost Sites screen into a table.

Usage: python3 -I parse_enzymes.py MOLECULE RUN_DIR
  RUN_DIR: a recording made by find_enzyme.sh (screens.txt)

Prints one line per listed site: molecule, kind (true or almost),
enzyme, recognition sequence, rank (order on screen), position, strand
(top, or complement for "C" or parentheses) and, for almost sites, the
match as shown (a dot per matching base, the molecule's base where it
differs; complement sites are compared with the reverse complement of
the recognition sequence). Checks that the number of lines equals the
"Number of Sites" shown.
"""
import re
import sys
from pathlib import Path

HEADER = re.compile(r'Search for: (\S+) +([A-Z]+)(?: +\(([A-Z]+)\))? +Number of Sites = (\d+)')
TRUE_SITE = re.compile(r'(\d+)(C?)')  # "C": a site on the complement strand
ALMOST_SITE = re.compile(r'(\d+) +(\()?([.ACGT]+)\)?')
SITE_LINE = re.compile(r'│ +\d')


def last_result_screen(screens):
    blocks = [block for block in screens.split('==== ') if HEADER.search(block)]
    if not blocks:
        raise SystemExit('no enzyme site screen')
    return blocks[-1]


def site_lines(block):
    """The lines between the "Search for" header and the molecule map"""
    lines = block.split('\n')
    first = next(index for index, line in enumerate(lines) if HEADER.search(line))
    for line in lines[first + 1:]:
        if line.startswith('┌'):  # the molecule map
            return
        if SITE_LINE.match(line):
            yield line.strip('│░▓ ')


def sites(block, kind):
    lines = list(site_lines(block))
    if kind == 'true':
        for line in lines:
            for position, complement in TRUE_SITE.findall(line):
                yield int(position), 'complement' if complement else 'top', '-'
    else:
        for line in lines:
            for position, parenthesis, match in ALMOST_SITE.findall(line):
                yield int(position), 'complement' if parenthesis else 'top', match


def main(molecule, run_dir):
    block = last_result_screen((Path(run_dir) / 'screens.txt').read_text(encoding='utf-8'))
    kind = 'almost' if 'ALMOST' in block else 'true'
    enzyme, site, _, count = HEADER.search(block).groups()
    listed = list(sites(block, kind))
    if len(listed) != int(count):
        raise SystemExit(f'{len(listed)} sites read, {count} shown (the list may need scrolling)')
    for rank, (position, strand, match) in enumerate(listed, start=1):
        print(molecule, kind, enzyme, site, rank, position, strand, match, sep='\t')


if __name__ == '__main__':
    main(*sys.argv[1:])
