"""Turn a recorded True Sites or Almost Sites screen into a table.

Usage: python3 -I parse_enzymes.py MOLECULE RUN_DIR
  RUN_DIR: a recording made by find_enzyme.sh (screens.txt)

Prints one line per listed site: molecule, kind (true or almost),
enzyme, recognition sequence, rank (order on screen), position, strand
(top, or complement for "C" or parentheses) and, for almost sites, the
match as shown (a dot per matching base, the molecule's base where it
differs; complement sites are compared with the reverse complement of
the recognition sequence). Long lists scroll: the list is rebuilt from
the result screens (screens shown for a single video frame are redraws
and are left out), from which site follows which. Checks that the
number of sites equals the "Number of Sites" shown and that positions
do not decrease.
"""
import itertools
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HEADER = re.compile(r'Search for: (\S+) +([A-Z]+)(?: +\(([A-Z]+)\))? +Number of Sites = (\d+)')
TRUE_SITE = re.compile(r'(\d+)(C?)')  # "C": a site on the complement strand
ALMOST_SITE = re.compile(r'(\d+) +(\()?([.ACGT]+)\)?')
SITE_LINE = re.compile(r'│ +\d')
FRAME = re.compile(r'f(\d+)\.png')


def result_screens(screens):
    """The result screens, without those shown for a single video frame"""
    blocks = [block.split('\n', 1) for block in screens.split('==== ')[1:]]
    numbers = [int(FRAME.match(header).group(1)) for header, _ in blocks]
    for index, (_, body) in enumerate(blocks):
        if index + 1 < len(blocks) and numbers[index + 1] - numbers[index] == 1:
            continue
        if HEADER.search(body):
            yield body


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


def merge(lists):
    """The whole list, from which site follows which on each screen"""
    following = defaultdict(Counter)
    for listed in lists:
        for site, next_site in itertools.pairwise(listed):
            following[site][next_site] += 1
    merged = [lists[0][0]]
    while merged[-1] in following:
        candidates = following[merged[-1]].most_common()
        if len(candidates) > 1:
            raise SystemExit(f'screens disagree on the site after {merged[-1]}: {candidates}')
        merged.append(candidates[0][0])
    return merged


def main(molecule, run_dir):
    blocks = list(result_screens((Path(run_dir) / 'screens.txt').read_text(encoding='utf-8')))
    if not blocks:
        raise SystemExit('no enzyme site screen')
    kind = 'almost' if 'ALMOST' in blocks[-1] else 'true'
    enzyme, site, _, count = HEADER.search(blocks[-1]).groups()
    lists = [listed for listed in (list(sites(block, kind)) for block in blocks) if listed]
    merged = merge(lists)
    if len(merged) != int(count):
        raise SystemExit(f'{len(merged)} sites read, {count} shown')
    if any(later[0] < earlier[0] for earlier, later in itertools.pairwise(merged)):
        raise SystemExit('positions decrease')
    for rank, (position, strand, match) in enumerate(merged, start=1):
        print(molecule, kind, enzyme, site, rank, position, strand, match, sep='\t')


if __name__ == '__main__':
    main(*sys.argv[1:])
