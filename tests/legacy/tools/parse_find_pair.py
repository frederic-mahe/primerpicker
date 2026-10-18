"""Turn a recorded Find Primer Pair list into a table.

Usage: python3 -I parse_find_pair.py WEIGHTS RUN_DIR
  WEIGHTS: the five weights of the run, for instance "0 5 0 0 1";
  RUN_DIR holds screens.txt (a recording made by find_pair.sh).

The list shows 14 pairs at a time and scrolls while the recording
moves down: the list is rebuilt from the screens, from which pair is
shown directly below which (the highlighted pair is not decoded, so
each pair must appear unhighlighted on some screen).
Prints one line per pair: the weights, the rank (order on screen),
then the values shown: positions of A and B, GC % of A and B, Tm of
A and B, max 3' dimers, max dimers (any), false-priming Tm of A and
B ("-" when none). Prints the counts of possible primers and pairs
on stderr. Fails when the weights form of the recording shows other
weights than WEIGHTS (a key was lost); the last weight cannot be
checked this way, as the form closes before showing it.
"""
import itertools
import re
from collections import Counter, defaultdict
import sys
from pathlib import Path

ROW = re.compile(r'│\s+(\d+)\s+(\d+)\s+(\d+) / (\d+)\s+(\d+) / (\d+)\s+(\d+)\s+(\d+)\s+(\d+|-) /\s+(\d+|-)\s')
WEIGHT_LABELS = ('GC Mismatch', 'Tm mismatch', "3'dimers", 'Any dimers', 'False Priming')
WEIGHT = re.compile(r"(GC Mismatch|Tm mismatch|3'dimers|Any dimers|False Priming) penalty:\s+(\S)")
FRAME = re.compile(r'f(\d+)\.png')
LIST_LINES = 14  # pairs shown at a time
COUNTS = re.compile(r'Primer A:\s+(\d+) possible.*?Primer B:\s+(\d+) possible.*?Pairs:\s+(\d+) possible', re.S)


def list_screens(screens):
    """Each list screen as its lines: a pair, or None (a blank line, or the
    highlighted pair, which screen_text.py does not decode). Screens shown
    for a single video frame are left out: they are captured while the
    list is redrawn, half scrolled."""
    blocks = [block.split('\n', 1) for block in screens.split('==== ')[1:]]
    numbers = [int(FRAME.match(header).group(1)) for header, _ in blocks]
    for index, (_, body) in enumerate(blocks):
        if index + 1 < len(blocks) and numbers[index + 1] - numbers[index] == 1:
            continue
        if 'FIND PRIMER PAIR' not in body:
            continue
        lines = body.split('\n')
        first = next((number for number, line in enumerate(lines) if '-A-' in line), None)
        if first is None:
            continue
        rows = [match.groups() if match is not None else None
                for match in map(ROW.search, lines[first + 2:first + 2 + LIST_LINES])]
        if any(row is not None for row in rows):
            yield rows


def merge(screens):
    """The list, from which pair is shown directly below which"""
    votes = defaultdict(Counter)
    for rows in screens:
        for above, below in itertools.pairwise(rows):
            if above is not None and below is not None:
                votes[above][below] += 1
    pairs = [screens[0][0]]
    while pairs[-1] in votes:
        ranked = votes[pairs[-1]].most_common()
        if len(ranked) > 1:
            raise ValueError(f'screens disagree on the pair after rank {len(pairs)}: {ranked}')
        if ranked[0][0] in pairs:
            raise ValueError(f'pair listed twice: {ranked[0][0]}')
        pairs.append(ranked[0][0])
    return pairs


def shown_weights(screens):
    """The weights on the last screen of the weights form that shows all five"""
    shown = None
    for block in screens.split('==== '):
        found = dict(WEIGHT.findall(block))
        if len(found) == len(WEIGHT_LABELS):
            shown = ' '.join(found[label] for label in WEIGHT_LABELS)
    return shown


def main(weights, run_dir):
    screens = (Path(run_dir) / 'screens.txt').read_text(encoding='utf-8')
    counts = COUNTS.findall(screens)
    if counts:
        print('possible A, B, pairs:', *counts[-1], file=sys.stderr)  # earlier ones may be half drawn
    shown = shown_weights(screens)
    # the form closes when the last digit is typed, before showing it
    if shown is not None and shown.split()[:-1] != weights.split()[:-1]:
        raise SystemExit(f'weights shown: {shown}, expected: {weights}')
    lists = list(list_screens(screens))
    if not lists:
        raise SystemExit('no Find Primer Pair list in ' + run_dir)
    for rank, row in enumerate(merge(lists), start=1):
        print(weights, rank, *row, sep='\t')


if __name__ == '__main__':
    main(*sys.argv[1:])
