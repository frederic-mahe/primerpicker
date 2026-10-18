"""Turn recorded Primer Pair Analysis screens into a table.

Usage: python3 -I parse_pair.py RUN_DIR... > table.tsv
(each RUN_DIR holds screens.txt; see pair_product.sh and
pair_direct.sh)

Columns: both primers; for A and B, length, GC % and Tm (with "!"
flags as shown) and the Comment column of each balance line; the
cautions (G/C runs at the 3' end, runs, hairpins: value and flag, and
the comment); the dimer matrix (3' end matches and adjacent
homologous bases for A:A, B:B and A:B, with flags).
"""
import re
import sys
from pathlib import Path

COLUMNS = ('primer_a', 'primer_b',
           'length_a', 'length_b', 'length_comment',
           'gc_a', 'gc_b', 'gc_comment',
           'tm_a', 'tm_b', 'tm_comment',
           'gc_run_3prime_a', 'gc_run_3prime_b', 'gc_run_comment',
           'runs_a', 'runs_b', 'runs_comment',
           'hairpins_a', 'hairpins_b', 'hairpins_comment',
           'dimer_3prime_aa', 'dimer_3prime_bb', 'dimer_3prime_ab',
           'dimer_any_aa', 'dimer_any_bb', 'dimer_any_ab')
FRAME = '│'
VALUE = r'(\S+)'


def inner(line):
    return line.split(FRAME)[1].rstrip() if line.count(FRAME) >= 2 else line.rstrip()


def balance_line(lines, label):
    """Two values and the comment of a balance or caution line."""
    for line in lines:
        text = inner(line)
        if text.strip().startswith(label):
            rest = text.strip()[len(label):]
            match = re.match(r'\s*' + VALUE + r'\s+' + VALUE + r'\s*(.*)$', rest)
            if match and '?' not in match.group(0):
                return match.group(1), match.group(2), match.group(3).strip() or '-'
    raise ValueError(label)


def dimer_line(lines, label):
    for line in lines:
        match = re.search(re.escape(label) + r'\s+' + VALUE + r'\s+' + VALUE + r'\s+' + VALUE, inner(line))
        if match and '?' not in match.group(0):
            return match.groups()
    raise ValueError(label)


def primers(lines):
    found = {}
    for line in lines:
        match = re.search(r"Primer ([AB]) \(bp.*?5'-([ACGT]+)-3'", line)
        if match:
            found[match.group(1)] = match.group(2)
    return found['A'], found['B']


def parse(screens):
    for block in screens.split('==== '):
        if '* * *   PRIMER PAIR ANALYSIS  * * *' not in block or 'Any adj homol bases' not in block:
            continue
        lines = block.split('\n')
        try:
            values = list(primers(lines))
            for label in ('Length', '% GC', 'Tm°C', "Runs of G/C at 3' end", 'Runs of bases (any)', 'Hairpins'):
                values += balance_line(lines, label)
            three_prime = dimer_line(lines, "3' end matches")
            adjacent = dimer_line(lines, 'Any adj homol bases')
            return values + list(three_prime) + list(adjacent)
        except (ValueError, KeyError):
            continue
    raise ValueError('no complete Primer Pair Analysis screen')


def main(run_dirs):
    print('# ' + '\t'.join(COLUMNS))
    seen = set()
    for run_dir in run_dirs:
        path = Path(run_dir) / 'screens.txt'
        try:
            values = parse(path.read_text(encoding='utf-8'))
        except (OSError, ValueError) as error:
            print(f'{run_dir}: {error}', file=sys.stderr)
            continue
        if tuple(values[:2]) in seen:
            continue
        seen.add(tuple(values[:2]))
        print('\t'.join(values))


if __name__ == '__main__':
    main(sorted(sys.argv[1:]))
