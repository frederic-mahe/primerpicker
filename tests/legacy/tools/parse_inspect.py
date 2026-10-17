"""Turn screens recorded by inspect_primer.sh into a table.

Usage: python3 -I parse_inspect.py RUN_DIR... > table.tsv
(each RUN_DIR is named after its primer and holds screens.txt)

Columns: the primer; the "Actual" column of the Primer Analysis
screen (3' dimers, any dimers, runs, 3' GC runs) and its hairpin
verdict (YES = meets the criterion, no hairpin); then the Hairpins
screen: the energy (kcal) and the three lines of the drawing (5' part,
pairing marks, 3' part reversed), with their leading spaces kept so
that columns line up, or "-" when no hairpin is found.
"""
import re
import sys
from pathlib import Path

COLUMNS = ('primer', 'dimers_3prime', 'dimers_any', 'runs', 'gc_run_3prime',
           'meets_no_hairpins', 'hairpin_kcal', 'hairpin_5prime', 'hairpin_pairs', 'hairpin_3prime')
FRAME = '│'


def blocks(screens):
    return screens.split('==== ')


def actual(block, label):
    match = re.search(re.escape(label) + r'.*?\s(-?\d+|--)\s+(YES|NO)\s', block)
    if match is None:
        raise ValueError(label)
    return match.group(1), match.group(2)


def analysis(screens):
    for block in blocks(screens):
        if 'PRIMER ANALYSIS' not in block or "No 3' GC runs" not in block:
            continue
        try:
            return (actual(block, "No 3' Dimers")[0], actual(block, 'No Dimers ')[0],
                    actual(block, 'No Runs ')[0], actual(block, "No 3' GC runs")[0],
                    actual(block, 'No Hairpins')[1])
        except ValueError:
            continue
    raise ValueError('no complete Primer Analysis screen')


def inner(line):
    """A screen line without the window frame, columns kept."""
    return line.split(FRAME)[1].rstrip() if line.count(FRAME) >= 2 else line.rstrip()


def hairpin(screens):
    for block in blocks(screens):
        if 'PRIMER HAIRPINS' not in block:
            continue
        if 'No Hairpins Found' in block:
            return ('-', '-', '-', '-')
        energy = re.search(r'Hairpin Energy = (-?[\d.]+) kcals', block)
        lines = block.split('\n')
        top = next((k for k, line in enumerate(lines) if "5'-" in line), None)
        if energy is None or top is None:
            continue
        drawing = tuple(inner(line) for line in lines[top:top + 3])
        if any('?' in line for line in drawing):  # the cursor hides a cell
            continue
        return (energy.group(1),) + drawing
    raise ValueError('no complete Hairpins screen')


def main(run_dirs):
    print('# ' + '\t'.join(COLUMNS))
    for run_dir in run_dirs:
        path = Path(run_dir)
        primer = path.name.removeprefix('insp_')
        screens = (path / 'screens.txt').read_text(encoding='utf-8')
        try:
            print('\t'.join((primer,) + analysis(screens) + hairpin(screens)))
        except ValueError as error:
            print(f'{primer}: {error}', file=sys.stderr)


if __name__ == '__main__':
    main(sorted(sys.argv[1:], key=lambda name: (len(Path(name).name), Path(name).name)))
