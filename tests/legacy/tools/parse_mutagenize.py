"""Turn Mutagenize Primer recordings (mutagenize.sh) into a table.

Usage: python3 -I parse_mutagenize.py RUN_DIR...

Each RUN_DIR is named MOLECULE_POSITION[_KEYS] and holds screens.txt;
the last Mutagenize Primer screen without a popup (codon table, enzyme
list) and without the blinking cursor in the lines read is used. Prints
one row per screen line of interest, in screen order: molecule,
position, keys, line (original, O, R, I, molecule, marker, primer, M,
U, T, mutated) and the text of the line without the border column
(right-trimmed); for "original" and "mutated", the values after the
label.
"""
import re
import sys
from pathlib import Path

ROWS = ('O', 'R', 'I', 'molecule', 'marker', 'primer', 'M', 'U', 'T')
POPUP = re.compile('[║╔╚]')


BOX_RIGHT = 71  # the column of the box's right border (the amino-acid legend follows)


def text(line):
    """The line inside the box: without the border columns, right-trimmed."""
    return line[1:BOX_RIGHT].rstrip()


def screen_rows(lines):
    original = next(index for index, line in enumerate(lines) if 'ORIGINAL:' in line)
    mutated = next(index for index, line in enumerate(lines) if 'MUTATED:' in line)
    rows = {'original': text(lines[original]).split('ORIGINAL:', 1)[1].strip(),
            'mutated': text(lines[mutated]).split('MUTATED:', 1)[1].strip()}
    # the frames start three lines after ORIGINAL; nine lines in all
    for offset, name in enumerate(ROWS):
        rows[name] = text(lines[original + 3 + offset])
    return rows


def usable(lines):
    if not any('MUTAGENIZE PRIMER' in line for line in lines):
        return False
    if not any('ORIGINAL:' in line for line in lines) or not any('MUTATED:' in line for line in lines):
        return False
    body = [line[:72] for line in lines]
    return not any(POPUP.search(line) or '?' in line[1:71] for line in body[3:22])


def main(run_dirs):
    print('\t'.join(['# molecule', 'position', 'keys', 'line', 'text']))
    for run_dir in run_dirs:
        molecule, position, *keys = Path(run_dir).name.split('_')
        screens = [screen.splitlines() for screen in (Path(run_dir) / 'screens.txt').read_text().split('==== ')]
        candidates = [lines for lines in screens if usable(lines)]
        if not candidates:
            raise SystemExit(f'{run_dir}: no usable Mutagenize screen')
        rows = screen_rows(candidates[-1])
        for name in ('original', *ROWS, 'mutated'):
            print('\t'.join([molecule, position, '_'.join(keys) or '-', name, rows[name]]))


if __name__ == '__main__':
    main(sys.argv[1:])
