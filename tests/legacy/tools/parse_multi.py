"""Turn Multi Primer Analysis recordings (multi_primer.sh) into tables.

Usage: python3 -I parse_multi.py balance|dimers|hairpins|duplexes RUN_DIR...

Each RUN_DIR is named SELECTION_PRIMER1_PRIMER2... and holds
screens.txt; the last Multi Primer Analysis screen is read.
- balance: one row per filled field (A to D): position (bp), primer,
  then the Primer Balance and Cautions column (values as shown, with
  their '!' flags);
- dimers: one row per cell of the upper triangle of the dimer
  matrices (row <= column): 3' end matches, then any adjacent bases;
- hairpins: from the hairpin pages (H, then PgDn), one row per primer:
  the energy shown ("-" for "No Hairpins Found") and the three lines of
  the drawing, without the border column (then in the columns of
  inspect/*.tsv);
- duplexes: from the dimer page (D, then the arrow keys), one row per
  pair shown ("Interactions between primers X and Y"): the three lines
  of its drawing, without the border column (from the first screen
  where they are drawn and hold no '?', the blinking cursor: the
  header changes after the next drawing).
"""
import re
import sys
from pathlib import Path

FIELDS = 'ABCD'
BALANCE_ROWS = ('Length', '% GC', 'Tm°C', "Runs of G/C at 3' end", 'Runs of bases (any)', 'Hairpins')
PRIMER_LINE = re.compile(r"Primer ([A-D]) \(bp +(\d+) \): +5'-([ACGT]+)-3'")
CELL_HALF_WIDTH = 3


def last_screen(path):
    screens = Path(path).read_text().split('==== ')
    matching = [screen for screen in screens if 'Primer Balance and Cautions' in screen]
    if not matching:
        raise ValueError(f'{path}: no Multi Primer Analysis screen')
    return matching[-1].splitlines()


def cell(line, column):
    """The value written around a column of the screen, '' if none."""
    return line[max(0, column - CELL_HALF_WIDTH):column + CELL_HALF_WIDTH + 1].strip(' │∙')


def balance(lines):
    header = next(line for line in lines if 'Primer Balance and Cautions' in line)
    columns = {field: header.index(f'-{field}-') + 1 for field in FIELDS}
    primers = {}
    for line in lines:
        match = PRIMER_LINE.search(line)
        if match:
            primers[match.group(1)] = (match.group(2), match.group(3))
    values = {}
    for name in BALANCE_ROWS:
        line = next(line for line in lines if line.lstrip('│ ').startswith(name))
        values[name] = {field: cell(line, columns[field]) for field in FIELDS}
    for field in FIELDS:
        if field in primers:
            yield [field, *primers[field], *(values[name][field] for name in BALANCE_ROWS)]


def dimers(lines):
    start = next(index for index, line in enumerate(lines) if 'Primer Dimers' in line)
    header_index = next(index for index in range(start, len(lines)) if re.search(r'A +B +C +D', lines[index]))
    header = lines[header_index]
    letters = [match.start() for match in re.finditer('[ABCD]', header)]
    three_prime, any_bases = letters[:4], letters[4:]
    for offset, row in enumerate(FIELDS):
        line = lines[header_index + 1 + offset]
        for index, column in enumerate(FIELDS):
            if index < offset:
                continue
            values = cell(line, three_prime[index]), cell(line, any_bases[index])
            if any(values):
                yield [row, column, *values]


HAIRPIN_HEADER = re.compile(r'PRIMER ([A-D])\s+(?:Hairpin Energy = (-?[\d.]+) kcals)?')


def drawing_line(line):
    """A drawing line without the border column: the hairpin page draws
    one column further right than the Primer Analysis screen, so this
    gives the columns of inspect/*.tsv."""
    return line[1:].rstrip(' │')


def hairpins(run_dir):
    screens = (Path(run_dir) / 'screens.txt').read_text().split('==== ')
    found = {}
    for screen in screens:
        if 'PRIMER HAIRPINS' not in screen:
            continue
        lines = screen.splitlines()
        for index, line in enumerate(lines):
            match = HAIRPIN_HEADER.search(line)
            if not match or match.group(1) in found:
                continue
            if match.group(2) is None:
                if 'No Hairpins Found' in lines[index + 1]:
                    found[match.group(1)] = ['-', '-', '-', '-']
                continue
            # the drawing follows two blank lines
            drawing = [drawing_line(lines[index + offset]) for offset in (3, 4, 5)]
            found[match.group(1)] = [match.group(2), *drawing]
    for field in FIELDS:
        if field in found:
            yield [field, *found[field]]


INTERACTION = re.compile(r'Interactions between primers ([A-D]) and ([A-D])')


def duplexes(run_dir):
    screens = (Path(run_dir) / 'screens.txt').read_text().split('==== ')
    found = {}
    for screen in screens:
        lines = screen.splitlines()
        for index, line in enumerate(lines):
            match = INTERACTION.search(line)
            if not match:
                continue
            # the drawing follows one blank line
            drawing = [drawing_line(lines[index + offset]) for offset in (2, 3, 4)]
            # skip the screens being redrawn: a blinking cursor, or a
            # duplex not drawn yet (empty top or bottom line); keep the
            # first complete one (when the next cell is drawn, the header
            # changes last)
            complete = not any('?' in text for text in drawing) and drawing[0] and drawing[2]
            key = (match.group(1), match.group(2))
            if complete and key not in found:
                found[key] = drawing
    for key in sorted(found):
        yield [*key, *found[key]]


def main(kind, run_dirs):
    columns = {
        'balance': ['field', 'position', 'primer', 'length', 'gc_percent', 'tm_celsius', 'gc_run_3prime', 'runs',
                    'hairpins'],
        'dimers': ['row', 'column', 'dimer_3prime', 'dimer_any'],
        'hairpins': ['field', 'hairpin_kcal', 'hairpin_5prime', 'hairpin_pairs', 'hairpin_3prime'],
        'duplexes': ['row', 'column', 'duplex_5prime', 'duplex_pairs', 'duplex_3prime'],
    }[kind]
    print('\t'.join(['# selection', 'list', *columns]))
    for run_dir in run_dirs:
        selection, *listed = Path(run_dir).name.split('_')
        # '_' in a selection (an empty field) also splits the name
        while listed and not re.fullmatch('[ACGT]+', listed[0]):
            selection += '_' + listed.pop(0)
        # the name ends with the EXTRA keys, if any
        while listed and not re.fullmatch('[ACGT]+', listed[-1]):
            listed.pop()
        if kind == 'hairpins':
            rows = hairpins(run_dir)
        elif kind == 'duplexes':
            rows = duplexes(run_dir)
        else:
            lines = last_screen(Path(run_dir) / 'screens.txt')
            rows = balance(lines) if kind == 'balance' else dimers(lines)
        for row in rows:
            print('\t'.join([selection, ','.join(listed), *row]))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
