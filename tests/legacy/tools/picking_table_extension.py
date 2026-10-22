"""Tabulate the settled states of picking screens of primers with a 5'
extension (probe_picking.sh runs).

Usage: python3 -I picking_table_extension.py MOLECULE RUN_DIR... > TABLE.tsv

As picking_table.py, for primers with a 5' extension: the header line
reads "5'extension: SEQUENCE", "Primer length" counts the extension, and
the highlighted window covers the bases of the molecule only (length
minus extension), from the shown position (top strand). Columns: those
of picking_table.py, with the extension after the length.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python3 -I leaves it out

import parse_picking  # noqa: E402
import picking_table  # noqa: E402


def state(text, colours):
    """(strand, length, extension, position, meets, value line, colour)
    of a complete picking screen of a top-strand primer with a 5'
    extension, or None"""
    whole = '\n'.join(text)
    extension = re.search(r"5'extension: ([ACGT]+)", whole)
    if not extension or 'Complementary Strand' in whole:
        return None
    if not any('PRIMER POSITION' in line for line in text):
        return None
    if '?' in ''.join(line for line in text if re.match(r'│\s+\d+\s{3}', line)):
        return None
    length = re.search(r'Primer length: (\d+)', whole)
    meets = re.search(r'MEETS CRITERIA SELECTED: (Yes|No)', whole)
    value_line = next((line for line in text if 'Pos:' in line), None)
    if not (length and meets and value_line) or '?' in value_line.replace('Tm°C: ?', ''):
        return None
    position = int(re.search(r'Pos: (\d+)', value_line).group(1))
    windows = []
    for first, last, colour in parse_picking.highlighted(text, colours):
        if colour not in picking_table.WINDOW_COLOURS:
            continue
        # the cell of its own colour (picking_table.py): one window
        if windows and windows[-1][2] == colour and windows[-1][1] + 2 == first:
            windows[-1] = (windows[-1][0], last, colour)
            continue
        windows.append((first, last, colour))
    if len(windows) != 1:
        return None
    first, last, colour = windows[0]
    size = int(length.group(1))
    if first != position or last - first + 1 != size - len(extension.group(1)):
        return None
    shown = value_line[picking_table.VALUE_LINE_COLUMN:].rstrip(' │').rstrip()
    return 'top', size, extension.group(1), position, meets.group(1), shown, colour


def main(molecule, run_dirs):
    print('# run\tmolecule\tframe\tstrand\tlength\textension\tposition\tmeets\tvalue_line\twindow_colour')
    picking_table.state = state  # the settled states of picking_table.py, with this state
    for run_dir in run_dirs:
        previous = None
        name = os.path.basename(os.path.normpath(run_dir))
        for frame, current in picking_table.states(run_dir):
            if current == previous:
                continue
            previous = current
            print('\t'.join(str(field) for field in [name, molecule, frame, *current]))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
