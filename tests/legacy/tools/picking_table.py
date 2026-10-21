"""Tabulate the settled states of picking screens (probe_picking.sh runs).

Usage: python3 -I picking_table.py MOLECULE RUN_DIR... > picking.tsv
MOLECULE: the molecule file of the runs, as written in the table (a
path from tests/legacy, such as pairs/TEST600.GB). For each run (a directory with colours.txt), the screens of MODIFY
PRIMER POSITION from parse_picking.py, keeping a screen only when it
is complete: the primer length, MEETS CRITERIA SELECTED and the value
line are drawn, and the highlighted window is one stretch of bases that
starts at the shown position (top strand) or ends there (complement
strand, "*** Complementary Strand ***"), as long as the primer. A state
equal to the one before it is printed once. Columns: run, molecule, frame,
strand, primer length, position, MEETS CRITERIA SELECTED, the value
line (spaces kept, as shown from column 2 of the screen), the
window's background colour.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python3 -I leaves it out
import parse_picking  # noqa: E402

VALUE_LINE_COLUMN = 2  # after the frame and one space
WINDOW_COLOURS = ('lightgrey', 'green')  # fails, meets the criteria
SETTLED_FRAMES = 30  # about 0.4 s of capture; keys come every 2.5 s or more


def frame_number(frame):
    return int(re.search(r'\d+', frame).group(0))


def state(text, colours):
    """(strand, length, position, meets, value line, colour) of a
    complete picking screen, or None"""
    if not any('PRIMER POSITION' in line for line in text):
        return None
    whole = '\n'.join(text)
    if '?' in ''.join(line for line in text if re.match(r'│\s+\d+\s{3}', line)):
        return None
    length = re.search(r'Primer length: (\d+)', whole)
    meets = re.search(r'MEETS CRITERIA SELECTED: (Yes|No)', whole)
    value_line = next((line for line in text if 'Pos:' in line), None)
    if not (length and meets and value_line) or '?' in value_line.replace('Tm°C: ?', ''):
        return None
    position = int(re.search(r'Pos: (\d+)', value_line).group(1))
    strand = 'complement' if 'Complementary Strand' in whole else 'top'
    windows = []
    for first, last, colour in parse_picking.highlighted(text, colours):
        if colour not in WINDOW_COLOURS:
            continue
        # one cell of the screen (row 12, column 40) keeps a colour of
        # its own (brown on the blue screen, magenta on green, black on
        # light grey): a window split by one base is one window
        if windows and windows[-1][2] == colour and windows[-1][1] + 2 == first:
            windows[-1] = (windows[-1][0], last, colour)
            continue
        windows.append((first, last, colour))
    if len(windows) != 1:
        return None
    first, last, colour = windows[0]
    size = int(length.group(1))
    expected_first = position if strand == 'top' else position - size + 1
    if first != expected_first or last - first + 1 != size:
        return None
    shown = value_line[VALUE_LINE_COLUMN:].rstrip(' │').rstrip()
    return strand, size, position, meets.group(1), shown, colour


def states(run_dir):
    """(first frame, state) of the settled states of a run"""
    groups = []  # [first frame, state], consecutive screens of one state
    for frame, text, colours in parse_picking.screens(os.path.join(run_dir, 'colours.txt')):
        current = state(text, colours)
        if groups and groups[-1][1] == current:
            continue
        groups.append([frame, current])
    for index, (frame, current) in enumerate(groups):
        if current is None:
            continue
        last = index + 1 == len(groups)
        if last or frame_number(groups[index + 1][0]) - frame_number(frame) > SETTLED_FRAMES:
            yield frame, current


def main(molecule, run_dirs):
    print('# run\tmolecule\tframe\tstrand\tlength\tposition\tmeets\tvalue_line\twindow_colour')
    for run_dir in run_dirs:
        previous = None
        name = os.path.basename(os.path.normpath(run_dir))
        for frame, current in states(run_dir):
            if current == previous:
                continue
            previous = current
            print('\t'.join(str(field) for field in [name, molecule, frame, *current]))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
