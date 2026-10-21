"""Summarise the picking screens of a probe_picking.sh recording.

Usage: python3 -I parse_picking.py COLOURS_TXT
prints one tab-separated line per screen of MODIFY PRIMER POSITION
(From Molecule, Modify -> Position) in colours.txt:
  frame, primer length, MEETS CRITERIA SELECTED, the value line
  ("Pos: ... Runs: ..."), the highlighted bases ("FIRST-LAST colour",
  molecule positions, from the background of the sequence rows), and
  the other cells whose background is not the screen's blue (row,
  columns, foreground/background; outside the value line)
Screens drawn halfway (with '?' cells in the sequence rows) are left
out.
"""
import re
import sys

SCREEN_BACKGROUND = 'blue'
FIRST_BASE_COLUMN = 11  # column of the first base of a sequence row
BLOCK = 10  # bases per block, blocks are separated by a space
VALUE_LINE_ROW = 22


def screens(path):
    for block in open(path, encoding='utf-8').read().split('==== ')[1:]:
        lines = block.split('\n')
        text = lines[1:26]
        colours = {}
        for line in lines[26:]:
            match = re.match(r'row (\d+): (.*)', line)
            if match:
                colours[int(match.group(1))] = [
                    (int(first), int(last), fore, back)
                    for first, last, fore, back in re.findall(r'(\d+)-(\d+) (\S+)/(\S+)', match.group(2))]
        yield lines[0], text, colours


def base_position(line_start, column):
    offset = column - FIRST_BASE_COLUMN
    block, within = divmod(offset, BLOCK + 1)
    if offset < 0 or within == BLOCK:
        return None
    return line_start + block * BLOCK + within


def highlighted(text, colours):
    """[(first, last, colour)] of the bases on another background"""
    found = []
    for row, line in enumerate(text):
        match = re.match(r'│\s+(\d+)\s{3}[ACGTN]', line)
        if not match:
            continue
        start = int(match.group(1))
        for first, last, _, back in colours.get(row, []):
            if back == SCREEN_BACKGROUND:
                continue
            positions = [base_position(start, column) for column in range(first, last + 1)]
            positions = [position for position in positions if position is not None]
            if not positions:
                continue
            if found and found[-1][2] == back and found[-1][1] + 1 == positions[0]:
                found[-1] = (found[-1][0], positions[-1], back)
            else:
                found.append((positions[0], positions[-1], back))
    return found


def other_cells(text, colours):
    sequence_rows = {row for row, line in enumerate(text) if re.match(r'│\s+\d+\s{3}[ACGTN]', line)}
    cells = []
    for row, runs in sorted(colours.items()):
        if row in sequence_rows or row == VALUE_LINE_ROW:
            continue
        for first, last, fore, back in runs:
            if back not in (SCREEN_BACKGROUND, '?'):
                cells.append(f'{row}:{first}-{last} {fore}/{back} "{text[row][first:last + 1]}"')
    return cells


def main(path):
    for frame, text, colours in screens(path):
        if not any('PRIMER POSITION' in line for line in text):
            continue
        sequence_lines = [line for line in text if re.match(r'│\s+\d+\s{3}', line)]
        if any('?' in line for line in sequence_lines):
            continue
        length = re.search(r'Primer length: (\d+)', '\n'.join(text))
        meets = re.search(r'MEETS CRITERIA SELECTED: (\S+)', '\n'.join(text))
        values = next((line.strip('│ ') for line in text if 'Pos:' in line), '')
        marks = ' '.join(f'{first}-{last} {colour}' for first, last, colour in highlighted(text, colours))
        print('\t'.join([frame, length.group(1) if length else '', meets.group(1) if meets else '',
                         re.sub(r'\s+', ' ', values), marks, '; '.join(other_cells(text, colours))]))


if __name__ == '__main__':
    main(sys.argv[1])
