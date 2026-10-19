"""Write a FASTA file whose sequence text has an exact size in bytes.

Usage: python3 -I make_sized_fasta.py NAME SECTION_BYTES > NAME.FA

Random bases (seed 3), DOS line ends, lines of 60 bases; the last line
is longer or shorter so that the text after the header line has
exactly SECTION_BYTES bytes (tests/legacy/README.md, "Molecule size").
"""
import random
import sys

LINE_BASES = 60
LINE_END = '\r\n'
LINE_BYTES = LINE_BASES + len(LINE_END)


def lines_for(section_bytes):
    full, rest = divmod(section_bytes, LINE_BYTES)
    # a last line needs at least one base: borrow a full line if needed
    if 0 < rest <= len(LINE_END):
        full -= 1
        rest += LINE_BYTES
    return full, rest


def main(name, section_bytes):
    rnd = random.Random(3)
    full, rest = lines_for(section_bytes)
    lines = [''.join(rnd.choice('ACGT') for _ in range(LINE_BASES)) for _ in range(full)]
    if rest:
        lines.append(''.join(rnd.choice('ACGT') for _ in range(rest - len(LINE_END))))
    body = ''.join(line + LINE_END for line in lines)
    assert len(body) == section_bytes
    sys.stdout.buffer.write(f'>{name}{LINE_END}{body}'.encode())


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
