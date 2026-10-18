"""Print the restriction enzyme table stored in PRIMER.exe.

Usage: python3 -I extract_enzymes.py PRIMER.EXE > enzymes/enzyme_table.tsv

PRIMER.exe is packed with Microsoft EXEPACK: the load image is
unpacked first (commands read backwards from the end of the packed
data: B0 = fill, B2 = copy, the low bit marks the last command). The
table then holds one 28-byte record per enzyme, in alphabetical
order:

    name                  8 bytes, NUL padded
    recognition sequence 14 bytes, IUPAC letters, NUL padded
    cut                   1 byte, signed (bases after the site start)
    (flags)               1 byte (2 when the site has non-ACGT letters;
                          1 for some enzymes that cut outside it)
    overhang              1 byte, signed (> 0: 5' overhang)
    length                1 byte (of the recognition sequence)
    flags                 1 byte (0x80: in the default enzyme menu)
    (unknown)             1 byte, always 0

The meaning of cut, overhang and both flag bytes is inferred from
known enzymes and from the default menu, not from documentation.
Prints one line per enzyme: name, recognition sequence, cut, overhang,
and whether the enzyme is in the default menu.
"""
import re
import struct
import sys

PARAGRAPH = 16
EXEPACK_SIGNATURE = b'RB'
RECORD_SIZE = 28
NAME_SIZE = 8
SEQUENCE_SIZE = 14
FIRST_RECORD = re.compile(rb'AatI\x00{4}AGGCCT\x00')
NAME = re.compile(rb'[A-Z][A-Za-z0-9-]{1,7}\x00*')
SEQUENCE = re.compile(rb'[ACGTMKRYSWBDHVN]{3,14}\x00*')
IN_DEFAULT_MENU = 0x80


def unpack_exepack(data):
    """The load image of an EXEPACK-packed MZ executable"""
    header_size = struct.unpack_from('<H', data, 8)[0] * PARAGRAPH
    initial_cs = struct.unpack_from('<H', data, 22)[0]
    exepack = header_size + initial_cs * PARAGRAPH
    fields = struct.unpack_from('<7H2s', data, exepack)
    destination_paragraphs, signature = fields[6], fields[7]
    if signature != EXEPACK_SIGNATURE:
        raise SystemExit('not an EXEPACK executable')
    packed = data[header_size:exepack]
    image = bytearray(destination_paragraphs * PARAGRAPH)
    image[:len(packed)] = packed
    source = len(packed)
    while packed[source - 1] == 0xff:  # padding
        source -= 1
    destination = len(image)
    while True:
        command = packed[source - 1]
        length = struct.unpack_from('<H', packed, source - 3)[0]
        source -= 3
        if command & 0xfe == 0xb0:  # fill
            fill = packed[source - 1]
            source -= 1
            image[destination - length:destination] = bytes([fill]) * length
        elif command & 0xfe == 0xb2:  # copy
            image[destination - length:destination] = packed[source - length:source]
            source -= length
        else:
            raise SystemExit(f'unknown EXEPACK command {command:#x}')
        destination -= length
        if command & 1:
            return bytes(image)


def records(image):
    start = FIRST_RECORD.search(image)
    if start is None:
        raise SystemExit('enzyme table not found')
    offset = start.start()
    while True:
        record = image[offset:offset + RECORD_SIZE]
        name, sequence = record[:NAME_SIZE], record[NAME_SIZE:NAME_SIZE + SEQUENCE_SIZE]
        if NAME.fullmatch(name) is None or SEQUENCE.fullmatch(sequence) is None:
            return
        cut, _, overhang, length, flags, _ = struct.unpack_from('<bBbBBB', record, NAME_SIZE + SEQUENCE_SIZE)
        sequence = sequence.rstrip(b'\x00').decode('ascii')
        if length != len(sequence):
            raise SystemExit(f'length field {length} for {sequence}')
        yield name.rstrip(b'\x00').decode('ascii'), sequence, cut, overhang, flags & IN_DEFAULT_MENU != 0
        offset += RECORD_SIZE


def main(executable):
    image = unpack_exepack(open(executable, 'rb').read())
    table = list(records(image))
    print('# name\trecognition_sequence\tcut\toverhang\tdefault_menu')
    for name, sequence, cut, overhang, in_menu in table:
        print(name, sequence, cut, overhang, 'Y' if in_menu else 'N', sep='\t')
    print(f'{len(table)} enzymes, {sum(entry[4] for entry in table)} in the default menu', file=sys.stderr)


if __name__ == '__main__':
    main(*sys.argv[1:])
