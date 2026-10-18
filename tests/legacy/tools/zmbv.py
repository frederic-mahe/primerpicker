"""Decode a ZMBV video (DOSBox-X's dx-capture /V) into PNG images.

Usage: python3 -I zmbv.py CAPTURE.avi OUTPUT_DIR

Writes one PNG per distinct frame (a frame identical to the previous
one is skipped), named fNNNNN.png after its frame number. Only what
DOSBox-X writes is supported: an AVI file whose video chunks are
ZMBV-compressed, 8-bit palettised or 32-bit colour frames.
"""
import hashlib
import struct
import sys
import zlib
from pathlib import Path

from PIL import Image

PALETTE_SIZE = 768  # 256 colours x 3 bytes
KEYFRAME = 0x01
DELTA_PALETTE = 0x02
BYTES_PER_PIXEL = {4: 1, 8: 4}  # ZMBV format code -> bytes per pixel
CONTAINERS = (b'RIFF', b'LIST')
CHUNK_HEADER_SIZE = 8
LIST_HEADER_SIZE = 12


def chunks(data):
    """The (id, payload) of every leaf chunk of an AVI (RIFF) file."""
    found = []

    def walk(start, end):
        position = start
        end = min(end, len(data))
        while position + CHUNK_HEADER_SIZE <= end:
            chunk_id = data[position:position + 4]
            size = struct.unpack('<I', data[position + 4:position + 8])[0]
            if chunk_id in CONTAINERS:
                list_end = position + CHUNK_HEADER_SIZE + size if size else len(data)
                walk(position + LIST_HEADER_SIZE, list_end)
            else:
                payload_start = position + CHUNK_HEADER_SIZE
                found.append((chunk_id, data[payload_start:payload_start + size]))
            position += CHUNK_HEADER_SIZE + size + (size & 1)  # chunks are word-aligned

    walk(LIST_HEADER_SIZE, len(data))
    return found


class Decoder:
    """Keeps the current frame, palette and zlib stream between frames."""

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.palette = bytearray(PALETTE_SIZE)
        self.bytes_per_pixel = 1
        self.frame = bytearray(width * height)
        self.block_width = self.block_height = 16
        self.stream = None

    def keyframe(self, payload):
        # header: version (2 bytes), compression, format, block width and height
        format_code, self.block_width, self.block_height = payload[3], payload[4], payload[5]
        self.bytes_per_pixel = BYTES_PER_PIXEL[format_code]
        self.stream = zlib.decompressobj()
        raw = self.stream.decompress(payload[6:])
        offset = 0
        if self.bytes_per_pixel == 1:
            self.palette[:] = raw[:PALETTE_SIZE]
            offset = PALETTE_SIZE
        size = self.width * self.height * self.bytes_per_pixel
        self.frame = bytearray(raw[offset:offset + size])

    def delta(self, flags, payload):
        raw = self.stream.decompress(payload)
        offset = 0
        if flags & DELTA_PALETTE and self.bytes_per_pixel == 1:
            for index in range(PALETTE_SIZE):
                self.palette[index] ^= raw[index]
            offset = PALETTE_SIZE
        blocks_across = (self.width + self.block_width - 1) // self.block_width
        blocks_down = (self.height + self.block_height - 1) // self.block_height
        block_count = blocks_across * blocks_down
        vectors = raw[offset:offset + 2 * block_count]
        offset += (2 * block_count + 3) & ~3  # vectors are padded to 4 bytes
        previous = bytes(self.frame)
        for block in range(block_count):
            offset = self.apply_block(raw, offset, previous, vectors, block, blocks_across)

    def apply_block(self, raw, offset, previous, vectors, block, blocks_across):
        """Copy one block from the previous frame (moved, maybe XORed)."""
        dx, dy = struct.unpack('bb', vectors[2 * block:2 * block + 2])
        xored = dx & 1
        dx >>= 1
        dy >>= 1
        if not xored and dx == 0 and dy == 0:
            return offset  # unchanged block
        left = (block % blocks_across) * self.block_width
        top = (block // blocks_across) * self.block_height
        block_width = min(self.block_width, self.width - left)
        block_height = min(self.block_height, self.height - top)
        size = self.bytes_per_pixel
        for row in range(top, top + block_height):
            source_row = row + dy
            for column in range(left, left + block_width):
                source_column = column + dx
                inside = 0 <= source_column < self.width and 0 <= source_row < self.height
                source = (source_row * self.width + source_column) * size
                target = (row * self.width + column) * size
                for channel in range(size):
                    value = previous[source + channel] if inside else 0
                    if xored:
                        value ^= raw[offset]
                        offset += 1
                    self.frame[target + channel] = value
        return offset

    def image(self):
        if self.bytes_per_pixel == 1:
            picture = Image.frombytes('P', (self.width, self.height), bytes(self.frame))
            picture.putpalette(bytes(self.palette))
            return picture
        return Image.frombytes('RGB', (self.width, self.height), bytes(self.frame), 'raw', 'BGRX')

    def fingerprint(self):
        return hashlib.md5(bytes(self.frame) + bytes(self.palette)).digest()


def decode(avi_path, output_dir):
    data = Path(avi_path).read_bytes()
    found = chunks(data)
    main_header = next(payload for chunk_id, payload in found if chunk_id == b'avih')
    width, height = struct.unpack('<II', main_header[32:40])
    frames = [payload for chunk_id, payload in found if chunk_id == b'00dc']
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    decoder = Decoder(width, height)
    last = None
    saved = 0
    for number, frame in enumerate(frames):
        if not frame:
            continue  # a dropped frame repeats the previous one
        flags = frame[0]
        if flags & KEYFRAME:
            decoder.keyframe(frame[1:])
        else:
            decoder.delta(flags, frame[1:])
        fingerprint = decoder.fingerprint()
        if fingerprint != last:
            decoder.image().save(output / f'f{number:05d}.png')
            saved += 1
            last = fingerprint
    print(f'{len(frames)} frames, {width}x{height}, saved {saved} distinct')


if __name__ == '__main__':
    decode(sys.argv[1], sys.argv[2])
