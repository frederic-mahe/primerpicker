"""Minimal ZMBV AVI decoder: writes PNGs of distinct frames."""
import struct, sys, zlib, hashlib
from pathlib import Path
from PIL import Image

def chunks(data):
    pos = 12
    stack = [(12, len(data))]
    out = []
    def walk(start, end):
        p = start; end = min(end, len(data))
        while p + 8 <= end:
            cid, size = data[p:p+4], struct.unpack('<I', data[p+4:p+8])[0]
            if cid in (b'RIFF', b'LIST'):
                walk(p + 12, p + 8 + size if size else len(data))
            else:
                out.append((cid, data[p+8:p+8+size]))
            p += 8 + size + (size & 1)
    walk(12, len(data))
    return out

def decode(path, outdir, every=1):
    data = Path(path).read_bytes()
    cs = chunks(data)
    avih = next(c for i, c in cs if i == b'avih')
    width, height = struct.unpack('<II', avih[32:40])
    frames = [c for i, c in cs if i == b'00dc']
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    pal = bytearray(768); cur = bytearray(width * height); dec = None
    bw = bh = 16; last = None; saved = 0
    for n, fr in enumerate(frames):
        if not fr: continue
        flags = fr[0]; p = 1
        if flags & 1:
            hi, lo, comp, fmt, bw, bh = fr[1:7]; p = 7
            bpp = {4: 1, 8: 4}[fmt]
            if len(cur) != width * height * bpp: cur = bytearray(width * height * bpp)
            dec = zlib.decompressobj()
            raw = dec.decompress(fr[p:])
            off = 768 if bpp == 1 else 0
            if bpp == 1: pal[:] = raw[:768]
            cur[:] = raw[off:off + width * height * bpp]
        else:
            raw = dec.decompress(fr[p:]); q = 0
            if flags & 2 and bpp == 1:
                for k in range(768): pal[k] ^= raw[k]
                q = 768
            bx = (width + bw - 1) // bw; by = (height + bh - 1) // bh
            nvec = bx * by
            vecs = raw[q:q + 2 * nvec]; q += (2 * nvec + 3) & ~3
            prev = bytes(cur)
            for b in range(nvec):
                dx = struct.unpack('b', vecs[2*b:2*b+1])[0]; dy = struct.unpack('b', vecs[2*b+1:2*b+2])[0]
                xor = dx & 1; dx >>= 1; dy >>= 1
                x0 = (b % bx) * bw; y0 = (b // bx) * bh
                w = min(bw, width - x0); h = min(bh, height - y0)
                if not xor and dx == 0 and dy == 0: continue
                for yy in range(h):
                    sy = y0 + yy + dy; row = (y0 + yy) * width
                    for xx in range(w):
                        sx = x0 + xx + dx
                        for c in range(bpp):
                            v = prev[(sy * width + sx) * bpp + c] if 0 <= sx < width and 0 <= sy < height else 0
                            if xor: v ^= raw[q]; q += 1
                            cur[(row + x0 + xx) * bpp + c] = v
        h = hashlib.md5(bytes(cur) + bytes(pal)).digest()
        if h != last and n % every == 0:
            if bpp == 1:
                img = Image.frombytes('P', (width, height), bytes(cur)); img.putpalette(bytes(pal))
            else:
                img = Image.frombytes('RGB', (width, height), bytes(cur), 'raw', 'BGRX')
            img.save(outdir / f'f{n:05d}.png'); saved += 1; last = h
    print(f'{len(frames)} frames, {width}x{height}, saved {saved} distinct')

if __name__ == '__main__':
    decode(sys.argv[1], sys.argv[2])
