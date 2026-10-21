"""Remove the splash screens (licence serial number) from a screens.txt.

screen_text.py prints every distinct screen of a recording, from the
DOS prompt on, including the splash screen of PRIMER.exe, which shows
the licence serial number: it must not be published. This keeps the
screens that follow the last splash screen (the first screen of the
program proper, usually its main menu, and everything after it).

Usage: python3 -I strip_splash.py SCREENS_TXT > PUBLISHED_SCREENS_TXT
"""
import sys

SPLASH_MARKERS = ('Serial Number', 'D E S I G N E R')


def screens(text):
    """Split screen_text.py output into (header, body) pairs."""
    current = None
    for line in text.splitlines():
        if line.startswith('==== '):
            if current is not None:
                yield current
            current = (line, [])
        elif current is not None:
            current[1].append(line)
    if current is not None:
        yield current


def main(path):
    with open(path, encoding='utf-8', errors='surrogateescape') as file:
        parsed = list(screens(file.read()))
    last_splash = -1
    for index, (_, body) in enumerate(parsed):
        if any(marker in line for line in body for marker in SPLASH_MARKERS):
            last_splash = index
    out = sys.stdout.buffer
    for header, body in parsed[last_splash + 1:]:
        for line in [header, *body]:
            out.write(line.encode('utf-8', errors='surrogateescape') + b'\n')


if __name__ == '__main__':
    main(sys.argv[1])
