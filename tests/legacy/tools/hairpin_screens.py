"""Keep only the Hairpins screen of inspect_primer.sh recordings.

Usage: python3 -I hairpin_screens.py RUN_DIR... > hairpins_screens.txt
(each RUN_DIR is named after its primer and holds screens.txt, with
the splash screens already removed by strip_splash.py)

Prints, for each primer, a line "#### PRIMER" and the first complete
PRIMER HAIRPINS screen (energy and drawing, or "No Hairpins Found")
of its recording, as screen_text.py wrote it.
"""
import sys
from pathlib import Path


def hairpin_screen(text):
    for block in text.split('==== ')[1:]:
        if 'PRIMER HAIRPINS' not in block:
            continue
        if 'No Hairpins Found' in block:
            return '==== ' + block.rstrip('\n')
        lines = block.split('\n')
        top = next((k for k, line in enumerate(lines) if "5'-" in line), None)
        if 'Hairpin Energy' not in block or top is None:
            continue
        if any('?' in line for line in lines[top:top + 3]):  # the cursor hides a cell
            continue
        return '==== ' + block.rstrip('\n')
    raise ValueError('no complete Hairpins screen')


def main(run_dirs):
    for run_dir in run_dirs:
        path = Path(run_dir)
        text = (path / 'screens.txt').read_text(encoding='utf-8')
        try:
            screen = hairpin_screen(text)
        except ValueError as error:
            print(f'{path.name}: {error}', file=sys.stderr)
            continue
        print(f'#### {path.name}')
        print(screen)


if __name__ == '__main__':
    main(sys.argv[1:])
