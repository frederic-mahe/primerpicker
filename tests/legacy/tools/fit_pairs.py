"""Check the product Tm and annealing temperature models of PRIMER.exe.

Models (see ../README.md, "Product Tm and annealing temperature"):
    product Tm = floor(59.9 + 41 x GC fraction - 675 / length)
    Ta = floor(0.3 x (lower primer Tm - 21.6) + 0.7 x product Tm - 14.9)
with the displayed (integer) primer and product Tm values.

Also checks that the recorded primers and product size agree with the
molecule, and that the primer Tm values follow fit_tm.py's model.

Usage: python3 -I fit_pairs.py ../pairs/pairs.tsv
(molecule files are read from the directory of the table)
"""
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from fit_tm import shown_tm  # noqa: E402

ENTROPY_CONSTANT = 60.699  # fit_tm.py


def molecule_sequence(path):
    text = path.read_text(encoding='ascii')
    origin = text.split('ORIGIN', 1)[1].split('//', 1)[0]
    return re.sub('[^acgt]', '', origin).upper()


def reverse_complement(sequence):
    return sequence[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def product_tm(product):
    gc_fraction = sum(base in 'GC' for base in product) / len(product)
    return math.floor(59.9 + 41 * gc_fraction - 675 / len(product))


def annealing_tm(primer_a_tm, primer_b_tm, product_tm_shown):
    lower = min(primer_a_tm, primer_b_tm)
    return math.floor(0.3 * (lower - 21.6) + 0.7 * product_tm_shown - 14.9)


def check(fields, directory):
    (molecule, position_a, position_b, primer_a, primer_b,
     tm_a, tm_b, size, shown_product_tm, shown_annealing_tm) = fields
    sequence = molecule_sequence(directory / molecule)
    product = sequence[int(position_a) - 1:int(position_b)]
    tm_a, tm_b = int(tm_a.rstrip('!')), int(tm_b.rstrip('!'))
    problems = []
    if len(product) != int(size):
        problems.append('size')
    if not product.startswith(primer_a) or not product.endswith(reverse_complement(primer_b)):
        problems.append('primer positions')
    if shown_tm(primer_a, ENTROPY_CONSTANT) != str(tm_a) or shown_tm(primer_b, ENTROPY_CONSTANT) != str(tm_b):
        problems.append('primer Tm')
    if product_tm(product) != int(shown_product_tm):
        problems.append('product Tm')
    if annealing_tm(tm_a, tm_b, int(shown_product_tm)) != int(shown_annealing_tm):
        problems.append('annealing temperature')
    return problems


def main(path):
    table = Path(path)
    records = 0
    failures = 0
    with table.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            records += 1
            problems = check(fields, table.parent)
            if problems:
                failures += 1
                print('mismatch:', ' '.join(fields[:3]), ', '.join(problems))
    print(f'{records} records, {failures} mismatches')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
