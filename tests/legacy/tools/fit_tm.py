"""Check the primer Tm model of PRIMER.exe against recorded values.

Model (see ../README.md, "Primer Tm"):
    Tm = floor(1000 * dH / (dS + b) - 273)    shown as '?' if Tm <= 0
with dH (kcal/mol) and dS (cal/K/mol) summed over the nearest-neighbour
pairs of the primer (5' extension excluded), Breslauer et al. 1986
parameters, and a constant b (cal/K/mol) fitted here.

Usage: python3 -I fit_tm.py ../tm/direct_entry_tm.tsv
"""
import math
import sys

# Breslauer KJ, Frank R, Blöcker H, Marky LA (1986) Predicting DNA
# duplex stability from the base sequence. PNAS 83:3746-3750, Table 2.
# Pair: (dH kcal/mol, dS cal/K/mol), as positive dissociation values.
# Checked against Primer3 2.6.1 (oligotm -tp 0), whose sums agree for
# all 16 pairs (Primer3 adds a constant initiation entropy).
BRESLAUER_1986 = {
    'AA': (9.1, 24.0), 'TT': (9.1, 24.0),
    'AT': (8.6, 23.9),
    'TA': (6.0, 16.9),
    'CA': (5.8, 12.9), 'TG': (5.8, 12.9),
    'GT': (6.5, 17.3), 'AC': (6.5, 17.3),
    'CT': (7.8, 20.8), 'AG': (7.8, 20.8),
    'GA': (5.6, 13.5), 'TC': (5.6, 13.5),
    'CG': (11.9, 27.8),
    'GC': (11.1, 26.7),
    'GG': (11.0, 26.6), 'CC': (11.0, 26.6),
}
MAX_PRIMER_LENGTH = 50  # longer entries are cut by the entry field
KELVIN_OFFSET = 273  # not 273.15: see README


def annealing_part(typed):
    """The primer as the program keeps it: A, C, G, T only, 50 at most.

    Other letters are rejected, and the key typed after them is lost
    (see README).
    """
    kept = []
    skip_next = False
    for base in typed:
        if skip_next:
            skip_next = False
        elif base in 'ACGT':
            kept.append(base)
        else:
            skip_next = True
    return ''.join(kept)[:MAX_PRIMER_LENGTH]


def sums(primer):
    pairs = [primer[i:i + 2] for i in range(len(primer) - 1)]
    return (sum(BRESLAUER_1986[pair][0] for pair in pairs),
            sum(BRESLAUER_1986[pair][1] for pair in pairs))


def shown_tm(primer, entropy_constant):
    enthalpy, entropy = sums(primer)
    tm = math.floor(1000 * enthalpy / (entropy + entropy_constant) - KELVIN_OFFSET)
    return str(tm) if tm > 0 else '?'


def read_records(path):
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            yield annealing_part(fields[1]), fields[5]


def main(path):
    records = list(read_records(path))
    candidates = [step / 10000 for step in range(600000, 620000)]
    fitting = [b for b in candidates
               if all(shown_tm(primer, b) == tm for primer, tm in records)]
    print(f'{len(records)} records')
    if not fitting:
        print('no value of b fits all records')
        return 1
    print(f'all records match for b in [{fitting[0]:.4f}, {fitting[-1]:.4f}] cal/K/mol')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
