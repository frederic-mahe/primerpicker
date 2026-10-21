"""Check recorded hairpins against the legacy hairpin model, and find
which G.T step values (or sums of them) the recordings determine.

Usage:
  python3 -I fit_gt_steps.py check [--old-zip] [--no-tandem] TSV...
  python3 -I fit_gt_steps.py fit TSV...          (needs scipy)

TSV: tables in the format of parse_inspect.py (primer in column 1, the
shown energy in column 7, the marks line in column 9).

The model is a Python copy of primerpicker's legacy_hairpin.hpp
(tests/legacy/README.md, "Primer criteria"), with the changes proposed
in PROBES_20261021_hairpin.md: the tandem step 5'-TG-3'/3'-GT-5'
(T.G then G.T) at 0.60 kcal, and stems zipped inwards through G.T
pairs. --old-zip zips over Watson-Crick pairs only (the rule of the
README before these probes); --no-tandem leaves the tandem step out
(stems stop there, as primerpicker did).

check: prints one line per recording whose shown energy or number of
pairs (':' marks) differs from the model, then a count.

fit: treats the eleven steps that hold a G.T pair as unknowns
(hundredths of kcal), turns each recording into linear inequalities
(the shown energy is truncated toward zero to tenths; the shown hairpin
is the most stable one), and prints the feasible range of the sums the
recordings can determine: the steps on both sides of a single G.T pair
(L(X) + R(Y)), the tandem steps with their neighbours, and GG/TT.
L(X) is the step from the Watson-Crick pair X to a G.T pair (G on the
5' strand), R(Y) the step from that G.T pair to Y.
"""
import sys

WATSON_CRICK = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
WOBBLE_PAIRS = {('G', 'T'), ('T', 'G')}
# Freier et al. 1986, Watson-Crick stacks (hundredths of kcal/mol)
STACK = {'AA': 90, 'AC': 210, 'AG': 170, 'AT': 90, 'CA': 180, 'CC': 290, 'CG': 200, 'CT': 170,
         'GA': 230, 'GC': 340, 'GG': 290, 'GT': 210, 'TA': 110, 'TC': 230, 'TG': 180, 'TT': 90}
# steps with a G.T pair, 5'-ab-3'/3'-cd-5' as (a, b, c, d): freier_1986.hpp,
# plus the proposed tandem step TG/GT
WOBBLE = {('A', 'G', 'T', 'T'): 50, ('A', 'T', 'T', 'G'): 70, ('C', 'G', 'G', 'T'): 150,
          ('C', 'T', 'G', 'G'): 150, ('G', 'A', 'T', 'T'): 50, ('G', 'C', 'T', 'G'): 190,
          ('G', 'G', 'C', 'T'): 130, ('G', 'G', 'T', 'T'): 50, ('G', 'T', 'T', 'G'): 50,
          ('T', 'A', 'G', 'T'): 70, ('T', 'G', 'G', 'T'): 60}
TANDEM = ('T', 'G', 'G', 'T')
LOOP = [(3, 740), (4, 590), (5, 440), (6, 430), (7, 410), (8, 410), (9, 420), (10, 430),
        (12, 490), (14, 560), (16, 610), (18, 670), (20, 710), (30, 890)]
MIN_LOOP = 3
PARAMS = list(WOBBLE)
NAMES = [f'{a}{b}/{c}{d}' for a, b, c, d in PARAMS]


def loop_penalty(length):
    if length < MIN_LOOP:
        return None
    for (low, low_penalty), (high, high_penalty) in zip(LOOP, LOOP[1:]):
        if length <= high:
            return low_penalty + (high_penalty - low_penalty) * (length - low) // (high - low)
    (low, low_penalty), (high, high_penalty) = LOOP[-2], LOOP[-1]
    return low_penalty + (high_penalty - low_penalty) * (length - low) // (high - low)


def pairs(first, second):
    return (first, second) in WATSON_CRICK


def wobble(first, second):
    return (first, second) in WOBBLE_PAIRS


def canonical(step, table):
    """the key of a step in the table (a step equals its reverse), or None"""
    a, b, c, d = step
    for key in (step, (d, c, b, a)):
        if key in table:
            return key
    return None


class Model:
    def __init__(self, zip_through_wobble=True, tandem=True, values=None):
        self.zip_through_wobble = zip_through_wobble
        self.table = dict(values or WOBBLE)
        if not tandem:
            self.table.pop(TANDEM, None)

    def step(self, primer, outer_5, outer_3):
        """the step between (outer_5, outer_3) and the next pair inwards:
        ('wc', value), ('wobble', key), or None if it cannot be in a stem"""
        a, b, c, d = primer[outer_5], primer[outer_5 + 1], primer[outer_3], primer[outer_3 - 1]
        if pairs(a, c) and pairs(b, d):
            return ('wc', STACK[a + b])
        if not ((pairs(a, c) or wobble(a, c)) and (pairs(b, d) or wobble(b, d))):
            return None
        key = canonical((a, b, c, d), self.table)
        return None if key is None else ('wobble', key)

    def zips_inward(self, primer, inner_5, inner_3):
        depth = 1
        while True:
            top, bottom = inner_5 + depth, inner_3 - depth
            if bottom - top - 1 < MIN_LOOP:
                return False
            if pairs(primer[top], primer[bottom]):
                return True
            if not (self.zip_through_wobble and wobble(primer[top], primer[bottom])):
                return False
            depth += 1

    def candidates(self, primer):
        """every scored stem: (start, pairs, loop, constant, wobble steps)"""
        found = []
        for inner_5 in range(len(primer)):
            for inner_3 in range(inner_5 + MIN_LOOP + 1, len(primer)):
                if not pairs(primer[inner_5], primer[inner_3]) or self.zips_inward(primer, inner_5, inner_3):
                    continue
                loop = inner_3 - inner_5 - 1
                length = 1
                while (length <= inner_5 and inner_3 + length < len(primer)
                       and self.step(primer, inner_5 - length, inner_3 + length) is not None):
                    length += 1
                while length > 1 and not pairs(primer[inner_5 + 1 - length], primer[inner_3 + length - 1]):
                    length -= 1
                constant, steps = loop_penalty(loop), []
                for offset in range(1, length):
                    kind, value = self.step(primer, inner_5 - offset, inner_3 + offset)
                    if kind == 'wc':
                        constant -= value
                    else:
                        steps.append(value)
                found.append((inner_5 + 1 - length, length, loop, constant, steps))
        return found

    def energy(self, candidate):
        return candidate[3] - sum(self.table[key] for key in candidate[4])

    def hairpin(self, primer):
        best = min(self.candidates(primer), key=self.energy, default=None)
        if best is None or self.energy(best) >= 0:
            return None
        return best


def shown(energy):
    """truncated toward zero to tenths, written as PRIMER.exe does (-.4)"""
    if energy is None:
        return '-'
    tenths = -((-energy) // 10)
    text = f'{tenths / 10:.1f}'
    return text.replace('-0.', '-.') if -10 < tenths < 0 else text


def recordings(paths):
    for path in paths:
        with open(path, encoding='utf-8') as table:
            for line in table:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.rstrip('\n').split('\t')
                marks = fields[8].count(':') if fields[6] != '-' else 0
                yield path, fields[0], fields[6], marks


def check(arguments):
    model = Model(zip_through_wobble='--old-zip' not in arguments, tandem='--no-tandem' not in arguments)
    paths = [argument for argument in arguments if not argument.startswith('--')]
    count = misfits = 0
    for path, primer, energy, marks in recordings(paths):
        count += 1
        hairpin = model.hairpin(primer)
        predicted = shown(model.energy(hairpin) if hairpin else None)
        length = hairpin[1] if hairpin else 0
        if predicted != energy or length != marks:
            misfits += 1
            print(f'misfit\t{path}\t{primer}\tshown {energy} ({marks} pairs)\tmodel {predicted} ({length} pairs)')
    print(f'{count} recordings, {misfits} misfits')


def fit(paths):
    import numpy
    from scipy.optimize import linprog
    model = Model()
    rows, bounds = [], []
    index = {key: position for position, key in enumerate(PARAMS)}

    def coefficients(candidate):
        vector = numpy.zeros(len(PARAMS))
        for key in candidate[4]:
            vector[index[key]] += 1
        return vector

    for _, primer, energy, _ in recordings(paths):
        candidates = model.candidates(primer)
        if energy == '-':
            for candidate in candidates:  # no stable hairpin: every energy >= 0
                rows.append(coefficients(candidate))
                bounds.append(candidate[3])
            continue
        high = round(float(energy) * 100)
        low = high - 10                       # low < energy <= high
        best = model.hairpin(primer)
        vector = coefficients(best)
        rows.append(vector)
        bounds.append(best[3] - low - 1)      # integer hundredths
        rows.append(-vector)
        bounds.append(high - best[3])
        for candidate in candidates:          # no other hairpin below the shown one
            if candidate[:3] != best[:3]:
                rows.append(coefficients(candidate))
                bounds.append(candidate[3] - low - 1)
    matrix, limits = numpy.array(rows), numpy.array(bounds, dtype=float)

    def named(*names):
        vector = numpy.zeros(len(PARAMS))
        for name in names:
            vector[NAMES.index(name)] += 1
        return vector

    left = {'A:T': 'AG/TT', 'C:G': 'CG/GT', 'G:C': 'GG/CT', 'T:A': 'TA/GT'}
    right = {'A:T': 'GA/TT', 'C:G': 'GC/TG', 'G:C': 'CT/GG', 'T:A': 'AT/TG'}
    sums = [('GG/TT', named('GG/TT'))]
    sums += [(f'L({x}) + R({y})', named(left[x], right[y])) for x in left for y in right]
    sums += [(f'L({x}) + GT/TG + L({y})', named(left[x], 'GT/TG', left[y])) for x in left for y in left if x <= y]
    sums += [(f'R({x}) + TG/GT + R({y})', named(right[x], 'TG/GT', right[y])) for x in right for y in right if x <= y]
    sums += [('GT/TG + TG/GT', named('GT/TG', 'TG/GT'))]
    print(f'{len(rows)} inequalities')
    for name, vector in sums:
        ends = []
        for sign in (1, -1):
            result = linprog(sign * vector, A_ub=matrix, b_ub=limits, bounds=[(-1000, 1000)] * len(PARAMS),
                             method='highs')
            ends.append(None if result.status != 0 else round(sign * result.fun, 2))
        value = sum(WOBBLE[PARAMS[k]] * vector[k] for k in range(len(PARAMS)))
        print(f'{name}\t[{ends[0]}, {ends[1]}]\tproposed table: {value:.0f}')


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('check', 'fit'):
        sys.exit(__doc__)
    if sys.argv[1] == 'check':
        check(sys.argv[2:])
    else:
        fit(sys.argv[2:])
