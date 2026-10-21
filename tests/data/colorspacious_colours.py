"""Write reference values for core/colour: the WCAG 2.2 relative
luminance and contrast ratio, the colour-vision simulation of Machado,
Oliveira and Fernandes (2009, IEEE TVCG 15(6):1291-1298; severity 1)
and the CIE 1976 colour difference, as colorspacious 1.1.2 computes the
last two (sRGB1+CVD with severity 100, deltaE in CIELab, D65).

usage: ~/.venvs/colour/bin/python -I colorspacious_colours.py > colorspacious_colours.tsv
(a virtual environment with colorspacious==1.1.2 and numpy)

Each line: the colour (#rrggbb), a second colour, the relative
luminance of the first, the contrast ratio of the two (WCAG 2.2
formulas, written here: colorspacious has no contrast), the
protanopia, deuteranopia and tritanopia simulations of the first
(sRGB components 0 to 1, clipped to [0, 1]), and the CIE76 difference
of the two. The colours are a fixed list: greys, primaries, and the
picking colours of the GUI (gui/picking_colours.hpp) at the time of
writing.
"""
import colorspacious
import numpy as np

assert colorspacious.__version__ == '1.1.2'

PAIRS = [
    ('#000000', '#ffffff'), ('#767676', '#ffffff'), ('#777777', '#ffffff'),
    ('#ff0000', '#00ff00'), ('#0000ff', '#ffff00'), ('#808080', '#404040'),
    ('#d4d4d4', '#8ccff0'), ('#555555', '#005f8f'), ('#1e1e1e', '#ffffff'),
    ('#c86432', '#3264c8'), ('#20a080', '#a02080'),
]


def rgb(colour):
    return np.array([int(colour[index:index + 2], 16) for index in (1, 3, 5)]) / 255


def luminance(colour):
    linear = np.where(colour <= 0.04045, colour / 12.92, ((colour + 0.055) / 1.055) ** 2.4)
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast(first, second):
    darker, lighter = sorted((luminance(first), luminance(second)))
    return (lighter + 0.05) / (darker + 0.05)


def simulate(colour, kind):
    space = {'name': 'sRGB1+CVD', 'cvd_type': kind, 'severity': 100}
    return np.clip(colorspacious.cspace_convert(colour, space, 'sRGB1'), 0, 1)


print('# colour\tsecond\tluminance\tcontrast\tprotanopia\tdeuteranopia\ttritanopia\tdelta_e76'
      ' (colorspacious 1.1.2; WCAG 2.2 formulas)')
for first_name, second_name in PAIRS:
    first, second = rgb(first_name), rgb(second_name)
    simulations = [','.join(f'{value:.6f}' for value in simulate(first, kind))
                   for kind in ('protanomaly', 'deuteranomaly', 'tritanomaly')]
    difference = colorspacious.deltaE(first, second, input_space='sRGB1', uniform_space='CIELab')
    print('\t'.join([first_name, second_name, f'{luminance(first):.6f}', f'{contrast(first, second):.6f}',
                     *simulations, f'{difference:.6f}']))
