# Recorded behaviour of the legacy program

Reference observations of *Primer Designer v2.0* (Scientific &
Educational Software, 1990-91), `PRIMER.exe`. The executable itself
is proprietary and is not part of this repository.

Recordings are made with DOSBox-X 2026.01.02, running headless:
keystrokes are injected with the built-in `AUTOTYPE` command, and
the screen is recorded with `dx-capture /V` (ZMBV video). The program
runs in 80x25 text mode (720x400 pixels).

- `tools/zmbv.py`: decodes a ZMBV AVI file into PNG images of the
  distinct frames (needs Python 3 and Pillow).
- `tools/system_defaults.conf`: DOSBox-X configuration that records
  the System Defaults screens (usage in the file header).


## Built-in defaults (recorded 2026-10-17)

Fresh start, no `pd2.cfg` file present. Screens: Modify -> System
Defaults.

| Screen | Setting | Default |
|---|---|---|
| `screens/system_defaults.png` | Primer length (incl. 5' extension) | 20 bases |
| | False priming and uniqueness search: exclude Tm <= 0 °C | Y |
| | Sort homologies by Tm | Y |
| | Match primer GC% to within ± | 5 % |
| | Match primer Tm to within ± | 10 °C |
| | Match Active Primer search range | start -50, end +50 bases from target |
| | Find Pair search range | normal start -100, complement end +100 bases from target |
| `screens/criteria_defaults.png` | GC% range | 50-60 % |
| | Tm range | 55-80 °C |
| | Hairpins: energy cutoff | 0.0 kcal |
| | Dimers: reject if >= | 3 matches at 3' end |
| | Dimers: reject if >= | 7 adjacent homologous bases |
| | Runs: reject runs of >= | 3 bases (any) |
| | Runs: reject runs of >= | 3 G or C at 3' end |
| `screens/weights_defaults.png` | Find-Pair weighting factors (0 = ignore, 1-5): GC mismatch, Tm mismatch, 3' dimers, any dimers, false priming | 1 each |
