#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of molecule files in PRIMER.exe's own format (written
## by Molecule-Ops -> Save/Remove -> S). Expected values come from
## PRIMER.exe (tests/legacy/molecule_save/, recorded 2026-10-22;
## tests/legacy/PROBES_20261022_molecule_save.md).
## Usage: test_native_molecule.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test native molecule files"
LINE=$(printf '%.0s-' {1..76})
printf "# %s %s\n" "${LINE:${#SCRIPT_NAME}}" "${SCRIPT_NAME}"

## Declare a color code for test results
RED="\033[1;31m"
GREEN="\033[1;32m"
NO_COLOR="\033[0m"

failure () {
    printf "%bFAIL%b: %s\n" "${RED}" "${NO_COLOR}" "${1}"
    exit 1
}

success () {
    printf "%bPASS%b: %s\n" "${GREEN}" "${NO_COLOR}" "${1}"
}

PRIMERPICKER="${1}"
LEGACY="${2}"
WRITTEN="${LEGACY}/molecule_save/written"
INPUTS="${LEGACY}/molecule_save/inputs"

## PRIMER.exe reads the files it writes: ZZ.GB (TEST600.GB saved) gives
## the same results as TEST600.GB
DESCRIPTION="a molecule saved by PRIMER.exe gives the results of the original file"
diff <("${PRIMERPICKER}" find-enzymes "${WRITTEN}/save_new.ZZ.GB" TaqI 2>&1) \
     <("${PRIMERPICKER}" find-enzymes "${LEGACY}/pairs/TEST600.GB" TaqI 2>&1) > /dev/null && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## a malformed file: PRIMER.exe's "Illegal File Format" (hand9)
DESCRIPTION="a malformed file of PRIMER.exe's format is an Illegal File Format"
"${PRIMERPICKER}" find-enzymes "${INPUTS}/ZH9" EcoRI 2>&1 > /dev/null | \
    grep -q "Illegal File Format" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="a malformed file of PRIMER.exe's format fails"
"${PRIMERPICKER}" find-enzymes "${INPUTS}/ZH9" EcoRI > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

exit 0
