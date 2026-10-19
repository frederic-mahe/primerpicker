#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "--config FILE" (a pd2.cfg settings file) and of
## "find-pair --offsets". Expected values come from PRIMER.exe
## (tests/legacy/search/find_pair_limits.tsv; tests/legacy/README.md,
## "Settings file pd2.cfg").
## Usage: test_options.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test options"
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
TLONG="${LEGACY}/enzymes/TLONG.GB"
RELAXED="${LEGACY}/config/PD2_RELAXED.CFG"

## errors
DESCRIPTION="--config needs a readable file"
"${PRIMERPICKER}" --config "${LEGACY}/config/MISSING.CFG" pair GTAAAACGACGGCCAGT AACAGCTATGACCATG \
                  > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="--config rejects a file that is not a pd2.cfg"
"${PRIMERPICKER}" --config "${TLONG}" pair GTAAAACGACGGCCAGT AACAGCTATGACCATG > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="--config without a command is an error"
"${PRIMERPICKER}" --config "${RELAXED}" > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="--offsets needs two numbers"
"${PRIMERPICKER}" find-pair --offsets 40 "${TLONG}" 400 800 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="--offsets rejects letters"
"${PRIMERPICKER}" find-pair --offsets a,40 "${TLONG}" 400 800 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## counts of a search, as "A B pairs"
counts () {
    "${PRIMERPICKER}" "$@" | \
        sed -n 's/^Primer A: *\([0-9]*\).*/\1/p; s/^Primer B: *\([0-9]*\).*/\1/p; s/^Pairs: *\([0-9]*\).*/\1/p' | \
        tr '\n' ' '
}

## built-in criteria, 111 positions for primer A: 9, 16 and 48
DESCRIPTION="find-pair --offsets: windows of the form"
[[ $(counts find-pair --offsets 110,50 "${TLONG}" 200 800) == "9 16 48 " ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## relaxed settings, 99 possible primers A: 99, 6 and 220
DESCRIPTION="--config: the settings file's criteria apply"
[[ $(counts --config "${RELAXED}" find-pair --offsets 98,5 "${TLONG}" 400 800) == "99 6 220 " ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="--offsets after --weights"
[[ $(counts --config "${RELAXED}" find-pair --weights 11111 --offsets 98,5 "${TLONG}" 400 800) == "99 6 220 " ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## relaxed settings, 100 possible primers A
DESCRIPTION="find-pair: 100 possible primers stop the search"
"${PRIMERPICKER}" --config "${RELAXED}" find-pair --offsets 99,5 "${TLONG}" 400 800 | \
    grep -qxF "More than 100 sites found.  Suggest you redefine matching limits." && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## relaxed settings, the file's windows (100 and 100): 101 positions
DESCRIPTION="--config: the settings file's search windows apply"
"${PRIMERPICKER}" --config "${RELAXED}" find-pair "${TLONG}" 400 800 | \
    grep -qxF "More than 100 sites found.  Suggest you redefine matching limits." && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## relaxed settings, 61 x 61 positions: the count stops at 1000 pairs
## (version 1 then warns, where PRIMER.exe hangs)
DESCRIPTION="find-pair: the pairs stop at 1000"
[[ $(counts --config "${RELAXED}" find-pair --offsets 60,60 "${TLONG}" 400 800) == "61 61 1000 " ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="find-pair: a warning follows the capped count"
"${PRIMERPICKER}" --config "${RELAXED}" find-pair --offsets 60,60 "${TLONG}" 400 800 | \
    grep -qxF "More than 1000 possible pairs found!" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## PD2_CRITERIA.CFG: GC 45-65 %, Tm 50-75 °C
DESCRIPTION="--config: analyze shows the settings file's criteria"
"${PRIMERPICKER}" --config "${LEGACY}/config/PD2_CRITERIA.CFG" analyze GTAAAACGACGGCCAGT | \
    grep -qE "^% GC +min 45, max 65 " && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
