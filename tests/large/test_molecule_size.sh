#!/bin/bash -
# shellcheck disable=SC2015

## Opt-in test of the largest molecules: a FASTA file of 4,300,000,100
## bases (more than 2^32) with a NotI site planted at 4,300,000,001;
## "find-enzymes" must report that position. Needs about 4.4 GB of disk
## and 13 GB of memory, and a minute or more: enabled with the CMake
## option PRIMERPICKER_LARGE_TESTS (ctest -L large).
## Usage: test_molecule_size.sh /path/to/primerpicker

## Print a header
SCRIPT_NAME="Test molecule size"
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
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT
MOLECULE="${WORKDIR}/large.fa"

## random bases, 60 per line, with GCGGCCGC (NotI) at 4,300,000,001
python3 -I - "${MOLECULE}" <<'EOP'
import random, sys
total, site = 4_300_000_100, 4_300_000_001
table = bytes.maketrans(bytes(range(256)), b'acgt' * 64)
rnd = random.Random(9)
written, block, line = 0, 6_000_000, 60
with open(sys.argv[1], 'wb') as out:
    out.write(b'>LARGE random\n')
    while written < total:
        size = min(block, total - written)
        data = bytearray(rnd.randbytes(size).translate(table))
        if written < site <= written + size - 8:
            data[site - written - 1:site - written + 7] = b'gcggccgc'
        out.write(b'\n'.join(bytes(data[i:i + line]) for i in range(0, size, line)) + b'\n')
        written += size
EOP

DESCRIPTION="find-enzymes: a site beyond 2^32 keeps its position"
"${PRIMERPICKER}" find-enzymes "${MOLECULE}" NotI | \
    grep -qx "4300000001" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
