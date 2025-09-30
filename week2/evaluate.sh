#!/bin/bash
set -euxo pipefail
# PATH=${PATH}:${HOME}/.codon/bin
export PATH="${HOME}/.codon/bin:$PATH"
echo "week2"
tests_motif=$(codon run week2/code/tests_motif.py)
echo $tests_motif