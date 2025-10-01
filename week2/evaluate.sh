#!/bin/bash
set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week2"
echo "running Codon tests"
tests_motif=$(codon run week2/code/tests_motif.py)
echo "running Python tests"
tests_motif=$(python3 run week2/code/tests_motif.py)