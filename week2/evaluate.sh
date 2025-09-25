#!/bin/bash
set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week2"
matrix=$(codon run week2/code/matrix.py)
echo $matrix