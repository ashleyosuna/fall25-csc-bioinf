#!/bin/bash
set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week2"
codon run -release week2/code/tests_motif.py
python3 week2/code/tests_python.py