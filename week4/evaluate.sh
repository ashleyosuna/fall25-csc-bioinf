#!/bin/bash
set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

printf "%-16s %-16s %-16s\n" "Method" "Language" "Runtime"
echo "--------------------------------------------------------"
python3 week4/code/tests.py

codon run -release week4/code/tests.py