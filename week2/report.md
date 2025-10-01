Steps (rough outline from what i can remember rn)
1. Port matrix functions first - as other files depend on it - and test individual functions
2. Port minimal.py and thresholds.py, and test functions
3. Convert _pwm.c to codon
4. Port motif class
5. Write tests in one file

Gotchas
* we skip the python test for test_format_clusterbuster (codon test passes); it tries to access Motif.weight which is not possible (problem with biopython?)

Time Spent
3 - 6 hours every day since Monday Sep 22 (9 days)