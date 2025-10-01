Steps (rough outline from what i can remember rn)

1. Port matrix functions first - as other files depend on it - and test individual functions
2. Port minimal.py and thresholds.py, and test functions
3. Convert \_pwm.c to codon
4. Port motif class
5. Write tests in one file

# Gotchas

## Tests

- We skip the python test for test_format_clusterbuster (codon test passes); it tries to access Motif.weight which is not possible (problem with biopython?)

- When testing test_pwm_getitem, we were unable to execute some tests that the original biopython tests included, for example, the following code
  ```
    t = counts[i2, i3:i12:i2]
    tc.assertAlmostEqual(t[i0], 0.0)
  ```
  (where i2, i3, i12, i2, and i0 are all integers)
  since t is of type dict[str, list[float]] and codon does not allow to index dictionaries of this type with integers.

Time Spent
3 - 6 hours every day since Monday Sep 22 (9 days)
