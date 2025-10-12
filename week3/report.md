# GOTCHAS

## Defining Custom Exception Types

When re-writing the TreeNode class into Codon, in the original code, custom exceptions called TreeException are raised in certain cases,
however, when I tried to create the class TreeException such that it was a subclass of the Exception class, the program would not compile. To work around this, I switched these exceptions to ValueError exceptions.

## Parametrize f string's precision for float numbers

In TreeNode's class, the **str** function takes a round distance parameter, to specify how many positions after the decimal point should be included in the string, however, codon does not seem to support writing f"{self.\_distance:.{<variable>}f}", so instead I hardcoded the arbitrary value of 3 here.

## numpy load text

The following line did not work for me, so I instead had to read the data file myself.
distances = np.loadtxt("week3/data/distances.txt", dtype=float)

## using frozenset in **eq** function

Since Codon does not have frozenset, I tried using a regular set instead (and extending the class as suggested in Piazza), however, I kept getting errors and could not make it work; instead, I modified the code in **eq** to check if two nodes have the same number of children, and if all of their children are the same (regardless of order).

## pytest not supported

Since I believe pytest is not supported in Codon, I used unittest for both codon and python tests instead, but kept the same tests as indicated.
