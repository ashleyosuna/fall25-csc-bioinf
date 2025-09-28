import math
# import unittest
from python import Bio.Seq as Seq
# import __init__ as motifs
import unittest
import __init__ as motifs

class TestMotif(unittest.TestCase):
    passed: int
    total: int

    def __init__(self):
        super().__init__()
        self.passed = 0
        self.total = 0
    
    def test_format(self):
        res = self.assertEqual(True, True)

        m = motifs.create([Seq.Seq("ATATA")])
        m.name = "Foo"
        s1 = m.format(format_spec="pfm")

        expected_pfm = """  1.00   0.00   1.00   0.00   1.00
  0.00   0.00   0.00   0.00   0.00
  0.00   0.00   0.00   0.00   0.00
  0.00   1.00   0.00   1.00   0.00
"""

        s2 = m.format(format_spec="jaspar")
        expected_jaspar = """>None Foo
A [  1.00   0.00   1.00   0.00   1.00]
C [  0.00   0.00   0.00   0.00   0.00]
G [  0.00   0.00   0.00   0.00   0.00]
T [  0.00   1.00   0.00   1.00   0.00]
"""
        self.assertEqual(s1, expected_pfm)
        self.assertEqual(s2, expected_jaspar)
        self.assertRaises(ValueError, lambda : m.format(format_spec="foo_bar"))

        if not res:
            self.passed += 1
        self.total += 1

    def test_relative_entropy(self):
        m = motifs.create([Seq.Seq("ATATA"), Seq.Seq("ATCTA"), Seq.Seq("TTGTA")])
        # print(m.alignment.length)
        self.assertEqual(len(m.alignment), 3)
        # self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    
    def __str__(self):
        return f"{self.passed} tests passed out of {self.total} tests"
    
    # def test_format(self):
    #     m = motifs.create([Seq("ATATA")])
        # m = matrix.FrequencyPositionMatrix(alphabet="ACGT", values={'A': [1.0, 2.0, 3.0]})

# runner = unittest.TextTestRunner(verbosity=2)
# unittest.main(testRunner=runner)
tests = TestMotif()
tests.test_format()
tests.test_relative_entropy()
print(tests)