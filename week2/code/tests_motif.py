import math
# import unittest
from python import Bio.Seq as Seq
# import __init__ as motifs
import unittest
import __init__ as motifs
import numpy as np

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
        self.assertEqual(len(m.alignment), 3)
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        self.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})

        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
            )
        )

        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.8186697601117167,
                        1.7369655941662063,
                        0.5419780939258206,
                        1.7369655941662063,
                        1.7369655941662063,
                    ]
                ),
            )
        )

        m.background = None
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        pseudocounts = math.sqrt(len(m.alignment))
        m.pseudocounts = {
            letter: m.background[letter] * pseudocounts for letter in "ACGT"
        }
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.3532586861097656,
                        0.7170228827697498,
                        0.11859369972847714,
                        0.7170228827697498,
                        0.7170228827697499,
                    ]
                ),
            )
        )

        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.19727984803857979,
                        0.561044044698564,
                        0.20984910512125132,
                        0.561044044698564,
                        0.5610440446985638,
                    ]
                ),
            )
        )
    
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