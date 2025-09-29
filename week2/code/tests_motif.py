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
        self.assertRaises(ValueError, lambda : m.format(format_spec="foo_bar"))

    def test_relative_entropy_alignment(self):
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
    
    def test_relative_entropy_counts(self):
        m = motifs.Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.08170417,
                        0.02904941,
                        0.01477186,
                    ]
                ),
            )
        )
    
    def test_pwm(self):
        m = motifs.create([Seq.Seq(("ATATA"))])
        expected = """        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
T:   0.00   1.00   0.00   1.00   0.00
"""
        self.assertEqual(expected, m.pwm.__str__())

    def test_pssm(self):
        m = motifs.create([Seq.Seq(("ATATA"))])
        expected="""        0      1      2      3      4
A:   2.00   -inf   2.00   -inf   2.00
C:   -inf   -inf   -inf   -inf   -inf
G:   -inf   -inf   -inf   -inf   -inf
T:   -inf   2.00   -inf   2.00   -inf
"""
        self.assertEqual(expected, m.pssm.__str__())

    def test_str(self):
        m = motifs.create([Seq.Seq(("ATATA"))])
        self.assertEqual("ATATA", m.__str__())
    
    def __str__(self):
        return f"{self.passed} tests passed out of {self.total} tests"

tests = TestMotif()
tests.test_format()
tests.test_relative_entropy_alignment()
tests.test_relative_entropy_counts()
tests.test_pwm()
tests.test_pssm()
tests.test_str()
# print(tests)