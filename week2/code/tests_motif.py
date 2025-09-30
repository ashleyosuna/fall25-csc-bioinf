import math
from python import Bio.Seq as Seq
import unittest
import __init__ as motifs
import numpy as np

class TestMotif(unittest.TestCase):
    def __init__(self):
        super().__init__()
    
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
    
    def test_format_transfac(self):
        m = motifs.create([Seq.Seq("ATATA")])
        m.name = "Foo"
        s = m.format(format_spec="transfac")
        expected_transfac = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(s, expected_transfac)

    def test_format_clusterbuster(self):
        m = motifs.create([Seq.Seq("ATATA")])
        m.name = "Foo"
        s = m.format(format_spec="clusterbuster")

        expected = """>Foo
1	0	0	0
0	0	0	1
1	0	0	0
0	0	0	1
1	0	0	0
"""

        self.assertEqual(s, expected)


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

        m.mask = "* * *"
        self.assertEqual('ATATA* * *\n', m.__str__(masked=True))

    def test_mask(self):
        m = motifs.create([Seq.Seq(("ATATA"))])
        self.assertEqual([1] * m.length, m.mask)

        m.mask = "* * *"
        self.assertEqual([1, 0, 1, 0, 1], m.mask)

        m.mask = [2, 0, 3, 0, 1]
        self.assertEqual([1, 0, 1, 0, 1], m.mask)

        def exception(): m.mask = "abcab"
        self.assertRaises(ValueError, exception)

        def exception(): m.mask = [1,2]
        self.assertRaises(ValueError, exception)

    def test_reverse_complement(self):
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        pseudocounts = 0.5
        m = motifs.create([Seq.Seq(("ATATA"))])
        m.background = background
        m.pseudocounts = pseudocounts

        received_forward = m.format(format_spec="transfac")
        expected_forward = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(received_forward, expected_forward)

        expected_forward_pwm = """        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        
        rc = m.reverse_complement()

        received_reverse = rc.format(format_spec="transfac")
        expected_reverse = """P0      A      C      G      T
01      0      0      0      1      T
02      1      0      0      0      A
03      0      0      0      1      T
04      1      0      0      0      A
05      0      0      0      1      T
XX
//
"""
        self.assertEqual(received_reverse, expected_reverse)

        expected_reverse_pwm = """        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(expected_reverse_pwm, str(rc.pwm))


        background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
        pseudocounts = 0.5
        m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        expected_forward_rna_counts = """        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   0.00   1.00   0.00   1.00   0.00
"""
        self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        
        expected_forward_rna_pwm = """        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
        expected_reverse_rna_counts = """        0      1      2      3      4
A:   0.00   1.00   0.00   1.00   0.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   1.00   0.00   1.00   0.00   1.00
"""
        self.assertEqual(
            str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
        )
        expected_reverse_rna_pwm = """        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

        m = motifs.create([Seq.Seq("ATATA")])
        counts = m.counts
        m = motifs.Motif(counts=counts)
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = m.format(format_spec="transfac")
        self.assertEqual(received_forward, expected_forward)
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = m.format(format_spec="transfac")
        self.assertEqual(received_reverse, expected_reverse)
        self.assertEqual(str(m.pwm), expected_reverse_pwm)
        # Same, but for RNA count matrix
        m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
        counts = m_rna.counts
        m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
        self.assertEqual(
            str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
        )
        self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

    @test
    def test_consensus(self):
        m = motifs.create([
            Seq.Seq("ATATA"),
            Seq.Seq("ATCTA"),
            Seq.Seq("TTGTA"),
            Seq.Seq("ATATA"),
        ])

        expected_consensus = "ATATA"
        self.assertEqual(str(m.consensus), expected_consensus)

    @test
    def test_anticonsensus(self):
        m = motifs.create([
            Seq.Seq("ATCGA"),
            Seq.Seq("ATCGA"),
            Seq.Seq("GGGTG"),
            Seq.Seq("GGGTG"),
            Seq.Seq("CCACC"),
            Seq.Seq("CCACC"),
            Seq.Seq("TATAT")
        ])

        expected_anticonsensus = "TATAT"
        self.assertEqual(str(m.anticonsensus), expected_anticonsensus)

    def test_degenerate_consensus(self):
        m = motifs.create([
            Seq.Seq("ATATA"),
            Seq.Seq("ATCTA"),
            Seq.Seq("TTGTA"),
            Seq.Seq("ATGTA")
        ])
        # Position-wise breakdown:
        # 1: A (3), T (1) → A
        # 2: T (4) → T
        # 3: A (1), C (1), G (2) → V (A/C/G)
        # 4: T (4) → T
        # 5: A (4) → A

        expected_degenerate_consensus = "ATVTA"
        self.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)

    def test_degenerate_consensus_with_ties(self):
        m = motifs.create([
            Seq.Seq("A"),
            Seq.Seq("C"),
            Seq.Seq("G"),
            Seq.Seq("T"),
        ])
        expected_degenerate_consensus = "N"  # all bases equally represented
        self.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)

    def test_degenerate_consensus_rna(self):
        m_rna = motifs.create([
            Seq.Seq("AUAUA"),
            Seq.Seq("AUCUA"),
            Seq.Seq("UUGUA"),
            Seq.Seq("AUGUA")
        ], alphabet="ACGU")

        expected_degenerate_consensus = "AUVUA"
        self.assertEqual(str(m_rna.degenerate_consensus), expected_degenerate_consensus)

    def test_getitem(self):
        m = motifs.create(["AACGCCA", "ACCGCCC", "AACTCCG"])
        expected="""AACGCC\nACCGCC\nAACTCC"""
        self.assertEqual(str(m[:-1]), expected)

        expected=""
        self.assertEqual(str(m[0:0]), expected)

        expected="A\nA\nA"
        self.assertEqual(str(m[0:1]), expected)
        
        # EXCEPTIONS
        self.assertRaises(TypeError, lambda: m[0])
        self.assertRaises(TypeError, lambda: m['A'])

tests = TestMotif()

# shitty test runner
def run_test(name, func):
    try:
        func()
        print(f"{name} passed.")
    except AssertionError as e:
        print(f"{name} failed: {e}")

run_test("test_format", tests.test_format)
run_test("test_relative_entropy_alignment", tests.test_relative_entropy_alignment)
run_test("test_relative_entropy_counts", tests.test_relative_entropy_counts)
run_test("test_pwm", tests.test_pwm)
run_test("test_pssm", tests.test_pssm)
run_test("test_str", tests.test_str)
run_test("test_mask", tests.test_mask)
run_test("test_reverse_complement", tests.test_reverse_complement)
run_test("test_consensus", tests.test_consensus)
run_test("test_anticonsensus", tests.test_anticonsensus)
run_test("test_degenerate_consensus", tests.test_degenerate_consensus)
run_test("test_degenerate_consensus_with_ties", tests.test_degenerate_consensus_with_ties)
run_test("test_degenerate_consensus_rna", tests.test_degenerate_consensus_rna)
run_test("test_getitem", tests.test_getitem)
run_test("test_format_transfac", tests.test_format_transfac)
run_test("test_format_transfac", tests.test_format_clusterbuster)