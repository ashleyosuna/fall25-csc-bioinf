import math
import unittest
import numpy as np
import __init__ as motifs
from python import Bio.Seq as Seq

@test
def test_format():
        tc = unittest.TestCase()
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
        tc.assertEqual(s1, expected_pfm)
        tc.assertEqual(s2, expected_jaspar)
        tc.assertRaises(ValueError, lambda : m.format(format_spec="foo_bar"))
test_format()

@test
def test_format_transfac():
    tc = unittest.TestCase()
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
    tc.assertEqual(s, expected_transfac)
test_format_transfac()

@test
def test_format_clusterbuster():
    tc = unittest.TestCase()
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
    tc.assertEqual(s, expected)
test_format_clusterbuster()

@test
def test_relative_entropy_alignment():
    m = motifs.create([Seq.Seq("ATATA"), Seq.Seq("ATCTA"), Seq.Seq("TTGTA")])
    tc = unittest.TestCase()
    tc.assertEqual(len(m.alignment), 3)
    tc.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    tc.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})

    tc.assertTrue(
        np.allclose(
            m.relative_entropy,
            np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
        )
    )

    m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
    tc.assertTrue(
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
    tc.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    pseudocounts = math.sqrt(len(m.alignment))
    m.pseudocounts = {
        letter: m.background[letter] * pseudocounts for letter in "ACGT"
    }
    tc.assertTrue(
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
    tc.assertTrue(
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
test_relative_entropy_alignment()

@test
def test_relative_entropy_counts():
    m = motifs.Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
    tc = unittest.TestCase()
    tc.assertTrue(
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
test_relative_entropy_counts()

@test
def test_pwm():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    expected = """        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
T:   0.00   1.00   0.00   1.00   0.00
"""
    tc.assertEqual(expected, m.pwm.__str__())
test_pwm()

@test
def test_pssm():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    expected="""        0      1      2      3      4
A:   2.00   -inf   2.00   -inf   2.00
C:   -inf   -inf   -inf   -inf   -inf
G:   -inf   -inf   -inf   -inf   -inf
T:   -inf   2.00   -inf   2.00   -inf
"""
    tc.assertEqual(expected, m.pssm.__str__())
test_pssm()

@test
def test_str():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    tc.assertEqual("ATATA", m.__str__())

    m.mask = "* * *"
    tc.assertEqual('ATATA* * *\n', m.__str__(masked=True))
test_str()

@test
def test_mask():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    tc.assertEqual([1] * m.length, m.mask)

    m.mask = "* * *"
    tc.assertEqual([1, 0, 1, 0, 1], m.mask)

    m.mask = [2, 0, 3, 0, 1]
    tc.assertEqual([1, 0, 1, 0, 1], m.mask)

    def exception(): m.mask = "abcab"
    tc.assertRaises(ValueError, exception)

    def exception(): m.mask = [1,2]
    tc.assertRaises(ValueError, exception)
test_mask()

@test
def test_reverse_complement():
    tc = unittest.TestCase()
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
    tc.assertEqual(received_forward, expected_forward)

    expected_forward_pwm = """        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
    tc.assertEqual(str(m.pwm), expected_forward_pwm)
    
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
    tc.assertEqual(received_reverse, expected_reverse)

    expected_reverse_pwm = """        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
    tc.assertEqual(expected_reverse_pwm, str(rc.pwm))


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
    tc.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
    
    expected_forward_rna_pwm = """        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.17   0.50   0.17   0.50   0.17
"""
    tc.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
    expected_reverse_rna_counts = """        0      1      2      3      4
A:   0.00   1.00   0.00   1.00   0.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   1.00   0.00   1.00   0.00   1.00
"""
    tc.assertEqual(
        str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
    )
    expected_reverse_rna_pwm = """        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.50   0.17   0.50   0.17   0.50
"""
    tc.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

    m = motifs.create([Seq.Seq("ATATA")])
    counts = m.counts
    m = motifs.Motif(counts=counts)
    m.background = background
    m.pseudocounts = pseudocounts
    received_forward = m.format(format_spec="transfac")
    tc.assertEqual(received_forward, expected_forward)
    tc.assertEqual(str(m.pwm), expected_forward_pwm)
    m = m.reverse_complement()
    received_reverse = m.format(format_spec="transfac")
    tc.assertEqual(received_reverse, expected_reverse)
    tc.assertEqual(str(m.pwm), expected_reverse_pwm)
    # Same, but for RNA count matrix
    m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
    counts = m_rna.counts
    m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
    m_rna.background = background_rna
    m_rna.pseudocounts = pseudocounts
    tc.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
    tc.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
    tc.assertEqual(
        str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
    )
    tc.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)
test_reverse_complement()

@test
def test_consensus():
    tc = unittest.TestCase()
    m = motifs.create([
        Seq.Seq("ATATA"),
        Seq.Seq("ATCTA"),
        Seq.Seq("TTGTA"),
        Seq.Seq("ATATA"),
    ])

    expected_consensus = "ATATA"
    tc.assertEqual(str(m.consensus), expected_consensus)
test_consensus()

@test
def test_anticonsensus():
    tc = unittest.TestCase()
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
    tc.assertEqual(str(m.anticonsensus), expected_anticonsensus)
test_anticonsensus()

@test
def test_degenerate_consensus():
    tc = unittest.TestCase()
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
    tc.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)
test_degenerate_consensus()

@test
def test_degenerate_consensus_with_ties():
    tc = unittest.TestCase()
    m = motifs.create([
        Seq.Seq("A"),
        Seq.Seq("C"),
        Seq.Seq("G"),
        Seq.Seq("T"),
    ])
    expected_degenerate_consensus = "N"  # all bases equally represented
    tc.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)
test_degenerate_consensus_with_ties()

@test
def test_degenerate_consensus_rna():
    tc = unittest.TestCase()
    m_rna = motifs.create([
        Seq.Seq("AUAUA"),
        Seq.Seq("AUCUA"),
        Seq.Seq("UUGUA"),
        Seq.Seq("AUGUA")
    ], alphabet="ACGU")

    expected_degenerate_consensus = "AUVUA"
    tc.assertEqual(str(m_rna.degenerate_consensus), expected_degenerate_consensus)
test_degenerate_consensus_rna()

@test
def test_getitem():
    tc = unittest.TestCase()
    m = motifs.create(["AACGCCA", "ACCGCCC", "AACTCCG"])
    expected="""AACGCC\nACCGCC\nAACTCC"""
    tc.assertEqual(str(m[:-1]), expected)

    expected="\n\n"
    tc.assertEqual(str(m[0:0]), expected)

    expected="A\nA\nA"
    tc.assertEqual(str(m[0:1]), expected)
    
    # EXCEPTIONS
    tc.assertRaises(TypeError, lambda: m[0])
    tc.assertRaises(TypeError, lambda: m['A'])
test_getitem()

@test
def test_minimal_parser_1():
    tc = unittest.TestCase()
    """Parse motifs/minimal_test.meme file."""
    with open("week2/tests/minimal_test.meme") as stream:
        record = motifs.parse(stream, "minimal")
    tc.assertEqual(record.version, "4")
    tc.assertEqual(record.alphabet, "ACGT")
    tc.assertEqual(len(record.sequences), 0)
    tc.assertEqual(record.command, "")
    tc.assertEqual(len(record), 3)
    motif = record[0]
    tc.assertEqual(motif.name, "KRP")
    # tc.assertEqual(record["KRP"], motif)
    # tc.assertEqual(motif.num_occurrences, 17)
    tc.assertEqual(motif.length, 19)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 4.1e-09, places=10)
    tc.assertEqual(motif.alphabet, "ACGT")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "TGTGATCGAGGTCACACTT")
    tc.assertEqual(motif.degenerate_consensus, "TGTGANNNWGNTCACAYWW")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    1.1684297174927525,
                    0.9432809925744818,
                    1.4307101633876265,
                    1.1549413780465179,
                    0.9308256303218774,
                    0.009164393966550805,
                    0.20124190687894253,
                    0.17618542656995528,
                    0.36777933103380855,
                    0.6635834532368525,
                    0.07729943368061855,
                    0.9838293592717438,
                    1.72489868427398,
                    0.8397561713453014,
                    1.72489868427398,
                    0.8455332015343343,
                    0.3106481207768122,
                    0.7382733641762232,
                    0.537435993300495,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "TGATCGA")
    motif = record[1]
    tc.assertEqual(motif.name, "IFXA")
    # tc.assertEqual(record["IFXA"], motif)
    # tc.assertEqual(motif.num_occurrences, 14)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 3.2e-35, places=36)
    tc.assertEqual(motif.alphabet, "ACGT")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
    tc.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.9632889858595118,
                    1.02677956765017,
                    2.451526420551951,
                    1.7098384161433415,
                    2.2598671267551107,
                    1.7098384161433415,
                    1.02677956765017,
                    1.391583804103081,
                    1.02677956765017,
                    1.1201961888781142,
                    0.27822438781180836,
                    0.36915366971717867,
                    1.7240522753630425,
                    0.3802185945622609,
                    0.790937683007783,
                    2.451526420551951,
                    1.7240522753630425,
                    1.3924085743645374,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CTGTATA")
test_minimal_parser_1()

@test
def test_minimal_parser_2():
    tc = unittest.TestCase()
    with open("week2/tests/minimal_test.meme") as stream:
        record = motifs.parse(stream, "minimal")
    motif = record[2]
    tc.assertEqual(motif.name, "IFXA_no_nsites_no_evalue")
    # tc.assertEqual(record["IFXA_no_nsites_no_evalue"], motif)
    # tc.assertEqual(motif.num_occurrences, 20)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 0.0, places=36)
    tc.assertEqual(motif.alphabet, "ACGT")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
    tc.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.99075309,
                    1.16078104,
                    2.45152642,
                    1.70983842,
                    2.25986713,
                    1.70983842,
                    1.16078104,
                    1.46052586,
                    1.16078104,
                    1.10213019,
                    0.29911041,
                    0.36915367,
                    1.72405228,
                    0.37696488,
                    0.85258086,
                    2.45152642,
                    1.72405228,
                    1.42793329,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CTGTATA")
test_minimal_parser_2()

@test
def test_minimal_parser_rna():
    tc = unittest.TestCase()
    """Test if Bio.motifs can parse MEME output files using RNA."""
    with open("week2/tests/minimal_test_rna.meme") as stream:
        record = motifs.parse(stream, "minimal")
    tc.assertEqual(record.version, "4")
    tc.assertEqual(record.alphabet, "ACGU")
    tc.assertEqual(len(record.sequences), 0)
    tc.assertEqual(record.command, "")
    tc.assertEqual(len(record), 3)
    motif = record[0]
    tc.assertEqual(motif.name, "KRP_fake_RNA")
    # tc.assertEqual(record["KRP_fake_RNA"], motif)
    # tc.assertEqual(motif.num_occurrences, 17)
    tc.assertEqual(motif.length, 19)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 4.1e-09, places=10)
    tc.assertEqual(motif.alphabet, "ACGU")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "UGUGAUCGAGGUCACACUU")
    tc.assertEqual(motif.degenerate_consensus, "UGUGANNNWGNUCACAYWW")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    1.1684297174927525,
                    0.9432809925744818,
                    1.4307101633876265,
                    1.1549413780465179,
                    0.9308256303218774,
                    0.009164393966550805,
                    0.20124190687894253,
                    0.17618542656995528,
                    0.36777933103380855,
                    0.6635834532368525,
                    0.07729943368061855,
                    0.9838293592717438,
                    1.72489868427398,
                    0.8397561713453014,
                    1.72489868427398,
                    0.8455332015343343,
                    0.3106481207768122,
                    0.7382733641762232,
                    0.537435993300495,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "UGAUCGA")
    motif = record[1]
    tc.assertEqual(motif.name, "IFXA_fake_RNA")
    # tc.assertEqual(record["IFXA_fake_RNA"], motif)
    # tc.assertEqual(motif.num_occurrences, 14)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 3.2e-35, places=36)
    tc.assertEqual(motif.alphabet, "ACGU")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "UACUGUAUAUAUAUCCAG")
    tc.assertEqual(motif.degenerate_consensus, "UACUGUAUAUAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.9632889858595118,
                    1.02677956765017,
                    2.451526420551951,
                    1.7098384161433415,
                    2.2598671267551107,
                    1.7098384161433415,
                    1.02677956765017,
                    1.391583804103081,
                    1.02677956765017,
                    1.1201961888781142,
                    0.27822438781180836,
                    0.36915366971717867,
                    1.7240522753630425,
                    0.3802185945622609,
                    0.790937683007783,
                    2.451526420551951,
                    1.7240522753630425,
                    1.3924085743645374,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CUGUAUA")

    motif = record[2]
    tc.assertEqual(motif.name, "IFXA_no_nsites_no_evalue_fake_RNA")
    # tc.assertEqual(record["IFXA_no_nsites_no_evalue_fake_RNA"], motif)
    # tc.assertEqual(motif.num_occurrences, 20)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 0.0, places=36)
    tc.assertEqual(motif.alphabet, "ACGU")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "UACUGUAUAUAUAUCCAG")
    tc.assertEqual(motif.degenerate_consensus, "UACUGUAUAUAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.99075309,
                    1.16078104,
                    2.45152642,
                    1.70983842,
                    2.25986713,
                    1.70983842,
                    1.16078104,
                    1.46052586,
                    1.16078104,
                    1.10213019,
                    0.29911041,
                    0.36915367,
                    1.72405228,
                    0.37696488,
                    0.85258086,
                    2.45152642,
                    1.72405228,
                    1.42793329,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CUGUAUA")
test_minimal_parser_rna()


# def test_pwm_getitem(self):
#     counts_ = {'A': [2.0, 9.0, 0.0, 1.0, 32.0, 3.0, 46.0, 1.0, 43.0, 15.0, 2.0, 2.0], 'C': [1.0, 33.0, 45.0, 45.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 'G': [39.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 43.0], 'T': [4.0, 2.0, 0.0, 0.0, 13.0, 42.0, 0.0, 45.0, 3.0, 30.0, 0.0, 0.0]}
#     m = motifs.Motif(counts=counts_)
#     counts = m.counts
#     python_integers = range(13)
#     numpy_integers = np.array(python_integers)
#     integers = {"python": python_integers, "numpy": numpy_integers}
#     for int_type in ("python", "numpy"):
#         i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12 = integers[int_type]
#         msg = f"using {int_type} integers as indices"
#         # slice, slice
#         d = counts[i1::i2, i2:i12:i3]
#         self.assertIsInstance(d, dict, msg=msg)
#         self.assertEqual(len(d), 2, msg=msg)
#         self.assertEqual(len(d["C"]), 4, msg=msg)
#         self.assertEqual(len(d["T"]), 4, msg=msg)
#         self.assertAlmostEqual(d["C"][i0], 45.0, msg=msg)
#         self.assertAlmostEqual(d["C"][i1], 1.0, msg=msg)
#         self.assertAlmostEqual(d["C"][i2], 0.0, msg=msg)
#         self.assertAlmostEqual(d["C"][i3], 1.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i0], 0.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i1], 42.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i2], 3.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i3], 0.0, msg=msg)
#         # slice, int
#         d = counts[i1::i2, i4]
#         self.assertIsInstance(d, dict, msg=msg)
#         self.assertEqual(len(d), 2, msg=msg)
#         self.assertAlmostEqual(d["C"], 1.0, msg=msg)
#         self.assertAlmostEqual(d["T"], 13.0, msg=msg)
#         # int, slice
#         t = counts[i2, i3:i12:i2]
#         self.assertIsInstance(t, tuple, msg=msg)
#         self.assertAlmostEqual(t[i0], 0.0, msg=msg)
#         self.assertAlmostEqual(t[i1], 0.0, msg=msg)
#         self.assertAlmostEqual(t[i2], 0.0, msg=msg)
#         self.assertAlmostEqual(t[i3], 0.0, msg=msg)
#         self.assertAlmostEqual(t[i4], 43.0, msg=msg)
#         # int, int
#         v = counts[i1, i5]
#         self.assertAlmostEqual(v, 1.0, msg=msg)
#         # tuple, slice
#         d = counts[(i0, i3), i3:i12:i2]
#         self.assertIsInstance(d, dict, msg=msg)
#         self.assertEqual(len(d), 2, msg=msg)
#         self.assertEqual(len(d["A"]), 5, msg=msg)
#         self.assertEqual(len(d["T"]), 5, msg=msg)
#         self.assertAlmostEqual(d["A"][i0], 1.0, msg=msg)
#         self.assertAlmostEqual(d["A"][i1], 3.0, msg=msg)
#         self.assertAlmostEqual(d["A"][i2], 1.0, msg=msg)
#         self.assertAlmostEqual(d["A"][i3], 15.0, msg=msg)
#         self.assertAlmostEqual(d["A"][i4], 2.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i0], 0.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i1], 42.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i2], 45.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i3], 30.0, msg=msg)
#         self.assertAlmostEqual(d["T"][i4], 0.0, msg=msg)
#         # tuple, int
#         d = counts[(i0, i3), i5]
#         self.assertIsInstance(d, dict, msg=msg)
#         self.assertEqual(len(d), 2, msg=msg)
#         self.assertAlmostEqual(d["A"], 3.0, msg=msg)
#         self.assertAlmostEqual(d["T"], 42.0, msg=msg)
#         # str, slice
#         t = counts["C", i2:i12:i4]
#         self.assertIsInstance(t, tuple, msg=msg)
#         self.assertAlmostEqual(t[i0], 45.0, msg=msg)
#         self.assertAlmostEqual(t[i1], 0.0, msg=msg)
#         self.assertAlmostEqual(t[i2], 0.0, msg=msg)
#         # str, int
#         self.assertAlmostEqual(counts["T", i4], 13.0, msg=msg)
    
@test
def test_pwm_mixed():
    tc = unittest.TestCase()
    counts_ = {'A': [2.0, 9.0, 0.0, 1.0, 32.0, 3.0, 46.0, 1.0, 43.0, 15.0, 2.0, 2.0], 'C': [1.0, 33.0, 45.0, 45.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 'G': [39.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 43.0], 'T': [4.0, 2.0, 0.0, 0.0, 13.0, 42.0, 0.0, 45.0, 3.0, 30.0, 0.0, 0.0]}
    m = motifs.Motif(counts=counts_)
    counts = m.counts
    pwm = counts.normalize(pseudocounts=0.25)
    pssm = pwm.log_odds()
    result = pssm.calculate(str(Seq.Seq("AcGTgTGCGtaGTGCGT")))
    tc.assertEqual(6, len(result))
    tc.assertAlmostEqual(float(result[0]), -29.18363571, places=5)
    tc.assertAlmostEqual(float(result[1]), -38.3365097, places=5)
    tc.assertAlmostEqual(float(result[2]), -29.17756271, places=5)
    tc.assertAlmostEqual(float(result[3]), -38.04542542, places=5)
    tc.assertAlmostEqual(float(result[4]), -20.3014183, places=5)
    tc.assertAlmostEqual(float(result[5]), -25.18009186, places=5)
test_pwm_mixed()

@test
def test_pwm_simple():
    tc = unittest.TestCase()
    counts = {'A': [2.0, 9.0, 0.0, 1.0, 32.0, 3.0, 46.0, 1.0, 43.0, 15.0, 2.0, 2.0], 'C': [1.0, 33.0, 45.0, 45.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 'G': [39.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 43.0], 'T': [4.0, 2.0, 0.0, 0.0, 13.0, 42.0, 0.0, 45.0, 3.0, 30.0, 0.0, 0.0]}
    m = motifs.Motif(counts=counts)
    s = str(Seq.Seq("ACGTGTGCGTAGTGCGT"))
    pwm = m.counts.normalize(pseudocounts=0.25)
    pssm = pwm.log_odds()
    result = pssm.calculate(s)
    tc.assertEqual(6, len(result))
    # The fast C-code in Bio/motifs/_pwm.c stores all results as 32-bit
    # floats; the slower Python code in Bio/motifs/__init__.py uses 64-bit
    # doubles. The C-code and Python code results will therefore not be
    # exactly equal. Test the first 5 decimal places only to avoid either
    # the C-code or the Python code to inadvertently fail this test.
    tc.assertAlmostEqual(float(result[0]), -29.18363571, places=5)
    tc.assertAlmostEqual(float(result[1]), -38.3365097, places=5)
    tc.assertAlmostEqual(float(result[2]), -29.17756271, places=5)
    tc.assertAlmostEqual(float(result[3]), -38.04542542, places=5)
    tc.assertAlmostEqual(float(result[4]), -20.3014183, places=5)
    tc.assertAlmostEqual(float(result[5]), -25.18009186, places=5)
test_pwm_simple()




# ------------------------------------------- old code ------------------------------------------- #

# class TestMotif(unittest.TestCase):
#     def __init__(self):
#         super().__init__()
    
    #     def test_format(self):
    #         m = motifs.create([Seq.Seq("ATATA")])
    #         m.name = "Foo"
    #         s1 = m.format(format_spec="pfm")

    #         expected_pfm = """  1.00   0.00   1.00   0.00   1.00
    # 0.00   0.00   0.00   0.00   0.00
    # 0.00   0.00   0.00   0.00   0.00
    # 0.00   1.00   0.00   1.00   0.00
    # """

    #         s2 = m.format(format_spec="jaspar")
    #         expected_jaspar = """>None Foo
    # A [  1.00   0.00   1.00   0.00   1.00]
    # C [  0.00   0.00   0.00   0.00   0.00]
    # G [  0.00   0.00   0.00   0.00   0.00]
    # T [  0.00   1.00   0.00   1.00   0.00]
    # """
    #         self.assertEqual(s1, expected_pfm)
    #         self.assertEqual(s2, expected_jaspar)
    #         self.assertRaises(ValueError, lambda : m.format(format_spec="foo_bar"))
        
    #     def test_format_transfac(self):
    #         m = motifs.create([Seq.Seq("ATATA")])
    #         m.name = "Foo"
    #         s = m.format(format_spec="transfac")
    #         expected_transfac = """P0      A      C      G      T
    # 01      1      0      0      0      A
    # 02      0      0      0      1      T
    # 03      1      0      0      0      A
    # 04      0      0      0      1      T
    # 05      1      0      0      0      A
    # XX
    # //
    # """
    #         self.assertEqual(s, expected_transfac)

#     def test_format_clusterbuster(self):
#         m = motifs.create([Seq.Seq("ATATA")])
#         m.name = "Foo"
#         s = m.format(format_spec="clusterbuster")

#         expected = """>Foo
# 1	0	0	0
# 0	0	0	1
# 1	0	0	0
# 0	0	0	1
# 1	0	0	0
# """

#         self.assertEqual(s, expected)


#     def test_relative_entropy_alignment(self):
#         m = motifs.create([Seq.Seq("ATATA"), Seq.Seq("ATCTA"), Seq.Seq("TTGTA")])
#         self.assertEqual(len(m.alignment), 3)
#         self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
#         self.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})

#         self.assertTrue(
#             np.allclose(
#                 m.relative_entropy,
#                 np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
#             )
#         )

#         m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
#         self.assertTrue(
#             np.allclose(
#                 m.relative_entropy,
#                 np.array(
#                     [
#                         0.8186697601117167,
#                         1.7369655941662063,
#                         0.5419780939258206,
#                         1.7369655941662063,
#                         1.7369655941662063,
#                     ]
#                 ),
#             )
#         )

#         m.background = None
#         self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
#         pseudocounts = math.sqrt(len(m.alignment))
#         m.pseudocounts = {
#             letter: m.background[letter] * pseudocounts for letter in "ACGT"
#         }
#         self.assertTrue(
#             np.allclose(
#                 m.relative_entropy,
#                 np.array(
#                     [
#                         0.3532586861097656,
#                         0.7170228827697498,
#                         0.11859369972847714,
#                         0.7170228827697498,
#                         0.7170228827697499,
#                     ]
#                 ),
#             )
#         )

#         m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
#         self.assertTrue(
#             np.allclose(
#                 m.relative_entropy,
#                 np.array(
#                     [
#                         0.19727984803857979,
#                         0.561044044698564,
#                         0.20984910512125132,
#                         0.561044044698564,
#                         0.5610440446985638,
#                     ]
#                 ),
#             )
#         )
    
#     def test_relative_entropy_counts(self):
#         m = motifs.Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
#         self.assertTrue(
#             np.allclose(
#                 m.relative_entropy,
#                 np.array(
#                     [
#                         0.08170417,
#                         0.02904941,
#                         0.01477186,
#                     ]
#                 ),
#             )
#         )
    
#     def test_pwm(self):
#         m = motifs.create([Seq.Seq(("ATATA"))])
#         expected = """        0      1      2      3      4
# A:   1.00   0.00   1.00   0.00   1.00
# C:   0.00   0.00   0.00   0.00   0.00
# G:   0.00   0.00   0.00   0.00   0.00
# T:   0.00   1.00   0.00   1.00   0.00
# """
#         self.assertEqual(expected, m.pwm.__str__())

#     def test_pssm(self):
#         m = motifs.create([Seq.Seq(("ATATA"))])
#         expected="""        0      1      2      3      4
# A:   2.00   -inf   2.00   -inf   2.00
# C:   -inf   -inf   -inf   -inf   -inf
# G:   -inf   -inf   -inf   -inf   -inf
# T:   -inf   2.00   -inf   2.00   -inf
# """
#         self.assertEqual(expected, m.pssm.__str__())

#     def test_str(self):
#         m = motifs.create([Seq.Seq(("ATATA"))])
#         self.assertEqual("ATATA", m.__str__())

#         m.mask = "* * *"
#         self.assertEqual('ATATA* * *\n', m.__str__(masked=True))

#     def test_mask(self):
#         m = motifs.create([Seq.Seq(("ATATA"))])
#         self.assertEqual([1] * m.length, m.mask)

#         m.mask = "* * *"
#         self.assertEqual([1, 0, 1, 0, 1], m.mask)

#         m.mask = [2, 0, 3, 0, 1]
#         self.assertEqual([1, 0, 1, 0, 1], m.mask)

#         def exception(): m.mask = "abcab"
#         self.assertRaises(ValueError, exception)

#         def exception(): m.mask = [1,2]
#         self.assertRaises(ValueError, exception)

#     def test_reverse_complement(self):
#         background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
#         pseudocounts = 0.5
#         m = motifs.create([Seq.Seq(("ATATA"))])
#         m.background = background
#         m.pseudocounts = pseudocounts

#         received_forward = m.format(format_spec="transfac")
#         expected_forward = """P0      A      C      G      T
# 01      1      0      0      0      A
# 02      0      0      0      1      T
# 03      1      0      0      0      A
# 04      0      0      0      1      T
# 05      1      0      0      0      A
# XX
# //
# """
#         self.assertEqual(received_forward, expected_forward)

#         expected_forward_pwm = """        0      1      2      3      4
# A:   0.50   0.17   0.50   0.17   0.50
# C:   0.17   0.17   0.17   0.17   0.17
# G:   0.17   0.17   0.17   0.17   0.17
# T:   0.17   0.50   0.17   0.50   0.17
# """
#         self.assertEqual(str(m.pwm), expected_forward_pwm)
        
#         rc = m.reverse_complement()

#         received_reverse = rc.format(format_spec="transfac")
#         expected_reverse = """P0      A      C      G      T
# 01      0      0      0      1      T
# 02      1      0      0      0      A
# 03      0      0      0      1      T
# 04      1      0      0      0      A
# 05      0      0      0      1      T
# XX
# //
# """
#         self.assertEqual(received_reverse, expected_reverse)

#         expected_reverse_pwm = """        0      1      2      3      4
# A:   0.17   0.50   0.17   0.50   0.17
# C:   0.17   0.17   0.17   0.17   0.17
# G:   0.17   0.17   0.17   0.17   0.17
# T:   0.50   0.17   0.50   0.17   0.50
# """
#         self.assertEqual(expected_reverse_pwm, str(rc.pwm))


#         background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
#         pseudocounts = 0.5
#         m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
#         m_rna.background = background_rna
#         m_rna.pseudocounts = pseudocounts
#         expected_forward_rna_counts = """        0      1      2      3      4
# A:   1.00   0.00   1.00   0.00   1.00
# C:   0.00   0.00   0.00   0.00   0.00
# G:   0.00   0.00   0.00   0.00   0.00
# U:   0.00   1.00   0.00   1.00   0.00
# """
#         self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        
#         expected_forward_rna_pwm = """        0      1      2      3      4
# A:   0.50   0.17   0.50   0.17   0.50
# C:   0.17   0.17   0.17   0.17   0.17
# G:   0.17   0.17   0.17   0.17   0.17
# U:   0.17   0.50   0.17   0.50   0.17
# """
#         self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
#         expected_reverse_rna_counts = """        0      1      2      3      4
# A:   0.00   1.00   0.00   1.00   0.00
# C:   0.00   0.00   0.00   0.00   0.00
# G:   0.00   0.00   0.00   0.00   0.00
# U:   1.00   0.00   1.00   0.00   1.00
# """
#         self.assertEqual(
#             str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
#         )
#         expected_reverse_rna_pwm = """        0      1      2      3      4
# A:   0.17   0.50   0.17   0.50   0.17
# C:   0.17   0.17   0.17   0.17   0.17
# G:   0.17   0.17   0.17   0.17   0.17
# U:   0.50   0.17   0.50   0.17   0.50
# """
#         self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

#         m = motifs.create([Seq.Seq("ATATA")])
#         counts = m.counts
#         m = motifs.Motif(counts=counts)
#         m.background = background
#         m.pseudocounts = pseudocounts
#         received_forward = m.format(format_spec="transfac")
#         self.assertEqual(received_forward, expected_forward)
#         self.assertEqual(str(m.pwm), expected_forward_pwm)
#         m = m.reverse_complement()
#         received_reverse = m.format(format_spec="transfac")
#         self.assertEqual(received_reverse, expected_reverse)
#         self.assertEqual(str(m.pwm), expected_reverse_pwm)
#         # Same, but for RNA count matrix
#         m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
#         counts = m_rna.counts
#         m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
#         m_rna.background = background_rna
#         m_rna.pseudocounts = pseudocounts
#         self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
#         self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
#         self.assertEqual(
#             str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
#         )
#         self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

#     def test_consensus(self):
#         m = motifs.create([
#             Seq.Seq("ATATA"),
#             Seq.Seq("ATCTA"),
#             Seq.Seq("TTGTA"),
#             Seq.Seq("ATATA"),
#         ])

#         expected_consensus = "ATATA"
#         self.assertEqual(str(m.consensus), expected_consensus)

#     def test_anticonsensus(self):
#         m = motifs.create([
#             Seq.Seq("ATCGA"),
#             Seq.Seq("ATCGA"),
#             Seq.Seq("GGGTG"),
#             Seq.Seq("GGGTG"),
#             Seq.Seq("CCACC"),
#             Seq.Seq("CCACC"),
#             Seq.Seq("TATAT")
#         ])

#         expected_anticonsensus = "TATAT"
#         self.assertEqual(str(m.anticonsensus), expected_anticonsensus)

#     def test_degenerate_consensus(self):
#         m = motifs.create([
#             Seq.Seq("ATATA"),
#             Seq.Seq("ATCTA"),
#             Seq.Seq("TTGTA"),
#             Seq.Seq("ATGTA")
#         ])
#         # Position-wise breakdown:
#         # 1: A (3), T (1) → A
#         # 2: T (4) → T
#         # 3: A (1), C (1), G (2) → V (A/C/G)
#         # 4: T (4) → T
#         # 5: A (4) → A

#         expected_degenerate_consensus = "ATVTA"
#         self.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)

#     def test_degenerate_consensus_with_ties(self):
#         m = motifs.create([
#             Seq.Seq("A"),
#             Seq.Seq("C"),
#             Seq.Seq("G"),
#             Seq.Seq("T"),
#         ])
#         expected_degenerate_consensus = "N"  # all bases equally represented
#         self.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)

#     def test_degenerate_consensus_rna(self):
#         m_rna = motifs.create([
#             Seq.Seq("AUAUA"),
#             Seq.Seq("AUCUA"),
#             Seq.Seq("UUGUA"),
#             Seq.Seq("AUGUA")
#         ], alphabet="ACGU")

#         expected_degenerate_consensus = "AUVUA"
#         self.assertEqual(str(m_rna.degenerate_consensus), expected_degenerate_consensus)

#     def test_getitem(self):
#         m = motifs.create(["AACGCCA", "ACCGCCC", "AACTCCG"])
#         expected="""AACGCC\nACCGCC\nAACTCC"""
#         self.assertEqual(str(m[:-1]), expected)

#         expected=""
#         self.assertEqual(str(m[0:0]), expected)

#         expected="A\nA\nA"
#         self.assertEqual(str(m[0:1]), expected)
        
#         # EXCEPTIONS
#         self.assertRaises(TypeError, lambda: m[0])
#         self.assertRaises(TypeError, lambda: m['A'])


# tests = TestMotif()

# # shitty test runner
# def run_test(name, func):
#     try:
#         func()
#         print(f"{name} passed.")
#     except:
#         print(f"{name} failed.")

# run_test("test_format", tests.test_format)
# run_test("test_relative_entropy_alignment", tests.test_relative_entropy_alignment)
# run_test("test_relative_entropy_counts", tests.test_relative_entropy_counts)
# run_test("test_pwm", tests.test_pwm)
# run_test("test_pssm", tests.test_pssm)
# run_test("test_str", tests.test_str)
# run_test("test_mask", tests.test_mask)
# run_test("test_reverse_complement", tests.test_reverse_complement)
# run_test("test_consensus", tests.test_consensus)
# run_test("test_anticonsensus", tests.test_anticonsensus)
# run_test("test_degenerate_consensus", tests.test_degenerate_consensus)
# run_test("test_degenerate_consensus_with_ties", tests.test_degenerate_consensus_with_ties)
# run_test("test_degenerate_consensus_rna", tests.test_degenerate_consensus_rna)
# run_test("test_getitem", tests.test_getitem)
# run_test("test_format_transfac", tests.test_format_transfac)
# run_test("test_format_clusterbuster", tests.test_format_clusterbuster)
# run_test("test_minimal_parser_1", tests.test_minimal_parser_1)
# run_test("test_minimal_parser_2", tests.test_minimal_parser_2)
# run_test("test_minimal_parser_rna", tests.test_minimal_parser_rna)
# run_test("test_pwm_getitem", tests.test_pwm_getitem)
# run_test("test_pwm_simple", tests.test_pwm_simple)
# run_test("test_pwm_mixed", tests.test_pwm_mixed)
