import math
# import unittest
from python import Bio.Seq as Seq
# import __init__ as motifs
import unittest
import __init__ as motifs

# if __codon__:
#    import __init__ as motifs  # depending on your naming
# else:
#    from Bio import motifs  # python path; you do not have to include Python files into the repository

class TestMotif(unittest.TestCase):
    passed: int
    total: int

    def __init__(self):
        super().__init__()
        self.passed = 0
        self.total = 0
    
    def testConstruction(self):
        res = self.assertEqual(True, True)

        m = motifs.create([Seq.Seq("ATATA")])
        m.name = "Foo"

        if not res:
            self.passed += 1
        self.total += 1
    
    def __str__(self):
        return f"{self.passed} tests passed out of {self.total} tests"
    
    # def test_format(self):
    #     m = motifs.create([Seq("ATATA")])
        # m = matrix.FrequencyPositionMatrix(alphabet="ACGT", values={'A': [1.0, 2.0, 3.0]})

# runner = unittest.TextTestRunner(verbosity=2)
# unittest.main(testRunner=runner)
tests = TestMotif()
tests.testConstruction()
print(tests)