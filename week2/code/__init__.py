from python import Bio.Align as Align
import numpy as np
import matrix
from typing import Optional, Dict, List
from python import Bio.Seq as Seq

def create(instances, alphabet="ACGT"):
    alignment = Align.Alignment(instances)
    return Motif(alphabet=alphabet, alignment=alignment)

class Motif:
    name: str
    counts: Optional[matrix.FrequencyPositionMatrix]
    length: Optional[int]
    alignment: Optional[pyobj]
    alphabet: str
    pseudocounts: Optional[int]
    # background: Optional[]
    # mask: Optional[]

    def __init__(self, alphabet="ACGT", alignment: Optional[pyobj]=None, counts: Optional[Dict[str, List[float]]]=None):
        self.name = ""

        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        elif counts is not None:
            self.alignment = None
            self.counts = matrix.FrequencyPositionMatrix(alphabet=alphabet, values=counts)
            self.length = self.counts.length
        elif alignment is not None:
            length = alignment.length
            frequencies: Dict[str, List[float]] = {}
            for letter in alphabet:
                if letter not in alignment.frequencies:
                    frequencies[letter] = [0.0 for _ in range(length)]
                else:
                    frequencies[letter] = [float(alignment.frequencies[letter][i]) for i in range(length)]
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
            self.length = length
        else:
            self.counts = None
            self.alignment = None
            self.length = None
        self.alphabet = alphabet
        self.pseudocounts = None
        # self.background = None
        # self.mask = None
    
    def __len__(self):
        return 0 if self.length is None else self.length
    
# motif = Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
# print('-->', len(motif))