# from python import Bio.Alignment as Alignment
import numpy as np
import matrix
from typing import Optional, Dict, List

class Motif:
    name: str
    counts: Optional[matrix.FrequencyPositionMatrix]
    length: Optional[int]
    # alignment: Optional[]
    alphabet: str
    pseudocounts: Optional[int]
    # background: Optional[]
    # mask: Optional[]

    def __init__(self, alphabet="ACGT", alignment=None, counts: Optional[Dict[str, List[float]]]=None):
        self.name = ""

        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        elif counts is not None:
            # self.alignment = None
            self.counts = matrix.FrequencyPositionMatrix(alphabet=alphabet, values=counts)
            self.length = self.counts.length
        # elif alignment is not None:
        #     length = alignment.length
        #     frequencies = alignment.frequencies
        #     for letter in alphabet:
        #         if letter not in frequencies:
        #             frequencies[letter] = np.zeros(length, int)
        #     self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
        #     self.alignment = alignment
        #     self.length = length
        else:
            self.counts = None
            # self.alignment = None
            self.length = None
        self.alphabet = alphabet
        self.pseudocounts = None
        # self.background = None
        # self.mask = None
    
    def __len__(self):
        return 0 if self.length is None else self.length
    
motif = Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
print('-->', len(motif))