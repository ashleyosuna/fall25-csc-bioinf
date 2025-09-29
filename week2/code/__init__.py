from python import Bio.Align as Align
import numpy as np
import matrix
from typing import Optional, Dict, List
from python import Bio.Seq as Seq
import utilities

def create(instances, alphabet="ACGT"):
    alignment = Align.Alignment(instances)
    return Motif(alphabet=alphabet, alignment=alignment)

class Motif:
    name: str
    counts: Optional[matrix.FrequencyPositionMatrix]
    length: int
    alignment: Optional[pyobj]
    alphabet: str
    _pseudocounts: Optional[Dict[str, float]]
    # background: Optional[Dict[str, float]]
    _background: Optional[Dict[str, float]]

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
                if letter not in list(alignment.frequencies.keys()):
                    frequencies[letter] = [0.0 for _ in range(length)]
                else:
                    frequencies[letter] = [float(alignment.frequencies[letter][i]) for i in range(length)]
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
            self.length = length
        else:
            self.counts = None
            self.alignment = None
            self.length = 0
        self.alphabet = alphabet
        self._pseudocounts = dict.fromkeys(self.alphabet, 0.0)
        self._background = None
        self.__set_background(value=None)
        # self.mask = None
    
    def __len__(self):
        return 0 if self.length is None else self.length
    
    def __format__(self, format_spec: Optional[str] = "default", **kwargs):
        if format_spec in ["pfm", "jaspar"]:
            motifs = [self]
            return utilities.jaspar_write(motifs, format_spec)
        # elif format_spec == "transfac":
        #     # from Bio.motifs import transfac
        #     from python import Bio.motifs.transfac as transfac

        #     motifs = [self]
        #     return transfac.write(motifs)
        # elif format_spec == "clusterbuster":
        #     # from Bio.motifs import clusterbuster
        #     from python import Bio.motifs.clustebuster as clusterbuster

        #     motifs = [self]
        #     return clusterbuster.write(motifs, **kwargs)
        elif format_spec == "default":
            # Follow python convention and default to using __str__
            return str(self)
        else:
            raise ValueError(f"Unknown format type {format_spec}")
    
    def format(self, format_spec):
        return self.__format__(format_spec)
    
    @property
    def pseudocounts(self):
        return self._pseudocounts
    
    @pseudocounts.setter
    def pseudocounts(self, value):
        self._pseudocounts = {}
        if isinstance(value, Dict[str, float]):
            self._pseudocounts = {letter: value[letter] for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = dict.fromkeys(self.alphabet, value)
    
    @property
    def background(self):
        return self._background
    
    def __set_background(self, value = None):
        if isinstance(value, Dict[str, float]): self._background = {letter: value[letter] for letter in self.alphabet}
        elif value is None: self._background = dict.fromkeys(self.alphabet, 1.0)
        elif isinstance(value, float) or isinstance(value, int):
            if not self._has_dna_alphabet() and not self._has_rna_alphabet():
                raise ValueError(
                    "Setting the background to a single value only works for DNA and RNA"
                    "motifs (in which case the value is interpreted as the GC content)"
                )
            value = float(value)
            T_or_U = "T" if self._has_dna_alphabet() else "U"
            d = {}
            d['A'] = (1.0 - value) / 2.0
            d["C"] = value / 2.0
            d["G"] = value / 2.0
            d[T_or_U] = (1.0 - value) / 2.0
            self._background = d
        
        # normalize
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total
    
    @background.setter
    def background(self, value = None):
        self.__set_background(value)

    def _has_dna_alphabet(self):
        return sorted(self.alphabet) == ["A", "C", "G", "T"]
    
    def _has_rna_alphabet(self):
        return sorted(self.alphabet) == ["A", "C", "G", "U"]
    
    @property
    def relative_entropy(self):
        background = self.background
        pseudocounts = self.pseudocounts
        alphabet = self.alphabet
        counts = self.counts
        length = self.length
        values = np.zeros(length)
        # if self.alignment is None:
        #     total = []
        #     # total = np.array(
        #     #     [
        #     #         sum(counts[c][i] + pseudocounts[c] for c in alphabet)
        #     #         for i in range(length)
        #     #     ]
        #     # )
        #     for letter, frequencies in counts.items():
        #         frequencies = np.array(frequencies) + pseudocounts[letter]
        #         mask = frequencies > 0
        #         frequencies = frequencies[mask] / total[mask]
        #         values[mask] += frequencies * np.log2(frequencies / background[letter])
        if self.alignment is not None:
            total = np.zeros(length)
            for letter in alphabet:
                frequencies = []
                for i in range(length):
                    frequencies.append(counts[letter][i])
                total += np.array(frequencies) + pseudocounts[letter]
            for letter in alphabet:
                frequencies = []
                for i in range(length):
                    frequencies.append(counts[letter][i])
                frequencies = np.array(frequencies) + pseudocounts[letter]
                mask = frequencies > 0
                frequencies = frequencies[mask] / total[mask]
                values[mask] += frequencies * np.log2(frequencies / background[letter])
        return values
    
# motif = Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
# print('-->', len(motif))