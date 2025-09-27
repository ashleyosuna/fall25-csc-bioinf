from typing import Dict, List, Optional
from python import Bio.Seq as Seq
import math
import numpy as np
import _pwm

class GenericPositionMatrix:
    alphabet: str
    length: int
    data: Dict[str, List[float]]

    def __init__(self, alphabet: str, values: Dict[str, List[int]]):
        self.data: Dict[str, List[float]] = Dict[str, List[float]]()
        self.alphabet = alphabet

        length = None
        for letter in alphabet:
            vals: List[float] = [float(v) for v in values[letter]]
            if length is None:
                self.length = len(vals)
            elif length != len(vals):
                raise Exception("data has inconsistent lengths")
            self.data[letter] = vals
    
    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        self.data: Dict[str, List[float]] = Dict[str, List[float]]()
        self.alphabet = alphabet

        length = None
        for letter in alphabet:
            vals: List[float] = [float(v) for v in values[letter]]  # normalize here
            if length is None:
                self.length = len(vals)
            elif length != len(vals):
                raise Exception("data has inconsistent lengths")
            self.data[letter] = vals

    def __str__(self):
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = [f"{i:6d}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = [f"{val:6.2f}" for val in self.data[letter]]
            line = f"{letter}: " + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text
    
#     # TODO: finish
    def __getitem__(self, key) -> List[float] | Dict[str, List[float]]:
        if isinstance(key, str):
            return self.data[key]

        elif isinstance(key, int):
            return self.data[self.alphabet[key]]

        # elif isinstance(key, slice):
        #     start, stop, stride = key.indices(len(self.alphabet))
        #     letters = [self.alphabet[i] for i in range(start, stop, stride)]
        #     d: Dict[str, List[float]] = Dict[str, List[float]]()
        #     for letter in letters:
        #         d[letter] = Dict[str, List[float]].__getitem__(self, letter)
        #     return d

        else:
            raise KeyError(f"Unsupported key type: {key}")

    @property
    def consensus(self):
        sequence = ""
        for i in range(self.length):
            max_count = float('-inf')
            for letter in self.alphabet:
                count = self[letter][i]
                if count > max_count:
                    max_count = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq.Seq(sequence)
    
    @property
    def anticonsensus(self):
        sequence = ""
        for i in range(self.length):
            min_count = float('inf')
            for letter in self.alphabet:
                count = self[letter][i]
                if count < min_count:
                    min_count = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq.Seq(sequence)

    @property
    def degenerate_consensus(self):
        """Return the degenerate consensus sequence."""
        degenerate_nucleotide: dict[str, str] = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "U": "U",
            "AC": "M",
            "AG": "R",
            "AT": "W",
            "AU": "W",
            "CG": "S",
            "CT": "Y",
            "CU": "Y",
            "GT": "K",
            "GU": "K",
            "ACG": "V",
            "ACT": "H",
            "ACU": "H",
            "AGT": "D",
            "AGU": "D",
            "CGT": "B",
            "CGU": "B",
            "ACGT": "N",
            "ACGU": "N",
        }
        sequence = ""
        for i in range(self.length):

            def get(nucleotide):
                return self[nucleotide][i]  # noqa: B023

            nucleotides = sorted(self.data, key=get, reverse=True)

            counts = [self[c][i] for c in nucleotides]
            # Follow the Cavener rules:
            if counts[0] > sum(counts[1:]) and counts[0] > 2 * counts[1]:
                key = nucleotides[0]
            elif 4 * sum(counts[:2]) > 3 * sum(counts):
                key = "".join(sorted(nucleotides[:2]))
            elif counts[3] == 0:
                key = "".join(sorted(nucleotides[:3]))
            else:
                key = "ACGT"
            
            # nucleotide = degenerate_nucleotide.get(key, key)
            if key in degenerate_nucleotide:
                nucleotide = degenerate_nucleotide[key]
            else:
                nucleotide = key

            sequence += nucleotide
        return Seq.Seq(sequence)

    def calculate_consensus(self, substitution_matrix=None, plurality=None, identity=0, setcase=None):
        alphabet = self.alphabet
        if set(alphabet).union(set("ACGTUN-")) == set("ACGTUN-"):
            undefined = "N"
        else:
            undefined = "X"
        if substitution_matrix is None:
            if plurality is not None:
                raise ValueError(
                    "plurality must be None if substitution_matrix is None"
                )
            sequence = ""
            for i in range(self.length):
                maximum: float = 0.0
                total: float = 0.0
                for letter in alphabet:
                    count = self[letter][i]
                    total += count
                    if count > maximum:
                        maximum = count
                        consensus_letter = letter
                if maximum < identity * total:
                    consensus_letter = undefined
                else:
                    if setcase is None:
                        setcase_threshold = total / 2
                    else:
                        setcase_threshold = setcase * total
                    if maximum <= setcase_threshold:
                        consensus_letter = consensus_letter.lower()
                sequence += consensus_letter
        else:
            raise NotImplementedError(
                "calculate_consensus currently only supports substitution_matrix=None"
            )
        return sequence

    @property
    def gc_content(self):
        """Compute the fraction GC content."""
        alphabet = self.alphabet
        gc_total = 0.0
        total = 0.0
        for i in range(self.length):
            for letter in alphabet:
                if letter in "CG":
                    gc_total += self[letter][i]
                total += self[letter][i]
        return gc_total / total
    
    def reverse_complement(self):
        values = {}
        if self.alphabet == "ACGU":
            values["A"] = self["U"][::-1]
            values["U"] = self["A"][::-1]
        else:
            values["A"] = self["T"][::-1]
            values["T"] = self["A"][::-1]
        values["G"] = self["C"][::-1]
        values["C"] = self["G"][::-1]
        alphabet = self.alphabet
        return self.__class__(alphabet, values)

    def __getalphabet__(self):
        return self.alphabet
    
    def __getlength__(self):
        return self.length
    
    def __getdata__(self):
        return self.data
    
class FrequencyPositionMatrix(GenericPositionMatrix):
    alphabet: str
    length: int

    def __init__(self, alphabet: str, values: Dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=values)
        self.length = super().__getlength__()
        self.alphabet = super().__getalphabet__()
    
    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        super().__init__(alphabet=alphabet, values=values)
        self.length = super().__getlength__()
        self.alphabet = super().__getalphabet__()

    def normalize(self, pseudocounts: int = 0):
        counts: Dict[str, List[float]] = {}
        pseudocounts = float(pseudocounts)
        for letter in self.alphabet:
            counts[letter] = [pseudocounts] * self.length
        for i in range(self.length):
            for letter in self.alphabet:
                counts[letter][i] += self[letter][i]
        return PositionWeightMatrix(alphabet=self.alphabet, counts=counts)
    
class PositionWeightMatrix(GenericPositionMatrix):
    length: int
    alphabet: str

    def __init__(self, alphabet: str, counts: Dict[str, List[float]]):
        super().__init__(alphabet=alphabet, values=counts)
        self.length = super().__getlength__()
        self.alphabet = alphabet

        for i in range(self.length):
            total = sum(self[letter][i] for letter in alphabet)
            for letter in alphabet:
                self[letter][i] /= total
    
    def __init__(self, alphabet: str, counts: Dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=counts)
        self.length = super().__getlength__()
        self.alphabet = alphabet

        for i in range(self.length):
            total = sum(self[letter][i] for letter in alphabet)
            for letter in alphabet:
                self[letter][i] /= total
    
    def log_odds(self, background: Optional[Dict[str, float]]=None):
        values: Dict[str, List[float]] = {}
        alphabet = self.alphabet

        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())

        for letter in alphabet:
            background[letter] /= total
            values[letter] = []
        for i in range(self.length):
            for letter in alphabet:
                b = background[letter]

                if b > 0:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.log(p / b, 2)
                    else:
                        logodds = -math.inf
                else:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.inf
                    else:
                        logodds = math.nan
                values[letter].append(logodds)
        pssm = PositionSpecificScoringMatrix(alphabet=alphabet, values=values)
        return pssm

class PositionSpecificScoringMatrix(GenericPositionMatrix):
    alphabet: str
    length: int

    def __init__(self, alphabet: str, values: Dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=values)
        self.alphabet = alphabet
        self.length = super().__getlength__()

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        super().__init__(alphabet=alphabet, values=values)
        self.alphabet = alphabet
        self.length = super().__getlength__()    
    
    def calculate(self, sequence: str):
        if sorted(self.alphabet) != ['A', 'C', 'G', 'T']:
            raise ValueError(f"PSSM has wrong alphabet: {self.alphabet} - Use only with DNA motifs")

        n = len(sequence)
        m = self.length

        scores = np.empty(n - m + 1, np.float32)
        logodds = np.array(
            [[self[letter][i] for letter in "ACGT"] for i in range(m)], float
        )
        _pwm.calculate(sequence, logodds, scores)
        return scores

#   # TODO: this whole thing
    # def search(self, sequence, threshold=0.0, both=True, chunksize=10**6):

    @property
    def max(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        score = 0.0
        letters = self.alphabet
        for position in range(self.length):
            score += max(self[letter][position] for letter in letters)
        return score

    @property
    def min(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        score = 0.0
        letters = self.alphabet
        for position in range(self.length):
            score += min(self[letter][position] for letter in letters)
        return score

    @property
    def gc_content(self):
        """Compute the GC-ratio."""
        return super().gc_content

values: Dict[str, List[int]] = {"A": [1, 2, 3], "C": [2, 3, 4], "G": [1, 2, 3], "T": [2, 4, 5]}
# positionMatrix = GenericPositionMatrix(alphabet="ACGT", values=values)
# print(positionMatrix)
# print(positionMatrix.consensus)
# print(positionMatrix.anticonsensus)
# print(positionMatrix.gc_content)
# TODO: degenerate_consensus causing issues?
# print(positionMatrix.degenerate_consensus)
# print(positionMatrix.calculate_consensus())

# rc = positionMatrix.reverse_complement()
# print(rc)

# TODO: fix whatever the hell is wrong here
# freqMatrix = FrequencyPositionMatrix(alphabet="ACGT", values={"A": [1, 2, 3], "C": [2, 3, 4], "G": [1, 2, 3], "T": [2, 4, 5]})
# print(freqMatrix)
# print(freqMatrix.normalize(pseudocounts=1))

# positionMatrix2 = PositionWeightMatrix(alphabet="ACGT", counts=values)
# print(positionMatrix2)
# print('logodds', positionMatrix2.log_odds())

# positionSpecific = PositionSpecificScoringMatrix(alphabet="ACGT", values=values)
# print(positionSpecific)
# print(positionSpecific.calculate("CGTA"))
# print(positionSpecific.max)
# print(positionSpecific.min)
# print(positionSpecific.gc_content)
