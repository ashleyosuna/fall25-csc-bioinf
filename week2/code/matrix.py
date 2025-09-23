from python import Bio.Seq as Seq
# import Bio.Seq as Seq

class GenericPositionMatrix(dict[str, List[float]]):
    alphabet: str
    length: int

    def __init__(self, alphabet: str, values: dict[str, List[int]]):
        super().__init__()
        length = None
        self.alphabet = alphabet
        for letter in alphabet:
            if length is None:
                self.length = len(values[letter])
            elif self.length != len(values[letter]):
                raise Exception("data has inconsistent lengths")
            self[letter] = [float(val) for val in values[letter]]
        print('->', super())
    
    def __init__(self, alphabet: str, values: dict[str, List[float]]):
        super().__init__()
        length = None
        self.alphabet = alphabet
        for letter in alphabet:
            if length is None:
                self.length = len(values[letter])
            elif self.length != len(values[letter]):
                raise Exception("data has inconsistent lengths")
            self[letter] = values[letter]

    def __str__(self):
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = [f"{i:6}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = [f"{val:6.2f}" for val in self[letter]]
            line = f"{letter}: " + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text
    
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
    
#     # TODO: degenerate consensus

#     # TODO: calculate consensus
    
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
    
class FrequencyPositionMatrix(GenericPositionMatrix):
    # alphabet: str
    # length: int

    def __init__(self, alphabet: str, values: dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=values)
        self.alphabet = super().__getalphabet__()
        self.length = super().__getlength__()
        print('=>', super())
    
    # def normalize(self):
    #     return

    def normalize(self, pseudocounts: int = 0):
        print('here', self.alphabet, self.length)
        counts: dict[str, List[float]] = {}
        pseudocounts = float(pseudocounts)
        for letter in self.alphabet:
            counts[letter] = [pseudocounts] * self.length
        print(counts)
        for i in range(self.length):
            for letter in self.alphabet:
                counts[letter][i] += self[letter][i]
        return counts



positionMatrix = GenericPositionMatrix(alphabet="ACGT", values={"A": [1, 2, 3], "C": [2, 3, 4], "G": [1, 2, 3], "T": [2, 4, 5]})
print(positionMatrix)
# print(positionMatrix.consensus)
# print(positionMatrix.anticonsensus)
# print(positionMatrix.gc_content)

rc = positionMatrix.reverse_complement()
print(rc)

# TODO: fix whatever the hell is wrong here
freqMatrix = FrequencyPositionMatrix(alphabet="ACGT", values={"A": [1, 2, 3], "C": [2, 3, 4], "G": [1, 2, 3], "T": [2, 4, 5]})
# print(freqMatrix)
print(freqMatrix.normalize(pseudocounts=1))