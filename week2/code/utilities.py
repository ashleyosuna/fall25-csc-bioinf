from matrix import FrequencyPositionMatrix

def jaspar_write(motifs, format):
    """Return the representation of motifs in "pfm" or "jaspar" format."""
    letters = "ACGT"
    lines: List[str] = []
    if format == "pfm":
        motif = motifs[0]
        counts: FrequencyPositionMatrix = motif.counts
        for letter in letters:
            terms = []
            for i in range(motif.length):
                terms.append(f"{counts[letter][i]:6.2f}")
            line = f"{' '.join(terms)}\n"
            lines.append(line)
    elif format == "jaspar":
        for m in motifs:
            counts: FrequencyPositionMatrix = m.counts
            # try:
            #     matrix_id = m.matrix_id
            # except AttributeError:
            #     matrix_id = None
            line = f">None {m.name}\n"
            lines.append(line)
            for letter in letters:
                terms = []
                for i in range(m.length):
                    terms.append(f"{counts[letter][i]:6.2f}")
                line = f"{letter} [{' '.join(terms)}]\n"
                lines.append(line)
    else:
        raise ValueError(f"Unknown JASPAR format {format}")

    # Finished; glue the lines together
    text = "".join(lines)
    return text

# def transfac_write(motifs):
#     blocks = []
#     try:
#         version = motifs.version
#     except AttributeError:
#         pass
#     else:
#         if version is not None:
#             block = (
#                 """\
# VV  %s
# XX
# //
# """
#                 % version
#             )
#             blocks.append(block)
#     multiple_value_keys = Motif.multiple_value_keys
#     sections = (
#         ("AC", "AS"),  # Accession
#         ("ID",),  # ID
#         ("DT", "CO"),  # Date, copyright
#         ("NA",),  # Name
#         ("DE",),  # Short factor description
#         ("TY",),  # Type
#         ("OS", "OC"),  # Organism
#         ("HP", "HC"),  # Superfamilies, subfamilies
#         ("BF",),  # Binding factors
#         ("P0",),  # Frequency matrix
#         ("BA",),  # Statistical basis
#         ("BS",),  # Factor binding sites
#         ("CC",),  # Comments
#         ("DR",),  # External databases
#         ("OV", "PV"),  # Versions
#     )
#     for motif in motifs:
#         lines = []
#         for section in sections:
#             blank = False
#             for key in section:
#                 if key == "P0":
#                     # Frequency matrix
#                     length = motif.length
#                     if length == 0:
#                         continue
#                     sequence = motif.degenerate_consensus
#                     letters = sorted(motif.alphabet)
#                     line = "      ".join(["P0"] + letters)

#                     lines.append(line)
#                     for i in range(length):
#                         line = (
#                             " ".join(["%02.d"] + ["%6.20g" for _ in letters])
#                             + "      %s"
#                         )
#                         line = line % tuple(
#                             [i + 1]
#                             + [motif.counts[_][i] for _ in letters]
#                             + [sequence[i]]
#                         )
#                         lines.append(line)
#                     blank = True
#                 else:
#                     try:
#                         value = motif.get(key)
#                     except AttributeError:
#                         value = None
#                     if value is not None:
#                         if key in multiple_value_keys:
#                             for v in value:
#                                 line = f"{key}  {v}"
#                                 lines.append(line)
#                         else:
#                             line = f"{key}  {value}"
#                             lines.append(line)
#                         blank = True
#                 if key == "PV":
#                     # References
#                     try:
#                         references = motif.references
#                     except AttributeError:
#                         pass
#                     else:
#                         keys = ("RN", "RX", "RA", "RT", "RL")
#                         for reference in references:
#                             for key in keys:
#                                 value = reference.get(key)
#                                 if value is None:
#                                     continue
#                                 line = f"{key}  {value}"
#                                 lines.append(line)
#                                 blank = True
#             if blank:
#                 line = "XX"
#                 lines.append(line)
#         # Finished this motif; glue the lines together
#         line = "//"
#         lines.append(line)
#         block = "\n".join(lines) + "\n"
#         blocks.append(block)
#     # Finished all motifs; glue the blocks together
#     text = "".join(blocks)
#     return text