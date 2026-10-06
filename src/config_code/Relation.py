from Bio import SeqIO
import pandas as pd

def find_matches(hpv_sequence, human_sequence, match_length=10):
    matches = []

    hpv_sequence = str(hpv_sequence).upper()
    human_sequence = str(human_sequence).upper()

    for i in range(len(hpv_sequence) - match_length + 1):
        hpv_piece = hpv_sequence[i:i + match_length]

        complement = str(Seq(hpv_piece).reverse_complement())

        if complement in human_sequence:
            matches.append({
                "hpv_position": i,
                "hpv_sequence": hpv_piece,
                "complement": complement
            })

    return matches

