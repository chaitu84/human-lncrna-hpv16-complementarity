from Bio import SeqIO


def read_fasta(filename):  # filename is a placeholder
    sequences = []

    for record in SeqIO.parse(filename, "fasta"):
        sequences.append(record)

    return sequences


if __name__ == "__main__":
    hpv_sequences = read_fasta(
        "data/HPV16_E6/ncbi_dataset/data/gene.fna"
    )

    for sequence in hpv_sequences:
        print("Name:", sequence.id)
        print("Sequence:", sequence.seq)
        print()
