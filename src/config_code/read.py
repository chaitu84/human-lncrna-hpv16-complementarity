
import gzip
from pathlib import Path
from Bio import SeqIO
import pandas as pd
#gencode.v50.lncRNA_transcripts.fa.gz
#gencode.v50.transcripts.fa.gz
#.fna from hpv no .zip or .gff3.gz



def check_fasta(file_path: str | Path):
    unique_ids = set()
    db = []

    # Supports both .gz and extracted FASTA files.
    path = Path(file_path)
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt", encoding="utf-8") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence_id = record.id.split("|")[0]
            sequence = str(record.seq).upper()

            duplicate = sequence_id in unique_ids
            unique_ids.add(sequence_id)

            db.append((
                sequence_id,
                sequence,
                len(sequence),
                sequence.count("N"),
                duplicate
            ))

    return db

# Check duplicate sequences separately.
# ☐ Flag empty sequences and other ambiguous/unexpected characters.
# ☐ Summarize minimum, maximum, and average length.
# ☐ Run on each intended file and save the results.


def process_fasta(records: list):
    db = pd.DataFrame(
        records,
        columns=["seq_id", "seq", "length", "N_count", "duplicate_ID"]
    )
    db["empty"] = db["length"] == 0
    db["dupe"] = db["seq"].duplicated()
    print(db["length"].agg(["mean", "max", "min"]))

    
    db["other_ambiguity_count"] = db["seq"].str.count("[RYSWKMBDHV]")
    db["unexpected_count"] = db["seq"].str.count("[^ACGTRYSWKMBDHVN]")
    



    return db





def main():
    data_folder = Path(
        r"C:\Users\chait\Downloads\BINF\BINF_2111\data"
    )

    output_folder = data_folder.parent / "results"
    output_folder.mkdir(parents=True, exist_ok=True)

    # Each input file has its own output filename.
    file_pairs = [
        (
            data_folder / "gencode.v50.lncRNA_transcripts.fa.gz",
            "lncrna_checks.csv"
        ),
        (
            data_folder / "gencode.v50.transcripts.fa.gz",
            "all_transcript_checks.csv"
        ),
        (
            data_folder / "HPV16_E7/ncbi_dataset/data/gene.fna",
            "HPV16_E7_checks.csv"
        ),
        (
            data_folder / "HPV16_E6/ncbi_dataset/data/gene.fna",
            "HPV16_E6_checks.csv"
        ),
    ]

    for input_path, output_name in file_pairs:
        print(f"\nReading: {input_path}")

        records = check_fasta(input_path)
        db = process_fasta(records)

        output_path = output_folder / output_name
        db.to_csv(output_path, index=False)

        print(f"Saved {len(db):,} records to {output_path}")

        # Free memory before reading the next file.
        del records, db


if __name__ == "__main__":
    main()

