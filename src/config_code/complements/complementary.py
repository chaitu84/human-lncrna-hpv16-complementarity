from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


WINDOW_SIZE = 10
TEST_LIMIT = None  # Change to None to process all human transcripts


def rev_comp(strand: str) -> str:
    return str(Seq(strand).reverse_complement_rna())


def read_fasta(path: Path) -> tuple[str, str]:
    record = SeqIO.read(path, "fasta")
    sequence = str(record.seq).upper().replace("T", "U")

    return record.id, sequence


def find_matches(query: str, human_sequence: str) -> list[int]:
    positions = []
    start = 0

    while True:
        position = human_sequence.find(query, start)

        if position == -1:
            break

        positions.append(position)

        # Move one base forward to include overlapping occurrences.
        start = position + 1

    return positions


def main():
    project = Path(
        "C:/Users/chait/Downloads/BINF/BINF_2111"
    )

    viral_path = (
        project
        / "src/config_code/hpv_step2/targets"
        / "HPV16_E6starI_E7_spliced_region.rna.fasta"
    )
    human_path = project / "results/lncrna_checks.csv"
    output_path = project / "results/exact_complementarity_matches.csv"

    viral_id, viral_sequence = read_fasta(viral_path)

    human_rna_csv = pd.read_csv(
        human_path,
        usecols=["seq_id", "seq"],
        dtype=str
    )

    if TEST_LIMIT is not None:
        human_rna_csv = human_rna_csv.head(TEST_LIMIT)

    # Prepare overlapping viral windows and their reverse complements.
    windows = []

    for viral_start in range(len(viral_sequence) - WINDOW_SIZE + 1):
        viral_end = viral_start + WINDOW_SIZE
        chunk = viral_sequence[viral_start:viral_end]

        # Only search windows containing unambiguous RNA bases.
        if any(base not in "ACGU" for base in chunk):
            continue

        windows.append(
            (viral_start, viral_end, chunk, rev_comp(chunk))
        )

    columns = [
        "hum_trspt_id",
        "viral_rna_id",
        "human_start",
        "human_end",
        "viral_start",
        "viral_end",
        "viral_window",
        "viral_reverse_complement",
        "human_matched_sequence",
        "matches",
        "match_length",
        "mismatches",
        "gaps"
    ]

    rows = []
    total = len(human_rna_csv)

    for number, row in enumerate(
        human_rna_csv.itertuples(index=False),
        start=1
    ):
        if number == 1 or number % 100 == 0 or number == total:
            print(
                f"Processing human sequence {number}/{total}",
                flush=True
            )

        if pd.isna(row.seq_id) or pd.isna(row.seq):
            continue

        hseq = row.seq.strip().upper().replace("T", "U")

        for viral_start, viral_end, chunk, comp in windows:
            positions = find_matches(comp, hseq)

            for human_start in positions:
                human_end = human_start + WINDOW_SIZE

                rows.append({
                    "hum_trspt_id": row.seq_id,
                    "viral_rna_id": viral_id,
                    "human_start": human_start,
                    "human_end": human_end,
                    "viral_start": viral_start,
                    "viral_end": viral_end,
                    "viral_window": chunk,
                    "viral_reverse_complement": comp,
                    "human_matched_sequence": hseq[human_start:human_end],
                    "matches": WINDOW_SIZE,
                    "match_length": WINDOW_SIZE,
                    "mismatches": 0,
                    "gaps": 0
                })

    sequences = pd.DataFrame(rows, columns=columns)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sequences.to_csv(output_path, index=False)

    print(sequences.head().to_string(index=False))
    print(f"\nTotal exact matches: {len(sequences)}")
    print(f"Saved to: {output_path}")

    return sequences


if __name__ == "__main__":
    results = main()