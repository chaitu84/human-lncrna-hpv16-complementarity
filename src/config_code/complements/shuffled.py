import csv
import random
from pathlib import Path

import pandas as pd

from complementary import rev_comp, read_fasta, find_matches


WINDOW_SIZE = 10
WINDOW_STEP = 10      # Use 1 if Step 3 uses overlapping windows.
NUM_SHUFFLES = 100
RANDOM_SEED = 42
TEST_LIMIT = None     # Set to 100 for a smaller test run.


def shuffled(sequence: str, rng: random.Random) -> str:
    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def main():
    project_path = Path(
        "C:/Users/chait/Downloads/BINF/BINF_2111"
    )

    viral_path = (
        project_path
        / "src/config_code/hpv_step2/targets"
        / "HPV16_E6starI_E7_spliced_region.rna.fasta"
    )

    human_csv_path = project_path / "results/lncrna_checks.csv"

    output_dir = project_path / "results/step4"
    output_dir.mkdir(parents=True, exist_ok=True)

    matches_path = output_dir / "shuffled_matches.csv"
    summary_path = output_dir / "shuffle_summary.csv"

    # Read the full viral RNA once.
    # This uses the same read_fasta return format as your working code.
    record, sequence = read_fasta(viral_path)
    viral_sequence = str(sequence).strip().upper().replace("T", "U")
    viral_id = getattr(record, "id", str(record))

    if len(viral_sequence) < WINDOW_SIZE:
        raise ValueError("The viral sequence is shorter than WINDOW_SIZE.")

    # Read and normalize the human sequences once.
    human_data = pd.read_csv(
        human_csv_path,
        usecols=["seq_id", "seq"],
        dtype=str
    )

    human_data = human_data.dropna(subset=["seq_id", "seq"]).copy()
    human_data["seq"] = (
        human_data["seq"]
        .str.strip()
        .str.upper()
        .str.replace("T", "U", regex=False)
    )

    human_data = human_data[human_data["seq"].str.len() > 0]

    if TEST_LIMIT is not None:
        human_data = human_data.head(TEST_LIMIT)

    if human_data.empty:
        raise ValueError("No valid human sequences were found.")

    total = len(human_data)

    columns = [
        "shuffle_id",
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

    summary_columns = [
        "shuffle_id",
        "total_exact_matches",
        "human_sequences_with_matches",
        "human_sequences_searched",
        "viral_windows_searched",
        "window_size",
        "window_step",
        "random_seed"
    ]

    # Initialize once so successive shuffles differ,
    # while the complete run remains reproducible.
    rng = random.Random(RANDOM_SEED)

    # Stream results to disk instead of keeping every match in memory.
    with (
        matches_path.open("w", newline="", encoding="utf-8") as matches_file,
        summary_path.open("w", newline="", encoding="utf-8") as summary_file
    ):
        match_writer = csv.DictWriter(
            matches_file, fieldnames=columns
        )
        summary_writer = csv.DictWriter(
            summary_file, fieldnames=summary_columns
        )

        match_writer.writeheader()
        summary_writer.writeheader()

        for shuffle_id in range(1, NUM_SHUFFLES + 1):
            shuffled_sequence = shuffled(viral_sequence, rng)

            windows = []

            for start in range(
                0,
                len(shuffled_sequence) - WINDOW_SIZE + 1,
                WINDOW_STEP
            ):
                end = start + WINDOW_SIZE
                chunk = shuffled_sequence[start:end]

                windows.append(
                    (start, end, chunk, rev_comp(chunk))
                )

            match_count = 0
            human_sequences_with_matches = 0

            for number, row in enumerate(
                human_data.itertuples(index=False),
                start=1
            ):
                if number == 1 or number % 1000 == 0 or number == total:
                    print(
                        f"Shuffle {shuffle_id}/{NUM_SHUFFLES}: "
                        f"human sequence {number}/{total}",
                        flush=True
                    )

                hseq = row.seq
                has_match = False

                for viral_start, viral_end, chunk, comp in windows:
                    positions = find_matches(comp, hseq)

                    for human_start in positions:
                        human_end = human_start + WINDOW_SIZE

                        match_writer.writerow({
                            "shuffle_id": shuffle_id,
                            "hum_trspt_id": row.seq_id,
                            "viral_rna_id": viral_id,
                            "human_start": human_start,
                            "human_end": human_end,
                            "viral_start": viral_start,
                            "viral_end": viral_end,
                            "viral_window": chunk,
                            "viral_reverse_complement": comp,
                            "human_matched_sequence":
                                hseq[human_start:human_end],
                            "matches": WINDOW_SIZE,
                            "match_length": WINDOW_SIZE,
                            "mismatches": 0,
                            "gaps": 0
                        })

                        match_count += 1
                        has_match = True

                if has_match:
                    human_sequences_with_matches += 1

            summary_writer.writerow({
                "shuffle_id": shuffle_id,
                "total_exact_matches": match_count,
                "human_sequences_with_matches":
                    human_sequences_with_matches,
                "human_sequences_searched": total,
                "viral_windows_searched": len(windows),
                "window_size": WINDOW_SIZE,
                "window_step": WINDOW_STEP,
                "random_seed": RANDOM_SEED
            })

            matches_file.flush()
            summary_file.flush()

            print(
                f"Finished shuffle {shuffle_id}: "
                f"{match_count} exact matches",
                flush=True
            )

    print(f"\nMatches saved to: {matches_path}")
    print(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()