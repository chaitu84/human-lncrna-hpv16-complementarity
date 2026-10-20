# HPV16 Step 2: verified partial viral RNA target regions

Step 2 is complete: two reference-derived partial E6/E7 regions have been built and independently verified. These are computational sequence targets, not measured RNA molecules or complete transcripts. Step 3 was not performed. Constructing these targets does not demonstrate RNA binding or HPV suppression.

## Evidence and coordinate mapping

All genomic coordinates are **1-based, inclusive, on the plus strand**. Target positions are also 1-based, inclusive. The reference FASTAs are genomic DNA; the RNA FASTAs represent the same retained bases with T replaced by U.

The entire 776-base sequence **K02718.1:83-858** equals **NC_001526.4:7125-7900**, base for base, with no gaps. Exhaustive searches of both strands, including overlapping occurrences and the circular origin, find one plus-strand occurrence at 7125 and no reverse-complement occurrence. Therefore, within this aligned interval only:

```text
NC_001526.4 position = K02718.1 position + 7042
```

Do not extrapolate this offset genome-wide. K02718.1 has 7904 bases; NC_001526.4 has 7906, different sequence origins, and corrections elsewhere. Even the GenBank comment about identity to K02718 is not a substitute for this direct comparison.

[Tang et al. (2006)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1472016/), Materials and Methods, **“Plasmid construction and 293 transfection”**, describes the pTMF58 E6*I construct as positions 104-855 with intron bases **227-408 removed**. Thus **226 and 409 are retained**, and 408 - 227 + 1 = **182 bases** are deleted. The article also identifies the 226-to-409 splice in its text. Its graphical Figure 8 coordinates were not visually audited; this verification relies on the written methods and splice description.

Applying the verified local mapping gives retained bases **7268 and 7451**, removing **7269-7450**, again 182 bases. The included GenBank record independently annotates the plus-strand **CDS** `join(7125..7268,7451..7459)`, with gene E6, product protein E6*, and a note describing it as spliced. This coding feature ends at the E6* stop codon; it does **not** specify the entire E6*I/E7 RNA or its transcript boundaries.

[Paget-Bailly et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6459911/), Materials and Methods, **“Plasmids, transient and stable transfections”**, provides another sequence check: the upstream construct anchor's reverse complement ends at 7268, and the downstream anchor starts at 7451. The builder checks these exact published sequence anchors. They corroborate the junction but do not establish full transcript boundaries.

## Deliberate analysis boundaries

The targets use conventional **104-858**, equivalent to **NC_001526.4:7146-7900**. Tang's construct starts at 104; [Stacey et al. (2000)](https://pmc.ncbi.nlm.nih.gov/articles/PMC112249/), Results, **“The short E6 5′UTR is fully efficient in supporting initiation at the E6 AUG”**, explicitly identifies the E6 start codon at 104. The reference sequence also contains ATG at 104 / 7146.

Your original E6 FASTA instead follows NCBI's E6 gene/CDS annotation at **7125-7601** (477 bases, conventional 83-559). Selecting 7146 intentionally omits its first **21 bases**. This is an analysis choice matching the paper's E6 start, not a newly established transcription start site and not preservation of the entire NCBI E6 CDS.

The targets extend to **7900 / 858**, preserving the complete native E7 coding sequence and its stop codon TAA (UAA in RNA). Tang's FLAG-tagged construct ends at 855; these targets include the next three native bases and contain no tags or vector sequence. They are therefore reference-derived regions, not literal reconstructions of the experimental plasmid.

| Feature | Conventional / K02718.1 | NC_001526.4 |
|---|---|---|
| Original E6 gene/CDS | 83-559 | 7125-7601 |
| Selected analysis start | 104 | 7146 |
| Last retained base before intron | 226 | 7268 |
| Intron removed | 227-408 | 7269-7450 |
| First retained base after intron | 409 | 7451 |
| Full E7 CDS including stop | 562-858 | 7604-7900 |
| Selected analysis end | 858 | 7900 |

## Outputs and verified results

| Target | Retained NC_001526.4 intervals | Length | E7 positions within target |
|---|---|---|---|
| Intron-containing region | 7146-7900 | 755 nt | 459-755 |
| E6*I/E7 spliced region | 7146-7268 joined to 7451-7900 | 573 nt | 277-573 |

The spliced left segment has **123** bases and the right segment **450**. Its junction lies between target positions **123 and 124**:

```text
DNA: CAGTTACTGCGACGTGAG|GTGTATTAACTGTCAAAAGC
                     7268|7451 in NC_001526.4
```

Both targets preserve the full **297-base E7 sequence**, exactly matching your original E7 FASTA, plus the two native TC bases between the E6 and E7 gene annotations. The intron-containing region retains this E6 intron; it does not assert that the rest of a complete transcript is unspliced.

The four FASTA files are **two biological targets in two alphabets**, not four independent targets:

| Target | DNA representation | RNA representation |
|---|---|---|
| Intron-containing | [HPV16_E6E7_intron_retained_region.dna.fasta](../targets/HPV16_E6E7_intron_retained_region.dna.fasta) | [HPV16_E6E7_intron_retained_region.rna.fasta](../targets/HPV16_E6E7_intron_retained_region.rna.fasta) |
| Spliced | [HPV16_E6starI_E7_spliced_region.dna.fasta](../targets/HPV16_E6starI_E7_spliced_region.dna.fasta) | [HPV16_E6starI_E7_spliced_region.rna.fasta](../targets/HPV16_E6starI_E7_spliced_region.rna.fasta) |

Choose the alphabet required by a future tool and do not count both representations as separate observations.

- [targets/target_metadata.csv](../targets/target_metadata.csv): intervals in both coordinate systems, plus-strand orientation, lengths, junction and E7 positions, source/boundary notes, output paths, limitations, and DNA/RNA sequence hashes.
- [audit/validation.json](../audit/validation.json): checks performed by the builder, feature evidence, exact-match locations, and reference/target hashes.
- [audit/independent_validation.json](../audit/independent_validation.json): standalone independent sequence verification, including the original E6/E7 files.
- [audit/reference_retrieval.json](../audit/reference_retrieval.json): fresh NCBI download comparison, URLs, times, HTTP status, and raw file hashes.
- [audit/literature_evidence.json](../audit/literature_evidence.json): direct primary-source article checks and the precise sections supporting the splice, start, and anchors.
- [audit/preservation_baseline.json](../audit/preservation_baseline.json) and [audit/preservation_validation.json](../audit/preservation_validation.json): integrity records for original inputs and unrelated repository files.

The builder and independent audit both pass. The verifier constructs expected sequences by enumerating 1-based positions and does not import the builder or use its target slicing functions. It checks the generated FASTAs against both references and checks the GenBank ORIGIN independently. Validation also runs under Python `-O`; a deliberately corrupted target was rejected in memory without changing the input files.

## Reproduce within this repository

The package remains under `BINF_2111/src/config_code/hpv_step2/`; original inputs remain under `data/HPV16_E6/` and `data/HPV16_E7/`. Both scripts use Python's standard library (Python 3.10 or newer) and require no network or third-party packages. The builder writes its own four target FASTAs, metadata, and validation report. The verifier reads these and writes only its independent audit report.

Exact commands tested in PowerShell from `C:\Users\chait\Downloads\BINF\BINF_2111`:

```powershell
& 'C:\Users\chait\AppData\Local\Programs\Python\Python313\python.exe' src/config_code/hpv_step2/scripts/build_targets.py
& 'C:\Users\chait\AppData\Local\Programs\Python\Python313\python.exe' src/config_code/hpv_step2/scripts/verify_targets.py
```

This environment has Python **3.13.5** at that path; `python` is absent from PATH and `py -3` does not find the installed interpreter. On another machine with Python configured, the equivalent commands from the repository root are:

```bash
python src/config_code/hpv_step2/scripts/build_targets.py
python src/config_code/hpv_step2/scripts/verify_targets.py
```

For an inclusive biological interval `a..b`, Python uses `genome[a-1:b]`: subtract one from the start and keep the inclusive end as the exclusive slice endpoint. The actual target slices are:

```python
intron_containing = genome[7145:7900]
spliced = genome[7145:7268] + genome[7450:7900]
```

The independent audit additionally needs the original gene files in their existing repository paths; copying just this folder is sufficient for building, but not for verifying the original files.

## Provenance, corrections, and limitations

Fresh NCBI EFetch responses verified on **October 5, 2026, America/New_York** (October 6 UTC) exactly match all three included raw reference files byte for byte. Raw references were not overwritten. Endpoints and retrieval timestamps are recorded in `audit/reference_retrieval.json`:

- [NC_001526.4](https://www.ncbi.nlm.nih.gov/nuccore/NC_001526.4), FASTA and GenBank.
- [K02718.1](https://www.ncbi.nlm.nih.gov/nuccore/K02718.1), FASTA.
- [NCBI EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi), with `db=nuccore`, the accession version, `rettype=fasta` or `gbwithparts`, and `retmode=text`.
- Tang et al. (2006), DOI [10.1128/JVI.80.9.4249-4263.2006](https://doi.org/10.1128/JVI.80.9.4249-4263.2006).
- Paget-Bailly et al. (2019), DOI [10.1038/s41598-019-42393-6](https://doi.org/10.1038/s41598-019-42393-6).

No nucleotide or coordinate correction was needed: the existing target sequences were correct. The script was strengthened to check accession identity, actual GenBank sequence and CDS qualifiers, range bounds, and exhaustive orientation/uniqueness; explicit exceptions replace assertions that Python `-O` could skip. Metadata, provenance, reproduction instructions, and the independent audit were added or expanded.

The prior README's claim that a separately pasted whole genome matched cannot be rechecked because that pasted input is not present here; the claim was removed. The original E6 and E7 files present in this repository were independently checked and are exact matches to their annotated reference intervals. All five requested mapping/target claims were verified; no required coordinate mapping remains unresolved.

The targets represent one HPV16 reference and one published splice event. They do not cover all variants, alternative splice forms, UTRs, complete transcript ends, or integration-derived RNAs. The offline scripts validate included data; the separate retrieval and literature checks establish external provenance. No human-lncRNA matching, human-transcript BLAST, binding prediction, ranking, shuffled controls, or Step 3 implementation was performed during this verification. Unrelated repository work and original inputs were preserved.

## Folder organization

The package is located in `src/config_code/hpv_step2/`. Python scripts live in `scripts/`, target metadata lives beside the FASTAs in `targets/`, and the builder validation report lives with the other reports in `audit/`. This document preserves the detailed Step 2 evidence. The scripts and reproduction commands use the new paths; the scientific target definitions and sequence files are unchanged. Original inputs and unrelated repository files were not moved.

The Step 2 folder was moved out of `data/` into `src/config_code/`. Package references and outputs resolve relative to the Step 2 folder; the verifier locates the repository root and reads the original E6/E7 files from `data/HPV16_E6/` and `data/HPV16_E7/`. These file paths were updated; Python imports and neighboring reader scripts did not need changes.
