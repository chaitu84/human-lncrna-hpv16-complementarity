"""Independently audit Step 2 offline, without importing build_targets.py.

Run from any directory with Python 3. Expected sequences are constructed by
enumerating 1-based reference coordinates, independently of the builder's slices.
Only audit/independent_validation.json is written; inputs and targets are read-only.
"""

from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import re
import sys


ROOT = Path(__file__).resolve().parent
REPORT = ROOT / "audit" / "independent_validation.json"
CHECKS = []
RAW_HASHES = {}


def check(condition, description):
    """Checks run even with Python -O; failures cannot silently become PASS."""
    CHECKS.append({"check": description, "passed": bool(condition)})
    if not condition:
        raise ValueError(description)


def raw_hash(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    RAW_HASHES[path.relative_to(ROOT.parent).as_posix()] = digest
    return digest


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_fasta(path, alphabet):
    lines = path.read_text(encoding="utf-8").splitlines()
    check(bool(lines) and lines[0].startswith(">"), f"FASTA header: {path.name}")
    check(sum(line.startswith(">") for line in lines) == 1,
          f"One FASTA record: {path.name}")
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    check(bool(sequence) and set(sequence) <= set(alphabet),
          f"FASTA alphabet {alphabet}: {path.name}")
    return lines[0][1:], sequence


def by_coordinates(sequence, first, last):
    if not 1 <= first <= last <= len(sequence):
        raise ValueError(f"Invalid 1-based interval {first}..{last}")
    return "".join(base for coordinate, base in enumerate(sequence, 1)
                   if first <= coordinate <= last)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def occurrences(genome, query, circular):
    """Test every possible starting base, including overlapping occurrences."""
    size = len(genome)
    search_sequence = genome + genome[:len(query) - 1] if circular else genome
    count = size if circular else size - len(query) + 1
    return [start + 1 for start in range(count)
            if search_sequence.startswith(query, start)]


def genbank_features(text):
    """Parse feature keys/locations separately from their qualifier text."""
    features = []
    in_features = False
    for line in text.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            break
        if not in_features:
            continue
        key = line[5:21].strip()
        content = line[21:].strip()
        if key:
            features.append({"type": key, "location": content, "body": []})
        elif features:
            features[-1]["body"].append(content)
    for feature in features:
        body = "\n".join(feature.pop("body"))
        feature["qualifiers"] = {
            key: " ".join((quoted if quoted else unquoted).split())
            for key, quoted, unquoted in re.findall(
                r'/([a-zA-Z_]+)=(?:"([^"]*)"|([^\s]+))', body, re.DOTALL)
        }
    return features


def annotated_feature(features, feature_type, location, gene, product=None):
    matches = [feature for feature in features
               if feature["type"] == feature_type
               and feature["location"] == location
               and feature["qualifiers"].get("gene") == gene]
    check(len(matches) == 1, f"GenBank {gene} {feature_type} plus-strand {location}")
    feature = matches[0]
    if product is not None:
        check(feature["qualifiers"].get("product") == product,
              f"GenBank {gene} product {product}")
        check(feature["qualifiers"].get("codon_start") == "1",
              f"GenBank {gene} CDS codon_start=1: {location}")
    return feature


def audit():
    reference_path = ROOT / "reference" / "NC_001526.4.fasta"
    conventional_path = ROOT / "reference" / "K02718.1.fasta"
    genbank_path = ROOT / "reference" / "NC_001526.4.gb"
    reference_header, reference = read_fasta(reference_path, "ACGT")
    conventional_header, conventional = read_fasta(conventional_path, "ACGT")
    check(reference_header.split()[0] == "NC_001526.4",
          "NC FASTA accession.version NC_001526.4")
    check(conventional_header.split()[0] == "K02718.1",
          "Conventional FASTA accession.version K02718.1")
    check("Human papillomavirus type 16" in reference_header
          and "Human papillomavirus type 16" in conventional_header,
          "Both FASTA headers identify HPV16")
    check(len(reference) == 7906 and len(conventional) == 7904,
          "Reference lengths: NC 7906 nt and K 7904 nt")
    text = genbank_path.read_text(encoding="utf-8")
    check(re.search(r"^VERSION\s+NC_001526\.4\s*$", text, re.MULTILINE)
          is not None, "GenBank VERSION NC_001526.4")
    check(re.search(r"^LOCUS\s+NC_001526\s+7906 bp\s+DNA\s+circular\b",
                    text, re.MULTILINE) is not None,
          "GenBank LOCUS: 7906-bp circular DNA")
    origin_match = re.search(r"^ORIGIN[^\n]*\n(.*?)^//", text,
                             re.MULTILINE | re.DOTALL)
    check(origin_match is not None, "GenBank ORIGIN present")
    origin = "".join(re.findall(r"[ACGTacgt]+", origin_match.group(1))).upper()
    check(origin == reference, "GenBank ORIGIN equals NC FASTA base for base")
    features = genbank_features(text)
    source = [feature for feature in features
              if feature["type"] == "source" and feature["location"] == "1..7906"]
    check(len(source) == 1
          and source[0]["qualifiers"].get("organism") == "Human papillomavirus 16"
          and source[0]["qualifiers"].get("mol_type") == "genomic DNA",
          "GenBank source identifies HPV16 genomic DNA")
    annotated_feature(features, "gene", "7125..7601", "E6")
    annotated_feature(features, "CDS", "7125..7601", "E6", "transforming protein E6")
    spliced_cds = annotated_feature(features, "CDS", "join(7125..7268,7451..7459)",
                                   "E6", "protein E6*")
    check("Spliced. E6*." in spliced_cds["qualifiers"].get("note", ""),
          "GenBank E6* feature explicitly described as spliced CDS")
    annotated_feature(features, "gene", "7604..7900", "E7")
    annotated_feature(features, "CDS", "7604..7900", "E7", "transforming protein E7")

    query = by_coordinates(conventional, 83, 858)
    check(len(query) == 776, "Mapping query length: 776 nt")
    mapping_hits = {}
    for circular in (False, True):
        mode = "circular" if circular else "linear"
        mapping_hits[mode] = {
            "+": occurrences(reference, query, circular),
            "-": occurrences(reference, reverse_complement(query), circular),
        }
        check(mapping_hits[mode] == {"+": [7125], "-": []},
              f"Exhaustive {mode} search: unique plus match at 7125, no minus match")
    check(query == by_coordinates(reference, 7125, 7900),
          "K02718.1:83..858 equals NC_001526.4:7125..7900 exactly")
    offset = mapping_hits["linear"]["+"][0] - 83
    check(offset == 7042, "Derived local offset +7042")
    check([position + offset for position in (104, 226, 227, 408, 409, 562, 858)]
          == [7146, 7268, 7269, 7450, 7451, 7604, 7900],
          "Mapped analysis, splice and E7 coordinates within verified local alignment")
    intron = by_coordinates(reference, 7269, 7450)
    check(len(intron) == 182 and intron.startswith("GT") and intron.endswith("AG"),
          "Removed interval: 182 nt with GT/AG reference intron boundaries")
    check(by_coordinates(reference, 7146, 7148) == "ATG"
          and by_coordinates(conventional, 104, 106) == "ATG",
          "Analysis starts at conventional 104 / NC 7146 ATG")
    check(by_coordinates(reference, 7457, 7459) == "TAA",
          "Annotated E6* CDS terminates at 7459 with TAA")

    originals = {}
    for gene, first, last in (("E6", 7125, 7601), ("E7", 7604, 7900)):
        path = ROOT.parent / f"HPV16_{gene}" / "ncbi_dataset" / "data" / "gene.fna"
        header, sequence = read_fasta(path, "ACGT")
        check(header.split()[0] == f"NC_001526.4:{first}-{last}"
              and header.split()[1] == gene,
              f"Original {gene} accession, gene name and coordinates")
        check(sequence == by_coordinates(reference, first, last),
              f"Original {gene} equals independently enumerated reference interval")
        originals[gene] = {"path": path.relative_to(ROOT.parent).as_posix(),
                           "length_nt": len(sequence), "raw_file_sha256": raw_hash(path)}

    e7 = by_coordinates(reference, 7604, 7900)
    check(len(e7) == 297 and e7.startswith("ATG") and e7.endswith("TAA"),
          "Native full E7: 297 nt with ATG start and TAA stop")
    check(by_coordinates(reference, 7602, 7603) == "TC",
          "Two native bases TC between annotated E6 and E7 retained")
    expected = {
        "HPV16_E6E7_intron_retained_region": {
            "sequence": by_coordinates(reference, 7146, 7900),
            "intervals": "7146-7900", "removed": "none", "length": 755,
            "conventional": by_coordinates(conventional, 104, 858),
            "e7_first": 459,
        },
        "HPV16_E6starI_E7_spliced_region": {
            "sequence": (by_coordinates(reference, 7146, 7268)
                         + by_coordinates(reference, 7451, 7900)),
            "intervals": "7146-7268;7451-7900", "removed": "7269-7450", "length": 573,
            "conventional": (by_coordinates(conventional, 104, 226)
                             + by_coordinates(conventional, 409, 858)),
            "e7_first": 277,
        },
    }
    targets = {}
    for name, specification in expected.items():
        dna = specification["sequence"]
        rna = "".join("U" if base == "T" else base for base in dna)
        check(len(dna) == specification["length"], f"Expected target length: {name}")
        check(dna == specification["conventional"],
              f"Independently reconstructed target agrees in both reference coordinate systems: {name}")
        for alphabet, sequence in (("DNA", dna), ("RNA", rna)):
            path = ROOT / "targets" / f"{name}.{alphabet.lower()}.fasta"
            header, actual = read_fasta(path, "ACGT" if alphabet == "DNA" else "ACGU")
            fields = header.split()
            attributes = dict(field.split("=", 1) for field in fields[1:] if "=" in field)
            check(fields[0] == name, f"Target identifier: {path.name}")
            for key, value in {
                "reference": "NC_001526.4", "strand": "+",
                "retained_1based": specification["intervals"],
                "partial_transcript_region": "true", "alphabet": alphabet,
            }.items():
                check(attributes.get(key) == value, f"Target header {key}: {path.name}")
            check(actual == sequence,
                  f"Target {alphabet} sequence equals independent coordinate enumeration: {name}")
        check(occurrences(dna, e7, False) == [specification["e7_first"]]
              and by_coordinates(dna, specification["e7_first"], len(dna)) == e7,
              f"Full E7 preserved at expected target position: {name}")
        check(dna.endswith("TAA") and rna.endswith("UAA"),
              f"Native E7 stop included in DNA/RNA representations: {name}")
        targets[name] = {
            "length_nt": len(dna), "retained_intervals_1based_inclusive": specification["intervals"],
            "dna_sequence_sha256": sequence_hash(dna), "rna_sequence_sha256": sequence_hash(rna),
            "e7_target_interval_1based_inclusive": f"{specification['e7_first']}-{len(dna)}",
            "e7_complete": True, "dna_rna_same_target": True,
        }

    left = by_coordinates(reference, 7146, 7268)
    right = by_coordinates(reference, 7451, 7900)
    spliced = expected["HPV16_E6starI_E7_spliced_region"]["sequence"]
    check(len(left) == 123 and by_coordinates(spliced, 1, 123) == left
          and by_coordinates(spliced, 124, len(spliced)) == right,
          "Spliced junction lies between target bases 123 and 124")
    context = by_coordinates(reference, 7251, 7268) + "|" + by_coordinates(reference, 7451, 7470)
    check(context == "CAGTTACTGCGACGTGAG|GTGTATTAACTGTCAAAAGC",
          "Spliced junction sequence context")

    with (ROOT / "target_metadata.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    check(len(rows) == 2 and {row["target"] for row in rows} == set(expected),
          "Metadata contains exactly the two biological targets")
    for row in rows:
        name = row["target"]
        specification = expected[name]
        values = {
            "reference": "NC_001526.4", "strand": "+",
            "retained_intervals_1based_inclusive": specification["intervals"],
            "length_nt": str(specification["length"]), "removed_interval": specification["removed"],
            "completeness": "partial transcript target region",
            "dna_sequence_sha256": sequence_hash(specification["sequence"]),
            "conventional_reference": "K02718.1",
            "conventional_retained_intervals_1based_inclusive":
                "104-858" if specification["removed"] == "none" else "104-226;409-858",
            "local_offset": "7042", "offset_valid_only_in": "K02718.1:83-858",
            "rna_sequence_sha256": targets[name]["rna_sequence_sha256"],
            "dna_fasta": f"targets/{name}.dna.fasta",
            "rna_fasta": f"targets/{name}.rna.fasta",
            "representation_relationship": "RNA is DNA with T replaced by U",
            "junction_target_1based": "" if specification["removed"] == "none" else "123|124",
            "E7_reference_interval": "7604-7900",
            "E7_target_interval": f"{specification['e7_first']}-{specification['length']}",
            "boundary_basis": "paper start 104; native E7 stop through 858",
            "splice_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC1472016/",
            "limitation": "reference-derived partial region; transcript ends and binding not established",
        }
        for key, value in values.items():
            check(row.get(key) == value, f"Metadata {key}: {name}")

    builder_validation = json.loads((ROOT / "validation.json").read_text(encoding="utf-8"))
    validation_values = {
        "status": "PASS", "reference_lengths": [7906, 7904],
        "alignment": "K02718.1:83-858 equals NC_001526.4:7125-7900 exactly (776/776 bases)",
        "literature_junction": "226|409", "mapped_junction": "7268|7451",
        "intron_length": 182, "spliced_length": 573, "unspliced_length": 755,
        "junction_after_target_base": 123, "junction_context_DNA": context,
        "E7_preserved": True,
        "reference_accessions": ["NC_001526.4", "K02718.1"],
        "genbank_origin_equals_fasta": True,
        "genbank_cds_features": [["7125..7601", "transforming protein E6"],
                                 ["join(7125..7268,7451..7459)", "protein E6*"],
                                 ["7604..7900", "transforming protein E7"]],
        "alignment_length_nt": 776, "alignment_plus_hits_1based": [7125],
        "alignment_reverse_complement_hits_1based": [],
        "local_offset": 7042, "offset_scope": "K02718.1:83-858 ONLY", "strand": "+",
        "removed_interval_1based_inclusive": "7269-7450",
        "retained_exon_lengths": [123, 450], "junction_before_target_base": 124,
        "E7_length_nt": 297, "E7_stop_codon_DNA": "TAA",
        "DNA_RNA_relationship": "same two targets; RNA alphabet replaces T with U",
        "target_sequence_sha256": {
            name: {"DNA": result["dna_sequence_sha256"], "RNA": result["rna_sequence_sha256"]}
            for name, result in targets.items()
        },
    }
    for key, value in validation_values.items():
        check(builder_validation.get(key) == value, f"Builder validation independently confirmed: {key}")
    for path in (reference_path, conventional_path, genbank_path):
        check(builder_validation.get("reference_file_sha256", {}).get(path.name) == raw_hash(path),
              f"Builder reference-file hash matches independently computed digest: {path.name}")

    return {
        "reference_identity": {"NC_FASTA": reference_header, "K_FASTA": conventional_header,
                               "GenBank_VERSION": "NC_001526.4", "GenBank_ORIGIN_matches_FASTA": True},
        "mapping": {"query": "K02718.1:83-858", "length_nt": 776,
                    "matched_interval": "NC_001526.4:7125-7900", "offset": offset,
                    "offset_scope": "Only the verified E6/E7 alignment; not genome-wide",
                    "exhaustive_exact_matches_by_strand": mapping_hits},
        "genbank_spliced_feature": spliced_cds,
        "original_inputs": originals,
        "targets": targets,
        "splice": {"retained_reference_bases": [7268, 7451], "removed_interval": "7269-7450",
                   "removed_length_nt": len(intron), "target_junction_between_bases": [123, 124],
                   "junction_context_DNA": context},
    }


def main():
    report = {
        "status": "FAIL", "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "offline": True, "independent_of_builder": True,
        "method": "Standalone standard-library verifier; 1-based coordinate enumeration; no builder import",
        "scope": "HPV Step 2 only; reference and target sequence verification",
        "limitations": [
            "Offline verification checks the included references and outputs; it does not fetch NCBI or reread papers.",
            "The published junction and analysis boundaries are supplied evidence; literature provenance is audited separately.",
            "The same source references support both scripts; independent computation does not replace external accession verification.",
            "These partial target regions establish neither complete transcript ends nor RNA binding or HPV suppression.",
        ],
    }
    try:
        report.update(audit())
        report["status"] = "PASS"
        exit_code = 0
    except (OSError, ValueError, KeyError, IndexError, TypeError) as error:
        report["error"] = str(error)
        exit_code = 1
    report["checks"] = CHECKS
    report["raw_input_file_sha256"] = RAW_HASHES
    report["checks_passed"] = sum(item["passed"] for item in CHECKS)
    report["checks_total"] = len(CHECKS)
    REPORT.parent.mkdir(exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Independent Step 2 audit: {report['status']} ({report['checks_passed']}/{report['checks_total']} checks).")
    print(f"Report: {REPORT}")
    if exit_code:
        print(report["error"], file=sys.stderr)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
