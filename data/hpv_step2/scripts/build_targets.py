"""Build only the two partial HPV16 Step 2 targets; offline, standard library."""
from pathlib import Path
import csv
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent
TANG_URL = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC1472016/'
ANCHOR_URL = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6459911/'


def require(condition, message):
    """Validation must also run when Python is invoked with -O."""
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fasta(path, accession):
    lines = path.read_text(encoding='utf-8').splitlines()
    require(bool(lines) and lines[0].startswith('>' + accession + ' '),
            f'{path.name}: expected accession {accession}')
    require(sum(line.startswith('>') for line in lines) == 1,
            f'{path.name}: expected one FASTA record')
    seq = ''.join(line.strip() for line in lines[1:]).upper()
    require(bool(seq) and set(seq) <= set('ACGT'), f'{path.name}: invalid DNA')
    return seq


def region(seq, start, end):
    """Convert 1-based inclusive coordinates to a 0-based half-open slice."""
    require(1 <= start <= end <= len(seq), f'Out-of-range interval {start}..{end}')
    return seq[start - 1:end]


def all_hits(query, genome):
    """Exhaustive overlapping matches, including across the circular origin."""
    circular = genome + genome[:len(query) - 1]
    return [i + 1 for i in range(len(genome)) if circular.startswith(query, i)]


def read_genbank(path, reference):
    gb = path.read_text(encoding='utf-8')
    require(re.search(r'^VERSION\s+NC_001526\.4\s*$', gb, re.M),
            'Unexpected GenBank version')
    require(re.search(r'^LOCUS\s+NC_001526\s+7906 bp\s+DNA\s+circular', gb, re.M),
            'Unexpected GenBank length, molecule or topology')
    require('Human papillomavirus type 16, complete genome.' in gb,
            'Unexpected GenBank identity')
    require(gb.count('\nORIGIN') == 1, 'Expected one GenBank ORIGIN')
    origin = gb.split('\nORIGIN', 1)[1].split('//', 1)[0]
    sequence = ''.join(re.findall(r'[A-Za-z]+', origin)).upper()
    require(sequence == reference, 'GenBank ORIGIN differs from NC FASTA')
    feature_text = gb.split('\nFEATURES', 1)[1].split('\nORIGIN', 1)[0]
    features = []
    for block in re.split(r'\n(?=     \S)', feature_text)[1:]:
        lines = block.splitlines()
        key, location = lines[0].strip().split(None, 1)
        qualifiers = '\n'.join(lines[1:])
        features.append((key, location, qualifiers))
    expected = [('7125..7601', 'transforming protein E6'),
                ('join(7125..7268,7451..7459)', 'protein E6*'),
                ('7604..7900', 'transforming protein E7')]
    for location, product in expected:
        matches = [q for key, loc, q in features if key == 'CDS' and loc == location]
        require(len(matches) == 1 and f'/product="{product}"' in matches[0],
                f'Missing or ambiguous plus-strand CDS {location}: {product}')
        require('/codon_start=1' in matches[0], f'Unexpected CDS frame: {location}')
        if product == 'protein E6*':
            require('/gene="E6"' in matches[0] and 'Spliced. E6*.' in matches[0],
                    'Spliced E6* feature qualifiers do not match')
    return expected


def write_fasta(path, name, seq, description):
    path.write_text('>' + name + ' ' + description + '\n' +
                    '\n'.join(seq[i:i + 70] for i in range(0, len(seq), 70)) + '\n',
                    encoding='utf-8', newline='\n')


def main():
    ref = fasta(ROOT / 'reference/NC_001526.4.fasta', 'NC_001526.4')
    lit = fasta(ROOT / 'reference/K02718.1.fasta', 'K02718.1')
    require(len(ref) == 7906 and len(lit) == 7904, 'Unexpected reference lengths')
    cds_features = read_genbank(ROOT / 'reference/NC_001526.4.gb', ref)
    query = region(lit, 83, 858)
    hits = all_hits(query, ref)
    reverse_query = query.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    reverse_hits = all_hits(reverse_query, ref)
    require(hits == [7125] and not reverse_hits,
            f'Region not unique in the expected orientation: {hits}, {reverse_hits}')
    offset = hits[0] - 83
    require(offset == 7042, 'Unexpected local coordinate offset')
    # This offset is established ONLY for K02718.1:83..858.
    start, donor, acceptor, end = [p + offset for p in (104, 226, 409, 858)]
    require((start, donor, acceptor, end) == (7146, 7268, 7451, 7900),
            'Unexpected mapped analysis/splice boundaries')
    intron = region(ref, donor + 1, acceptor - 1)
    require(len(intron) == 182 and intron.startswith('GT') and intron.endswith('AG'),
            'Unexpected intron length or splice-site motifs')
    # Published construct sequence anchors are evidence checks, not a lab protocol.
    upstream = 'CTCACGTCGCAGTAACTG'.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    downstream = 'GTGTATTAACTGTCAAAAGCC'
    require(region(ref, donor - len(upstream) + 1, donor) == upstream,
            'Published upstream junction anchor mismatch')
    require(region(ref, acceptor, acceptor + len(downstream) - 1) == downstream,
            'Published downstream junction anchor mismatch')
    left, right = region(ref, start, donor), region(ref, acceptor, end)
    unspliced, spliced = region(ref, start, end), left + right
    require((len(left), len(right), len(unspliced), len(spliced)) == (123, 450, 755, 573),
            'Unexpected exon/target lengths')
    require(region(lit, 104, 858) == unspliced and
            region(lit, 104, 226) + region(lit, 409, 858) == spliced,
            'Target differs between reference coordinate systems')
    e7 = region(ref, 7604, 7900)
    require(len(e7) == 297 and spliced.endswith(e7) and unspliced.endswith(e7),
            'Full E7 sequence not preserved')
    require(e7.startswith('ATG') and e7.endswith('TAA'), 'E7 start/stop mismatch')
    require(region(ref, 7602, 7603) == 'TC', 'Unexpected intergenic bases')
    require(left[-1] == ref[7267] and right[0] == ref[7450],
            'Splice junction retained-base mismatch')
    out = ROOT / 'targets'
    out.mkdir(exist_ok=True)
    rows = []
    target_hashes = {}
    for name, seq, kept, conventional, e7_start, junction in [
        ('HPV16_E6E7_intron_retained_region', unspliced, '7146-7900', '104-858', 459, ''),
        ('HPV16_E6starI_E7_spliced_region', spliced, '7146-7268;7451-7900',
         '104-226;409-858', 277, '123|124')
    ]:
        description = (f'reference=NC_001526.4 strand=+ retained_1based={kept} '
                       'partial_transcript_region=true')
        paths, hashes = {}, {}
        for alphabet, sequence in [('DNA', seq), ('RNA', seq.replace('T', 'U'))]:
            path = out / f'{name}.{alphabet.lower()}.fasta'
            write_fasta(path, name, sequence, description + f' alphabet={alphabet}')
            paths[alphabet] = path.relative_to(ROOT).as_posix()
            hashes[alphabet] = sha256(sequence.encode('ascii'))
        target_hashes[name] = hashes
        rows.append(dict(target=name, reference='NC_001526.4', strand='+',
                         retained_intervals_1based_inclusive=kept, length_nt=len(seq),
                         removed_interval='7269-7450' if junction else 'none',
                         completeness='partial transcript target region',
                         dna_sequence_sha256=hashes['DNA'],
                         conventional_reference='K02718.1',
                         conventional_retained_intervals_1based_inclusive=conventional,
                         local_offset=offset, offset_valid_only_in='K02718.1:83-858',
                         rna_sequence_sha256=hashes['RNA'], dna_fasta=paths['DNA'],
                         rna_fasta=paths['RNA'], representation_relationship='RNA is DNA with T replaced by U',
                         junction_target_1based=junction,
                         E7_reference_interval='7604-7900',
                         E7_target_interval=f'{e7_start}-{len(seq)}',
                         boundary_basis='paper start 104; native E7 stop through 858',
                         splice_source=TANG_URL,
                         limitation='reference-derived partial region; transcript ends and binding not established'))
    with (ROOT / 'target_metadata.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    checks = dict(status='PASS', validation_kind='builder checks; run verify_targets.py for independent audit',
                  reference_lengths=[len(ref), len(lit)],
                  reference_accessions=['NC_001526.4', 'K02718.1'],
                  genbank_origin_equals_fasta=True, genbank_cds_features=cds_features,
                  alignment='K02718.1:83-858 equals NC_001526.4:7125-7900 exactly (776/776 bases)',
                  alignment_length_nt=776, alignment_plus_hits_1based=hits,
                  alignment_reverse_complement_hits_1based=reverse_hits,
                  uniqueness_search='all overlapping starts on both strands, including circular origin',
                  local_offset=offset, offset_scope='K02718.1:83-858 ONLY', strand='+',
                  literature_junction='226|409', mapped_junction='7268|7451',
                  removed_interval_1based_inclusive='7269-7450', intron_length=182,
                  spliced_length=573, unspliced_length=755, retained_exon_lengths=[123, 450],
                  junction_after_target_base=123, junction_before_target_base=124,
                  junction_context_DNA=left[-18:] + '|' + right[:20],
                  splice_sources=[TANG_URL, ANCHOR_URL],
                  independent_support='NCBI spliced CDS annotation and published deletion-construct sequence anchors',
                  E7_preserved=True, E7_length_nt=297, E7_stop_codon_DNA='TAA',
                  DNA_RNA_relationship='same two targets; RNA alphabet replaces T with U',
                  target_sequence_sha256=target_hashes,
                  reference_file_sha256={p.name: sha256(p.read_bytes()) for p in
                                         sorted((ROOT / 'reference').iterdir()) if p.is_file()},
                  limitations=['partial regions, not complete transcript boundaries',
                               'no binding or suppression evidence', 'no Step 3 analysis'])
    (ROOT / 'validation.json').write_text(json.dumps(checks, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
