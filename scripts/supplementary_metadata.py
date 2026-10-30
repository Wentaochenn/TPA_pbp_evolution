#!/usr/bin/env python3
"""Export per-codon amino-acid states for the manuscript metadata table."""

import argparse
import collections
import csv
import json
import re
from pathlib import Path


METADATA_COLUMNS = (
    'sample', 'country', 'continent', 'year', 'lineage', 'bioproject',
    'biosample', 'study', 'source', 'clade', 'raw_lineage',
)
EXCLUDED_COLUMNS = (
    'included_primary_phylogeny', 'included_mrcA_selection_alignment',
    'included_pbp1_selection_alignment', 'included_pbp2_selection_alignment',
    'included_tp47_selection_alignment', 'included_strict_dating',
    'included_relaxed_dating', 'complete_mrcA_sequence', 'gene_haplotype',
)
SITE_GROUPS = (
    ('mrcA', ('G708', 'M625', 'A506', 'P44', 'I487', 'S759', 'R167',
              'V676', 'Q292', 'V451', 'T631', 'S339', 'A661', 'H209',
              'P448', 'D537', 'R641')),
    ('pbp1', ('P564', 'I16', 'Q145', 'T98')),
    ('tp47', ('S394', 'R107', 'R312', 'T27', 'S99', 'E184', 'R428')),
    ('lytB', ('A78', 'R147', 'V284', 'S303')),
    ('pbp2', ('I415', 'A366', 'E351', 'A549', 'Y291', 'G361', 'H510',
              'R340', 'K490', 'D520')),
)
EXPECTED_ALTERNATIVES = {
    'mrcA:G708': 'S', 'mrcA:M625': 'V', 'mrcA:A506': 'MTV',
    'mrcA:P44': 'L', 'mrcA:I487': 'L', 'mrcA:S759': 'P',
    'mrcA:R167': 'H', 'mrcA:V676': 'F', 'mrcA:Q292': 'R',
    'mrcA:V451': 'A', 'mrcA:T631': 'P', 'mrcA:S339': 'G',
    'mrcA:A661': 'V', 'mrcA:H209': 'R', 'mrcA:P448': 'L',
    'mrcA:D537': 'G', 'mrcA:R641': 'W',
    'pbp1:P564': 'L', 'pbp1:I16': 'M', 'pbp1:Q145': 'H', 'pbp1:T98': 'M',
    'tp47:S394': 'R', 'tp47:R107': 'G', 'tp47:R312': 'C',
    'tp47:T27': 'A', 'tp47:S99': 'P', 'tp47:E184': '*', 'tp47:R428': 'C',
    'lytB:A78': 'S', 'lytB:R147': 'H', 'lytB:V284': 'I', 'lytB:S303': 'G',
    'pbp2:I415': 'FM', 'pbp2:A366': 'TV', 'pbp2:E351': 'G',
    'pbp2:A549': 'D', 'pbp2:Y291': 'H', 'pbp2:G361': 'S',
    'pbp2:H510': 'N', 'pbp2:R340': 'C', 'pbp2:K490': 'E', 'pbp2:D520': 'N',
}
AA_ONE = dict(zip(
    ('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
     'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'Ter'),
    'ARNDCQEGHILKMFPSTWYV*',
))
MISSING = {'', 'NA', 'NAN', 'NULL', 'NONE', 'UNCALLABLE'}
FOOTNOTES = [
    'Amino-acid columns identify Nichols-reference residues and positions; '
    'each cell gives the observed one-letter residue. An asterisk denotes a '
    'termination codon. Reference and alternative residues are both shown. '
    'All alternatives at a multiallelic codon are evaluated jointly.',
    'For read-derived calls, a codon is callable when at least five accepted '
    'complete-codon observations support an unambiguous nucleotide codon at a '
    'frequency of at least 90%. Accepted alignments have mapping quality '
    '\u226520; complete-codon observations are excluded if any recorded '
    'codon-base quality is <20. Unmapped, secondary, '
    'supplementary, duplicate and quality-control-failed alignments are excluded. '
    'Observations sharing a query name count once when concordant; conflicting '
    'observations are discarded. Assembly-derived calls require one unambiguous '
    'aligned codon with mapping quality \u226520. Codons failing these criteria, '
    'or showing an unsupported amino-acid state or discordance between the '
    'sequence-derived variant and the BAM-supported state, are labelled '
    'uncallable. Missingness is assigned per codon, not to an entire gene.',
    'mrcA_protein_type reports the amino-acid sequence type assigned from an '
    'unambiguous complete mrcA coding sequence with a valid terminal stop codon. '
    'Unassigned protein types are shown as untyped; this designation is '
    'independent of callability at an individual codon. WT denotes the '
    'wild-type protein and Mutation #1\u2013#15 identify alternative protein types.',
    'Carrier totals count callable genomes whose residue differs from the '
    'reference residue. Callable totals exclude cells labelled uncallable. '
    'Summary counts describe individual positions, so alternative residues at '
    'the same position are counted together.',
]
FOOTNOTES_ZH = [
    '\u6c28\u57fa\u9178\u5217\u540d\u6807\u793aNichols\u53c2\u8003\u5e8f\u5217\u7684\u6b8b\u57fa\u53ca\u4f4d\u7f6e\uff1b\u5355\u5143\u683c\u4ee5\u5355\u5b57\u6bcd\u663e\u793a\u5b9e\u9645\u6c28\u57fa\u9178\u72b6\u6001\uff0c*\u8868\u793a\u7ec8\u6b62\u5bc6\u7801\u5b50\u3002\u540c\u4e00\u5bc6\u7801\u5b50\u7684\u591a\u79cd\u66ff\u4ee3\u72b6\u6001\u8054\u5408\u5224\u5b9a\u3002',
    '\u8bfb\u6bb5\u6570\u636e\u4e2d\uff0c\u53ef\u5224\u5b9a\u5bc6\u7801\u5b50\u9700\u81f3\u5c115\u4e2a\u5408\u683c\u7684\u5b8c\u6574\u5bc6\u7801\u5b50\u89c2\u6d4b\uff0c\u4e14\u4f18\u52bf\u6838\u82f7\u9178\u5bc6\u7801\u5b50\u6bd4\u4f8b\u226590%\u3002\u5408\u683c\u6bd4\u5bf9\u7684\u6bd4\u5bf9\u8d28\u91cf\u226520\uff1b\u5982\u5bc6\u7801\u5b50\u4e2d\u4efb\u4e00\u5df2\u8bb0\u5f55\u7684\u78b1\u57fa\u8d28\u91cf<20\uff0c\u5219\u5254\u9664\u8be5\u5b8c\u6574\u5bc6\u7801\u5b50\u89c2\u6d4b\u3002\u5254\u9664\u672a\u6bd4\u5bf9\u3001\u6b21\u8981\u3001\u8865\u5145\u3001\u91cd\u590d\u53ca\u8d28\u63a7\u672a\u901a\u8fc7\u7684\u6bd4\u5bf9\u3002\u540c\u4e00query name\u7684\u4e00\u81f4\u89c2\u6d4b\u53ea\u8ba1\u4e00\u6b21\uff0c\u51b2\u7a81\u89c2\u6d4b\u4e88\u4ee5\u5254\u9664\u3002\u7ec4\u88c5\u5e8f\u5217\u9700\u6709\u4e00\u4e2a\u65e0\u6b67\u4e49\u3001\u6bd4\u5bf9\u8d28\u91cf\u226520\u7684\u5b8c\u6574\u6bd4\u5bf9\u5bc6\u7801\u5b50\u3002\u4e0d\u6ee1\u8db3\u4e0a\u8ff0\u6761\u4ef6\u3001\u6c28\u57fa\u9178\u72b6\u6001\u4e0d\u53d7\u652f\u6301\uff0c\u6216\u5e8f\u5217\u53d8\u5f02\u4e0eBAM\u652f\u6301\u72b6\u6001\u4e0d\u4e00\u81f4\u8005\u6807\u8bb0\u4e3auncallable\u3002\u8be5\u6807\u8bb0\u5bf9\u5e94\u5355\u4e2a\u5bc6\u7801\u5b50\uff0c\u4e0d\u8868\u793a\u6574\u4e2a\u57fa\u56e0\u5747\u65e0\u6cd5\u5224\u5b9a\u3002',
    'mrcA_protein_type\u6839\u636e\u65e0\u6b67\u4e49\u4e14\u5177\u6709\u6709\u6548\u672b\u7aef\u7ec8\u6b62\u5bc6\u7801\u5b50\u7684\u5b8c\u6574mrcA\u7f16\u7801\u5e8f\u5217\u8fdb\u884c\u6c28\u57fa\u9178\u5e8f\u5217\u5206\u578b\u3002\u672a\u80fd\u5b8c\u6210\u86cb\u767d\u578b\u5224\u5b9a\u8005\u6807\u8bb0\u4e3auntyped\uff0c\u4e0e\u5355\u4e2a\u5bc6\u7801\u5b50\u7684\u53ef\u5224\u5b9a\u6027\u5206\u5f00\u8868\u793a\u3002WT\u8868\u793a\u91ce\u751f\u578b\uff0cMutation #1\u2013#15\u8868\u793a\u4e0d\u540c\u53d8\u5f02\u86cb\u767d\u578b\u3002',
    '\u643a\u5e26\u8005\u603b\u6570\u7edf\u8ba1\u53ef\u5224\u5b9a\u6837\u672c\u4e2d\u6c28\u57fa\u9178\u4e0d\u540c\u4e8e\u53c2\u8003\u6b8b\u57fa\u7684\u6837\u672c\u6570\uff1b\u53ef\u5224\u5b9a\u603b\u6570\u4e0d\u5305\u62ecuncallable\u3002\u540c\u4e00\u4f4d\u7f6e\u7684\u4e0d\u540c\u66ff\u4ee3\u6b8b\u57fa\u5408\u5e76\u8ba1\u6570\u3002',
]


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)):
            raise ValueError('Missing or duplicate table column names: ' + str(path))
        result = list(reader)
    if not result or any(None in row or any(v is None for v in row.values()) for row in result):
        raise ValueError('Empty or malformed table: ' + str(path))
    return fields, result


def keyed_rows(rows, label):
    result = {}
    for row in rows:
        sample = row.get('sample', '').strip()
        if not sample or sample in result:
            raise ValueError('Missing or duplicate sample in ' + label + ': ' + sample)
        result[sample] = row
    return result


def indicator_groups(fields):
    groups = collections.defaultdict(list)
    for field in fields:
        match = re.fullmatch(r'([^:]+):p\.([A-Za-z]+)(\d+)([A-Za-z]+)', field)
        if not match:
            continue
        gene, reference, position, alternative = match.groups()
        if gene not in dict(SITE_GROUPS):
            continue
        if reference not in AA_ONE or alternative not in AA_ONE:
            raise ValueError('Unsupported amino-acid name: ' + field)
        if reference == alternative:
            raise ValueError('Reference equals alternative: ' + field)
        groups[(gene, AA_ONE[reference] + position)].append((field, AA_ONE[alternative]))
    expected = {(gene, site) for gene, sites in SITE_GROUPS for site in sites}
    if set(groups) != expected:
        raise ValueError('Audited codon set differs from Table S1: missing=%s extra=%s' %
                         (sorted(expected - set(groups)), sorted(set(groups) - expected)))
    for key, indicators in groups.items():
        alternatives = [alternative for _, alternative in indicators]
        if len(set(alternatives)) != len(alternatives):
            raise ValueError('Duplicate alternative indicators at %s' % (key,))
        field = ':'.join(key)
        if set(alternatives) != set(EXPECTED_ALTERNATIVES[field]):
            raise ValueError('Audited alternatives differ at ' + field)
    return groups


def amino_acid_state(row, indicators, reference):
    values = []
    for field, alternative in indicators:
        value = str(row[field]).strip().upper()
        if value in MISSING:
            value = None
        elif value not in ('0', '1'):
            raise ValueError('Invalid ternary state for %s/%s: %r' %
                             (row.get('sample', '?'), field, row[field]))
        values.append((value, alternative))
    present = [alternative for value, alternative in values if value == '1']
    if len(present) > 1:
        raise ValueError('Mutually exclusive alternatives are both present in ' +
                         row.get('sample', '?'))
    if any(value is None for value, _ in values):
        return 'uncallable'
    return present[0] if present else reference


def build_table(metadata_rows, matrix_rows, matrix_fields):
    metadata = keyed_rows(metadata_rows, 'metadata')
    matrix = keyed_rows(matrix_rows, 'ternary matrix')
    if set(metadata) != set(matrix):
        raise ValueError('Metadata/matrix sample sets differ: metadata-only=%d matrix-only=%d' %
                         (len(set(metadata) - set(matrix)), len(set(matrix) - set(metadata))))
    required = set(METADATA_COLUMNS) | {'mrcA_protein_type'}
    if not all(required.issubset(row) for row in metadata.values()):
        raise ValueError('Metadata is missing a required publication field')
    groups = indicator_groups(matrix_fields)
    site_columns = [gene + ':' + site for gene, sites in SITE_GROUPS for site in sites]
    columns = list(METADATA_COLUMNS) + site_columns + ['mrcA_protein_type']
    rows = []
    for sample in sorted(metadata):
        source = metadata[sample]
        row = {field: source[field] for field in METADATA_COLUMNS}
        for gene, sites in SITE_GROUPS:
            for site in sites:
                row[gene + ':' + site] = amino_acid_state(matrix[sample], groups[(gene, site)], site[0])
        protein_type = str(source['mrcA_protein_type']).strip()
        row['mrcA_protein_type'] = 'untyped' if protein_type.upper() in MISSING else protein_type
        rows.append(row)
    summaries = []
    for gene, sites in SITE_GROUPS:
        for site in sites:
            field = gene + ':' + site
            states = collections.Counter(row[field] for row in rows)
            callable_count = len(rows) - states['uncallable']
            summaries.append({
                'field': field, 'gene': gene, 'site': site, 'reference': site[0],
                'carriers': callable_count - states[site[0]],
                'callable': callable_count, 'uncallable': states['uncallable'],
                'residue_counts': dict(sorted(states.items())),
            })
    display_labels = {'sample': 'Sample ID', 'year': 'isolation_year'}
    specs = []
    for field in columns:
        if ':' in field:
            gene, site = field.split(':')
            specs.append({'field': field, 'group': gene, 'label': site, 'kind': 'amino_acid'})
        else:
            specs.append({'field': field, 'group': '', 'label': display_labels.get(field, field),
                          'kind': 'protein_type' if field == 'mrcA_protein_type' else 'metadata'})
    spec = {
        'table': 1,
        'title': 'Genome metadata, amino-acid states and mrcA protein types',
        'title_zh': '\u57fa\u56e0\u7ec4\u5143\u6570\u636e\u3001\u6c28\u57fa\u9178\u72b6\u6001\u53camrcA\u86cb\u767d\u578b',
        'row_count': len(rows), 'column_count': len(columns),
        'amino_acid_site_count': len(site_columns), 'columns': specs,
        'gene_groups': [{'gene': gene, 'sites': list(sites)} for gene, sites in SITE_GROUPS],
        'footnotes': FOOTNOTES, 'footnotes_zh': FOOTNOTES_ZH, 'summaries': summaries,
        'protein_types_assigned': sum(row['mrcA_protein_type'] != 'untyped' for row in rows),
        'protein_types_untyped': sum(row['mrcA_protein_type'] == 'untyped' for row in rows),
    }
    return columns, rows, spec


def export(metadata_path, matrix_path, outdir, spec_path):
    _, metadata_rows = read_tsv(metadata_path)
    fields, matrix_rows = read_tsv(matrix_path)
    columns, rows, spec = build_table(metadata_rows, matrix_rows, fields)
    outdir, spec_path = Path(outdir), Path(spec_path)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / 'TableS1_sample_metadata.tsv'
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    spec_path.parent.mkdir(parents=True, exist_ok=True)
    spec_path.write_text(json.dumps(spec, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    notes = outdir / 'TableS1_footnotes.txt'
    notes.write_text('\n\n'.join([spec['title']] + FOOTNOTES) + '\n', encoding='utf-8')
    return spec


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--matrix', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--spec', required=True, type=Path)
    args = parser.parse_args(argv)
    spec = export(args.metadata, args.matrix, args.outdir, args.spec)
    print('Table S1: %d genomes, %d columns, %d audited amino-acid sites.' %
          (spec['row_count'], spec['column_count'], spec['amino_acid_site_count']))


if __name__ == '__main__':
    main()
