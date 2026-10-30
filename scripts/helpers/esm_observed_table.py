#!/usr/bin/env python3
"""Export four-gene scores for replacements observed in the study matrix."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

GENES = {'mrcA': 'TP0705', 'pbp1': 'TP0500', 'pbp2': 'TP0760', 'tp47': 'TP0574'}
AA3 = dict(zip(('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter').split(),
               'ARNDCQEGHILKMFPSTWYV*'))
MISSING = {'', 'NA', 'NaN', 'nan', '.', '-', 'None'}


def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def observed_counts(path):
    rows = read(path)
    if not rows or 'sample' not in rows[0]:
        raise ValueError('Expected a nonempty sample-by-variant matrix')
    if len({r['sample'] for r in rows}) != len(rows):
        raise ValueError('Duplicate sample identifiers')
    result = []
    for locus in rows[0]:
        match = re.fullmatch(r'([a-zA-Z0-9]+):p\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})', locus)
        if not match or match[1] not in GENES:
            continue
        states = []
        for row in rows:
            value = row[locus]
            if value in MISSING:
                states.append(None)
            elif value in {'0', '0.0'}:
                states.append(0)
            elif value in {'1', '1.0'}:
                states.append(1)
            else:
                raise ValueError('Invalid matrix state at ' + locus)
        carriers = states.count(1)
        if carriers:
            result.append({'gene': match[1], 'locus_tag': GENES[match[1]], 'locus': locus,
                           'position': int(match[3]), 'reference_aa': AA3[match[2]],
                           'alternate_aa': AA3[match[4]], 'carriers': carriers,
                           'callable_genomes': len(states) - states.count(None)})
    return sorted(result, key=lambda r: (list(GENES).index(r['gene']), r['position'], r['alternate_aa']))


def build(matrix, scores):
    observed = observed_counts(matrix)
    full = {}
    for row in read(scores):
        if row['protein'] not in GENES or row.get('scan_mode', 'full_dms') != 'full_dms':
            continue
        key = row['protein'], int(row['position']), row['wt_aa'], row['mut_aa']
        value = float(row['llr'])
        if key in full or not math.isfinite(value):
            raise ValueError('Duplicate or non-finite score')
        full[key] = value
    for row in observed:
        key = row['gene'], row['position'], row['reference_aa'], row['alternate_aa']
        if key not in full:
            raise ValueError('No archived score for observed replacement: ' + row['locus'])
        row['replacement'] = row['reference_aa'] + str(row['position']) + row['alternate_aa']
        row['carrier_percent_callable'] = 100 * row['carriers'] / row['callable_genomes']
        row['ESM2_LLR'] = full[key]
    return observed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a new supplementary output directory')
    rows = build(args.matrix, args.scores)
    if not rows:
        raise ValueError('No observed replacements found')
    args.out.mkdir(parents=True)
    target = args.out / 'observed_PBP_ESM2_scores.tsv'
    fields = ['gene', 'locus_tag', 'locus', 'replacement', 'position', 'reference_aa', 'alternate_aa',
              'carriers', 'callable_genomes', 'carrier_percent_callable', 'ESM2_LLR']
    with target.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    report = {'state': 'complete', 'rows': len(rows), 'genes': list(GENES),
              'rows_by_gene': {g: sum(r['gene'] == g for r in rows) for g in GENES},
              'matrix_sha256': digest(args.matrix), 'scores_sha256': digest(args.scores),
              'table_sha256': digest(target), 'all_position_scores_in_supplementary_table': False}
    (args.out / 'ESM_SUPPLEMENT_STATUS.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
