"""Project deposited selection inputs onto the manuscript's FUBAR data schema."""
import csv
import json
import math
from collections import Counter
from pathlib import Path

from integrity import sha

GENES = {'mrcA': 884, 'pbp1': 624, 'pbp2': 652, 'tp47': 434}
SELECTION_FIELDS = (
    'gene', 'site', 'FUBAR_Ppositive', 'FUBAR_Pnegative',
    'FUBAR_BayesFactor_positive', 'FUBAR_Ppositive_above_09',
    'callable_codon_sequences', 'observed_codon_types',
    'amino_acid_counts_in_selection_alignment', 'matrix_indicators',
)
FOCAL_FIELDS = ('locus',) + SELECTION_FIELDS + (
    'n_mutant', 'n_wildtype', 'n_missing', 'simmap_gain_mean',
    'simmap_gain_q025', 'simmap_gain_q975', 'minimum_gains_among_optima',
    'maximum_gains_among_optima', 'optimal_root_states',
    'delta_AIC_ARD_ER', 'AIC_evidence', 'ARD_root_Pmut',
)
TOPOLOGY_FIELDS = (
    'gene', 'codon', 'planned_bootstrap_trees', 'complete_bootstrap_trees',
    'primary_fubar_posterior', 'fubar_gt_09_count', 'fubar_gt_09_fraction',
    'stability_fraction_status', 'fubar_posterior_positive_min',
    'fubar_posterior_positive_q025', 'fubar_posterior_positive_median',
    'fubar_posterior_positive_q975', 'fubar_posterior_positive_max',
)


def read_table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t', escapechar='\\')
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError('Malformed selection input: ' + str(path))
        return reader.fieldnames, rows


def keep_source_key(key):
    if key.endswith('.json') and key.split('_', 1)[0] in GENES:
        return key in {gene + '_fubar.json' for gene in GENES}
    return True


def project_table(path, allowed, required):
    fields, rows = read_table(path)
    if not set(required) <= set(fields):
        raise ValueError('Required FUBAR fields are missing from ' + str(path))
    selected = [field for field in allowed if field in fields]
    temporary = Path(str(path) + '.partial')
    with temporary.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=selected, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows({field: row[field] for field in selected} for row in rows)
    temporary.replace(path)
    return rows


def project(run):
    """Preserve all FUBAR values while restricting newly copied plotting inputs."""
    run = Path(run).resolve()
    data = run / 'work/figures/data'
    fields, rows = read_table(data / 'selection.tsv')
    if Counter(row['gene'] for row in rows) != GENES:
        raise ValueError('Selection data must contain 2,594 codons across the four genes')
    if len({(row['gene'], int(row['site'])) for row in rows}) != 2594:
        raise ValueError('Selection data contain duplicate codons')
    for row in rows:
        value = float(row['FUBAR_Ppositive'])
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('Invalid FUBAR posterior probability')
    project_table(data / 'selection.tsv', SELECTION_FIELDS, ('gene', 'site', 'FUBAR_Ppositive'))
    project_table(data / 'focal.tsv', FOCAL_FIELDS, ('locus', 'gene', 'site', 'FUBAR_Ppositive'))
    topology = data / 'topology_stability.tsv'
    if topology.is_file():
        project_table(topology, TOPOLOGY_FIELDS,
                      ('gene', 'codon', 'primary_fubar_posterior', 'complete_bootstrap_trees',
                       'fubar_gt_09_count', 'fubar_posterior_positive_min',
                       'fubar_posterior_positive_median', 'fubar_posterior_positive_max'))
    source_path = run / 'SOURCE_MANIFEST.json'
    sources = {key: value for key, value in json.loads(source_path.read_text()).items()
               if keep_source_key(key)}
    if topology.is_file() and 'TOPOLOGY' in sources:
        sources['TOPOLOGY'] = dict(path=str(topology), sha256=sha(topology), bytes=topology.stat().st_size)
    source_path.write_text(json.dumps(sources, indent=2) + '\n')
    resolution = run / 'INPUT_RESOLUTION.tsv'
    if resolution.is_file():
        fields, records = read_table(resolution)
        with resolution.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for row in records:
                if not keep_source_key(row['input']):
                    continue
                if row['input'] == 'TOPOLOGY' and topology.is_file():
                    row['path'] = str(topology)
                writer.writerow(row)
    manifest = {path.name: sha(path) for path in data.iterdir() if path.is_file()}
    (run / 'work/INPUT_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    scope = dict(method='FUBAR', codons=2594, posterior_threshold=0.90,
                 selection_models_refitted=False)
    (run / 'FUBAR_INPUT_SCOPE.json').write_text(json.dumps(scope, indent=2) + '\n')
    return scope
