"""Validate the model comparisons and stochastic-map records used for plotting."""
from collections import defaultdict
import math
from pathlib import Path


def validate_model_mappings(run):
    from prepare import read
    data = Path(run)/'work/figures/data'
    summary = read(data/'er_ard.tsv')
    if len(summary) != 23 or len({r['locus'] for r in summary}) != 23:
        raise ValueError('Expected 23 distinct ER/ARD comparisons')
    if any(r.get('analysis_status') != 'complete' for r in summary):
        raise ValueError('All ER/ARD comparisons must be complete')
    draws = read(data/'simmap_draws.tsv')
    groups = defaultdict(list)
    for row in draws:
        groups[row['locus'], row['model']].append(row)
    expected = {(r['locus'], m): int(r['mapping_n']) for r in summary for m in ['ER','ARD']}
    if set(groups) != set(expected):
        raise ValueError('Model comparison and stochastic-map scopes differ')
    for key, n in expected.items():
        rows = groups[key]
        if n < 2 or len(rows) != n or {int(r['map_index']) for r in rows} != set(range(1,n+1)):
            raise ValueError('Invalid stochastic-map count or indices: '+str(key))
        for row in rows:
            for name in ['gain','loss']:
                value = float(row[name])
                if not math.isfinite(value) or value < 0 or value != int(value):
                    raise ValueError('Invalid transition count: '+str(key))
    return {'state': 'PASS', 'variants': len(summary), 'models': len(groups), 'maps': len(draws)}
