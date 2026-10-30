"""Read site probabilities from a complete FUBAR result."""
import json
import math
from pathlib import Path

def json_sites(path,method,ncodons,ntips=None):
    d=json.loads(Path(path).read_text());inp=d['input']
    if int(inp['number of sites'])!=ncodons:raise ValueError('HyPhy site count mismatch')
    if ntips is not None and int(inp['number of sequences'])!=ntips:raise ValueError('HyPhy sample count mismatch')
    headers=[str(x[0]) if isinstance(x,list) else str(x) for x in d['MLE']['headers']]
    if method!='fubar':raise ValueError('Only FUBAR result files are supported')
    key='Prob[alpha<beta]'
    if headers.count(key)!=1:raise ValueError('Missing/ambiguous HyPhy header '+key)
    idx=headers.index(key);out={};content=d['MLE']['content']
    for part,rows in content.items():
        coverage=d['data partitions'][str(part)]['coverage']
        positions=[int(i) for group in coverage for i in group]
        if len(rows)!=len(positions):raise ValueError('Coverage/row count mismatch')
        for pos,row in zip(positions,rows):
            if not isinstance(row,list) or len(row)<len(headers):raise ValueError('Short HyPhy result row')
            v=float(row[idx])
            if not math.isfinite(v) or not 0<=v<=1:raise ValueError('Invalid HyPhy probability')
            if pos in out:raise ValueError('Duplicate result site')
            out[pos]=v
    if set(out)!=set(range(ncodons)):raise ValueError('Missing/out-of-range result sites')
    return out
