#!/usr/bin/env python3
"""Fixed-sample, fixed-root temporal diagnostic. Not a BEAST randomized-date test.
Single-date monophyletic clusters follow Murray et al. (2016),
DOI 10.1111/2041-210X.12466. Independent implementation; see README.
"""
import argparse, collections, copy, csv, hashlib, json, math, sys
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET
try:
    import numpy as np
    import Bio
    from Bio import Phylo
except ImportError as e:
    raise SystemExit('Use the existing tp_wgs environment with numpy and biopython. Missing dependency: '+str(e))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.common import tsv as write_tsv, sha

def xml_dates(path):
    root=ET.parse(path).getroot()
    traits=[x for x in root.iter() if x.get('traitname') in ('date','date-forward') and x.get('value')]
    if len(traits)!=1:raise ValueError('Expected exactly one forward calendar-date trait')
    out={}
    for item in traits[0].get('value').split(','):
        k,v=item.strip().rsplit('=',1);k=k.strip()
        if k in out:raise ValueError('Duplicate tip date: '+k)
        val=float(v)
        if not math.isfinite(val):raise ValueError('Non-finite date')
        out[k]=val
    taxa={x.get('taxon') for x in root.iter() if x.tag.split('}')[-1]=='sequence' and x.get('taxon')}
    if taxa and taxa!=set(out):raise ValueError('XML sequence taxon set differs from date trait set')
    return out

def validate_tree(tree):
    tips=tree.get_terminals();names=[t.name for t in tips]
    if None in names or len(names)!=len(set(names)):raise ValueError('Missing/duplicate tip labels')
    for c in tree.find_clades():
        b=c.branch_length
        if c is tree.root and b is None:continue
        if b is None or not math.isfinite(b) or b<0:raise ValueError('Missing, negative or non-finite edge length')
    if not any((t.branch_length or 0)>0 for t in tree.find_clades()):raise ValueError('Zero-length tree')

def cluster_indices(tree,dates):
    """Maximal rooted clades whose tips share exactly one recorded date."""
    names=[t.name for t in tree.get_terminals()];index={s:i for i,s in enumerate(names)}
    desc={};same={}
    for c in tree.find_clades(order='postorder'):
        if c.is_terminal():desc[c]=[index[c.name]];same[c]=dates[c.name]
        else:
            desc[c]=[i for ch in c.clades for i in desc[ch]]
            vals=[same[ch] for ch in c.clades]
            same[c]=vals[0] if vals and vals[0] is not None and all(v==vals[0] for v in vals) else None
    groups=[]
    stack=[tree.root]
    while stack:
        c=stack.pop()
        if same[c] is not None:groups.append(desc[c])
        else:stack.extend(reversed(c.clades))
    assert sorted(i for g in groups for i in g)==list(range(len(names)))
    return names,groups

def correlation(x,y):
    xc=x-x.mean();yc=y-y.mean();den=np.linalg.norm(xc)*np.linalg.norm(yc)
    if den==0:raise ValueError('Constant dates or root-to-tip distances: correlation undefined')
    return float(np.dot(xc,yc)/den)

def unique_permutations(values):
    counts=collections.Counter(values);keys=sorted(counts);out=[]
    def rec():
        if len(out)==len(values):yield tuple(out);return
        for k in keys:
            if counts[k]:
                counts[k]-=1;out.append(k);yield from rec();out.pop();counts[k]+=1
    yield from rec()

def permutation_test(x,year,groups,B,seed):
    vals=np.array([year[g[0]] for g in groups]);membership=np.empty(len(year),dtype=int)
    for j,g in enumerate(groups):
        if not np.all(year[g]==vals[j]):raise ValueError('Cluster contains mixed recorded dates')
        membership[g]=j
    observed=correlation(x,year)
    counts=collections.Counter(vals.tolist())
    logn=math.lgamma(len(vals)+1)-sum(math.lgamma(v+1) for v in counts.values())
    mode='exact' if logn<=math.log(20000)+1e-8 else 'monte_carlo'
    rng=np.random.default_rng(seed)
    assignments=unique_permutations(vals.tolist()) if mode=='exact' else (rng.permutation(vals) for _ in range(B))
    null=np.array([correlation(x,np.asarray(v)[membership]) for v in assignments])
    exceed=int(np.sum(null>=observed-1e-12))
    p=exceed/len(null) if mode=='exact' else (exceed+1)/(len(null)+1)
    return {'r':observed,'r_squared':observed**2,'p_one_sided':p,'permutation_mode':mode,
            'n_permutations':len(null),'n_units':len(groups),'log_unique_assignments':logn,
            'positive_slope':bool(observed>0),'minimum_attainable_p':1/len(null) if mode=='exact' else 1/(len(null)+1)},null
