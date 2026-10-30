"""Describe binary-state histories on the supplied, unchanged rooted tree.

No rate fitting, sequence prediction, model selection or phenotypic inference.
Bounds are exact marginal bounds among globally minimum-change histories.
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from treeio import read_tree

BASE = Path(__file__).resolve().parent
INF = 10**9


def fit(tree, observations, edge_regions, root_state=None):
    costs = {}
    for n in tree.postorder:
        node = tree.nodes[n]
        v = observations.get(node.label) if not node.children else None
        costs[n] = [
            (0 if v is None or v == s else INF) +
            sum(min(costs[c][t] + (s != t) for t in (0, 1)) for c in node.children)
            for s in (0, 1)
        ]
    allowed = (0, 1) if root_state is None else (root_state,)
    score = min(costs[tree.root][s] for s in allowed)
    roots = [s for s in allowed if costs[tree.root][s] == score]
    possible = {tree.root: set(roots)}
    pairs = {}
    for n in tree.preorder:
        for c in tree.nodes[n].children:
            pairs[c] = set()
            for s in possible[n]:
                best = min(costs[c][t] + (s != t) for t in (0, 1))
                pairs[c].update((s, t) for t in (0, 1)
                                if costs[c][t] + (s != t) == best)
            possible[c] = {t for s, t in pairs[c]}
    bounds = {}
    for region in ['ALL'] + sorted(set(edge_regions.values())):
        for transition, key in [((0, 1), 'gains'), ((1, 0), 'losses')]:
            lo, hi = {}, {}
            for n in tree.postorder:
                lo[n], hi[n] = [0, 0], [0, 0]
                for s in (0, 1):
                    for c in tree.nodes[n].children:
                        best = min(costs[c][t] + (s != t) for t in (0, 1))
                        ts = [t for t in (0, 1) if costs[c][t] + (s != t) == best]
                        def inc(t):
                            return int((s, t) == transition and
                                       (region == 'ALL' or edge_regions[c] == region))
                        lo[n][s] += min(lo[c][t] + inc(t) for t in ts)
                        hi[n][s] += max(hi[c][t] + inc(t) for t in ts)
            bounds.setdefault(region, {})[key] = [min(lo[tree.root][s] for s in roots),
                                                  max(hi[tree.root][s] for s in roots)]
    return {'score': score, 'roots': roots, 'bounds': bounds}, possible, pairs


def table(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--input-dir",type=Path,default=BASE/"input")
    parser.add_argument("--outdir",type=Path,default=BASE/"results")
    args=parser.parse_args()
    inp,out=args.input_dir,args.outdir
    out.mkdir(parents=True,exist_ok=True)
    def rows(name):
        with (inp/name).open() as f:
            return list(csv.DictReader(f, delimiter='\t'))
    matrix, meta = rows('matrix.tsv'), rows('metadata.tsv')
    assert len({r['sample'] for r in matrix}) == len(matrix)
    assert len({r['sample'] for r in meta}) == len(meta)
    tree = read_tree(inp/'primary.treefile')
    state = {r['sample']: int(r['mrcA:p.Ala506Thr'])
             if r['mrcA:p.Ala506Thr'] in ('0','1') else None for r in matrix}
    lineage = {r['sample']:r['lineage'] for r in meta}
    assert set(tree.tips) == set(state) == set(lineage)
    assert set(lineage.values()) == {'Nichols','SS14'}
    conflicts = [r['sample'] for r in matrix if sum(r[k]=='1' for k in
                 ['mrcA:p.Ala506Thr','mrcA:p.Ala506Val','mrcA:p.Ala506Met']) > 1]
    assert not conflicts, conflicts
    counts, groups, descendants = {}, {}, {}
    for n in tree.postorder:
        node = tree.nodes[n]
        if node.children:
            descendants[n] = [s for c in node.children for s in descendants[c]]
        else:
            descendants[n] = [node.label]
        counts[n] = Counter(state[s] for s in descendants[n])
        groups[n] = set(lineage[s] for s in descendants[n])
    regions = {n: next(iter(groups[n])) if len(groups[n]) == 1 else 'mixed_backbone'
               for n in tree.nodes if n != tree.root}
    fits = {}
    for root in [None, 0, 1]:
        key = 'unrestricted' if root is None else f'root_{root}'
        fits[key], _, _ = fit(tree, state, regions, root)
    result, possible, pairs = fit(tree, state, regions)
    lineage_summary = []
    for group in ['Nichols','SS14']:
        samples = [s for s in tree.tips if lineage[s] == group]
        candidates = [n for n in tree.nodes if set(samples) <= set(descendants[n])]
        mrca = min(candidates, key=lambda n:len(descendants[n]))
        assert set(descendants[mrca]) == set(samples), 'Lineage is not monophyletic'
        count = Counter(state[s] for s in samples)
        # A carrier-compatible clade has no observed noncarrier tips; missing tips remain unknown.
        pure = [n for n in tree.nodes if groups[n] == {group} and counts[n][1] and not counts[n][0]]
        largest = max(pure, key=lambda n: counts[n][1])
        lineage_summary.append(dict(lineage=group, tips=len(samples), carriers=count[1],
             noncarriers=count[0], missing=count[None], callable=count[0]+count[1],
             carrier_percent=100*count[1]/(count[0]+count[1]), lineage_mrca=mrca,
             mrca_possible_states=','.join(map(str,sorted(possible[mrca]))),
             gains_min=result['bounds'][group]['gains'][0],
             gains_max=result['bounds'][group]['gains'][1],
             losses_min=result['bounds'][group]['losses'][0],
             losses_max=result['bounds'][group]['losses'][1],
             largest_carrier_compatible_clade=largest,
             largest_clade_carriers=counts[largest][1], largest_clade_missing=counts[largest][None],
             largest_clade_fraction_of_carriers=counts[largest][1]/count[1]))
    edge_rows=[]
    for n in tree.preorder:
        if n == tree.root:continue
        ps=pairs[n]
        edge_rows.append(dict(node=n,parent=tree.nodes[n].parent,lineage=regions[n],
            terminal=not tree.nodes[n].children,branch_length=tree.nodes[n].length,
            possible_transitions=';'.join(f'{a}>{b}' for a,b in sorted(ps)),
            gain_possible=(0,1) in ps,gain_forced=ps=={(0,1)},
            loss_possible=(1,0) in ps,loss_forced=ps=={(1,0)},
            descendant_carriers=counts[n][1],descendant_noncarriers=counts[n][0],
            descendant_missing=counts[n][None]))
    table(out/'lineage_summary.tsv',lineage_summary)
    table(out/'all_edges.tsv',edge_rows)
    table(out/'possible_gain_edges.tsv',[r for r in edge_rows if r['gain_possible']])
    table(out/'tip_states.tsv',[dict(sample=s,lineage=lineage[s],A506T=state[s] if state[s] is not None else 'NA') for s in tree.tips])
    document={'analysis':'Binary-state equal-cost parsimony on unchanged full rooted tree',
        'edge_assignment':'Lineage if all descendant tips belong to that lineage; includes lineage stem',
        'bounds':'Exact separate marginal extrema over globally optimal histories, not confidence intervals',
        'zero':'Not carrying A506T; may carry another residue at codon 506',
        'missing':'Unknown, not forced to zero', 'fits':fits,
        'lineage_summary':lineage_summary,
        'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inp.iterdir()},
        'analysis_completed':True}
    (out/'analysis.json').write_text(json.dumps(document,indent=2)+'\n')
    print(json.dumps(document,indent=2))

if __name__=='__main__':main()
