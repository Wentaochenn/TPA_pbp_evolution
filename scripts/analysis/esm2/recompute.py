"""Translate the annotated reference and compute its four-protein ESM-2 scan."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace

from esm2_stage import GENES, audit_scores, dump, fasta, run, sha, write


def load_extractor():
    path = Path(__file__).resolve().parents[2] / 'upstream/script/p5_01_extract_reference_pbp_proteins.py'
    spec = importlib.util.spec_from_file_location('pbp_reference_translation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def legacy_sequences(path):
    """Recover the original reference at every scored site without assuming initiation."""
    import csv
    residues = {g: {} for g in GENES}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            gene = row.get('protein') or row.get('gene')
            if gene not in GENES or row.get('scan_mode', 'full_dms') != 'full_dms':
                continue
            site = int(row.get('position') or row.get('site'))
            aa = row.get('wt_aa') or row.get('wt')
            if site in residues[gene] and residues[gene][site] != aa:
                raise ValueError('Previous scores have conflicting reference residues')
            residues[gene][site] = aa
    sequences = {}
    for gene, length in GENES.items():
        if set(residues[gene]) != set(range(1, length + 1)):
            raise ValueError('Previous score table has incomplete positions: ' + gene)
        sequences[gene] = ''.join(residues[gene][i] for i in range(1, length+1))
    return sequences


def compare_scores(previous, current, reference, out):
    old_reference = legacy_sequences(previous)
    old_summary, old, _ = audit_scores(previous, old_reference)
    new_summary, new, _ = audit_scores(current, reference)
    if not old_summary['coverage_complete'] or not new_summary['coverage_complete']:
        raise ValueError('A compared full substitution scan is incomplete')
    differences, summaries = [], []
    for gene in GENES:
        altered = [i for i, (a, b) in enumerate(zip(old_reference[gene], reference[gene]), 1) if a != b]
        keys = sorted(k for k in old.keys() | new.keys() if k[0] == gene)
        deltas = []
        for key in keys:
            a, b = old.get(key), new.get(key)
            # The alternatives at a changed reference site belong to different comparisons.
            comparable = a is not None and b is not None and a['wt_aa'] == b['wt_aa']
            delta = b['llr'] - a['llr'] if comparable else None
            if comparable:
                deltas.append(delta)
            differences.append(dict(gene=gene, position=key[1], alternate=key[2],
                previous_reference=a['wt_aa'] if a else '', current_reference=b['wt_aa'] if b else '',
                previous_llr=a['llr'] if a else '', current_llr=b['llr'] if b else '',
                same_reference_comparison=comparable, delta_llr=delta if delta is not None else ''))
        unchanged = not altered
        summaries.append(dict(gene=gene, reference_changed_positions=altered,
            compared_substitutions=len(deltas), changed_scores_above_1e_6=sum(abs(x)>1e-6 for x in deltas),
            max_absolute_difference=max(map(abs,deltas),default=None),
            unchanged_reference_control=unchanged,
            unchanged_control_pass=(all(abs(x)<=1e-6 for x in deltas) if unchanged else None)))
    write(out/'score_comparison.tsv', differences, list(differences[0]))
    result = dict(state='PASS' if all(r['unchanged_control_pass'] is not False for r in summaries) else 'REVIEW_REQUIRED',
        previous_scores_sha256=sha(previous), new_scores_sha256=sha(current), absolute_tolerance=1e-6,
        genes=summaries, reference_changed_site_scores_not_directly_compared=True,
        historical_inputs_modified=False, release_adoption='pending review of regenerated figures and tables')
    dump(out/'SCORE_COMPARISON.json', result)
    return result


def recompute(args):
    root, out = args.project.resolve(), args.out.resolve()
    if out.exists():
        raise ValueError('Use a new recomputation output directory')
    previous = args.previous_scores.resolve() if args.previous_scores else None
    old_ref = None
    if previous is not None:
        if not previous.is_file():
            raise ValueError('Previous complete four-protein score table is missing')
        old_ref = legacy_sequences(previous)
        if not audit_scores(previous, old_ref)[0]['coverage_complete']:
            raise ValueError('Previous four-protein scan is incomplete')
    out.mkdir(parents=True)
    try:
        dump(out/'RUN_STATUS.json', dict(state='running', stage='reference_translation'))
        reference = args.reference.resolve() if args.reference else root/'refs/TPA_Nichols.fasta'
        genbank = args.genbank.resolve() if args.genbank else root/'refs/TPA_Nichols.gbff'
        proteins = out/'reference_proteins.fasta'
        load_extractor().extract(reference, genbank, proteins, out/'REFERENCE_AUDIT.json')
        if old_ref is not None:
            for gene, sequence in fasta(proteins).items():
                changes = [i for i, (a,b) in enumerate(zip(old_ref[gene],sequence),1) if a!=b]
                if any(i != 1 for i in changes):
                    raise ValueError('Unexpected reference change outside initiation residue: '+gene)
        dump(out/'RUN_STATUS.json',dict(state='running',stage='ESM2_inference'))
        run(SimpleNamespace(action='run', protein_fasta=proteins, out=out/'inference',
                            model=args.model, device=args.device, frozen_release=None))
        comparison = compare_scores(previous, out/'inference/esm2_all_substitutions.tsv', fasta(proteins), out) if previous else None
        dump(out/'RUN_STATUS.json', dict(state='complete' if comparison is None or comparison['state']=='PASS' else 'review_required',
            stage='complete',exit_code=0, score_comparison=comparison['state'] if comparison else 'not_requested',
            score_comparison_requested=previous is not None,
            expected_scores=49286, reference_fasta_sha256=sha(proteins), historical_inputs_modified=False))
        print('RECOMPUTATION_OUTPUT:',out,flush=True)
    except Exception as error:
        dump(out/'RUN_STATUS.json',dict(state='failed',exit_code=1,error=str(error)))
        raise
