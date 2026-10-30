"""Bind a completed annotated-reference ESM scan into new figure and table outputs."""
import csv
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
from integrity import sha


def completed_scan(inference):
    """Validate the reference and optional comparison recorded by recompute."""
    inference = Path(inference).resolve()
    required = [inference/'RUN_STATUS.json', inference/'INPUT_CONTEXT.json',
                inference/'MODEL_CONTEXT.json', inference/'ESM2_COVERAGE.json',
                inference/'esm2_all_substitutions.tsv', inference.parent/'RUN_STATUS.json',
                inference.parent/'REFERENCE_AUDIT.json', inference.parent/'reference_proteins.fasta']
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError('Use a completed annotated-reference recomputation; missing: '+', '.join(missing))
    status = json.loads((inference/'RUN_STATUS.json').read_text())
    if status.get('state') != 'complete' or status.get('exit_code') != 0:
        raise ValueError('ESM-2 inference is not complete')
    parent_status=json.loads((inference.parent/'RUN_STATUS.json').read_text())
    if parent_status.get('state')!='complete' or parent_status.get('exit_code')!=0:
        raise ValueError('Annotated-reference recomputation is not complete')
    comparison = parent_status.get('score_comparison')
    if comparison not in ('PASS', 'not_requested'):
        raise ValueError('Requested historical-score comparison did not pass')
    requested = parent_status.get('score_comparison_requested', comparison!='not_requested')
    if not isinstance(requested, bool) or requested != (comparison=='PASS'):
        raise ValueError('Historical-score comparison status is inconsistent')
    audit_path = inference.parent/'REFERENCE_AUDIT.json'
    reference_path = inference.parent/'reference_proteins.fasta'
    audit = json.loads(audit_path.read_text())
    if audit.get('state') != 'PASS' or audit['protein_fasta_sha256'] != sha(reference_path):
        raise ValueError('Annotated reference translation audit is missing or mismatched')
    context = json.loads((inference/'INPUT_CONTEXT.json').read_text())
    if context['protein_fasta_sha256'] != sha(reference_path):
        raise ValueError('Scored reference FASTA differs from translation audit')
    if parent_status.get('reference_fasta_sha256') != sha(reference_path):
        raise ValueError('Recomputation reference identity is missing or mismatched')
    if requested:
        report_path = inference.parent/'SCORE_COMPARISON.json'
        if not report_path.is_file():
            raise ValueError('Requested historical-score comparison report is missing')
        report = json.loads(report_path.read_text())
        if report.get('state')!='PASS' or report.get('new_scores_sha256')!=sha(inference/'esm2_all_substitutions.tsv'):
            raise ValueError('Historical-score comparison report does not match the completed scores')
    return reference_path, audit_path, comparison


def refresh(run, inference):
    run, inference = Path(run).resolve(), Path(inference).resolve()
    reference_path, audit_path, comparison = completed_scan(inference)
    module_path = Path(__file__).resolve().parent/'analysis/esm2/esm2_stage.py'
    spec = importlib.util.spec_from_file_location('esm2_binding_validation',module_path)
    stage = importlib.util.module_from_spec(spec); spec.loader.exec_module(stage)
    sequences = stage.fasta(reference_path)
    scores = inference/'esm2_all_substitutions.tsv'
    coverage, _, _ = stage.audit_scores(scores,sequences)
    if not coverage['coverage_complete']:
        raise ValueError('ESM-2 substitution scan is incomplete')
    source_path = run/'SOURCE_MANIFEST.json'
    sources = json.loads(source_path.read_text())
    matrix = Path(sources['matrix.tsv']['path'])
    if not matrix.is_absolute(): matrix = (run/matrix).resolve()
    if sha(matrix) != sources['matrix.tsv']['sha256']:
        raise ValueError('Callability matrix differs from figure input manifest')
    data = run/'work/figures/data'
    before = {p.name:sha(p) for p in data.iterdir() if p.is_file()}
    from prepare import esm
    esm(SimpleNamespace(used={'ESM':scores,'matrix.tsv':matrix}),data)
    allowed = {'esm_site_means.tsv','esm_specific.tsv','esm_observed.tsv'}
    after = {p.name:sha(p) for p in data.iterdir() if p.is_file()}
    if any(before[k] != after.get(k) for k in before if k not in allowed):
        raise ValueError('An unrelated prepared figure input changed')
    sources['ESM'] = dict(path=str(scores),sha256=sha(scores),bytes=scores.stat().st_size)
    source_path.write_text(json.dumps(sources,indent=2)+'\n')
    (run/'work/INPUT_MANIFEST.json').write_text(json.dumps(after,indent=2)+'\n')
    resolution = run/'INPUT_RESOLUTION.tsv'
    if resolution.exists():
        with resolution.open(newline='') as f:
            reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
        for row in rows:
            if row.get('input')=='ESM': row.update(path=str(scores),status='FOUND')
        with resolution.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
    import shutil
    provenance=run/'ESM2_Inputs';provenance.mkdir()
    for source in [reference_path,audit_path,inference/'INPUT_CONTEXT.json',inference/'MODEL_CONTEXT.json',inference/'ESM2_COVERAGE.json']:
        shutil.copy2(source,provenance/source.name)
    result=dict(state='PASS',scores_sha256=sha(scores),reference_sha256=sha(reference_path),
                historical_score_comparison=comparison,
                expected_scores=49286,changed_prepared_inputs=[k for k in sorted(allowed) if before.get(k)!=after[k]],
                source_data_and_formal_table='regenerated by the subsequent R export')
    (run/'ESM2_INPUT_UPDATE.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
