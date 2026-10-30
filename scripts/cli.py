"""Reproduce study figures or run the included analysis modules."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MODULES = {
    'clustered-dates': 'clustered_dates/temporal.py',
    'selection-topologies': 'selection_topologies/ensemble.py',
    'protein-types': 'protein_types/run.py',
    'clock-histories': 'clock_histories/run.py',
}


def contained(path, root):
    return path == root or root in path.parents


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    for name in ['reproduce', 'publish']:
        q = sub.add_parser(name, help='Generate figures and tables from ' +
                           ('the Zenodo dataset' if name=='reproduce' else 'completed analysis results'))
        q.add_argument('--project', type=Path, default=Path.cwd())
        q.add_argument('--out', type=Path, required=name=='reproduce')
        if name=='reproduce':
            q.add_argument('--data', type=Path, required=True)
            q.add_argument('--compare', action='store_true', help='Compare outputs with the deposited figures and tables')
            q.add_argument('--rerender-structure', action='store_true', help='Re-render the PDB instead of using its deposited rendering')
        q.add_argument('--inputs', type=Path, help='JSON mapping of analytical input files')
        q.add_argument('--esm2-results', type=Path, help='Completed annotated-reference ESM recomputation inference directory')
        q.add_argument('--pymol-python', type=Path)
        q.add_argument('--structure-render', type=Path)
    for name in ['upstream', 'esm2', 'inputs']:
        sub.add_parser(name).add_argument('arguments', nargs=argparse.REMAINDER)
    q = sub.add_parser('analysis')
    q.add_argument('module', choices=MODULES)
    q.add_argument('arguments', nargs=argparse.REMAINDER)
    q = sub.add_parser('verify-data')
    q.add_argument('--data', type=Path, required=True)
    q = sub.add_parser('compare')
    q.add_argument('--reference', type=Path, required=True)
    q.add_argument('--candidate', type=Path, required=True)
    a = p.parse_args(argv)
    if a.action=='verify-data':
        from integrity import verify_data
        print(json.dumps(verify_data(a.data), indent=2))
        return
    if a.action=='compare':
        subprocess.run([sys.executable, str(HERE/'compare_results.py'), '--reference',
                        str(a.reference), '--candidate', str(a.candidate)], check=True)
        return
    if a.action in ['reproduce', 'publish']:
        if a.action=='reproduce' and a.compare and a.esm2_results:
            raise ValueError('Use --compare only when reproducing the same score release; omit it for corrected ESM inputs')
        project = a.project.resolve()
        if a.out:
            out = a.out.resolve()
            if contained(out, HERE.parent):
                raise ValueError('Choose an output directory outside the code directory')
            if out.exists():
                raise ValueError('Use a new output directory: ' + str(out))
        cmd = [sys.executable, str(HERE/'run.py'), '--project', str(project)]
        if a.out:
            cmd += ['--out', str(a.out.resolve())]
        render = a.structure_render
        if a.action=='reproduce':
            data = a.data.resolve()
            if contained(out, data) or contained(data, out):
                raise ValueError('Keep the output directory separate from the deposited data')
            from integrity import verify_data
            print('[VERIFY]', json.dumps(verify_data(data)), flush=True)
            cmd += ['--reuse-run', str(data/'frozen_run')]
            if not render and not a.rerender_structure and (data/'Structure_Render').is_dir():
                render = data/'Structure_Render'
        for flag, value in [('inputs', a.inputs), ('esm2-results', a.esm2_results), ('pymol-python', a.pymol_python), ('structure-render', render)]:
            if value:
                cmd += ['--'+flag, str(value.resolve())]
    else:
        script = HERE/('upstream.py' if a.action=='upstream' else
                       'study_inputs.py' if a.action=='inputs' else
                       'analysis/esm2/esm2_stage.py' if a.action=='esm2' else 'analysis/'+MODULES[a.module])
        arguments = a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
        cmd = [sys.executable, str(script)] + arguments
    subprocess.run(cmd, check=True)
    if a.action=='reproduce' and a.compare:
        subprocess.run([sys.executable, str(HERE/'compare_results.py'), '--reference',
                        str(data), '--candidate', str(out)], check=True)


if __name__=='__main__':
    main()
