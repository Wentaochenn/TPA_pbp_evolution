#!/usr/bin/env python3
"""Export manuscript figures and tables directly from an existing cluster project."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from prepare import prepare

PACKAGE=Path(__file__).resolve().parent


def execute(args,log):
    with log.open('w') as f:
        result=subprocess.run([str(x) for x in args],stdout=f,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError('Command failed; inspect '+str(log))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path.cwd())
    p.add_argument('--out',type=Path)
    p.add_argument('--extra-root',type=Path,action='append',default=[])
    p.add_argument('--inputs',type=Path,help='JSON overrides for named source paths')
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--esm2-results',type=Path,help='Completed annotated-reference ESM inference directory')
    p.add_argument('--reuse-run',type=Path,help='Read prepared inputs from a previous successful figure run')
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--pymol-python',type=Path)
    p.add_argument('--pymol-environment',default='pbp_structure_render')
    p.add_argument('--structure-render',type=Path,help='Reuse a completed render with matching model, script, camera and resolution')
    p.add_argument('--structure-width',type=int,default=2400)
    a=p.parse_args();project=a.project.resolve()
    stamp=time.strftime('%Y%m%d_%H%M%S')
    out=a.out.resolve() if a.out else project/'results/figures'/f'run_{os.environ.get("JOB_ID","manual")}_{stamp}'
    if out.exists():raise ValueError('Use a new output directory')
    status={'state':'preparing','project':str(project),'output':str(out),'scientific_models_rerun':False}
    try:
        if a.reuse_run:
            previous=a.reuse_run.resolve()
            if not (previous/'work/INPUT_MANIFEST.json').is_file():
                raise ValueError('Previous run has no prepared input manifest')
            out.mkdir(parents=True)
            shutil.copytree(previous/'work',out/'work')
            for name in ('figures/R','figures/qa'):
                shutil.rmtree(out/'work'/name,ignore_errors=True)
                shutil.copytree(PACKAGE/name,out/'work'/name)
            for name in ('SOURCE_MANIFEST.json','INPUT_RESOLUTION.tsv','PREPARE_STATUS.json'):
                shutil.copyfile(previous/name,out/name)
            source_manifest=json.loads((out/'SOURCE_MANIFEST.json').read_text())
            for record in source_manifest.values():
                if not Path(record['path']).is_absolute():
                    record['path']=str((previous/record['path']).resolve())
            (out/'SOURCE_MANIFEST.json').write_text(json.dumps(source_manifest,indent=2)+'\n')
            shutil.copytree(previous/'lineage_history',out/'lineage_history')
            shutil.copyfile(out/'lineage_history/lineage_summary.tsv',out/'work/figures/data/lineage_summary.tsv')
            from prepare import sha
            data=out/'work/figures/data'
            # Verify original prepared bytes before adding the lineage source-data component.
            frozen=json.loads((previous/'work/INPUT_MANIFEST.json').read_text())
            if not all(sha(data/name)==value for name,value in frozen.items()):
                raise ValueError('Previous prepared input checksums differ')
            (out/'work/INPUT_MANIFEST.json').write_text(json.dumps({p.name:sha(p) for p in data.iterdir() if p.is_file()},indent=2)+'\n')
            status['prepared_inputs_reused_from']=str(previous)
        else:
            prepare(project,out,a.extra_root,a.inputs)
        from publication_selection import project as project_selection
        status['selection_input_scope']=project_selection(out)
        if a.esm2_results:
            from esm2_inputs import refresh
            status['esm2_input_update']=refresh(out,a.esm2_results)
        if not a.prepare_only:
            print('[RUN] Publication runtime and R syntax preflight',flush=True)
            execute([sys.executable,PACKAGE/'preflight.py','--package',PACKAGE,
                     '--rscript',a.rscript,'--out',out/'ENVIRONMENT_CHECK.json'],
                    out/'R_syntax_preflight.log')
            from input_bindings import validate_model_mappings
            status['model_mapping_validation']=validate_model_mappings(out)
            from structure_stage import ensure_structure
            status['structure_render'] = ensure_structure(out,a.pymol_python,a.pymol_environment,a.structure_render,a.structure_width)
            from prepare import sha
            data=out/'work/figures/data'
            from prepare import sasa
            sasa(data/'mrcA_model.pdb',data)
            (out/'work/INPUT_MANIFEST.json').write_text(json.dumps({p.name:sha(p) for p in data.iterdir() if p.is_file()},indent=2)+'\n')
        (out/'RUN_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
        if a.prepare_only:
            status.update(state='prepared',figures_rendered=False);return
        diagnostics = ['TableS12A_clustered_date_tests.tsv',
                       'TableS12B_subset_stability.tsv',
                       'TableS13B_strict_ESS.tsv', 'TableS13C_relaxed_ESS.tsv']
        missing = [name for name in diagnostics
                   if not (out/'work/figures/data'/name).is_file()]
        if missing:
            raise ValueError('Retained diagnostic source files are missing: '+', '.join(missing))
        rscript=shutil.which(a.rscript)
        if not rscript:raise RuntimeError('Rscript not found; activate the figure_environment.yml environment')
        status['state']='rendering'
        (out/'RUN_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
        print('[RUN] R figures; progress in R_render.console.log',flush=True)
        execute([rscript,PACKAGE/'run_figures.R',out/'work',out/'R_export'],out/'R_render.console.log')
        print('[RUN] PDF text and collision checks',flush=True)
        execute([sys.executable,PACKAGE/'figures/qa/audit_all.py',out/'R_export'],out/'R_render_QA.console.log')
        status['state']='exporting_tables'
        (out/'RUN_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
        print('[RUN] supplementary tables',flush=True)
        execute([rscript,PACKAGE/'export_tables.R',out,PACKAGE],out/'table_export.console.log')
        for language in ('EN','ZH'):
            shutil.copyfile(PACKAGE/f'Figure_legends_and_table_list_{language}.md',
                            out/f'Figure_legends_and_table_list_{language}.md')
        shutil.copytree(out/'R_export/figures',out/'Figures')
        shutil.copytree(out/'R_export/panel_data',out/'Figure_Source_Data')
        shutil.copytree(out/'R_export/rendered_layers',out/'Figure_Source_Data/rendered_layers')
        shutil.copyfile(PACKAGE/'SOURCE_DATA_README.md',out/'Figure_Source_Data/README.md')
        shutil.copyfile(out/'R_export/FIGURE_INDEX.json',out/'FIGURE_INDEX.json')
        shutil.copytree(out/'R_export/qa',out/'QA')
        exported=json.loads((out/'FIGURE_INDEX.json').read_text())
        if len(exported)!=9:raise ValueError('Expected five main and four supplementary figures')
        table_status=json.loads((out/'TABLE_EXPORT_STATUS.json').read_text())
        from check_outputs import validate
        validation=validate(out)
        status.update(state='complete_with_optional_table_gaps' if table_status['optional_table_gaps'] else 'complete',
                      figures=9,panels=validation['panels'],output_validation='PASS',
                      figure_source_data='exported for every panel',**table_status,exit_code=0)
    except Exception as error:
        status.update(state='failed',error=f'{type(error).__name__}: {error}',exit_code=1)
        raise
    finally:
        if out.exists():(out/'RUN_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')
        print('OUTPUT:',out,flush=True)


if __name__=='__main__':main()
