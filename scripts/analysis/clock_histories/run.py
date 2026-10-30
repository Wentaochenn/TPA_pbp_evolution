"""Summarize existing strict/relaxed chains and reconstruct matched acquisition dates."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
from prepare import read, table, sha, sample_posterior, audit_chain
from treeannotator_support import command, require_tree
from compare import compare_events

HERE=Path(__file__).resolve().parent


def execute(args,log,env):
    print('[RUN]',log.name,flush=True)
    log.with_suffix('.command.json').write_text(json.dumps([str(x) for x in args],indent=2)+'\n')
    with log.open('w') as handle:
        subprocess.run([str(x) for x in args],stdout=handle,stderr=subprocess.STDOUT,env=env,check=True)


def extract_ess(report, target):
    header = None
    row = None
    for line in report.read_text().splitlines():
        tokens = line.split()
        if tokens[:2] == ['sample','filename']:
            header = tokens
        elif header and len(tokens)==len(header) and tokens[0].isdigit():
            row=tokens
    if header is None or row is None:
        raise ValueError('No complete ESS row; inspect '+str(report))
    records=[]
    for i,key in enumerate(header):
        if key.endswith('.ESS'):
            try:value=float(row[i])
            except ValueError:value=float('nan')
            records.append({'parameter':key[:-4],'ESS':value,
                            'status':'PASS' if value>=200 else 'CHECK'})
    if not records:raise ValueError('No ESS columns')
    table(target,records)
    if any(r['status']!='PASS' for r in records):
        raise ValueError('ESS gate not met; inspect '+str(target))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True);p.add_argument('--rscript',default='Rscript')
    p.add_argument('--treeannotator',default='treeannotator');p.add_argument('--loganalyser',default='loganalyser')
    p.add_argument('--java-heap',default='64g')
    a=p.parse_args();root=a.root.resolve();cfg=json.loads(a.config.read_text());out=a.outdir.resolve()
    if out.exists():raise ValueError('Use a new clock-history output directory')
    def path(value):
        q=Path(value);return q.resolve() if q.is_absolute() else (root/q).resolve()
    matrix=path(cfg['matrix']);dates=path(cfg['dates'])
    chains={m:{k:path(cfg[m][k]) for k in ['xml','trace','posterior','target_mcc']} for m in ['strict','relaxed']}
    inputs=[matrix,dates,*[q for model in chains.values() for q in model.values()]]
    for q in inputs:
        if not q.is_file() or not q.stat().st_size:raise ValueError('Missing clock input: '+str(q))
    before={str(q):sha(q) for q in inputs};out.mkdir(parents=True)
    for folder in ['audit','logs','sampled','clock_sensitivity/strict','clock_sensitivity/relaxed']:(out/folder).mkdir(parents=True)
    state={'state':'running','stage':'start','new_MCMC_run':False,'started':datetime.now(timezone.utc).isoformat()}
    env=dict(os.environ);env.pop('DISPLAY',None)
    env['_JAVA_OPTIONS']='-Xmx'+a.java_heap
    env['JAVA_TOOL_OPTIONS']=env.get('JAVA_TOOL_OPTIONS','')+' -Djava.awt.headless=true'
    (out/'INPUT_CONTEXT.json').write_text(json.dumps({'inputs':before,'config':cfg,'posterior_trees':500,
         'posterior_burnin':0.10,'posterior_thin':10,'maps_per_tree':1,'gain_root_prior':[0.999,0.001]},indent=2)+'\n')
    try:
        for model,source in chains.items():
            state['stage']=model+'_ESS';report=out/'audit'/(model+'_loganalyser.txt')
            execute([a.loganalyser,'-oneline','-b','10',source['trace']],report,env)
            extract_ess(report,out/'audit'/(model+'_ESS.tsv'))
            state['stage']=model+'_sample_posterior';sample=out/'sampled'/(model+'.trees')
            info=sample_posterior(source['posterior'],sample)
            audit_chain(source['trace'],info,source['xml'],out/'audit'/(model+'_chain_audit.json'))
            folder=out/'clock_sensitivity'/model;ca=folder/('mcc_'+model+'_CA.tree')
            state['stage']=model+'_CA_MCC';log=out/'logs'/('treeannotator_'+model+'.log')
            execute(command(a.treeannotator,source['posterior'],source['target_mcc'],ca),log,env)
            require_tree(ca,log)
            execute([a.rscript,HERE/'check_tree.R',ca,dates,folder/'chronogram_check.tsv',source['target_mcc']],out/'logs'/('chronogram_'+model+'.log'),env)
            mapping_env=dict(os.environ,PBP_MK_HELPER=str(HERE.parents[1]/'upstream/script/pbp_mk_missing_helpers.R'),
                PBP_MATRIX_FILE=str(matrix),PBP_DATES_FILE=str(dates),PBP_MCC_TREE_FILE=str(ca),
                PBP_POSTERIOR_TREES_FILE=str(sample),PBP_OUTDIR=str(folder),PBP_POSTERIOR_BURNIN='0',
                PBP_POSTERIOR_THIN='1',PBP_SIMMAPS_PER_TREE='1',PBP_MAX_POSTERIOR_TREES='500',
                PBP_MATCH_JACCARD='0.50',PBP_MIN_MATCH_SUPPORT='0.50',PBP_STOP_IF_NO_POSTERIOR='TRUE')
            state['stage']=model+'_PBP_dates'
            execute([a.rscript,HERE/'payload/date_origins.R'],out/'logs'/('dates_'+model+'.log'),mapping_env)
            maps=read(folder/'pbp_mapping_execution_audit.tsv');counts={}
            for row in maps:counts[row['locus']]=counts.get(row['locus'],0)+1
            if len(counts)!=6 or any(n!=500 for n in counts.values()):raise ValueError('Incomplete posterior-tree/locus accounting')
        state['stage']='event_matching';compare_events(out)
        if any(sha(Path(q))!=h for q,h in before.items()):raise ValueError('Original chain input changed during summarization')
        state.update(state='complete',stage='complete',exit_code=0)
    except Exception as error:
        state.update(state='failed',error=type(error).__name__+': '+str(error),exit_code=1);raise
    finally:(out/'RUN_STATUS.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__':main()
