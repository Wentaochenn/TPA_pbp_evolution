"""Genome-tree ensemble sensitivity. Original alignments; descriptive stability only."""
import argparse, collections, copy, fcntl, gzip, io, json, math, os, random, re
import shutil, subprocess, sys, time, traceback
from contextlib import contextmanager
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.common import alignment, check_known_hashes, context, dump, fasta, fingerprint, prune, read_tsv, resolve, sha, Tasks, tsv
from analysis.selection_topologies.fubar_sites import json_sites

PACKAGE=Path(__file__).resolve().parent

def load(p):return json.loads(Path(p).read_text())

def names(tr):
    n=[t.name for t in tr.get_terminals()]
    if not n or any(not x for x in n) or len(set(n))!=len(n):raise ValueError('Missing or duplicate tree tip IDs')
    return set(n)

def splits(tr):
    ids=sorted(names(tr));idx={n:i for i,n in enumerate(ids)};full=(1<<len(ids))-1;found=set()
    def rec(c):
        mask=1<<idx[c.name] if c.is_terminal() else 0
        for child in c.clades:mask|=rec(child)
        size=mask.bit_count()
        if 1<size<len(ids)-1:found.add(min(mask,full^mask))
        return mask
    rec(tr.root)
    return found

def topo_hash(tr):
    import hashlib
    text=json.dumps([sorted(names(tr)),sorted(splits(tr))],separators=(',',':'))
    return hashlib.sha256(text.encode()).hexdigest()

def topology_text(tr):
    def rec(c):
        if c.is_terminal():
            if not re.fullmatch(r'S[0-9]+',c.name):raise ValueError('Tree output requires safe tip IDs')
            return c.name
        return '('+','.join(rec(x) for x in c.clades)+')'
    return rec(tr.root)+';\n'

def prepared_tree(tr,ids,aliases):
    t=copy.deepcopy(tr)
    for c in t.find_clades():c.branch_length=0.0 if c.branch_length is None else c.branch_length
    t=prune(t,ids)
    while len(t.root.clades)==1:t.root=t.root.clades[0]
    for c in t.find_clades():
        if c.is_terminal():c.name=aliases[c.name]
        else:c.name=None
        c.confidence=None;c.comment=None
    return t

def iter_trees(path):
    from Bio import Phylo
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:yield from Phylo.parse(f,'newick')

def version(exe):
    r=subprocess.run([exe,'--version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60)
    if not r.stdout.strip():raise ValueError('Empty HyPhy version output')
    return r.stdout.strip()

def hyphy_cmd(exe,ver,threads,seed,method,aln,tree,out):
    match=re.search(r'(\d+)\.(\d+)\.(\d+)',ver);cmd=[exe,f'CPU={threads}']
    if match and tuple(map(int,match.groups()))>=(2,5,101):cmd.append(f'ENV=RANDOM_SEED={seed}')
    return cmd+[method,'--alignment',str(aln),'--tree',str(tree),'--output',str(out),'--code','Universal']

@contextmanager
def lock(path):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,'a') as f:
        try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('Another process holds '+str(path))
        yield

def candidate_source(path,cfg,root,explicit):
    s=str(path);s=s[:-3] if s.endswith('.gz') else s
    prefix=re.sub(r'\.(ufboot|boottrees|bootstrap|trees)$','',s)
    ml=resolve(root,cfg['source_ml_tree']) if explicit and cfg.get('source_ml_tree') else Path(prefix+'.treefile')
    report=resolve(root,cfg['provenance_report']) if explicit and cfg.get('provenance_report') else Path(prefix+'.iqtree')
    if not report.exists() and not (explicit and cfg.get('provenance_report')):report=Path(prefix+'.log')
    return ml,report

def inventory(root,cfg,out,primary):
    from Bio import Phylo
    explicit=bool(cfg.get('bootstrap_trees'))
    files=[resolve(root,cfg['bootstrap_trees'])] if explicit else sorted({p for pat in ['*.ufboot','*.ufboot.gz','*.boottrees','*.boottrees.gz'] for p in resolve(root,cfg['bootstrap_search_dir']).rglob(pat)})
    rows=[];valid=[];tipset=names(primary)
    for p in files:
        ml,report=candidate_source(p,cfg,root,explicit)
        row=dict(bootstrap_trees=str(p),source_ml_tree=str(ml),provenance_report=str(report),status='BLOCKED',reason='')
        try:
            if not p.is_file() or not ml.is_file() or not report.is_file():
                raise ValueError('Bootstrap trees and associated ML tree/report are required')
            source=Phylo.read(str(ml),'newick')
            if names(source)!=tipset:raise ValueError('Source ML tips differ from the primary tree')
            text=report.read_text(errors='replace')
            if not re.search(r'IQ-TREE',text,re.I) or not re.search(r'bootstrap|ufboot',text,re.I):
                raise ValueError('The report must document a genome-wide IQ-TREE bootstrap analysis')
            n=0
            for tr in iter_trees(p):
                n+=1
                if names(tr)!=tipset:raise ValueError('Bootstrap taxon set differs at replicate '+str(n))
            if n<cfg['n_bootstrap_trees']:raise ValueError('Insufficient bootstrap replicates')
            row.update(status='ELIGIBLE',n_trees=n,n_tips=len(tipset))
            valid.append((p,ml,report,n))
        except Exception as error:row['reason']=str(error)
        rows.append(row)
    tsv(out/'BOOTSTRAP_CANDIDATES.tsv',rows,['bootstrap_trees','source_ml_tree','provenance_report','status','reason','n_trees','n_tips'])
    if len(valid)!=1:
        raise ValueError('Specify one genome-wide bootstrap source and its associated ML tree and report in the configuration')
    return valid[0]


def check_cfg(cfg):
    for key in ['n_bootstrap_trees','threads','max_parallel']:
        if type(cfg[key]) is not int or cfg[key]<1:raise ValueError('Invalid '+key)
    if sum(cfg['genes'].values())!=2594:raise ValueError('This study uses 2,594 codons across four genes')
    if cfg['genes']!={'mrcA':884,'pbp1':624,'pbp2':652,'tp47':434}:raise ValueError('Gene scope changed; requires a new analysis design')
    if cfg['bootstrap_kind'] not in ['auto','ultrafast','standard']:raise ValueError('Invalid bootstrap_kind')
    if cfg['fubar_posterior_threshold']!=0.9:raise ValueError('Predeclared posterior threshold changed')

def validate_manifest(run):
    import Bio
    ctx=load(run/'CONTEXT.json')
    for category in ['inputs','code','executables','prepared']:
        for p,h in ctx[category].items():
            if not Path(p).is_file() or sha(p)!=h:raise ValueError('Immutable run input changed: '+p)
    if Bio.__version__!=ctx['biopython']:raise ValueError('Biopython version changed')
    if version(ctx['hyphy'])!=ctx['hyphy_version']:raise ValueError('HyPhy version changed')
    return ctx


def prepare(root,run):
    import Bio
    from Bio import Phylo
    cfg=load(run/'config.json');check_cfg(cfg)
    if (run/'PREPARE_STATUS.json').exists() and load(run/'PREPARE_STATUS.json').get('state')=='complete':
        validate_manifest(run);print('REUSE preparation',flush=True);return
    dump(run/'PREPARE_STATUS.json',{'state':'running'})
    try:
        primary_path=resolve(root,cfg['primary_tree']);primary=Phylo.read(str(primary_path),'newick')
        p,ml,report,total=inventory(root,cfg,run,primary)
        hy=shutil.which(cfg['hyphy'])
        if not hy:raise ValueError('HyPhy not found in active environment')
        hy=str(Path(hy).resolve());ver=version(hy)
        inputs=[primary_path,p,ml,report,run/'config.json']+[resolve(root,cfg['selection_dir'])/(g+'_codon.fasta') for g in cfg['genes']]
        baseline_paths={g:{m:str(resolve(root,cfg['selection_dir'])/(g+'_'+m+'.json')) for m in ['fubar']} for g in cfg['genes']}
        for gene,paths in baseline_paths.items():
            for method,path in paths.items():
                json_sites(path,method,cfg['genes'][gene]);inputs.append(Path(path))
        check_known_hashes(root,cfg,inputs)
        original_hashes=fingerprint(inputs)
        base_context={'inputs':original_hashes,'code':fingerprint(list(PACKAGE.glob('*.py'))+list(PACKAGE.glob('*.sh'))),'executables':fingerprint([hy]),'hyphy':hy,'hyphy_version':ver,'biopython':Bio.__version__,'baseline_paths':baseline_paths}
        if (run/'CONTEXT.json').exists():
            old=load(run/'CONTEXT.json')
            if any(old[k]!=v for k,v in base_context.items()):raise ValueError('Run context changed; start a new run')
        sample_indices=sorted(random.Random(cfg['tree_selection_seed']).sample(range(1,total+1),cfg['n_bootstrap_trees']))
        selected={};wanted=set(sample_indices)
        for idx,tr in enumerate(iter_trees(p),1):
            if idx in wanted:selected[idx]=tr
        plan=[{'task':1,'tree_id':'primary','source_index':0,'role':'control'}]+[{'task':j+2,'tree_id':f'bootstrap_{j+1:03d}','source_index':i,'role':'bootstrap'} for j,i in enumerate(sample_indices)]
        # Membership is frozen before any selection likelihood is evaluated.
        dump(run/'SAMPLING_PLAN.json',{'seed':cfg['tree_selection_seed'],'source_trees':total,'sampled_indices':sample_indices,'duplicates':'retain replicate multiplicity; no topology-based selection','plan':plan})
        prepared=[];genes={};stats=[]
        for gene,ncodons in cfg['genes'].items():
            seqs=alignment(resolve(root,cfg['selection_dir'])/(gene+'_codon.fasta'),ncodons)
            if not set(seqs)<=names(primary):raise ValueError('Gene samples absent from primary tree')
            if any(set(s)<=set('-') for s in seqs.values()):raise ValueError('All-missing sample in '+gene)
            aliases={s:f'S{i+1:05d}' for i,s in enumerate(sorted(seqs))};folder=run/'prepared'/gene;folder.mkdir(parents=True,exist_ok=True)
            a=folder/'alignment.fasta';fasta(a,{aliases[k]:v for k,v in seqs.items()});prepared.append(a)
            cross=folder/'TIP_CROSSWALK.tsv';tsv(cross,[{'sample':k,'safe_id':v} for k,v in aliases.items()]);prepared.append(cross)
            genes[gene]={'ncodons':ncodons,'ntips':len(seqs),'alignment':str(a),'trees':{}}
            for entry in plan:
                tr=primary if entry['role']=='control' else selected[entry['source_index']]
                t=prepared_tree(tr,seqs,aliases);f=folder/(entry['tree_id']+'.nwk');f.write_text(topology_text(t));prepared.append(f);genes[gene]['trees'][entry['tree_id']]=str(f)
        ps=splits(primary)
        for entry in plan:
            tr=primary if entry['role']=='control' else selected[entry['source_index']];sp=splits(tr);rf=len(ps^sp);den=len(ps)+len(sp)
            stats.append({**entry,'topology_sha256':topo_hash(tr),'rf_distance_to_primary':rf,'normalized_rf':rf/den if den else 0,'n_nontrivial_splits':len(sp)})
        tsv(run/'TREE_INVENTORY.tsv',stats);prepared.append(run/'TREE_INVENTORY.tsv')
        dump(run/'GENES.json',genes);prepared.extend([run/'GENES.json',run/'SAMPLING_PLAN.json'])
        kind=cfg['bootstrap_kind']
        if kind=='auto':kind='ultrafast' if '.ufboot' in p.name else 'standard'
        dump(run/'SOURCE_PROVENANCE.json',{'bootstrap_file':str(p),'source_ml_tree':str(ml),'report':str(report),'kind':kind,'selection_scope':'FUBAR site-level support across genome-wide bootstrap topologies','branch_lengths':'input trees contain topology only; HyPhy estimates branch lengths from each unchanged gene alignment','unique_selected_topologies':len(set(x['topology_sha256'] for x in stats if x['role']=='bootstrap')),'source_check':'Declared genome-wide bootstrap source, associated source ML/report and matching taxon sets','hyphy_seed':'CLI seed only effective >=2.5.101; older version default stochastic seed, raw outputs retained'})
        prepared.append(run/'SOURCE_PROVENANCE.json')
        audit=run/'source_audit';audit.mkdir(exist_ok=True)
        for source,target in [(report,audit/'bootstrap_source_report.txt'),(ml,audit/'source_ml_tree.nwk'),(primary_path,audit/'primary_tree.nwk')]:
            shutil.copyfile(source,target);prepared.append(target)
        context(run,{**base_context,'prepared':fingerprint(prepared)})
        if fingerprint(inputs)!=original_hashes:raise ValueError('Source inputs changed during preparation')
        dump(run/'PREPARE_STATUS.json',{'state':'complete','planned_tasks':len(plan),'source_tree_count':total,'selected_bootstrap_trees':len(sample_indices)})
        print('PREPARE COMPLETE: '+str(len(plan))+' tree tasks',flush=True)
    except Exception as e:
        dump(run/'PREPARE_STATUS.json',{'state':'failed','error':str(e)});raise

def fit(run,task_id):
    if load(run/'PREPARE_STATUS.json')['state']!='complete':raise ValueError('Preparation did not complete; inspect prepare.console.log')
    ctx=validate_manifest(run);cfg=load(run/'config.json');plan=load(run/'SAMPLING_PLAN.json')['plan'];genes=load(run/'GENES.json')
    entries=[x for x in plan if x['task']==task_id]
    if len(entries)!=1:raise ValueError('Invalid array task ID')
    entry=entries[0];folder=run/'fits'/entry['tree_id'];folder.mkdir(parents=True,exist_ok=True)
    with lock(folder/'.lock'):
        tasks=Tasks(folder/'tasks');state={'state':'running','tree':entry,'completed':[]};dump(folder/'RUN_STATUS.json',state)
        try:
            for gene,g in genes.items():
                for method in ['fubar']:
                    d=folder/gene/method;d.mkdir(parents=True,exist_ok=True);a=d/'alignment.fasta'
                    if not a.exists():shutil.copyfile(g['alignment'],a)
                    if sha(a)!=sha(g['alignment']):raise ValueError('Fit alignment changed')
                    result=d/'result.json';cmd=hyphy_cmd(ctx['hyphy'],ctx['hyphy_version'],cfg['threads'],cfg['tree_selection_seed']+task_id,method,a,g['trees'][entry['tree_id']],result)
                    tasks.run(gene+'_'+method,cmd,[result],lambda:json_sites(result,method,g['ncodons'],g['ntips']))
                    state['completed'].append(gene+'/'+method);dump(folder/'RUN_STATUS.json',state)
            validate_manifest(run);state['state']='complete'
        except Exception as e:state.update(state='failed',error=str(e));raise
        finally:dump(folder/'RUN_STATUS.json',state)

def quantile(values,p):
    a=sorted(values)
    if not a:return None
    k=(len(a)-1)*p;i=int(k);return a[i]+(a[min(i+1,len(a)-1)]-a[i])*(k-i)

def collect(run):
    cfg=load(run/'config.json');out=run/'summary';out.mkdir(exist_ok=True)
    rows=[];statuses=[];per_tree={};planned=cfg['n_bootstrap_trees']
    if not (run/'PREPARE_STATUS.json').exists() or load(run/'PREPARE_STATUS.json')['state']!='complete':
        raise ValueError('No completed preparation')
    ctx=validate_manifest(run);genes=load(run/'GENES.json');plan=load(run/'SAMPLING_PLAN.json')['plan']
    for entry in plan:
        current=[];errors=[];folder=run/'fits'/entry['tree_id']
        for gene,settings in genes.items():
            try:
                result=folder/gene/'fubar/result.json'
                done=load(folder/'tasks'/(gene+'_fubar')/'DONE.json')
                if done['outputs'].get(str(result))!=sha(result):raise ValueError('Result checksum mismatch')
                probabilities=json_sites(result,'fubar',settings['ncodons'],settings['ntips'])
                current.extend(dict(tree_id=entry['tree_id'],role=entry['role'],source_index=entry['source_index'],
                    gene=gene,codon=i+1,fubar_posterior_positive=probabilities[i]) for i in range(settings['ncodons']))
            except Exception as error:errors.append(gene+': '+str(error))
        if not errors and len(current)==2594:
            per_tree[entry['tree_id']]=current;rows.extend(current)
        statuses.append({**entry,'status':'failed_or_incomplete' if errors else 'complete','error':' | '.join(errors)})
    tsv(out/'TREE_TASK_STATUS.tsv',statuses)
    tsv(out/'ALL_TREE_SITE_RESULTS.tsv',rows,['tree_id','role','source_index','gene','codon','fubar_posterior_positive'])
    boot=[entry for entry in statuses if entry['role']=='bootstrap' and entry['status']=='complete']
    all_complete=len(boot)==planned and 'primary' in per_tree
    bysite=collections.defaultdict(list);control={}
    for row in rows:
        if row['role']=='bootstrap':bysite[row['gene'],row['codon']].append(row)
        else:control[row['gene'],row['codon']]=row
    summary=[]
    for gene,ncodons in cfg['genes'].items():
        for codon in range(1,ncodons+1):
            sampled=bysite[gene,codon];baseline=control.get((gene,codon),{})
            count=sum(row['fubar_posterior_positive']>cfg['fubar_posterior_threshold'] for row in sampled)
            row=dict(gene=gene,codon=codon,planned_bootstrap_trees=planned,complete_bootstrap_trees=len(sampled),
                primary_fubar_posterior=baseline.get('fubar_posterior_positive'),fubar_gt_09_count=count,
                fubar_gt_09_fraction=count/planned if all_complete else None,
                stability_fraction_status='complete_descriptive_only' if all_complete else 'withheld_incomplete_plan')
            probabilities=[row['fubar_posterior_positive'] for row in sampled]
            for name,probability in [('min',0),('q025',.025),('median',.5),('q975',.975),('max',1)]:
                row['fubar_posterior_positive_'+name]=quantile(probabilities,probability)
            summary.append(row)
    tsv(out/'SITE_STABILITY.tsv',summary)
    original=[]
    for gene,ncodons in cfg['genes'].items():
        probabilities=json_sites(ctx['baseline_paths'][gene]['fubar'],'fubar',ncodons,genes[gene]['ntips'])
        for i in range(ncodons):
            baseline=control.get((gene,i+1),{});p=baseline.get('fubar_posterior_positive')
            original.append(dict(gene=gene,codon=i+1,archived_fubar_posterior=probabilities[i],
                control_fubar_posterior=p,fubar_threshold_agrees=(probabilities[i]>.9)==(p>.9) if p is not None else None))
    tsv(out/'PRIMARY_CONTROL_COMPARISON.tsv',original)
    dump(run/'RUN_STATUS.json',dict(state='complete' if all_complete else 'incomplete',
        planned_bootstrap_trees=planned,complete_bootstrap_trees=len(boot),primary_complete='primary' in per_tree,
        interpretation='Descriptive FUBAR support across sampled genome-wide topologies'))
    return 0 if all_complete else 1


def status(run):
    for f in ['PREPARE_STATUS.json','RUN_STATUS.json']:
        if (run/f).exists():print(f, json.dumps(load(run/f),ensure_ascii=False))
    for p in sorted((run/'fits').glob('*/RUN_STATUS.json')):
        d=load(p);print(p.parent.name,d['state'],str(len(d['completed']))+'/4')

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','fit','collect','status']);p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--run',type=Path,required=True);p.add_argument('--task',type=int);p.add_argument('--config',type=Path)
    a=p.parse_args();run=a.run.resolve();root=a.root.resolve();run.mkdir(parents=True,exist_ok=True)
    if a.action=='prepare':
        with lock(run/'.prepare.lock'):
            if a.config:
                source=a.config.resolve();check_cfg(load(source))
                if (run/'config.json').exists():
                    if sha(run/'config.json')!=sha(source):raise ValueError('Run configuration is already frozen and differs from --config')
                else:shutil.copy2(source,run/'config.json')
            if not (run/'config.json').is_file():raise ValueError('Provide --config with explicit bootstrap, source ML and provenance report paths. The template is analysis/selection_topologies/config.json.')
            prepare(root,run)
    elif a.action=='fit':fit(run,a.task)
    elif a.action=='collect':
        with lock(run/'.collect.lock'):
            try:code=collect(run)
            except Exception as e:
                dump(run/'RUN_STATUS.json',{'state':'failed','error':str(e)});traceback.print_exc();code=1
            raise SystemExit(code)
    else:status(run)
if __name__=='__main__':main()
