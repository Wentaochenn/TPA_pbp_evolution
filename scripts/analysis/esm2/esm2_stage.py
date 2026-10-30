#!/usr/bin/env python3
"""Compute and validate the four-protein reference-context ESM-2 scan."""
import argparse, csv, hashlib, importlib.util, json, math, os, platform
import subprocess, sys, traceback
from importlib.metadata import version
from pathlib import Path

PACKAGE=Path(__file__).resolve().parent
AA=set('ACDEFGHIKLMNPQRSTVWY')
GENES={'mrcA':884,'pbp1':624,'pbp2':652,'tp47':434}
FOCAL=[('mrcA',506,'T'),('mrcA',506,'V'),('mrcA',506,'M'),('mrcA',625,'V'),('mrcA',708,'S'),('pbp1',564,'L'),('tp47',394,'R')]
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def write(p,rows,fields):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta(path):
    found={};name=None
    for line in Path(path).read_text().splitlines():
        if not line.strip():continue
        if line.startswith('>'):
            name=line[1:].split()[0]
            if name in found:raise ValueError('Duplicate FASTA ID: '+name)
            found[name]=''
        elif name:found[name]+=line.strip().upper()
        else:raise ValueError('Missing FASTA header')
    for gene,length in GENES.items():
        if gene not in found or len(found[gene])!=length:raise ValueError('Reference protein length differs: '+gene)
        if set(found[gene])-AA:raise ValueError('Noncanonical amino acid or stop codon in '+gene)
        if not found[gene].startswith('M'):raise ValueError('Reference CDS must start with initiating methionine: '+gene)
    for gene,pos,wt in [('mrcA',506,'A'),('mrcA',625,'M'),('mrcA',708,'G'),('pbp1',564,'P'),('tp47',394,'S')]:
        if found[gene][pos-1]!=wt:raise ValueError('Nichols reference residue differs: '+gene+str(pos))
    return {g:found[g] for g in GENES}
def audit_scores(path,sequences):
    accepted={};outside=0
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        for r in reader:
            gene=r.get('protein') or r.get('gene')
            if gene not in sequences:outside+=1;continue
            mode=r.get('scan_mode','full_dms')
            if mode!='full_dms':outside+=1;continue
            pos=int(r.get('position') or r.get('site'));wt=r.get('wt_aa') or r.get('wt');alt=r.get('mut_aa') or r.get('mut')
            value=float(r.get('llr') or r.get('ESM2_LLR'))
            if not 1<=pos<=len(sequences[gene]) or sequences[gene][pos-1]!=wt:raise ValueError('Reference residue mismatch in score file')
            if alt not in AA or alt==wt or not math.isfinite(value):raise ValueError('Invalid substitution or non-finite score')
            key=(gene,pos,alt)
            if key in accepted:raise ValueError('Duplicate substitution row: '+str(key))
            accepted[key]={'protein':gene,'position':pos,'wt_aa':wt,'mut_aa':alt,'llr':value}
    missing=[];means=[]
    for gene,seq in sequences.items():
        for pos,wt in enumerate(seq,1):
            keys=[(gene,pos,alt) for alt in sorted(AA-{wt})]
            absent=[k for k in keys if k not in accepted]
            if absent:missing.append({'gene':gene,'site':pos,'available':19-len(absent)});continue
            means.append({'gene':gene,'site':pos,'score':sum(accepted[k]['llr'] for k in keys)/19,'n_substitutions_for_mean':19})
    summary={'expected_substitutions':sum(map(len,sequences.values()))*19,'observed_core_full_scan_rows':len(accepted),'expected_positions':sum(map(len,sequences.values())),'positions_with_all_19_alternatives':len(means),'missing_positions':missing,'rows_outside_core_full_scan':outside,'coverage_complete':not missing}
    return summary,accepted,means
def export(path,sequences,out):
    summary,accepted,means=audit_scores(path,sequences)
    dump(out/'ESM2_COVERAGE.json',summary)
    if not summary['coverage_complete']:raise ValueError('Incomplete full substitution matrix; see ESM2_COVERAGE.json')
    write(out/'esm2_all_substitutions.tsv',accepted.values(),['protein','position','wt_aa','mut_aa','llr'])
    write(out/'esm2_site_means.tsv',means,['gene','site','score','n_substitutions_for_mean'])
    focal=[]
    for gene,pos,alt in FOCAL:
        if gene in sequences:
            r=accepted[(gene,pos,alt)];focal.append({'gene':gene,'site':pos,'wt':r['wt_aa'],'mut':alt,'label':r['wt_aa']+str(pos)+alt,'llr':r['llr'],'n_rows':1,'scan_mode':'full_dms'})
    write(out/'esm2_focal_scores.tsv',focal,['gene','site','wt','mut','label','llr','n_rows','scan_mode'])
    return summary
def compare_frozen(out,freeze):
    specifications=[('esm2_site_means.tsv','esm_site_means.tsv',['gene','site'],'score'),('esm2_focal_scores.tsv','esm_specific.tsv',['gene','site','mut'],'llr')]
    report=[]
    for actual,expected,keys,column in specifications:
        def keyed(path):
            with path.open(encoding='utf-8-sig',newline='') as handle:
                data=list(csv.DictReader(handle,delimiter='\t'))
            result={tuple(r[k] for k in keys):float(r[column]) for r in data}
            if len(result)!=len(data):raise ValueError('Duplicate keys in ESM-2 comparison table: '+str(path))
            return result
        a=keyed(out/actual);b=keyed(freeze/'figures/data'/expected)
        common=a.keys()&b.keys();deltas=[abs(a[k]-b[k]) for k in common]
        report.append({'generated_table':actual,'frozen_table':expected,'matched_rows':len(common),'generated_rows':len(a),'frozen_rows':len(b),'keys_equal':a.keys()==b.keys(),'max_absolute_difference':max(deltas,default=None),'all_scores_match_within_1e-6':a.keys()==b.keys() and all(d<=1e-6 for d in deltas)})
    result={'frozen_release':str(freeze),'absolute_tolerance':1e-6,'tables':report,'matches_frozen':all(r['all_scores_match_within_1e-6'] for r in report),'frozen_files_modified':False}
    dump(out/'FROZEN_SCORE_COMPARISON.json',result)
    return result
def run(args):
    sequences=fasta(args.protein_fasta);out=args.out.resolve()
    if out.exists():raise ValueError('Use a new ESM-2 output directory')
    out.mkdir(parents=True)
    dump(out/'RUN_STATUS.json',{'state':'running','stage':'ESM2','scores_do_not_update_manuscript_freeze':True})
    try:
        if args.action=='audit':
            source=args.scores.resolve();summary=export(source,sequences,out)
            context={'score_source':str(source),'score_sha256':sha(source)}
        else:
            model=args.model.resolve()
            if not model.is_dir():raise ValueError('Local ESM-2 model folder is required')
            config=json.loads((model/'config.json').read_text())
            if config.get('model_type')!='esm' or config.get('num_hidden_layers')!=36 or config.get('hidden_size')!=2560:raise ValueError('Expected the existing ESM2 t36 3B model architecture')
            os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1'
            import torch
            if args.device=='cuda' and not torch.cuda.is_available():raise ValueError('No allocated CUDA GPU; this entry does not request a GPU from SGE')
            if args.device=='cuda' and not torch.cuda.is_bf16_supported():raise ValueError('The existing bfloat16 calculation requires a compatible GPU')
            spec=importlib.util.spec_from_file_location('esm2_reference_scorer',PACKAGE/'scorer.py')
            scorer=importlib.util.module_from_spec(spec);spec.loader.exec_module(scorer)
            files=[p for p in sorted(model.rglob('*')) if p.is_file()]
            context={'model_path':str(model),'model_files':[{'path':str(p.relative_to(model)),'sha256':sha(p),'bytes':p.stat().st_size} for p in files],'torch':torch.__version__,'python':platform.python_version(),'transformers':version('transformers'),'pandas':version('pandas'),'cuda_runtime':torch.version.cuda,'cudnn':torch.backends.cudnn.version(),'device':args.device,'dtype':'bfloat16','scoring_source_sha256':sha(PACKAGE/'scorer.py')}
            if args.device=='cuda':context['gpu']=torch.cuda.get_device_name(0)
            dump(out/'MODEL_CONTEXT.json',context)
            logger=scorer.setup_logging(out/'esm2.console.log')
            network,tokenizer=scorer.load_model(str(model),args.device,logger)
            for aa in AA:
                if tokenizer.convert_tokens_to_ids(aa)==tokenizer.unk_token_id:raise ValueError('Tokenizer lacks a canonical amino acid')
            table=scorer.run_full_dms_mode(network,tokenizer,sequences,args.device,logger,proteins_to_scan=list(GENES))
            source=out/'esm2_dms_results.tsv';table.to_csv(source,sep='\t',index=False)
            summary=export(source,sequences,out)
        dump(out/'INPUT_CONTEXT.json',dict(context,protein_fasta=str(args.protein_fasta.resolve()),protein_fasta_sha256=sha(args.protein_fasta)))
        comparison=compare_frozen(out,args.frozen_release.resolve()) if args.frozen_release else None
        dump(out/'RUN_STATUS.json',{'state':'complete','stage':'ESM2','exit_code':0,'summary':summary,'frozen_score_comparison':comparison,'manuscript_freeze_modified':False})
        print('ESM2_OUTPUT:',out)
    except Exception as e:
        dump(out/'RUN_STATUS.json',{'state':'failed','exit_code':1,'error':str(e),'traceback':traceback.format_exc()});raise
def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('recompute',help='Translate the annotated reference and compute scores; optionally compare previous results')
    q.add_argument('--project',type=Path,required=True);q.add_argument('--model',type=Path,required=True)
    q.add_argument('--out',type=Path,required=True);q.add_argument('--previous-scores',type=Path,help='Optional complete score table for a keyed comparison')
    q.add_argument('--reference',type=Path);q.add_argument('--genbank',type=Path)
    q.add_argument('--device',choices=['cuda','cpu'],default='cuda')
    for action in ['audit','run']:
        a=s.add_parser(action);a.add_argument('--protein-fasta',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--frozen-release',type=Path)
        if action=='audit':a.add_argument('--scores',type=Path,required=True)
        else:a.add_argument('--model',type=Path,required=True);a.add_argument('--device',choices=['cuda','cpu'],default='cuda')
    args=p.parse_args()
    if args.action=='recompute':
        from recompute import recompute
        recompute(args)
    else:run(args)
if __name__=='__main__':main()
