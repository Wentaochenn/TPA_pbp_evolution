"""Sampling sensitivity: complete dated lineage + 20 signal-blind random subsets."""
import argparse, collections, copy, csv, json, math, random, sys
from pathlib import Path
import numpy as np
import Bio
from Bio import Phylo
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.common import sha,dump,tsv,read_tsv,resolve,fingerprint,context,validate_tree,prune,write_tree,check_known_hashes,midpoint_tree
from analysis.clustered_dates.cluster_core import xml_dates,cluster_indices,permutation_test

def metadata_rows(path):
    if path.suffix.lower()=='.xlsx':
        import openpyxl
        wb=openpyxl.load_workbook(path,read_only=True,data_only=True)
        it=wb.worksheets[0].iter_rows(values_only=True);headers=[str(v).strip() if v is not None else '' for v in next(it)]
        if len(set(headers))!=len(headers):raise ValueError('Duplicate Excel header')
        rows=[dict(zip(headers,values)) for values in it];wb.close();return rows
    return read_tsv(path)

def load_dates(path,cfg):
    if cfg['dates_header']:
        rows=read_tsv(path);pairs=[(r[cfg['dates_sample_column']],r[cfg['dates_year_column']]) for r in rows]
    else:
        with path.open() as f:
            rows=list(csv.reader(f,delimiter='\t'))
        if any(len(r)!=2 for r in rows if r):raise ValueError('Expected two-column headerless dates file')
        pairs=[r for r in rows if r]
    dates={};excluded=[];seen=set()
    for key,value in pairs:
        key=str(key).strip()
        if not key or key in seen:raise ValueError('Missing/duplicate sample date identifier '+key)
        seen.add(key)
        try:v=float(value)
        except (TypeError,ValueError):excluded.append({'sample':key,'reason':'invalid_or_missing_date'});continue
        if not math.isfinite(v) or not 1800<=v<=2100:excluded.append({'sample':key,'reason':'invalid_or_out_of_range_date'});continue
        dates[key]=v
    return dates,excluded

def memberships(names,n,count,seed):
    names=sorted(names)
    if n>=len(names):raise ValueError('Random subset must be smaller than complete eligible cohort')
    return {f'random_{i+1:02d}':sorted(random.Random(seed+i).sample(names,n)) for i in range(count)}

def main(a):
    root=a.root.resolve();out=a.outdir.resolve();out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads(a.config.read_text());paths={k:resolve(root,cfg[k]) for k in ['primary_tree','metadata','dates','selected_xml']}
    check_known_hashes(root,cfg,paths.values())
    source=fingerprint(paths.values());code=fingerprint(Path(__file__).parent.glob('*.py'))
    context(out,{'settings':cfg,'inputs':source,'code':code,'numpy':np.__version__,'biopython':Bio.__version__,'design':'full_cohort_and_unconditional_random_subsets; no PBP rescue or rtd/date ranking'})
    status={'state':'running','complete_memberships':[]};dump(out/'RUN_STATUS.json',status)
    try:
        if cfg['random_subsets']<1 or cfg['permutations']<99:raise ValueError('Require positive subset count and at least 99 permutations')
        if cfg['tree_units']!='substitutions-per-site':raise ValueError('A genetic-distance tree is required')
        tree=Phylo.read(str(paths['primary_tree']),'newick');tips=validate_tree(tree)
        dates,exclusions=load_dates(paths['dates'],cfg);meta={};sc=cfg['metadata_sample_column'];lc=cfg['metadata_lineage_column']
        for row in metadata_rows(paths['metadata']):
            if sc not in row or lc not in row:raise ValueError('Configured metadata columns absent')
            if row[sc] is None or str(row[sc]).strip()=='':continue
            key=str(row[sc]).strip();lineage=str(row[lc] or '').strip()
            if key in meta:raise ValueError('Duplicate metadata sample '+key)
            meta[key]=lineage
        accepted={s.lower() for s in cfg['lineage_values']};eligible=[]
        for key in sorted(tips):
            if key not in meta:reason='no_metadata'
            elif meta[key].lower() not in accepted:reason='outside_explicit_lineage_labels'
            elif key not in dates:reason='no_valid_date'
            else:eligible.append(key);continue
            exclusions.append({'sample':key,'reason':reason})
        tsv(out/'EXCLUSIONS.tsv',exclusions,['sample','reason'])
        tsv(out/'ELIGIBLE_SAMPLES.tsv',[{'sample':k,'lineage':meta[k],'date':dates[k]} for k in eligible])
        selected=xml_dates(paths['selected_xml'])
        difference=[{'sample':k,'xml_date':v,'source_date':dates.get(k,''),'eligible':k in eligible} for k,v in selected.items() if k not in eligible or k not in dates or abs(dates[k]-v)>1e-8]
        tsv(out/'SELECTED_DATE_DIFFERENCES.tsv',difference,['sample','xml_date','source_date','eligible'])
        if difference:raise ValueError('Original selected XML dates/membership differ from eligible cohort; inspect SELECTED_DATE_DIFFERENCES.tsv before comparing')
        if len(selected)!=cfg['subset_n']:raise ValueError('Selected XML count != configured matched subset size')
        groups={'all_eligible':eligible,'original_selected':sorted(selected)}
        groups.update(memberships(eligible,cfg['subset_n'],cfg['random_subsets'],cfg['seed']))
        dump(out/'PREDEFINED_MEMBERSHIPS.json',groups)
        full=prune(tree,eligible);rooted={}
        for rooting in ('inherited','midpoint'):
            tr=copy.deepcopy(full)
            if rooting=='midpoint':tr=midpoint_tree(tr)
            validate_tree(tr)
            rooted[rooting]=tr
        allrows=[]
        for index,(name,keep) in enumerate(groups.items()):
            folder=out/name;folder.mkdir(exist_ok=True);done=folder/'DONE.json';table=folder/'TEMPORAL_RESULTS.tsv'
            if done.exists():
                d=json.loads(done.read_text())
                for f,h in d['outputs'].items():
                    if not Path(f).exists() or sha(f)!=h:raise ValueError('Temporal checkpoint modified '+name)
                rows=read_tsv(table);print('REUSE',name,flush=True)
            else:
                print('RUN',name,len(keep),flush=True);rows=[];saved=[]
                for rooting,full_tree in rooted.items():
                    tr=prune(full_tree,keep)
                    tr.collapse_all(lambda c:c is not tr.root and not c.is_terminal() and c.branch_length==0)
                    names,units=cluster_indices(tr,dates);depths=tr.depths()
                    x=np.array([depths[t] for t in tr.get_terminals()]);year=np.array([dates[k] for k in names])
                    points=[]
                    for j,g in enumerate(units,1):
                        points.extend({'sample':names[i],'date':float(year[i]),'root_to_tip':float(x[i]),'cluster':j,'cluster_n':len(g)} for i in g)
                    pointfile=folder/(rooting+'_points.tsv');tsv(pointfile,points);saved.append(pointfile)
                    treefile=folder/(rooting+'.nwk');write_tree(treefile,tr);saved.append(treefile)
                    for mode,permunits in [('individual',[[i] for i in range(len(names))]),('clustered',units)]:
                        row={'sample_set':name,'rooting':rooting,'randomization':mode,'n':len(keep),'date_min':float(year.min()),'date_max':float(year.max()),'n_date_values':len(set(year))}
                        try:
                            stat,null=permutation_test(x,year,permunits,cfg['permutations'],cfg['seed']+index*10+(mode=='clustered'))
                            slope=float(np.dot(year-year.mean(),x-x.mean())/np.dot(year-year.mean(),year-year.mean()))
                            row.update(stat,slope=slope,status='complete')
                            nullfile=folder/(rooting+'_'+mode+'_null.tsv');tsv(nullfile,({'permutation':i+1,'r':float(v)} for i,v in enumerate(null)));saved.append(nullfile)
                        except ValueError as e:
                            row.update(status='not_estimable',reason=str(e))
                        rows.append(row)
                fields=list(dict.fromkeys(k for r in rows for k in r));tsv(table,rows,fields);saved.append(table)
                dump(done,{'outputs':fingerprint(saved)})
            allrows.extend(rows);status['complete_memberships'].append(name);dump(out/'RUN_STATUS.json',status)
        fields=list(dict.fromkeys(k for r in allrows for k in r));tsv(out/'TEMPORAL_ALL.tsv',allrows,fields)
        summaries=[]
        for rooting in rooted:
            for mode in ('individual','clustered'):
                subset=[r for r in allrows if str(r['sample_set']).startswith('random_') and r['rooting']==rooting and r['randomization']==mode]
                ok=[r for r in subset if r['status']=='complete']
                summaries.append({'rooting':rooting,'randomization':mode,'planned_random_subsets':cfg['random_subsets'],'estimable_subsets':len(ok),'positive_and_p_lt_005':sum(float(r['r'])>0 and float(r['p_one_sided'])<.05 for r in ok),'median_r':float(np.median([float(r['r']) for r in ok])) if ok else None,'interpretation':'descriptive stability across overlapping subsets; not independent replications or a pooled P value'})
        tsv(out/'RANDOM_SUBSET_STABILITY.tsv',summaries)
        if fingerprint(paths.values())!=source:raise ValueError('Source changed during run')
        status.update(state='complete',scope='Sampling sensitivity of root-to-tip signal; no BEAST refit, no calibrated selection-adjusted P for original subset',roots='Both root positions fixed on full eligible genetic tree before subsetting',scientific_review='pending')
    except Exception as e:
        status.update(state='failed',error=str(e));raise
    finally:dump(out/'RUN_STATUS.json',status)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);main(p.parse_args())
