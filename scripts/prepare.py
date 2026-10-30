#!/usr/bin/env python3
"""Prepare plotting inputs from existing local study outputs, without model fitting."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys

PACKAGE=Path(__file__).resolve().parent
GENES={'mrcA':884,'pbp1':624,'pbp2':652,'tp47':434}
FOCAL=['mrcA:p.Ala506Thr','mrcA:p.Ala506Val','mrcA:p.Gly708Ser','mrcA:p.Met625Val','pbp1:p.Pro564Leu','tp47:p.Ser394Arg']
MISSING={'','NA','NaN','nan','None'}
AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODONS=dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),AA))


def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        records=list(csv.DictReader(f,delimiter='\t',escapechar='\\'))
    if any(None in r or any(v is None for v in r.values()) for r in records):raise ValueError('Malformed table: '+str(path))
    return records


def write(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);rows=list(rows)
    fields=fields or list(rows[0])
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def fasta(path):
    result={};name=None;parts=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:result[name]=''.join(parts).upper()
            name=line[1:].split()[0];parts=[]
            if name in result:raise ValueError('Duplicate alignment ID')
        else:parts.append(line.strip())
    if name is not None:result[name]=''.join(parts).upper()
    return result




def hyphy(path,n):
    obj=json.loads(path.read_text());headers=[h[0] if isinstance(h,list) else h for h in obj['MLE']['headers']];result={}
    for part,values in obj['MLE']['content'].items():
        coverage=obj['data partitions'][part]['coverage']
        while len(coverage)==1 and isinstance(coverage[0],list):coverage=coverage[0]
        if len(coverage)!=len(values):raise ValueError('HyPhy coverage and row counts differ')
        for site,row in zip(coverage,values):
            if len(row)<len(headers):raise ValueError('Short HyPhy row')
            site=int(site)+1
            if site in result:raise ValueError('Duplicate HyPhy site')
            result[site]=dict(zip(headers,row))
    if set(result)!=set(range(1,n+1)):raise ValueError('Incomplete codon results: '+str(path))
    return result


class Sources:
    def __init__(self,project,extra,overrides):
        self.project=project.resolve();self.roots=list(dict.fromkeys([self.project]+[p.resolve() for p in extra]))
        self.paths=json.loads((PACKAGE/'config/source_paths.json').read_text());self.used={};self.records=[]
        if overrides:self.paths.update(json.loads(overrides.read_text()))

    def find(self,key,required=True):
        candidates=self.paths.get(key,[])
        if isinstance(candidates,str):candidates=[candidates]
        found=None
        for candidate in candidates:
            p=Path(candidate)
            options=[p] if p.is_absolute() else [root/p for root in self.roots]
            for option in options:
                if option.is_file() and option.stat().st_size:
                    found=option.resolve();break
            if found:break
        self.records.append({'input':key,'required':required,'status':'FOUND' if found else 'MISSING',
                             'path':str(found) if found else '', 'expected_paths':';'.join(candidates)})
        if found:self.used[key]=found
        return found

    def manifest(self,out):
        write(out/'INPUT_RESOLUTION.tsv',self.records)
        (out/'SOURCE_MANIFEST.json').write_text(json.dumps({k:{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for k,p in self.used.items()},indent=2)+'\n')


def selection(sources,out,matrix):
    indicators=defaultdict(list)
    for name in matrix[0]:
        match=re.fullmatch(r'([^:]+):p\.[A-Za-z]+([0-9]+)[A-Za-z]+',name)
        if match:indicators[match[1],int(match[2])].append(name)
    rows=[];counts={}
    for gene,n in GENES.items():
        seqs=fasta(sources.used[gene+'_codon.fasta']);counts[gene]=len(seqs)
        if {len(x) for x in seqs.values()}!={n*3}:raise ValueError('Codon alignment length differs for '+gene)
        if set(sources.used[gene+'_ids.txt'].read_text().split())!=set(seqs):raise ValueError('Selection ID list and alignment membership differ for '+gene)
        fubar=hyphy(sources.used[gene+'_fubar.json'],n)
        for i in range(1,n+1):
            observed=Counter(seq[(i-1)*3:i*3] for seq in seqs.values() if set(seq[(i-1)*3:i*3])<=set('ACGT'))
            amino=Counter()
            for codon,count in observed.items():amino[CODONS[codon]]+=count
            f=fubar[i];positive=float(f['Prob[alpha<beta]'])
            rows.append(dict(gene=gene,site=i,FUBAR_Ppositive=positive,
                FUBAR_Pnegative=f['Prob[alpha>beta]'],FUBAR_BayesFactor_positive=f['BayesFactor[alpha<beta]'],
                FUBAR_Ppositive_above_09=positive>.9,
                callable_codon_sequences=sum(observed.values()),observed_codon_types=len(observed),
                amino_acid_counts_in_selection_alignment=json.dumps(dict(sorted(amino.items()))),
                matrix_indicators=';'.join(sorted(indicators[gene,i]))))
    write(out/'selection.tsv',rows)
    return rows,counts


def esm(sources,out):
    values={};by_site=defaultdict(list)
    for r in read(sources.used['ESM']):
        if r['protein'] not in GENES or r.get('scan_mode','full_dms')!='full_dms':continue
        gene=r['protein'];site=int(r['position']);key=gene,site,r['wt_aa'],r['mut_aa'];value=float(r['llr'])
        if key in values or not math.isfinite(value) or not 1<=site<=GENES[gene] or r['wt_aa']==r['mut_aa'] or not {r['wt_aa'],r['mut_aa']}<=set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('Invalid ESM key/score')
        values[key]=value;by_site[gene,site].append(value)
    if len(values)!=sum(GENES.values())*19 or len(by_site)!=sum(GENES.values()) or {len(x) for x in by_site.values()}!={19}:raise ValueError('Incomplete four-gene ESM scan')
    rows=[dict(gene=g,site=i,score=sum(by_site[g,i])/19,n_substitutions_for_mean=19) for g,n in GENES.items() for i in range(1,n+1)]
    write(out/'esm_site_means.tsv',rows)
    highlights=[('mrcA',506,'A','M'),('mrcA',506,'A','T'),('mrcA',506,'A','V'),('mrcA',625,'M','V'),('mrcA',708,'G','S'),('pbp1',564,'P','L'),('tp47',394,'S','R')]
    write(out/'esm_specific.tsv',[dict(gene=g,site=i,wt=w,mut=m,label=w+str(i)+m,llr=values[g,i,w,m],n_rows=1,scan_mode='full_dms') for g,i,w,m in highlights])
    sys.path.insert(0,str(PACKAGE/'helpers'));from esm_observed_table import build
    observed=build(sources.used['matrix.tsv'],sources.used['ESM'])
    write(out/'esm_observed.tsv',observed)
    return len(values)


def sasa(pdb,out):
    """Shrake–Rupley reference SASA; Fibonacci sphere and standard atomic radii."""
    import numpy as np
    from scipy.spatial import cKDTree
    atoms=[x for x in pdb.read_text().splitlines() if x.startswith('ATOM')]
    xyz=np.array([[float(x[30:38]),float(x[38:46]),float(x[46:54])] for x in atoms])
    radius=np.array([{'C':1.7,'N':1.55,'O':1.52,'S':1.8,'H':1.2,'P':1.8}[x[76:78].strip()]+1.4 for x in atoms])
    tree=cKDTree(xyz);rows=[]
    for n in [2880]:
        sphere=np.zeros((n,3),dtype=np.float32);dz=2./n;angle=0.;z=1-dz/2;increment=np.pi*(3-np.sqrt(5))
        for i in range(n):
            rr=np.sqrt(1-z*z);sphere[i]=[np.cos(angle)*rr,np.sin(angle)*rr,z];z-=dz;angle+=increment
        for pos in [506,625,708]:
            indices=[i for i,x in enumerate(atoms) if x[21]=='A' and int(x[22:26])==pos];area=0.
            if not indices:raise ValueError('Structural residue missing')
            for i in indices:
                pts=sphere*radius[i]+xyz[i];visible=np.ones(n,dtype=bool)
                for j in tree.query_ball_point(xyz[i],radius[i]+max(radius)):
                    if j!=i and np.sum((xyz[i]-xyz[j])**2)<(radius[i]+radius[j])**2:
                        visible &= np.sum((pts-xyz[j])**2,axis=1)>radius[j]**2
                area+=int(visible.sum())*radius[i]**2*4*np.pi/n
            rows.append(dict(chain='A',position=pos,residue=atoms[indices[0]][17:20],n_points=n,probe_radius_A=1.4,SASA_A2=area,rounded_A2=f'{area:.2f}',n_atoms=len(indices)))
    write(out/'sasa_resolution.tsv',rows)


def historical(sources,out):
    direct=sources.find('temporal_extended_historical.tsv',False)
    if direct:shutil.copyfile(direct,out/'temporal_extended_historical.tsv');return
    allrows=[]
    plans=[('01_reference_substitution/05_temporal',{'full_nicholsref','nichols_lineage_nicholsref'},'Stage 1: Reference substitution'),
           ('02_nichols_subsampling',{'profile_A_light_downsample','profile_B_aggressive_downsample','profile_C_era_stratified'},'Stage 2: Nichols-lineage subsampling'),
           ('04_crossref_validation',{'crossref_a_ss14_alone_nicholsref','crossref_b_nichols_alone_light_nicholsref'},'Stage 4: Cross-reference validation (mirror)')]
    for relative,keys,section in plans:
        options=[r/'results/p6_01_nichols_temporal_signal'/relative/'temporal_signal_summary.tsv' for r in sources.roots]
        source=next((p for p in options if p.is_file()),None)
        if not source:return
        selected=[r for r in read(source) if r.get('dataset',r.get('profile')) in keys]
        if len(selected)!=len(keys):return
        for r in selected:
            slope=float(r['slope']);positive=slope>0 and float(r['real_r2'])>float(r['perm_95pct'])
            allrows.append(dict(section=section,display_dataset=r['label'],n=r['n_samples'],R2=r['real_r2'],slope=r['slope'],p_value=r['p_value'],temporal_signal_class='Temporal signal supported' if positive else 'Negative clock slope' if slope<0 else 'No supported temporal signal'))
        sources.used['historical_'+relative]=source
    if len(allrows)==7:
        write(out/'temporal_extended_historical.tsv',allrows)
        for record in sources.records:
            if record['input']=='temporal_extended_historical.tsv':record.update(status='DERIVED',path=str(out/'temporal_extended_historical.tsv'))


def prepare(project,out,extra=(),overrides=None):
    out.mkdir(parents=True,exist_ok=False);data=out/'work/figures/data';data.mkdir(parents=True)
    sources=Sources(project,extra,overrides)
    required=['matrix.tsv','metadata.tsv','ASR','ESM','mrcA_model.pdb','clock_rates.tsv','dates.tsv','er_ard.tsv','haplotype_lineages.tsv','haplotype_samples.tsv','metadata_provenance.tsv','nucleotide_haplotypes.tsv','parsimony.tsv','primary.treefile','relaxed.tree','relaxed_calibration.tsv','relaxed_origins.tsv','rtt_points.tsv','rtt_summary.tsv','simmap_draws.tsv','source_project_crosswalk.tsv','strict.tree','strict_calibration.tsv','strict_origins.tsv','structural_metrics.tsv','temporal_signal.tsv']
    required += [g+'_'+ext for g in GENES for ext in ['fubar.json','codon.fasta','ids.txt']]
    missing=[k for k in required if sources.find(k) is None]
    for key in ['TOPOLOGY','TableS12A_clustered_date_tests.tsv','TableS12B_subset_stability.tsv','TableS13B_strict_ESS.tsv','TableS13C_relaxed_ESS.tsv']:sources.find(key,False)
    sources.manifest(out)
    if missing:raise ValueError('Required cluster inputs missing; see INPUT_RESOLUTION.tsv: '+', '.join(missing))
    for name,path in sources.used.items():
        if name.endswith(('.tsv','.tree','.treefile','.pdb','.txt')) and '/' not in name:shutil.copyfile(path,data/name)
    matrix=read(data/'matrix.tsv');meta=read(data/'metadata.tsv')
    if len(matrix)!=1861 or len(meta)!=1861 or len({r['sample'] for r in matrix})!=1861 or {r['sample'] for r in matrix}!={r['sample'] for r in meta}:raise ValueError('Study cohort mismatch')
    geo=read(PACKAGE/'config/geography.tsv');lookup={r['Geo_Country'].replace('_',' '):r for r in geo};counts=Counter()
    for row in meta:
        key=row['country'].replace('_',' ')
        if key not in lookup:raise ValueError('Country annotation missing: '+key)
        row['continent']=lookup[key]['continent'];counts[key]+=1
    write(data/'metadata.tsv',meta)
    fields=[k for k in geo[0] if k!='continent']
    write(data/'map_coordinates.tsv',[{**{k:r[k] for k in fields},'Count':counts[r['Geo_Country'].replace('_',' ')]} for r in geo if counts[r['Geo_Country'].replace('_',' ')]])
    selected,aln_counts=selection(sources,data,matrix)
    nr=esm(sources,data)
    # Merge existing estimates; do not refit transition or selection models.
    asr={r['locus']:r for r in read(sources.used['ASR'])};er={r['locus']:r for r in read(data/'er_ard.tsv')}
    if len(er)!=23 or any(r['analysis_status']!='complete' for r in er.values()):raise ValueError('Use the completed 23-row ER/ARD table, including R107G')
    par={r['locus']:r for r in read(data/'parsimony.tsv') if r['tree_scope']=='primary_full' and r['root_policy']=='unrestricted'}
    sites={(r['gene'],r['site']):r for r in selected};focal=[]
    for locus in FOCAL:
        gene=locus.split(':')[0];site=int(re.search(r'([0-9]+)',locus.split(':')[1])[1]);row={'locus':locus,**sites[gene,site]}
        row.update({k:asr[locus][k] for k in ['n_mutant','n_wildtype','n_missing','simmap_gain_mean','simmap_gain_q025','simmap_gain_q975']})
        row.update({k:par[locus][k] for k in ['minimum_gains_among_optima','maximum_gains_among_optima','optimal_root_states']})
        row.update({k:er[locus][k] for k in ['delta_AIC_ARD_ER','AIC_evidence','ARD_root_Pmut']});focal.append(row)
    write(data/'focal.tsv',focal)
    hs=read(data/'haplotype_samples.tsv')
    sys.path.insert(0,str(PACKAGE/'helpers/lineage_parsimony'));from treeio import read_tree
    dated=set(read_tree(data/'strict.tree').tips);relaxed=set(read_tree(data/'relaxed.tree').tips)
    if dated!=relaxed or len(dated)!=533 or not dated<={r['sample'] for r in matrix}:raise ValueError('Two-clock sample membership differs from the selected 533-genome cohort')
    ndate=len(dated)
    coverage=[dict(analysis='Whole cohort',n=len(matrix))]+[dict(analysis=g+' codon alignment',n=aln_counts[g]) for g in GENES]
    coverage += [dict(analysis='Complete mrcA sequence',n=sum(r['gene_complete_acgt'].lower()=='true' for r in hs)),dict(analysis='SS14 dating subset',n=ndate)]
    write(data/'coverage.tsv',coverage)
    sasa(data/'mrcA_model.pdb',data)
    if 'TOPOLOGY' in sources.used:shutil.copyfile(sources.used['TOPOLOGY'],data/'topology_stability.tsv')
    historical(sources,data);sources.manifest(out)
    lineage_out=out/'lineage_history'
    import subprocess
    with (out/'lineage_history.log').open('w') as f:subprocess.run([sys.executable,str(PACKAGE/'helpers/lineage_parsimony/analyse.py'),'--input-dir',str(data),'--outdir',str(lineage_out)],stdout=f,stderr=subprocess.STDOUT,check=True)
    shutil.copyfile(lineage_out/'lineage_summary.tsv',data/'lineage_summary.tsv')
    shutil.copytree(PACKAGE/'figures/R',out/'work/figures/R');shutil.copytree(PACKAGE/'figures/qa',out/'work/figures/qa')
    (out/'work/INPUT_MANIFEST.json').write_text(json.dumps({str(p.relative_to(data)):sha(p) for p in data.iterdir() if p.is_file()},indent=2)+'\n')
    optional=[r['input'] for r in sources.records if not r['required'] and r['status']=='MISSING']
    (out/'PREPARE_STATUS.json').write_text(json.dumps(dict(state='complete',cohort=len(matrix),codons=len(selected),ESM_scores=nr,optional_sources_missing=optional,scientific_models_rerun=False),indent=2)+'\n')
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--extra-root',action='append',type=Path,default=[]);p.add_argument('--inputs',type=Path);a=p.parse_args()
    print(prepare(a.project,a.out,a.extra_root,a.inputs))
