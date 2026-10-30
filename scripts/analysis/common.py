"""Validation, immutable provenance and resumable task utilities."""
import copy, csv, hashlib, json, math, os, subprocess, sys, time
from pathlib import Path
sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def dump(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n'); tmp.replace(path)

def tsv(path,rows,fields=None):
    rows=list(rows); fields=fields or (list(rows[0]) if rows else [])
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def read_tsv(path):
    with open(path,newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))

def resolve(root,value):
    p=Path(value);return p.resolve() if p.is_absolute() else (root/p).resolve()

def fingerprint(paths):
    return {str(Path(p).resolve()):sha(p) for p in paths}

def context(out,obj):
    f=out/'CONTEXT.json'
    if f.exists() and json.loads(f.read_text())!=obj:raise ValueError('Inputs, settings, code or tool versions changed. Use a NEW run directory; do not reuse these checkpoints.')
    dump(f,obj)

def validate_tree(tree):
    tips=tree.get_terminals();names=[t.name for t in tips]
    if None in names or len(set(names))!=len(names):raise ValueError('Missing/duplicate tree identifiers')
    for c in tree.find_clades():
        if c is tree.root and c.branch_length is None:continue
        if c.branch_length is None or not math.isfinite(c.branch_length) or c.branch_length<0:raise ValueError('Missing, negative or non-finite branch length')
    return set(names)

def prune(tree,keep):
    tr=copy.deepcopy(tree);keep=set(keep)
    if not keep<=validate_tree(tr):raise ValueError('Samples absent from tree')
    def rec(c):
        if c.is_terminal():return c if c.name in keep else None
        c.clades=[x for x in (rec(ch) for ch in c.clades) if x is not None]
        if not c.clades:return None
        if len(c.clades)==1 and c is not tr.root:
            ch=c.clades[0];ch.branch_length=(ch.branch_length or 0)+(c.branch_length or 0);return ch
        return c
    tr.root=rec(tr.root)
    if tr.root is None:raise ValueError('Empty retained tree')
    return tr

def write_tree(path,tr):
    from Bio import Phylo
    Phylo.write(tr,str(path),'newick',format_branch_length='%1.12g')

def alignment(path,expected=None):
    from Bio import SeqIO
    seqs={}
    for r in SeqIO.parse(str(path),'fasta'):
        if r.id in seqs:raise ValueError('Duplicate sequence '+r.id)
        s=str(r.seq).upper()
        if len(s)%3 or (expected is not None and len(s)!=3*expected):raise ValueError('Codon length mismatch: '+str(path))
        if set(s)-set('ACGT-'):raise ValueError('Unexpected character in audited codon alignment')
        for i in range(0,len(s),3):
            c=s[i:i+3]
            if '-' in c and c!='---':raise ValueError('Partial codon gap')
            if c in ('TAA','TAG','TGA'):raise ValueError('Internal stop in codon alignment')
        seqs[r.id]=s
    if not seqs or len({len(s) for s in seqs.values()})!=1:raise ValueError('Empty/ragged alignment')
    return seqs

def fasta(path,seqs):
    with open(path,'w') as f:
        for k in sorted(seqs):f.write('>'+k+'\n'+seqs[k]+'\n')



class Tasks:
    def __init__(self,out):self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True)
    def run(self,name,cmd,outputs,validator=None,timeout=None):
        folder=self.out/name;folder.mkdir(parents=True,exist_ok=True);done=folder/'DONE.json'
        outputs=[Path(p) for p in outputs]
        if done.exists():
            d=json.loads(done.read_text())
            if d['command']!=[str(x) for x in cmd] or any(not p.is_file() or sha(p)!=d['outputs'].get(str(p)) for p in outputs):raise ValueError('Checkpoint output changed: '+name)
            if validator:validator()
            print('REUSE',name,flush=True);return
        dump(folder/'COMMAND.json',[str(x) for x in cmd]);print('RUN',name,flush=True)
        start=time.time()
        with (folder/'console.log').open('a') as log:
            result=subprocess.run([str(x) for x in cmd],cwd=folder,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        if result.returncode:raise RuntimeError(f'{name}: exit={result.returncode}; see {folder}/console.log')
        if any(not p.is_file() or p.stat().st_size==0 for p in outputs):raise ValueError('Expected output missing: '+name)
        if validator:validator()
        dump(done,{'command':[str(x) for x in cmd],'outputs':fingerprint(outputs),'seconds':time.time()-start})


def check_known_hashes(root,cfg,paths):
    expected=cfg.get('expected_sha256',{})
    for path in paths:
        for rel,h in expected.items():
            if resolve(root,rel)==Path(path).resolve() and sha(path)!=h:
                raise ValueError('Input differs from the declared checksum: '+str(path)+'. Review the new source before changing expected_sha256.')


def midpoint_tree(tree):
    tr=copy.deepcopy(tree)
    while len(tr.root.clades)==1:
        tr.root=tr.root.clades[0];tr.root.branch_length=0.0
    tr.root_at_midpoint();validate_tree(tr)
    return tr
