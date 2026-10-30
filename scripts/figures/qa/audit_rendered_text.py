#!/usr/bin/env python3
"""Effective PDF font sizes (including Tm/CTM scaling used by R pdf())."""
import json,sys
from pathlib import Path
import fitz
p=Path(sys.argv[1]);rows=[]
with fitz.open(p) as doc:
 for i,page in enumerate(doc,1):
  for b in page.get_text('dict')['blocks']:
   if b['type']!=0:continue
   for line in b['lines']:
    for s in line['spans']:
     if s['text'].strip():rows.append({'page':i,'text':s['text'],'size_pt':s['size']})
failed=[r for r in rows if r['size_pt']<5-1e-5]
report={'pdf':str(p),'method':'PyMuPDF effective text size including PDF text matrix','minimum_pt':min(r['size_pt'] for r in rows),'n_spans':len(rows),'below_5pt':failed,'status':'FAIL' if failed else 'PASS'}
print(json.dumps(report,indent=2));sys.exit(bool(failed))
