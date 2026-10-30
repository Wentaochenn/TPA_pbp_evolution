#!/usr/bin/env python3
"""Read-only rendered geometry checks. Does not redraw any figure."""
from pathlib import Path
import sys,json,subprocess
base=Path(sys.argv[1]).resolve();here=Path(__file__).resolve().parent
try:import fitz
except ImportError:
 report={'status':'NOT_AUDITABLE','reason':'PyMuPDF is unavailable; rendered checks were not run.'}
 (base/'qa/rendered_audit_summary.json').write_text(json.dumps(report,indent=2))
 print(report['reason']);sys.exit(2)
reports=[]
index=json.loads((base/'FIGURE_INDEX.json').read_text())
for p in [base/'figures'/f'{name}.pdf' for name in sorted(index)]:
 cr=subprocess.run([sys.executable,str(here/'audit_figure_collisions.py'),str(p),'--json-out',str(base/'qa'/f'{p.stem}.collision.json')],capture_output=True,text=True)
 tr=subprocess.run([sys.executable,str(here/'audit_rendered_text.py'),str(p)],capture_output=True,text=True)
 (base/'qa'/f'{p.stem}.effective_text.json').write_text(tr.stdout or tr.stderr)
 reports.append({'figure':p.stem,'collision_exit':cr.returncode,'effective_text_exit':tr.returncode})
 if cr.returncode:print(cr.stdout+cr.stderr)
 if tr.returncode:print(tr.stdout+tr.stderr)
failed=any(r['collision_exit'] or r['effective_text_exit'] for r in reports)
report={'status':'FAIL' if failed else 'PASS','figures':reports}
(base/'qa/rendered_audit_summary.json').write_text(json.dumps(report,indent=2))
print('Rendered audits:',report['status']);sys.exit(bool(failed))
