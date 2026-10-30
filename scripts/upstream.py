"""Prepare an isolated upstream project or run one documented analysis stage."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
STAGES=['check','init','wgs','p0','p0_02','p1','p2','p3','selection','temporal',
        'treetime','beast_prepare','beast_run','pbp_temporal','esm2','structural','figures',
        'all_until_beauti','all_after_beauti','all','all_from_p0']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--source-project',type=Path,required=True);q.add_argument('--run-root',type=Path,required=True)
    q=sub.add_parser('run');q.add_argument('stage',choices=STAGES);q.add_argument('--project',type=Path,required=True)
    a=p.parse_args()
    if a.action=='prepare':
        command=[sys.executable,str(HERE/'upstream/bin/prepare_run.py'),'--source-project',str(a.source_project.resolve()),'--run-root',str(a.run_root.resolve())]
        subprocess.run(command,check=True)
    else:
        project=a.project.resolve()
        if not (project/'run_inputs.json').is_file() or (project/'script').resolve()!=(HERE/'upstream/script').resolve():
            raise ValueError('Use master.sh upstream prepare to create an isolated project linked to this code release')
        env=dict(os.environ,PROJECT_ROOT=str(project),SCRIPT_DIR=str(HERE/'upstream/script'),PBP_PUBLICATION_PACKAGE=str(HERE))
        subprocess.run(['bash',str(HERE/'upstream/master.sh'),a.stage],cwd=project,env=env,check=True)


if __name__=='__main__':main()
