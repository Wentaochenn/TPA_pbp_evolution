"""Check the installed publication runtime before structural rendering and export."""
import argparse
import importlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import os


def check(package, rscript):
    report={'python':sys.version,'python_packages':{},'R':None,'software_installed':False,
            'locale':{k:os.environ.get(k) for k in ['LANG','LC_CTYPE','LC_NUMERIC','LC_COLLATE']}}
    for module,distribution in [('numpy','numpy'),('scipy','scipy'),('PIL','Pillow'),('fitz','PyMuPDF')]:
        importlib.import_module(module)
        report['python_packages'][distribution]=importlib.metadata.version(distribution)
    exe=shutil.which(rscript)
    if not exe:raise RuntimeError('Rscript is unavailable; activate pbp_figures_release before publishing')
    expression='''p <- c("ggplot2","patchwork","dplyr","tidyr","ggrepel","ape","jsonlite","maps","svglite","ragg","systemfonts","textshaping","gtable","scales","openxlsx")
missing <- p[!vapply(p, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if(length(missing)) stop(paste("Missing R packages:",paste(missing,collapse=", ")))
cat(R.version.string,"\\n")
for(n in p) cat(n,as.character(packageVersion(n)),"\\n")
for(f in list.files(commandArgs(TRUE)[1],pattern="[.]R$",recursive=TRUE,full.names=TRUE)) parse(file=f)
'''
    result=subprocess.run([exe,'--vanilla','-e',expression,str(Path(package).resolve())],
                          capture_output=True,text=True)
    if result.returncode:
        raise RuntimeError('R runtime or syntax preflight failed:\n'+result.stdout+result.stderr)
    report['R']=result.stdout.strip().splitlines()
    report['state']='PASS'
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,required=True)
    parser.add_argument('--rscript',default='Rscript')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    try:
        record=check(args.package,args.rscript)
    except Exception as error:
        record={'state':'FAIL','error':type(error).__name__+': '+str(error),'software_installed':False}
        args.out.write_text(json.dumps(record,indent=2)+'\n')
        raise
    args.out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
