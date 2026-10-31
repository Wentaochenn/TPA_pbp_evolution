"""Generate the structural figure inputs from the frozen reference coordinates."""
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MODEL_SHA256 = '8a13293e93d580c94a59cddd6cb56587106a4947a18d5c7fe21ae6b6d1327587'


def sha(path):
    from prepare import sha as digest
    return digest(path)


def verify_render(folder, pdb, width):
    from PIL import Image
    record = json.loads((folder / 'RENDER_STATUS.json').read_text())
    if record.get('state') != 'complete' or record.get('pdb_sha256') != sha(pdb):
        raise ValueError('Incomplete or mismatched structural rendering')
    if record.get('script_sha256') != sha(HERE / 'structure/render_reference.py'):
        raise ValueError('Structural rendering script differs')
    if record.get('camera_preset_sha256') != sha(HERE / 'structure/camera_presets.json'):
        raise ValueError('Structural camera preset differs')
    if record.get('original_image_assets_used') is not False:
        raise ValueError('Structure inputs must be generated from PDB coordinates')
    for stem, height in [('MrcA_reference_structure', round(width * 0.625)),
                         ('MrcA_local_context', width)]:
        with Image.open(folder / (stem + '.png')) as image:
            if image.size != (width, height):
                raise ValueError('Unexpected structural image dimensions: ' + stem)
            image.verify()
        if not (folder / (stem + '.pse')).is_file():
            raise ValueError('Missing saved PyMOL session: ' + stem)
    record['vector_label_validation']=verify_label_overlay(folder,record,width)
    return record



def verify_label_overlay(folder,record,width):
    """Verify the vector labels against the actual saved coordinates and camera."""
    if record.get('raster_labels') is not False or record.get('vector_label_font')!='Arial':
        raise ValueError('Require Arial vector labels without raster text')
    path=folder/'MrcA_label_positions.tsv'
    if not path.is_file() or record.get('vector_label_sha256')!=sha(path):
        raise ValueError('Structural label coordinates missing or changed')
    spec=importlib.util.spec_from_file_location('_arial_renderer',HERE/'structure/render_reference.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    with (folder/'MrcA_marked_residue_coordinates.tsv').open(newline='') as f:
        marked=list(csv.DictReader(f,delimiter='\t'))
    views=json.loads((folder/'camera_views.json').read_text())
    presets=json.loads((HERE/'structure/camera_presets.json').read_text())['views']
    for name,wanted in presets.items():
        actual=views.get(name,[])
        if len(actual)!=18 or not all(math.isclose(x,y,abs_tol=1e-5) for x,y in zip(actual,wanted)):
            raise ValueError('Actual camera differs from fixed preset')
    expected=module._build_label_overlay(marked,views['MrcA_local_context'],width,width)
    with path.open(newline='') as f:actual=list(csv.DictReader(f,delimiter='\t'))
    if len(actual)!=3:raise ValueError('Expected three vector labels')
    for got,wanted in zip(actual,expected):
        if set(got)!=set(wanted):raise ValueError('Unexpected label columns')
        for key,value in wanted.items():
            equal=(math.isclose(float(got[key]),value,rel_tol=1e-10,abs_tol=1e-10)
                   if isinstance(value,(int,float)) else got[key]==value)
            if not equal:raise ValueError('Label geometry differs: '+key)
    return dict(sites=[506,625,708],font='Arial',geometry_verified=True)


def ensure_structure(run, python=None, environment='pbp_structure_render', source=None, width=2400):
    from PIL import Image
    run = Path(run)
    data = run / 'work/figures/data'
    pdb = data / 'mrcA_model.pdb'
    if sha(pdb) != MODEL_SHA256:
        raise ValueError('Reference model differs from the frozen publication model')
    dest = run / 'Structure_Render'
    if source:
        source = Path(source).resolve()
        verify_render(source, pdb, width)
        shutil.copytree(source, dest)
    else:
        if python:
            command = [str(python)]
        else:
            conda = shutil.which('conda') or os.environ.get('CONDA_EXE')
            if not conda:
                raise ValueError('Provide --pymol-python or make conda available')
            command = [conda, 'run', '--no-capture-output', '-n', environment, 'python']
        command += [str(HERE / 'structure/render_reference.py'), '--pdb', str(pdb),
                    '--outdir', str(dest), '--width', str(width)]
        print('[RUN] PDB structure rendering; progress in structure_render.console.log', flush=True)
        with (run / 'structure_render.console.log').open('w') as output:
            subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=True)
    record = verify_render(dest, pdb, width)
    for stem in ['MrcA_reference_structure', 'MrcA_local_context']:
        shutil.copy2(dest / (stem + '.png'), data / (stem + '.png'))
        with Image.open(dest / (stem + '.png')) as image:
            image.convert('RGB').save(data / (stem + '.ppm'))
    for src, name in [('camera_views.json', 'structure_camera_views.json'),
                      ('RENDER_STATUS.json', 'structure_render_provenance.json'),
                      ('MrcA_marked_residue_coordinates.tsv', 'structure_marked_residues.tsv'),
                      ('MrcA_label_positions.tsv', 'structure_label_overlay.tsv')]:
        shutil.copy2(dest / src, data / name)
    # Remove superseded assets only from this new working copy.
    for path in data.glob('original_structure_*'):
        if path.is_file():
            path.unlink()
    record['reused_render'] = str(source) if source else None
    return record
