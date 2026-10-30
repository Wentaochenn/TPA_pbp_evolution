# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = ["pymol-open-source-whl==3.2.0.2"]
# ///
"""Render the frozen MrcA reference model and three previously reported sites."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import traceback

EXPECTED_SHA256 = '8a13293e93d580c94a59cddd6cb56587106a4947a18d5c7fe21ae6b6d1327587'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pdb',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--width',type=int,default=2400)
    p.add_argument('--camera-presets',type=Path,default=Path(__file__).with_name('camera_presets.json'))
    a=p.parse_args()
    presets=json.loads(a.camera_presets.read_text())
    if presets.get('pdb_sha256')!=EXPECTED_SHA256:p.error('Camera preset model mismatch')
    for name in ['MrcA_reference_structure','MrcA_local_context']:
        view=presets['views'][name]
        if len(view)!=18 or not all(isinstance(v,(int,float)) and math.isfinite(v) for v in view):
            p.error('Invalid camera view: '+name)
    pdb=a.pdb.resolve(strict=True)
    if not pdb.is_file():p.error('PDB is not a regular file')
    sha=hashlib.sha256(pdb.read_bytes()).hexdigest()
    if sha!=EXPECTED_SHA256:p.error('PDB differs from the frozen MrcA reference model')
    out=a.outdir.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'RENDER_STATUS.json').exists():p.error('Use a new output directory')
    os.environ['PYOPENGL_PLATFORM']='osmesa'
    os.environ.pop('DISPLAY',None)
    import pymol
    pymol.pymol_argv=['pymol','-cq']
    pymol.finish_launching()
    from pymol import cmd
    failed=False
    try:
        cmd.load(str(pdb),'mrcA')
        atoms=cmd.count_atoms('all')
        if not atoms:raise RuntimeError('PyMOL loaded zero atoms')
        ca_count=cmd.count_atoms('mrcA and name CA')
        if ca_count!=884:raise ValueError('Expected 884 reference C-alpha atoms')
        if ca_count==atoms:cmd.set('cartoon_trace_atoms',1)
        cmd.hide('everything','all');cmd.show('cartoon','mrcA')
        cmd.set_color('reference_grey',[0.81,0.83,0.84])
        cmd.set_color('site_red',[0.698,0.094,0.169])
        cmd.set_color('site_orange',[0.902,0.49,0.09])
        cmd.color('reference_grey','mrcA')
        cmd.set('cartoon_transparency',0.15)
        cmd.set('cartoon_fancy_helices',1)
        cmd.set('cartoon_fancy_sheets',1)
        cmd.set('orthoscopic',1);cmd.set('depth_cue',0);cmd.set('ray_shadows',0)
        cmd.set('ray_opaque_background',1);cmd.set('opaque_background',1)
        cmd.set('antialias',2);cmd.set('specular',0.25)
        cmd.set('max_threads',max(1,int(os.environ.get('NSLOTS','1'))))
        cmd.bg_color('white')
        marked=[]
        expected={506:'ALA',625:'MET',708:'GLY'}
        for site,resn in expected.items():
            sel='mrcA and resi {} and name CA'.format(site)
            model=cmd.get_model(sel)
            if len(model.atom)!=1 or model.atom[0].resn!=resn:
                raise ValueError('Reference residue identity/uniqueness mismatch at '+str(site))
            atom=model.atom[0]
            marked.append({'residue':resn,'position':site,'chain':atom.chain,'x_A':atom.coord[0],
                           'y_A':atom.coord[1],'z_A':atom.coord[2],'pLDDT':atom.b,
                           'colour':'red' if site==506 else 'orange','marker':'C-alpha sphere'})
            cmd.color('site_red' if site==506 else 'site_orange','mrcA and resi '+str(site))
        cmd.select('marked_ca','mrcA and name CA and resi 506+625+708')
        cmd.show('spheres','marked_ca');cmd.set('sphere_scale',1.0,'marked_ca')
        cmd.select('marked_residues','mrcA and resi 506+625+708')
        cmd.show('sticks','marked_residues');cmd.set('stick_radius',0.18,'marked_residues')
        cmd.set('label_color','black');cmd.set('label_size',10);cmd.set('label_font_id',7)
        cmd.set('float_labels',1);cmd.set('label_outline_color','white')
        for label,site in [('A506',506),('M625',625),('G708',708)]:
            cmd.label('mrcA and name CA and resi '+str(site),repr(label))
        for site,position in [(506,[-5.5,-3.5,3.0]),(625,[4.0,1.0,3.0]),(708,[-4.0,3.0,3.0])]:
            cmd.set('label_position',position,'mrcA and name CA and resi '+str(site))
        views={}
        for name in ['MrcA_reference_structure','MrcA_local_context']:
            cmd.hide('labels','all')
            cmd.set('cartoon_transparency',0.20 if name=='MrcA_reference_structure' else 0.55)
            if name=='MrcA_local_context':cmd.show('labels','marked_ca')
            cmd.set_view(presets['views'][name])
            views[name]=list(cmd.get_view())
            cmd.png(str(out/(name+'.png')),width=a.width,height=(round(a.width*0.625) if name=='MrcA_reference_structure' else a.width),dpi=300,ray=1,quiet=0)
            if not (out/(name+'.png')).is_file():raise RuntimeError('PNG was not written')
            cmd.save(str(out/(name+'.pse')))
        with (out/'MrcA_marked_residue_coordinates.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(marked[0]),delimiter='\t');w.writeheader();w.writerows(marked)
        (out/'camera_views.json').write_text(json.dumps(views,indent=2)+'\n')
        provenance={'state':'complete','pdb':str(pdb),'pdb_sha256':sha,'pymol_version':list(cmd.get_version()),
                    'renderer':'PyMOL CPU ray tracer through cmd.png(ray=1)',
                    'original_image_assets_used':False,'mutant_model_generated':False,
                    'camera_preset_sha256':hashlib.sha256(a.camera_presets.read_bytes()).hexdigest(),
                    'camera_provenance':presets['view_provenance'],
                    'marked_positions':[506,625,708],'image_width_px':a.width,
                    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (out/'RENDER_STATUS.json').write_text(json.dumps(provenance,indent=2,default=str)+'\n')
        print('[DONE] '+str(out),flush=True)
    except BaseException as exc:
        failed=True;traceback.print_exc()
        (out/'RENDER_STATUS.json').write_text(json.dumps({'state':'failed','error':str(exc)},indent=2)+'\n')
    finally:
        cmd.quit(1 if failed else 0)


if __name__=='__main__':main()
