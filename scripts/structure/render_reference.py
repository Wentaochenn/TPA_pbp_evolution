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



def _build_label_overlay(marked, view, width, height):
    """Project C-alpha markers to bottom-left normalized orthographic image space.

    PyMOL get_view rotation is column-major: camera = R.T*(p-origin)+translation.
    Half-height = -view[11]*tan(field_of_view/2), half-width includes image aspect.
    Text and leaders are drawn by R; no labels are baked into the raster.
    """
    if len(view)!=18 or not all(math.isfinite(float(v)) for v in view):
        raise ValueError('Invalid label camera')
    if width<=0 or height<=0 or view[11]>=0 or abs(abs(view[17])-20)>1e-5:
        raise ValueError('Expected fixed orthographic 20-degree camera')
    for i in range(3):
        for j in range(3):
            dot=sum(view[3*k+i]*view[3*k+j] for k in range(3))
            if abs(dot-(1 if i==j else 0))>1e-4:raise ValueError('Nonorthogonal camera')
    expected={506:('ALA','A506'),625:('MET','M625'),708:('GLY','G708')}
    by_site={int(row['position']):row for row in marked}
    if len(marked)!=3 or set(by_site)!=set(expected):raise ValueError('Expected three marked sites')
    hh=-view[11]*math.tan(math.radians(20)/2);hw=hh*width/height
    offsets={506:(-.13,-.04),625:(.11,-.02),708:(-.08,.055)}
    rows=[]
    for site,(resn,label) in expected.items():
        atom=by_site[site]
        if atom['residue']!=resn:raise ValueError('Reference residue differs')
        p=[float(atom[a+'_A']) for a in ('x','y','z')]
        rel=[p[k]-view[12+k] for k in range(3)]
        cam=[sum(view[3*k+i]*rel[k] for k in range(3))+view[9+i] for i in range(3)]
        mx,my=.5+cam[0]/(2*hw),.5+cam[1]/(2*hh)
        dx,dy=offsets[site];lx,ly=mx+dx,my+dy
        if not all(.03<v<.97 for v in (mx,my,lx,ly)):raise ValueError('Label outside camera')
        length=math.hypot(dx*2*hw,dy*2*hh)
        rows.append(dict(view='MrcA_local_context',position=site,residue_label=label,
            marker_x=mx,marker_y=my,label_x=lx,label_y=ly,
            leader_x=lx+(.012 if dx<0 else -.012),leader_y=ly,
            leader_xend=mx+dx*1.7/length,leader_yend=my+dy*1.7/length,
            hjust=1 if dx<0 else 0,font_family='Arial',font_face='bold',
            coordinate_system='normalized_image_bottom_left',image_width_px=width,image_height_px=height))
    return rows


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
        # R draws editable Arial labels; the reference PNG has no text.
        cmd.hide('labels','all')
        views={}
        for name in ['MrcA_reference_structure','MrcA_local_context']:
            cmd.hide('labels','all')
            cmd.set('cartoon_transparency',0.20 if name=='MrcA_reference_structure' else 0.55)
            cmd.set_view(presets['views'][name])
            cmd.set('orthoscopic',1);cmd.set('field_of_view',20)
            views[name]=list(cmd.get_view())
            cmd.png(str(out/(name+'.png')),width=a.width,height=(round(a.width*0.625) if name=='MrcA_reference_structure' else a.width),dpi=300,ray=1,quiet=0)
            if not (out/(name+'.png')).is_file():raise RuntimeError('PNG was not written')
            cmd.save(str(out/(name+'.pse')))
        with (out/'MrcA_marked_residue_coordinates.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(marked[0]),delimiter='\t');w.writeheader();w.writerows(marked)
        (out/'camera_views.json').write_text(json.dumps(views,indent=2)+'\n')
        overlay=_build_label_overlay(marked,views['MrcA_local_context'],a.width,a.width)
        overlay_path = out/'MrcA_label_positions.tsv'
        with overlay_path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(overlay[0]),delimiter='\t');w.writeheader();w.writerows(overlay)
        provenance={'state':'complete','pdb':str(pdb),'pdb_sha256':sha,'pymol_version':list(cmd.get_version()),
                    'renderer':'PyMOL CPU ray tracer through cmd.png(ray=1)',
                    'original_image_assets_used':False,'mutant_model_generated':False,
                    'camera_preset_sha256':hashlib.sha256(a.camera_presets.read_bytes()).hexdigest(),
                    'camera_provenance':presets['view_provenance'],
                    'marked_positions':[506,625,708],'image_width_px':a.width,
                    'raster_labels':False,'vector_label_font':'Arial',
                    'vector_label_sha256':hashlib.sha256(overlay_path.read_bytes()).hexdigest(),
                    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (out/'RENDER_STATUS.json').write_text(json.dumps(provenance,indent=2,default=str)+'\n')
        print('[DONE] '+str(out),flush=True)
    except BaseException as exc:
        failed=True;traceback.print_exc()
        (out/'RENDER_STATUS.json').write_text(json.dumps({'state':'failed','error':str(exc)},indent=2)+'\n')
    finally:
        cmd.quit(1 if failed else 0)


if __name__=='__main__':main()
