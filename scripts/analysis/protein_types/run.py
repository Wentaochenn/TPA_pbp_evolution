"""Prepare reference-coordinate MrcA sequences and exact nucleotide haplotypes."""
import argparse
import json
from pathlib import Path
from common import read_matrix, normalize_metadata, write_json
from prepare_existing_sequences import prepare
from sequence_types import read_fasta, _coordinates, _sequence_outputs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True,help='Normalized metadata TSV or original workbook')
    p.add_argument('--settings',type=Path,help='Optional JSON with sequence and metadata adapters')
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();settings=json.loads(a.settings.read_text()) if a.settings else {}
    if a.outdir.exists():raise ValueError('Use a new sequence-typing output directory')
    matrix=read_matrix(a.matrix);metadata,audit=normalize_metadata(a.metadata,settings.get('metadata',{}))
    a.outdir.mkdir(parents=True)
    prepared=prepare(a.root,a.outdir/'prepared',matrix,settings.get('sequences',{}))
    if prepared['status']!='complete':raise ValueError('Reference-coordinate sequence preparation did not complete')
    cfg=prepared['settings'];sequences=read_fasta(cfg['alignment']);coords=_coordinates(cfg,sequences)
    sequences={sid:seq for sid,seq in sequences.items() if sid in matrix}
    if len(sequences)<2:raise ValueError('Too few matrix-matched regional sequences')
    _,haplotypes,complete=_sequence_outputs(sequences,matrix,metadata,coords,a.outdir)
    write_json(a.outdir/'RUN_STATUS.json',{'state':'complete','n_complete':complete,
               'n_nucleotide_haplotypes':len(haplotypes),'metadata_adapter':audit,
               'protein_type_translation':'Performed by the publication R code from the exact gene haplotypes'})


if __name__=='__main__':main()
