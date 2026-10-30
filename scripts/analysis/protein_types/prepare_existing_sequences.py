#!/usr/bin/env python3
"""Extract the observed study region from existing reference-coordinate inputs."""
import ast
from collections import Counter
import json
from pathlib import Path
import re

from common import sha256, table, write_json
ALLOWED = frozenset("ACGTRYSWKMBDHVN?-")

# Corrected coordinates in the project's p3_05_hyphy_selection_tests.py.
# Do not use the legacy gene_coordinates.tsv mrcA label: that row is pbp2.
SS14_COORDINATES = (772313, 774967, '-', 884)
NICHOLS_COORDINATES = (772433, 775087, '-', 884)
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN?-', 'TGCAYRSWMKVHDBN?-')
BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a+b+c: AA[i*16+j*4+k] for i,a in enumerate(BASES)
        for j,b in enumerate(BASES) for k,c in enumerate(BASES)}
AA3 = {'Ala':'A', 'Thr':'T', 'Val':'V', 'Met':'M', 'Gly':'G', 'Ser':'S',
       'Arg':'R', 'His':'H', 'Gln':'Q', 'Pro':'P', 'Leu':'L', 'Ile':'I',
       'Asp':'D', 'Glu':'E', 'Phe':'F', 'Tyr':'Y', 'Trp':'W', 'Cys':'C',
       'Asn':'N', 'Lys':'K', 'Ter':'*'}


def records(path):
    """One full genome at a time, not the multi-gigabyte alignment in memory."""
    seen, sid, pieces = set(), None, []
    with Path(path).open(encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sid is not None:
                    sequence = ''.join(pieces)
                    if not sequence:
                        raise ValueError('Empty FASTA sequence: ' + sid)
                    yield sid, sequence
                fields = line[1:].split()
                if not fields or fields[0] in seen:
                    raise ValueError('Empty or duplicate FASTA identifier')
                sid, pieces = fields[0], []
                seen.add(sid)
            else:
                if sid is None or set(line.upper()) - ALLOWED:
                    raise ValueError('Malformed nucleotide FASTA')
                pieces.append(line.upper())
    if sid is not None:
        sequence = ''.join(pieces)
        if not sequence:
            raise ValueError('Empty FASTA sequence: ' + sid)
        yield sid, sequence


def revcomp(value):
    return value.translate(COMPLEMENT)[::-1]


def protein(value):
    return ''.join(CODE.get(value[i:i+3], 'X') for i in range(0, len(value)-2, 3))


def _path(root, raw):
    path = Path(raw).expanduser()
    return (root/path).resolve() if not path.is_absolute() else path.resolve()


def _yaml_scalar(path, key):
    if not path.is_file():
        return None
    found = []
    for line in path.read_text(encoding='utf-8').splitlines():
        match = re.match(r'^' + re.escape(key) + r'\s*:\s*(.+?)\s*$', line)
        if match:
            value = match[1].split(' #', 1)[0].strip()
            if value[:1] in ('"', "'") and value[-1:] == value[:1]:
                value = value[1:-1]
            if not value or value.startswith(('[', '{', '&', '*', '!')):
                raise ValueError('Unsupported YAML scalar for ' + key)
            found.append(value)
    if len(found) > 1:
        raise ValueError('Duplicate YAML input key: ' + key)
    return found[0] if found else None


def validate_settings(settings):
    if not isinstance(settings, dict) or set(settings) - {'alignment', 'reference', 'nichols_reference', 'flank', 'max_discordance_fraction'}:
        raise ValueError('Invalid existing-sequence settings')
    if type(settings.get('flank', 500)) is not int or settings.get('flank', 500) < 0:
        raise ValueError('flank must be a nonnegative integer')
    fraction = settings.get('max_discordance_fraction', 0.05)
    if isinstance(fraction, bool) or not isinstance(fraction, (float,int)) or not 0 <= fraction <= 1:
        raise ValueError('max_discordance_fraction must be between 0 and 1')
    for field in ('alignment', 'reference', 'nichols_reference'):
        if field in settings and (not isinstance(settings[field], str) or not settings[field]):
            raise ValueError(field + ' must be a nonempty path')


def prepare(root, outdir, matrix, settings):
    """Return local-module settings and the complete input provenance chain."""
    validate_settings(settings)
    root, outdir = Path(root).resolve(), Path(outdir)
    if outdir.exists():
        raise ValueError('Preparation output directory must be new')
    config = root/'config.yaml'
    if settings.get('alignment'):
        alignment = _path(root, settings['alignment'])
        authority = 'explicit_preparation_setting'
    elif _yaml_scalar(config, 'input_msa'):
        alignment = _path(root, _yaml_scalar(config, 'input_msa'))
        authority = 'project_config_yaml_input_msa'
    else:
        choices = [root/p for p in ('consensus/TPA_ss14ref_msa.fasta', 'msa.fasta',
                                   'consensus/msa.fasta', 'results/00_msa/msa.fasta') if (root/p).is_file()]
        if len(choices) > 1:
            raise ValueError('Ambiguous full alignments; declare input_msa in the existing project config')
        alignment = choices[0] if choices else root/'consensus/TPA_ss14ref_msa.fasta'
        authority = 'established_pipeline_path'
    reference = _path(root, settings.get('reference') or _yaml_scalar(config, 'reference') or 'refs/TPA_SS14.fasta')
    if not alignment.is_file() or not reference.is_file():
        return {'status': 'blocked', 'settings': None, 'inputs': [str(p) for p in (config,alignment,reference) if p.is_file()],
                'audit': {'reason': 'Existing project inputs unavailable in this environment',
                          'alignment_expected': str(alignment), 'reference_expected': str(reference)}}
    annotation = root/'script/p3_05_hyphy_selection_tests.py'
    nichols = _path(root, settings.get('nichols_reference', 'refs/TPA_Nichols.fasta'))
    if not nichols.is_file():
        raise FileNotFoundError('Nichols reference required for cross-reference numbering audit: ' + str(nichols))
    nichols_annotation = root/'script/p2_01_extract_codon_corrected_variants.py'
    source_inputs = [alignment,reference,nichols]+([config] if config.is_file() else [])+[
        p for p in (annotation,nichols_annotation) if p.is_file()]
    initial_stats = {str(p): (p.stat().st_size,p.stat().st_mtime_ns) for p in source_inputs}
    initial_hashes = {str(p):sha256(p) for p in source_inputs}
    forbidden = ('filtered_polymorphic', '1294sites', 'variable_only', 'snp_only.fasta', 'masked.fasta',
                 'mrcA_seq_data', 'theoretical', 'mrcA_codon')
    if any(token.lower() in str(alignment).lower() for token in forbidden):
        raise ValueError('Do not use screened SNP columns, masked/codon output or theoretical constructs as raw region input')
    ref_records = list(records(reference))
    if len(ref_records) != 1:
        raise ValueError('Expected one SS14 chromosome reference')
    ref_id, genome = ref_records[0]
    # Unambiguous bases are required in the region actually used downstream.
    # IUPAC ambiguity elsewhere retains chromosome positions and is irrelevant
    # to this extraction; alignment gaps never define chromosome coordinates.
    if not genome or set(genome) - set('ACGTRYSWKMBDHVN'):
        raise ValueError('Chromosome reference must be nonempty, gap-free nucleotide sequence in forward coordinates')
    if 'NC_021490' in ref_id:
        raise ValueError('Nichols reference cannot anchor an SS14-coordinate MSA')
    coordinates = SS14_COORDINATES
    if annotation.is_file():
        tree = ast.parse(annotation.read_text(encoding='utf-8'))
        definitions = [node for node in tree.body if isinstance(node, ast.Assign) and
                       any(isinstance(t, ast.Name) and t.id == 'GENES' for t in node.targets)]
        if len(definitions) != 1:
            raise ValueError('Ambiguous or missing corrected project GENES annotation')
        parsed = ast.literal_eval(definitions[0].value)
        declared = tuple(parsed.get('mrcA', ()))
        if declared != coordinates:
            raise ValueError('Corrected project mrcA coordinates differ from verified release; stop for reference review')
    start, end, strand, expected = coordinates
    if not 1 <= start <= end <= len(genome):
        raise ValueError('Gene coordinates fall outside the reference chromosome')
    flank = settings.get('flank', 500)
    lo, hi = max(1,start-flank), min(len(genome),end+flank)
    reference_ambiguity = {base:count for base,count in Counter(genome).items() if base not in 'ACGT'}
    region_ambiguity = {base:count for base,count in Counter(genome[lo-1:hi]).items() if base not in 'ACGT'}
    if region_ambiguity:
        raise ValueError('Reference extraction interval '+str(lo)+'-'+str(hi)+
                         ' must be gap-free ACGT; ambiguous counts: '+str(region_ambiguity))
    cds = revcomp(genome[start-1:end])
    translated = protein(cds)
    if len(cds) != (expected+1)*3 or len(translated) != expected+1 or translated[-1] != '*' or '*' in translated[:-1]:
        raise ValueError('Reference mrcA length/frame/terminal-stop validation failed; reject coordinate mismatch')
    residues = {str(p): translated[p-1] for p in (506,625,708)}
    if residues['506'] not in 'ATVM' or residues['625'] not in 'MV' or residues['708'] not in 'GS':
        raise ValueError('Reference key-residue validation failed; reject the legacy pbp2/mrcA label mismatch')
    if nichols_annotation.is_file():
        definitions = [node for node in ast.parse(nichols_annotation.read_text(encoding='utf-8')).body
                       if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and
                       t.id == 'RESISTANCE_GENES' for t in node.targets)]
        if len(definitions) != 1:
            raise ValueError('Ambiguous or missing Nichols RESISTANCE_GENES annotation')
        declared = ast.literal_eval(definitions[0].value).get('mrcA', ())
        if tuple(declared[1:4]) != NICHOLS_COORDINATES[:3]:
            raise ValueError('Nichols project mrcA coordinates differ; stop for reference review')
    from reference_numbering import reference_audit
    numbering, proteins = reference_audit(reference, nichols, outdir/'reference_numbering',
        {'mrcA': {'SS14': SS14_COORDINATES, 'Nichols': NICHOLS_COORDINATES}})
    rows, concordance, seqs = [], [], {}
    patterns = []
    for locus in next(iter(matrix.values())):
        match = re.fullmatch(r'mrcA:p\.([A-Za-z]+)(\d+)([A-Za-z]+)', locus)
        if locus.startswith('mrcA:') and (not match or match[1] not in AA3 or match[3] not in AA3):
            raise ValueError('Unsupported mrcA protein indicator: ' + locus)
        if match:
            position = int(match[2])
            if not 1 <= position <= expected or proteins['Nichols']['mrcA'][position-1] != AA3[match[1]]:
                raise ValueError('Matrix reference amino acid differs from Nichols reference: ' + locus)
            patterns.append((locus, int(match[2]), AA3[match[3]]))
    n_checked = n_discordant = 0
    for sid, seq in records(alignment):
        if len(seq) != len(genome):
            raise ValueError('Alignment not in full reference coordinates: ' + sid)
        if sid not in matrix:
            rows.append({'sample':sid, 'scope':'alignment_outside_current_matrix'})
            continue
        region = list(seq[lo-1:hi])
        observed_cds = revcomp(seq[start-1:end])
        aa = protein(observed_cds)
        mask_positions = set()
        for locus, position, alt in patterns:
            if not 1 <= position <= expected:
                raise ValueError('Matrix mrcA codon outside corrected reference gene')
            state = matrix[sid][locus]
            observed = aa[position-1]
            reason = ''
            if state is None:
                mask_positions.add(position)
                reason = 'matrix_full_codon_not_callable'
            elif observed == 'X':
                reason = 'source_sequence_uncertain'
            else:
                n_checked += 1
                if int(observed == alt) != state:
                    n_discordant += 1
                    mask_positions.add(position)
                    reason = 'source_sequence_matrix_discordance_masked'
            if reason:
                concordance.append({'sample':sid,'locus':locus,'observed_aa':observed,
                                    'matrix_state':state,'action':reason})
        for position in mask_positions:
            # Negative-strand codon [end-3p+1, end-3p+3], in forward region.
            codon_lo = end-3*position+1
            for genomic in range(codon_lo, codon_lo+3):
                region[genomic-lo] = 'N'
        seqs[sid] = ''.join(region)
    for sid in sorted(set(matrix)-set(seqs)):
        rows.append({'sample':sid,'scope':'matrix_without_source_sequence'})
    if len(seqs) < 2:
        raise ValueError('Insufficient exact matrix/whole-alignment sample overlap')
    fraction = n_discordant/n_checked if n_checked else 0
    table(outdir/'source_alignment_sample_audit.tsv', rows, ['sample','scope'])
    table(outdir/'source_matrix_codon_concordance.tsv',concordance,
          ['sample','locus','observed_aa','matrix_state','action'])
    if fraction > settings.get('max_discordance_fraction',0.05):
        raise ValueError('Excessive source/matrix discordance; review coordinates, consensus source and sample IDs')
    ref_name = 'REFERENCE_NC_021508_1'
    if ref_name in matrix:
        raise ValueError('Reserved reference identifier collides with study sample')
    with (outdir/'region.fasta').open('w') as handle:
        handle.write('>'+ref_name+'\n'+genome[lo-1:hi]+'\n')
        for sid, seq in sorted(seqs.items()):
            handle.write('>'+sid+'\n'+seq+'\n')
    if any((p.stat().st_size,p.stat().st_mtime_ns) != initial_stats[str(p)] or sha256(p) != initial_hashes[str(p)]
           for p in source_inputs):
        raise ValueError('Existing sequence source changed during preparation')
    local = {'alignment':str((outdir/'region.fasta').resolve()),'reference_sample':ref_name,
             'coordinate_system':'reference_1based_inclusive_ungapped','region_start':lo,
             'gene_start':start,'gene_end':end,'gene_strand':strand}
    audit = {'status':'complete','alignment':str(alignment),'selection_authority':authority,
             'reference':str(reference),'reference_id':ref_id,'reference_residues':residues,
             'cross_reference_numbering':numbering,
             'reference_ambiguity_outside_used_interval':reference_ambiguity,
             'reference_used_interval_ambiguity':region_ambiguity,
             'reference_validation_scope':'Gap-free chromosome positions retained; extracted gene and flanks require ACGT; outside-region IUPAC ambiguity allowed without imputation',
             'gene_coordinates':coordinates,'n_selected_samples':len(seqs),'region_start':lo,'region_end':hi,
             'source_mask_policy':'Preserve all input N/IUPAC/gaps; upstream consensus can be DP3-masked',
             'tested_codon_policy':'Matrix noncallable/discordant codons masked N; no other bases imputed',
             'n_compared_indicator_calls':n_checked,'n_discordant_indicator_calls':n_discordant,
             'discordance_fraction':fraction,
             'sequence_scope':'reference-coordinate SNP-consensus substitution backgrounds; indels not represented',
             'inputs_sha256':initial_hashes}
    write_json(outdir/'SEQUENCE_PREPARATION.json',audit)
    write_json(outdir/'local_module_settings.json',local)
    return {'status':'complete','settings':local,'inputs':[str(p) for p in source_inputs],'audit':audit}
