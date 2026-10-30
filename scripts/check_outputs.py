#!/usr/bin/env python3
"""Check exported figures, panel data and study tables without drawing figures."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

FIGURES = {
    'Fig1_dataset_and_variants', 'Fig2_lineages_and_recurrence',
    'Fig3_selection', 'Fig4_mrcA_backgrounds', 'Fig5_time_calibrated_histories',
    'Supp1_reconstruction_sensitivity', 'Supp3_ESM2',
    'Supp4_temporal_regressions', 'Supp2_structural_context',
}
PANELS = {('Figure' + str(i), p) for i, letters in
          [(1, 'ABCDEF'), (2, 'ABCDE'), (3, 'A'), (4, 'ABC'), (5, 'ABCD')]
          for p in letters}
PANELS |= {('FigureS' + str(i), p) for i, letters in
           [(1, 'A'), (2, 'ABCD'), (3, 'AB'), (4, 'A')] for p in letters}
GENES = {'mrcA': 884, 'pbp1': 624, 'pbp2': 652, 'tp47': 434}


def substituted_glyphs(text):
    return '..' in text or '\ufffd' in text


def table(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t', escapechar='\\')
        rows = list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError('Malformed TSV: ' + str(path))
        return reader.fieldnames, rows


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def validate_carrier_year_rows(rows):
    """Reject assigned dates for zero-count layout rows and missing carrier dates."""
    for row in rows:
        count = int(row['n'])
        value = row.get('year')
        missing = value is None or str(value).strip().lower() in ('', 'na', 'nan')
        if count < 0:
            raise ValueError('Carrier counts cannot be negative')
        if count == 0 and not missing:
            raise ValueError('A zero-carrier category has an assigned calendar year')
        if count > 0 and (missing or not math.isfinite(float(value))):
            raise ValueError('A dated-carrier count has no finite calendar year')


def validate_metadata_exports(run):
    """Cross-check the formal metadata table against the frozen audited states."""
    from supplementary_metadata import (
        EXCLUDED_COLUMNS, build_table, read_tsv,
    )
    run = Path(run)
    full = run / 'Reproducibility_Data/Analysis_Tables/sample_metadata.tsv'
    _, metadata = read_tsv(full)
    matrix_fields, matrix = read_tsv(run / 'work/figures/data/matrix.tsv')
    expected_fields, expected_rows, spec = build_table(metadata, matrix, matrix_fields)
    fields, actual = read_tsv(run / 'Supplementary_Tables/TableS1_sample_metadata.tsv')
    if fields != expected_fields or actual != expected_rows:
        raise ValueError('Supplementary Table 1 differs from the audited metadata and codon states')
    if set(fields).intersection(EXCLUDED_COLUMNS):
        raise ValueError('Removed analysis-membership fields remain in Supplementary Table 1')
    if len(actual) != 1861 or len(fields) != 54:
        raise ValueError('Expected 1861 genomes and 54 formal metadata columns')

    path = run / 'Supplementary_Tables/TableS1_sample_metadata.xlsx'
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with ZipFile(path) as archive:
        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            root = ET.fromstring(archive.read('xl/sharedStrings.xml'))
            shared = [''.join(n.text or '' for n in item.findall('.//m:t', ns))
                      for item in root.findall('m:si', ns)]
        sheet = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
    cells = {}
    for cell in sheet.findall('.//m:sheetData/m:row/m:c', ns):
        value = cell.find('m:v', ns)
        if cell.get('t') == 's':
            text = shared[int(value.text)] if value is not None else ''
        elif cell.get('t') == 'inlineStr':
            text = ''.join(n.text or '' for n in cell.findall('.//m:t', ns))
        else:
            text = (value.text or '') if value is not None else ''
        cells[cell.attrib['r']] = text

    def column_name(number):
        result = ''
        while number:
            number, digit = divmod(number-1, 26)
            result = chr(65+digit) + result
        return result

    def equivalent(a, b):
        if a == b or (a in ('', 'NA', 'NaN', 'nan') and b in ('', 'NA', 'NaN', 'nan')):
            return True
        try:
            return math.isclose(float(a), float(b), rel_tol=0, abs_tol=1e-10)
        except (ValueError, TypeError):
            return False

    for row_number, record in enumerate(actual, 3):
        for j, field in enumerate(fields, 1):
            coordinate = column_name(j) + str(row_number)
            if not equivalent(cells.get(coordinate, ''), record[field]):
                raise ValueError('Supplementary Table 1 XLSX/TSV mismatch: ' + coordinate)
    last = len(actual)+2
    for j, column in enumerate(spec['columns'], 1):
        header_row = 2 if column['kind'] == 'amino_acid' else 1
        if cells.get(column_name(j)+str(header_row)) != column['label']:
            raise ValueError('Supplementary Table 1 display header differs: ' + column['field'])
    for summary in spec['summaries']:
        j = fields.index(summary['field'])+1
        for offset, key in ((2, 'carriers'), (3, 'callable')):
            coordinate = column_name(j)+str(last+offset)
            if not equivalent(cells.get(coordinate, ''), str(summary[key])):
                raise ValueError('Supplementary Table 1 summary differs: ' + coordinate)
    filter_node = sheet.find('m:autoFilter', ns)
    if filter_node is None or filter_node.get('ref') != 'A2:' + column_name(len(fields)) + str(last):
        raise ValueError('Supplementary Table 1 filter must include only headers and sample rows')
    all_text = '\n'.join(cells.values())
    if not all(note in all_text for note in spec['footnotes']):
        raise ValueError('Supplementary Table 1 lacks the full callability footnotes')
    return {'genomes': len(actual), 'amino_acid_positions': spec['amino_acid_site_count'],
            'uncallable_cells': sum(s['uncallable'] for s in spec['summaries']),
            'protein_types_assigned': spec['protein_types_assigned'],
            'workbook_cells_checked': len(actual)*len(fields)}


def validate(run):
    from PIL import Image
    import fitz
    run = Path(run).resolve()
    checks = []

    def require(condition, message):
        if not condition:
            raise ValueError(message)
        checks.append(message)

    index = json.loads((run / 'FIGURE_INDEX.json').read_text())
    require(set(index) == FIGURES, 'Five main and four supplementary figures')
    for name in sorted(FIGURES):
        for ext in ('pdf', 'svg', 'png'):
            path = run / 'Figures' / (name + '.' + ext)
            require(path.is_file() and path.stat().st_size > 0, 'Nonempty ' + path.name)
        with fitz.open(run / 'Figures' / (name + '.pdf')) as doc:
            require(len(doc) == 1, name + ': one PDF page')
            text=doc[0].get_text()
            require(not substituted_glyphs(text),
                    name + ': no substituted glyph sequences in PDF text')
        with Image.open(run / 'Figures' / (name + '.png')) as image:
            image.verify()

    panel_root = run / 'Figure_Source_Data'
    _, panels = table(panel_root / 'PANEL_DATA_INDEX.tsv')
    require({(r['figure'], r['panel']) for r in panels} == PANELS,
            'All 27 manuscript panels have indexed source data')
    for row in panels:
        path = (panel_root / row['file']).resolve()
        require(path.is_relative_to(panel_root) and path.is_file() and path.stat().st_size > 0,
                'Panel data: ' + row['file'])

    _, carrier_years = table(panel_root / 'Figure1/F/carrier_counts_by_year.tsv')
    validate_carrier_year_rows(carrier_years)
    require(True, 'Carrier-year source data have no dated zero-count placeholders')

    from structure_stage import verify_render
    render = json.loads((run/'Structure_Render/RENDER_STATUS.json').read_text())
    verify_render(run/'Structure_Render', run/'work/figures/data/mrcA_model.pdb', render['image_width_px'])
    require(render.get('original_image_assets_used') is False, 'Structure panels regenerated from the reference PDB')

    data = run / 'work/figures/data'
    frozen = json.loads((run / 'work/INPUT_MANIFEST.json').read_text())
    require(all((data / name).is_file() and digest(data / name) == expected
                for name, expected in frozen.items()),
            'Prepared input checksums unchanged through export')

    tables = run / 'Supplementary_Tables'
    full_tables = run / 'Reproducibility_Data/Analysis_Tables'
    _, metadata = table(full_tables / 'sample_metadata.tsv')
    require(len(metadata) == 1861 and len({r['sample'] for r in metadata}) == 1861,
            'Metadata: 1861 unique genomes')
    for model in ('strict', 'relaxed'):
        ids = {r['sample'] for r in metadata
               if r['included_' + model + '_dating'].upper() == 'TRUE'}
        require(len(ids) == 533 and all(r['lineage'] == 'SS14' for r in metadata
                                       if r['sample'] in ids), model + ': 533 SS14 genomes')
        if model == 'strict':
            strict_ids = ids
        else:
            require(ids == strict_ids, 'Both clocks use the same genomes')
    assigned = [r for r in metadata if r['mrcA_protein_type'] not in ('', 'NA')]
    require(len(assigned) == 1328 and len({r['mrcA_protein_type'] for r in assigned}) == 16,
            '1328 assigned genomes and 16 protein types')
    metadata_validation = validate_metadata_exports(run)
    require(True, 'Formal metadata TSV and XLSX match all audited codon states and footnotes')

    _, models = table(full_tables / 'ER_ARD_model_sensitivity.tsv')
    require(len(models) == 23 and len({r['locus'] for r in models}) == 23
            and all(r['analysis_status'] == 'complete' for r in models),
            '23 completed ER-ARD comparisons, including R107G')
    require(any(r['locus'] == 'tp47:p.Arg107Gly' for r in models), 'R107G result present')
    differences = [float(r['delta_AIC_ARD_ER']) for r in models]
    require(all(math.isfinite(x) for x in differences)
            and sum(abs(x) < 2 for x in differences) == 22,
            '22 comparisons have absolute AIC differences below two')

    _, selection = table(full_tables / 'all_codon_selection_results.tsv')
    require(Counter(r['gene'] for r in selection) == GENES
            and len({(r['gene'], r['site']) for r in selection}) == 2594,
            '2594 unique codons across four genes')
    _, observed = table(tables / 'TableS2_PBP_variants_and_ESM2_scores.tsv')
    require(len(observed) == 40 and len({r['locus'] for r in observed}) == 40
            and {r['gene'] for r in observed} <= set(GENES)
            and all(int(r['carriers']) > 0 and math.isfinite(float(r['ESM2_LLR']))
                    for r in observed), 'ESM supplementary table: 40 observed replacements only')
    _, residue_rows = table(tables / 'TableS1_sample_metadata.tsv')
    for variant in observed:
        position = str(int(variant['position']))
        field = variant['gene'] + ':' + variant['reference_aa'] + position
        states = [r[field] for r in residue_rows]
        require(sum(x == variant['alternate_aa'] for x in states) == int(variant['carriers'])
                and sum(x != 'uncallable' for x in states) == int(variant['callable_genomes']),
                'Supplementary tables 1 and 2 agree: ' + variant['locus'])

    _, structural = table(full_tables / 'reference_structural_metrics.tsv')
    require(len(structural) == 3 and {r['gene'] for r in structural} == {'mrcA'}
            and {int(r['position']) for r in structural} == {506,625,708}
            and all(int(r['SASA_sphere_points']) == 2880 for r in structural),
            'Structural table contains only three MrcA residues, each at 2880-point resolution')

    _, table_index = table(tables / 'TABLE_INDEX.tsv')
    require(len(table_index) == 3 and {r['table'] for r in table_index} == set('123'),
            'Exactly three manuscript supplementary tables')
    for row in table_index:
        fields, rows = table(tables / row['part'])
        require(len(rows) == int(row['rows']) and len(fields) == int(row['columns']),
                'Supplementary table index: ' + row['part'])
    public_files = [p for folder in ('Figures', 'Figure_Source_Data', 'Supplementary_Tables')
                    for p in (run / folder).rglob('*') if p.is_file()]
    require(not any('pagel' in p.name.lower() or 'cooccurrence' in p.name.lower()
                    for p in public_files), 'Publication exports exclude Pagel and co-occurrence files')
    files = sorted(public_files + [run / 'FIGURE_INDEX.json'])
    manifest = {str(p.relative_to(run)): {'bytes': p.stat().st_size, 'sha256': digest(p)}
                for p in files}
    (run / 'PUBLICATION_OUTPUT_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    result = {'state': 'PASS', 'figures': 9, 'panels': len(PANELS),
              'metadata_export': metadata_validation,
              'indexed_panel_files': len(panels), 'supplementary_table_files': len(table_index),
              'output_checksums': len(manifest), 'checks_passed': len(checks), 'checks': checks}
    (run / 'OUTPUT_VALIDATION.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    try:
        report = validate(args.run)
        print(json.dumps({k: v for k, v in report.items() if k != 'checks'}, indent=2))
    except Exception as error:
        (args.run / 'OUTPUT_VALIDATION.json').write_text(json.dumps(
            {'state': 'FAIL', 'error': str(error)}, indent=2) + '\n')
        raise
