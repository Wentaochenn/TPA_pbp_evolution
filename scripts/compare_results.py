"""Compare regenerated figures and tables with the deposited reference outputs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read_table(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.reader(handle, delimiter='\t', escapechar='\\'))
    if not rows or len(set(rows[0])) != len(rows[0]):
        raise ValueError('Missing or duplicate table headers: ' + str(path))
    if any(len(row) != len(rows[0]) for row in rows):
        raise ValueError('Malformed table: ' + str(path))
    return rows


def same_cell(a, b, atol=1e-12, rtol=1e-10):
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    return (math.isnan(x) and math.isnan(y)) or math.isclose(x, y, abs_tol=atol, rel_tol=rtol)


def compare_table(reference, candidate):
    a, b = read_table(reference), read_table(candidate)
    if a[0] != b[0] or len(a) != len(b):
        return False
    return all(len(x) == len(y) and all(same_cell(i, j) for i, j in zip(x, y))
               for x, y in zip(a[1:], b[1:]))


def pdf_signature(path):
    import fitz
    pages = []
    with fitz.open(path) as doc:
        for page in doc:
            spans = []
            for block in page.get_text('dict')['blocks']:
                for line in block.get('lines', []):
                    for s in line.get('spans', []):
                        spans.append([s['text'], s['font'], round(s['size'], 4),
                                      [round(x, 3) for x in s['bbox']]])
            pages.append([round(page.rect.width, 3), round(page.rect.height, 3), spans])
    return pages


def image_signature(path):
    from PIL import Image
    with Image.open(path) as image:
        image = image.convert('RGB')
        return image.size, hashlib.sha256(image.tobytes()).hexdigest()


def compare(reference, candidate):
    reference, candidate = Path(reference).resolve(), Path(candidate).resolve()
    figures, tables = [], []
    from check_outputs import FIGURES
    for name in sorted(FIGURES):
        png_a, png_b = (x/'Figures'/(name+'.png') for x in [reference, candidate])
        pdf_a, pdf_b = (x/'Figures'/(name+'.pdf') for x in [reference, candidate])
        if not all(p.is_file() for p in [png_a, png_b, pdf_a, pdf_b]):
            raise ValueError('Missing reference or candidate figure: ' + name)
        figures.append({'figure': name, 'pixels_equal': image_signature(png_a)==image_signature(png_b),
                        'PDF_text_and_geometry_equal': pdf_signature(pdf_a)==pdf_signature(pdf_b)})
    # Exclude rasterisation internals; compare every indexed scientific TSV and all analysis tables.
    for folder in ['Figure_Source_Data', 'Supplementary_Tables', 'Reproducibility_Data/Analysis_Tables']:
        def names(root):
            return {p.relative_to(root/folder).as_posix() for p in (root/folder).rglob('*.tsv')
                    if 'rendered_layers' not in p.relative_to(root/folder).parts}
        a, b = names(reference), names(candidate)
        if not a or a != b:
            raise ValueError('Reference/candidate table sets differ in ' + folder +
                             '; reference-only=' + str(sorted(a-b)) + '; candidate-only=' + str(sorted(b-a)))
        for name in sorted(a):
            tables.append({'file': folder+'/'+name,
                           'values_equal': compare_table(reference/folder/name, candidate/folder/name)})
    numeric = all(r['values_equal'] for r in tables)
    graphics = all(r['pixels_equal'] and r['PDF_text_and_geometry_equal'] for r in figures)
    return {'state': 'PASS' if numeric and graphics else 'DIFFERENT',
            'table_values_equal': numeric, 'figure_rendering_equal': graphics,
            'figures_compared': len(figures), 'tables_compared': len(tables),
            'numeric_absolute_tolerance': 1e-12, 'numeric_relative_tolerance': 1e-10,
            'figures': figures, 'tables': tables}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--out', type=Path)
    a = p.parse_args()
    result = compare(a.reference, a.candidate)
    target = a.out or a.candidate/'REPRODUCTION_CHECK.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['figures','tables']}, indent=2))
    return 0 if result['state']=='PASS' else 1


if __name__=='__main__':
    raise SystemExit(main())
