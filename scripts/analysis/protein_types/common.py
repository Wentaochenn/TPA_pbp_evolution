"""Strict dependency-free table and metadata adapters. No sample ID guessing."""
import csv
import hashlib
import json
import math
import posixpath
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

FOCAL = ['mrcA:p.Ala506Thr', 'mrcA:p.Ala506Val', 'mrcA:p.Gly708Ser',
         'mrcA:p.Met625Val', 'pbp1:p.Pro564Leu', 'tp47:p.Ser394Arg']
MISSING = {'', 'NA', 'N/A', 'NAN', 'NONE', 'NULL', '?', '.', '-', 'UNKNOWN',
           'NOT AVAILABLE', 'NOT KNOWN', 'NOT PROVIDED', 'NOT APPLICABLE'}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def strict_json(path):
    def pairs(items):
        result = {}
        for k, v in items:
            if k in result:
                raise ValueError('Duplicate JSON key: ' + k)
            result[k] = v
        return result
    def bad_constant(value):
        raise ValueError('Non-finite JSON value: ' + value)
    def finite_float(value):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError('Non-finite JSON number: ' + value)
        return number
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),
                      object_pairs_hook=pairs, parse_constant=bad_constant,
                      parse_float=finite_float)


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                    allow_nan=False) + '\n', encoding='utf-8')


def table(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        out = csv.DictWriter(handle, fieldnames=fields, delimiter='\t',
                             extrasaction='raise', lineterminator='\n')
        out.writeheader()
        for row in rows:
            out.writerow({k: 'NA' if row.get(k) is None else row.get(k, '') for k in fields})


def _rows_to_dicts(rows):
    rows = iter(rows)
    try:
        header = [str(x).strip() for x in next(rows)]
    except StopIteration:
        raise ValueError('Empty table') from None
    if not header or any(not x for x in header) or len(set(header)) != len(header):
        raise ValueError('Blank or duplicate table header')
    result = []
    for n, row in enumerate(rows, 2):
        if not row or all(str(x).strip() == '' for x in row):
            continue
        if len(row) != len(header):
            raise ValueError(f'Table row {n}: expected {len(header)} fields, got {len(row)}')
        result.append(dict(zip(header, [str(x).strip() for x in row])))
    return header, result


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return _rows_to_dicts(csv.reader(handle, delimiter='\t'))


def read_matrix(path):
    header, rows = read_tsv(path)
    if len(header) < 2 or header[0] != 'sample':
        raise ValueError('Matrix must start with sample and contain locus columns')
    result = {}
    for row in rows:
        sid = row['sample']
        if not sid or sid in result:
            raise ValueError('Empty or duplicate matrix sample: ' + sid)
        values = {}
        for key in header[1:]:
            v = row[key]
            if v.upper() in MISSING:
                values[key] = None
            elif v in ('0', '1'):
                values[key] = int(v)
            else:
                raise ValueError(f'Invalid ternary matrix value {sid}/{key}: {v!r}')
        result[sid] = values
    if not result:
        raise ValueError('Matrix contains no samples')
    return result


def read_xlsx(path, sheet=None):
    """Read stored cell values; reject formulas and merged cells rather than guessing.

    Sheet selection is explicit or first visible sheet. This adapter reads years,
    not Excel serial dates. Numeric/date ambiguity is rejected by metadata checks.
    """
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
          'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read('xl/workbook.xml'))
        rel = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        links = {r.get('Id'): r.get('Target') for r in rel}
        sheets = list(wb.find('s:sheets', ns))
        candidates = [s for s in sheets if s.get('name') == sheet] if sheet else [
            s for s in sheets if s.get('state', 'visible') == 'visible'][:1]
        if len(candidates) != 1:
            raise ValueError('XLSX sheet missing/ambiguous: ' + str(sheet))
        selected = candidates[0]
        target = links[selected.get('{' + ns['r'] + '}id')]
        member = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)
        shared = []
        if 'xl/sharedStrings.xml' in z.namelist():
            for node in ET.fromstring(z.read('xl/sharedStrings.xml')):
                shared.append(''.join(t.text or '' for t in node.iter('{' + ns['s'] + '}t')))
        doc = ET.fromstring(z.read(member))
        if doc.find('s:mergeCells', ns) is not None:
            raise ValueError('Merged cells in metadata sheet; export a flat TSV')
        stored_rows = []
        width = None
        for row in doc.findall('s:sheetData/s:row', ns):
            vals = {}
            for cell in row:
                reference = cell.get('r', '')
                match = re.fullmatch(r'([A-Z]+)[0-9]+', reference)
                if not match:
                    raise ValueError('Invalid XLSX cell reference')
                idx = 0
                for letter in match[1]:
                    idx = idx * 26 + ord(letter) - 64
                if idx in vals:
                    raise ValueError('Duplicate XLSX cell')
                if cell.find('s:f', ns) is not None:
                    raise ValueError(f'Formula in metadata cell {reference}; export values-only TSV')
                if cell.get('t') == 'e':
                    raise ValueError('Excel error cell: ' + reference)
                if cell.get('t') == 'inlineStr':
                    value = ''.join(t.text or '' for t in cell.findall('.//s:t', ns))
                else:
                    value = cell.findtext('s:v', default='', namespaces=ns)
                    if cell.get('t') == 's':
                        index = int(value)
                        if index < 0 or index >= len(shared):
                            raise ValueError('Invalid XLSX shared-string index: ' + value)
                        value = shared[index]
                vals[idx] = value
            if not any(str(v).strip() for v in vals.values()):
                continue
            if width is None:
                width = max(k for k, v in vals.items() if str(v).strip())
            if any(k > width and str(v).strip() for k, v in vals.items()):
                raise ValueError('Metadata data extends beyond header')
            stored_rows.append([vals.get(i, '') for i in range(1, width + 1)])
        header, records = _rows_to_dicts(stored_rows)
        return header, records, selected.get('name')


ALIASES = {
    'sample': ['sample', 'SRA_ID (Key)', 'sample_id', 'taxon'],
    'year': ['year', 'Isolation Year', 'Sample_Year', 'sampling_year'],
    'lineage': ['lineage', 'Clade (Lineage)', 'TPA_Lineage'],
    'country': ['country', 'Country', 'Geo_Country'],
    'source': ['source', 'study', 'BioProject', 'bioproject', 'study_accession'],
    'clade': ['clade'],
}


def normalize_metadata(path, settings):
    if Path(path).suffix.lower() == '.xlsx':
        header, rows, sheet = read_xlsx(path, settings.get('sheet'))
    else:
        header, rows = read_tsv(path)
        sheet = None
    mapping = settings.get('columns', {})
    unknown = set(mapping) - set(ALIASES)
    if unknown:
        raise ValueError('Unknown normalized metadata column: ' + str(sorted(unknown)))
    chosen = {}
    for key, aliases in ALIASES.items():
        if key in mapping:
            name = mapping[key]
            if name is not None and name not in header:
                raise ValueError('Configured metadata column missing: ' + str(name))
        else:
            matches = [x for x in aliases if x in header]
            if len(matches) > 1:
                raise ValueError(f'Ambiguous metadata columns {key}: {matches}; configure columns')
            name = matches[0] if matches else None
        chosen[key] = name
    if chosen['sample'] is None:
        raise ValueError('No metadata sample column')
    data = {}
    lineage_map = settings.get('lineage_map', {})
    for row in rows:
        entry = {key: row[name] if name else '' for key, name in chosen.items()}
        sid = entry['sample']
        if not sid or sid in data:
            raise ValueError('Empty or duplicate metadata sample: ' + sid)
        entry['raw_lineage'] = entry['lineage']
        raw = entry['lineage']
        if raw in lineage_map:
            entry['lineage'] = lineage_map[raw]
        elif re.search('Nichols', raw, re.I) and not re.search('SS14', raw, re.I):
            entry['lineage'] = 'Nichols'
        elif re.search('SS14', raw, re.I) and not re.search('Nichols', raw, re.I):
            entry['lineage'] = 'SS14'
        elif raw.upper() in MISSING:
            entry['lineage'] = ''
        else:
            entry['lineage'] = raw
        year = entry['year']
        if year.upper() not in MISSING:
            try:
                numeric = float(year)
            except ValueError:
                raise ValueError(f'Unparsed collection year for {sid}: {year!r}; supply explicit numeric year') from None
            if not math.isfinite(numeric) or not 1800 <= numeric <= 2100:
                raise ValueError(f'Invalid year or Excel serial date for {sid}: {year!r}')
            entry['year'] = str(numeric)
        else:
            entry['year'] = ''
        for key in ('country', 'source', 'clade'):
            if entry[key].upper() in MISSING:
                entry[key] = ''
        data[sid] = entry
    return data, {'sheet': sheet, 'columns': chosen, 'n_metadata': len(data)}
