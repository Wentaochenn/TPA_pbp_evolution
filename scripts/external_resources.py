"""Validate and export the identities of external analysis databases."""
import json
from pathlib import PurePosixPath
import re

REQUIRED = ('strainseeker', 'snpeff_ss14', 'snpeff_nichols')
NAMES = {
    'strainseeker': 'StrainSeeker database',
    'snpeff_ss14': 'snpEff TPA_SS14 database',
    'snpeff_nichols': 'snpEff TPA_Nichols database',
}


def informative(value):
    text = str(value or '').strip()
    return bool(text) and text.lower() not in ('unknown', 'unavailable', 'not recorded', 'pending', 'na', 'n/a') and not text.lower().startswith('must match')


def relative(value):
    path = PurePosixPath(str(value))
    if not str(value) or path.is_absolute() or '..' in path.parts or '\\' in str(value):
        raise ValueError('Resource file names must be relative: ' + str(value))
    return str(path)


def public_provenance(value):
    if not isinstance(value, dict) or value.get('schema_version') != 1 or not isinstance(value.get('resources'), list):
        raise ValueError('Expected resource_provenance.json schema_version 1 with resources')
    output, seen, missing = [], set(), []
    for record in value['resources']:
        resource_id = record.get('resource_id')
        if not isinstance(resource_id, str) or resource_id in seen:
            raise ValueError('Resource identifiers must be nonempty and unique')
        seen.add(resource_id)
        clean = {key: record[key] for key in ('resource_id', 'name', 'software_version', 'source_description',
                 'release', 'build_date', 'recipe', 'download_url') if key in record}
        files = []
        for row in record.get('files', []):
            if not re.fullmatch(r'[a-fA-F0-9]{64}', str(row.get('sha256', ''))):
                raise ValueError('Resource file SHA256 is invalid: ' + resource_id)
            size = row.get('size', row.get('size_bytes'))
            if not isinstance(size, int) or isinstance(size, bool) or size < 0:
                raise ValueError('Resource file size is invalid: ' + resource_id)
            files.append({'path': relative(row.get('path', '')), 'size': size, 'sha256': row['sha256'].lower()})
        if len({r['path'] for r in files}) != len(files):
            raise ValueError('Duplicate resource file name: ' + resource_id)
        clean['files'] = files
        if record.get('archived_path'):
            clean['archived_path'] = relative(record['archived_path'])
        # A local path or a manually set boolean alone does not establish reproducible identity.
        identity = record.get('identity_complete') is True and bool(files)
        source = informative(clean.get('source_description'))
        access = informative(clean.get('archived_path')) or informative(clean.get('recipe')) or (
            informative(clean.get('download_url')) and informative(clean.get('release')))
        retrieval = record.get('retrieval_complete') is True and source and access
        clean['identity_complete'] = identity
        clean['retrieval_complete'] = retrieval
        if resource_id in REQUIRED and not (identity and retrieval):
            missing.append(resource_id)
        output.append(clean)
    missing.extend(resource_id for resource_id in REQUIRED if resource_id not in seen)
    return {'schema_version': 1, 'resources': output,
            'upstream_database_versions_complete': not missing,
            'unresolved_resources': sorted(missing)}


def load_provenance(path):
    return public_provenance(json.loads(path.read_text()))


def resource_rows(provenance, model_id):
    by_id = {r['resource_id']: r for r in (provenance or {}).get('resources', [])}
    rows = [dict(resource='Human host reference', version='hg19', purpose='BBMap host-read removal', included='no'),
            dict(resource='Rabbit host reference', version='oryCun2', purpose='BBMap host-read removal', included='no')]
    for resource_id in REQUIRED:
        record = by_id.get(resource_id, {})
        version = record.get('release') or record.get('build_date') or ('exact file identities in resource_provenance.json' if record.get('identity_complete') else 'not recorded')
        rows.append(dict(resource=NAMES[resource_id], version=version,
                         purpose='Read classification' if resource_id == 'strainseeker' else 'Variant annotation',
                         included='yes' if record.get('archived_path') else 'no'))
    rows.append(dict(resource=model_id, version='file checksums in esm2_model.json', purpose='ESM-2 inference', included='no'))
    return rows
