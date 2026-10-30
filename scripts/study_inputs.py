"""Prepare public study inputs, retrieve paired reads, or package a cluster release."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import tarfile
import urllib.parse
import urllib.request

from integrity import file_manifest, sha, verify_data, verify_manifest
from external_resources import load_provenance, public_provenance, resource_rows

CODE = Path(__file__).resolve().parent.parent
UPSTREAM = CODE / 'scripts/upstream'
RUN_ID = re.compile(r'(?:SRR|ERR|DRR)[0-9]+\Z')


def rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def cohort(metadata):
    records = rows(metadata)
    required = {'sample', 'year', 'lineage', 'country'}
    if not records or not required.issubset(records[0]):
        raise ValueError('Use Supplementary Table 1 with sample, year, lineage and country columns')
    by_sample = {r['sample']: r for r in records}
    if len(by_sample) != len(records) or any(not s for s in by_sample):
        raise ValueError('Metadata sample identifiers must be nonempty and unique')
    excluded = set((UPSTREAM / 'remove_samples.txt').read_text().split())
    groups = {}
    for group in ('ss14', 'nichols'):
        selected = [s for s in (UPSTREAM / 'list' / ('tpa.' + group + '.sampleList.txt')).read_text().split() if s not in excluded]
        if len(selected) != len(set(selected)):
            raise ValueError('Duplicate read identifier in the study list')
        for s in selected:
            if not RUN_ID.fullmatch(s) or s in groups:
                raise ValueError('Invalid or repeated run identifier: ' + s)
            groups[s] = group
    assemblies = rows(UPSTREAM / 'config/complete_genome_manifest.tsv')
    expected = set(groups) | {r['sample_id'] for r in assemblies}
    if set(by_sample) != expected:
        raise ValueError('Supplementary Table 1 does not match the deposited cohort: missing=%s; extra=%s' %
                         (sorted(expected-set(by_sample)), sorted(set(by_sample)-expected)))
    if len(records) != 1861 or len(groups) != 1856 or len(assemblies) != 5:
        raise ValueError('Expected 1861 genomes, 1856 read runs and five assemblies')
    for r in records:
        if r['lineage'] not in ('Nichols', 'SS14'):
            raise ValueError('Unexpected lineage: ' + r['sample'])
        value = r['year'].strip()
        if value and value.lower() not in ('na', 'nan'):
            year = float(value)
            if not math.isfinite(year) or not 1800 <= year <= 2100:
                raise ValueError('Invalid sampling year: ' + r['sample'])
    return records, groups, assemblies


def export_metadata(metadata, out, dates=None):
    from openpyxl import Workbook
    records, groups, assemblies = cohort(metadata)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(metadata, out / 'sample_metadata.tsv')
    book = Workbook()
    sheet = book.active
    sheet.title = 'Metadata'
    sheet.append(['SRA_ID (Key)', 'Isolation Year', 'Clade (Lineage)', 'Country', 'Source', 'BioProject', 'BioSample'])
    date_values = {}
    for r in records:
        value = r['year'].strip()
        year = float(value) if value and value.lower() not in ('na', 'nan') else None
        if year is not None:
            date_values[r['sample']] = year
        sheet.append([r['sample'], year, r['lineage'], r['country'], r.get('source', ''), r.get('bioproject', ''), r.get('biosample', '')])
    if dates:
        original = {}
        with Path(dates).open() as handle:
            for line in handle:
                values = line.strip().split()
                if not values:
                    continue
                if len(values) != 2:
                    raise ValueError('Sampling dates must have two columns without a header')
                sample, value = values
                year = float(value)
                if not math.isfinite(year) or sample in original:
                    raise ValueError('Invalid or duplicated sampling date: ' + sample)
                original[sample] = year
        original = {s: v for s, v in original.items() if s in {r['sample'] for r in records}}
        if set(original) != set(date_values):
            raise ValueError('Sampling-date coverage differs from Supplementary Table 1')
        for sample, year in original.items():
            if int(year) != int(date_values[sample]):
                raise ValueError('Sampling year differs from Supplementary Table 1: ' + sample)
        date_values = original
    book.save(out / 'TPA_master_metadata.xlsx')
    (out / 'sample_dates.tsv').write_text(''.join('%s\t%s\n' % (s, format(v, '.17g')) for s, v in sorted(date_values.items())))
    lookup = {r['sample']: r for r in records}
    manifest = [{'run_accession': s, 'reference_group': group, 'biosample': lookup[s].get('biosample', ''),
                 'bioproject': lookup[s].get('bioproject', '')} for s, group in sorted(groups.items())]
    write_tsv(out / 'read_accessions.tsv', manifest, ['run_accession', 'reference_group', 'biosample', 'bioproject'])
    shutil.copyfile(UPSTREAM / 'config/complete_genome_manifest.tsv', out / 'complete_genome_manifest.tsv')
    for group in ('ss14', 'nichols'):
        target = out / 'list' / ('tpa.' + group + '.sampleList.txt')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(''.join(s+'\n' for s, g in groups.items() if g == group))
    return {'genomes': len(records), 'read_runs': len(groups), 'assemblies': len(assemblies), 'dated_genomes': len(date_values)}


def prepare(data, out):
    data, out = data.resolve(), out.resolve()
    verify_data(data)
    if out.exists() or out.is_relative_to(data) or data.is_relative_to(out) or out.is_relative_to(CODE):
        raise ValueError('Use a new source directory outside the deposited data and code')
    inputs = data / 'Upstream_Inputs'
    if not inputs.is_dir():
        raise ValueError('This data release lacks Upstream_Inputs; create the updated cluster package first')
    shutil.copytree(inputs, out)
    summary = export_metadata(out / 'sample_metadata.tsv', out / 'metadata_adapter', out / 'sample_dates.tsv')
    for name in ('TPA_master_metadata.xlsx', 'sample_dates.tsv'):
        shutil.copyfile(out / 'metadata_adapter' / name, out / name)
    print(json.dumps(dict(state='prepared', **summary, source=str(out)), indent=2))


def md5(path):
    h = hashlib.md5()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def download_reads(source, selected=None):
    source = source.resolve()
    records = rows(source / 'read_accessions.tsv')
    if selected:
        wanted = set(selected)
        unknown = wanted - {r['run_accession'] for r in records}
        if unknown:
            raise ValueError('Requested run outside the study: ' + ','.join(sorted(unknown)))
        records = [r for r in records if r['run_accession'] in wanted]
    for r in records:
        run, group = r['run_accession'], r['reference_group']
        if not RUN_ID.fullmatch(run) or group not in ('ss14', 'nichols'):
            raise ValueError('Invalid read manifest entry')
        query = urllib.parse.urlencode(dict(accession=run, result='read_run', fields='run_accession,library_layout,fastq_ftp,fastq_md5', format='json'))
        with urllib.request.urlopen('https://www.ebi.ac.uk/ena/portal/api/filereport?' + query, timeout=60) as response:
            report = json.load(response)
        if len(report) != 1 or report[0]['run_accession'] != run or report[0]['library_layout'] != 'PAIRED':
            raise ValueError('ENA does not report one paired-end run for ' + run)
        links, checksums = report[0]['fastq_ftp'].split(';'), report[0]['fastq_md5'].split(';')
        pairs = {}
        for link, digest in zip(links, checksums):
            match = re.search(r'_([12])\.fastq\.gz$', link)
            if not match or not re.fullmatch('[0-9a-fA-F]{32}', digest):
                raise ValueError('Unexpected ENA file naming or checksum for ' + run)
            mate = match.group(1)
            if mate in pairs:
                raise ValueError('Multiple ENA files for the same mate: ' + run)
            url = urllib.parse.urlsplit(link if '://' in link else 'https://' + link)
            if url.scheme not in ('http', 'https', 'ftp') or url.hostname != 'ftp.sra.ebi.ac.uk':
                raise ValueError('Unexpected ENA read-file host')
            pairs[mate] = ('https://' + url.netloc + url.path, digest.lower())
        if len(links) != len(checksums) or set(pairs) != {'1', '2'}:
            raise ValueError('Expected exactly two paired ENA files for ' + run)
        for mate, (url, digest) in sorted(pairs.items()):
            dest = source / 'raw_fastq' / group / (run + '_' + mate + '.fastq.gz')
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                if md5(dest) != digest:
                    raise ValueError('Existing FASTQ differs from ENA checksum: ' + str(dest))
                continue
            partial = dest.with_suffix(dest.suffix + '.partial')
            with urllib.request.urlopen(url, timeout=120) as response, partial.open('wb') as handle:
                shutil.copyfileobj(response, handle)
            if md5(partial) != digest:
                raise ValueError('Downloaded FASTQ differs from ENA checksum: ' + str(partial))
            partial.replace(dest)
        print('Downloaded and verified ' + run, flush=True)


def copy_required(source, target):
    if not source.is_file():
        raise FileNotFoundError('Required cluster input missing: ' + str(source))
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def copy_study_input(source, target):
    if target.is_file() and source.is_file() and sha(source) != sha(target):
        raise ValueError('Current project input differs from the deposited file: ' + str(source))
    copy_required(source, target)


def clean_model_context(value):
    fields = ['model_files', 'torch', 'python', 'transformers', 'pandas', 'cuda_runtime', 'cudnn', 'device', 'dtype', 'scoring_source_sha256', 'gpu']
    result = {key: value[key] for key in fields if key in value}
    result['model_id'] = 'facebook/esm2_t36_3B_UR50D'
    for record in result.get('model_files', []):
        path = Path(record['path'])
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Model file names must be relative')
    return result


def download_model(specification, out, revision=None):
    from huggingface_hub import snapshot_download
    spec = clean_model_context(json.loads(specification.read_text()))
    records = spec.get('model_files', [])
    if not records:
        raise ValueError('Model specification has no archived file checksums')
    if out.exists():
        raise ValueError('Use a new local model directory')
    snapshot_download(repo_id=spec['model_id'], revision=revision,
                      allow_patterns=[r['path'] for r in records], local_dir=str(out))
    for record in records:
        path = out / record['path']
        if not path.is_file() or sha(path) != record['sha256']:
            raise ValueError('Model download does not match archived weights/tokenizer: ' + record['path'] +
                             '. Select a repository revision matching the archived checksums.')
    print('Model files match the archived ESM-2 checksums:', out)


def archive(path, entries):
    with path.open('wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w', format=tarfile.GNU_FORMAT) as tar:
            for source, name in sorted(entries, key=lambda pair: pair[1]):
                info = tar.gettarinfo(str(source), arcname=name)
                info.uid = info.gid = 0
                info.uname = info.gname = ''
                info.mtime = 0
                with source.open('rb') as handle:
                    tar.addfile(info, handle)


def package(project, data, out, resources=None):
    project, data, out = project.resolve(), data.resolve(), out.resolve()
    code_entries = json.loads((CODE / 'PACKAGE_MANIFEST.json').read_text())['files']
    verify_manifest(CODE, code_entries)
    verified = verify_data(data)
    metadata = data / 'Supplementary_Tables/TableS1_sample_metadata.tsv'
    cohort(metadata)
    if out.exists() or any(out.is_relative_to(r) or r.is_relative_to(out) for r in (data, CODE)):
        raise ValueError('Use a new package output directory outside code and deposited data')
    out.mkdir(parents=True)
    stage = out / '_package_inputs'
    inputs = stage / 'Upstream_Inputs'
    original = json.loads((data / 'DATA_MANIFEST.json').read_text())
    for relative in original:
        if relative.startswith('Upstream_Inputs/'):
            copy_required(data / relative, stage / relative)
    summary = export_metadata(metadata, inputs, project / 'sample_dates.tsv')
    for name in ('TPA_SS14.fasta', 'TPA_Nichols.fasta'):
        copy_study_input(project / 'refs' / name, inputs / 'refs' / name)
    for name in ('*.gff', '*.gff3', '*.gb', '*.gbk', '*.gbff'):
        for source in sorted((project / 'refs').glob(name)):
            if source.stat().st_size > 50*1024*1024:
                raise ValueError('Unexpectedly large TPA annotation: ' + str(source))
            copy_study_input(source, inputs / 'refs' / source.name)
    for r in rows(UPSTREAM / 'config/complete_genome_manifest.tsv'):
        copy_study_input(project / r['fasta'], inputs / r['fasta'])
    for relative in ('results/01_masked_alignment/premask_regions.bed', 'nichols_premask_regions.bed'):
        source = project / relative
        if source.is_file():
            copy_study_input(source, inputs / 'masks' / source.name)
    protein = data / 'Upstream_Inputs/reference_proteins.fasta'
    if not protein.is_file():
        protein = data / 'Model_Inputs/ESM2_reference_proteins.fasta'
    copy_required(protein, inputs / 'reference_proteins.fasta')
    model_spec = data / 'Upstream_Inputs/esm2_model.json'
    if model_spec.is_file():
        context = clean_model_context(json.loads(model_spec.read_text()))
    else:
        contexts = list((data / 'Model_Inputs').rglob('MODEL_CONTEXT.json'))
        if len(contexts) != 1:
            raise ValueError('Expected one archived ESM-2 MODEL_CONTEXT.json')
        context = clean_model_context(json.loads(contexts[0].read_text()))
    if not context.get('model_files'):
        raise ValueError('Archived ESM-2 file checksums are missing')
    write_json(inputs / 'esm2_model.json', context)
    provenance_path = resources.resolve() if resources else inputs / 'resource_provenance.json'
    provenance = load_provenance(provenance_path) if provenance_path.is_file() else public_provenance({'schema_version': 1, 'resources': []})
    if resources and not provenance_path.is_file():
        raise FileNotFoundError('Resource provenance file is missing: ' + str(provenance_path))
    write_json(inputs / 'resource_provenance.json', provenance)
    for record in provenance['resources']:
        if not record.get('archived_path'):
            continue
        source_dir = provenance_path.parent / record['archived_path']
        target_dir = inputs / record['archived_path']
        for row in record['files']:
            source = source_dir / row['path']
            target = target_dir / row['path']
            if source.resolve() != target.resolve():
                copy_required(source, target)
            size = row.get('size', row.get('size_bytes'))
            if not target.is_file() or target.stat().st_size != size or sha(target) != row['sha256']:
                raise ValueError('Archived resource differs from its declared identity: ' + str(target))
    external = resource_rows(provenance, context['model_id'])
    write_tsv(inputs / 'external_resources.tsv', external, ['resource', 'version', 'purpose', 'included'])
    adapter = project / 'refs/TruSeq3-PE.fa'
    if adapter.is_file():
        copy_study_input(adapter, inputs / 'refs/TruSeq3-PE.fa')
    instructions = ('# Study input files\n\n'
        'Sample metadata are Supplementary Table 1. References and five assembly FASTAs are the files used by the study. '
        'Sampling dates retain the original numerical values.\n\n'
        'Prepare inputs with `bash master.sh inputs prepare --data /path/to/Data_for_Zenodo --out /path/to/source_data`. '
        'Download paired reads with `bash master.sh inputs download-reads --source /path/to/source_data`; '
        'downloads use ENA file reports and MD5 verification.\n\n'
        'External databases and model weights are not redistributed. Set the database paths documented in the code README. '
        'ESM-2 model-file checksums are in `esm2_model.json`. '
        'The supplied files support reconstruction of study inputs; a complete raw-read rerun still requires the original '
        'StrainSeeker and snpEff database versions.\n')
    (inputs / 'README.md').write_text(instructions)
    input_files = file_manifest(inputs)
    input_files.pop('INPUT_MANIFEST.json', None)
    write_json(inputs / 'INPUT_MANIFEST.json', input_files)
    extra = file_manifest(stage)
    overlap = set(original) & set(extra)
    if any(not name.startswith('Upstream_Inputs/') for name in overlap):
        raise ValueError('Unexpected overlap outside Upstream_Inputs')
    entries = dict(original, **extra)
    write_json(stage / 'DATA_MANIFEST.json', entries)
    code_files = [(CODE / name, name) for name in sorted(code_entries)]
    code_files.append((CODE / 'PACKAGE_MANIFEST.json', 'PACKAGE_MANIFEST.json'))
    archive(out / 'TPA_pbp_evolution_code.tar.gz', code_files)
    data_files = [((stage if name in extra else data) / name, 'Data_for_Zenodo/' + name) for name in entries]
    data_files.append((stage / 'DATA_MANIFEST.json', 'Data_for_Zenodo/DATA_MANIFEST.json'))
    archive(out / 'TPA_pbp_evolution_data.tar.gz', data_files)
    archives = [out / 'TPA_pbp_evolution_code.tar.gz', out / 'TPA_pbp_evolution_data.tar.gz']
    (out / 'SHA256SUMS').write_text(''.join(sha(path) + '  ' + path.name + '\n' for path in archives))
    status = dict(state='packaged', **summary, deposited_files=verified['files_checked'], added_input_files=len(extra),
                  upstream_database_versions_complete=provenance['upstream_database_versions_complete'],
                  unresolved_resources=provenance['unresolved_resources'], original_data_modified=False)
    write_json(out / 'PACKAGE_STATUS.json', status)
    with tarfile.open(out / 'input_check.tar.gz', 'w:gz', format=tarfile.GNU_FORMAT) as tar:
        for name in ('PACKAGE_STATUS.json', 'SHA256SUMS'):
            tar.add(out / name, arcname=name)
        for path in sorted(inputs.rglob('*')):
            if path.is_file():
                tar.add(path, arcname=str(path.relative_to(stage)))
    print(json.dumps(status, indent=2))
    print('CODE:', archives[0])
    print('ZENODO_DATA:', archives[1])
    print('INPUT_CHECK:', out / 'input_check.tar.gz')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    q = sub.add_parser('prepare')
    q.add_argument('--data', type=Path, required=True)
    q.add_argument('--out', type=Path, required=True)
    q = sub.add_parser('download-reads')
    q.add_argument('--source', type=Path, required=True)
    q.add_argument('--runs', nargs='+')
    q = sub.add_parser('package')
    q.add_argument('--project', type=Path, required=True)
    q.add_argument('--data', type=Path, required=True)
    q.add_argument('--out', type=Path, required=True)
    q.add_argument('--resources', type=Path, help='Resource provenance JSON; defaults to the deposited input metadata')
    q = sub.add_parser('download-model')
    q.add_argument('--spec', type=Path, required=True)
    q.add_argument('--out', type=Path, required=True)
    q.add_argument('--revision', help='Hugging Face repository revision; all files must match archived checksums')
    args = p.parse_args()
    if args.action == 'prepare':
        prepare(args.data, args.out)
    elif args.action == 'download-reads':
        download_reads(args.source, args.runs)
    elif args.action == 'download-model':
        download_model(args.spec, args.out, args.revision)
    else:
        package(args.project, args.data, args.out, args.resources)


if __name__ == '__main__':
    main()
