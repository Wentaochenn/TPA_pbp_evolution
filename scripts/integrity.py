"""Portable checksum validation for code, deposited data and prepared inputs."""
import hashlib
import json
from pathlib import Path


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def checked_path(root, relative):
    root = Path(root).resolve()
    name = Path(relative)
    if name.is_absolute() or not name.parts or '..' in name.parts:
        raise ValueError('Manifest path must be relative and contained: ' + str(relative))
    path = root / name
    if any(parent.is_symlink() for parent in [path, *path.parents] if parent != root and parent.is_relative_to(root)):
        raise ValueError('Manifest entries must not traverse symlinks: ' + str(relative))
    if not path.resolve().is_relative_to(root) or not path.is_file():
        raise ValueError('Manifest file is missing or outside the data directory: ' + str(relative))
    return path


def file_manifest(root, exclude=()):
    root = Path(root)
    result = {}
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if '__pycache__' in relative.parts or path.suffix == '.pyc' or str(relative) in exclude:
            continue
        if path.is_symlink():
            raise ValueError('Release files must be materialized, not symlinks: ' + str(relative))
        if path.is_file():
            result[str(relative)] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
    return result


def verify_manifest(root, entries):
    if not isinstance(entries, dict) or not entries:
        raise ValueError('The checksum manifest is empty or invalid')
    for relative, record in entries.items():
        path = checked_path(root, relative)
        expected = record if isinstance(record, str) else record['sha256']
        if isinstance(record, dict) and path.stat().st_size != record['bytes']:
            raise ValueError('File size differs from the frozen manifest: ' + relative)
        if sha(path) != expected:
            raise ValueError('File checksum differs from the frozen manifest: ' + relative)
    return {'state': 'PASS', 'files_checked': len(entries)}


def write_code_manifest(root):
    root = Path(root)
    record = {'files': file_manifest(root, exclude=['PACKAGE_MANIFEST.json'])}
    (root / 'PACKAGE_MANIFEST.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def verify_data(root):
    root = Path(root).resolve()
    entries = json.loads((root / 'DATA_MANIFEST.json').read_text())
    result = verify_manifest(root, entries)
    frozen = root / 'frozen_run'
    prepared = json.loads((frozen / 'work/INPUT_MANIFEST.json').read_text())
    verify_manifest(frozen / 'work/figures/data', prepared)
    sources = json.loads((frozen / 'SOURCE_MANIFEST.json').read_text())
    for key, record in sources.items():
        if Path(record['path']).is_absolute():
            raise ValueError('Frozen source paths must be relative: ' + key)
        path = (frozen / record['path']).resolve()
        if not path.is_relative_to(root):
            raise ValueError('Frozen source escapes the deposited data: ' + key)
        relative = str(path.relative_to(root))
        checked_path(root, relative)
        if relative not in entries or sha(path) != record['sha256']:
            raise ValueError('Frozen source differs from the deposited manifest: ' + key)
    result.update(prepared_files=len(prepared), source_files=len(sources))
    return result
