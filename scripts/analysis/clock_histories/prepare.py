#!/usr/bin/env python3
"""Read-only input audit and two-pass posterior-tree sampling; Python stdlib."""
import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

TREE_LINE = re.compile(rb'^\s*tree\s+\S+\s*=', re.I)
STATE = re.compile(rb'^\s*tree\s+STATE_(\d+)\b', re.I)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def table(path, rows, columns=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = columns or (list(rows[0]) if rows else [])
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def trace_summary(path):
    """Audit the observed trace, independently of the configured XML duration."""
    path = Path(path)
    stat = path.stat()
    header = None
    first = previous = None
    intervals = set()
    records = 0
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for line_number, line in enumerate(handle, 1):
            digest.update(line)
            tokens = line.decode('utf-8-sig').strip().split()
            if not tokens or tokens[0].startswith('#'):
                continue
            if header is None:
                header = [t.lower() for t in tokens]
                if header[0] not in ('sample', 'state') or 'posterior' not in header:
                    raise ValueError(f'Invalid trace header: {path}:{line_number}')
                continue
            if len(tokens) != len(header) or not tokens[0].isdigit():
                raise ValueError(f'Invalid/incomplete trace record: {path}:{line_number}')
            state = int(tokens[0])
            try:
                posterior = float(tokens[header.index('posterior')])
            except ValueError:
                raise ValueError(f'Non-numeric posterior: {path}:{line_number}') from None
            if not math.isfinite(posterior):
                raise ValueError(f'Non-finite posterior: {path}:{line_number}')
            if previous is not None:
                if state <= previous:
                    raise ValueError(f'Duplicate/unordered trace states: {path}:{line_number}')
                intervals.add(state - previous)
            else:
                first = state
            previous = state
            records += 1
    if records < 2:
        raise ValueError('Too few trace records: ' + str(path))
    if (path.stat().st_size, path.stat().st_mtime_ns) != (stat.st_size, stat.st_mtime_ns):
        raise RuntimeError('Trace changed while being audited: ' + str(path))
    return {'source': str(path), 'sha256': digest.hexdigest(),
            'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
            'records': records, 'first_state': first, 'last_state': previous,
            'logging_intervals': sorted(intervals)}


def audit_chain(trace, posterior, xml, out):
    """Gate on observed endpoints; XML length alone cannot establish run history."""
    trace = trace_summary(trace)
    lengths = [int(n.get('chainLength')) for n in ET.parse(xml).getroot().iter()
               if n.get('chainLength') is not None]
    configured = lengths[0] if len(lengths) == 1 else None
    matches = trace['last_state'] == posterior['last_state']
    result = {'trace': trace, 'posterior': posterior,
              'current_xml': {'source': str(xml), 'sha256': sha(xml),
                              'chain_lengths': lengths},
              'endpoint_agreement': matches,
              'observed_last_state': trace['last_state'] if matches else None,
              'xml_length_status': ('matches_observed_endpoint' if configured == trace['last_state']
                                    else 'run_history_review_required'),
              'interpretation': 'Agreement verifies the available chain endpoints, not historical completion. '
                                'Current XML chainLength may differ after resumes or edits; '
                                'review launch history before reporting the intended chain duration.'}
    Path(out).write_text(json.dumps(result, indent=2))
    if not matches:
        raise ValueError(f'Trace/posterior endpoint mismatch: trace={trace["last_state"]}, '
                         f'trees={posterior["last_state"]}; inspect {out}. '
                         'Do not rename states or change expected values to bypass this check.')
    return result


def sample_posterior(source, target, count=500, burnin=.10, thin=10):
    source, target = Path(source), Path(target)
    if count < 2 or thin < 1 or not 0 <= burnin < 1:
        raise ValueError('Invalid sampling settings')
    n = 0
    previous = None
    digest = hashlib.sha256()
    stat = source.stat()
    with source.open('rb') as handle:
        for line in handle:
            digest.update(line)
            if TREE_LINE.match(line):
                label = STATE.match(line)
                if label is None or not line.rstrip().endswith(b';'):
                    raise ValueError(f'Expected one complete STATE tree per line: {source}')
                state = int(label.group(1))
                if previous is not None and state <= previous:
                    raise ValueError(f'Duplicate/unordered posterior states: {source}')
                previous = state
                n += 1
    burn = math.floor(n * burnin)
    eligible = range(burn, n, thin)
    k = min(count, len(eligible))
    if k < 2:
        raise ValueError('Too few post-burn-in trees')
    positions = [round(i * (len(eligible) - 1) / (k - 1)) for i in range(k)]
    selected = {eligible[i] for i in positions}
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix('.tmp')
    mapping = []
    seen_tree = False
    tree_i = 0
    with source.open('rb') as handle, tmp.open('wb') as dest:
        for line in handle:
            if TREE_LINE.match(line):
                seen_tree = True
                if tree_i in selected:
                    # Keep topology/branch lengths/STATE labels. Annotations are
                    # unnecessary for stochastic mapping and can be very large.
                    dest.write(re.sub(rb'\[[^\]]*\]', b'', line))
                    mapping.append({'sample_index': len(mapping)+1,
                                    'original_tree_index': tree_i+1,
                                    'mcmc_state': int(STATE.match(line).group(1))})
                tree_i += 1
            elif not seen_tree:
                dest.write(line)
        dest.write(b'\nEnd;\n')
    if source.stat().st_size != stat.st_size or source.stat().st_mtime_ns != stat.st_mtime_ns:
        tmp.unlink(missing_ok=True)
        raise RuntimeError('Posterior changed while being sampled')
    if len(mapping) != k:
        raise ValueError('Posterior sampling count mismatch')
    tmp.replace(target)
    table(target.with_suffix('.indices.tsv'), mapping)
    metadata = {'source': str(source), 'sha256': digest.hexdigest(),
                'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                'total_trees': n, 'last_state': previous,
                'burnin_fraction': burnin, 'burnin_trees': burn, 'thin': thin,
                'eligible_trees': len(eligible), 'selected_trees': k,
                'first_selected_state': mapping[0]['mcmc_state'],
                'last_selected_state': mapping[-1]['mcmc_state'],
                'selection': 'evenly spaced over the entire post-burn-in thinned chain'}
    target.with_suffix('.sampling.json').write_text(json.dumps(metadata, indent=2))
    return metadata


def xml_summary(path):
    root = ET.parse(path).getroot()
    by_id={n.get('id'):n for n in root.iter() if n.get('id')}

    def weight_value(value):
        if value.startswith('@') and value[1:] in by_id:
            node=by_id[value[1:]]
            value=node.get('value') or node.text or ''
        tokens=re.split(r'[,\s]+',value.strip())
        try:return ' '.join(str(int(float(t))) for t in tokens)
        except ValueError:return ' '.join(value.split())
    sequences = {}
    for node in root.iter('sequence'):
        taxon = node.get('taxon')
        child = node.find('taxon')
        if not taxon and child is not None:
            taxon = child.get('idref') or child.get('id')
        seq = re.sub(r'\s+', '', node.get('value') or (node.text or '')).upper()
        if taxon and seq:
            if taxon in sequences:
                raise ValueError(f'Duplicate XML sequence taxon: {taxon}')
            sequences[taxon] = seq
    date_traits = []
    for node in root.iter():
        if str(node.get('traitname', '')).startswith('date'):
            values = {}
            for token in (node.get('value') or node.text or '').split(','):
                if '=' in token:
                    a, b = token.split('=', 1)
                    values[a.strip()] = float(b.strip())
            date_traits.append({'traitname': node.get('traitname'),
                                'units': node.get('units', 'year').rstrip('s'),
                                'values': values})
    constants = []
    specs = []
    model_nodes = []
    for node in root.iter():
        spec = node.get('spec', '')
        if spec:
            specs.append(spec)
        for key, value in node.attrib.items():
            if key == 'constantSiteWeights':
                constants.append(weight_value(value))
        if node.tag=='constantSiteWeights' or node.get('name')=='constantSiteWeights':
            constants.append(weight_value(node.get('value') or node.text or ''))
        if node.tag in ('siteModel', 'substModel') or 'Coalescent' in spec:
            model_nodes.append(ET.tostring(node, encoding='unicode'))
    return {'path': str(path), 'sha256': sha(path), 'n_sequences': len(sequences),
            'configured_chain_lengths': [int(n.get('chainLength')) for n in root.iter()
                                         if n.get('chainLength') is not None],
            'sequence_lengths': sorted({len(s) for s in sequences.values()}),
            'alignment_hash': hashlib.sha256(json.dumps(sequences, sort_keys=True).encode()).hexdigest(),
            'date_traits': date_traits, 'constant_site_weights': sorted(set(constants)),
            'clock_specs': sorted({s for s in specs if 'clock' in s.lower()}),
            'site_and_tree_model_xml': model_nodes,
            'required_packages': root.get('required')}
