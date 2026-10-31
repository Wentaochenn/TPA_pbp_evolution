#!/usr/bin/env python3
"""Extract existing CA-height HPD annotations for time-tree plotting.

Reads completed annotated trees; does not estimate a tree or posterior. Clades
are matched by complete descendant sample sets, independent of ladderisation.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass, field
import hashlib
import io
import json
import math
from pathlib import Path
import re


@dataclass(eq=False)
class Node:
    label: str | None = None
    length: float | None = None
    annotations: dict[str, str] = field(default_factory=dict)
    children: list['Node'] = field(default_factory=list)
    parent: 'Node | None' = None


def _split_top(text, separator=','):
    items, start, braces, brackets, quote = [], 0, 0, 0, None
    i = 0
    while i < len(text):
        char = text[i]
        if quote:
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    i += 2
                    continue
                quote = None
        elif char in "'\"": quote = char
        elif char == '{': braces += 1
        elif char == '}': braces -= 1
        elif char == '[': brackets += 1
        elif char == ']': brackets -= 1
        elif char == separator and braces == 0 and brackets == 0:
            items.append(text[start:i].strip()); start = i + 1
        if braces < 0 or brackets < 0: raise ValueError('Unbalanced annotation delimiters')
        i += 1
    if quote or braces or brackets: raise ValueError('Unclosed annotation/translation field')
    items.append(text[start:].strip())
    return items


def _unquote(text):
    text = text.strip()
    if len(text) > 1 and text[0] in "'\"" and text[-1] == text[0]:
        return text[1:-1].replace(text[0] * 2, text[0])
    return text


def _statement_end(text, start):
    quote, brackets, i = None, 0, start
    while i < len(text):
        char = text[i]
        if quote:
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    i += 2; continue
                quote = None
        elif char in "'\"": quote = char
        elif char == '[': brackets += 1
        elif char == ']': brackets -= 1
        elif char == ';' and brackets == 0: return i
        i += 1
    raise ValueError('NEXUS statement has no terminating semicolon')


def _tokens(text):
    i = 0
    while i < len(text):
        char = text[i]
        if char.isspace(): i += 1
        elif char in '(),:;':
            yield char, char; i += 1
        elif char == '[':
            start, depth, quote = i + 1, 1, None; i += 1
            while i < len(text) and depth:
                char = text[i]
                if quote:
                    if char == quote:
                        if i + 1 < len(text) and text[i + 1] == quote:
                            i += 2; continue
                        quote = None
                elif char in "'\"": quote = char
                elif char == '[': depth += 1
                elif char == ']': depth -= 1
                i += 1
            if depth: raise ValueError('Unclosed Newick comment')
            yield 'comment', text[start:i - 1]
        elif char in "'\"":
            start, quote = i, char; i += 1
            while i < len(text):
                if text[i] == quote:
                    if i + 1 < len(text) and text[i + 1] == quote:
                        i += 2; continue
                    i += 1; break
                i += 1
            else: raise ValueError('Unclosed quoted Newick label')
            yield 'label', _unquote(text[start:i])
        else:
            start = i
            while i < len(text) and not text[i].isspace() and text[i] not in '(),:;[]': i += 1
            yield 'label', text[start:i]


def read_tree(path):
    text = Path(path).read_text(encoding='utf-8-sig')
    matches = list(re.finditer(r'(?im)^\s*tree\s+[^=\r\n]+=', text))
    if len(matches) != 1: raise ValueError('Expected one annotated summary tree: ' + str(path))
    translations = {}
    translation = re.search(r'(?im)^\s*translate\b', text)
    if translation:
        end = _statement_end(text, translation.end())
        for item in _split_top(text[translation.end():end]):
            pair = re.fullmatch(r'\s*(\S+)\s+(.+?)\s*', item, re.S)
            if not pair: raise ValueError('Malformed NEXUS translation: ' + item[:80])
            key, name = pair.group(1), _unquote(pair.group(2))
            if key in translations: raise ValueError('Duplicate translation key: ' + key)
            translations[key] = name
    start = matches[0].end(); tree_text = text[start:_statement_end(text, start) + 1]
    root, current, stack, need_length, ended = None, None, [], False, False
    for kind, value in _tokens(tree_text):
        if ended: raise ValueError('Text after completed Newick tree')
        if need_length:
            if kind != 'label' or current is None: raise ValueError('Branch length missing')
            current.length = float(value)
            if not math.isfinite(current.length): raise ValueError('Nonfinite branch length')
            need_length = False; continue
        if kind == '(':
            node = Node()
            if stack:
                node.parent = stack[-1]; stack[-1].children.append(node)
            elif root is not None: raise ValueError('Multiple Newick roots')
            else: root = node
            stack.append(node); current = None
        elif kind == ')':
            if not stack: raise ValueError('Unmatched closing parenthesis')
            current = stack.pop()
            if len(current.children) != 2: raise ValueError('Expected bifurcating summary tree')
        elif kind == ',':
            if not stack or current is None: raise ValueError('Unexpected comma')
            current = None
        elif kind == 'label':
            if current is not None:
                if not current.children or current.label is not None: raise ValueError('Unexpected label: ' + value)
                current.label = value
            else:
                if not stack: raise ValueError('Expected rooted bifurcating tree')
                current = Node(label=translations.get(value, value), parent=stack[-1])
                stack[-1].children.append(current)
        elif kind == 'comment':
            if current is not None and value.startswith('&'):
                for item in _split_top(value[1:]):
                    if '=' not in item: continue
                    key, raw = item.split('=', 1); key = key.strip()
                    if key in current.annotations and current.annotations[key] != raw.strip():
                        raise ValueError('Conflicting node annotation: ' + key)
                    current.annotations[key] = raw.strip()
        elif kind == ':':
            if current is None or current.length is not None: raise ValueError('Unexpected branch length')
            need_length = True
        elif kind == ';':
            if stack or root is None or current is not root: raise ValueError('Incomplete Newick tree')
            ended = True
    if not ended or need_length: raise ValueError('Incomplete annotated summary tree')
    return root


def _float(node, name):
    if name not in node.annotations: raise ValueError('Node annotation absent: ' + name)
    value = float(node.annotations[name])
    if not math.isfinite(value): raise ValueError('Nonfinite annotation: ' + name)
    return value


def _interval(node, name):
    raw = node.annotations.get(name, '')
    if not (raw.startswith('{') and raw.endswith('}')): raise ValueError('HPD interval absent: ' + name)
    values = tuple(float(x) for x in _split_top(raw[1:-1]))
    if len(values) != 2 or not all(math.isfinite(x) for x in values) or values[0] > values[1]:
        raise ValueError('Invalid node HPD interval: ' + name)
    return values


def extract_clock(data_dir, clock):
    data = Path(data_dir); tree = read_tree(data / (clock + '.tree'))
    with (data / (clock + '_calibration.tsv')).open(newline='') as handle:
        calibrations = list(csv.DictReader(handle, delimiter='\t'))
    if len(calibrations) != 1: raise ValueError('Expected one calibration record')
    cal = calibrations[0]; root_date = float(cal['mcc_root_date_median']); reference_year = float(cal['latest_tip_date'])
    if not all(math.isfinite(x) for x in (root_date, reference_year)): raise ValueError('Nonfinite calibration')
    with (data / 'metadata.tsv').open(newline='') as handle:
        metadata_rows = list(csv.DictReader(handle, delimiter='\t'))
    metadata = {row['sample']: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows): raise ValueError('Duplicate metadata sample identifiers')
    nodes, depths, stack = [], {tree: 0.0}, [tree]
    while stack:
        node = stack.pop(); nodes.append(node)
        for child in node.children:
            if child.length is None or child.length < 0: raise ValueError('Missing or negative dated branch length')
            depths[child] = depths[node] + child.length
        stack.extend(reversed(node.children))
    tips = [n for n in nodes if not n.children]; internal = [n for n in nodes if n.children]
    names = [n.label for n in tips]
    if len(tips) != 533 or len(internal) != 532 or len(set(names)) != 533:
        raise ValueError('Expected 533 distinct tips and 532 internal nodes')
    if int(cal['n_mcc_tips']) != 533 or int(cal['n_mcc_tips_with_dates']) != 533:
        raise ValueError('Calibration counts differ from tree')
    tip_dates = []
    for tip in tips:
        if not tip.label or ';' in tip.label or not tip.label.isascii() or tip.label not in metadata:
            raise ValueError('Unusable/unmatched dated sample identifier: ' + str(tip.label))
        year = float(metadata[tip.label]['year'])
        if not math.isfinite(year) or not math.isclose(root_date + depths[tip], year, abs_tol=1e-6, rel_tol=0):
            raise ValueError('Terminal tree date differs from metadata: ' + tip.label)
        tip_dates.append(year)
    if not math.isclose(max(tip_dates), reference_year, abs_tol=1e-8, rel_tol=0):
        raise ValueError('Latest-tip reference differs from dated samples')
    descendants, rows = {}, []
    for node in reversed(nodes):
        if not node.children: descendants[node] = (node.label,); continue
        labels = tuple(sorted(label for child in node.children for label in descendants[child]))
        descendants[node] = labels
        mean = _float(node, 'CAheight_mean'); low, high = _interval(node, 'CAheight_95%_HPD')
        point = root_date + depths[node]
        if not math.isclose(point, reference_year - mean, abs_tol=1e-6, rel_tol=0):
            raise ValueError('CA mean height and drawn node date disagree')
        posterior = _float(node, 'posterior')
        if not 0 <= posterior <= 1: raise ValueError('Node posterior outside [0,1]')
        key_bytes = json.dumps(labels, ensure_ascii=True, separators=(',', ':')).encode('ascii')
        rows.append(dict(clock=clock, node_clade_key=hashlib.sha256(key_bytes).hexdigest(),
            descendant_tip_set=';'.join(labels), n_descendants=len(labels), posterior=posterior,
            point_date=point, hpd_lower_date=reference_year-high, hpd_upper_date=reference_year-low,
            ca_height_mean=mean, ca_height_hpd_lower=low, ca_height_hpd_upper=high,
            reference_year=reference_year, height_summary_method='common_ancestor_mean',
            interval_annotation='CAheight_95%_HPD'))
    if len({r['node_clade_key'] for r in rows}) != len(rows): raise ValueError('Duplicate clade identities')
    rows.sort(key=lambda r: (r['n_descendants'], r['node_clade_key']))
    root_row = next(r for r in rows if r['n_descendants'] == 533)
    return rows, dict(clock=clock, tips=533, internal_nodes=532, terminal_dates_match_metadata=True,
        negative_branches=0, root_date=root_date, root_hpd_lower_date=root_row['hpd_lower_date'],
        root_hpd_upper_date=root_row['hpd_upper_date'], reference_year=reference_year)


def prepare_clock_hpd(data_dir, output_dir=None):
    """Write node tables into NEW working data, or an explicit read-only-review output."""
    data = Path(data_dir).resolve(strict=True)
    output = Path(output_dir).resolve() if output_dir is not None else data
    extracted = {clock: extract_clock(data, clock) for clock in ('strict', 'relaxed')}
    output.mkdir(parents=True, exist_ok=True)
    reports = []
    for clock, (rows, report) in extracted.items():
        path = output / (clock + '_node_hpd.tsv')
        text = io.StringIO(newline='')
        writer = csv.DictWriter(text, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
        expected = text.getvalue().encode('utf-8')
        if path.exists():
            if path.is_symlink() or path.read_bytes() != expected:
                raise ValueError('Existing node-HPD table differs from its annotated tree: ' + str(path))
        else:
            path.write_bytes(expected)
        report['file'] = str(path); report['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        reports.append(report)
    return dict(state='PASS', scientific_models_recomputed=False,
        node_identity='sorted_complete_descendant_sample_set', clocks=reports)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    a = p.parse_args(); print(json.dumps(prepare_clock_hpd(a.data, a.out), ensure_ascii=False, indent=2))


if __name__ == '__main__': main()
