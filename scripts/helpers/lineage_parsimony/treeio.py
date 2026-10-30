"""Strict, iterative reader for one rooted topology in Newick or NEXUS.

Root placement is supplied by the author. A syntactically rooted Newick is not
independent evidence that the root is biologically correct.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
import re


@dataclass
class Node:
    id: int
    parent: int | None = None
    children: list[int] = field(default_factory=list)
    label: str | None = None
    length: float | None = None


@dataclass
class Tree:
    nodes: dict[int, Node]
    root: int
    tips: dict[str, int]
    preorder: list[int]
    postorder: list[int]


def _tokens(text: str):
    """Tokenize without recursion; comments and escaped quoted labels supported."""
    i = 0
    while i < len(text):
        if text[i].isspace():
            i += 1
            continue
        if text[i] == '[':
            comment_start = i + 1
            depth = 1
            i += 1
            while i < len(text) and depth:
                if text[i] == '[':
                    depth += 1
                elif text[i] == ']':
                    depth -= 1
                i += 1
            if depth:
                raise ValueError('Unterminated Newick/NEXUS comment')
            if re.match(r'^\s*&U(?:\s*,|\s*$)', text[comment_start:i - 1], flags=re.I):
                raise ValueError('Explicit unrooted [&U] annotation is not accepted')
            continue
        if text[i] in "'\"":
            quote, value = text[i], []
            i += 1
            while i < len(text):
                if text[i] == quote:
                    if i + 1 < len(text) and text[i + 1] == quote:
                        value.append(quote)
                        i += 2
                        continue
                    i += 1
                    break
                value.append(text[i])
                i += 1
            else:
                raise ValueError('Unterminated quoted tree label')
            yield ('label', ''.join(value))
            continue
        if text[i] in '(),:;=':
            yield (text[i], text[i])
            i += 1
            continue
        start = i
        while i < len(text) and not text[i].isspace() and text[i] not in "[](),:;='\"":
            i += 1
        if i == start:
            raise ValueError('Invalid tree character at position ' + str(i))
        yield ('label', text[start:i])


def parse_newick(text: str) -> Tree:
    tokens = list(_tokens(text))
    nodes: dict[int, Node] = {}
    stack: list[int] = []
    root = None
    current = None
    expect_subtree = True
    expect_length = False
    ended = False
    just_closed = False
    label_seen = False

    def new_node(label=None):
        nonlocal root
        idx = len(nodes)
        parent = stack[-1] if stack else None
        nodes[idx] = Node(idx, parent=parent, label=label)
        if parent is not None:
            nodes[parent].children.append(idx)
        elif root is None:
            root = idx
        else:
            raise ValueError('More than one tree root')
        return idx

    for kind, value in tokens:
        if ended:
            raise ValueError('More than one tree or trailing tokens after semicolon')
        if expect_length:
            if kind != 'label':
                raise ValueError('Expected a numeric branch length')
            try:
                length = float(value)
            except ValueError as exc:
                raise ValueError('Invalid branch length: ' + value) from exc
            if not math.isfinite(length) or length < 0:
                raise ValueError('Non-finite or negative branch length: ' + value)
            nodes[current].length = length
            expect_length = False
            just_closed = False
            continue
        if kind == '(':
            if not expect_subtree:
                raise ValueError('Unexpected opening parenthesis')
            current = new_node()
            stack.append(current)
            expect_subtree = True
            just_closed = False
            label_seen = False
        elif kind == 'label':
            if expect_subtree:
                if not value:
                    raise ValueError('Empty tip label')
                current = new_node(value)
                expect_subtree = False
                label_seen = True
                just_closed = False
            elif just_closed and not label_seen:
                nodes[current].label = value
                label_seen = True
            else:
                raise ValueError('Unexpected extra tree label: ' + value)
        elif kind == ')':
            if expect_subtree or not stack:
                raise ValueError('Unbalanced or empty tree parentheses')
            current = stack.pop()
            expect_subtree = False
            just_closed = True
            label_seen = False
        elif kind == ',':
            if expect_subtree or not stack:
                raise ValueError('Unexpected comma in tree')
            expect_subtree = True
            just_closed = False
        elif kind == ':':
            if expect_subtree or current is None or nodes[current].length is not None:
                raise ValueError('Unexpected or duplicate branch length')
            expect_length = True
        elif kind == ';':
            if expect_subtree or stack:
                raise ValueError('Incomplete tree at semicolon')
            ended = True
        else:
            raise ValueError('Unexpected Newick token: ' + value)
    if root is None or stack or expect_subtree or expect_length or not ended:
        raise ValueError('Incomplete tree; exactly one semicolon-terminated tree required')
    preorder, todo = [], [root]
    while todo:
        nid = todo.pop()
        preorder.append(nid)
        todo.extend(reversed(nodes[nid].children))
    tips = {}
    for nid in preorder:
        node = nodes[nid]
        if not node.children:
            if not node.label or node.label in tips:
                raise ValueError('Missing or duplicate tree tip: ' + str(node.label))
            tips[node.label] = nid
    return Tree(nodes, root, tips, preorder, list(reversed(preorder)))


def read_tree(path: str | Path) -> Tree:
    text = Path(path).read_text(encoding='utf-8-sig')
    if not text.lstrip().upper().startswith('#NEXUS'):
        return parse_newick(text)
    # Split at tokenized semicolons, so semicolons inside quoted labels are safe.
    commands, command = [], []
    for token in _tokens(text):
        if token[0] == ';':
            commands.append(command)
            command = []
        else:
            command.append(token)
    if command:
        raise ValueError('Unterminated final NEXUS command')
    translate, tree_commands = {}, []
    in_trees = False
    for cmd in commands:
        # #NEXUS and BEGIN TREES may be in the same initial command.
        words = [v.upper() for k, v in cmd if k == 'label']
        if 'BEGIN' in words and 'TREES' in words:
            in_trees = True
            continue
        if not in_trees or not cmd:
            continue
        name = cmd[0][1].upper()
        if name in ('END', 'ENDBLOCK'):
            in_trees = False
        elif name == 'TRANSLATE':
            entries, entry = [], []
            for tok in cmd[1:] + [(',', ',')]:
                if tok[0] == ',':
                    entries.append(entry)
                    entry = []
                else:
                    entry.append(tok)
            for entry in entries:
                if len(entry) != 2 or any(k != 'label' for k, _ in entry):
                    raise ValueError('Invalid NEXUS TRANSLATE entry')
                key, label = entry[0][1], entry[1][1]
                if key in translate:
                    raise ValueError('Duplicate NEXUS TRANSLATE key: ' + key)
                translate[key] = label
        elif name in ('TREE', 'UTREE'):
            if name == 'UTREE':
                raise ValueError('Explicitly unrooted NEXUS tree is not accepted for event ordering')
            eqs = [i for i, tok in enumerate(cmd) if tok[0] == '=']
            if len(eqs) != 1:
                raise ValueError('Invalid NEXUS TREE declaration')
            tree_commands.append(cmd[eqs[0] + 1:])
    if len(tree_commands) != 1:
        raise ValueError('Expected exactly one NEXUS tree; found ' + str(len(tree_commands)))
    def quote_label(label):
        return "'" + label.replace("'", "''") + "'"
    rendered = ''.join(quote_label(v) if k == 'label' else v for k, v in tree_commands[0])
    tree = parse_newick(rendered + ';')
    if translate:
        tips = {}
        for old, nid in tree.tips.items():
            if old not in translate:
                raise ValueError('Tree tip absent from NEXUS TRANSLATE: ' + old)
            label = translate[old]
            if label in tips:
                raise ValueError('Duplicate translated tree tip: ' + label)
            tree.nodes[nid].label = label
            tips[label] = nid
        tree.tips = tips
    return tree
