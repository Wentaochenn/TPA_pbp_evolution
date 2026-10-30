#!/usr/bin/env python3
"""Binary Mk likelihood fitting on a supplied known-state tree with a flat root prior."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy.linalg import expm
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

try:
    from numba import njit
except ImportError:
    def njit(func):
        return func


@njit
def loglik(a, b, starts, ends, children, lengths, tip):
    n = len(starts)
    L = np.ones((n, 2))
    logscale = 0.
    tot = a + b
    for i in range(n):
        if starts[i] == ends[i]:
            if tip[i] >= 0:
                L[i, 1-tip[i]] = 0.
        else:
            v0 = 1.
            v1 = 1.
            for k in range(starts[i], ends[i]):
                j = children[k]
                z = -math.expm1(-tot*lengths[k]) if tot > 0 else 0.
                p01 = a/tot*z if tot > 0 else 0.
                p10 = b/tot*z if tot > 0 else 0.
                v0 *= (1.-p01)*L[j, 0] + p01*L[j, 1]
                v1 *= p10*L[j, 0] + (1.-p10)*L[j, 1]
            m = max(v0, v1)
            if m <= 0:
                return -math.inf
            L[i, 0] = v0/m
            L[i, 1] = v1/m
            logscale += math.log(m)
    return logscale + math.log(.5*(L[-1, 0]+L[-1, 1]))


class BinaryLikelihood:
    def __init__(self, tree_path, states_path):
        self.tree = Phylo.read(str(tree_path), 'newick')
        with Path(states_path).open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if not rows or len({r['sample'] for r in rows}) != len(rows):
            raise ValueError('Missing or duplicate sample states')
        self.states = {r['sample']: int(r['state']) for r in rows}
        tips = [n.name for n in self.tree.get_terminals()]
        if len(set(tips)) != len(tips) or set(tips) != set(self.states):
            raise ValueError('Tree and known-state sample sets differ')
        if set(self.states.values()) != {0, 1}:
            raise ValueError('Both binary states are required')
        nodes = list(self.tree.find_clades(order='postorder'))
        ids = {id(n): i for i, n in enumerate(nodes)}
        starts, ends, children, lengths, tip = [], [], [], [], []
        for n in nodes:
            starts.append(len(children))
            for child in n.clades:
                t = child.branch_length
                if t is None or not math.isfinite(t) or t < 0:
                    raise ValueError('Invalid branch length')
                children.append(ids[id(child)])
                lengths.append(t)
            ends.append(len(children))
            tip.append(self.states[n.name] if n.is_terminal() else -1)
        self.arrays = tuple(np.array(v, dtype=float if i==3 else np.int64)
                            for i, v in enumerate([starts, ends, children, lengths, tip]))

    def __call__(self, a, b):
        return loglik(float(a), float(b), *self.arrays)

    def objective(self, p):
        if max(p) > 35 or min(p) < -35:
            return 1e100
        return -self(*np.exp(p))

    def independent_loglik(self, a, b):
        Q = np.array([[-a, a], [b, -b]])
        conditional = {}
        for node in self.tree.find_clades(order='postorder'):
            if node.is_terminal():
                v = np.array([-np.inf, -np.inf])
                v[self.states[node.name]] = 0.
            else:
                v = np.zeros(2)
                for child in node.clades:
                    transition = expm(Q*child.branch_length)
                    if not np.isfinite(transition).all() or np.min(transition) < 0:
                        raise ValueError('Invalid matrix exponential')
                    with np.errstate(divide='ignore'):
                        v += logsumexp(np.log(transition)+conditional[id(child)][None, :], axis=1)
            conditional[id(node)] = v
        return float(logsumexp(conditional[id(self.tree.root)]+np.log([.5, .5])))


def fit_pair(likelihood):
    er = minimize_scalar(lambda z: -likelihood(math.exp(z), math.exp(z)),
                         bounds=(-20, 25), method='bounded', options={'xatol': 1e-11})
    if not er.success or not math.isfinite(er.fun):
        raise ValueError('ER optimization failed')
    attempts = []
    for d0 in [-7, -3, 0, 3, 7]:
        for d1 in [-7, -3, 0, 3, 7]:
            start = np.array([er.x+d0, er.x+d1])
            r = minimize(likelihood.objective, start, method='Nelder-Mead',
                         options={'maxiter': 2500, 'xatol': 1e-9, 'fatol': 1e-11})
            attempts.append({'start': start.tolist(), 'log_rates': r.x.tolist(),
                             'logL': float(-r.fun), 'converged': bool(r.success)})
    passing = [r for r in attempts if r['converged'] and math.isfinite(r['logL'])]
    if not passing:
        raise ValueError('ARD optimization failed')
    best = max(passing, key=lambda r: r['logL'])
    if sum(abs(r['logL']-best['logL']) < 1e-6 for r in passing) < 2:
        raise ValueError('ARD optimum was not recovered from multiple starts')
    finite = [r['logL'] for r in attempts if math.isfinite(r['logL'])]
    if max(finite)-best['logL'] > 1e-6 or best['logL']+1e-8 < -er.fun:
        raise ValueError('ARD optimum failed the finite-likelihood or nested-model check')
    boundary = []
    for direction in [0, 1]:
        def objective(z):
            q = [0., 0.]
            q[direction] = math.exp(z)
            return -likelihood(*q)
        r = minimize_scalar(objective, bounds=(-20, 25), method='bounded')
        value = float(-r.fun)
        boundary.append({'variable_rate': direction, 'logL': value if math.isfinite(value) else None})
        if math.isfinite(value) and value > best['logL']+1e-5:
            raise ValueError('A boundary fit exceeds the interior ARD solution')
    qer = math.exp(er.x)
    a, b = (float(x) for x in np.exp(best['log_rates']))
    rows = [{'model': 'ER', 'q01': qer, 'q10': qer, 'loglikelihood': float(-er.fun), 'AIC': float(2+2*er.fun)},
            {'model': 'ARD', 'q01': a, 'q10': b, 'loglikelihood': best['logL'], 'AIC': 4-2*best['logL']}]
    for row in rows:
        checked = likelihood.independent_loglik(row['q01'], row['q10'])
        row['independent_loglikelihood'] = checked
        if abs(checked-row['loglikelihood']) > 1e-8:
            raise ValueError('Independent matrix-exponential likelihood check failed')
    return rows, {'root_prior': [.5, .5], 'attempts': attempts, 'boundaries': boundary,
                  'delta_AIC_ARD_ER': rows[1]['AIC']-rows[0]['AIC']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tree', type=Path, required=True)
    p.add_argument('--states', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    rows, record = fit_pair(BinaryLikelihood(a.tree, a.states))
    with a.out.open('w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=rows[0], delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    a.out.with_suffix('.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')


if __name__=='__main__':
    main()
