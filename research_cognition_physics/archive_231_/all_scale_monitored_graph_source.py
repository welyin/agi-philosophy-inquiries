"""Round 517: all-finite-size source theorem for the original monitored flip model.

The infinite-time assertion is proved in the note, not inferred from sampled
spectra. Integer graph certificates and actual unknown-reference evolution are
independent diagnostics. No uniform-in-size mixing budget is claimed.
"""
import argparse
from collections import deque
from fractions import Fraction as Q
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'all_scale_monitored_graph_source_results.json'
OBS = {}


def words(counts, prefix=()):
    if not any(counts):
        yield prefix
    for v, count in enumerate(counts):
        if count:
            nxt = list(counts)
            nxt[v] -= 1
            yield from words(nxt, prefix + (v,))


def family(i):
    n = 2*i + 2
    trees = sorted((old.prufer_tree(n, w) for w in words([2]*i)),
                   key=lambda t: tuple(sorted(t)))
    lookup = {g:k for k,g in enumerate(trees)}
    transitions = []
    for g in trees:
        transitions.append({lookup[h]:w for h,w in old.tree_flips(g,n).items()})
    return n, trees, transitions


def connected(adj):
    seen = {0}
    todo = deque(seen)
    while todo:
        for w in adj[todo.popleft()]:
            if w not in seen:
                seen.add(w)
                todo.append(w)
    return len(seen) == len(adj)


def coefficient_dimension(g, h, n):
    """Exact rank by equalities/zero constraints; no floating rank tolerance."""
    adj = old.adjacency(g & h, n)
    forced = {v for e in g ^ h for v in e}
    seen = set()
    free = 0
    for start in range(n):
        if start in seen:
            continue
        component = {start}
        todo = [start]
        while todo:
            for v in adj[todo.pop()]:
                if v not in component:
                    component.add(v)
                    todo.append(v)
        seen |= component
        free += not bool(component & forced)
    return int(free)


def sparse_h(trees, transitions, n):
    """One-excitation restriction of sum SWAP plus original graph flips."""
    m = len(trees)
    rows = []
    cols = []
    values = []
    for g, tree in enumerate(trees):
        adj = old.adjacency(tree,n)
        for v in range(n):
            col = v*m + g
            rows.append(col);cols.append(col);values.append(len(tree)-len(adj[v]))
            for w in adj[v]:
                rows.append(w*m+g);cols.append(col);values.append(1)
            for gp, weight in transitions[g].items():
                rows.append(v*m+gp);cols.append(col);values.append(weight)
    return tuple(np.asarray(x,dtype=np.int64) for x in (rows,cols,values))


def apply_h(data, entries):
    rows, cols, values = entries
    out = np.zeros_like(data,dtype=complex)
    np.add.at(out, rows, values[:,None]*data[cols])
    return out


def evolve_columns(data, entries, time, order=18):
    term = data.astype(complex)
    out = term.copy()
    for k in range(1,order+1):
        term = (-1j*float(time)/k)*apply_h(term,entries)
        out += term
    return out


def unitary(h, time):
    ev,vec = np.linalg.eigh(h)
    return (vec*np.exp(-1j*time*ev))@vec.conj().T


def channel(x, u, n, m):
    y = u@x@u.conj().T
    out = np.zeros_like(y)
    for v in range(n):
        sl = slice(v*m,(v+1)*m)
        out[sl,sl] = y[sl,sl]
    return out


class Audit(unittest.TestCase):
    def test_01_original_all_scale_configuration_support(self):
        rows = []
        for i in (1,2,3,4):
            n,trees,f = family(i)
            m = len(trees)
            self.assertEqual(m,math.factorial(2*i)//2**i)
            self.assertEqual(len(set(trees)),m)
            self.assertTrue(connected(f))
            for g,t in enumerate(trees):
                adj = old.adjacency(t,n)
                self.assertTrue(connected(adj))
                self.assertEqual([len(a) for a in adj],[3]*i+[1]*(i+2))
                self.assertEqual(sum(f[g].values()),4*(i-1))
                for h,w in f[g].items():
                    self.assertEqual(w,1)
                    self.assertEqual(f[h][g],w)
            if i == 3:
                # Independent old bit-mask implementation, not only tree_flips.
                masks = [old.mask(t,n) for t in trees]
                for g,mask in enumerate(masks):
                    self.assertEqual(old.previous.flips(mask,n),
                                     {masks[h]:w for h,w in f[g].items()})
            rows.append(dict(internal=i,vertices=n,graphs=m,
                             directed_flips=sum(map(len,f)),
                             flip_degree=4*(i-1),configuration_connected=True))
        OBS['finite_configuration_certificates'] = rows

    def test_02_all_pair_intertwiner_constraints_without_common_edge(self):
        n,trees,f = family(3)
        common = set.intersection(*(set(t) for t in trees))
        self.assertEqual(common,set())
        for g in trees:
            self.assertEqual(coefficient_dimension(g,g,n),1)
        count = 0
        for g,h in itertools.combinations(trees,2):
            self.assertEqual(coefficient_dimension(g,h,n),0)
            count += 1
        OBS['exact_pair_constraints'] = dict(vertices=n,graphs=len(trees),
            unordered_distinct_pairs=count,common_to_all_edges=0,
            each_diagonal_pair_free_coefficients=1,
            each_offdiagonal_pair_free_coefficients=0)

    def test_03_actual_H_and_passive_reference_in_larger_family(self):
        n,trees,f = family(3)
        m = len(trees);d = n*m
        entries = sparse_h(trees,f,n)
        rows,cols,values = entries
        actual = {}
        for r,c,w in zip(rows,cols,values):
            if w:
                actual[(int(r),int(c))] = actual.get((int(r),int(c)),0)+int(w)
        independent = {}
        for g,t in enumerate(trees):
            for v in range(n):
                col = v*m+g
                bits = 1 << v
                for a,b in t:
                    switched = bits ^ ((1<<a)|(1<<b)) if v in (a,b) else bits
                    w = switched.bit_length()-1
                    key = (w*m+g,col)
                    independent[key] = independent.get(key,0)+1
                for gp,weight in f[g].items():
                    key = (v*m+gp,col)
                    independent[key] = independent.get(key,0)+weight
        self.assertEqual(actual,independent)
        self.assertTrue(all(actual.get((c,r))==w for (r,c),w in actual.items()))
        psi = np.zeros((d,2),complex)
        psi[0,0] = 1/np.sqrt(2)
        psi[m-1,1] = 1/np.sqrt(2)
        h_star = (n-1)+4*(3-1)
        tau = Q(1,8*h_star)
        evolved = evolve_columns(psi,entries,tau)
        blocks = evolved.reshape(n,m,2)
        probabilities = np.einsum('vgr,vgr->v',blocks.conj(),blocks).real
        reference = sum(b.T@b.conj() for b in blocks)
        self.assertLess(float(np.max(np.abs(reference-np.eye(2)/2))),2e-13)
        self.assertLess(abs(float(probabilities.sum())-1),2e-13)
        self.assertGreater(float(probabilities[1:].sum()),0)
        # Whole isometry retains coherent outcome labels; dephasing is only
        # the declared active-state marginal, not deletion of the record.
        tail = 2*Q(1,8)**19/math.factorial(19)
        self.assertLess(tail,Q(1,10**30))
        OBS['larger_actual_reference_diagnostic'] = dict(vertices=n,graphs=m,
            active_dimension=d,passive_reference_dimension=2,
            sparse_integer_H_equal_to_original_SWAPS=True,
            h_norm_upper=h_star,time=str(tau),taylor_order=18,
            analytic_operator_remainder_upper=str(tail),
            source_port_probabilities=[float(f'{p:.10g}') for p in probabilities],
            reference_unchanged_to_tolerance=True,
            numerical_mixing_time_or_diamond_certificate_claimed=False)

    def test_04_necessary_contract_boundaries(self):
        x = np.array([[0,1],[1,0]],complex)
        z = np.diag([1,-1]).astype(complex)
        eye = np.eye(2)
        tau = 0.125
        h = np.kron(x,eye)+np.kron(eye,x)
        u = unitary(h,tau)
        plus = np.array([1,1])/np.sqrt(2)
        minus = np.array([1,-1])/np.sqrt(2)
        witness = np.kron(eye,np.outer(plus,minus))
        residual = channel(witness,u,2,2)-np.exp(-2j*tau)*witness
        self.assertLess(float(np.max(np.abs(residual))),2e-13)
        # Connected graph-register flips cannot reveal duplicate graph labels.
        self.assertEqual(coefficient_dimension(frozenset({(0,1)}),
                                               frozenset({(0,1)}),2),1)
        no_j = unitary(np.kron(eye,x),tau)
        self.assertLess(float(np.max(np.abs(channel(np.kron(z,eye),no_j,2,2)
                                            -np.kron(z,eye)))),2e-13)
        disconnected_h = np.kron(eye,x)
        q = np.diag([1,1,0,0]).astype(complex)
        ud = unitary(disconnected_h,tau)
        self.assertLess(float(np.max(np.abs(channel(q,ud,4,1)-q))),2e-13)
        alias = unitary(x,math.pi)
        self.assertLess(float(np.max(np.abs(channel(z,alias,2,1)-z))),2e-13)
        OBS['boundary_counterexamples'] = dict(
            duplicate_graph_labels_have_nontrivial_peripheral_operator=True,
            no_data_coupling_preserves_port_information=True,
            disconnected_physical_graph_preserves_component_information=True,
           _aliased_sampling_preserves_nontrivial_diagonal_information=True)

    def test_05_finite_source_error_and_resource_quantifiers(self):
        rows = []
        for n in (6,8,10,100):
            eta = Q(1,100*n)
            p_min = Q(1,n)-eta
            conditional = eta/p_min
            self.assertEqual(conditional,Q(1,99))
            self.assertGreater(p_min,0)
            attempts = 100*n
            failure = (1-p_min)**attempts
            self.assertLess(failure,Q(1,1000))
            rows.append(dict(vertices=n,hypothetical_half_diamond_eta=str(eta),
                success_probability_lower=str(p_min),
                conditional_reference_error_upper=str(conditional),
                attempts=attempts,abort_bound_below_one_per_thousand=True,
                mixing_block_length_not_numerically_claimed=True))
        OBS['source_contract_examples'] = rows


def run():
    OBS.clear();stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28',round=517,scientific_base_through_round=516,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('branching_tree_distance_audit.py',
                         'simultaneous_edge_exchange_audit.py','research_note_492.md',
                         'research_note_493.md')},
        scope=dict(all_finite_distinct_connected_graph_families_theorem=True,
            unchanged_model_with_declared_port_instrument=True,
            arbitrary_initial_active_state_and_passive_reference=True,
            all_binary_tree_sizes_configuration_connected=True,
            no_data_or_graph_reset_between_monitoring_steps=True,
            uniform_in_size_mixing_budget_proved=False,
            full_record_environment_decoupled=False,
            coherent_uniform_graph_state_generated=False,
            autonomous_measurement_and_clock_generated=False,
            stable_macroscopic_position_generated=False,
            dimension_three_generated=False,full_GR_goal_completed=False),
        observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args();result = run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
