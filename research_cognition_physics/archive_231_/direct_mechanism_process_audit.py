"""Round 392: infer wiring from a specified intervention process, with error bounds.

The laboratory slots and their port factorization are inputs. A recovered edge
is an operational mechanism, not a derivation of spatial contact or dimension.
Only NumPy; --write-results creates a new result and never overwrites evidence.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'direct_mechanism_process_audit_results.json'
TOL = 2e-10
I2 = np.eye(2, dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
# Axes: A.out.a, A.out.b, B.in, B.out.c, C.in. Trivial ports omitted.
DIMS = (2, 2, 2, 2, 2)
INPUTS = {'B': 2, 'C': 4}
OUTPUTS = {'a': (0, 'A'), 'b': (1, 'A'), 'c': (3, 'B')}


def trace_norm(a):
    return float(np.sum(np.linalg.svd(a, compute_uv=False)))


def marginal(a, dims, keep):
    """Partial trace with output axes in precisely the requested order."""
    keep = tuple(keep)
    rest = tuple(i for i in range(len(dims)) if i not in keep)
    order = keep + rest
    n = len(dims)
    dk = int(np.prod([dims[i] for i in keep]))
    dr = int(np.prod([dims[i] for i in rest]))
    t = a.reshape(tuple(dims)*2).transpose(order + tuple(n+i for i in order))
    return np.trace(t.reshape(dk, dr, dk, dr), axis1=1, axis2=3)


def product_factors(factors, dims):
    """Assemble disjoint matrix factors, then restore canonical axis order."""
    result = np.ones((1, 1), dtype=complex)
    order = []
    for axes, factor in factors:
        order.extend(axes)
        result = np.kron(result, factor)
    if sorted(order) != list(range(len(dims))):
        raise ValueError('Factors must partition all axes.')
    n = len(dims)
    permutation = tuple(order.index(i) for i in range(n))
    expanded = tuple(dims[i] for i in order)
    return result.reshape(expanded*2).transpose(
        permutation + tuple(n+i for i in permutation)).reshape(result.shape)


def depolarize(a, dims, axis):
    rest = tuple(i for i in range(len(dims)) if i != axis)
    return product_factors([((axis,), np.eye(dims[axis])/dims[axis]),
                            (rest, marginal(a, dims, rest))], dims)


def choi(kraus):
    """Unnormalized Choi, input axis first, without the lab-map transpose."""
    din = kraus[0].shape[1]
    dout = kraus[0].shape[0]
    out = np.zeros((din*dout, din*dout), dtype=complex)
    for k in kraus:
        v = k.T.reshape(-1)
        out += np.outer(v, v.conj())
    return out


def bell(sign=1):
    v = np.array([1, 0, 0, sign], dtype=complex)/np.sqrt(2)
    return np.outer(v, v.conj())


def xor_choi():
    """Normalized channel from (b,c) to v: measure Z, write b XOR c."""
    ks = []
    for b in range(2):
        for c in range(2):
            k = np.zeros((2, 4), dtype=complex)
            k[b ^ c, 2*b+c] = 1
            ks.append(k)
    return choi(ks)/4


def process(kind, strength=1.):
    if kind == 'chain':
        link = strength*bell() + (1-strength)*np.eye(4)/4
        return product_factors([((0, 2), link), ((3, 4), bell()),
                                ((1,), I2/2)], DIMS)
    if kind == 'shortcut':
        return product_factors([((0, 2), bell()),
                                ((1, 3, 4), xor_choi())], DIMS)
    if kind == 'memory':
        return sum(product_factors([((0, 2), bell(s)), ((3, 4), bell(s)),
                                    ((1,), I2/2)], DIMS)/2 for s in (1, -1))
    if kind == 'white':
        return np.eye(32, dtype=complex)/32
    raise ValueError(kind)


def infer(omega, dims=DIMS, inputs=INPUTS, outputs=OUTPUTS, eta=0.):
    """Threshold candidates; exact Markov status requires separate checks.

    With eta>0 the graph is certified only under the theorem's promised ideal
    tensor-Markov model and gap. A small residual never certifies exact Markov.
    """
    threshold = 2*eta + TOL
    parents = {name: [] for name in inputs}
    edges, open_ports, diagnostics = [], [], {}
    for port, (axis, owner) in outputs.items():
        delta = omega-depolarize(omega, dims, axis)
        gamma = trace_norm(delta)
        residuals = {name: trace_norm(marginal(delta, dims,
                      tuple(i for i in range(len(dims)) if i != target)))
                     for name, target in inputs.items()}
        if gamma <= threshold:
            open_ports.append(port)
            recipient = None
        else:
            matches = [name for name, r in residuals.items() if r <= threshold]
            if len(matches) != 1:
                raise ValueError('No unique recipient; the promised gap/model may fail.')
            recipient = matches[0]
            if recipient == owner:
                raise ValueError('Self-loop is outside the DAG model.')
            parents[recipient].append(port)
            edges.append((owner, recipient))
        diagnostics[port] = dict(gamma=gamma, after_input_trace=residuals,
                                 recipient=recipient)
    edges = sorted(set(edges))
    if any(a == b for a, b in closure(edges)):
        raise ValueError('Cycle is outside the DAG model.')
    blocks = []
    tp_errors = {}
    for name, target in inputs.items():
        axes = tuple(outputs[p][0] for p in parents[name]) + (target,)
        factor = marginal(omega, dims, axes)
        blocks.append((axes, factor))
        pd = int(np.prod([dims[i] for i in axes[:-1]]))
        tp_errors[name] = trace_norm(marginal(factor, tuple(dims[i] for i in axes),
                                   tuple(range(len(axes)-1)))-np.eye(pd)/pd)
    for port in open_ports:
        axis = outputs[port][0]
        blocks.append(((axis,), marginal(omega, dims, (axis,))))
    reconstructed = product_factors(blocks, dims)
    residual = trace_norm(omega-reconstructed)
    open_errors = [trace_norm(marginal(omega, dims, (outputs[p][0],))-
                             np.eye(dims[outputs[p][0]])/dims[outputs[p][0]])
                   for p in open_ports]
    return dict(edges=[list(e) for e in edges], open_ports=open_ports,
                parents=parents, diagnostics=diagnostics,
                blocks=[list(a) for a, _ in blocks], factorization_residual=residual,
                channel_trace_errors=tp_errors,
                numeric_exact_model_pass=(residual < TOL and max(tp_errors.values(), default=0.) < TOL
                                          and max(open_errors, default=0.) < TOL))


def closure(edges):
    found = set(tuple(e) for e in edges)
    while True:
        new = found | {(a, d) for a, b in found for c, d in found if b == c}
        if new == found:
            return sorted(found)
        found = new


def random_unitary(n, rng):
    q, r = np.linalg.qr(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
    phases = np.diag(r)
    return q @ np.diag(phases/np.abs(phases))


def experiment_probability(omega, rho_a, unitary_b, effect_c):
    # Process normalization is product of laboratory output dimensions 4*2.
    labs = np.kron(np.kron(rho_a.T, choi([unitary_b]).T), effect_c)
    return float(np.trace(8*omega@labs).real)


def circuit_probability(kind, rho_a, unitary_b, effect_c):
    out = np.zeros((2, 2), dtype=complex)
    if kind in ('chain', 'memory'):
        rho = marginal(rho_a, (2, 2), (0,))
        signs = (0,) if kind == 'chain' else (0, 1)
        for z in signs:
            zz = I2 if z == 0 else Z
            u = zz@unitary_b@zz
            out += u@rho@u.conj().T/len(signs)
    elif kind == 'shortcut':
        u = np.kron(unitary_b, I2)
        after = u@rho_a@u.conj().T  # Axes c,b.
        for c in range(2):
            for b in range(2):
                out[b ^ c, b ^ c] += after[2*c+b, 2*c+b]
    else:
        raise ValueError(kind)
    return float(np.trace(out@effect_c).real)


def report():
    models = {k: infer(process(k)) for k in ('chain', 'shortcut', 'memory')}
    calibrated = {}
    for kind in ('chain', 'shortcut'):
        omega = process(kind)
        estimate = .99*omega+.01*process('white')
        eta = trace_norm(estimate-omega)
        found = infer(estimate, eta=eta)
        gap = min(d['gamma'] for d in models[kind]['diagnostics'].values() if d['recipient'])
        calibrated[kind] = dict(eta=eta, threshold=2*eta, gap=gap,
                                gap_exceeds_4eta=gap > 4*eta, inferred=found)
    memory = .99*process('memory')+.01*process('white')
    eta = trace_norm(memory-process('memory'))
    found = infer(memory, eta=eta)
    rng = np.random.default_rng(392)
    probability_error = 0.
    for _ in range(8):
        vec = random_unitary(4, rng)[:, 0]
        rho = np.outer(vec, vec.conj())
        ub = random_unitary(2, rng)
        evec = random_unitary(2, rng)[:, 0]
        effect = np.outer(evec, evec.conj())
        for kind in ('chain', 'shortcut', 'memory'):
            probability_error = max(probability_error, abs(
                experiment_probability(process(kind), rho, ub, effect)-
                circuit_probability(kind, rho, ub, effect)))
    return dict(round=392, scope='Given laboratory slots, port factors, complete process and calibration; infer minimal tensor-Markov mechanisms. No spatial contact, 3D, autonomous tomography or GR derivation.',
                assumptions=['fixed finite input/output slots and port factorization',
                             'valid normalized complete intervention process',
                             'ideal tensor-Markov DAG for the recovery guarantee',
                             'trace-norm calibration and active-port margin for robustness'],
                normalized_process_dimension=32, process_trace=8,
                models=models, equal_transitive_closure=closure(models['chain']['edges']) == closure(models['shortcut']['edges']),
                calibration=calibrated,
                memory_rejection=dict(eta=eta, residual=found['factorization_residual'],
                                      m=len(found['blocks']), exclusion_threshold=(len(found['blocks'])+1)*eta,
                                      rejected=found['factorization_residual'] > (len(found['blocks'])+1)*eta),
                weak_edge=dict(strength=.01, distance_to_absent=trace_norm(process('chain', .01)-process('chain', 0.)),
                               active_residual=infer(process('chain', .01))['diagnostics']['a']['gamma']),
                independent_probability_error=probability_error,
                source='https://arxiv.org/html/1704.00800')


class Checks(unittest.TestCase):
    def test_01_partial_trace_and_depolarization(self):
        np.testing.assert_allclose(marginal(bell(), (2, 2), (0,)), I2/2, atol=TOL)
        np.testing.assert_allclose(depolarize(bell(), (2, 2), 0), np.eye(4)/4, atol=TOL)
        for kind in ('chain', 'shortcut', 'memory', 'white'):
            omega = process(kind)
            self.assertAlmostEqual(np.trace(omega).real, 1.)
            self.assertGreater(np.linalg.eigvalsh(omega).min(), -TOL)
            np.testing.assert_allclose(marginal(omega, DIMS, (0, 1, 3)), np.eye(8)/8, atol=TOL)

    def test_02_process_probabilities_against_independent_circuit(self):
        self.assertLess(report()['independent_probability_error'], TOL)
        rng = np.random.default_rng(393)
        for kind in ('chain', 'shortcut', 'memory'):
            for _ in range(4):
                v = random_unitary(4, rng)[:, 0]
                self.assertAlmostEqual(experiment_probability(process(kind), np.outer(v, v.conj()),
                                                             random_unitary(2, rng), I2), 1.)

    def test_03_distinct_minimal_graphs_same_reachability(self):
        chain, shortcut = infer(process('chain')), infer(process('shortcut'))
        self.assertEqual(chain['edges'], [['A', 'B'], ['B', 'C']])
        self.assertEqual(shortcut['edges'], [['A', 'B'], ['A', 'C'], ['B', 'C']])
        self.assertTrue(chain['numeric_exact_model_pass'] and shortcut['numeric_exact_model_pass'])
        self.assertEqual(closure(chain['edges']), closure(shortcut['edges']))

    def test_04_axis_permutation_does_not_supply_wiring(self):
        permutation = (4, 1, 3, 0, 2)
        inv = {old: new for new, old in enumerate(permutation)}
        for kind in ('chain', 'shortcut'):
            p = process(kind).reshape(DIMS*2).transpose(permutation+tuple(5+i for i in permutation)).reshape(32, 32)
            found = infer(p, inputs={k: inv[v] for k, v in INPUTS.items()},
                          outputs={k: (inv[v], owner) for k, (v, owner) in OUTPUTS.items()})
            self.assertEqual(found['edges'], infer(process(kind))['edges'])
            self.assertTrue(found['numeric_exact_model_pass'])

    def test_05_xor_parents_invisible_in_pair_marginals(self):
        c = xor_choi()
        for axes in ((0, 2), (1, 2)):
            np.testing.assert_allclose(marginal(c, (2, 2, 2), axes), np.eye(4)/4, atol=TOL)
        # Fix the other input to |0>: either parent can select the output bit.
        self.assertEqual([int(np.argmax(np.diag(c)[2*j:2*j+2])) for j in range(4)], [0, 1, 1, 0])
        model = infer(process('shortcut'))
        self.assertEqual(model['parents']['C'], ['b', 'c'])
        self.assertAlmostEqual(model['diagnostics']['b']['gamma'], 1.)

    def test_06_shared_internal_noise_requires_markov_rejection(self):
        found = infer(process('memory'))
        self.assertEqual(found['edges'], [['A', 'B'], ['B', 'C']])
        self.assertFalse(found['numeric_exact_model_pass'])
        self.assertAlmostEqual(found['factorization_residual'], 1.)
        dephased = (bell(1)+bell(-1))/2
        independent = product_factors([((0, 2), dephased), ((3, 4), dephased), ((1,), I2/2)], DIMS)
        plus = (I2+X)/2
        rho = np.kron(plus, I2/2)
        self.assertAlmostEqual(experiment_probability(process('memory'), rho, I2, plus), 1.)
        self.assertAlmostEqual(experiment_probability(independent, rho, I2, plus), .5)

    def test_07_exact_unique_recipient_gap(self):
        for kind in ('chain', 'shortcut'):
            found = infer(process(kind))
            for d in found['diagnostics'].values():
                if d['recipient'] is None:
                    continue
                for recipient, residual in d['after_input_trace'].items():
                    self.assertAlmostEqual(residual, 0. if recipient == d['recipient'] else d['gamma'])

    def test_08_calibration_bounds_and_promised_graph_recovery(self):
        for kind in ('chain', 'shortcut'):
            ideal = process(kind)
            for alternative in ('chain', 'shortcut', 'memory', 'white'):
                estimate = .99*ideal+.01*process(alternative)
                eta = trace_norm(estimate-ideal)
                exact, noisy = infer(ideal), infer(estimate, eta=eta)
                self.assertEqual(noisy['edges'], exact['edges'])
                for port, data in exact['diagnostics'].items():
                    measured = noisy['diagnostics'][port]
                    self.assertLessEqual(abs(data['gamma']-measured['gamma']), 2*eta+TOL)
                    for name, value in data['after_input_trace'].items():
                        self.assertLessEqual(abs(value-measured['after_input_trace'][name]), 2*eta+TOL)
                self.assertLessEqual(noisy['factorization_residual'], (len(noisy['blocks'])+1)*eta+TOL)

    def test_09_non_markov_exclusion_survives_calibration(self):
        witness = report()['memory_rejection']
        self.assertTrue(witness['rejected'])
        self.assertGreater(witness['residual'], 10*witness['exclusion_threshold'])

    def test_10_no_uniform_absence_certificate_without_gap(self):
        absent = process('chain', 0.)
        for strength in (.1, .01, .001):
            active = process('chain', strength)
            distance = trace_norm(active-absent)
            self.assertAlmostEqual(distance, 1.5*strength)
            self.assertIn(['A', 'B'], infer(active)['edges'])
            self.assertNotIn(['A', 'B'], infer(absent)['edges'])
            self.assertAlmostEqual(infer(active)['diagnostics']['a']['gamma'], distance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['checks'] = dict(run=checked.testsRun, failures=len(checked.failures), errors=len(checked.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
