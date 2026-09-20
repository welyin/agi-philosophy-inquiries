"""Round 225: triple orientation consistency and cross-system CP certificates.

Finite parity enumeration is an exact combinatorial check. Matrix examples
check the analytic argument, not a numerical proof for all theories.
"""
import argparse
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np

from reversible_dynamics_bridge import dagger, hermitian_basis, sample_generator, sample_state, unitary
from tensor_process_bridge import apply_kraus, choi, extend_left, kraus_from_choi, partial_transpose


def edge_bits(n, bits):
    out = np.zeros((n, n), dtype=int)
    for (i, j), bit in zip(itertools.combinations(range(n), 2), bits):
        out[i, j] = out[j, i] = bit
    return out


def solve_gauge(tau, reference=0):
    if np.any(np.diag(tau)) or np.any(tau != tau.T):
        raise ValueError('Relative orientations must be symmetric with zero diagonal')
    s = tau[:, reference].copy()
    if np.any((s[:, None] ^ s[None, :]) != tau):
        raise ValueError('An odd triangle obstructs a consistent orientation')
    return s


def parity_certificate(n):
    total = 1 << (n*(n-1)//2)
    accepted = 0
    for bits in itertools.product((0, 1), repeat=n*(n-1)//2):
        tau = edge_bits(n, bits)
        triangles = all((tau[i, j] ^ tau[j, k] ^ tau[i, k]) == 0
                        for i, j, k in itertools.combinations(range(n), 3))
        try:
            s = solve_gauge(tau)
            solved = True
            assert np.array_equal(tau, s[:, None] ^ s[None, :])
        except ValueError:
            solved = False
        assert triangles == solved
        accepted += int(solved)
    return {'objects': n, 'all_edge_assignments': total, 'consistent': accepted,
            'expected_consistent': 2**(n-1), 'rejected': total-accepted}


def embed_pair(x, dims, pair):
    """Insert identity on omitted factors, preserving original factor order."""
    rest = [i for i in range(len(dims)) if i not in pair]
    order = list(pair)+rest
    ordered_dims = [dims[i] for i in order]
    large = np.kron(x, np.eye(int(np.prod([dims[i] for i in rest]))))
    positions = [order.index(i) for i in range(len(dims))]
    axes = positions+[i+len(dims) for i in positions]
    total = int(np.prod(dims))
    return large.reshape(ordered_dims*2).transpose(axes).reshape(total, total)


def rank_two_bell(a, b):
    v = np.zeros(a*b, dtype=complex)
    v[0], v[b+1] = 1/np.sqrt(2), 1/np.sqrt(2)
    return np.outer(v, v.conj())


def marginal_sign_certificate():
    dims = (2, 3, 2)
    rows = []
    for alpha in itertools.product((0, 1), repeat=3):
        for pair in itertools.combinations(range(3), 2):
            a, b = [dims[i] for i in pair]
            for e, f in itertools.product((0, 1), repeat=2):
                mismatch = (alpha[pair[0]] ^ e, alpha[pair[1]] ^ f)
                image = embed_pair(partial_transpose(rank_two_bell(a, b), a, b, *mismatch), dims, pair)
                low = float(np.linalg.eigvalsh(image).min())
                parity_ok = (alpha[pair[0]] ^ alpha[pair[1]]) == (e ^ f)
                assert abs(low-(0 if parity_ok else -.5)) < 1e-13
                rows.append({'triple_bits': list(alpha), 'pair': list(pair),
                             'pair_bits': [e, f], 'positive_on_witness': parity_ok,
                             'min_eigenvalue': low})
    return {'dimensions': list(dims), 'cases': len(rows),
            'consistent_cases': sum(r['positive_on_witness'] for r in rows),
            'inconsistent_cases': sum(not r['positive_on_witness'] for r in rows),
            'inconsistent_min_eigenvalue': -.5,
            'nonadjacent_pair_cases': sum(r['pair'] == [0, 2] for r in rows),
            'all_cases_match_analytic_parity': True}


class PairChart:
    def __init__(self, a, b, sa, sb, sab):
        self.a, self.b = a, b
        self.sa, self.sb, self.sab = sa, sb, sab
        self.u = unitary(sample_generator(a*b, seed=2250+a), .219)

    def encode(self, x):
        z = self.u @ partial_transpose(x, self.a, self.b, self.sa, self.sb) @ dagger(self.u)
        return z.T if self.sab else z

    def decode(self, x):
        x = x.T if self.sab else x
        return partial_transpose(dagger(self.u)@x@self.u, self.a, self.b, self.sa, self.sb)


def random_rectangular_channel(a, b, count=4):
    rng = np.random.default_rng(22500+10*a+b)
    raw = [rng.normal(size=(b, a))+1j*rng.normal(size=(b, a)) for _ in range(count)]
    vals, vecs = np.linalg.eigh(sum(dagger(k)@k for k in raw))
    normalizer = (vecs*(1/np.sqrt(vals)))@dagger(vecs)
    return [k@normalizer for k in raw]


def cross_channel_certificate(a, b, bits):
    sa, sb, sr, sar, sbr = bits
    r = a
    incoming, outgoing = PairChart(a, r, sa, sr, sar), PairChart(b, r, sb, sr, sbr)
    ks = random_rectangular_channel(a, b)
    standard = lambda x: apply_kraus(ks, x)
    def old(x):
        z = standard(x.T if sa else x)
        return z.T if sb else z
    def fixed(x):
        z = old(x.T if sa else x)
        return z.T if sb else z
    w = sample_state(a*r, seed=2251)
    global_input = incoming.u@w@dagger(incoming.u)
    if sar:
        global_input = global_input.T
    actual = outgoing.encode(extend_left(old, incoming.decode(global_input), a, r))
    expected = outgoing.u@extend_left(standard, w, a, r)@dagger(outgoing.u)
    if sbr:
        expected = expected.T
    c = choi(fixed, a)
    rebuilt = kraus_from_choi(c, a)
    errors = [np.linalg.norm(fixed(x)-apply_kraus(rebuilt, x)) for x in hermitian_basis(a)]
    return {'input_dim': a, 'output_dim': b, 'chart_bits': list(bits),
            'extension_error': float(np.linalg.norm(actual-expected)),
            'choi_min_eigenvalue': float(np.linalg.eigvalsh(c).min()),
            'kraus_reconstruction_error': float(max(errors)),
            'tp_error': float(np.linalg.norm(sum(dagger(k)@k for k in rebuilt)-np.eye(a))),
            'output_min_eigenvalue': float(np.linalg.eigvalsh(actual).min())}


def flatten(tree):
    if isinstance(tree, str):
        return tree
    return flatten(tree[0])+flatten(tree[1])


def all_trees(word):
    if len(word) == 1:
        return [word]
    return [(a, b) for k in range(1, len(word))
            for a in all_trees(word[:k]) for b in all_trees(word[k:])]


def object_unitary(word):
    return unitary(sample_generator(2**len(word), seed=2252+sum(map(ord, word))), .117)


def tree_unitary(tree):
    if isinstance(tree, str):
        return np.eye(2, dtype=complex)
    left, right = tree
    wl, wr = flatten(left), flatten(right)
    # Arbitrary implementer phases deliberately differ across parenthesizations.
    phase = np.exp(1j*.13*(len(wl)+2*len(wr)))
    join = phase*object_unitary(wl+wr) @ np.kron(dagger(object_unitary(wl)), dagger(object_unitary(wr)))
    return join @ np.kron(tree_unitary(left), tree_unitary(right))


def coherence_certificate():
    us = [tree_unitary(t) for t in all_trees('ABCD')]
    n = us[0].shape[0]
    errors = []
    phases = []
    for u in us:
        ratio = u@dagger(us[0])
        phase = np.trace(ratio)/n
        phases.append(float(np.angle(phase)))
        errors.append(float(np.linalg.norm(ratio-phase*np.eye(n))))
    return {'four_factor_parenthesizations': len(us),
            'implementer_relative_phases': phases,
            'max_noncentral_difference': max(errors),
            'meaning': 'Constructed coherent model; scalar phase differences leave all density-operator maps identical.'}


class OrientationTests(unittest.TestCase):
    def test_triangle_constraints_exactly_match_global_gauges(self):
        for n in range(2, 6):
            c = parity_certificate(n)
            self.assertEqual(c['consistent'], c['expected_consistent'])

    def test_odd_triangle_rejected(self):
        with self.assertRaises(ValueError):
            solve_gauge(edge_bits(3, (1, 1, 1)))

    def test_reference_change_is_one_global_flip(self):
        s = np.array([1, 0, 1, 1, 0])
        tau = s[:, None] ^ s[None, :]
        for r in range(5):
            self.assertTrue(np.array_equal(solve_gauge(tau, r), s ^ s[r]))

    def test_all_pair_marginals_including_middle_discard_detect_wrong_parity(self):
        c = marginal_sign_certificate()
        self.assertEqual(c['cases'], 96)
        self.assertEqual(c['consistent_cases'], 48)
        self.assertEqual(c['nonadjacent_pair_cases'], 32)

    def test_embedding_signs_telescope_through_composite_objects(self):
        for sa, sb, sab, sr, sabr in itertools.product((0, 1), repeat=5):
            ea, eb = sa ^ sab, sb ^ sab
            outer, reference = sab ^ sabr, sr ^ sabr
            self.assertEqual(ea ^ outer ^ reference, sa ^ sr)
            self.assertEqual(eb ^ outer ^ reference, sb ^ sr)

    def test_regauging_removes_both_local_transposes(self):
        x = sample_state(6, 2253)
        for sa, sb, sab in itertools.product((0, 1), repeat=3):
            chart = PairChart(2, 3, sa, sb, sab)
            result = chart.encode(partial_transpose(x, 2, 3, sa, sb))
            if sab:
                result = result.T
            np.testing.assert_allclose(result, chart.u@x@dagger(chart.u), atol=2e-14)

    def test_rectangular_channels_in_all_input_output_and_partner_charts(self):
        for a, b in ((2, 3), (3, 2), (2, 4)):
            for bits in itertools.product((0, 1), repeat=5):
                c = cross_channel_certificate(a, b, bits)
                self.assertLess(c['extension_error'], 3e-13)
                self.assertLess(c['kraus_reconstruction_error'], 3e-13)
                self.assertLess(c['tp_error'], 3e-13)
                self.assertGreaterEqual(c['choi_min_eigenvalue'], -3e-13)
                self.assertGreaterEqual(c['output_min_eigenvalue'], -3e-13)

    def test_bad_cross_type_chart_can_make_identity_look_like_transpose(self):
        old = lambda x: x.T
        fixed = lambda x: old(x).T
        self.assertAlmostEqual(np.linalg.eigvalsh(choi(old, 2)/2).min(), -.5)
        self.assertGreaterEqual(np.linalg.eigvalsh(choi(fixed, 2)/2).min(), -1e-14)

    def test_four_factor_coherence_is_at_the_operator_level(self):
        c = coherence_certificate()
        self.assertEqual(c['four_factor_parenthesizations'], 5)
        self.assertLess(c['max_noncentral_difference'], 3e-13)
        self.assertGreater(max(abs(p) for p in c['implementer_relative_phases']), .1)


def report():
    cross = [cross_channel_certificate(a, b, bits)
             for a, b in ((2, 3), (3, 2), (2, 4))
             for bits in itertools.product((0, 1), repeat=5)]
    return {'round': 225, 'date': '2026-09-20',
            'analytic_result': 'Positive three-body marginals force even relative-orientation parity. A fixed nontrivial reference type gives one coherent orientation for all objects, including composites. All actual A-to-B events are CP in these charts.',
            'framework_used': ['associative parallel composition', 'discarding any factor including the middle one', 'compatibility of marginal probabilities and product effects', 'closure under required finite composites'],
            'additional_physical_axioms': [],
            'parity_enumeration': [parity_certificate(n) for n in range(2, 6)],
            'marginal_sign_witness': marginal_sign_certificate(),
            'cross_system_certificates': {'cases': len(cross),
                'max_extension_error': max(c['extension_error'] for c in cross),
                'max_kraus_reconstruction_error': max(c['kraus_reconstruction_error'] for c in cross),
                'max_tp_error': max(c['tp_error'] for c in cross),
                'min_choi_eigenvalue': min(c['choi_min_eigenvalue'] for c in cross),
                'min_output_eigenvalue': min(c['output_min_eigenvalue'] for c in cross),
                'dimensions': [[2, 3], [3, 2], [2, 4]], 'chart_assignments_per_dimension': 32},
            'coherence': coherence_certificate(),
            'scope_limits': ['no proof all CP maps are physically allowed', 'no requirement all integer Hilbert dimensions occur', 'operator-level coherence, not a chosen phase-coherent lift to pure vectors', 'no particular Hamiltonian or gravity derived'],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'next': 'Use U in the coherent tensor representation to audit steering-state purity and measurement permissions.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(OrientationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('global_orientation_bridge_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
