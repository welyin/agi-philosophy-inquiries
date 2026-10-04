"""Round 256: CHSH witness in typed interaction records, with no discarded trials.

This transfers supplied Bell correlations to local partner choices. It does
not derive quantum theory or geometry, and does not test all classical models.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
import mutual_proposal_interaction_audit as protocol

X, Z = protocol.X, protocol.Z
GRAPH = ((1,), (0, 2), (1,), (4,), (3, 5), (4,))
ALICE = (Z, X)
BOB = ((Z+X)/np.sqrt(2), (Z-X)/np.sqrt(2))
SIGNS = (1, -1)


def source(dephased=False):
    # Tensor order: leaf, center A, leaf, leaf, center B, leaf.
    vectors = [np.zeros(64), np.zeros(64)]
    for index, bits in enumerate(itertools.product((0, 1), repeat=6)):
        if bits[1] == bits[4]:
            vectors[bits[1]][index] = .25
    if dephased:
        return tuple((.5, psi) for psi in vectors)
    return ((1., (vectors[0]+vectors[1])/np.sqrt(2)),)


@lru_cache(maxsize=None)
def probabilities(eta, theta=np.pi/4, dephased=False):
    values = np.zeros((2, 2, 2, 2))
    for x, y in itertools.product(range(2), repeat=2):
        for ia, ib in itertools.product(range(2), repeat=2):
            record = (1, (0, 2)[ia], 1, 4, (3, 5)[ib], 4)
            h = protocol.branch(GRAPH, record, eta, theta, {1: ALICE[x], 4: BOB[y]})
            values[x, y, ia, ib] = sum(weight*float(np.vdot(h@psi, h@psi).real)
                                       for weight, psi in source(dephased))
    return values


def correlations(values):
    return np.einsum('xyab,a,b->xy', values, SIGNS, SIGNS)


def chsh(values):
    e = correlations(values)
    return float(e[0, 0]+e[0, 1]+e[1, 0]-e[1, 1])


def analytic(eta):
    expected = np.zeros((2, 2, 2, 2))
    for x, y, a, b in itertools.product(range(2), repeat=4):
        dot = float(np.trace(ALICE[x]@BOB[y]).real/2)
        expected[x, y, a, b] = (1+SIGNS[a]*SIGNS[b]*eta**2*dot)/4
    return expected


def classical_scores():
    return [a0*b0+a0*b1+a1*b0-a1*b1
            for a0, a1, b0, b1 in itertools.product(SIGNS, repeat=4)]


def entanglement_fidelity(eta):
    return float(sum(abs(np.trace(protocol.weak(Z, g, eta)))**2 for g in SIGNS)/4)


def report():
    cases = []
    for eta in (.6, .9, 1.):
        values = probabilities(eta)
        cases.append({'eta': eta, 'graph_CHSH': chsh(values),
                      'predicted_CHSH': float(2*np.sqrt(2)*eta**2),
                      'local_proposal_entanglement_fidelity': entanglement_fidelity(eta),
                      'probabilities_x_y_a_b': values.tolist(),
                      'analytic_probability_max_error': float(np.max(np.abs(values-analytic(eta))))})
    values = probabilities(.9)
    return {'round': 256, 'candidate_modules': 2, 'ports_per_module': 3,
            'retained_trial_fraction': 1., 'accepted_pair_events_per_trial': 2,
            'classical_deterministic_strategies': len(classical_scores()),
            'classical_CHSH_bound': max(classical_scores()),
            'violation_eta_threshold': float(2**(-.25)),
            'fidelity_at_threshold': float((1+np.sqrt(1-1/np.sqrt(2)))/2),
            'cases': cases,
            'normalization_max_error': float(np.max(np.abs(values.sum(axis=(2, 3))-1))),
            'no_signaling_max_error': float(max(
                np.max(np.abs(values.sum(axis=3)[:, 0]-values.sum(axis=3)[:, 1])),
                np.max(np.abs(values.sum(axis=2)[0]-values.sum(axis=2)[1])))),
            'dephased_source_CHSH_eta_09': chsh(probabilities(.9, dephased=True)),
            'scope': 'Supplied entanglement and quantum instruments produce CHSH violation in marked partner-choice records. The excluded class is setting-independent local classical shared-variable arbitration, with no intermodule communication during each trial. This is a finite simulation, not a spacetime derivation or experiment.'}


class Audit(unittest.TestCase):
    def test_01_every_trial_retained_and_both_modules_interact(self):
        rs = protocol.records(GRAPH)
        self.assertEqual(len(rs), 4)
        for r in rs:
            edges = protocol.matching(r)
            self.assertEqual(len(edges), 2)
            self.assertEqual(sum(1 in e for e in edges), 1)
            self.assertEqual(sum(4 in e for e in edges), 1)

    def test_02_full_six_qubits_agree_with_two_qubit_formula(self):
        for eta in (0., .6, .9, 1.):
            np.testing.assert_allclose(probabilities(eta), analytic(eta), atol=3e-14)

    def test_03_normalization_and_no_signaling(self):
        for eta in (0., .6, .9, 1.):
            p = probabilities(eta)
            self.assertGreaterEqual(float(p.min()), -1e-14)
            np.testing.assert_allclose(p.sum(axis=(2, 3)), 1., atol=3e-14)
            np.testing.assert_allclose(p.sum(axis=3), .5, atol=3e-14)
            np.testing.assert_allclose(p.sum(axis=2), .5, atol=3e-14)

    def test_04_exhaustive_classical_extreme_points(self):
        self.assertEqual(len(classical_scores()), 16)
        self.assertEqual(set(classical_scores()), {-2, 2})

    def test_05_violation_strength_and_nonviolation_control(self):
        self.assertLess(chsh(probabilities(.6)), 2)
        self.assertGreater(chsh(probabilities(.9)), 2)
        self.assertAlmostEqual(chsh(probabilities(1.)), 2*np.sqrt(2))
        for eta in (2**(-.25)-1e-5, 2**(-.25), 2**(-.25)+1e-5):
            self.assertAlmostEqual(chsh(probabilities(eta)), 2*np.sqrt(2)*eta**2)

    def test_06_post_record_pair_gates_do_not_select_samples(self):
        for theta in (0., .137, np.pi/4, 1.2):
            np.testing.assert_allclose(probabilities(.9, theta), probabilities(.9), atol=3e-14)

    def test_07_dephased_source_is_classical_control(self):
        value = chsh(probabilities(.9, dephased=True))
        self.assertAlmostEqual(value, np.sqrt(2)*.9**2)
        self.assertLess(value, 2)

    def test_08_measurement_disturbance_computed_from_kraus(self):
        for eta in (0., .6, .9, 1.):
            gamma = np.sqrt(1-eta**2)
            self.assertAlmostEqual(entanglement_fidelity(eta), (1+gamma)/2)
            channel_x = sum(protocol.weak(Z, g, eta)@X@protocol.weak(Z, g, eta) for g in SIGNS)
            np.testing.assert_allclose(channel_x, gamma*X, atol=3e-14)

    def test_09_settings_are_local_operator_choices(self):
        # All cross-module Kraus operators commute, even on entangled inputs.
        for x, y, a, b in itertools.product(range(2), repeat=4):
            left = protocol.embed(6, {1: protocol.weak(ALICE[x], SIGNS[a], .9)})
            right = protocol.embed(6, {4: protocol.weak(BOB[y], SIGNS[b], .9)})
            np.testing.assert_allclose(left@right, right@left, atol=3e-14)

    def test_10_source_resource_and_unbiased_settings_family(self):
        for dephased in (False, True):
            self.assertAlmostEqual(sum(w for w, _ in source(dephased)), 1.)
            for _, psi in source(dephased):
                self.assertAlmostEqual(float(np.vdot(psi, psi).real), 1.)
        for q in ALICE+BOB:
            np.testing.assert_allclose(q@q, np.eye(2), atol=3e-14)
            self.assertAlmostEqual(float(np.trace(q)), 0.)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python': platform.python_version(), 'numpy': np.__version__}
    result['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    target = Path(__file__).with_name('interaction_record_bell_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve existing results; use a new result version.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
