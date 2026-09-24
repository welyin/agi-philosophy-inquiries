"""Round 420: physical local hardware extensions, not a fixed processor.

Python/NumPy only.  Run --check to avoid writing; default exclusively creates JSON.
Every datum is in an original physical tensor factor.  No circuit-dependent data
encoding is used.  Head alphabet, couplings, and extension hardware are resources.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
SWAP = np.eye(4, dtype=complex)[[0, 2, 1, 3]]
OBS: dict[str, object] = {}


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * w)) @ v.conj().T


def apply_local(a, sites, dims, columns):
    """Apply a literal one/two-cell matrix to columns in full physical space."""
    sites = tuple(sites)
    other = tuple(i for i in range(len(dims)) if i not in sites)
    order = sites + other + (len(dims),)
    shaped = columns.reshape(tuple(dims) + (columns.shape[1],))
    moved = shaped.transpose(order)
    out = a @ moved.reshape(a.shape[0], -1)
    out = out.reshape(tuple(dims[i] for i in sites + other) + (columns.shape[1],))
    return out.transpose(np.argsort(order)).reshape(columns.shape)


def embed_data(a, sites, data_dims):
    return apply_local(a, sites, data_dims, np.eye(math.prod(data_dims), dtype=complex))


class Device:
    """Each cell is enable E (2), head Q (power of 2), original D."""

    def __init__(self, data_dims, old_terms, head_positions, gates, gmax=1.):
        self.data_dims = tuple(data_dims)
        self.n = len(data_dims)
        self.L = len(gates)
        assert self.L >= 1 and len(head_positions) == self.L + 1
        self.Q = 1 << (self.L + 1).bit_length()
        self.dims = tuple(2 * self.Q * d for d in data_dims)
        self.size = math.prod(self.dims)
        self.D = math.prod(data_dims)
        self.pos = tuple(head_positions)
        mL = math.sqrt(((self.L + 1) ** 2) // 4)
        self.omega = 2 * gmax / mL
        self.time = math.pi / self.omega
        self.weights = [self.omega * math.sqrt((t + 1) * (self.L - t)) / 2
                        for t in range(self.L)]
        self.terms = []
        self.old_terms = old_terms
        self.raw_old = []
        self.compensation = []
        self.transitions = []
        self.target = np.eye(self.D, dtype=complex)
        self.prefix = [self.target]
        for sites, a in old_terms:
            raw = self.lift(a, sites)
            correction = -self.lift(a, sites, enabled=True)
            self.raw_old.append((tuple(sites), raw))
            self.compensation.append((tuple(sites), correction))
            self.terms.extend([(tuple(sites), raw), (tuple(sites), correction)])
        for t, (sites, gate) in enumerate(gates):
            sites = tuple(sites)
            x, y = self.pos[t:t + 2]
            if x != y:
                assert abs(x - y) == 1 and sites == tuple(sorted((x, y)))
                assert np.allclose(gate, np.eye(gate.shape[0]))
            else:
                assert x in sites
            assert len(sites) <= 2
            assert len(sites) == 1 or sites[1] - sites[0] == 1
            source = [0] * len(sites)
            target = [0] * len(sites)
            source[sites.index(x)] = t + 1
            target[sites.index(y)] = t + 2
            forward = self.lift(gate, sites, head_pair=(source, target), enabled=True)
            h = self.weights[t] * (forward + forward.conj().T)
            self.transitions.append((sites, h))
            self.terms.append((sites, h))
            self.target = embed_data(gate, sites, data_dims) @ self.target
            self.prefix.append(self.target)

    def lift(self, a, sites, head_pair=None, enabled=False):
        sites = tuple(sites)
        dd = tuple(self.data_dims[i] for i in sites)
        physical = tuple(self.dims[i] for i in sites)
        out = np.zeros((math.prod(physical), math.prod(physical)), complex)
        data_states = list(itertools.product(*(range(d) for d in dd)))
        enables = [tuple([1] * len(sites))] if enabled else itertools.product((0, 1), repeat=len(sites))
        for es in enables:
            heads = [(head_pair[0], head_pair[1])] if head_pair is not None else [
                (q, q) for q in itertools.product(range(self.Q), repeat=len(sites))]
            for qs, qt in heads:
                source = [np.ravel_multi_index(tuple((es[k] * self.Q + qs[k]) * dd[k] + b[k]
                                                     for k in range(len(sites))), physical)
                          for b in data_states]
                target = [np.ravel_multi_index(tuple((es[k] * self.Q + qt[k]) * dd[k] + b[k]
                                                     for k in range(len(sites))), physical)
                          for b in data_states]
                out[np.ix_(target, source)] += a
        return out

    def apply(self, columns, terms=None):
        out = np.zeros_like(columns, dtype=complex)
        for sites, a in self.terms if terms is None else terms:
            out += apply_local(a, sites, self.dims, columns)
        return out

    def inject(self, t, enable=1, prefix=True):
        out = np.zeros((self.size, self.D), complex)
        for b in itertools.product(*(range(d) for d in self.data_dims)):
            physical = tuple((enable * self.Q + (t + 1 if i == self.pos[t] else 0))
                             * self.data_dims[i] + b[i] for i in range(self.n))
            out[np.ravel_multi_index(physical, self.dims), np.ravel_multi_index(b, self.data_dims)] = 1
        return out @ self.prefix[t] if prefix else out

    def history(self):
        return np.column_stack([self.inject(t) for t in range(self.L + 1)])

    def clock(self):
        f = np.zeros((self.L + 1, self.L + 1), complex)
        for t, g in enumerate(self.weights):
            f[t + 1, t] = f[t, t + 1] = g
        return f


class Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        theta = math.sqrt(2.) / 7
        gate = exp_h(np.kron(X, Y) + .3 * np.kron(Z, I2), theta)
        old = [((0,), .31 * Z), ((1,), -.23 * Y), ((0, 1), .47 * np.kron(X, X))]
        cls.base = Device((2, 2), old, (0, 1, 1), [((0, 1), np.eye(4)), ((0, 1), gate)])
        cls.H = cls.base.apply(np.eye(cls.base.size, dtype=complex))
        cls.U = exp_h(cls.H, cls.base.time)

    def test_01_literal_local_rule_and_compensation(self):
        d = self.base
        self.assertLess(opnorm(self.H - self.H.conj().T), 1e-12)
        b = d.inject(0, prefix=False)
        cancellation = d.apply(b, d.raw_old + d.compensation)
        self.assertLess(opnorm(cancellation), 1e-12)
        self.assertGreater(opnorm(d.apply(b, d.raw_old)), .2)
        OBS['literal_physical_dimension'] = d.size
        OBS['drift_cancel_residual'] = opnorm(cancellation)
        OBS['uncancelled_drift_norm'] = opnorm(d.apply(b, d.raw_old))

    def test_02_disabled_preserves_original_natural_h(self):
        d = self.base
        b = d.inject(0, enable=0, prefix=False)
        old = sum((embed_data(a, s, d.data_dims) for s, a in d.old_terms), np.zeros((d.D, d.D), complex))
        error = opnorm(d.apply(b) - b @ old)
        self.assertLess(error, 1e-12)
        evolved = exp_h(self.H, .83) @ b
        self.assertLess(opnorm(evolved - b @ exp_h(old, .83)), 1e-11)
        OBS['disabled_original_h_error'] = error

    def test_03_full_unknown_reference_intertwining(self):
        d = self.base
        j = d.history()
        error = opnorm(d.apply(j) - j @ np.kron(d.clock(), np.eye(d.D)))
        self.assertLess(opnorm(j.conj().T @ j - np.eye(j.shape[1])), 1e-12)
        self.assertLess(error, 1e-11)
        endpoint = self.U @ d.inject(0)
        expected = (-1j) ** d.L * d.inject(d.L)
        exact_error = opnorm(endpoint - expected)
        self.assertLess(exact_error, 1e-11)
        # Vectorization is a maximally entangled input with a untouched reference.
        bell_error = float(np.linalg.norm((endpoint - expected).reshape(-1) / math.sqrt(d.D)))
        OBS['all_input_isometry_error'] = exact_error
        OBS['reference_vector_error'] = bell_error
        OBS['intertwining_error'] = error

    def test_04_no_compensation_negative_control(self):
        d = self.base
        bad = d.apply(np.eye(d.size, dtype=complex), d.raw_old + d.transitions)
        error = opnorm(exp_h(bad, d.time) @ d.inject(0) - (-1j) ** d.L * d.inject(d.L))
        self.assertGreater(error, .1)
        OBS['uncompensated_endpoint_error'] = error

    def test_05_unknown_new_subject_actual_swap(self):
        # A at cell 0, empty work register at 1, genuinely unknown B at 2.
        old = [((0,), .2 * X), ((1, 2), .3 * np.kron(Z, Y))]
        d = Device((2, 2, 2), old, (1, 1), [((1, 2), SWAP)])
        j = d.history()
        self.assertLess(opnorm(d.apply(j) - j @ np.kron(d.clock(), np.eye(8))), 1e-11)
        incoming = np.zeros((8, 4), complex)
        expected = np.zeros((8, 4), complex)
        for a, b in itertools.product((0, 1), repeat=2):
            incoming[4 * a + b, 2 * a + b] = 1
            expected[4 * a + 2 * b, 2 * a + b] = 1
        self.assertLess(opnorm(d.target @ incoming - expected), 1e-12)
        # The full isometry preserves *all* A-B-reference inputs, including entangled A,B.
        end = d.inject(1) @ incoming
        direct = d.inject(1, prefix=False) @ expected
        self.assertLess(opnorm(end - direct), 1e-12)
        OBS['new_subject_physical_factor_isometry_error'] = opnorm(end - direct)
        OBS['new_subject_old_factor_preserved'] = True

    def test_06_next_device_controls_entire_previous_device(self):
        # First device has a complete 16-D cell, including unused head/enable sectors.
        first = Device((2,), [((0,), .17 * X)], (0, 0), [((0,), exp_h(Y, .37))])
        old_h = first.apply(np.eye(first.size, dtype=complex))
        rng = np.random.default_rng(420)
        q, r = np.linalg.qr(rng.normal(size=(first.size, first.size)) +
                            1j * rng.normal(size=(first.size, first.size)))
        target = q @ np.diag(np.diag(r) / np.abs(np.diag(r)))
        outer = Device((first.size,), [((0,), old_h)], (0, 0), [((0,), target)])
        all_inputs = outer.history()
        error = opnorm(outer.apply(all_inputs) - all_inputs @ np.kron(outer.clock(), np.eye(first.size)))
        self.assertLess(error, 1e-11)
        h = outer.apply(np.eye(outer.size, dtype=complex))
        end_error = opnorm(exp_h(h, outer.time) @ outer.inject(0) + 1j * outer.inject(1))
        self.assertLess(end_error, 1e-10)
        OBS['complete_old_device_dimension'] = first.size
        OBS['next_level_dimension'] = outer.size
        OBS['all_old_sectors_control_error'] = end_error

    def test_07_independent_composition(self):
        ha, hb, hc = .23 * X, .19 * Y, -.41 * Z
        hab = np.kron(ha, I2) + np.kron(I2, hb)
        left = np.kron(hab, I2) + np.kron(np.eye(4), hc)
        hbc = np.kron(hb, I2) + np.kron(I2, hc)
        right = np.kron(ha, np.eye(4)) + np.kron(I2, hbc)
        self.assertLess(opnorm(left - right), 1e-12)
        error = opnorm(exp_h(left, .79) - np.kron(np.kron(exp_h(ha, .79), exp_h(hb, .79)), exp_h(hc, .79)))
        self.assertLess(error, 1e-12)
        OBS['independent_three_type_composition_error'] = error

    def test_08_internal_record_and_feedback(self):
        # One complete physical cell consists of data D, log M, environment E, reader R.
        dims = (2, 2, 2, 2)
        cnot = np.eye(4, dtype=complex)[[0, 1, 3, 2]]
        def cg(a, b):
            return embed_data(cnot, (a, b), dims)
        target = cg(1, 0) @ cg(1, 3) @ cg(1, 2) @ cg(0, 1)
        old = .07 * embed_data(X, (0,), dims) + .03 * embed_data(Y, (2,), dims)
        d = Device((16,), [((0,), old)], (0, 0), [((0,), target)])
        j = d.history()
        self.assertLess(opnorm(d.apply(j) - j @ np.kron(d.clock(), np.eye(16))), 1e-12)
        start = np.zeros((16, 2), complex)
        start[0, 0], start[8, 1] = 1, 1
        expected = np.zeros((16, 2), complex)
        expected[0, 0], expected[7, 1] = 1, 1
        self.assertLess(opnorm(target @ start - expected), 1e-12)
        # After tracing E, outcome-labelled branches send |b> to |0,b,b> on D,M,R.
        v = (target @ start).reshape(2, 2, 2, 2, 2)
        kraus = [v[:, :, e, :, :].reshape(8, 2) for e in range(2)]
        comp = sum((k.conj().T @ k for k in kraus), np.zeros((2, 2), complex))
        self.assertLess(opnorm(comp - I2), 1e-12)
        OBS['record_reader_feedback_isometry_error'] = opnorm(target @ start - expected)
        OBS['record_environment_physically_retained'] = True

    def test_09_local_strength_and_initial_independence(self):
        d = self.base
        by_support = {}
        for sites, a in d.transitions:
            by_support[sites] = by_support.get(sites, np.zeros_like(a)) + a
        maximum = max(opnorm(a) for a in by_support.values())
        self.assertLessEqual(maximum, 2 * max(d.weights) + 1e-12)
        self.assertTrue(all(len(s) <= 2 and (len(s) == 1 or s[1] == s[0] + 1) for s, _ in d.terms))
        b = d.inject(0)
        self.assertLess(opnorm(b.conj().T @ b - np.eye(d.D)), 1e-12)
        OBS['max_compiled_head_term_norm'] = maximum
        OBS['max_transition_weight'] = max(d.weights)
        OBS['runtime'] = d.time
        OBS['initial_auxiliary_does_not_depend_on_unknown_input'] = True

    def test_10_endpoint_is_not_release_of_original_device(self):
        d = self.base
        endpoint = d.inject(d.L, prefix=False)
        old = sum((embed_data(a, s, d.data_dims) for s, a in d.old_terms),
                  np.zeros((d.D, d.D), complex))
        release_defect = opnorm(d.apply(endpoint) - endpoint @ old)
        self.assertGreater(release_defect, .1)
        reduced_generator = endpoint.conj().T @ d.apply(endpoint)
        rho = np.zeros((d.D, d.D), complex)
        rho[0, 0] = 1
        actual_derivative = -1j * (reduced_generator @ rho - rho @ reduced_generator)
        original_derivative = -1j * (old @ rho - rho @ old)
        derivative_defect = opnorm(actual_derivative - original_derivative)
        self.assertLess(opnorm(actual_derivative), 1e-12)
        self.assertGreater(derivative_defect, .1)
        # The final controller remains on: the old drift remains cancelled,
        # and the final head has a backward transition.  No switch-off occurs.
        for i in range(d.n):
            enabled = np.diag(np.repeat([0., 1.], d.Q * d.data_dims[i]))
            p = apply_local(enabled, (i,), d.dims, np.eye(d.size, dtype=complex))
            self.assertLess(opnorm(self.H @ p - p @ self.H), 1e-12)
            self.assertLess(opnorm(p @ endpoint - endpoint), 1e-12)
        OBS['endpoint_original_free_evolution_defect'] = release_defect
        OBS['endpoint_reduced_data_generator_norm'] = opnorm(reduced_generator)
        OBS['original_data_derivative_defect'] = derivative_defect
        OBS['original_device_automatically_released'] = False


    def test_11_same_endpoint_channel_distinct_future_processes(self):
        # Same literal apparatus and unknown-input factor, different independent
        # enable preparation.  All-input endpoint maps coincide, later maps do not.
        d = Device((2,), [((0,), Z)], (0, 0), [((0,), I2)], gmax=.5)
        h = d.apply(np.eye(d.size, dtype=complex))
        self.assertAlmostEqual(d.time, math.pi)
        bell = I2.reshape(-1) / math.sqrt(2)
        identity_choi = np.outer(bell, bell.conj())
        def kraus_at(t, enable):
            v = exp_h(h, t) @ d.inject(0, enable=enable, prefix=False)
            return v.reshape(2 * d.Q, 2, 2)
        endpoint_errors = []
        for enable in (0, 1):
            ks = kraus_at(d.time, enable)
            vectors = ks.reshape(2 * d.Q, 4)
            choi = np.einsum('ai,aj->ij', vectors, vectors.conj()) / 2
            endpoint_errors.append(opnorm(choi - identity_choi))
        self.assertLess(max(endpoint_errors), 1e-12)
        plus = np.ones(2, complex) / math.sqrt(2)
        rho = np.outer(plus, plus.conj())
        later = []
        for enable in (0, 1):
            ks = kraus_at(d.time + math.pi / 2, enable)
            later.append(sum((k @ rho @ k.conj().T for k in ks), np.zeros((2, 2), complex)))
        distance = float(np.sum(np.abs(np.linalg.eigvalsh(later[0] - later[1]))) / 2)
        self.assertAlmostEqual(distance, 1.)
        self.assertLess(opnorm(later[1] - rho), 1e-12)
        self.assertLess(opnorm(later[0] - Z @ rho @ Z), 1e-12)
        OBS['two_preparations_identity_endpoint_choi_errors'] = endpoint_errors
        OBS['same_endpoint_later_trace_distance'] = distance
        OBS['endpoint_channel_is_complete_future_process'] = False


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Audit)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    payload = {
        'round': 420, 'baseline_round': 418, 'status': 'verified_conditional_result',
        'python': platform.python_version(), 'numpy': np.__version__,
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'observations': OBS,
        'scope': {
            'fixed_local_rule_for_all_tasks': False,
            'exact_continuous_hardware_types_are_model_input': True,
            'declared_original_factors_and_arbitrary_reference': True,
            'all_complete_old_device_states_supported': True,
            'unknown_remote_arrival_generated': False,
            'continuous_fabrication_or_attachment_process_generated': False,
            'permanent_output_guarantee': False,
            'black_box_devices_physically_serial_composable': False,
            'endpoint_restores_original_free_hamiltonian': False,
            'three_dimensional_space_derived': False,
            'full_GR_goal_completed': False,
            'complete_local_process_FUCP_implementation_proved': False,
            'full_cognitive_countermodel_completed': False,
            'phase_closure_triggered': False,
        },
    }
    if not args.check:
        with Path(__file__).with_name('local_hardware_family_audit_results.json').open('x', encoding='utf-8') as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
            f.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    run()
