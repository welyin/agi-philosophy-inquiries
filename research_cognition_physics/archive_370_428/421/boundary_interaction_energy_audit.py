"""Round 421: exact all-input gates and the interaction-energy boundary.

Only Python/NumPy.  --check is read-only; default creates a new result file.
The declared data Hamiltonian and the actual round-420 tensor factors are fixed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
from local_hardware_family_audit import Device, I2, X, Y, Z, embed_data, exp_h, opnorm

OBS = {}


def beta(a):
    values = np.linalg.eigvalsh(a)
    return float((values[-1] - values[0]) / 2)


def trace_distance(a, b):
    return float(np.linalg.svd(a - b, compute_uv=False).sum() / 2)


def pure_choi(u):
    vector = u.reshape(-1) / math.sqrt(u.shape[0])
    return np.outer(vector, vector.conj())


def isometry_choi(v, da, db):
    # v has row order (output A, environment B), column order input A.
    vectors = v.reshape(da, db, da).transpose(0, 2, 1).reshape(da * da, db)
    return vectors @ vectors.conj().T / da


def ladder(m):
    assert m >= 2
    bdim = m + 2
    w = np.eye(2 * bdim, dtype=complex)
    for n in range(1, bdim):
        j, k = n, bdim + n - 1
        w[j, j] = w[k, k] = 0
        w[j, k] = w[k, j] = 1
    charge = np.kron(np.diag([0., 1.]), np.eye(bdim)) + np.kron(I2, np.diag(np.arange(bdim)))
    battery = np.zeros(bdim, complex)
    battery[1:m + 1] = 1 / math.sqrt(m)
    injection = np.kron(I2, battery[:, None])
    return w, charge, battery, injection


class Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        gate = exp_h(np.kron(X, Y) + .3 * np.kron(Z, I2), math.sqrt(2) / 7)
        old = [((0,), .31 * Z), ((1,), -.23 * Y), ((0, 1), .47 * np.kron(X, X))]
        cls.device = Device((2, 2), old, (0, 1, 1), [((0, 1), np.eye(4)), ((0, 1), gate)])
        d = cls.device
        cls.h = d.apply(np.eye(d.size, dtype=complex))
        cls.ha_full = d.apply(np.eye(d.size, dtype=complex), d.raw_old)
        cls.ha = sum((embed_data(a, s, d.data_dims) for s, a in old), np.zeros((d.D, d.D), complex))
        cls.v = cls.h - cls.ha_full
        cls.b0, cls.b1 = d.inject(0, prefix=False), d.inject(d.L, prefix=False)

    def test_01_actual_endpoint_with_complete_reference(self):
        d = self.device
        endpoint = exp_h(self.h, d.time) @ self.b0
        wanted = (-1j) ** d.L * self.b1 @ d.target
        residual = opnorm(endpoint - wanted)
        self.assertLess(residual, 1e-12)
        bell = np.eye(d.D).reshape(-1) / math.sqrt(d.D)
        reference_error = np.linalg.norm(np.kron(endpoint - wanted, np.eye(d.D)) @ bell)
        self.assertLess(reference_error, 1e-12)
        OBS['actual_h_dimension'] = d.size
        OBS['all_input_endpoint_error'] = residual
        OBS['complete_reference_error'] = float(reference_error)

    def test_02_boundary_energy_operator_identity(self):
        u = self.device.target
        a0, a1 = self.b0.conj().T @ self.v @ self.b0, self.b1.conj().T @ self.v @ self.b1
        delta = u.conj().T @ self.ha @ u - self.ha
        residual = opnorm(delta - (a0 - u.conj().T @ a1 @ u))
        self.assertLess(residual, 1e-12)
        self.assertLess(opnorm(a0 + self.ha), 1e-12)
        self.assertLess(opnorm(a1 + self.ha), 1e-12)
        self.assertGreater(beta(delta), .1)
        self.assertLessEqual(beta(delta), beta(a0) + beta(a1) + 1e-12)
        OBS['boundary_identity_error'] = residual
        OBS['data_energy_nonscalar_change'] = beta(delta)
        OBS['sum_boundary_energy_variations'] = beta(a0) + beta(a1)

    def test_03_sharp_boundary_witness(self):
        d = Device((2,), [((0,), Z)], (0, 0), [((0,), X)])
        h = d.apply(np.eye(d.size, dtype=complex))
        hf = d.apply(np.eye(d.size, dtype=complex), d.raw_old)
        v = h - hf
        b0, b1 = d.inject(0, prefix=False), d.inject(1, prefix=False)
        a0, a1 = b0.conj().T @ v @ b0, b1.conj().T @ v @ b1
        delta = X @ Z @ X - Z
        self.assertAlmostEqual(beta(delta), beta(a0) + beta(a1))
        actual = exp_h(h, d.time) @ b0
        data_change = np.real(np.diag(actual.conj().T @ hf @ actual - b0.conj().T @ hf @ b0))
        interaction_change = np.real(np.diag(actual.conj().T @ v @ actual - b0.conj().T @ v @ b0))
        self.assertTrue(np.allclose(data_change, [-2, 2]))
        self.assertTrue(np.allclose(data_change + interaction_change, 0))
        OBS['sharp_beta_delta'] = beta(delta)
        OBS['sharp_beta_boundary_sum'] = beta(a0) + beta(a1)
        OBS['opposite_data_energy_changes'] = data_change.tolist()
        OBS['interaction_energy_balance_error'] = float(np.max(np.abs(data_change + interaction_change)))

    def test_04_mixed_consumable_auxiliary_needs_no_return(self):
        # Append an independent mixed spectator to an actual gate device.
        d = Device((2,), [((0,), Z)], (0, 0), [((0,), X)])
        h = d.apply(np.eye(d.size, dtype=complex))
        hs = .37 * Y
        sigma = np.diag([.8, .2])
        vs = exp_h(hs, d.time)
        tau = vs @ sigma @ vs.conj().T
        rho = np.array([[.4, .12 + .17j], [.12 - .17j, .6]])
        b0, b1 = d.inject(0, prefix=False), d.inject(1, prefix=False)
        w = exp_h(np.kron(h, I2) + np.kron(np.eye(d.size), hs), d.time)
        initial = np.kron(b0 @ rho @ b0.conj().T, sigma)
        expected = np.kron(b1 @ X @ rho @ X @ b1.conj().T, tau)
        residual = opnorm(w @ initial @ w.conj().T - expected)
        self.assertLess(residual, 1e-12)
        self.assertGreater(trace_distance(sigma, tau), .2)
        self.assertAlmostEqual(float(np.trace(hs @ sigma).real), float(np.trace(hs @ tau).real))
        OBS['mixed_auxiliary_factorization_error'] = residual
        OBS['mixed_auxiliary_nonreturn_distance'] = trace_distance(sigma, tau)

    def test_05_ladder_charge_and_autonomous_witness(self):
        worst = 0.
        for m in (2, 3, 8):
            w, q, _, _ = ladder(m)
            self.assertLess(opnorm(w @ w - np.eye(len(w))), 1e-12)
            self.assertLess(opnorm(w @ q - q @ w), 1e-12)
            # Q is a conserved charge, not automatically the generator.
            h = q + (np.eye(len(w)) - w) / 4
            self.assertGreaterEqual(float(np.linalg.eigvalsh(h).min()), -1e-12)
            worst = max(worst, opnorm(exp_h(h, 2 * math.pi) - w))
        self.assertLess(worst, 1e-12)
        OBS['ladder_autonomous_endpoint_error'] = worst

    def test_06_ladder_full_channel_and_reference_error(self):
        table = []
        for m in (2, 3, 8, 32):
            w, _, _, b = ladder(m)
            choi = isometry_choi(w @ b, 2, m + 2)
            ideal = pure_choi(X)
            predicted = (1 - 1 / m) * ideal + pure_choi(X @ Z) / m
            self.assertLess(opnorm(choi - predicted), 1e-12)
            bell_lower_bound = trace_distance(choi, ideal)
            self.assertAlmostEqual(bell_lower_bound, 1 / m)
            # Analytic mixture upper bound is 1/M, including all references.
            table.append(dict(M=m, half_diamond_error=1 / m, bell_witness=bell_lower_bound))
        OBS['ladder_reference_error_table'] = table

    def test_07_ladder_resources_and_input_dependent_battery(self):
        table = []
        for m in (2, 3, 8, 32):
            w, q, battery, b = ladder(m)
            n = np.arange(m + 2)
            mean = float(np.vdot(battery, n * battery).real)
            var = float(np.vdot(battery, (n - mean) ** 2 * battery).real)
            out = (w @ b).reshape(2, m + 2, 2)
            final0, final1 = out[1, :, 0], out[0, :, 1]
            overlap = float(np.vdot(final0, final1).real)
            self.assertAlmostEqual(overlap, (m - 2) / m)
            self.assertAlmostEqual(mean, (m + 1) / 2)
            self.assertAlmostEqual(var, (m * m - 1) / 12)
            self.assertAlmostEqual(float(np.vdot(final0, n * final0).real) - mean, -1)
            self.assertAlmostEqual(float(np.vdot(final1, n * final1).real) - mean, 1)
            self.assertLess(opnorm(b.conj().T @ (w.conj().T @ q @ w - q) @ b), 1e-12)
            table.append(dict(M=m, mean_charge=mean, charge_variance=var, conditional_overlap=overlap))
        OBS['ladder_resource_table'] = table

    def test_08_incoherent_uniform_battery_negative_control(self):
        m = 8
        w, _, _, _ = ladder(m)
        choi = np.zeros((4, 4), complex)
        for n in range(1, m + 1):
            ket = np.eye(m + 2)[:, n:n + 1]
            choi += isometry_choi(w @ np.kron(I2, ket), 2, m + 2) / m
        expected = (pure_choi(X) + pure_choi(X @ Z)) / 2
        self.assertLess(opnorm(choi - expected), 1e-12)
        self.assertAlmostEqual(trace_distance(choi, pure_choi(X)), .5)
        OBS['incoherent_battery_half_diamond_error'] = .5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise SystemExit(1)
    saved = dict(round=421, baseline_round=420, status='verified_conditional_result',
                 python=platform.python_version(), numpy=np.__version__, tests_run=result.testsRun,
                 failures=len(result.failures), errors=len(result.errors), observations=OBS,
                 scope=dict(fixed_energy_decomposition_is_model_input=True,
                            scalar_interaction_boundaries_are_extra_contract=True,
                            actual_round420_hamiltonian_checked=True,
                            auxiliary_catalytic_return_required=False,
                            permanent_output_required=False,
                            approximate_control_excluded=False,
                            optimal_resource_cost_proved=False,
                            three_dimensional_space_derived=False,
                            full_cognitive_countermodel_completed=False,
                            full_GR_goal_completed=False, phase_closure_triggered=False))
    if not args.check:
        with Path(__file__).with_name('boundary_interaction_energy_audit_results.json').open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(saved, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(saved, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
