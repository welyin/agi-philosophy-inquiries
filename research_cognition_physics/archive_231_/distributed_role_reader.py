"""Round 507: distributed cat reader with the original graph dynamics left on.

Local occupation/readout coupling, fresh GHZ resources and timing are declared
inputs. Exact code preservation differs from approximation of ideal Luders read.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old
import joint_orbit_port_interface as orbit

HERE = Path(__file__).resolve().parent
TARGET = HERE/'distributed_role_reader_results.json'
OBS = {}


def exponential(a):
    values, vectors = np.linalg.eigh(a)
    return (vectors*np.exp(-1j*values)) @ vectors.conj().T


def operators(h, q, tau, perturbation=None):
    u0 = exponential(float(tau)*h)
    argument = float(tau)*h+np.pi*q
    if perturbation is not None:
        argument += float(tau)*perturbation
    u1 = exponential(argument)
    ks = [(u0+u1)/2, (u0-u1)/2]
    ideal = [(np.eye(len(q))-q) @ u0, q @ u0]
    return u0, u1, ks, ideal


def opnorm(a):
    return float(np.sqrt(max(0, np.linalg.eigvalsh(a.conj().T @ a)[-1])))


def choi_branch(k):
    v = k.ravel()/np.sqrt(len(k))
    return np.outer(v, v.conj())


class Audit(unittest.TestCase):
    def test_01_local_meter_tensor_restricts_to_two_cat_branches(self):
        h = old.model(2)[4]
        n, g, d, meter = 6, 6, 36, 64
        # Basis ordering: physical DG then all six meter bits.
        embedding = np.zeros((d*meter, d*2), dtype=np.int64)
        for x in range(d):
            embedding[x*meter, 2*x] = 1
            embedding[x*meter+meter-1, 2*x+1] = 1
        interaction_diagonal = np.zeros(d*meter, dtype=np.int64)
        for v in range(n):
            for j in range(g):
                for bits in range(meter):
                    interaction_diagonal[(v*g+j)*meter+bits] = int(v < 2 and bool(bits & (1 << (n-1-v))))
        local_action = interaction_diagonal[:, None]*embedding
        q = np.diag(np.repeat([1, 1, 0, 0, 0, 0], g))
        logical = np.kron(q, np.diag([0, 1]))
        self.assertTrue(np.array_equal(local_action, embedding @ logical))
        h_action = np.einsum('ab,bmc->amc', h, embedding.reshape(d, meter, d*2)).reshape(d*meter, d*2)
        self.assertTrue(np.array_equal(h_action, embedding @ np.kron(h, np.eye(2, dtype=np.int64))))
        OBS['local_tensor'] = dict(DG_dimension=d, meter_qubits=n, expanded_dimension=d*meter,
            cat_sector_dimension=d*2, integer_intertwining_exact=True,
            original_h_left_on=True, interaction_two_body_given_local_role=True)

    def test_02_exact_uniform_role_commutator(self):
        records = []
        for i in (2, 3, 4):
            trees = old.previous.ensemble(i)[0]
            n, leaves = 2*i+2, i+2
            for tree in trees:
                b = np.zeros((i, leaves), dtype=np.int64)
                for u, v in tree:
                    if u >= i and v < i:
                        u, v = v, u
                    if u < i <= v:
                        b[u, v-i] = 1
                self.assertTrue(np.array_equal(b.sum(axis=0), np.ones(leaves)))
                m = b.sum(axis=1)
                self.assertTrue(np.array_equal(b @ b.T, np.diag(m)))
                self.assertEqual(int(max(m)), 2)
            records.append(dict(I=i, trees=len(trees), incidence_gram_diagonal=True,
                commutator_norm_squared_exact=2))
        OBS['uniform_commutator'] = records

    def test_03_all_fine_meter_results_and_code_preservation(self):
        d = orbit.quotient(3)
        h, q = d['k'], np.diag(d['role'])
        _, _, ks, _ = operators(h, q, Q(1, 16))
        n = d['n']
        summed = [np.zeros((d['q']**2, d['q']**2), dtype=complex) for _ in range(2)]
        for bits in range(2**n):
            r = bits.bit_count() % 2
            fine = ks[r]/np.sqrt(2**(n-1))
            summed[r] += choi_branch(fine)
        for r in range(2):
            self.assertLess(np.linalg.norm(summed[r]-choi_branch(ks[r])), 1e-12)
        # Independently compare the full 720-dimensional action on every code column.
        full_h = old.model(3)[4]
        full_q = np.diag(np.repeat(np.arange(8) < 3, d['g']).astype(float))
        e = orbit.embedding(d)
        tau = Q(1, 16)
        u0e = old.evolve(float(tau)*full_h, e, 1, order=48)
        u1e = old.evolve(float(tau)*full_h+np.pi*full_q, e, 1, order=48)
        residual = max(np.linalg.norm((u0e+u1e)/2-e @ ks[0]), np.linalg.norm((u0e-u1e)/2-e @ ks[1]))
        self.assertLess(residual, 1e-11)
        self.assertLess(np.linalg.norm(sum(k.conj().T @ k for k in ks)-np.eye(d['q'])), 1e-12)
        OBS['complete_meter_instrument'] = dict(meter_strings=2**n,
            all_fine_results_proportional_to_parity_branch=True,
            full_unknown_code_columns_checked=True, full_code_intertwining_residual=old.short(residual),
            graph_measured=False, postselected=False)

    def test_04_uniform_stinespring_error_and_finite_histories(self):
        records = []
        for i in (2, 3, 4):
            d = orbit.quotient(i)
            h, q = d['k'], np.diag(d['role'])
            for tau in (Q(1, 8), Q(1, 32), Q(1, 128)):
                u0, u1, ks, ideal = operators(h, q, tau)
                delta = opnorm(u1-(np.eye(len(q))-2*q) @ u0)/np.sqrt(2)
                self.assertLessEqual(delta, np.pi*float(tau)/2+1e-12)
                observed = sum(float(np.sum(np.abs(np.linalg.eigvalsh(choi_branch(a)-choi_branch(b))))) for a, b in zip(ks, ideal))
                self.assertLessEqual(observed, 2*delta+1e-12)
                records.append(dict(I=i, tau=str(tau), dilation_error=old.short(delta),
                    ordinary_diamond_upper=old.short(np.pi*float(tau)),
                    normalized_choi_trace_distance=old.short(observed)))
        # Adaptive waits based on old recorded parity, retaining all eight histories.
        d = orbit.quotient(3)
        h, q = d['k'], np.diag(d['role'])
        paths = [([], np.eye(d['q'], dtype=complex), np.eye(d['q'], dtype=complex))]
        for step in range(3):
            new = []
            for history, actual, target in paths:
                tau = Q(1, 32) if not history or history[-1] == 0 else Q(1, 64)
                _, _, ks, ideal = operators(h, q, tau)
                for r in range(2):
                    new.append((history+[r], ks[r] @ actual, ideal[r] @ target))
            paths = new
        complete = sum(a.conj().T @ a for _, a, _ in paths)
        self.assertLess(np.linalg.norm(complete-np.eye(d['q'])), 1e-12)
        difference = sum(float(np.sum(np.abs(np.linalg.eigvalsh(choi_branch(a)-choi_branch(b))))) for _, a, b in paths)
        self.assertLess(difference, 3*np.pi/32)
        OBS['error_and_histories'] = dict(cases=records, adaptive_history_count=len(paths),
            full_history_choi_trace_distance=old.short(difference), full_history_diamond_upper=old.short(3*np.pi/32),
            unknown_reference_covered_by_dilation_proof=True)

    def test_05_classical_cat_populations_cannot_supply_role_record(self):
        d = orbit.quotient(3)
        h, q = d['k'], np.diag(d['role'])
        u0, u1, _, _ = operators(h, q, Q(1, 16))
        # In each parity branch the classical mixture has Kraus u0/2 and u1/2.
        effect = (u0.conj().T @ u0+u1.conj().T @ u1)/4
        self.assertLess(np.linalg.norm(effect-np.eye(d['q'])/2), 1e-12)
        coherent_cat = np.ones((2, 2))/2
        classical_cat = np.eye(2)/2
        self.assertTrue(np.array_equal(coherent_cat.diagonal(), classical_cat.diagonal()))
        self.assertAlmostEqual(float(np.sum(np.abs(np.linalg.eigvalsh(coherent_cat-classical_cat)))), 1)
        OBS['source_boundary'] = dict(same_meter_computational_populations=True,
            classical_cat_parity_probability_all_inputs='1/2', classical_cat_role_information=False,
            full_meter_string_probability_all_inputs='2^-N', no_fault_tolerance_claim=True)

    def test_06_finite_noise_budget_and_resources(self):
        h = old.model(2)[4]
        q = np.diag(np.repeat([1, 1, 0, 0, 0, 0], 6)).astype(float)
        tau, p = Q(1, 32), Q(1, 100)
        drift = np.diag(np.repeat([1/8, -1/8, 0, 0, 0, 0], 6))
        _, _, actual, ideal = operators(h, q, tau, drift)
        # A GHZ phase-error mixture interchanges parity branches with probability p.
        distance = 0.0
        for r in range(2):
            vectors = np.column_stack([actual[r].ravel(), actual[1-r].ravel(), ideal[r].ravel()])/np.sqrt(len(h))
            _, triangular = np.linalg.qr(vectors)
            reduced = triangular @ np.diag([1-float(p), float(p), -1]) @ triangular.conj().T
            distance += float(np.sum(np.abs(np.linalg.eigvalsh(reduced))))
        bound = np.pi*float(tau)+2*float(tau)/8+2*float(p)
        self.assertLess(distance, bound)
        budgets = []
        for i, m, eps in ((2, 1, Q(1, 100)), (20, 10, Q(1, 1000))):
            duration = eps/(4*m)
            self.assertLess(Q(22, 7)*m*duration, eps)
            budgets.append(dict(I=i, reads=m, target_error=str(eps), duration_per_read=str(duration),
                coupling_g_over_pi=str(1/duration), fresh_meter_qubits=m*(2*i+2),
                record_bits=m*(2*i+2), parity_xor_count=m*(2*i+1)))
        OBS['resources_and_noise'] = dict(noisy_full_DG_choi_trace_distance=old.short(distance),
            complete_noise_diamond_upper=old.short(bound), global_cat_trace_error=str(2*p),
            maximum_local_coupling_error='1/8', budgets=budgets,
            source_preparation_and_timing_still_inputs=True, energy_dissipation_lower_bound_claimed=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=507, scientific_baseline_round=506,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py', 'joint_orbit_port_interface.py',
            'role_covariant_graph_obstruction.py', 'branching_tree_distance_audit.py', 'research_note_506.md')},
        scope=dict(original_h_continues_during_read=True, local_occupation_reader_model_explicit=True,
            exact_joint_code_preservation_at_finite_duration=True, size_uniform_full_instrument_error=True,
            arbitrary_reference_in_single_excitation_sector=True, finite_history_error_control=True,
            global_cat_source_derived=False, autonomous_timing_derived=False,
            pure_original_swap_only_implementation=False, exact_ideal_read_claimed_at_finite_duration=False,
            dimension_three_generated=False, full_GR_goal_completed=False, phase_closure_triggered=False),
        observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(data, ensure_ascii=False, indent=2))
