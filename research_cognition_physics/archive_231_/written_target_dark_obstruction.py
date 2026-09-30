"""Round 515: an actually written target need not be reached by target-only monitoring.

The unchanged 510 Hamiltonian and the exact 514 source are used throughout.
Integer invariant-subspace identities and rational Taylor bounds prove the claim;
finite numerical histories are diagnostics, not a proof of infinite-time failure.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import historical_membership_interface as old
import written_record_persistence as formation
from distributed_role_reader import exponential

HERE = Path(__file__).resolve().parent
TARGET = HERE/'written_target_dark_obstruction_results.json'
OBS = {}
T = Q(1, 2048)


def coefficients():
    state = formation.initial(2, np.int64)[:, :, 0:1]
    out = [state[1, 12:18, 0].copy()]
    for _ in range(4):
        state = old.apply_historical(state)
        out.append(state[1, 12:18, 0].copy())
    return out


def rational_bounds():
    r = formation.tail(11*T, 5)/T**3
    low = Q(1, 6)-2*r
    up = Q(1, 6)+Q(5, 12)*T+r
    one = Q(1, 6)-r
    return dict(tail_over_T3=r, difference_lower_over_T3=low,
                norm_upper_over_T3=up, source_probability_lower=one**2*T**6,
                source_probability_upper=up**2*T**6,
                source_and_never_detected_lower=low**2*T**6/8,
                conditional_never_detected_lower=low**2/(8*up**2),
                initial_target_adjacency_lower=one**2/up**2)


def actual_source():
    assert formation.FORMATION == T
    k = formation.selected(2, formation.source(2))[0][:, 0]/float(T**3)
    probability_scaled = float(np.vdot(k, k).real)
    psi = np.zeros(36, complex)
    psi[12:18] = k/math.sqrt(probability_scaled)
    return k, psi, probability_scaled


def invariant_integer_basis():
    # w=2z, with graph indices exactly those returned by old.model().
    w = np.zeros(36, dtype=np.int64)
    w[[12, 17, 18, 23]] = [1, -1, -1, 1]
    representatives = [m for m in range(64) if m < (m ^ 12)]
    index = {m: j for j, m in enumerate(representatives)}
    basis = np.zeros((64, 36, 32), dtype=np.int64)
    x2 = np.zeros((32, 32), dtype=np.int64)
    for j, m in enumerate(representatives):
        basis[m, :, j] = w
        basis[m ^ 12, :, j] = w
        flipped = m ^ 4
        x2[index[min(flipped, flipped ^ 12)], j] = 1
    return w, basis, x2


class Audit(unittest.TestCase):
    def test_01_exact_same_H_dark_subspace(self):
        a = old.model()
        w, basis, x2 = invariant_integer_basis()
        self.assertTrue(np.array_equal(a['h']@w, -w))
        self.assertEqual(int(w@w), 4)
        self.assertFalse(np.any(basis[:, :6, :]))
        flat = basis.reshape(2304, 32)
        self.assertTrue(np.array_equal(flat.T@flat, 8*np.eye(32, dtype=int)))
        for coupling in (-2, 0, 1, 3):
            lhs = old.apply_historical(basis, coupling)
            rhs = basis@(-np.eye(32, dtype=int)+coupling*x2)
            self.assertTrue(np.array_equal(lhs, rhs))
        self.assertTrue(np.array_equal(x2@x2, np.eye(32, dtype=int)))
        self.assertEqual(int(abs(a['h']).sum(axis=1).max())+1, 11)
        OBS['integer_dark_certificate'] = dict(
            graph_edges=[sorted(map(list, tree)) for tree in a['trees']],
            unnormalized_DG_vector=w.tolist(), vector_norm_squared=4,
            memory_pair_flip_mask=12, witness_dimension=32,
            integer_basis_gram=8, exact_restricted_H='-I + g X_2',
            target=0, target_projection_annihilates_basis=True,
            complete_dark_space_classification_claimed=False,
            H_abs_row_bound_at_g1=11)

    def test_02_exact_actual_source_coefficients(self):
        cs = coefficients()
        expected = [[0]*6]*3+[[1, 0, 0, 0, 0, 0], [-8, 4, 4, 1, 1, 0]]
        self.assertEqual([v.tolist() for v in cs], expected)
        self.assertEqual(int(cs[4]@cs[4]), 98)
        self.assertIn((0, 2), old.model()['trees'][0])
        OBS['actual_source_coefficients'] = dict(anchor=2, written_target=0,
            input_graph=0, memory_initial_Z_string=0, selected_memory_Z_string=1,
            selected_integer_H_powers=[v.tolist() for v in cs],
            fourth_coefficient_norm_squared=98, first_nonzero_order=3)

    def test_03_rational_source_and_conditional_obstruction(self):
        b = rational_bounds()
        low, up = b['difference_lower_over_T3'], b['norm_upper_over_T3']
        self.assertGreater(low, 0)
        self.assertGreater(9*low**2, 8*up**2)
        self.assertGreater(b['source_probability_lower'], T**6/49)
        self.assertGreater(b['conditional_never_detected_lower'], Q(1, 9))
        self.assertGreater(b['initial_target_adjacency_lower'], Q(9, 10))
        OBS['rational_certificate'] = dict(duration=str(T), coupling=1,
            **{key: str(value) for key, value in b.items()},
            simplified_conditional_timeout_lower='1/9',
            all_waiting_schedules_and_all_finite_detection_counts=True)

    def test_04_full_actual_source_diagnostic(self):
        k, psi, p_scaled = actual_source()
        b = rational_bounds()
        state = formation.source(2)[:, :, 0]
        completeness = abs(float(np.vdot(state, state).real)-1)
        dark = abs(k[0]-k[5])**2/(8*p_scaled)
        self.assertLess(completeness, 2e-13)
        self.assertGreater(dark, float(b['conditional_never_detected_lower']))
        self.assertGreater(p_scaled, float(b['source_probability_lower']/T**6))
        self.assertLess(p_scaled, float(b['source_probability_upper']/T**6))
        self.assertAlmostEqual(float(np.vdot(psi, psi).real), 1, places=13)
        OBS['source_diagnostic'] = dict(
            branch_probability_divided_by_T6=old.old.short(p_scaled),
            conditional_dark_witness_weight=old.old.short(float(dark)),
            complete_source_instrument_error=old.old.short(completeness),
            rare_branch_not_assumed_free=True, prior_graph_not_reset=True,
            selected_graph_vector_over_T3_real=[old.old.short(float(v.real)) for v in k],
            selected_graph_vector_over_T3_imag=[old.old.short(float(v.imag)) for v in k])

    def test_05_complete_monitored_histories_same_H(self):
        a = old.model()
        k, psi, p_scaled = actual_source()
        dark = float(abs(k[0]-k[5])**2/(8*p_scaled))
        cases = []
        for name, intervals in [('uniform', [Q(1, 8)]*32),
                                ('nonuniform', [Q(1, 16), Q(1, 4), Q(1, 8)]*8)]:
            # The actual source is M=e_0, including its relative X-sector phases.
            waiting = a['hadamard'][:, 1, None]*np.tile(psi, (64, 1))/8
            unitaries = {t: np.array([exponential(float(t)*(a['h']+
                np.diag(np.repeat(s, 6)))) for s in a['signs']]) for t in set(intervals)}
            arrivals = []
            survival = []
            for t in intervals:
                advanced = np.einsum('sij,sj->si', unitaries[t], waiting)
                arrivals.append(float(np.vdot(advanced[:, :6], advanced[:, :6]).real))
                advanced[:, :6] = 0
                waiting = advanced
                mass = float(np.vdot(waiting, waiting).real)
                survival.append(mass)
                self.assertGreaterEqual(mass+2e-12, dark)
            completeness = abs(sum(arrivals)+survival[-1]-1)
            self.assertLess(completeness, 3e-12)
            cases.append(dict(schedule=name, intervals=[str(t) for t in intervals],
                all_first_arrival_probabilities=[old.old.short(v) for v in arrivals],
                timeout_probabilities=[old.old.short(v) for v in survival],
                full_instrument_probability_error=old.old.short(completeness)))
        OBS['monitored_histories'] = dict(cases=cases,
            all_64_memory_X_sectors_retained=True,
            no_memory_monitoring_or_source_reset=True,
            infinite_time_claim_uses_integer_certificate_not_numerical_plateau=True)

    def test_06_nonzero_reception_and_finite_instrument_error_scope(self):
        b = rational_bounds()
        wait = Q(1, 1024)
        remainder = formation.tail(11*wait, 2)
        self.assertLess(remainder, wait/4)
        self.assertGreater(b['initial_target_adjacency_lower'], Q(9, 10))
        # sqrt(9/10)>3/4; leading arrival >=3u/4 and remainder <=u/4.
        lower = wait**2/4
        a = old.model()
        _, psi, _ = actual_source()
        arrival = sum(float(np.vdot((exponential(float(wait)*(a['h']+
            np.diag(np.repeat(s, 6))))@psi)[:6],
            (exponential(float(wait)*(a['h']+np.diag(np.repeat(s, 6))))@psi)[:6]).real)
            for s in a['signs'])/64
        self.assertGreater(arrival, float(lower))
        # This eta is a hypothetical complete-protocol error, not achieved hardware.
        eta = b['source_and_never_detected_lower']/10
        conditional = (b['source_and_never_detected_lower']-eta/2)/(
            b['source_probability_upper']+eta/2)
        self.assertGreater(conditional, Q(1, 10))
        OBS['scope_and_error'] = dict(first_wait=str(wait),
            first_detection_probability=old.old.short(arrival),
            first_detection_rational_lower=str(lower),
            nonzero_reception_compatible_with_permanent_undetected_component=True,
            hypothetical_complete_protocol_diamond_eta=str(eta),
            perturbed_conditional_lower=str(conditional),
            implemented_precision_or_infinite_error_accumulation_not_claimed=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(date='2026-09-28', round=515, scientific_base_through_round=514,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('historical_membership_interface.py','written_record_persistence.py',
                         'distributed_role_reader.py','research_note_510.md',
                         'research_note_512.md','research_note_514.md')},
        scope=dict(same_H_as_actual_record_formation=True,
            actual_514_source_with_strictly_positive_probability=True,
            target_only_arbitrary_schedule_obstruction_proved=True,
            all_sources_impossible=False, arbitrary_active_routing_excluded=False,
            one_failed_spatial_model_refutes_cognitive_axioms=False,
            autonomous_readout_or_timing_generated=False,
            physical_position_or_dimension_generated=False,
            full_GR_goal_completed=False), observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('round','tests_run','failures','errors')},
                     ensure_ascii=False))
