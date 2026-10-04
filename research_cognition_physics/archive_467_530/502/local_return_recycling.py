"""Round 502: local monitored return restores a reusable query source.

Uses known finite-dimensional monitored subspace recurrence. Only fresh
one-step replies are accepted; later replies are recorded and recycled locally.
"""
import argparse
from fractions import Fraction as F
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import local_neighbor_receipt as frozen501
import current_neighbor_detection as frozen500

HERE = Path(__file__).resolve().parent
TARGET = HERE/'local_return_recycling_results.json'
OBS = {}


def cycle_blocks(u, source, graph_count, colours, cap):
    """Root is port zero. Return time and colour are distinct retained records."""
    away = np.ones(len(u), bool)
    away[:graph_count*colours] = False
    survivor = np.eye(len(u), dtype=complex)[:, source]
    fresh, reset = [], []
    for step in range(1, cap+1):
        next_state = u@survivor
        for colour in range(colours):
            block = next_state[np.arange(graph_count)*colours+colour]
            (fresh if step == 1 and colour != 0 else reset).append((step, colour, block))
        survivor = next_state*away[:, None]
    return fresh, reset, survivor


def reference_after(k, rho):
    """rho indexed graph, reference, graph, reference; K may be rectangular."""
    return np.einsum('ag,grhs,ah->rs', k, rho, k.conj())


def outward(value, denominator, upper=False):
    x = value*denominator
    integer = -((-x.numerator)//x.denominator) if upper else x.numerator//x.denominator
    return F(integer, denominator)


class Checks(unittest.TestCase):
    def test_01_exact_recurrence_identity(self):
        d = 5
        u = np.roll(np.eye(d, dtype=np.int64), 1, axis=0)
        p = np.diag([1, 0, 1, 0, 0])
        ell = (np.eye(d, dtype=np.int64)-p)@u
        self.assertTrue(np.array_equal(ell@ell.T, np.eye(d, dtype=np.int64)-p))
        power = np.eye(d, dtype=np.int64)
        total = np.zeros((d, d), dtype=np.int64)
        expected = total.copy()
        survival = []
        rho = np.diag([1, 0, 0, 0, 0])
        for k in range(7):
            total += power@p@power.T
            expected += p@power.T@power@p
            survival.append(int(np.trace(power@rho@power.T)))
            power = ell@power
            self.assertTrue(np.array_equal(total, np.eye(d, dtype=np.int64)-power@power.T))
        self.assertEqual(sum(survival), 2)
        self.assertEqual(int(np.trace(expected)), 5)
        self.assertLessEqual(np.linalg.eigvalsh(expected)[-1], d)
        OBS['exact_recurrence'] = dict(dimension=d, observed_subspace_rank=2,
            partial_sum_trace=int(np.trace(total)), test_state_mean_return=2,
            exact_telescoping=True, unknown_reference_not_in_dimension_bound=True)

    def test_02_actual_two_cycle_instrument(self):
        trees, _, (h, labels, source, _, _) = frozen501.six()
        u = frozen500.evolution(h, F(1, 64))
        fresh, reset, timeout = cycle_blocks(u, source, 6, 6, 8)
        complete = timeout.conj().T@timeout
        for _, _, k in fresh+reset:
            complete += k.conj().T@k
        self.assertLess(np.linalg.norm(complete-np.eye(6)), 1e-12)
        psi = np.array([[1, 1j], [2j, -1], [1+1j, 2], [-1j, 1], [2, -2j], [1, 3j]], complex)
        psi /= np.linalg.norm(psi)
        initial = np.outer(psi.ravel(), psi.ravel().conj())
        active = initial.copy()
        reference = np.zeros((2, 2), complex)
        accepted = bad_mass = timeout_mass = 0.0
        for cycle in range(2):
            rho4 = active.reshape(6, 2, 6, 2)
            for _, colour, k in fresh:
                kr = np.kron(k, np.eye(2))
                out = kr@active@kr.conj().T
                accepted += float(np.trace(out).real)
                bad = np.repeat([(0, labels[colour]) not in tree for tree in trees], 2)
                bad_mass += float(np.trace(out[np.ix_(bad, bad)]).real)
                reference += reference_after(k, rho4)
            ref_timeout = reference_after(timeout, rho4)
            timeout_mass += float(np.trace(ref_timeout).real)
            reference += ref_timeout
            following = np.zeros_like(active)
            for _, _, k in reset:
                kr = np.kron(k, np.eye(2))
                following += kr@active@kr.conj().T
            active = following
        unfinished = float(np.trace(active).real)
        reference += np.einsum('grgs->rs', active.reshape(6, 2, 6, 2))
        self.assertAlmostEqual(accepted+timeout_mass+unfinished, 1, places=12)
        self.assertLess(np.linalg.norm(reference-psi.T@psi.conj()), 2e-12)
        self.assertLess(bad_mass/accepted, 1/10000)
        OBS['actual_cycle_diagnostics'] = dict(active_dimension=len(h), cycles=2,
            maximum_measurements_per_cycle=8, natural_step='1/64',
            fresh_success_branches=len(fresh), local_reset_branches=len(reset),
            accepted=accepted, timeout=timeout_mass, final_cycle_exhaustion=unfinished,
            current_joint_error_given_acceptance=bad_mass/accepted,
            one_cycle_completeness_residual=float(np.linalg.norm(complete-np.eye(6))),
            full_reference_residual=float(np.linalg.norm(reference-psi.T@psi.conj())),
            graph_copies=1, packet_copies=1, high_success_protocol_actually_run=False)

    def test_03_local_reset_retains_colour(self):
        colours = 6
        # |a> -> |q>_packet |a>_record: the colour is transferred, not erased.
        w = np.zeros((colours*colours, colours), dtype=np.int64)
        for a in range(colours):
            w[a, a] = 1
        self.assertTrue(np.array_equal(w.T@w, np.eye(colours, dtype=np.int64)))
        psi = np.arange(1, 13).reshape(6, 2).astype(complex)
        psi[:, 1] *= 1j
        psi /= np.linalg.norm(psi)
        out = (w@psi).reshape(colours, colours, 2)
        self.assertFalse(np.any(out[1:]))
        self.assertLess(np.linalg.norm(out[0]-psi), 1e-15)
        OBS['local_reset'] = dict(colours=colours, isometry_exact=True,
            packet_restored_to_query=True, colour_reference_state_retained_in_record=True,
            root_only_feedback_is_extra_input=True, no_unrecorded_erasure=True)

    def test_04_late_reply_is_not_current_certificate(self):
        # Static three-vertex path: a late label 2 is certainly not adjacent to 0.
        trees = [frozenset({(0, 1), (1, 2)})]
        h, _, source, reply, bad = frozen501.model(trees, np.zeros((1, 1), dtype=np.int64), 3)
        real, imag, scale, tail, norm = frozen501.source_taylor(h, np.arange(9), 8, 12)
        vr = np.eye(9, dtype=object)[:, source]
        vi = np.zeros_like(vr)
        sum_reply = sum_bad = 0
        denominator = 1
        cap = 64
        away = np.arange(9) >= 3
        for step in range(1, cap+1):
            zr, zi = real@vr-imag@vi, real@vi+imag@vr
            sum_reply *= scale*scale
            sum_bad *= scale*scale
            if step > 1:
                sum_reply += sum(int(x*x+y*y) for x,y in zip(zr[reply].ravel(), zi[reply].ravel()))
                sum_bad += sum(int(x*x+y*y) for x,y in zip(zr[bad].ravel(), zi[bad].ravel()))
            denominator *= scale*scale
            vr, vi = zr*away[:, None], zi*away[:, None]
        # Each product has <=cap contractions; truncations have norm <=1+tail.
        e = cap*tail/(1-cap*tail)
        total_error = cap*(2*e+e*e)
        lower_bad = outward(F(sum_bad, denominator)-total_error, 10**9)
        upper_reply = outward(F(sum_reply, denominator)+total_error, 10**9, upper=True)
        self.assertGreater(lower_bad, upper_reply/4)
        OBS['late_reply_counterexample'] = dict(N=3, graph='0-1-2', kappa=0,
            step_time='1/8', last_return_step=cap, accepted_steps_in_counterexample='2..64',
            norm_bound=norm, Taylor_degree=12, exact_unitary_tail=str(tail),
            late_nonedge_mass_lower=str(lower_bad), late_reply_mass_upper=str(upper_reply),
            conditional_nonedge_error_strict_lower='1/4',
            proof_arithmetic='integer repeated products, rational tail and outward rounding',
            actual_protocol_accepts_these_late_replies=False)

    def test_05_fixed_step_scope(self):
        # H=X: first pi/2 step is X up to phase, later pi steps are I up to phase.
        p = np.diag([1, 0]); q = np.eye(2, dtype=np.int64)-p
        first = np.array([[0, 1], [1, 0]], dtype=np.int64)
        later = np.eye(2, dtype=np.int64)
        state = np.array([1, 0], dtype=np.int64)
        state = q@first@state
        self.assertTrue(np.array_equal(state, [0, 1]))
        self.assertTrue(np.array_equal(q@later@state, state))
        OBS['fixed_step_scope'] = dict(same_H='Pauli X', first_interval='pi/2',
            subsequent_interval='pi', nonreturn_probability_exact=1,
            arbitrary_changed_sampling_schedule_supported=False,
            theorem_requires_initial_state_in_entire_observed_root_subspace=True)

    def test_06_finite_high_success_budget(self):
        s0, eps, dim = F(1, 10**13), F(1, 10000), 216
        m = int(4/s0)
        cap = 50*m*dim
        self.assertGreater(sum(F(4**k, math.factorial(k)) for k in range(9)), 50)
        self.assertLess(F(m*dim, cap+1), F(1, 50))
        failure = F(1, 25)
        total_delta = F(1, 10000)
        delay = F(1, 4096)
        delayed_error = (F(1, 100)+24*delay)**2
        actual_error = (delayed_error+total_delta)/(1-failure-total_delta)
        self.assertLess(actual_error, F(1, 1000))
        # Each of m cycles uses <=cap complete measurement/control slots.
        slots = m*cap
        each_error = total_delta/(2*slots)
        source_error = total_delta/2
        self.assertEqual(source_error+slots*each_error, total_delta)
        OBS['finite_resource_certificate'] = dict(s0=str(s0), fresh_error_upper=str(eps),
            active_dimension=dim, maximum_cycles=m, measurements_per_cycle=cap,
            complete_slots=slots, ideal_failure_strict_upper='1/25',
            delay=str(delay), actual_success_strict_lower=str(1-failure-total_delta),
            actual_joint_error_upper=str(actual_error), actual_joint_error_strict_upper='1/1000',
            source_half_trace_error=str(source_error),
            per_complete_slot_half_diamond_error=str(each_error),
            root_raw_record_bits=3*slots, natural_waiting_sum_upper=str(F(slots, 64)),
            graph_copies=1, packet_register_sets=1, controller_clock_isolation_extra=True,
            massive_protocol_executed=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    dep = HERE/'local_neighbor_receipt_results.json'
    return dict(round=502, scientific_baseline_round=501,
        reused_frozen_rounds=[493, 500, 501], tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_results_sha256={dep.name: hashlib.sha256(dep.read_bytes()).hexdigest()},
        primary_source=dict(url='https://arxiv.org/abs/1302.7286', version='1302.7286v1',
            consulted='Sections 2 and 3, Theorem 3.2, Appendix A; finite monitored subspace recurrence',
            tool_theory_claimed_original=False),
        scope=dict(same_unknown_graph_and_packet_recycled=True, root_only_instrument=True,
            first_step_fresh_replies_only=True, late_replies_recorded_and_locally_reset=True,
            finite_fixed_U_and_known_dimension_bound_required=True,
            unknown_reference_dimension_not_charged_to_return_bound=True,
            new_local_feedback_and_finer_failure_readout_are_inputs=True,
            autonomous_complete_controller_generated=False, natural_H_suspended=False,
            common_deadline_graph_freeze=False, three_dimensional_positions_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False), observations=OBS.copy())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
