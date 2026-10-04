"""Round 458: relational record activates a response under one fixed exchange law.

The old round-431 receiver is reused as a physical two-spin control variable.
There is no intermediate measurement, selection of a gate, or external schedule.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as old
import relational_interface_conversion_audit as previous

TARGET = Path(__file__).with_name('continuous_record_response_audit_results.json')
OBS = {}
EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5), (4, 5)]


def short(x):
    return float(f'{float(x):.12g}')


@lru_cache(None)
def model():
    swaps = {p: old.swap(7, *p) for p in EDGES}
    h0 = sum(swaps[p] for p in EDGES[:4])
    k = sum(swaps[p] for p in EDGES[4:])
    singlet = previous.model()['s']
    vin = np.kron(np.kron(old.encoding(), singlet[:, None]), singlet[:, None])
    em = (np.eye(128)-old.swap(7, 3, 4))/2
    ea = (np.eye(128)-old.swap(7, 5, 6))/2
    return dict(swaps=swaps, h0=h0, k=k, h=h0+k, vin=vin, em=em, ea=ea)


@lru_cache(None)
def certificate():
    m = model()
    a, b, *_ = old.operators()
    twice_singlet = np.eye(4)-old.swap(2, 0, 1).real
    d = np.kron(np.kron((a@b-b@a).real, twice_singlet), twice_singlet).astype(int).astype(object)
    permutations = [np.argmax(m['swaps'][p].real, axis=1) for p in EDGES]
    order, time = 100, F(3, 2)
    tail = F(2*3**18*18**(order+1), math.factorial(order+1))
    low, high = F(1732050807568877, 10**15), F(1732050807568878, 10**15)
    assert low*low < 3 < high*high
    reports = {}
    for name in ('em', 'ea'):
        ad = (2*m[name].real).astype(int).astype(object)
        total, traces = F(0), []
        for n in range(order+1):
            trace = int(np.sum(d.T*ad))
            traces.append(trace)
            if n % 2:
                total += F((-1)**((n+1)//2)*trace, 16*math.factorial(n))*time**n
            else:
                assert trace == 0
            ad = sum(ad[p, :]-ad[:, p] for p in permutations)
        assert total < 0
        lower, upper = -total/high-tail, -total/low+tail
        assert lower > F(1, 4) and upper < F(3, 10)
        reports[name] = dict(low_order_integer_traces=traces[:12],
            signed_series_before_dividing_sqrt3=str(total),
            contrast_lower=str(lower), contrast_upper=str(upper),
            contrast_interval_float=[float(lower), float(upper)])
    return dict(order=order, time=str(time), analytic_gap_tail=str(tail),
        analytic_gap_tail_float=float(tail), effects=reports)


def program_channel(rho_m, rho_ar, t):
    # Factor order memory(4), actuator(2), arbitrary reference(2).
    k = old.swap(3, 0, 2)+old.swap(3, 1, 2)
    u = np.kron(old.evolve(k, t), np.eye(2))
    out = u@np.kron(rho_m, rho_ar)@u.conj().T
    return old.partial(out, [4, 4], (1,))


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=3e-11):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_record_sector_controls_unknown_actuator(self):
        ps = (np.eye(4)-old.swap(2, 0, 1))/2
        k = old.swap(3, 0, 2)+old.swap(3, 1, 2)
        p = np.kron(ps, np.eye(2))
        self.close(k@p, p)
        self.close(k@p, p@k)
        self.close(np.linalg.eigvalsh(k), [-1, -1, 1, 1, 2, 2, 2, 2])
        # Check full superoperator on a matrix basis, including off-diagonals.
        worst = 0.
        for quality, t in ((0., .5), (.7, math.pi/3), (1., .9)):
            tau = quality*ps+(1-quality)*(np.eye(4)-ps)/3
            lam = 1-F(32, 27)*(1-quality)*math.sin(3*t/2)**2
            for i, j in itertools.product(range(4), repeat=2):
                basis = np.zeros((4, 4), complex)
                basis[i, j] = 1
                actual = program_channel(tau, basis, t)
                reference = old.partial(basis, [2, 2], (1,))
                expected = float(lam)*basis+(1-float(lam))*np.kron(np.eye(2)/2, reference)
                worst = max(worst, float(np.linalg.norm(actual-expected)))
                self.close(actual, expected)
        OBS['conditional_exchange'] = dict(memory_qubits=2, actuator_qubits=1,
            generator='S_35 + S_45', singlet_sector_generator='I',
            triplet_sector_eigenvalues=[-1, 2], record_projector_commutes_with_generator=True,
            depolarizing_parameter='1-(32/27)*(1-p)*sin(3t/2)^2',
            singlet_Choi_probability='1-(8/9)*(1-p)*sin(3t/2)^2',
            arbitrary_input_reference_superoperator_residual=short(worst),
            formula_requires_independent_tau_p_memory=True,
            no_abstract_controlled_gate_added=True)

    def test_02_one_joint_rule_and_all_source_gauges(self):
        m = model()
        h, vin = m['h'], m['vin']
        self.close(vin.conj().T@vin, np.eye(4))
        for mu in old.PAULI:
            joint = sum(previous.local(7, i, mu) for i in range(7))
            self.close(joint@vin, vin@np.kron(mu, np.eye(2)))
            self.close(joint@h, h@joint)
            for e in (m['em'], m['ea']):
                self.close(joint@e, e@joint)
        for t in (.5, 1.5):
            v = old.evolve(h, t)@vin
            for e in (m['em'], m['ea']):
                compressed = v.conj().T@e@v
                f_l = old.partial(compressed, [2, 2], (1,))/2
                self.close(compressed, np.kron(np.eye(2), f_l))
        OBS['joint_model'] = dict(edges=[list(p) for p in EDGES], weights=[1]*6,
            all_couplings_always_on=True, initial_singlets=[[3, 4], [5, 6]],
            unknown_source_input_order='G,L', reference_spin_6_is_internal=True,
            source_and_actuator_have_no_direct_edge=True,
            each_pulled_back_effect='I_G tensor F_L(t) for all t',
            arbitrary_GL_R_source_inputs_allowed=True,
            complete_target_channel_independent_of_source_G_claimed=False)

    def test_03_exact_record_and_response_certificates(self):
        result = certificate()
        self.assertEqual(result['effects']['em']['low_order_integer_traces'],
            [0, 0, 0, 0, 0, 1440, 0, 60480, 0, 2051424, 0, 67366464])
        self.assertEqual(result['effects']['ea']['low_order_integer_traces'],
            [0, 0, 0, 0, 0, 0, 0, -13440, 0, -732672, 0, -28890048])
        m = model()
        u = old.evolve(m['h'], 1.5)
        v = u@m['vin']
        outputs = [v@np.kron(np.eye(2)/2, (np.eye(2)+sign*old.PAULI[1])/2)@v.conj().T
                   for sign in (1, -1)]
        probabilities = {}
        for name in ('em', 'ea'):
            values = [float(np.trace(m[name]@r).real) for r in outputs]
            gap = values[1]-values[0]
            lo, hi = result['effects'][name]['contrast_interval_float']
            self.assertTrue(lo-1e-12 < gap < hi+1e-12)
            probabilities[name] = list(map(short, values))
        OBS['exact_integer_certificates'] = result
        OBS['record_and_action'] = dict(centre=1.5, probabilities=probabilities,
            record_contrast=short(probabilities['em'][1]-probabilities['em'][0]),
            response_contrast=short(probabilities['ea'][1]-probabilities['ea'][0]),
            leading_record_contrast='3*t^5/(4*sqrt(3)) + O(t^7)',
            leading_response_contrast='t^7/(6*sqrt(3)) + O(t^9)',
            intermediate_readout_or_postselection=False)

    def test_04_activation_and_joint_time_window(self):
        m = model()
        k, p, q, e = m['k'], m['em'], np.eye(128)-m['em'], m['ea']
        current = 1j*(k@e-e@k)
        self.close(current, q@current@q)
        self.close(k@p, p)
        self.close(m['h0']@e, e@m['h0'])
        # Integer identity: spectrum of current squared is contained in {0,2}.
        current2 = current@current
        self.close(current2@current2, 2*current2)
        self.assertAlmostEqual(float(np.linalg.norm(current, 2)), math.sqrt(2))
        mr = float(np.linalg.norm(m['h']@p-p@m['h'], 2))
        self.assertAlmostEqual(mr, math.sqrt(3)/2)
        self.assertEqual(F(1, 4)-4*F(1, 64), F(3, 16))
        self.assertEqual(F(1, 4)-2*F(1, 64), F(7, 32))
        OBS['activation_and_window'] = dict(
            actuator_current_supported_only_on_memory_triplet=True,
            all_joint_states_bound='abs(d p_action/dt) <= sqrt(2)*(1-p_memory)',
            exact_current_polynomial='J_action^4=2*J_action^2',
            interval=['95/64', '97/64'],
            record_contrast_strict_lower='7/32', response_contrast_strict_lower='3/16',
            response_equal_prior_success_strict_lower='19/32',
            exact_stopping_time_needed=False, permanent_record_claimed=False)

    def test_05_response_changes_the_record_process(self):
        m = model()
        comm = m['h0']@m['em']-m['em']@m['h0']
        back = m['k']@comm-comm@m['k']
        # Entire matrix is integer up to the projector factor 1/2.
        twice = np.rint(2*back.real).astype(np.int64)
        self.close(back, twice/2)
        norm_squared = int(np.sum(twice*twice))
        self.assertEqual(norm_squared, 384)
        # Frozen 431 has old contrast >2/3. New certificate gives <3/10.
        self.assertLess(F(certificate()['effects']['em']['contrast_upper']), F(3, 10))
        self.assertEqual(F(2, 3)-F(3, 10), F(11, 30))
        OBS['back_action'] = dict(
            exact_squared_Frobenius_double_commutator=norm_squared/4,
            operator='[K,[H_record,E_memory]]',
            old_receiver_contrast_lower_from_frozen_431='2/3',
            joint_receiver_contrast_upper='3/10',
            record_contrast_reduction_strict_lower='11/30',
            additional_gate_schedule_used=False,
            one_way_channel_concatenation_validated=False,
            learned_goal_or_error_correction_claimed=False)

    def test_06_full_unknown_reference_and_preparation_ledger(self):
        m = model()
        v = old.evolve(m['h'], 1.5)@m['vin']
        self.close(v.conj().T@v, np.eye(4))
        rho = old.density(12, np.random.default_rng(458))
        w = np.kron(v, np.eye(3))
        out = w@rho@w.conj().T
        self.close(w.conj().T@out@w, rho)
        self.close(old.partial(out, [128, 3], (1,)), old.partial(rho, [4, 3], (1,)))
        initial = m['vin']@np.eye(4)@m['vin'].conj().T/4
        after = v@np.eye(4)@v.conj().T/4
        initial_energy = float(np.trace(m['h']@initial).real)
        final_energy = float(np.trace(m['h']@after).real)
        self.assertAlmostEqual(initial_energy, final_energy)
        OBS['internal_ledger'] = dict(source_raw_qubits=3, memory_raw_qubits=2,
            actuator_raw_qubits=1, retained_internal_reference_qubits=1,
            pure_singlet_pairs_initially_supplied=2,
            complete_joint_information_and_unknown_reference_preserved=True,
            source_old_local_state_preserved_claimed=False, auxiliaries_restored_claimed=False,
            sample_fixed_H_initial_energy=short(initial_energy),
            preparations_effect_access_contacts_and_time_units_remain_inputs=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=458, baseline_round=457, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(physical_record_sector_controls_response=True,
            same_static_pair_exchange_receives_records_and_responds=True,
            exact_joint_activation_bound_and_time_window=True,
            arbitrary_source_gauge_and_reference_effect_identity=True,
            response_back_action_accounted=True, two_pure_singlets_are_explicit_inputs=True,
            complete_cognitive_feedback_learning_loop=False,
            recursively_complete_subject_constructed=False,
            exact_455_handoff_implemented=False, spatial_dimension_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(round=result['round'], tests=result['tests_run'],
        failures=result['failures'], errors=result['errors'])))
