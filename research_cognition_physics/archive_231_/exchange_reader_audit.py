"""Round 431: a two-qubit internal reader under continuous exchange only.

Uses the frozen round 430 operators. An exact integer commutator series and
rational remainder certify transfer to an internal pair, without a gate schedule.
"""
import argparse
from fractions import Fraction
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('exchange_reader_audit_results.json')
OBS = {}


def setup():
    a, b, c, p, x, y, z = core.operators()
    meter = (np.eye(4)-core.swap(2, 0, 1))/2
    plus, minus = np.kron((p+y)/4, meter), np.kron((p-y)/4, meter)
    exchanges = [core.swap(5, j, j+1) for j in range(4)]
    h = sum(exchanges)
    effect = np.kron(np.eye(8), meter)
    return h, effect, plus, minus, exchanges


def exact_certificate(order=80):
    a, b, *_ = core.operators()
    difference_scaled = np.kron((a@b-b@a).real, np.eye(4)-core.swap(2, 0, 1).real).astype(object)
    difference_scaled = np.vectorize(int)(difference_scaled).astype(object)
    effect_scaled = np.kron(np.eye(8), np.eye(4)-core.swap(2, 0, 1).real)
    ad = np.vectorize(int)(effect_scaled).astype(object)
    permutations = [np.argmax(core.swap(5, j, j+1).real, axis=1) for j in range(4)]
    total = Fraction(0)
    traces = []
    for k in range(order+1):
        trace = int(np.sum(difference_scaled.T*ad))
        traces.append(trace)
        if k % 2:
            total += Fraction((-1)**((k+1)//2)*trace, 8*math.factorial(k))*Fraction(3, 2)**k
        else:
            assert trace == 0
        ad = sum(ad[perm, :]-ad[:, perm] for perm in permutations)
    # Delta probability = total / sqrt(3) plus this analytic remainder.
    # ||H|| <= 4, ||E|| = 1, ||rho_plus-rho_minus||_1 = 2, t = 3/2.
    tail = Fraction(2*3**12*12**(order+1), math.factorial(order+1))
    low_sqrt = Fraction(1732050807568877, 10**15)
    high_sqrt = Fraction(1732050807568878, 10**15)
    assert low_sqrt**2 < 3 < high_sqrt**2 and total < 0
    lower, upper = -total/high_sqrt-tail, -total/low_sqrt+tail
    assert lower > Fraction(2, 3) and upper < 1
    return dict(order=order, model_time='3/2', low_order_integer_traces=traces[:8],
        signed_series_before_dividing_sqrt3=str(total),
        analytic_tail_upper=str(tail), analytic_tail_upper_float=float(tail),
        contrast_lower_rational=str(lower), contrast_upper_rational=str(upper),
        contrast_interval_float=[float(lower), float(upper)],
        contrast_exceeds_two_thirds_by_exact_arithmetic=True), lower, upper


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=5e-12):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_common_rule_and_initial_reader(self):
        h, effect, plus, minus, exchanges = setup()
        self.close(h, h.conj().T)
        for s in exchanges:
            self.close(s@s, np.eye(32))
        self.assertAlmostEqual(float(np.linalg.norm(h, 2)), 4)
        self.close(effect@effect, effect)
        for rho in (plus, minus):
            self.close(np.trace(rho), 1)
            self.assertGreater(np.linalg.eigvalsh(rho).min(), -1e-12)
            self.assertAlmostEqual(float(np.trace(effect@rho).real), 1)
        self.close(core.partial(plus, [2]*5, (3, 4)), core.partial(minus, [2]*5, (3, 4)))
        OBS['setup'] = dict(subject_qubits=3, reader_qubits=2,
            constant_exchange_pairs=['12', '23', '34', '45'], coupling_strengths=[1, 1, 1, 1],
            reader_initial_state='singlet_45', all_couplings_permanently_on=True,
            coupling_pattern_and_singlet_preparation_are_inputs=True)

    def test_02_exact_commutator_and_interval_certificate(self):
        certificate, lower, upper = exact_certificate()
        self.assertEqual(certificate['low_order_integer_traces'], [0, 0, 0, 0, 0, 720, 0, 19488])
        self.assertGreater(lower, Fraction(2, 3))
        self.assertLess(upper, 1)
        OBS['exact_transfer_certificate'] = certificate

    def test_03_full_unitary_and_reduced_reader(self):
        h, effect, plus, minus, exchanges = setup()
        u = core.evolve(h, 1.5)
        outputs = [u@rho@u.conj().T for rho in (plus, minus)]
        readers = [core.partial(rho, [2]*5, (3, 4)) for rho in outputs]
        probabilities = [float(np.trace(effect@rho).real) for rho in outputs]
        gap = probabilities[1]-probabilities[0]
        self.assertGreater(gap, 2/3)
        certificate = OBS['exact_transfer_certificate']
        low, high = certificate['contrast_interval_float']
        self.assertGreater(gap, low-1e-13)
        self.assertLess(gap, high+1e-13)
        self.assertAlmostEqual(core.distance(*readers), gap)
        singlet = (np.eye(4)-core.swap(2, 0, 1))/2
        for probability, rho in zip(probabilities, readers):
            expected = probability*singlet+(1-probability)*(np.eye(4)-singlet)/3
            self.close(rho, expected)
        self.assertAlmostEqual(core.distance(*outputs), 1)
        OBS['reader_transfer'] = dict(time=1.5, singlet_probabilities=probabilities,
            reader_trace_distance=gap, complete_system_trace_distance=1.,
            equal_prior_single_read_success_probability=(1+gap)/2,
            no_postselection=True)

    def test_04_finite_continuous_read_window(self):
        h, effect, plus, minus, exchanges = setup()
        commutator = 1j*(h@effect-effect@h)
        norm = float(np.linalg.norm(commutator, 2))
        self.assertAlmostEqual(norm, math.sqrt(3)/2)
        rows = []
        for t in (17/12, 1.5, 19/12):
            u = core.evolve(h, t)
            gap = float(np.trace(effect@u@(minus-plus)@u.conj().T).real)
            self.assertGreater(gap, .5)
            rows.append(dict(time=t, contrast=gap))
        # Certified gap at centre >2/3; |gap'| <=sqrt(3)<2; |dt|<=1/12.
        self.assertEqual(Fraction(2, 3)-2*Fraction(1, 12), Fraction(1, 2))
        OBS['finite_read_window'] = dict(interval=['17/12', '19/12'],
            effect_derivative_norm=norm, contrast_strict_lower_bound='1/2',
            equal_prior_success_strict_lower_bound='3/4', checks_at_three_times=rows,
            exact_stop_or_exact_read_instant_required=False,
            permanent_stable_record_claimed=False)

    def test_05_unknown_reference_and_back_action(self):
        h, effect, plus, minus, exchanges = setup()
        rng = np.random.default_rng(43105)
        rho_dr = core.density(16, rng)
        meter = (np.eye(4)-core.swap(2, 0, 1))/2
        initial = np.kron(rho_dr, meter).reshape(8, 2, 4, 8, 2, 4).transpose(0, 2, 1, 3, 5, 4).reshape(64, 64)
        u = np.kron(core.evolve(h, 1.5), np.eye(2))
        after = u@initial@u.conj().T
        inverse_error = float(np.linalg.norm(u.conj().T@after@u-initial))
        self.close(u.conj().T@after@u, initial)
        self.close(core.partial(initial, [8, 4, 2], (2,)), core.partial(after, [8, 4, 2], (2,)))
        initial_d = core.partial(initial, [8, 4, 2], (0,))
        final_d = core.partial(after, [8, 4, 2], (0,))
        disturbance = core.distance(initial_d, final_d)
        self.assertGreater(disturbance, .01)
        OBS['reference_and_disturbance'] = dict(arbitrary_input_dimension=8,
            reference_dimension=2, inverse_recovery_error=inverse_error,
            sample_subject_trace_distance_change=disturbance,
            local_subject_state_preservation_claimed=False,
            full_joint_unknown_information_preserved=True)

    def test_06_internal_axis_independence_and_reader_minimum(self):
        h, effect, plus, minus, exchanges = setup()
        one = core.evolve(.31*core.PAULI[0]-.24*core.PAULI[1]+.17*core.PAULI[2], 1)
        common = one
        for _ in range(4):
            common = np.kron(common, one)
        self.close(common@h@common.conj().T, h)
        self.close(common@effect@common.conj().T, effect)
        for rho in (plus, minus):
            self.close(common@rho@common.conj().T, rho)
            for t in (.4, 1.5):
                u = core.evolve(h, t)
                after = u@rho@u.conj().T
                for j in range(5):
                    self.close(core.partial(after, [2]*5, (j,)), np.eye(2)/2)
        OBS['axis_and_reader_size'] = dict(no_absolute_spin_axis_in_state_rule_or_effect=True,
            every_single_qubit_marginal_remains_maximally_mixed=True,
            two_qubit_joint_reader_carries_distinguishable_relation=True,
            minimality_limited_to_collective_invariant_qubit_readers=True,
            common_identification_derived=False)


def run():
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=431, baseline_round=430, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(internal_pair_record_transfer_constructed=True,
            exact_integer_and_rational_signal_certificate=True,
            finite_read_window_proved=True, same_continuous_exchange_rule_only=True,
            unknown_full_joint_reference_information_preserved=True,
            singlet_source_and_coupling_pattern_are_inputs=True,
            exact_programmable_processor_required=False, permanent_record_proved=False,
            complete_readout_amplifier_constructed=False, physical_clock_derived=False,
            full_cognitive_implementation_completed=False,
            three_dimensional_space_unconditionally_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
