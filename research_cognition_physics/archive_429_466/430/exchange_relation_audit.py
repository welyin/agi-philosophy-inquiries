"""Round 430: relational dynamics of three qubits under fixed exchange.

The exchange-only encoding is established literature. This audit connects it
to round 429 and checks a complete static-pair/dynamic-pair distinction.
Only NumPy and the standard library are used; output creation is exclusive.
"""
import argparse
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

TARGET = Path(__file__).with_name('exchange_relation_audit_results.json')
OBS = {}
PAULI = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex), np.diag([1., -1.])]


def swap(n, first, second):
    out = np.zeros((2**n, 2**n), complex)
    for col in range(2**n):
        bits = [(col >> (n-1-j)) & 1 for j in range(n)]
        bits[first], bits[second] = bits[second], bits[first]
        row = sum(value << (n-1-j) for j, value in enumerate(bits))
        out[row, col] = 1
    return out


def evolve(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values)) @ vectors.conj().T


def partial(rho, dimensions, keep):
    keep = tuple(keep)
    drop = tuple(j for j in range(len(dimensions)) if j not in keep)
    order = keep+drop
    tensor = rho.reshape(tuple(dimensions)*2)
    tensor = tensor.transpose(order+tuple(j+len(dimensions) for j in order))
    left = math.prod(dimensions[j] for j in keep)
    right = math.prod(dimensions[j] for j in drop)
    return np.trace(tensor.reshape(left, right, left, right), axis1=1, axis2=3)


def density(d, rng):
    a = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
    out = a@a.conj().T
    return out/np.trace(out)


def distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)


def operators():
    a, b, c = swap(3, 0, 1), swap(3, 1, 2), swap(3, 0, 2)
    p = np.eye(8)-(a+b+c)/3
    x, z = (b-c)/math.sqrt(3), -p@a
    y = 1j*(a@b-b@a)/math.sqrt(3)
    return a, b, c, p, x, y, z


def encoding():
    _, _, _, _, x, _, _ = operators()
    singlet = np.array([0, 1, -1, 0], complex)/math.sqrt(2)
    first, second = np.kron(singlet, [1, 0]), np.kron(singlet, [0, 1])
    # Columns ordered as gauge0/logical0, gauge0/logical1, gauge1/logical0, ...
    return np.column_stack([first, x@first, second, x@second])


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=3e-12):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_minimum_noncommuting_relation(self):
        a, b, c, p, x, y, z = operators()
        basis = [np.eye(8)-p, p, x, y, z]
        self.assertEqual(np.linalg.matrix_rank(np.column_stack([q.ravel() for q in basis])), 5)
        s = swap(2, 0, 1)
        two = [np.eye(4), s]
        self.assertEqual(np.linalg.matrix_rank(np.column_stack([q.ravel() for q in two])), 2)
        self.close(a+b+c, 3*(np.eye(8)-p))
        self.close(a@b+b@a, a+b+c-np.eye(8))
        self.assertGreater(np.linalg.norm(a@b-b@a), 1)
        OBS['minimal_relation_algebra'] = dict(two_qubit_algebra='C direct_sum C',
            three_qubit_algebra='C direct_sum M_2', dimensions=[2, 5],
            minimality_uses_analytic_representation_decomposition=True)

    def test_02_intrinsic_pauli_algebra(self):
        a, b, c, p, x, y, z = operators()
        self.close(p@p, p)
        self.assertAlmostEqual(float(np.trace(p).real), 4)
        errors = []
        for q in (x, y, z):
            self.close(q@q, p)
            self.close(q, q.conj().T)
        for left, right, last in ((x, y, z), (y, z, x), (z, x, y)):
            errors.append(float(np.linalg.norm(left@right-1j*last)))
            self.close(left@right, 1j*last)
        v = encoding()
        self.close(v.conj().T@v, np.eye(4))
        self.close(v@v.conj().T, p)
        for q, logical in zip((x, y, z), PAULI):
            self.close(v.conj().T@q@v, np.kron(np.eye(2), logical))
        OBS['intrinsic_qubit'] = dict(support_rank=4, gauge_dimension=2,
            logical_dimension=2, maximum_pauli_error=max(errors))

    def test_03_identical_full_pair_marginals(self):
        a, b, c, p, x, y, z = operators()
        plus, minus = (p+y)/4, (p-y)/4
        for rho in (plus, minus):
            self.close(np.trace(rho), 1)
            self.assertGreater(np.linalg.eigvalsh(rho).min(), -1e-12)
        pair_expected = np.eye(4)/3-swap(2, 0, 1)/6
        errors = []
        for keep in ((0, 1), (0, 2), (1, 2)):
            left, right = partial(plus, [2]*3, keep), partial(minus, [2]*3, keep)
            errors.append(float(np.linalg.norm(left-right)))
            self.close(left, pair_expected)
            self.close(right, pair_expected)
        for keep in ((0,), (1,), (2,)):
            self.close(partial(plus, [2]*3, keep), np.eye(2)/2)
            self.close(partial(minus, [2]*3, keep), np.eye(2)/2)
        self.assertAlmostEqual(distance(plus, minus), 1)
        self.assertAlmostEqual(np.trace(y@plus).real, 1)
        self.assertAlmostEqual(np.trace(y@minus).real, -1)
        OBS['static_pair_equality'] = dict(maximum_pair_difference=max(errors),
            complete_joint_trace_distance=1., pair_marginal='I_4/3 - SWAP/6',
            all_pair_swap_expectations=0., logical_chiral_expectations=[1., -1.])

    def test_04_fixed_continuous_dynamics_reveals_relation(self):
        a, b, c, p, x, y, z = operators()
        h = a+b
        effect = (np.eye(8)-b)/2
        records, errors = [], []
        for t in (0., .13, .37, math.pi/4, 1.12):
            u = evolve(h, t)
            observed = []
            for sign in (1, -1):
                rho = (p+sign*y)/4
                out = u@rho@u.conj().T
                probability = float(np.trace(effect@out).real)
                predicted = .5-sign*math.sqrt(3)/4*math.sin(2*t)
                errors.append(abs(probability-predicted))
                self.assertAlmostEqual(probability, predicted)
                observed.append(probability)
            records.append(dict(time=t, pair_singlet_probabilities=observed))
        self.close(1j*(h@b-b@h), math.sqrt(3)*y)
        gap = abs(records[3]['pair_singlet_probabilities'][1]-records[3]['pair_singlet_probabilities'][0])
        self.assertAlmostEqual(gap, math.sqrt(3)/2)
        OBS['continuous_pair_readout'] = dict(hamiltonian='SWAP_12 + SWAP_23',
            maximum_formula_error=max(errors), records=records,
            contrast_at_pi_over_four=gap, target_dependent_switching_used=False,
            pair_effect_apparatus_and_time_calibration_derived=False)

    def test_05_continuous_closure_and_homogeneous_boundary(self):
        a, b, c, p, x, y, z = operators()
        v = encoding()
        rng = np.random.default_rng(43005)
        errors = []
        for weights in ((1., 1., 0.), (1., 1., 1.), (.2, -.4, .7)):
            ja, jb, jc = weights
            h = ja*a+jb*b+jc*c
            vector = np.array([math.sqrt(3)*(jb-jc)/2, 0., -ja+(jb+jc)/2])
            logical_h = sum(value*q for value, q in zip(vector, PAULI))
            self.close(v.conj().T@h@v, np.kron(np.eye(2), logical_h))
            frequency_squared = ja*ja+jb*jb+jc*jc-ja*jb-jb*jc-jc*ja
            self.assertAlmostEqual(float(vector@vector), frequency_squared)
            r = rng.normal(size=3)
            r *= .8/np.linalg.norm(r)
            rho = (p+sum(value*q for value, q in zip(r, (x, y, z))))/4
            derivative = -1j*(h@rho-rho@h)
            observed = np.array([np.trace(q@derivative).real for q in (x, y, z)])
            errors.append(float(np.linalg.norm(observed-2*np.cross(vector, r))))
            self.close(observed, 2*np.cross(vector, r))
        self.close((a+b+c)@p, np.zeros((8, 8)))
        OBS['closed_relational_dynamics'] = dict(maximum_generator_error=max(errors),
            all_equal_couplings_freeze_relational_qubit=True,
            unequal_couplings_are_explicit_model_input=True,
            bloch_components_are_position_coordinates=False)

    def test_06_reference_and_common_frame(self):
        a, b, c, p, x, y, z = operators()
        v = encoding()
        one = evolve(.21*PAULI[0]-.17*PAULI[1]+.31*PAULI[2], 1.)
        collective = np.kron(np.kron(one, one), one)
        self.close(collective@v, v@np.kron(one, np.eye(2)))
        for q in (p, x, y, z):
            self.close(collective@q@collective.conj().T, q)
        rng = np.random.default_rng(43006)
        # Arbitrary unknown state on gauge, logical qubit and a three-level reference.
        rho = density(12, rng)
        embed = np.kron(v, np.eye(3))
        physical = embed@rho@embed.conj().T
        t = .47
        whole = np.kron(evolve(a+b, t), np.eye(3))
        h_logical = math.sqrt(3)/2*PAULI[0]-.5*PAULI[2]
        expected_u = np.kron(np.kron(np.eye(2), evolve(h_logical, t)), np.eye(3))
        actual = whole@physical@whole.conj().T
        expected = embed@expected_u@rho@expected_u.conj().T@embed.conj().T
        error = float(np.linalg.norm(actual-expected))
        self.close(actual, expected)
        self.close(whole.conj().T@actual@whole, physical)
        OBS['reference_and_frame'] = dict(arbitrary_reference_dimension=3,
            full_reference_intertwining_error=error,
            collective_rotation_acts_only_on_gauge=True,
            independent_local_identifications_removed=False,
            unknown_physical_input_automatically_encoded=False)


def run():
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=430, baseline_round=429, status='exchange_selected_relational_qubit_and_continuous_readout',
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(known_exchange_encoding_explicitly_attributed=True,
            three_is_minimum_for_noncommuting_qubit_relations=True,
            identical_complete_pair_marginals_can_evolve_differently=True,
            fixed_continuous_hamiltonian_used=True,
            exact_programmable_processor_required=False,
            common_identification_and_couplings_derived=False,
            internal_readout_apparatus_completed=False,
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
