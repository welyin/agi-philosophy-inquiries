"""Round 435: composition of continuous singlet-block exchange interactions.

Second-order cross-edge terms vanish even at a shared block. A specified
three-edge loop produces a genuine third-order three-block interaction.
The graph and singlet preparation remain model inputs.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as core
from encoded_exchange_response_audit import code_block, norm

TARGET = Path(__file__).with_name('exchange_block_composition_audit_results.json')
OBS = {}


def kron_all(items):
    out = np.array([[1.]])
    for item in items:
        out = np.kron(out, item)
    return out


class Blocks:
    def __init__(self, count=3):
        self.count = count
        self.dimension = 16**count
        self.code = kron_all([code_block()]*count)
        self.c4 = sum(core.swap(4, a, b).real for a, b in itertools.combinations(range(4), 2))
        energies, self.basis = np.linalg.eigh(self.c4)
        self.energies = np.rint(energies)
        energy_sum = np.zeros([16]*count)
        for j in range(count):
            shape = [1]*count
            shape[j] = 16
            energy_sum += self.energies.reshape(shape)
        self.inverse_energies = np.zeros_like(energy_sum)
        np.divide(1., energy_sum, out=self.inverse_energies, where=energy_sum > 0)
        self.indices = np.arange(self.dimension, dtype=np.int64)
        self.permutations = {}

    def swap(self, a, b, columns):
        pair = tuple(sorted((a, b)))
        if pair not in self.permutations:
            i, j = 4*self.count-1-a, 4*self.count-1-b
            different = ((self.indices >> i) ^ (self.indices >> j)) & 1
            self.permutations[pair] = self.indices ^ (different << i) ^ (different << j)
        return columns[self.permutations[pair]]

    def interaction(self, links, columns):
        out = np.zeros_like(columns, dtype=np.result_type(columns.dtype, float))
        for a, b, weight in links:
            out += weight*(self.swap(a, b, columns)-columns/2)
        return out

    def block_op(self, operator, index, columns):
        tensor = columns.reshape([16]*self.count+[columns.shape[1]])
        tensor = np.moveaxis(np.tensordot(operator, tensor, axes=(1, index)), 0, index)
        return tensor.reshape(columns.shape)

    def h0(self, columns):
        return sum(self.block_op(self.c4, j, columns) for j in range(self.count))

    def inverse(self, columns):
        out = columns
        for j in range(self.count):
            out = self.block_op(self.basis.T, j, out)
        out = (out.reshape([16]*self.count+[columns.shape[1]])*self.inverse_energies[..., None]).reshape(columns.shape)
        for j in range(self.count):
            out = self.block_op(self.basis, j, out)
        return out


def block_matrices():
    v = code_block()
    return {(i, j): (3*np.eye(2) if i == j else 2*v.conj().T@core.swap(4, i, j)@v-np.eye(2))
            for i in range(4) for j in range(4)}


def pair_formula(weights):
    m = block_matrices()
    out = np.zeros((4, 4), complex)
    for i, j, k, l in itertools.product(range(4), repeat=4):
        out -= weights[i, j]*weights[k, l]*np.kron(m[i, k], m[j, l])/48
    return out


def triangle_links():
    # AB uses A1 B1, BC uses B2 C1, CA uses C2 A2 (one-based within blocks).
    return [[(0, 4, 1.)], [(5, 8, 1.)], [(9, 1, 1.)]]


def third_order_bound(epsilon, tau):
    e, s = abs(epsilon), abs(tau)
    return ((3/4+2781*s/1024)*e+(27/32+1701*s/4096)*e*e
            +(135/128+1215*s/8192)*e**3)


class Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blocks = Blocks(3)

    def close(self, a, b, tolerance=3e-11):
        self.assertLess(norm(a-b), tolerance)

    def test_01_local_scalar_projection(self):
        v, m = code_block(), block_matrices()
        sigma = {}
        for i in range(4):
            for alpha, pauli in enumerate(core.PAULI):
                sigma[i, alpha] = kron_all([pauli if j == i else np.eye(2) for j in range(4)])
                self.close(v.conj().T@sigma[i, alpha]@v, np.zeros((2, 2)))
        errors = []
        for i, j, a, b in itertools.product(range(4), range(4), range(3), range(3)):
            actual = v.conj().T@sigma[i, a]@sigma[j, b]@v
            expected = m[i, j]/3 if a == b else np.zeros((2, 2))
            self.close(actual, expected)
            errors.append(norm(actual-expected))
        OBS['local_identity'] = dict(vector_projection_zero=True,
            two_vector_identities_checked=len(errors), maximum_operator_error=max(errors),
            projected_pair_scalar='M_ii=3I; M_ij=2 S_ij(logical)-I for i != j')

    def test_02_second_order_shared_block_additivity(self):
        b, v = self.blocks, self.blocks.code
        edges = [([(0, 4, 1.), (1, 5, 1.)]),
                 ([(4, 8, 1.), (6, 9, 1.)]),
                 ([(0, 8, .3), (2, 11, -.7)])]
        images = [b.interaction(edge, v) for edge in edges]
        cross_errors = []
        for image in images:
            self.close(b.h0(image), 4*image)
            self.close(b.inverse(image), image/4)
        for i, j in itertools.permutations(range(3), 2):
            cross = images[i].conj().T@images[j]/4
            self.close(cross, np.zeros((8, 8)))
            cross_errors.append(norm(cross))
        total = sum(images)
        k2_direct = -total.conj().T@b.inverse(total)
        contributions = [-image.conj().T@image/4 for image in images]
        self.close(k2_direct, sum(contributions))
        commutator = norm(contributions[0]@contributions[1]-contributions[1]@contributions[0])
        self.assertGreater(commutator, .05)
        OBS['second_order_composition'] = dict(shared_block_cross_terms_checked=6,
            maximum_cross_term_norm=max(cross_errors),
            total_operator_formula_error=norm(k2_direct-sum(contributions)),
            two_contributions_commutator_norm=commutator,
            additive_does_not_mean_commuting=True,
            arbitrary_finite_block_graph_proof_in_note=True)

    def test_03_real_plane_boundary_of_second_order(self):
        b = Blocks(2)
        weights = np.array([[1, 2, 0, -1], [3, 0, 2, 1], [-2, 1, 1, 0], [0, 2, -1, 1]], float)/4
        links = [(i, 4+j, weights[i, j]) for i, j in itertools.product(range(4), repeat=2) if weights[i, j]]
        image = b.interaction(links, b.code)
        k2 = -image.conj().T@image/4
        self.close(k2, pair_formula(weights))
        paulis = [np.eye(2), *core.PAULI]
        y_coefficients = []
        for i, j in itertools.product(range(4), repeat=2):
            if 2 in (i, j):
                value = np.trace(np.kron(paulis[i], paulis[j])@k2)/4
                self.assertLess(abs(value), 1e-12)
                y_coefficients.append(float(abs(value)))
        yy_swap = float(np.trace(np.kron(core.PAULI[1], core.PAULI[1])@core.swap(2, 0, 1)).real/4)
        self.assertAlmostEqual(yy_swap, .5)
        uniform = [(i, 4+j, 1.) for i, j in itertools.product(range(4), repeat=2)]
        self.close(b.interaction(uniform, b.code), np.zeros_like(b.code))
        OBS['second_order_scope'] = dict(pair_formula_error=norm(k2-pair_formula(weights)),
            maximum_coefficient_with_a_logical_Y=max(y_coefficients),
            logical_swap_YY_coefficient=yy_swap,
            arbitrary_logical_SWAP_not_obtained_by_this_second_order_formula=True,
            fully_uniform_cross_block_coupling_annihilates_singlet_codes=True)

    def test_04_triangle_words_and_genuine_three_block_term(self):
        b, v, edges = self.blocks, self.blocks.code, triangle_links()
        images = [b.interaction(edge, v) for edge in edges]
        d = -np.eye(2)-2*core.PAULI[2]
        product = kron_all([d, d, d])
        k3 = np.zeros((8, 8), complex)
        counts = dict(single_edge_cubes=0, triangle_words=0, vanishing_other_words=0)
        for i, j, k in itertools.product(range(3), repeat=3):
            term = images[i].conj().T@b.interaction(edges[j], images[k])/16
            if i == j == k:
                expected = -3*np.eye(8)/64
                counts['single_edge_cubes'] += 1
            elif len({i, j, k}) == 3:
                expected = product/1152
                counts['triangle_words'] += 1
            else:
                expected = np.zeros((8, 8))
                counts['vanishing_other_words'] += 1
            self.close(term, expected)
            k3 += term
        expected = -9*np.eye(8)/64+product/192
        self.close(k3, expected)
        coefficient = float(np.trace(kron_all([core.PAULI[2]]*3)@k3).real/8)
        self.assertAlmostEqual(coefficient, -1/24)
        OBS['triangle'] = dict(ordered_words=counts,
            k3='-9 I/64 + (-I-2 Z_A)(-I-2 Z_B)(-I-2 Z_C)/192',
            k2='-9 I/16', genuine_ZZZ_coefficient=coefficient,
            formula_error=norm(k3-expected),
            three_block_interaction_present_at_order=3)

    def test_05_exact_integer_coefficient_certificate(self):
        b = self.blocks
        elementary = np.eye(16, dtype=np.int64)
        zero = elementary[5]-elementary[6]-elementary[9]+elementary[10]
        one = 2*elementary[3]+2*elementary[12]-elementary[5]-elementary[6]-elementary[9]-elementary[10]
        unnormalized = np.kron(np.kron(np.stack([zero, one], axis=1),
                                      np.stack([zero, one], axis=1)), np.stack([zero, one], axis=1))
        norm_squares = np.array([64*3**j.bit_count() for j in range(8)], dtype=np.int64)
        self.assertTrue(np.array_equal(unnormalized.T@unnormalized, np.diag(norm_squares)))
        links = [(a, c) for edge in triangle_links() for a, c, _ in edge]
        def doubled_t(columns):
            return sum(2*b.swap(a, c, columns)-columns for a, c in links)
        second = unnormalized.T@doubled_t(doubled_t(unnormalized))
        third = unnormalized.T@doubled_t(doubled_t(doubled_t(unnormalized)))
        self.assertTrue(np.array_equal(second, 9*np.diag(norm_squares)))
        self.assertTrue(np.array_equal(third, np.diag(np.diag(third))))
        exact_values = []
        for j in range(8):
            product = math.prod(1 if (j >> bit) & 1 else -3 for bit in range(3))
            expected = -Fraction(9, 64)+Fraction(product, 192)
            actual = Fraction(int(third[j, j]), 128*int(norm_squares[j]))
            self.assertEqual(actual, expected)
            exact_values.append(actual)
        zzz_coefficient = sum(((-1)**j.bit_count())*value for j, value in enumerate(exact_values))/8
        self.assertEqual(zzz_coefficient, -Fraction(1, 24))
        OBS['exact_certificate'] = dict(integer_arithmetic=True,
            code_norm_squares=norm_squares.tolist(),
            k3_diagonal=[str(value) for value in exact_values],
            exact_ZZZ_coefficient=str(zzz_coefficient),
            maximum_integer_magnitude=int(np.max(np.abs(third))),
            floating_fit_used=False)

    def test_06_third_order_uniform_time_bound(self):
        b, v = self.blocks, self.blocks.code
        links = sum(triangle_links(), [])
        action = lambda columns: b.interaction(links, columns)
        tv = action(v)
        k2 = -tv.conj().T@tv/4
        self.close(k2, -9*np.eye(8)/16)
        k3 = tv.conj().T@action(tv)/16
        w1 = -b.inverse(tv)
        w2 = -b.inverse(action(w1))
        w3 = -b.inverse(action(w2)+(9/16)*w1)
        self.close(b.h0(w1)+tv, np.zeros_like(v))
        self.close(b.h0(w2)+action(w1)+(9/16)*v, np.zeros_like(v))
        self.close(b.h0(w3)+action(w2)+(9/16)*w1, v@k3)
        a4 = action(w3)+(9/16)*w2-w1@k3
        a5 = (9/16)*w3-w2@k3
        a6 = -w3@k3
        constants = [('w1', w1, 3/8), ('w2', w2, 27/64), ('w3', w3, 135/256),
                     ('a4', a4, 2781/1024), ('a5', a5, 1701/4096), ('a6', a6, 1215/8192)]
        for _, matrix, upper in constants:
            self.assertLessEqual(norm(matrix), upper+1e-12)
        e = 1/13
        w = v+e*w1+e**2*w2+e**3*w3
        residual = b.h0(w)+e*action(w)+(9/16)*e**2*w-w@(e**3*k3)
        self.close(residual, e**4*a4+e**5*a5+e**6*a6)
        witness = third_order_bound(1/512, 3*math.pi)
        self.assertLess(witness, .052)
        OBS['third_order_dynamics'] = dict(
            scalar_energy_shift_removed='-9 epsilon^2 I/16',
            uniform_error_bound='(3/4+2781|tau|/1024)|epsilon| + (27/32+1701|tau|/4096)epsilon^2 + (135/128+1215|tau|/8192)|epsilon|^3',
            effective_time='tau=epsilon^3*t for epsilon>0',
            bound_at_epsilon_inverse_512_tau_3pi=witness,
            constants=[dict(name=name, numerical_norm=norm(matrix), analytic_upper_bound=upper)
                       for name, matrix, upper in constants],
            arbitrary_unknown_code_input_and_reference_covered=True,
            exact_polynomial_identity_error=norm(residual-(e**4*a4+e**5*a5+e**6*a6)),
            full_long_time_spectral_exponentiation_claimed=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=435, baseline_round=434, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(second_order_cross_edge_terms_vanish_for_singlet_blocks=True,
            simultaneous_noncommuting_pair_terms_allowed=True,
            third_order_triangle_generates_genuine_three_block_interaction=True,
            third_order_uniform_unknown_input_dynamical_bound_proved=True,
            integer_coefficient_certificate_verified=True,
            exact_433_network_implemented=False,
            selected_block_graph_and_port_weights_are_inputs=True,
            gap_initial_code_and_timescale_are_free=False,
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
