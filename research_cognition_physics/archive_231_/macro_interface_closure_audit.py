"""Round 452: exact output algebra and physical relation-coherence witness.
The subalgebra conclusion does not grant tomography or arbitrary control.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools as it
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

import protected_partner_processing_audit as previous
import agreement_exchange_audit as ranks
import whole_subject_exchange_audit as transport

HERE = Path(__file__).resolve().parent
TARGET = HERE / "macro_interface_closure_audit_results.json"
OBS = {}


@lru_cache(None)
def model():
    p, bits, words, active, record, local, packet, c, z, w1, w2, k2, d, b2, *_ = previous.model()
    match = previous.old.matching_words(4)
    a0 = z.T @ (active[:, None]*z)
    l0 = z.T @ local @ z
    effects = [np.array([b[i] for _ in match for b in bits], dtype=np.int64) for i in range(4)]
    permutations = [transport.permutation_matrix(transport.data_perm(m, 2)).astype(np.int64)
                    for m in match]
    return bits, words, active, record, packet, z, k2, a0, l0, b2, effects, permutations


def small_edges():
    result = []
    for m, n in it.combinations(range(3), 2):
        a = np.zeros((3, 3), dtype=np.int64)
        a[m, n] = a[n, m] = 1
        result.append((m, n, a))
    return result


def null_vector(matrix):
    a = [[F(int(v)) for v in row] for row in matrix]
    row = 0
    pivots = []
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        value = a[row][col]
        a[row] = [v/value for v in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                coefficient = a[i][col]
                a[i] = [v-coefficient*w for v, w in zip(a[i], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    free = next(i for i in range(len(a[0])) if i not in pivots)
    v = [F(0) for _ in a[0]]
    v[free] = F(1)
    for i, col in enumerate(pivots):
        v[col] = -a[i][free]
    scale = math.lcm(*(x.denominator for x in v))
    values = [int(x*scale) for x in v]
    divisor = math.gcd(*values)
    return np.array([x//divisor for x in values], dtype=np.int64), len(pivots)


@lru_cache(None)
def phase_certificate():
    *_, h2, effects, permutations = model()
    first, second = 1, 24
    order, time = 20, F(1, 8)
    power = np.eye(48, dtype=object)[:, [first, second]]
    re, im = np.zeros((48, 2), dtype=object), np.zeros((48, 2), dtype=object)
    for k in range(order+1):
        coefficient = F((-1)**(k//2), 16**k*math.factorial(k))
        if k % 2:
            im -= coefficient*power
        else:
            re += coefficient*power
        power = h2.astype(object)@power
    probabilities = []
    for sign in (1, -1):
        r = re[:, 0]-sign*im[:, 1]
        v = im[:, 0]+sign*re[:, 1]
        probabilities.append(sum((F(r[i]**2+v[i]**2, 2)
                                  for i in np.flatnonzero(effects[0])), F(0)))
    # ||B|| <= 7, t=1/8, exp(7/8)<3: uniform vector remainder, not entrywise.
    remainder = 3*F(7, 8)**21/math.factorial(21)
    lower = probabilities[0]-probabilities[1]-4*remainder-2*remainder**2
    delta = F(1, 1000)
    b = previous.error_bound(F(1, 1024), time+delta)
    actual_lower = F(3, 25)-6*delta-2*b
    return probabilities, remainder, lower, actual_lower


class Audit(unittest.TestCase):
    def test_01_existing_physical_frame_and_small_group_decomposition(self):
        bits, words, active, record, packet, z, k2, a0, l0, h2, effects, perms = model()
        reconstructed = -4*np.eye(48, dtype=np.int64)
        single = [i for i, b in enumerate(bits) if sum(b) == 1]
        coefficient_matrices = []
        for m, n, edge in small_edges():
            move = perms[m]@perms[n]
            reconstructed -= np.kron(edge, move)
            coefficient_matrices.append(move[np.ix_(single, single)])
        np.testing.assert_array_equal(reconstructed, k2)
        rank, pivots = ranks.modular_rank(np.column_stack([q.reshape(-1) for q in coefficient_matrices]))
        self.assertEqual(rank, 3)
        for site, effect in enumerate(effects):
            actual = np.array([int(w[2*site] == -1) for w in words], dtype=int)
            np.testing.assert_array_equal(np.diag(z.T@(actual[:, None]*z)), effect)
        one_indices = [16*m+x for m in range(3) for x in single]
        np.testing.assert_array_equal(a0[np.ix_(one_indices, one_indices)], 2*np.eye(12, dtype=int))
        # The single-bit flip connects the 16 bit labels; no K or h_A connects different weights.
        for x, y in it.product(range(16), repeat=2):
            distance = sum(a != b for a,b in zip(bits[x], bits[y]))
            block = l0[x::16, y::16]
            np.testing.assert_array_equal(block, np.eye(3, dtype=int)*(distance == 1))
            if sum(bits[x]) != sum(bits[y]):
                np.testing.assert_array_equal(k2[x::16, y::16], np.zeros((3,3),dtype=int))
        OBS["physical_interface"] = dict(code_dimension=48, local_blank_effects=4,
            controlled_frame_not_confused_with_receiver_frame=True,
            double_transposition_coefficient_rank=rank, modular_pivots=pivots,
            active_potential_on_one_excitation="2 I_12",
            internal_flip_connects_all_16_logic_labels=True)

    def test_02_exact_joint_commutant_and_parameter_independence(self):
        matrices = [a for _, _, a in small_edges()]
        rank, pivots = ranks.modular_rank(ranks.commutator_constraints(matrices))
        self.assertEqual(rank, 8)
        ident = np.eye(3, dtype=int)
        for a in matrices:
            np.testing.assert_array_equal(a@ident-ident@a, np.zeros_like(a))
        # Explicit matrix-unit generation complements the modular lower rank certificate.
        a01, a02, a12 = matrices
        e00 = (a01@a01+a02@a02-a12@a12)//2
        np.testing.assert_array_equal(e00, np.diag([1,0,0]))
        self.assertEqual(int((e00@a01)[0,1]), 1)
        self.assertEqual(int((e00@a02)[0,2]), 1)
        OBS["minimal_quantum_algebra"] = dict(matching_commutator_rank=rank,
            modulus=1009, pivots=pivots, scalar_kernel_dimension=1,
            finite_Cstar_invariant_algebra_lemma_used=True,
            minimal_algebra="M_48(C)", complex_vector_dimension=2304,
            iff_condition="kappa != 0 and ell != 0; j arbitrary real",
            classification_is_analytic_not_parameter_scan=True,
            arbitrary_affine_or_non_subalgebra_compression_excluded=False)

    def test_03_missing_process_boundaries(self):
        _, _, _, _, _, _, k2, a0, l0, _, effects, _ = model()
        ident = np.eye(16, dtype=int)
        match_projector = np.kron(np.diag([1,0,0]), ident)
        no_rematch_h = a0+l0
        np.testing.assert_array_equal(match_projector@no_rematch_h, no_rematch_h@match_projector)
        weight = np.diag(sum(effects))
        no_flip_h = k2+2*a0
        np.testing.assert_array_equal(weight@no_flip_h, no_flip_h@weight)
        for e in effects:
            np.testing.assert_array_equal(match_projector@np.diag(e), np.diag(e)@match_projector)
        self.assertGreater(np.ptp(np.diag(weight)), 0)
        OBS["boundary"] = dict(kappa_zero_has_matching_sector_projector=True,
            ell_zero_has_nonconstant_Hamming_weight_conservation=True,
            uniform_all_subject_readout_scope_fixed=True,
            arbitrary_small_coupling_finite_time_effect_not_inferred_from_exact_algebra=True)

    def test_04_coherent_relation_record_and_rational_window(self):
        bits, *rest = model()
        h2, effects, _ = rest[-3:]
        real2 = np.zeros((48,48), dtype=np.int64)
        real2[1,1] = real2[24,24] = 1
        imag2 = np.zeros_like(real2)
        imag2[1,24], imag2[24,1] = -1, 1
        # rho_+/-=(real2 +/- i imag2)/2, both pure.
        rho = (real2+1j*imag2)/2
        np.testing.assert_allclose(rho@rho, rho, atol=0)
        def marginal_data(q):
            return np.trace(q.reshape(3,16,3,16), axis1=0, axis2=2)
        def marginal_match(q):
            return np.trace(q.reshape(3,16,3,16), axis1=1, axis2=3)
        np.testing.assert_array_equal(marginal_data(imag2), np.zeros((16,16)))
        np.testing.assert_array_equal(marginal_match(imag2), np.zeros((3,3)))
        for m in range(3):
            np.testing.assert_array_equal(imag2[16*m:16*(m+1),16*m:16*(m+1)], np.zeros((16,16)))
        probs, remainder, lower, actual = phase_certificate()
        self.assertGreater(lower, F(3,25))
        self.assertGreater(actual, F(1,10))
        OBS["relation_phase_witness"] = dict(indices=[1,24],
            physical_bits=[list(bits[1]),list(bits[8])], relative_phases=["+i","-i"],
            full_matching_marginal_equal=True, full_data_marginal_equal=True,
            complete_matching_dephased_quantum_data_equal=True,
            tau_center="1/8", tau_window=["31/250","63/500"],
            taylor_order=20, uniform_remainder=previous.old.fraction_data(remainder),
            central_probability_numeric=[previous.old.short(float(v)) for v in probs],
            effective_gap_strict_lower="3/25",
            actual_gap_lower=previous.old.fraction_data(actual),
            actual_gap_strict_lower="1/10", effect="physical subject 0 A-register blank",
            derivative_gap_Lipschitz_upper=6, no_postselection=True)

    def test_05_full_raw_dynamics_transfer(self):
        _, words, active, record, packet, z, _, _, l0, h2, effects, _ = model()
        local = previous.model()[5]
        h = np.diag(1024*record)+packet+(np.diag(active)+local)/1024
        u = previous.old.unitary(h, 128)
        v = previous.old.unitary(h2/2, 1/8)
        psi = np.zeros((48,2),dtype=complex)
        psi[1] = 1/np.sqrt(2)
        psi[24] = np.array([1j,-1j])/np.sqrt(2)
        actual, ideal = u@z@psi, z@v@psi
        effect = np.array([int(w[0] == -1) for w in words])
        probabilities = (effect[:,None]*np.abs(actual)**2).sum(axis=0)
        errors = np.linalg.norm(actual-ideal, axis=0)
        self.assertGreater(probabilities[0]-probabilities[1], 1/10)
        self.assertLess(float(max(errors)), float(previous.error_bound(F(1,1024), F(1,8))))
        OBS["raw_transfer"] = dict(raw_single_marker_dimension=384, Delta=1024,
            g=1, J_A_and_lambda="1/1024", physical_center_time=128,
            numerical_probabilities=[previous.old.short(x) for x in probabilities],
            numerical_gap=previous.old.short(probabilities[0]-probabilities[1]),
            numerical_vector_errors=[previous.old.short(x) for x in errors],
            all_unknown_input_reference_bound_reused_from_round_449=True,
            actual_measurement_or_input_preparation_device_derived=False)

    def test_06_full_algebra_does_not_grant_single_wait_tomography(self):
        *_, h2, effects, _ = model()
        powers = [np.eye(48,dtype=np.int64)]
        for _ in range(5):
            powers.append(powers[-1]@h2)
        constraints = np.array([[np.trace(q) for q in powers]]+
            [[np.dot(e,np.diag(q)) for q in powers] for e in effects], dtype=np.int64)
        coeff, rank = null_vector(constraints)
        d = sum((int(c)*q for c,q in zip(coeff,powers)), np.zeros((48,48),dtype=np.int64))
        self.assertTrue(np.any(d))
        np.testing.assert_array_equal(d@h2, h2@d)
        self.assertEqual(int(np.trace(d)),0)
        for e in effects:
            self.assertEqual(int(np.dot(e,np.diag(d))),0)
        bound = int(np.max(np.sum(np.abs(d),axis=1)))
        # rho +/-= I/48 +/- d/(48 bound). Gershgorin proves positivity exactly.
        self.assertGreater(bound,0)
        gram = int(np.trace(d@d))
        self.assertGreater(gram,0)
        OBS["access_boundary"] = dict(stationary_polynomial_coefficients=coeff.tolist(),
            polynomial_constraint_rank=rank, integer_operator_row_norm_upper=bound,
            normalization_denominator=48*bound, difference_Hilbert_Schmidt_squared_numerator=4*gram,
            difference_Hilbert_Schmidt_squared_denominator=(48*bound)**2,
            positivity_by_integer_row_sum_bound=True,
            all_four_single_wait_output_histories_identical=True,
            full_state_tomography_from_these_histories=False,
            nontrivial_stationary_invisible_direction_dimension_at_least=43)


def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=452,baseline_round=451,date="2026-09-24",
        runtime=dict(python=platform.python_version(),numpy=np.__version__),
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(existing_physical_readout_interface_used=True,
            exact_minimal_invariant_Cstar_algebra_classified=True,
            matching_coherence_future_physical_signal_certified=True,
            actual_raw_finite_window_gap_transferred=True,
            linear_statistics_and_quantum_algebra_distinguished=True,
            all_arbitrary_quantum_compressions_excluded=False,
            complete_macroscopic_agent_or_SoCA_implemented=False,
            macro_full_tomography_or_control_derived=False,
            new_rule_derived_from_cognitive_principles=False,
            full_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run",action="store_true")
    args=parser.parse_args()
    report=run()
    if not args.dry_run:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2))

