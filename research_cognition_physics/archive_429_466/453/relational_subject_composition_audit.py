"""Round 453: exact autonomy of two fixed three-qubit relational interfaces.

Baseline 451; independent of round 452. No new cognitive axiom or no-go for
encoded, endpoint, prepared-gauge, enlarged-interface or low-energy schemes.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools as it
import json
from pathlib import Path
import platform
import unittest

import numpy as np

import agreement_exchange_audit as rank_tools
import exchange_relation_audit as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / "relational_subject_composition_audit_results.json"
OBS = {}


def short(value):
    return float(f"{float(value):.12g}")


def fraction_data(value):
    return dict(numerator=value.numerator, denominator=value.denominator,
                decimal=short(value))


def swap(n, i, j):
    return old.swap(n, i, j).real.astype(np.int64)


def site_pauli(n, site, mu):
    result = np.array([[1]], dtype=complex)
    for i in range(n):
        result = np.kron(result, old.PAULI[mu] if i == site else np.eye(2))
    return result


@lru_cache(None)
def model():
    local = [swap(3, i, j) for i, j in ((0, 1), (1, 2), (0, 2))]
    p3 = 3*np.eye(8, dtype=np.int64)-sum(local)
    first = np.kron(np.array([0, 1, -1, 0]), [1, 0])
    second = np.kron(np.array([0, 1, -1, 0]), [0, 1])
    diff = local[1]-local[2]
    integer_v = np.column_stack((first, diff@first, second, diff@second))
    norms = np.array([2, 6, 2, 6])
    v = integer_v/np.sqrt(norms)[None, :]
    # Original order GA,LA,GB,LB -> grouped GA,GB,LA,LB.
    order = [ga*8+la*4+gb*2+lb for ga, gb, la, lb in it.product(range(2), repeat=4)]
    w = np.kron(v, v)[:, order]
    p9 = np.kron(p3, p3)
    cross = [swap(6, a, b+3) for a in range(3) for b in range(3)]
    comms = np.column_stack([(s@p9-p9@s).ravel() for s in cross])
    eye, x, z = np.eye(2), old.PAULI[0], old.PAULI[2]
    aa = [(eye-z)/3-x/np.sqrt(3), (eye-z)/3+x/np.sqrt(3), (eye+2*z)/3]
    cg = sum((np.kron(q, q) for q in old.PAULI), np.zeros((4, 4), complex))
    return local, p3, integer_v, v, w, p9, cross, comms, aa, cg


class Audit(unittest.TestCase):
    def test_01_fixed_relation_code_and_vector_projection(self):
        local, p3, iv, v, w, p9, cross, comms, aa, cg = model()
        np.testing.assert_array_equal(p3@p3, 3*p3)
        np.testing.assert_array_equal(iv.T@iv, np.diag([2, 6, 2, 6]))
        np.testing.assert_array_equal(p3@iv, 3*iv)
        np.testing.assert_allclose(v, old.encoding(), atol=1e-15)
        np.testing.assert_allclose(w.conj().T@w, np.eye(16), atol=2e-15)
        np.testing.assert_allclose(w@w.conj().T, p9/9, atol=2e-15)
        worst = 0.0
        for a, mu in it.product(range(3), repeat=2):
            actual = v.conj().T@site_pauli(3, a, mu)@v
            expected = np.kron(old.PAULI[mu], aa[a])
            worst = max(worst, float(np.linalg.norm(actual-expected)))
        self.assertLess(worst, 4e-15)
        np.testing.assert_allclose(sum(aa), np.eye(2), atol=1e-15)
        OBS["fixed_interface"] = dict(raw_dimension=64, original_joint_code_dimension=16,
            retained_logic_dimension=4, discarded_gauge_dimension=4,
            integer_code_columns=iv.tolist(), code_column_norms_squared=[2, 6, 2, 6],
            vector_projection_max_error=short(worst),
            spin_not_assumed_to_be_physical_space_direction=True)

    def test_02_exact_all_weight_leakage_gram(self):
        *_, comms, aa, cg = model()
        gram = comms.T@comms
        expected = np.empty((9, 9), dtype=np.int64)
        pairs = list(it.product(range(3), repeat=2))
        for i, (a, b) in enumerate(pairs):
            for j, (c, d) in enumerate(pairs):
                expected[i, j] = 2592*(a == c and b == d)-720*((a == c)+(b == d))+192
        np.testing.assert_array_equal(gram, expected)
        ones, eye = np.ones((3, 3), dtype=np.int64), np.eye(3, dtype=np.int64)
        rowcol9 = 3*(np.kron(eye, ones)+np.kron(ones, eye))-2*np.kron(ones, ones)
        double9 = np.kron(3*eye-ones, 3*eye-ones)
        np.testing.assert_array_equal(9*gram, 432*rowcol9+2592*double9)
        np.testing.assert_array_equal(gram@np.ones(9, dtype=int), np.zeros(9, dtype=int))
        rank = rank_tools.modular_rank(comms.astype(complex))[0]
        self.assertEqual(rank, 8)
        # Check the arbitrary-weight formula with exact rational row/column decomposition.
        values = [F(v) for v in (1, -2, 3, 0, 4, -1, 2, 5, -3)]
        average = sum(values)/9
        u = [sum(values[3*a+b] for b in range(3))/3-average for a in range(3)]
        vv = [sum(values[3*a+b] for a in range(3))/3-average for b in range(3)]
        centered = [values[3*a+b]-average-u[a]-vv[b] for a, b in pairs]
        predicted = 16*(sum(x*x for x in u)+sum(x*x for x in vv))+32*sum(x*x for x in centered)
        actual = sum(values[i]*int(gram[i,j])*values[j] for i,j in it.product(range(9),repeat=2))/81
        self.assertEqual(actual, predicted)
        OBS["exact_code_retention"] = dict(integer_projection_scale=9,
            all_cross_weight_commutator_gram=gram.tolist(), constraint_rank_mod1009=rank,
            kernel="all nine J_ab equal", nonzero_gram_eigenvalues=[432, 2592],
            multiplicities=[4, 4],
            physical_commutator_norm_squared="16*(||u||^2+||v||^2)+32*||D||_F^2",
            mean_initial_leakage_t2_coefficient="(||u||^2+||v||^2)/2+||D||_F^2",
            arbitrary_real_weights_classified_not_sampled=True)

    def test_03_compressed_gauge_dependence_and_logical_kernel(self):
        local, p3, iv, v, w, p9, cross, comms, aa, cg = model()
        worst = 0.0
        for index, (a, b) in enumerate(it.product(range(3), repeat=2)):
            actual = w.conj().T@cross[index]@w
            expected = np.eye(16)/2+np.kron(cg, np.kron(aa[a], aa[b]))/2
            worst = max(worst, float(np.linalg.norm(actual-expected)))
        self.assertLess(worst, 8e-15)
        # Rows of R are coefficients of 3 A_a in basis I, X/sqrt(3), Z.
        coefficients = np.array([[1, -3, -1], [1, 3, -1], [1, 0, 2]], dtype=np.int64)
        expansion = np.kron(coefficients.T, coefficients.T)
        non_scalar = expansion[1:, :]
        rank = rank_tools.modular_rank(non_scalar.astype(complex))[0]
        self.assertEqual(rank, 8)
        np.testing.assert_array_equal(non_scalar@np.ones(9, dtype=int), np.zeros(8, dtype=int))
        np.testing.assert_array_equal((cg+3*np.eye(4))@(cg-np.eye(4)), np.zeros((4,4)))
        OBS["autonomous_logic"] = dict(cross_projection_max_error=short(worst),
            gauge_eigenvalues=[-3, 1], A_coefficient_matrix=coefficients.tolist(),
            non_scalar_K_constraint_rank_mod1009=rank,
            gauge_independent_reduced_flow_iff_K_scalar=True,
            K_scalar_iff_all_cross_weights_equal=True,
            compressed_H_automatically_equal_to_true_restricted_dynamics=False,
            product_gauge_inputs_already_detect_failure=True)

    def test_04_uniform_cross_weights_preserve_all_unknown_inputs_and_reference(self):
        local, p3, iv, v, w, p9, cross, comms, aa, cg = model()
        h_a = sum((c*s for c,s in zip((1,-2,3),local)), np.zeros((8,8),dtype=int))
        h_b = sum((c*s for c,s in zip((2,1,-1),local)), np.zeros((8,8),dtype=int))
        h = np.kron(h_a,np.eye(8))+np.kron(np.eye(8),h_b)+2*sum(cross)
        hla = (v.conj().T@h_a@v)[:2,:2]
        hlb = (v.conj().T@h_b@v)[:2,:2]
        hl = np.kron(hla,np.eye(2))+np.kron(np.eye(2),hlb)
        hg = 9*np.eye(4)+cg
        target = np.kron(hg,np.eye(4))+np.kron(np.eye(4),hl)
        np.testing.assert_allclose(h@w, w@target, atol=2e-14)
        full_u, ug, ul = old.evolve(h,.37), old.evolve(hg,.37), old.evolve(hl,.37)
        error = float(np.linalg.norm(full_u@w-w@np.kron(ug,ul),2))
        self.assertLess(error, 2e-14)
        rng = np.random.default_rng(453)
        rho = old.density(48,rng)  # arbitrary correlated gauge, logic, reference mixed input
        u = np.kron(np.kron(ug,ul),np.eye(3))
        actual = old.partial(u@rho@u.conj().T,[4,4,3],(1,2))
        logical_ref = old.partial(rho,[4,4,3],(1,2))
        expected = np.kron(ul,np.eye(3))@logical_ref@np.kron(ul.conj().T,np.eye(3))
        ref_error = float(np.linalg.norm(actual-expected))
        self.assertLess(ref_error, 2e-14)
        OBS["sufficient_branch"] = dict(uniform_cross_weight=2,
            nontrivial_gauge_interaction_allowed=True,
            logical_generator_is_sum_of_two_local_generators=True,
            full_input_intertwining_error=short(error), reference_dimension=3,
            correlated_mixed_reference_error=short(ref_error),
            exact_operator_proof_covers_any_reference=True)

    def test_05_actual_local_readout_sees_unknown_gauge_without_postselection(self):
        local, p3, iv, v, w, p9, cross, comms, aa, cg = model()
        s = swap(6,0,5)
        # O = Y_relation,A/sqrt(3) on A, identity on B; numerator is Gaussian integer.
        o3 = np.kron(1j*(local[0]@local[1]-local[1]@local[0]), np.eye(8,dtype=int))
        e6 = 3*np.eye(64,dtype=int)+o3
        np.testing.assert_array_equal(o3,o3.conj().T)
        np.testing.assert_array_equal(o3@o3,np.kron(p3,np.eye(8,dtype=int)))
        probabilities, leakage, derivatives = [], [], []
        for gauge_b in (0,1):
            raw = np.kron(iv[:,0],iv[:,2*gauge_b]) # normalized input = raw/2
            self.assertEqual(int(raw@raw),4)
            np.testing.assert_array_equal(p9@raw,9*raw)
            self.assertEqual(np.vdot(raw,o3@raw),0)
            self.assertEqual(np.vdot(s@raw,o3@(s@raw)),0)
            derivative_num = np.vdot(raw,1j*(s@o3-o3@s)@raw)
            self.assertEqual(derivative_num.imag,0)
            derivatives.append(F(int(derivative_num.real),12))
            # At pi/4, normalized output is (raw-i S raw)/sqrt(8).
            output = raw-1j*s@raw
            numerator = np.vdot(output,e6@output)
            self.assertEqual(numerator.imag,0)
            probabilities.append(F(int(numerator.real),48))
            survival = np.vdot(output,p9@output)
            self.assertEqual(survival.imag,0)
            leakage.append(1-F(int(survival.real),72))
        self.assertEqual(probabilities,[F(7,12),F(5,12)])
        self.assertEqual(leakage,[F(1,3),F(1,6)])
        self.assertEqual(derivatives,[F(1,3),-F(1,3)])
        gap_window = F(1,6)*(1-F(1,4)**2/2)
        self.assertEqual(gap_window,F(31,192))
        self.assertGreater(gap_window,F(1,8))
        OBS["physical_gauge_witness"] = dict(cross_edge=[0,5], time="pi/4",
            same_initial_logic="|0_L 0_L>", gauge_inputs=["|0_G 0_G>","|0_G 1_G>"],
            initial_complete_subjects_independent=True,
            fixed_local_A_effect="(I + Y_relation,A/sqrt(3))/2",
            probabilities=[fraction_data(q) for q in probabilities],
            exact_probability_difference=fraction_data(F(1,6)),
            original_joint_code_leakage=[fraction_data(q) for q in leakage],
            all_time_probability_difference="sin(2*t)/6",
            window="[pi/4-1/8,pi/4+1/8]", window_lower=fraction_data(gap_window),
            physical_effect_defined_on_all_raw_states=True,
            postselection_or_physical_PHP_projection_used=False,
            preparation_and_readout_still_inputs=True)

    def test_06_endpoint_exchange_does_not_obey_all_time_contract(self):
        local, p3, iv, v, w, p9, cross, comms, aa, cg = model()
        aligned = [swap(6,i,i+3) for i in range(3)]
        h = sum(aligned)
        self.assertTrue(np.any(h@p9-p9@h))
        block = aligned[0]@aligned[1]@aligned[2]
        np.testing.assert_array_equal(block@p9,p9@block)
        # Product exp[-i(pi/2)S_r]=i S_block.
        endpoint = 1j*block
        factorized = np.kron(rank_tools.swap(2),rank_tools.swap(2))
        np.testing.assert_allclose(endpoint@w,1j*w@factorized,atol=1e-15)
        numerical = old.evolve(h,np.pi/2)
        self.assertLess(float(np.linalg.norm(numerical-endpoint,2)),2e-14)
        OBS["boundary"] = dict(aligned_three_pair_H_leaks_during_evolution=True,
            endpoint_time="pi/2", endpoint="i * SWAP_gauge * SWAP_logic",
            endpoint_covers_full_unknown_code_and_reference=True,
            exact_endpoint_exchange_excluded_by_theorem=False,
            pulsed_encoded_entangling_gates_excluded=False,
            enlarged_interface_or_prepared_gauge_excluded=False,
            low_energy_434_436_results_preserved=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream,verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=453,baseline_round=451,date="2026-09-24",
        runtime=dict(python=platform.python_version(),numpy=np.__version__),
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(two_fixed_three_qubit_relational_interfaces_only=True,
            arbitrary_static_real_pair_exchange_weights=True,
            exact_all_time_code_retention_classified=True,
            unknown_gauge_independent_compressed_logic_classified=True,
            exact_code_retention_already_forces_no_cross_logical_interaction=True,
            nontrivial_macro_interaction_obtained=False,
            all_possible_relational_organizations_ruled_out=False,
            original_logic_interface_preservation_imposed_as_cognitive_axiom=False,
            new_subsystem_or_exchange_gate_theory_claimed=False,
            depends_on_round_452=False,
            full_soca_or_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run",action="store_true")
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2))
