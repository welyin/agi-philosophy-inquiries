"""1022: conditional chiral U(1) species selection by one Higgs and anomalies.

Exact charge arithmetic and finite matrix calibrations accompany an analytic
classification; finite examples do not prove the unbounded existence theorem.
No SM selection, physical Wilson instrument, or UV completion is claimed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import reduce
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "chiral_higgs_selection_results.json"
HISTORICAL = (
    "archive_531_553/research_note_531.md",
    "archive_585_628/research_note_627.md",
    "archive_1009_/research_note_1012.md",
    "archive_1009_/1021/research_round_1021_checks.json",
    "archive_1009_/1021/input_dependency_update_v0_10.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
)
TEMPLATE_8 = [(4, -2), (5, -3), (-1, -1), (4, -6)]
TEMPLATE_12 = [(4, -2), (5, -3), (13, -11), (-1, -1), (-1, -1), (12, -14)]


def norm(matrix):
    return float(np.linalg.norm(matrix))


def exact_determinant(matrix):
    """Independent rational elimination, not a product of the known pair blocks."""
    a = [[F(x) for x in row] for row in matrix]
    size = len(a); result = F(1)
    for j in range(size):
        pivot = next((i for i in range(j, size) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]; result = -result
        diagonal = a[j][j]; result *= diagonal
        for i in range(j + 1, size):
            multiplier = a[i][j] / diagonal
            if multiplier:
                for k in range(j + 1, size):
                    a[i][k] -= multiplier*a[j][k]
                a[i][j] = F(0)
    return result


def charge_data(charges):
    q = [F(x) for x in charges]
    return dict(charges=[int(x) if x.denominator == 1 else str(x) for x in q],
                linear_anomaly=str(sum(q, F(0))),
                cubic_anomaly=str(sum((x**3 for x in q), F(0))),
                quadratic_sum=str(sum((x*x for x in q), F(0))),
                no_zero_charge=all(x != 0 for x in q),
                no_opposite_charges=all(-x not in q for x in q),
                primitive_integer_charge_gcd=reduce(math.gcd, (abs(int(x)) for x in q))
                if all(x.denominator == 1 for x in q) else None)


def matrix_from_terms(size, terms):
    matrix = np.zeros((size, size), dtype=int)
    for i, j, mass in terms:
        matrix[i, j] += mass
        if i != j:
            matrix[j, i] += mass
    return matrix


def same_charge_unitary(charges, seed):
    """A complex unitary inside each repeated-charge subspace."""
    rng = np.random.default_rng(seed)
    U = np.eye(len(charges), dtype=complex)
    for charge in sorted(set(charges)):
        positions = [i for i, q in enumerate(charges) if q == charge]
        z = rng.normal(size=(len(positions), len(positions)))
        z = z + 1j*rng.normal(size=z.shape)
        block, _ = np.linalg.qr(z)
        U[np.ix_(positions, positions)] = block
    return U


def mass_calibration(charges, terms, higgs_charge, seed):
    q = np.array(charges, dtype=float); Q = np.diag(q); h = float(higgs_charge)
    M = matrix_from_terms(len(charges), terms)
    sums = q[:, None] + q[None, :]
    plus = M * (sums == h); minus = M * (sums == -h)
    assert np.array_equal(plus + minus, M)
    determinant = exact_determinant(M.tolist()); assert determinant != 0
    U = same_charge_unitary(charges, seed)
    transformed = U.T @ M @ U
    Qnew = U.conj().T @ Q @ U
    pnew, mnew = U.T@plus@U, U.T@minus@U
    singular = np.linalg.svd(M, compute_uv=False)
    alpha = .37
    D = np.diag(np.exp(1j*alpha*q))
    gauge_transformed_mass = np.exp(-1j*alpha*h)*plus + np.exp(1j*alpha*h)*minus
    errors = dict(symmetry=norm(M-M.T),
                  positive_sum_covariance=norm(Q.T@plus + plus@Q-h*plus),
                  negative_sum_covariance=norm(Q.T@minus + minus@Q+h*minus),
                  same_charge_unitarity=norm(U.conj().T@U-np.eye(len(q))),
                  same_charge_commutation=norm(Q@U-U@Q),
                  transformed_symmetry=norm(transformed-transformed.T),
                  transformed_positive_covariance=norm(Qnew.T@pnew+pnew@Qnew-h*pnew),
                  transformed_negative_covariance=norm(Qnew.T@mnew+mnew@Qnew+h*mnew),
                  Takagi_singular_values=norm(np.linalg.svd(transformed, compute_uv=False)-singular),
                  finite_gauge_covariance=norm(D.T@gauge_transformed_mass@D-M))
    assert max(errors.values()) < 3e-10
    return dict(higgs_charge=higgs_charge, mass_terms=[list(t) for t in terms],
                exact_determinant=str(determinant), numerical_rank=int(np.linalg.matrix_rank(M)),
                Takagi_singular_values=singular.tolist(), errors=errors,
                arbitrary_same_charge_unitary_covariance_is_analytic=True,
                arbitrary_coherent_mass_matrix_selection_is_analytic=True)


def template_case(number):
    assert number >= 8 and number % 4 == 0
    # k = 2a+3b; for even k use 8-blocks, for odd k use one 12-block.
    k = number//4; b = k % 2; a = (k-3*b)//2
    pairs = TEMPLATE_8*a + TEMPLATE_12*b
    charges = [q for pair in pairs for q in pair]
    assert len(charges) == number
    data = charge_data(charges)
    assert data["linear_anomaly"] == data["cubic_anomaly"] == "0"
    assert data["no_zero_charge"] and data["no_opposite_charges"]
    terms = [(2*j, 2*j+1, j+1) for j in range(number//2)]
    data.update(number=number, copies_of_8=a, copies_of_12=b,
                pair_sums=[sum(pair) for pair in pairs],
                mass=mass_calibration(charges, terms, 2, 20261022+number))
    return data


def odd_closed_walks():
    """Finite calibration of the loop-layer lemma, not its all-length proof."""
    records = []; total = 0
    for size in (1, 3, 5, 7, 9):
        count = 0
        for signs in itertools.product((-1, 1), repeat=size):
            q0 = sum((F(sign)*(-1)**j for j, sign in enumerate(signs)), F(0))/2
            q = [q0]
            for sign in signs[:-1]:
                q.append(F(sign)-q[-1])
            assert q[-1]+q0 == signs[-1]
            loop_positions = [i for i, x in enumerate(q) if abs(x) == F(1, 2)]
            assert loop_positions
            assert not (F(1, 2) in q and F(-1, 2) in q)
            pivot = loop_positions[0]
            remaining = [(pivot+j) % size for j in range(1, size)]
            pairs = [(remaining[j], remaining[j+1]) for j in range(0, size-1, 2)]
            assert all(q[i]+q[j] in (-1, 1) for i, j in pairs)
            assert sum((q[i]+q[j] for i, j in pairs), q[pivot]) == sum(q)
            count += 1
        total += count; records.append(dict(length=size, sign_words_checked=count))
    assert total == 682
    return dict(cases=records, total_odd_sign_words=total,
                every_walk_has_a_self_pairable_charge=True,
                both_loop_signs_in_one_orbit=False,
                all_lengths_proved_by_charge_orbit_graph_not_enumeration=True)


def counterexamples():
    q = [1, 1, 1, -4, -4, 5]
    terms = [(0, 1, 1), (2, 3, 2), (4, 5, 3)]
    powers = [-2, 3, -1]  # negative means H conjugate; h=1
    M = matrix_from_terms(6, terms)
    assert all(q[i]+q[j]+power == 0 for (i,j,_),power in zip(terms,powers))
    higher = dict(**charge_data(q), number=6, higgs_charge=1,
                  mass_terms=[list(t) for t in terms], signed_Higgs_powers=powers,
                  operator_dimensions=[3+abs(k) for k in powers],
                  exact_determinant=str(exact_determinant(M.tolist())),
                  numerical_rank=int(np.linalg.matrix_rank(M)), removed_assumption="dimension-four-only mass operators")

    q = [1, 5, -7, -8, 9]; terms = [(0, 4, 1), (1, 1, 2), (2, 3, 3)]
    signed_h = [-10, -10, 15]
    assert all(q[i]+q[j]+h == 0 for (i,j,_),h in zip(terms,signed_h))
    M = matrix_from_terms(5, terms)
    two = dict(**charge_data(q), number=5, higgs_charges=[10,15],
               mass_terms=[list(t) for t in terms], signed_Higgs_charges=signed_h,
               operator_dimensions=[4,4,4], exact_determinant=str(exact_determinant(M.tolist())),
               numerical_rank=int(np.linalg.matrix_rank(M)), removed_assumption="single Higgs")

    controls = []
    for name, q, h, terms in (
        ("mixed gravitational-U1 anomaly", [4,2,3,3,-1,-5], 6, [(0,1,1),(2,3,2),(4,5,3)]),
        ("cubic U1 anomaly", [4,-2,5,-7], 2, [(0,1,1),(2,3,2)]),
        ("complete chirality", [1,-1], 2, [(0,0,1),(1,1,2)]),
    ):
        controls.append(dict(**charge_data(q), number=len(q), removed_assumption=name,
                             no_bare_mass_term_used=True,
                             mass=mass_calibration(q,terms,h,20260000+len(q))))
    assert higher["linear_anomaly"] == higher["cubic_anomaly"] == "0"
    assert higher["numerical_rank"] == 6 and higher["no_opposite_charges"]
    assert two["linear_anomaly"] == two["cubic_anomaly"] == "0"
    assert two["numerical_rank"] == 5 and two["no_opposite_charges"]
    return dict(higher_dimension=higher, second_Higgs=two, other_controls=controls)


def equal_quadratic_inequivalent_spectra():
    A = [52,-26,65,-39,3,-29,50,-76]
    B = [52,-26,65,-39,43,-69,20,-46]
    terms = [(2*j,2*j+1,j+1) for j in range(4)]
    records = []
    for i, q in enumerate((A,B)):
        data = charge_data(q)
        assert data["linear_anomaly"] == data["cubic_anomaly"] == "0"
        assert data["quadratic_sum"] == "18252"
        assert data["no_zero_charge"] and data["no_opposite_charges"]
        assert data["primitive_integer_charge_gcd"] == 1
        data["mass"] = mass_calibration(q,terms,26,20269000+i)
        records.append(data)
    assert sorted(A) != sorted(B) and sorted(A) != sorted(-q for q in B)
    ar = [F(q,13) for q in A]; br = [F(q,13) for q in B]
    fourth = sum((q**4 for q in ar),F(0))-sum((q**4 for q in br),F(0))
    sixth = sum((abs(q)**6 for q in ar+br),F(0))
    x = F(1,10)
    fourth_signal = fourth*x**4/F(16*24)
    sixth_remainder = sixth*x**6/F(16*720)
    lower = fourth_signal-sixth_remainder
    assert lower > 0
    effects = [.5+sum(math.cos(float(q*x)) for q in qr)/16 for qr in (ar,br)]
    actual_gap = effects[0]-effects[1]
    assert 0 <= effects[0] <= 1 and 0 <= effects[1] <= 1
    assert actual_gap > float(lower)
    lipschitz = sum((abs(q) for q in ar+br),F(0))/16
    epsilon = F(1,10**6)
    packet_lower = lower-lipschitz*epsilon
    assert lipschitz == F(175,52) and packet_lower > 0
    return dict(candidates=records, same_higgs_charge=26, same_masses=[1,1,2,2,3,3,4,4],
                not_equivalent_by_permutation_or_U1_charge_inversion=True,
                same_quadratic_charge_sum=18252,
                charge_normalization="primitive integer compact-U1 charges; x=13*theta",
                holonomy_parameter_x=str(x), fourth_moment_difference=str(fourth),
                sixth_absolute_moment_sum=str(sixth),
                fourth_order_signal=str(fourth_signal), sixth_order_remainder_bound=str(sixth_remainder),
                strict_effect_gap_lower_exact=str(lower), strict_effect_gap_lower=float(lower),
                effect_probabilities=effects, numerical_effect_gap=actual_gap,
                difference_Lipschitz_bound_exact=str(lipschitz),
                smooth_packet_support_halfwidth_x=str(epsilon),
                packet_averaged_gap_lower_exact=str(packet_lower),
                packet_averaged_gap_lower=float(packet_lower),
                any_probability_distribution_in_support_satisfies_bound=True,
                compact_smooth_wavepacket_has_finite_electric_quadratic_form_by_analysis=True,
                packet_or_physical_instrument_constructed=False,
                shared_readout_is_representation_character_Wilson_effect=True,
                result_is_not_equal_all_physics_from_equal_quadratic_charge_sum=True)


def run():
    examples = [template_case(n) for n in (8,12,16,20,24,28,32)]
    union = [q for pair in TEMPLATE_8+TEMPLATE_12 for q in pair]
    assert all(q != 0 and -q not in union for q in union)
    walks = odd_closed_walks()
    controls = counterexamples()
    spectra = equal_quadratic_inequivalent_spectra()
    return dict(round=1022, date="2026-10-08", all_scientific_calibrations_passed=True,
                new_calibration_groups=1, cumulative_test_groups=3800, new_cognitive_axioms=0,
                scope="Four-dimensional Spin Weyl fields, one compact U(1), full chirality, one nonzero-charge Higgs VEV, dimension-four Yukawa-only full-rank symmetric mass, vanishing linear and cubic anomalies.",
                analytic_number_classification="N is a multiple of 4 and N >= 8; finite examples are calibrations only",
                all_mass_mixings_covered_by_analytic_cycle_cover_argument=True,
                arbitrary_N_existence_from_8_and_12_template_union=True,
                templates_have_no_cross_opposite_charges=True,
                charged_Higgs_VEV_assumed_not_derived=True,
                single_Higgs_renormalizable_mass_rule_is_extra_input=True,
                spatial_dimension_or_SM_group_selected=False,
                unique_chiral_charge_spectrum_selected=False,
                quantum_gravity_or_nonperturbative_UV_completion_proved=False,
                anomaly_free_effective_branch_not_full_FUCP_countermodel=True,
                actual_Wilson_instrument_or_preparation_built=False,
                mass_generation_does_not_claim_discrete_gauge_remnant_removed=True,
                existence_calibrations=examples, odd_cycle_lemma_calibration=walks,
                removed_assumption_controls=controls, residual_charge_freedom=spectra,
                historical_source_sha256={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
                                          for name in HISTORICAL})


def compare(fresh, saved, path="result"):
    if isinstance(fresh, dict):
        assert isinstance(saved, dict) and fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key],saved[key],path+"."+key)
    elif isinstance(fresh, list):
        assert isinstance(saved,list) and len(fresh)==len(saved),path
        for i,(a,b) in enumerate(zip(fresh,saved)):
            compare(a,b,path+f"[{i}]")
    elif isinstance(fresh,float):
        assert math.isclose(fresh,saved,rel_tol=4e-10,abs_tol=3e-11),(path,fresh,saved)
    else:
        assert fresh==saved,(path,fresh,saved)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true")
    args=parser.parse_args(); start=time.perf_counter(); result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    else:
        compare(result,json.loads(OUT.read_text("utf8")))
    print(json.dumps(dict(round=1022,passed=True,elapsed_seconds=time.perf_counter()-start,
        existence_cases=len(result["existence_calibrations"]),
        odd_cycle_words=result["odd_cycle_lemma_calibration"]["total_odd_sign_words"],
        Wilson_gap_lower=result["residual_charge_freedom"]["strict_effect_gap_lower_exact"],
        packet_gap_lower=result["residual_charge_freedom"]["packet_averaged_gap_lower_exact"],
        cumulative_test_groups=3800),ensure_ascii=False,indent=2))
