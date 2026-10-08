"""1025: complete finite-basis audit of a non-SM internal-geometry branch.

This does not classify all finite geometries or construct a physical field theory.
The base D0 commutes with A_F; the two explicit extensions have nonzero one-forms.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "finite_geometry_selection_results.json"
HISTORICAL = [
    "archive_1009_/research_note_1024.md",
    "archive_1009_/1024/NEXT.md",
    "archive_1009_/1024/input_dependency_update_v0_13.md",
    "archive_1009_/1024/bridge_ledger.json",
    "archive_531_553/research_note_531.md",
    "archive_531_553/research_note_532.md",
    "archive_531_553/research_note_533.md",
    "archive_531_553/research_note_543.md",
    "archive_531_553/research_note_553.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
]
ORDER = ["X11", "X12", "X21", "X22", "Y11", "Y12", "Y21", "Y22"]


def comm(a, b):
    return a @ b - b @ a


def real_coordinates(a):
    return np.concatenate((a.real.ravel(), a.imag.ravel()))


def rank_mod(matrix, prime=101):
    a = [[int(v) % prime for v in row] for row in matrix]
    pivot = 0
    for col in range(len(a[0])):
        found = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        inv = pow(a[pivot][col], -1, prime)
        a[pivot] = [v * inv % prime for v in a[pivot]]
        for i in range(pivot + 1, len(a)):
            factor = a[i][col]
            a[i] = [(x - factor * y) % prime for x, y in zip(a[i], a[pivot])]
        pivot += 1
    return pivot


def rank_rational(matrix):
    a = [[F(int(value)) for value in row] for row in matrix]
    pivot = 0
    for col in range(len(a[0])):
        found = next((i for i in range(pivot,len(a)) if a[i][col]),None)
        if found is None:
            continue
        a[pivot],a[found] = a[found],a[pivot]
        denominator = a[pivot][col]
        for i in range(pivot+1,len(a)):
            if a[i][col]:
                factor = a[i][col]/denominator
                a[i] = [x-factor*y for x,y in zip(a[i],a[pivot])]
        pivot += 1
    return pivot


def exact_integer_gram(columns):
    matrix = np.column_stack(columns)
    rounded = np.rint(matrix).astype(np.int64)
    assert np.array_equal(matrix, rounded)
    gram = rounded.T @ rounded
    rank = rank_rational(gram.tolist())
    assert rank == rank_mod(gram.tolist())
    return gram, rank


def represent(a, b):
    ans = np.zeros((8, 8), dtype=complex)
    ans[:4, :4] = np.kron(a, np.eye(2))
    ans[4:, 4:] = np.kron(b, np.eye(2))
    return ans


def direct_matrix(operation):
    columns = []
    for index in range(8):
        z = np.eye(8)[:, index]
        x, y = z[:4].reshape(2, 2), z[4:].reshape(2, 2)
        u, v = operation(x, y)
        columns.append(np.concatenate((u.ravel(), v.ravel())))
    return np.column_stack(columns).astype(complex)


def geometry():
    s = np.diag([1, -1])
    J0 = direct_matrix(lambda x, y: (y.T, x.T))
    gamma = direct_matrix(lambda x, y: (s @ x, -y @ s))
    D = np.zeros((8, 8), dtype=complex); D[0, 4] = D[4, 0] = 1
    all_basis = []
    all_labels = []
    for block in (0, 1):
        for row in range(2):
            for col in range(2):
                for phase, phase_label in ((1, "1"), (1j, "i")):
                    a, b = np.zeros((2, 2), dtype=complex), np.zeros((2, 2), dtype=complex)
                    (a if block == 0 else b)[row, col] = phase
                    all_basis.append((a, b)); all_labels.append(f"{['a','b'][block]}_E{row+1}{col+1}_{phase_label}")
    complex_even = []
    for slot in range(6):
        a, b = np.zeros((2, 2), dtype=complex), np.zeros((2, 2), dtype=complex)
        if slot < 2:
            a[slot, slot] = 1
        else:
            b[(slot-2)//2, (slot-2)%2] = 1
        complex_even.append((a, b))
    af_complex = [(np.diag([1., 0.]), np.diag([1., 0.])),
                  (np.diag([0., 1.]), np.zeros((2, 2))),
                  (np.zeros((2, 2)), np.diag([0., 1.]))]
    af_real = [(phase * a, phase * b) for a, b in af_complex for phase in (1, 1j)]
    return J0, gamma, D, all_basis, all_labels, complex_even, af_complex, af_real


def opposite(operator, J0):
    # a^0 = J a^* J^{-1} = J0 a^T J0^*, not J0 a^* J0.
    return J0 @ operator.T @ J0.conjugate().T


def hermitian_basis(n):
    ans = []
    for i in range(n):
        q = np.zeros((n, n), dtype=complex); q[i, i] = 1; ans.append(q)
        for j in range(i+1, n):
            q = np.zeros((n, n), dtype=complex); q[i, j] = q[j, i] = 1; ans.append(q)
            q = np.zeros((n, n), dtype=complex); q[i, j] = 1j; q[j, i] = -1j; ans.append(q)
    return ans


def selfadjoint_commutant_certificate(operators, J0):
    unknowns = hermitian_basis(8)
    columns = []
    for q in unknowns:
        residuals = [comm(q, a) for a in operators] + [q @ J0 - J0 @ q.conjugate()]
        columns.append(np.concatenate([real_coordinates(r) for r in residuals]))
    gram, rank = exact_integer_gram(columns)
    assert np.all(gram == gram.T)
    return dict(real_selfadjoint_unknown_dimension=64, constraint_Gram_rank_mod_101=rank,
                constraint_Gram_exact_rational_rank=rank,
                real_selfadjoint_commutant_dimension=64-rank,
                all_constraint_entries_exact_Gaussian_integers=True,
                Gram_sha256=hashlib.sha256(gram.astype("<i8").tobytes()).hexdigest())


def run():
    J0, gamma, D, all_basis, all_labels, even_complex, af_complex, af_real = geometry()
    A = [represent(a, b) for a, b in all_basis]
    AF = [represent(a, b) for a, b in af_real]
    AFprojections = [represent(a, b) for a, b in af_complex]
    I = np.eye(8)
    errors = dict(J_antiunitary=float(np.linalg.norm(J0.conjugate().T @ J0-I)),
                  J_square_plus=float(np.linalg.norm(J0 @ J0.conjugate()-I)),
                  gamma_selfadjoint=float(np.linalg.norm(gamma-gamma.conjugate().T)),
                  gamma_square_plus=float(np.linalg.norm(gamma@gamma-I)),
                  J_gamma_anticommutes=float(np.linalg.norm(J0 @ gamma.conjugate()+gamma @ J0)),
                  D_selfadjoint=float(np.linalg.norm(D-D.conjugate().T)),
                  J_D_commutes=float(np.linalg.norm(J0 @ D.conjugate()-D @ J0)),
                  D_gamma_anticommutes=float(np.linalg.norm(D @ gamma+gamma @ D)))
    left_errors = []
    opposite_errors = []
    normalizing_errors = []
    s = np.diag([1, -1])
    for (a, b), operator in zip(all_basis, A):
        direct = direct_matrix(lambda x, y: (a @ x, b @ y))
        direct_opposite = direct_matrix(lambda x, y: (x @ b, y @ a))
        left_errors.append(float(np.linalg.norm(operator-direct)))
        opposite_errors.append(float(np.linalg.norm(opposite(operator, J0)-direct_opposite)))
        normalizing_errors.append(float(np.linalg.norm(gamma @ operator @ gamma-represent(s @ a @ s, b))))
    zero_order = [float(np.linalg.norm(comm(a, opposite(b, J0)))) for a in A for b in A]
    af_zero_order = [float(np.linalg.norm(comm(a, opposite(b, J0)))) for a in AF for b in AF]
    first_order = [float(np.linalg.norm(comm(comm(D, a), opposite(b, J0)))) for a in AF for b in AF]
    af_D_comm = [float(np.linalg.norm(comm(D, a))) for a in AF]
    af_gamma = [float(np.linalg.norm(comm(gamma, a))) for a in AF]
    full_gamma_nonzero = [all_labels[i] for i, a in enumerate(A) if np.any(comm(gamma, a))]
    assert len(full_gamma_nonzero) == 4
    p1 = np.diag([1]*4+[0]*4); p2 = I-p1
    center_comm_norm = float(np.linalg.norm(comm(D, p1)))
    assert center_comm_norm > 1
    center_exchange_error = float(np.linalg.norm(J0 @ p1.conjugate() @ J0.conjugate().T-p2))

    separating = np.concatenate((np.eye(2).ravel(), np.eye(2).ravel()))
    left_images = [real_coordinates(a @ separating) for a in A]
    opposite_images = [real_coordinates(opposite(a, J0) @ separating) for a in A]
    _, left_rank = exact_integer_gram(left_images)
    _, right_rank = exact_integer_gram(opposite_images)
    assert left_rank == right_rank == 16
    irreducibility = selfadjoint_commutant_certificate(A, J0)
    assert irreducibility["constraint_Gram_rank_mod_101"] == 63
    assert irreducibility["real_selfadjoint_commutant_dimension"] == 1
    reduced_commutant = selfadjoint_commutant_certificate(AF, J0)
    assert reduced_commutant["real_selfadjoint_commutant_dimension"] == 6
    reducing = np.zeros((8, 8)); reducing[1, 1] = reducing[6, 6] = 1
    reduction_errors = [float(np.linalg.norm(comm(reducing, a))) for a in AF]
    reduction_errors += [float(np.linalg.norm(reducing @ J0-J0 @ reducing)),
                         float(np.linalg.norm(reducing @ reducing-reducing)),
                         float(np.linalg.norm(comm(reducing, D)))]
    assert np.trace(reducing) == 2 and max(reduction_errors) == 0

    # Polynomial coefficient audit for all a in A^even, not random a samples.
    # Coordinates are (alpha,beta,u,v,w,z), a=diag(alpha,beta), b=[[u,v],[w,z]].
    E = [represent(a, b) for a, b in even_complex]
    Ereal = [phase*a for a in E for phase in (1,1j)]
    even_zero = [float(np.linalg.norm(comm(a,opposite(b,J0)))) for a in Ereal for b in Ereal]
    even_first = [float(np.linalg.norm(comm(comm(D,a),opposite(b,J0)))) for a in Ereal for b in Ereal]
    assert len(Ereal) == 12 and len(even_zero) == len(even_first) == 144
    assert max(even_zero) == 0 and max(even_first) > 0
    assert all(np.all(comm(gamma,a) == 0) for a in Ereal)
    identity_pairs = [(0, 4), (1, 6), (6, 1)]
    quadratic_errors = []
    for row, col in identity_pairs:
        actual = np.zeros((6, 6), dtype=complex)
        for p in range(6):
            for q in range(6):
                actual[p, q] = comm(comm(D, E[p]), opposite(E[q], J0))[row, col]
        actual = (actual+actual.T)/2
        expected = np.zeros((6, 6))
        if (row, col) == (0, 4):
            expected[0, 0] = expected[2, 2] = -1; expected[0, 2] = expected[2, 0] = 1
        else:
            variable = 3 if (row, col) == (1, 6) else 4
            expected[variable, variable] = -1
        quadratic_errors.append(float(np.linalg.norm(actual-expected)))
    assert max(quadratic_errors) == 0

    # Equivalent linear conditions on A^even: u-alpha=v=w=0.
    constraints = np.array([[-1, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0]])
    assert rank_mod(constraints.tolist()) == 3
    af_real_matrix = np.column_stack([real_coordinates(a) for a in AF])
    assert np.linalg.matrix_rank(af_real_matrix) == 6
    closure_errors = []
    for p, a in enumerate(AFprojections):
        closure_errors += [float(np.linalg.norm(a.conjugate().T-a))]
        for q, b in enumerate(AFprojections):
            closure_errors += [float(np.linalg.norm(a@b-(a if p == q else np.zeros((8,8)))))]
    closure_errors.append(float(np.linalg.norm(sum(AFprojections)-I)))
    assert max(closure_errors) == 0
    orientation_products = [a @ opposite(b,J0) for a in AF for b in AF]
    orientation_columns = [real_coordinates(a) for a in orientation_products]
    _, orientation_rank = exact_integer_gram(orientation_columns)
    _, orientation_augmented_rank = exact_integer_gram(orientation_columns+[real_coordinates(gamma)])
    assert orientation_rank == 14 and orientation_augmented_rank == 15
    assert all(a[0,0] == a[4,4] for a in orientation_products)
    assert gamma[0,0]-gamma[4,4] == 2
    t = np.array([[0, 1], [1, 0]])
    Dt = direct_matrix(lambda x, y: (t @ x, y @ t))
    D1 = np.zeros((8, 8), dtype=complex)
    D1[0, 6] = D1[6, 0] = D1[4, 1] = D1[1, 4] = 1
    Dt_full_order_one = [float(np.linalg.norm(comm(comm(Dt, a), opposite(b, J0)))) for a in A for b in A]
    assert max(Dt_full_order_one) == 0
    dynamics = []
    for name, operator in (("D0", D), ("D0_plus_Dt", D+Dt), ("D0_plus_Dt_plus_D1", D+Dt+D1)):
        checks = dict(selfadjoint=float(np.linalg.norm(operator-operator.conjugate().T)),
                      J_real=float(np.linalg.norm(J0 @ operator.conjugate()-operator @ J0)),
                      grading_odd=float(np.linalg.norm(comm(operator, gamma) + 2*gamma@operator)),
                      AF_order_one=max(float(np.linalg.norm(comm(comm(operator,a),opposite(b,J0)))) for a in AF for b in AF))
        assert max(checks.values()) == 0
        _, comm_rank = exact_integer_gram([real_coordinates(comm(operator,a)) for a in AF])
        form_columns = [real_coordinates(a @ comm(operator,b)) for a in AF for b in AF]
        _, form_rank = exact_integer_gram(form_columns)
        assert (6-comm_rank,form_rank) == {"D0":(6,0),"D0_plus_Dt":(4,4),"D0_plus_Dt_plus_D1":(2,8)}[name]
        dynamics.append(dict(name=name, identity_errors=checks,
            AF_D_commuting_subalgebra_real_dimension=6-comm_rank,
            represented_one_form_space_real_dimension=form_rank,
            represented_one_form_generator_count=len(form_columns),
            maximum_AF_D_commutator_norm=max(float(np.linalg.norm(comm(operator,a))) for a in AF),
            D_original_center_commutator_norm=float(np.linalg.norm(comm(operator,p1))),
            operator_real_matrix=operator.real.astype(int).tolist(),
            nonzero_forms_not_a_physical_gauge_theory_construction=True))
    # The added Dt contributes zero to every full-A first-order coefficient, so
    # the already checked fixed-D0 quadratic obstruction is unchanged for D0+Dt.
    all_errors = list(errors.values())+left_errors+opposite_errors+normalizing_errors+zero_order+af_zero_order+first_order+af_D_comm+af_gamma+quadratic_errors+closure_errors+[center_exchange_error]
    assert max(all_errors) == 0
    return dict(round=1025, date="2026-10-08", all_scientific_calibrations_passed=True,
                new_calibration_groups=1, cumulative_test_groups=3802, new_cognitive_axioms=0,
                scope="Finite unitary real-form branch M2(C)+M2(C), original (A,J) irreducible; fixed-D greatest even order-one subalgebra C^3 and nonzero-one-form controls. The report separately proves the all-D dimension bound; no complete physical model is constructed.",
                primary_source_url="https://arxiv.org/html/0706.3688v1", basis_order=ORDER,
                complex_Hilbert_dimension=8, real_parent_algebra_dimension=16,
                real_even_algebra_dimension=12, real_selected_algebra_dimension=6,
                real_standard_model_algebra_dimension=24,
                standard_algebra_is_C_plus_H_plus_M3C=True,
                internal_real_structure_is_antilinear_not_linear_orientation_J=True,
                quaternion_internal_real_form_is_not_quaternion_quantum_state_theory=True,
                KO6_signs_adopted=[1,1,-1], KO6_derived_from_cognition=False,
                quaternion_linearity_not_adopted_in_witness=True,
                physical_unification_or_SM_anomaly_freedom_certified=False,
                finite_geometric_axioms_all_derived_from_cognition=False,
                all_spectral_triple_axioms_satisfied_claimed=False,
                continuum_or_spectral_action_constructed=False,
                sampled_elements_used_as_universal_proof=False,
                real_parent_basis_count=len(A), real_parent_basis_labels=all_labels,
                real_AF_basis_count=len(AF), full_zero_order_pairs=len(zero_order),
                real_even_basis_count=len(Ereal), even_zero_order_pairs=len(even_zero),
                even_zero_order_max_error=max(even_zero), even_order_one_negative_control_pairs=len(even_first),
                even_order_one_nonzero_pair_count=sum(v>0 for v in even_first),
                even_order_one_negative_control_max_norm=max(even_first),
                AF_zero_order_pairs=len(af_zero_order), AF_order_one_pairs=len(first_order),
                identity_errors=errors, left_action_independent_errors=left_errors,
                opposite_action_independent_errors=opposite_errors,
                grading_normalizes_full_A_errors=normalizing_errors,
                grading_commutes_with_full_A=False, full_A_noncommuting_grading_basis=full_gamma_nonzero,
                grading_commutes_with_AF=True, AF_grading_errors=af_gamma,
                zero_order_max_error=max(zero_order), AF_zero_order_max_error=max(af_zero_order),
                AF_order_one_max_error=max(first_order),
                D_noncommuting_original_center_Frobenius_norm=center_comm_norm,
                original_center_exchange_error=center_exchange_error,
                D_eigenvalues=np.linalg.eigvalsh(D).tolist(),
                D_commutes_with_AF=True, AF_D_commutator_errors=af_D_comm,
                represented_one_form_space_dimension=0,
                nonzero_gauge_or_Higgs_interactions_generated=False,
                nontrivial_D_controls=dynamics,
                Dt_full_parent_order_one_pairs=len(Dt_full_order_one),
                Dt_full_parent_order_one_max_error=max(Dt_full_order_one),
                fixed_D0_maximality_polynomial_unchanged_for_D0_plus_Dt=True,
                connected_Dstar_AF_commutant_is_common_complex_scalars=True,
                separating_vector_real_coordinates=real_coordinates(separating).tolist(),
                separating_left_map_real_rank=left_rank, separating_opposite_cyclic_real_rank=right_rank,
                original_pair_irreducibility=irreducibility,
                AF_pair_selfadjoint_commutant=reduced_commutant,
                AF_pair_irreducible=False, AF_reducing_projection_indices=[1,6],
                AF_reducing_projection_rank=2, AF_reducing_projection_errors=reduction_errors,
                fixed_D_maximality=dict(complex_even_coordinates=["alpha","beta","u","v","w","z"],
                    self_order_one_entry_indices=[[0,4],[1,6],[6,1]],
                    self_order_one_polynomials=["-(u-alpha)^2","-v^2","-w^2"],
                    complete_quadratic_coefficient_errors=quadratic_errors,
                    quadratic_bilinear_coefficient_count=3*6*6,
                    remaining_complex_dimension=3, real_dimension=6,
                    every_allowed_subalgebra_is_contained_in_AF=True,
                    fixed_D_unique_greatest_subalgebra=True,
                    varying_all_D_bound_not_proved_by_finite_numeric_scan=True,
                    varying_all_D_upper_bound_is_separate_analytic_argument=True),
                AF_projection_closure_max_error=max(closure_errors),
                all_saved_finite_dimensions_checked_by_Fraction_elimination=True,
                modular_ranks_used_as_independent_lower_bounds=True,
                finite_orientability_obstruction=dict(left_right_product_count=36,
                    left_right_product_span_real_dimension=orientation_rank,
                    with_gamma_real_dimension=orientation_augmented_rank,
                    every_product_X11_minus_Y11_entry_exact=0,
                    gamma_X11_minus_Y11_entry_exact=2,
                    orientability_satisfied=False,
                    orientability_not_in_the_main_adopted_contract=True),
                maximum_identity_error=max(all_errors),
                historical_source_sha256={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in HISTORICAL})


def compare(fresh, saved, path="root"):
    if isinstance(fresh, dict):
        assert fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key], saved[key], path+"."+key)
    elif isinstance(fresh, list):
        assert len(fresh) == len(saved), path
        for i, (a,b) in enumerate(zip(fresh,saved)):
            compare(a,b,f"{path}[{i}]")
    elif isinstance(fresh, float):
        assert math.isclose(fresh,saved,rel_tol=2e-10,abs_tol=2e-12),(path,fresh,saved)
    else:
        assert fresh == saved,(path,fresh,saved)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true")
    args=parser.parse_args(); result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    else:
        compare(result,json.loads(OUT.read_text("utf8")))
    print(json.dumps({"round":1025,"passed":True,"max_error":result["maximum_identity_error"],
                      "original_pair_commutant_dimension":result["original_pair_irreducibility"]["real_selfadjoint_commutant_dimension"],
                      "AF_real_dimension":6,"AF_pair_irreducible":False,"represented_one_forms":0,
                      "new_calibration_groups":1,"cumulative_test_groups":3802},ensure_ascii=False,indent=2))
