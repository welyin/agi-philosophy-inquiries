"""Limited round-1025 delivery verification; live navigation is excluded."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import finite_geometry_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / "research_note_1025.md"
OUT = HERE / "research_round_1025_checks.json"
OWN = ["finite_geometry_selection.py", "finite_geometry_selection_results.json", "review.md",
       "selection_audit.md", "input_dependency_update_v0_14.md", "NEXT.md", "verify_round1025.py"]
HISTORICAL = {"archive_1009_/research_note_1024.md", "archive_1009_/1024/NEXT.md",
              "archive_1009_/1024/input_dependency_update_v0_13.md", "archive_1009_/1024/bridge_ledger.json",
              "archive_1009_/1009/input_dependency_ledger_v0_1.md"} | {
              f"archive_531_553/research_note_{r}.md" for r in (531, 532, 533, 543, 553)}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_mod103(a):
    a = np.array(a, dtype=np.int64).copy() % 103
    r = 0
    for col in range(a.shape[1]):
        nonzero = np.flatnonzero(a[r:, col])
        if not len(nonzero):
            continue
        k = r + nonzero[0]
        a[[r, k]] = a[[k, r]]
        a[r] = a[r] * pow(int(a[r, col]), -1, 103) % 103
        for i in range(r+1, len(a)):
            a[i] = (a[i] - a[i, col]*a[r]) % 103
        r += 1
    return r


def exact_rank(matrix):
    rows = [[Fraction(int(value)) for value in row] for row in matrix]
    pivots = 0
    for column in range(len(rows[0])):
        candidates = [i for i in range(pivots,len(rows)) if rows[i][column]]
        if not candidates:
            continue
        selected = candidates[0]
        rows[pivots],rows[selected] = rows[selected],rows[pivots]
        scale = rows[pivots][column]
        rows[pivots] = [v/scale for v in rows[pivots]]
        for i in range(pivots+1,len(rows)):
            coefficient = rows[i][column]
            if coefficient:
                rows[i] = [u-coefficient*v for u,v in zip(rows[i],rows[pivots])]
        pivots += 1
    assert pivots == rank_mod103(matrix)
    return pivots


def independent_parent_basis():
    answer = []
    for block in (0, 1):
        for i in range(2):
            for j in range(2):
                for phase in (1, 1j):
                    a = np.zeros((8, 8), dtype=complex)
                    for column in range(2):
                        a[4*block+2*i+column, 4*block+2*j+column] = phase
                    answer.append(a)
    return answer


def realvec(a):
    return np.r_[a.real.ravel(), a.imag.ravel()]


def commutant_rank(operators, J0):
    hermitian = []
    for i in range(8):
        q = np.zeros((8, 8), complex); q[i, i] = 1; hermitian.append(q)
        for j in range(i+1, 8):
            for phase in (1, 1j):
                q = np.zeros((8, 8), complex); q[i, j] = phase; q[j, i] = np.conjugate(phase)
                hermitian.append(q)
    columns = []
    for q in hermitian:
        values = [realvec(q@a-a@q) for a in operators]
        values += [realvec(q@J0-J0@q.conjugate())]
        columns.append(np.concatenate(values))
    constraints = np.column_stack(columns)
    assert np.array_equal(constraints, constraints.astype(np.int64))
    return exact_rank(constraints.astype(np.int64).T @ constraints.astype(np.int64))


def local_links(text):
    for match in LINK.finditer(EXCLUDED.sub(lambda match: " "*len(match[0]), text)):
        stripped = (match[1] if match[1] is not None else match[2]).strip()
        target = stripped[1:stripped.index(">")] if stripped.startswith("<") and ">" in stripped else re.split(r"\s+[\"']", stripped, 1)[0]
        local = unquote(target.split("#", 1)[0])
        if local and not re.match(r"^[a-zA-Z]+:", local) and not local.startswith("//"):
            yield target, local


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text("utf8")))
    assert fresh["round"] == 1025 and fresh["all_scientific_calibrations_passed"]
    assert fresh["new_calibration_groups"] == 1 and fresh["cumulative_test_groups"] == 3802
    assert fresh["new_cognitive_axioms"] == 0
    assert fresh["basis_order"] == ["X11","X12","X21","X22","Y11","Y12","Y21","Y22"]
    assert [fresh[k] for k in ("complex_Hilbert_dimension", "real_parent_algebra_dimension", "real_even_algebra_dimension", "real_selected_algebra_dimension", "real_standard_model_algebra_dimension")] == [8,16,12,6,24]
    assert fresh["KO6_signs_adopted"] == [1,1,-1]
    for key in ("KO6_derived_from_cognition", "physical_unification_or_SM_anomaly_freedom_certified",
                "finite_geometric_axioms_all_derived_from_cognition", "all_spectral_triple_axioms_satisfied_claimed",
                "continuum_or_spectral_action_constructed", "sampled_elements_used_as_universal_proof",
                "nonzero_gauge_or_Higgs_interactions_generated"):
        assert fresh[key] is False, key
    assert fresh["internal_real_structure_is_antilinear_not_linear_orientation_J"]
    assert fresh["quaternion_internal_real_form_is_not_quaternion_quantum_state_theory"]
    assert fresh["quaternion_linearity_not_adopted_in_witness"]

    A = independent_parent_basis()
    permutation = [4,6,5,7,0,2,1,3]
    J0 = np.eye(8)[permutation]
    gamma = np.diag([1,1,-1,-1,-1,1,-1,1])
    D0 = np.zeros((8,8)); D0[0,4] = D0[4,0] = 1
    af_diagonals = [np.array(v,complex) for v in ([1,1,0,0,1,1,0,0], [0,0,1,1,0,0,0,0], [0,0,0,0,0,0,1,1])]
    af_diagonals = [phase*v for v in af_diagonals for phase in (1,1j)]
    AF = [np.diag(v) for v in af_diagonals]
    assert np.array_equal(J0@J0, np.eye(8)) and np.array_equal(J0@gamma, -gamma@J0)
    assert np.array_equal(J0@D0, D0@J0) and np.array_equal(gamma@D0, -D0@gamma)
    assert fresh["full_zero_order_pairs"] == 256 and fresh["AF_zero_order_pairs"] == fresh["AF_order_one_pairs"] == 36
    for a in A:
        for b in A:
            op = J0@b.T@J0
            assert np.all(a@op == op@a)
    assert len(fresh["full_A_noncommuting_grading_basis"]) == 4
    assert sum(bool(np.any(gamma@a-a@gamma)) for a in A) == 4
    even = [a for a in A if np.all(gamma@a == a@gamma)]
    assert len(even) == fresh["real_even_basis_count"] == 12
    even_residuals = []
    for a in even:
        for b in even:
            opposite = J0@b.T@J0
            first = D0@a-a@D0
            even_residuals.append(float(np.linalg.norm(first@opposite-opposite@first)))
    assert len(even_residuals) == fresh["even_order_one_negative_control_pairs"] == 144
    assert sum(v>0 for v in even_residuals) == fresh["even_order_one_nonzero_pair_count"] == 56
    assert math.isclose(max(even_residuals), fresh["even_order_one_negative_control_max_norm"])
    assert fresh["even_zero_order_pairs"] == 144 and fresh["even_zero_order_max_error"] == 0
    assert not fresh["grading_commutes_with_full_A"] and fresh["grading_commutes_with_AF"]
    assert all(np.array_equal(gamma@a, a@gamma) for a in AF)
    assert all(np.array_equal(D0@a,a@D0) for a in AF)
    assert fresh["D_commutes_with_AF"] and fresh["represented_one_form_space_dimension"] == 0
    p1 = np.diag([1]*4+[0]*4)
    assert np.linalg.matrix_rank(D0@p1-p1@D0) == 2
    assert math.isclose(fresh["D_noncommuting_original_center_Frobenius_norm"], math.sqrt(2))
    assert fresh["D_eigenvalues"] == [-1.,0.,0.,0.,0.,0.,0.,1.]
    assert fresh["separating_left_map_real_rank"] == fresh["separating_opposite_cyclic_real_rank"] == 16
    xi = np.array([1,0,0,1,1,0,0,1])
    left = np.column_stack([realvec(a@xi) for a in A])
    right = np.column_stack([realvec(J0@a.T@J0@xi) for a in A])
    assert rank_mod103(left) == rank_mod103(right) == 16
    assert commutant_rank(A, J0) == fresh["original_pair_irreducibility"]["constraint_Gram_rank_mod_101"] == 63
    assert commutant_rank(AF, J0) == fresh["AF_pair_selfadjoint_commutant"]["constraint_Gram_rank_mod_101"] == 58
    assert fresh["original_pair_irreducibility"]["real_selfadjoint_commutant_dimension"] == 1
    assert fresh["AF_pair_selfadjoint_commutant"]["real_selfadjoint_commutant_dimension"] == 6
    reducing = np.diag([0,1,0,0,0,0,1,0])
    assert np.trace(reducing) == fresh["AF_reducing_projection_rank"] == 2
    assert np.all(reducing@J0 == J0@reducing)
    assert all(np.all(reducing@a == a@reducing) for a in AF)
    assert not fresh["AF_pair_irreducible"]

    maximal = fresh["fixed_D_maximality"]
    assert maximal["self_order_one_polynomials"] == ["-(u-alpha)^2","-v^2","-w^2"]
    assert maximal["quadratic_bilinear_coefficient_count"] == 108
    assert maximal["complete_quadratic_coefficient_errors"] == [0.,0.,0.]
    assert maximal["remaining_complex_dimension"] == 3 and maximal["real_dimension"] == 6
    assert maximal["every_allowed_subalgebra_is_contained_in_AF"] and maximal["fixed_D_unique_greatest_subalgebra"]
    # Nonzero excluded elements independently demonstrate that all-zero residuals
    # on AF are not a disabled first-order test.
    excluded = [(np.diag([0,0]), np.diag([1,0]), (0,4)),
                (np.zeros((2,2)), np.array([[0,1],[0,0]]), (1,6)),
                (np.zeros((2,2)), np.array([[0,0],[1,0]]), (6,1))]
    obstruction_values = []
    for a,b,(i,j) in excluded:
        op = np.zeros((8,8)); op[:4,:4] = np.kron(a,np.eye(2)); op[4:,4:] = np.kron(b,np.eye(2))
        first = D0@op-op@D0; opposite = J0@op.T@J0
        obstruction = first@opposite-opposite@first
        assert obstruction[i,j] == -1
        obstruction_values.append(int(obstruction[i,j]))

    cases = fresh["nontrivial_D_controls"]
    assert len(cases) == 3
    dimensions = []
    form_dimensions = []
    for row in cases:
        D = np.array(row["operator_real_matrix"],float)
        assert np.array_equal(D,D.T) and np.array_equal(D@J0,J0@D)
        assert np.array_equal(D@gamma,-gamma@D)
        for a in af_diagonals:
            for b in af_diagonals:
                opposite = b[permutation]
                coefficient = D*(a[None,:]-a[:,None])*(opposite[None,:]-opposite[:,None])
                assert np.all(coefficient == 0)
        columns = np.column_stack([realvec(D@a-a@D) for a in AF])
        dim = 6-exact_rank(columns.astype(np.int64).T@columns.astype(np.int64))
        forms = np.column_stack([realvec(a@(D@b-b@D)) for a in AF for b in AF])
        form_rank = exact_rank(forms.astype(np.int64).T@forms.astype(np.int64))
        assert dim == row["AF_D_commuting_subalgebra_real_dimension"]
        assert form_rank == row["represented_one_form_space_real_dimension"]
        assert max(row["identity_errors"].values()) == 0
        dimensions.append(dim); form_dimensions.append(form_rank)
    assert dimensions == [6,4,2] and form_dimensions == [0,4,8]
    assert cases[-1]["maximum_AF_D_commutator_norm"] > 0
    assert fresh["Dt_full_parent_order_one_pairs"] == 256 and fresh["Dt_full_parent_order_one_max_error"] == 0
    assert fresh["fixed_D0_maximality_polynomial_unchanged_for_D0_plus_Dt"]
    assert fresh["connected_Dstar_AF_commutant_is_common_complex_scalars"]
    assert fresh["all_saved_finite_dimensions_checked_by_Fraction_elimination"]
    orientation = fresh["finite_orientability_obstruction"]
    products = [a@J0@b.T@J0 for a in AF for b in AF]
    assert len(products) == orientation["left_right_product_count"] == 36
    columns = np.column_stack([realvec(a) for a in products]).astype(np.int64)
    extended = np.column_stack([columns, realvec(gamma)]).astype(np.int64)
    assert exact_rank(columns.T@columns) == orientation["left_right_product_span_real_dimension"] == 14
    assert exact_rank(extended.T@extended) == orientation["with_gamma_real_dimension"] == 15
    assert all(a[0,0]-a[4,4] == 0 for a in products)
    assert gamma[0,0]-gamma[4,4] == orientation["gamma_X11_minus_Y11_entry_exact"] == 2
    assert not orientation["orientability_satisfied"] and orientation["orientability_not_in_the_main_adopted_contract"]
    assert fresh["maximum_identity_error"] == 0
    assert max(fresh["identity_errors"].values()) == 0
    assert set(fresh["historical_source_sha256"]) == HISTORICAL
    for name,digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name) == digest, name

    assets = [NOTE] + [HERE/name for name in OWN]
    assert len(assets) == 8
    assert all(p.name not in {"README.md","research_direction.md","RESEARCH_STATE.md"} for p in assets)
    links = 0
    for path in assets:
        assert path.is_file(),path
        if path.suffix == ".md":
            text = path.read_text("utf-8-sig"); assert text.count("$$") % 2 == 0,path
            for target,local in local_links(text):
                resolved = (path.parent/local.replace("\\","/")).resolve()
                assert resolved.exists() or (prospective and resolved == OUT),(path,target)
                links += 1
    assert all(f"## {n}." in NOTE.read_text("utf8") for n in range(1,11))
    assert "独立科学签审通过" in (HERE/"review.md").read_text("utf8")
    return dict(round=1025,date="2026-10-08",all_delivery_checks_passed=True,
                scientific_result_reproduced=True,new_calibration_groups=1,cumulative_research_groups=3802,
                calibration_kind="finite algebra only; no new physical prediction or experiment",
                new_cognitive_axioms=0,local_links_checked=links,frozen_current_files=8,historical_input_files=len(HISTORICAL),
                live_navigation_frozen=False,neighboring_round_frozen=False,
                Hilbert_complex_dimension=8,parent_real_algebra_dimension=16,AF_real_dimension=6,
                full_order_zero_pairs=256,AF_order_one_pairs_per_D=36,Dirac_controls=3,
                even_basis_order_one_negative_controls=144,even_nonzero_order_one_controls=56,
                polynomial_coefficient_checks=108,excluded_order_one_entries=obstruction_values,
                independent_commutant_Gram_prime=103,original_pair_commutant_dimension=1,
                rational_elimination_and_modular_lower_bounds_both_checked=True,
                AF_pair_commutant_dimension=6,AF_D_commuting_dimensions=dimensions,
                represented_one_form_real_dimensions=form_dimensions,
                initial_A_not_even_commuting_with_grading=True,AF_not_claimed_irreducible=True,
                fixed_D_greatest_subalgebra_certificate=True,
                orientability_span_real_dimension=14,with_grading_real_dimension=15,
                all_D_dimension_bound_requires_report_analytic_argument=True,
                all_spectral_axioms_or_physical_model_certified=False,
                maximum_exact_matrix_identity_error=0.,goal_completed=False,visual_checks_performed=False,
                historical_source_sha256=fresh["historical_source_sha256"],
                source_sha256={str(path.relative_to(ROOT)).replace("\\","/"):sha(path) for path in assets})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--write", "--write-receipt",dest="write",action="store_true")
    args = parser.parse_args(); result = verify(prospective=args.write)
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")},ensure_ascii=False,indent=2))
