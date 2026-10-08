"""Round 1024: finite certificates for a typed parent-object audit.

No new physical model, observable prediction, or implementation is certified.
Default execution only recomputes and compares the saved results.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "parent_bridge_audit_results.json"
HISTORICAL = [f"archive_1009_/research_note_{n}.md" for n in range(1009, 1024)] + [
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_990_1008/993/common_candidate_v1.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
]


def rref(rows, columns=None):
    a = [[F(x) for x in row] for row in rows]
    width = len(a[0]) if a else columns
    pivots = []
    for col in range(width):
        found = next((i for i in range(len(pivots), len(a)) if a[i][col]), None)
        if found is None:
            continue
        p = len(pivots)
        a[p], a[found] = a[found], a[p]
        factor = a[p][col]
        a[p] = [v / factor for v in a[p]]
        for i in range(len(a)):
            if i != p and a[i][col]:
                factor = a[i][col]
                a[i] = [u - factor * v for u, v in zip(a[i], a[p])]
        pivots.append(col)
    return a, pivots


def nullspace(rows, columns):
    reduced, pivots = rref(rows, columns)
    basis = []
    for free in range(columns):
        if free in pivots:
            continue
        vector = [F(0)] * columns
        vector[free] = F(1)
        for i, pivot in enumerate(pivots):
            vector[pivot] = -reduced[i][free]
        basis.append(vector)
    return pivots, basis


def exact_rank(matrix):
    return len(rref(matrix)[1])


def field_certificate():
    labels = [f"Qup_{c}" for c in range(3)] + [f"Qdown_{c}" for c in range(3)]
    labels += [f"uC_{c}" for c in range(3)] + [f"dC_{c}" for c in range(3)] + ["N", "D", "E"]
    per_y = [1] * 6 + [-4] * 3 + [2] * 3 + [-3] * 2 + [6]
    per_em = [2] * 3 + [-1] * 3 + [-2] * 3 + [1] * 3 + [0, -3, 3]
    y = np.array(per_y * 3, dtype=int)
    em = np.array(per_em * 3, dtype=int)
    n = len(y)
    h = 3
    adjacency = {i: [j for j in range(n) if abs(int(y[i] + y[j])) == h] for i in range(n)}
    unseen = set(range(n))
    components = []
    while unseen:
        seed = min(unseen)
        color = {seed: 0}
        pending = [seed]
        while pending:
            i = pending.pop()
            for j in adjacency[i]:
                if j not in color:
                    color[j] = 1 - color[i]
                    pending.append(j)
                assert color[i] != color[j]
        unseen.difference_update(color)
        parts = [[i for i in color if color[i] == c] for c in (0, 1)]
        sizes = sorted(map(len, parts))
        components.append(dict(charge_levels=sorted({int(y[i]) for i in color}),
                               partition_sizes=sizes, symmetric_mass_rank_upper_bound=2 * min(sizes)))

    one_h = np.zeros((n, n), dtype=int)
    one_h_conjugate = np.zeros((n, n), dtype=int)
    pair_list = []
    for generation in range(3):
        o = 15 * generation
        pairs = [(o + c, o + 6 + c, 3 * generation + 1) for c in range(3)]
        pairs += [(o + 3 + c, o + 9 + c, 3 * generation + 2) for c in range(3)]
        pairs += [(o + 13, o + 14, 3 * generation + 3)]
        for i, j, mass in pairs:
            target = one_h if y[i] + y[j] == -h else one_h_conjugate
            target[i, j] = target[j, i] = mass
            pair_list.append(dict(indices=[i, j], mass=mass, hypercharge_sum=int(y[i] + y[j])))
    charged = one_h + one_h_conjugate
    two_h = np.zeros((n, n), dtype=int)
    for generation in range(3):
        two_h[15 * generation + 12, 15 * generation + 12] = generation + 1
    complete = charged + two_h
    color_generators = []
    for i in range(3):
        for j in range(i + 1, 3):
            symmetric = np.zeros((3, 3), dtype=complex)
            symmetric[i, j] = symmetric[j, i] = 1
            antisymmetric = np.zeros((3, 3), dtype=complex)
            antisymmetric[i, j] = -1j; antisymmetric[j, i] = 1j
            color_generators.extend([symmetric, antisymmetric])
    color_generators += [np.diag([1., -1., 0.]), np.diag([1., 1., -2.])]
    color_ward_errors = []
    for generator in color_generators:
        T = np.zeros((n, n), dtype=complex)
        for generation in range(3):
            o = 15 * generation
            for offset in (0, 3):
                T[o+offset:o+offset+3, o+offset:o+offset+3] = generator
            for offset in (6, 9):
                T[o+offset:o+offset+3, o+offset:o+offset+3] = -generator.conjugate()
        color_ward_errors.append(float(np.linalg.norm(T.T @ complete + complete @ T)))
    assert len(color_ward_errors) == 8 and max(color_ward_errors) == 0
    qy = np.diag(y)
    qe = np.diag(em)
    identities = dict(
        H_coefficient=float(np.linalg.norm(qy @ one_h + one_h @ qy + h * one_h)),
        H_conjugate_coefficient=float(np.linalg.norm(qy @ one_h_conjugate + one_h_conjugate @ qy - h * one_h_conjugate)),
        HH_coefficient=float(np.linalg.norm(qy @ two_h + two_h @ qy + 2 * h * two_h)),
        unbroken_electric_mass=float(np.linalg.norm(qe @ complete + complete @ qe)),
    )
    covariance = []
    for alpha in (.13, -.29):
        transformed_h = np.exp(1j * h * alpha)
        transformed_mass = transformed_h * one_h + transformed_h.conjugate() * one_h_conjugate + transformed_h**2 * two_h
        field_phase = np.diag(np.exp(1j * y * alpha))
        covariance.append(float(np.linalg.norm(field_phase.T @ transformed_mass @ field_phase - complete)))
    weak_t3 = [F(1, 2)] * 3 + [-F(1, 2)] * 3 + [F(0)] * 6 + [F(1, 2), -F(1, 2), F(0)]
    assert [3 * t + F(q, 2) for t, q in zip(weak_t3, per_y)] == list(map(F, per_em))
    bound = sum(row["symmetric_mass_rank_upper_bound"] for row in components)
    rank_charged = exact_rank(charged.tolist())
    rank_complete = exact_rank(complete.tolist())
    assert bound == rank_charged == 42 and rank_complete == 45
    assert sum(map(int, y)) == sum(int(q)**3 for q in y) == 0
    assert all(q and -q not in y for q in y)
    assert max(identities.values()) == 0 and max(covariance) < 1e-11
    return dict(
        number_of_generations=3, number_of_left_handed_Weyl_fields=n,
        per_generation_labels=labels, per_generation_6Y=per_y, three_generation_6Y=y.tolist(),
        hypercharge_multiplicities={str(k): v for k, v in sorted(Counter(y.tolist()).items())},
        linear_anomaly_exact=str(sum(map(int, y))), cubic_anomaly_exact=str(sum(int(q)**3 for q in y)),
        fully_chiral_under_hypercharge=True, Higgs_6Y=h, dimension_four_graph_components=components,
        dimension_four_maximum_rank=bound, charged_matrix_exact_rank=rank_charged,
        charged_pairs=pair_list, H_mass_matrix=one_h.tolist(), H_conjugate_mass_matrix=one_h_conjugate.tolist(),
        color_ward_generator_count=8, color_ward_identity_errors=color_ward_errors,
        per_generation_up_down_electron_masses=[[3*g+1, 3*g+2, 3*g+3] for g in range(3)],
        each_quark_mass_shared_across_three_colors=True,
        HH_Weinberg_mass_matrix=two_h.tolist(), Weinberg_extended_exact_rank=rank_complete,
        Weinberg_neutral_slots=[12, 27, 42], polynomial_charge_identity_errors=identities,
        finite_phase_covariance_errors=covariance,
        per_generation_3Qem=per_em, electric_multiplicities={str(k): v for k, v in sorted(Counter(em.tolist()).items())},
        electric_zero_charge_multiplicity=int(np.count_nonzero(em == 0)), electric_opposite_charge_pairs=[1, 2, 3],
        fully_chiral_under_electric_charge=False, electric_Higgs_vev_charge=0,
        two_distinct_U1_generators_not_interchangeable=True,
        round1022_full_rank_dimension_four_hypothesis_failed=True,
        Weinberg_extension_changes_operator_dimension=True,
        anomaly_check_is_not_all_parent_consistency_conditions=True,
    )


def polynomial_hessian(kappa):
    """Differentiate the entire polynomial exactly; no sampled 'for all' claim."""
    poly = {}
    def add(exponent, value):
        key = tuple(exponent)
        poly[key] = poly.get(key, F(0)) + F(value)
    zero = [0] * 5
    add(zero, F(1, 4))
    for i in range(4):
        power = [0] * 5; power[i] = 2; add(power, -F(1, 2))
        power[i] = 4; add(power, F(1, 4))
        for j in range(i + 1, 4):
            power = [0] * 5; power[i] = power[j] = 2; add(power, F(1, 2))
        power = [0] * 5; power[i] = power[4] = 2; add(power, F(kappa, 4))
    power = [0] * 5; power[4] = 2; add(power, F(1, 2))
    power[4] = 4; add(power, F(1, 4))
    coefficients = {}
    for exponent, coefficient in poly.items():
        if not coefficient:
            continue
        for i in range(5):
            for j in range(5):
                power = list(exponent)
                multiplier = power[i]
                power[i] -= 1
                multiplier *= power[j]
                power[j] -= 1
                if not multiplier or min(power) < 0:
                    continue
                key = tuple(power)
                if key not in coefficients:
                    coefficients[key] = np.full((5, 5), F(0), dtype=object)
                coefficients[key][i, j] += coefficient * multiplier
    return coefficients


def hessian_certificate(kappa):
    coefficients = polynomial_hessian(kappa)
    pairs = [(i, j) for i in range(5) for j in range(i, 5)]
    unknown_basis = []
    for i, j in pairs:
        q = np.full((5, 5), F(0), dtype=object)
        q[i, j] = q[j, i] = F(1)
        unknown_basis.append(q)
    rows = []
    for exponent in sorted(coefficients):
        h = coefficients[exponent]
        commutators = [q @ h - h @ q for q in unknown_basis]
        for i in range(5):
            for j in range(i + 1, 5):
                row = [c[i, j] for c in commutators]
                if any(row):
                    rows.append(row)
    pivots, basis = nullspace(rows, len(pairs))
    matrices = [sum((x * q for x, q in zip(vector, unknown_basis)),
                    np.full((5, 5), F(0), dtype=object)) for vector in basis]
    assert len(basis) == (1 if kappa else 2)
    assert all(np.all(q @ h == h @ q) for q in matrices for h in coefficients.values())
    eval_errors = []
    for point in ([0, 0, 0, 0, 0], [1, 2, -1, 0, 2], [F(1, 2), F(1, 3), 0, F(2, 3), -F(1, 2)]):
        point = list(map(F, point))
        hpoly = np.full((5, 5), F(0), dtype=object)
        for exponent, coeff in coefficients.items():
            monomial = math.prod(p**e for p, e in zip(point, exponent))
            hpoly += monomial * coeff
        phi = np.array(point[:4], dtype=object); sigma = point[4]
        radius = sum(p*p for p in phi)
        direct = np.full((5, 5), F(0), dtype=object)
        for i in range(4):
            for j in range(4):
                direct[i, j] = 2 * phi[i] * phi[j]
                if i == j:
                    direct[i, j] += radius - 1 + F(kappa, 2) * sigma**2
            direct[i, 4] = direct[4, i] = kappa * phi[i] * sigma
        direct[4, 4] = 1 + 3 * sigma**2 + F(kappa, 2) * radius
        assert np.array_equal(hpoly, direct)
        eval_errors.append(float(np.linalg.norm(np.array(hpoly - direct, dtype=float))))
    return dict(kappa_exact=str(kappa), coefficient_monomial_count=len(coefficients),
                symmetric_unknown_dimension=len(pairs), nonzero_equation_count=len(rows),
                exact_constraint_rank=len(pivots), exact_commutant_dimension=len(basis),
                monomial_exponents=[list(e) for e in sorted(coefficients)],
                coefficient_matrices=[[[str(v) for v in row] for row in coefficients[e]] for e in sorted(coefficients)],
                nullspace_matrices=[[[str(v) for v in row] for row in q] for q in matrices],
                exact_closed_formula_evaluation_errors=eval_errors,
                all_polynomial_coefficients_tested=True,
                sampled_Hessians_not_used_to_prove_commutant=True,
                potential_algebra_only_not_full_gauge_dynamics=True)


def gauge_truncation_certificate():
    rows = []
    for t in (0., .3, 1.1, 2.7):
        phi = np.array([0., np.exp(-1j*t)], dtype=complex)
        dot = -1j * phi
        ddot = -phi
        norm = float(np.vdot(phi, phi).real)
        dV = (2 * norm - 1) * phi
        d_up_0 = -dot
        # With D=partial-i*g*y*A and signature -+++, this is dL_m/dA_0 / g.
        coefficient = -1j * .5 * (np.vdot(phi, d_up_0) - np.vdot(d_up_0, phi))
        rows.append(dict(t=t, scalar_equation_residual=float(np.linalg.norm(ddot + dV)),
                         sigma_equation_residual=0., hypercharge_action_source_over_g=float(coefficient.real),
                         source_imaginary_error=float(abs(coefficient.imag)),
                         Maxwell_left_hand_divergence_at_A_zero=0.,
                         Maxwell_Euler_residual_over_g=float(coefficient.real)))
    assert max(row["scalar_equation_residual"] for row in rows) < 1e-12
    assert all(abs(row["hypercharge_action_source_over_g"] - 1) < 1e-12 for row in rows)
    return dict(signature="-+++", covariant_derivative="D_mu=partial_mu-i*g*y*A_mu",
                source_definition="delta L_m / delta A_0 / g; Maxwell Euler equation partial_mu F^{mu0} + delta L_m/delta A_0 = 0",
                physical_Maxwell_RHS_current_has_opposite_sign=True, hypercharge_y_exact="1/2",
                solution="Phi=(0,exp(-it)), sigma=0, A=0, lambda_H=v^2=omega^2=1",
                exact_source_coefficient_over_g="1", time_samples=rows,
                scalar_only_equations_satisfied=True, full_gauge_equations_satisfied=False,
                A_zero_not_a_consistent_all_scalar_data_truncation=True,
                no_failure_of_full_parent_dynamics_inferred=True,
                no_universal_gravity_coupling_conclusion_inferred=True)


def run():
    fields = field_certificate()
    hessians = [hessian_certificate(1), hessian_certificate(0)]
    truncation = gauge_truncation_certificate()
    return dict(round=1024, date="2026-10-08", all_audit_certificates_passed=True,
                new_cognitive_axioms=0, new_physical_calibration_groups=0,
                new_audit_certificate_groups=1, cumulative_research_groups=3801,
                scope="Typed common-parent dependency audit: explicit finite premise and interface certificates; not new physical predictions, instruments, a full shared process, or a cosmic countermodel.",
                parent_hypotheses_derived_from_cognition=False, full_parent_process_certified=False,
                machine_checks_prove_human_applicability_classifications=False,
                manual_bridge_ledger_checked_by_delivery_verifier=True,
                chiral_mass_interface=fields, Higgs_portal_commutant=hessians,
                gauge_truncation=truncation,
                historical_source_sha256={name: hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in HISTORICAL})


def compare(fresh, saved, path="root"):
    if isinstance(fresh, dict):
        assert fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key], saved[key], path+"."+key)
    elif isinstance(fresh, list):
        assert len(fresh) == len(saved), path
        for i, (a, b) in enumerate(zip(fresh, saved)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(fresh, float):
        assert math.isclose(fresh, saved, rel_tol=2e-10, abs_tol=2e-12), (path, fresh, saved)
    else:
        assert fresh == saved, (path, fresh, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    else:
        compare(result, json.loads(OUT.read_text("utf8")))
    print(json.dumps({"round": 1024, "passed": True, "dimension_four_rank": result["chiral_mass_interface"]["dimension_four_maximum_rank"],
                      "Weinberg_rank": result["chiral_mass_interface"]["Weinberg_extended_exact_rank"],
                      "commutant_dimensions": [r["exact_commutant_dimension"] for r in result["Higgs_portal_commutant"]],
                      "scalar_truncation_full_gauge_solution": False,
                      "new_physical_calibration_groups": 0, "new_audit_certificate_groups": 1,
                      "cumulative_research_groups": 3801}, ensure_ascii=False, indent=2))
