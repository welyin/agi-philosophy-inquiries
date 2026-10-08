"""1018: finite exact calibrations of a conditional gravity-ideal obstruction.

The all-local-correction / regular-open-algebra argument is analytic, not a
numerical extrapolation. This file checks explicit Killing jets, actual
sequential scalar variations, the full transported ideal product, and the
initial-data identities used in the first-order recoil lemma.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "gravity_ideal_selection_results.json"


def zeros(m, n=None):
    return [[F(0) for _ in range(m if n is None else n)] for _ in range(m)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum((x*y for x, y in zip(row, col)), F(0))
             for col in zip(*b)] for row in a]


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    return value


def defect(kappa, q):
    return [[q[a]*q[b] - (kappa[a]*q[a] if a == b else 0)
             for b in range(len(q))] for a in range(len(q))]


def is_zero(a):
    return all(v == 0 for row in a for v in row)


def support_condition(kappa, q):
    active = [i for i, v in enumerate(q) if v]
    return not active or (len(active) == 1 and
                          q[active[0]] == kappa[active[0]])


def series_product(a, b, n):
    return [sum((a[j]*b[k-j] for j in range(k+1)
                 if j < len(a) and k-j < len(b)), F(0))
            for k in range(n+1)]


def potential_gradient_series(potential, series, n):
    result = [[F(0)]*(n+1) for _ in series]
    for powers, coefficient in potential.items():
        for axis, power in enumerate(powers):
            if not power:
                continue
            exponent = list(powers)
            exponent[axis] -= 1
            product = [F(1)] + [F(0)]*n
            for coordinate, count in enumerate(exponent):
                for _ in range(count):
                    product = series_product(product, series[coordinate], n)
            for order in range(n+1):
                result[axis][order] += coefficient*power*product[order]
    return result


def scalar_solution_jets(potential):
    """Taylor jets of the analytic local ODE f'' = grad V(f), not a solution
    claim for the finite truncated polynomial or for the coupled theory."""
    degree = 10
    series = [[F(1, 5), F(2, 7)] + [F(0)]*(degree-1),
              [F(-1, 10), F(1, 6)] + [F(0)]*(degree-1)]
    for n in range(degree-1):
        rhs = potential_gradient_series(potential, series, n)
        for i in range(2):
            series[i][n+2] = rhs[i][n] / ((n+2)*(n+1))
    rhs = potential_gradient_series(potential, series, degree-2)
    residual = [[(n+2)*(n+1)*series[i][n+2]-rhs[i][n]
                 for n in range(degree-1)] for i in range(2)]
    assert is_zero(residual)
    gradient = zeros(4, 2)
    gradient[2] = [row[1] for row in series]
    hessian = [[[F(0), F(0)] for _ in range(4)] for _ in range(4)]
    hessian[2][2] = [2*row[2] for row in series]
    return series, residual, gradient, hessian


def killing_checks(gradient, hessian):
    xi = [F(0), F(1), F(0), F(0)]
    eta = [F(0)]*4
    dxi = zeros(4)
    deta = zeros(4)
    deta[1][2] = F(1)
    deta[2][1] = F(-1)
    signs = [-1, 1, 1, 1]
    d_xi = [[signs[nu]*dxi[mu][nu] + signs[mu]*dxi[nu][mu]
             for nu in range(4)] for mu in range(4)]
    d_eta = [[signs[nu]*deta[mu][nu] + signs[mu]*deta[nu][mu]
              for nu in range(4)] for mu in range(4)]
    assert is_zero(d_xi) and is_zero(d_eta)
    bracket = [sum((xi[mu]*deta[mu][nu]-eta[mu]*dxi[mu][nu]
                    for mu in range(4)), F(0)) for nu in range(4)]
    assert bracket == [0, 0, 1, 0]

    def sequential(left, right, dright):
        return [sum((left[mu]*(dright[mu][nu]*gradient[nu][i] +
                     right[nu]*hessian[mu][nu][i])
                     for mu in range(4) for nu in range(4)), F(0))
                for i in range(2)]

    # Variations act on fields; external parameters do not vary. Thus the
    # ordered gauge difference below is L_xi L_eta - L_eta L_xi.
    eta_after_xi = sequential(xi, eta, deta)
    xi_after_eta = sequential(eta, xi, dxi)
    commutator = [a-b for a, b in zip(eta_after_xi, xi_after_eta)]
    bracket_action = [sum((bracket[mu]*gradient[mu][i]
                          for mu in range(4)), F(0)) for i in range(2)]
    assert commutator == bracket_action == [F(2, 7), F(1, 6)]
    return dict(xi=xi, eta_at_origin=eta, derivative_eta=deta,
                D_xi=d_xi, D_eta=d_eta,
                all_higher_parameter_derivatives_zero_by_affinity=True,
                all_derivatives_of_D_xi_D_eta_zero=True,
                bracket=bracket,
                delta_eta_after_delta_xi=eta_after_xi,
                delta_xi_after_delta_eta=xi_after_eta,
                leading_commutator=commutator,
                convention="delta_eta(delta_xi Phi)-delta_xi(delta_eta Phi)",
                parameter_jets_retained=True)


def exact_rotation(kappa, q, o):
    n = len(q)
    qt = [sum((o[a][i]*q[i] for i in range(n)), F(0)) for a in range(n)]
    at = [[[sum((kappa[i]*o[c][i]*o[a][i]*o[b][i]
                 for i in range(n)), F(0)) for b in range(n)]
           for a in range(n)] for c in range(n)]
    dt = [[qt[a]*qt[b]-sum((qt[c]*at[c][a][b] for c in range(n)), F(0))
           for b in range(n)] for a in range(n)]
    transported = multiply(multiply(o, defect(kappa, q)), transpose(o))
    assert dt == transported
    return dict(O=o, q_transformed=qt, transported_product=at,
                defect_transformed=dt, covariance_exact=True)


def alignment_control():
    s = [[F(1), F(1)], [F(-1), F(1)]]
    kappa = [F(2), F(2)]
    q = [F(1), F(1)]
    qhat = [sum((s[a][i]*q[i] for i in range(2)), F(0))/2
            for a in range(2)]  # q' / sqrt(2)
    ahat = [[[sum((kappa[i]*s[c][i]*s[a][i]*s[b][i]
                   for i in range(2)), F(0))/4 for b in range(2)]
             for a in range(2)] for c in range(2)]  # A' / sqrt(2)
    dt = [[2*(qhat[a]*qhat[b]-sum((qhat[c]*ahat[c][a][b]
                                 for c in range(2)), F(0)))
           for b in range(2)] for a in range(2)]
    expected = [[F(0), F(0)], [F(0), F(-2)]]
    covariance = [[v/2 for v in row] for row in
                  multiply(multiply(s, defect(kappa, q)), transpose(s))]
    assert dt == expected == covariance
    wrong = [[[ahat[c][a][b] if a == b == c else F(0)
               for b in range(2)] for a in range(2)] for c in range(2)]
    fake_d = [[2*(qhat[a]*qhat[b]-sum((qhat[c]*wrong[c][a][b]
                                     for c in range(2)), F(0)))
               for b in range(2)] for a in range(2)]
    assert is_zero(fake_d) and ahat[0][1][1] == 1
    o = np.array(s, float)/math.sqrt(2)
    at_float = np.einsum("i,ci,ai,bi->cab", np.array(kappa, float), o, o, o)
    qt_float = o @ np.array(q, float)
    dt_float = np.outer(qt_float, qt_float)-np.einsum("c,cab->ab", qt_float, at_float)
    error = max(float(np.max(np.abs(at_float-np.array(ahat, float)*math.sqrt(2)))),
                float(np.max(np.abs(dt_float-np.array(expected, float)))))
    assert error < 3e-14
    return dict(kappa=kappa, q=q, O_times_sqrt2=s,
                q_prime_div_sqrt2=qhat, A_prime_div_sqrt2=ahat,
                correct_defect_prime=dt,
                all_mixed_product_entries_wrongly_deleted_defect=fake_d,
                deleted_A_1_22_div_sqrt2=ahat[0][1][1],
                deleting_mixed_product_changes_the_theory=True,
                source_alignment_alone_sufficient=False,
                maximum_numpy_cross_check_residual=error)


def pf_initial_constraints():
    """Finite spatial jets of the proposed de Donder initial data.

    No PDE solve and no claim that these arbitrary declared source components
    already give a complete coupled scalar-gravity solution.
    """
    rows = []
    for sample in range(4):
        uh = [[F((sample+1)*(i+j+1), 5+sample) for j in range(3)]
              for i in range(3)]
        dw = [[F((sample+1)*(2*i-j), 7) for j in range(3)] for i in range(3)]
        ddw = [[[F((sample+1)*(l+i+1)*(j+2)-(l*i), 11)
                 for j in range(3)] for i in range(3)] for l in range(3)]
        assert all(ddw[l][i] == ddw[i][l] for i in range(3) for l in range(3))
        trace_dw = sum((dw[k][k] for k in range(3)), F(0))
        p = [[dw[i][j]+dw[j][i]-(trace_dw if i == j else 0)
              for j in range(3)] for i in range(3)]
        dp = [[[ddw[l][i][j]+ddw[l][j][i] -
                 (sum((ddw[l][k][k] for k in range(3)), F(0)) if i == j else 0)
                 for j in range(3)] for i in range(3)] for l in range(3)]
        divp = [sum((dp[i][i][j] for i in range(3)), F(0)) for j in range(3)]
        lapw = [sum((ddw[i][i][j] for i in range(3)), F(0)) for j in range(3)]
        lapu = sum((uh[i][i] for i in range(3)), F(0))
        source_0nu = [lapu]+[-v for v in lapw]
        # H=-p_0nu + partial_i barh_i_nu=0 on the declared initial slice.
        h_constraint = [F(0)]*4
        dot_h = [source_0nu[0]-lapu]+[source_0nu[j+1]+divp[j] for j in range(3)]
        assert divp == lapw and all(v == 0 for v in dot_h)
        rows.append(dict(u_hessian=uh, W_first_jet=dw, W_second_jet=ddw,
                         p_ij=p, spatial_divergence_p=divp, laplacian_W=lapw,
                         J_0nu=source_0nu, H_nu=h_constraint, dt_H_nu=dot_h))
    return dict(samples=rows, exact_initial_constraints_pass=True,
                wave_operator="-partial_t^2 + spatial_Laplacian",
                source_relation="Box bar_h_mu_nu = J_mu_nu",
                source_is_not_claimed_to_be_a_complete_scalar_backreaction=True,
                conservation_and_local_PDE_existence_provided_by_analytic_lemma=True,
                numerical_PDE_solution_claimed=False)


def projector_sanity(previous):
    source = previous["mixed_vertex_counterexample"]
    projectors = [[[F(v) for v in row] for row in source[name]]
                  for name in ("Q_plus", "Q_minus")]
    kappa = [F(2), F(3)]
    qs = [[[kappa[a]*v for v in row] for row in projectors[a]] for a in range(2)]
    residuals = []
    for a in range(2):
        for b in range(2):
            product = multiply(qs[a], qs[b])
            residual = [[product[i][j] - (kappa[a]*qs[a][i][j] if a == b else 0)
                         for j in range(2)] for i in range(2)]
            assert is_zero(residual)
            residuals.append(residual)
    return dict(projectors=projectors, kappa=kappa, matter_representation_matrices=qs,
                product_residuals=residuals, both_sources_nonzero=True,
                distinct_decoupled_matter_sectors=True,
                contradicted_by_single_irreducible_sector_claim=False)


def compare(fresh, saved, path="root"):
    if isinstance(fresh, float):
        assert isinstance(saved, (int, float)) and math.isclose(fresh, saved, rel_tol=2e-11, abs_tol=2e-13), path
    elif isinstance(fresh, dict):
        assert isinstance(saved, dict) and fresh.keys() == saved.keys(), path
        for key in fresh:
            compare(fresh[key], saved[key], path+"."+key)
    elif isinstance(fresh, list):
        assert isinstance(saved, list) and len(fresh) == len(saved), path
        for i, (a, b) in enumerate(zip(fresh, saved)):
            compare(a, b, path+f"[{i}]")
    else:
        assert fresh == saved, (path, fresh, saved)


def run():
    previous_path = BASE/"archive_1009_/1017/stress_source_commutant_results.json"
    previous = json.loads(previous_path.read_text(encoding="utf8"))
    family = next(row for row in previous["exact_polynomial_families"]
                  if row["name"] == "genuine_mixed_quartic")
    assert family["symmetric_commutant_dimension"] == 1
    potential = {tuple(map(int, powers.split(","))): F(value)
                 for powers, value in family["potential"].items()}
    series, ode_residual, grad, hess = scalar_solution_jets(potential)
    killing = killing_checks(grad, hess)
    cases = [("one_active", [2, 3], [2, 0]),
             ("two_active", [2, 3], [2, 3]),
             ("wrong_self_coefficient", [2, 3], [1, 0]),
             ("free_ideal_nonzero_coupling", [0, 3], [1, 0]),
             ("zero_source", [2, 3], [0, 0]),
             ("three_ideals_one_active", [F(3, 2), -2, 0], [0, -2, 0]),
             ("four_ideals_wrong_support", [1, 2, 3, 4], [0, 2, 0, 4])]
    rows = []
    for name, kk, qq in cases:
        kk, qq = list(map(F, kk)), list(map(F, qq))
        d = defect(kk, qq)
        residuals = []
        for a in range(len(qq)):
            for b in range(len(qq)):
                # Independently calculated sequential derivative, minus the
                # COMPLETE adopted first-order bracket action.
                sequential = [qq[a]*qq[b]*v for v in killing["leading_commutator"]]
                bracket_term = [(kk[a]*qq[a]*v if a == b else F(0))
                                for v in killing["leading_commutator"]]
                remaining = [v-w for v, w in zip(sequential, bracket_term)]
                assert remaining == [d[a][b]*v for v in grad[2]]
                residuals.append(dict(a=a, b=b, commutator=sequential,
                                      complete_B1_action=bracket_term,
                                      residual=remaining))
        assert is_zero(d) == support_condition(kk, qq)
        rows.append(dict(name=name, kappa=kk, q=qq, D=d,
                         obstruction_zero=is_zero(d), sequential_checks=residuals))
    enumerated = 0
    for kk in ([F(2), F(3)], [F(0), F(3, 2), F(-2)],
               [F(1), F(2), F(3), F(4)]):
        choices = [sorted({F(0), v, v/2, F(1)}) for v in kk]
        for qq in itertools.product(*choices):
            assert is_zero(defect(kk, qq)) == support_condition(kk, qq)
            enumerated += 1
    rotated = []
    for n in (3, 4):
        o = [[F(int(a == b)) for b in range(n)] for a in range(n)]
        o[0][0], o[0][1], o[1][0], o[1][1] = F(3, 5), F(4, 5), F(-4, 5), F(3, 5)
        assert multiply(o, transpose(o)) == [[int(a == b) for b in range(n)] for a in range(n)]
        rotated.append(exact_rotation([F(i+1) for i in range(n)],
                                      [F(2, i+1) for i in range(n)], o))
    sources = [BASE/"archive_301_341/research_note_326.md",
               BASE/"archive_342_369/research_note_358.md",
               BASE/"archive_1009_/research_note_1017.md", previous_path,
               BASE/"archive_1009_/1009/input_dependency_ledger_v0_1.md"]
    return encode(dict(round=1018, new_calibration_groups=1,
        cumulative_test_groups=3796, new_cognitive_axioms=0,
        all_scientific_calibrations_passed=True,
        scope="fixed complete minimal first-order Noether data, positive-kinetic spin2 ideals, canonical irreducible scalar sector; conditional necessary coupling obstruction",
        first_order_minimal_data_adopted=True,
        higher_local_corrections_arbitrarily_assumed_zero=False,
        regular_open_algebra_allowed=True, M1_assumed_zero=False,
        first_order_recoil_exists_by_analytic_lemma=True,
        nonlinear_completion_unique=False,
        arbitrary_R2_or_open_terms_numerically_verified=False,
        full_backreaction_solution_numerically_produced=False,
        Einstein_action_derived_by_this_code=False,
        gravity_generated=False, number_of_spin2_fields_selected=False,
        killing_jet_calibration=killing,
        scalar_on_shell_local_jet_calibration=dict(
            potential=family["potential"], exact_commutant_dimension=1,
            series_coefficients=series, equation_residual_coefficients=ode_residual,
            equation="f_second_derivative = gradient_V(f)",
            degree=10, verified_equation_orders=list(range(9)),
            finite_Taylor_polynomial_claimed_exact_solution=False,
            actual_local_analytic_solution_uses_ODE_existence=True),
        exact_coupling_cases=rows,
        finite_support_classification_samples=enumerated,
        all_parameter_classification_is_analytic_not_sampling=True,
        source_alignment_control=alignment_control(),
        rational_orthogonal_rotation_controls=rotated,
        first_order_recoil_initial_constraints=pf_initial_constraints(),
        separated_matter_boundary=projector_sanity(previous),
        analytic_obligations=[
            "Killing jets make delta0 R2 vanish for every regular finite-jet local R2",
            "actual free matter solution plus conserved-source first-order gravitational recoil gives E=O(t^2)",
            "regular open term t M1 E is then O(t^3), without assuming M1=0",
            "D_ab=0 implies q=0 or q_r=kappa_r for one nonfree ideal",
            "one-active branch has an adopted Einstein-scalar completion; its uniqueness is not claimed"],
        retained_freedom=["zero-source branch", "spectator free or self-interacting gravity ideals",
            "nonzero kappa magnitude and normalization", "vacuum constants and independent invariant improvements",
            "nonminimal first-order data outside the contract", "other matter and higher derivative completions"],
        historical_source_sha256={str(p.relative_to(BASE)).replace("\\", "/"):
                                  hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1018, passed=result["all_scientific_calibrations_passed"],
        coupling_cases=len(result["exact_coupling_cases"]),
        support_samples=result["finite_support_classification_samples"],
        rotated_defect=result["source_alignment_control"]["correct_defect_prime"],
        recoil_initial_data_samples=len(result["first_order_recoil_initial_constraints"]["samples"])), indent=2))
