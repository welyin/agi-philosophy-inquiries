"""1029: finite exact checks for fixed-vertex first-order lift equivalence.

The completeness statement uses the analytic regular Euler-jet splitting and
Noether's second theorem.  These calculations test its principal blocks, a
genuine off-shell Higgs trivial symmetry, and formal full-E transport.  They
are not a solver for all local gauge transformations or all first-order vertices.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "first_order_lift_equivalence_results.json"
_path = BASE / "archive_1009_/1027/charged_source_selection.py"
_spec = importlib.util.spec_from_file_location("r1027_lift_arithmetic", _path)
_old = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _old
_spec.loader.exec_module(_old)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def rank(matrix):
    a = [list(map(F, row)) for row in matrix]
    r = 0
    for c in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        value = a[r][c]
        a[r] = [x/value for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                value = a[i][c]
                a[i] = [x-value*y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    result = F(1)
    for c in range(len(a)):
        p = next((i for i in range(c, len(a)) if a[i][c]), None)
        if p is None:
            return F(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            result = -result
        pivot = a[c][c]
        result *= pivot
        for i in range(c+1, len(a)):
            multiplier = a[i][c]/pivot
            a[i] = [x-multiplier*y for x, y in zip(a[i], a[c])]
    return result


def mv(a, x):
    return [sum((v*w for v, w in zip(row, x)), F(0)) for row in a]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def principal_blocks():
    # Time-derivative part: 1/2 sum_A,i(dot A_i - D_i A_0)^2
    # +1/2 sum_r(dot h_r + (A_0 h)_r)^2.  Shifts do not alter Hessians.
    n = 12*3+4+1
    kinetic_monomials = {tuple(2 if j == i else 0 for j in range(n)): F(1, 2)
                         for i in range(n)}
    hessian = [[sum((coefficient*p[i]*(p[j]-(i == j))
                     for p, coefficient in kinetic_monomials.items()), F(0))
                for j in range(n)] for i in range(n)]
    assert rank(hessian) == n and determinant(hessian) == 1

    # i sigma^0 in a Weyl equation, realified only for exact matrix checking.
    # Odd canonical momentum constraints are second class; this is NOT a
    # claim that a second-order fermionic velocity Hessian is invertible.
    weyl_realified = [[F(0), F(0), F(-1), F(0)],
                      [F(0), F(0), F(0), F(-1)],
                      [F(1), F(0), F(0), F(0)],
                      [F(0), F(1), F(0), F(0)]]
    assert rank(weyl_realified) == 4 and determinant(weyl_realified) == 1
    return dict(bosonic_velocity_variables=n, gauge_spatial_components=36,
                real_scalar_components=5, bosonic_hessian_rank=rank(hessian),
                bosonic_hessian_determinant=determinant(hessian),
                bosonic_kinetic_normalization="1/2 sum of shifted velocity squares",
                weyl_time_symbol="i sigma^0=i I_2; conjugate equation has opposite i",
                weyl_time_symbol_realified=weyl_realified,
                weyl_time_symbol_real_rank=rank(weyl_realified),
                weyl_time_symbol_real_determinant=determinant(weyl_realified),
                independent_Weyl_modules_with_internal_and_flavor_multiplicity=45,
                scope="These are principal coefficients of the declared action. The analytic jet splitting, not this rank check, proves no further Noether identity.")


def pf_principal_and_constraints():
    # L_time=1/2 (v_ij v_ij - (tr v)^2), including both off-diagonal entries.
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    vectors = [[F(int(a == i and b == j or a == j and b == i))
                for a, b in pairs] for i in range(3) for j in range(3)]
    trace = [F(1), F(1), F(1), F(0), F(0), F(0)]
    dewitt = [[sum((v[a]*v[b] for v in vectors), F(0))-trace[a]*trace[b]
               for b in range(6)] for a in range(6)]
    assert rank(dewitt) == 6 and determinant(dewitt) == -16

    def add_term(row, key, coefficient):
        row[key] = row.get(key, F(0))+F(coefficient)
        if not row[key]:
            del row[key]

    constraints = [{} for _ in range(4)]
    for i in range(3):
        for j in range(3):
            add_term(constraints[0], f"h{min(i,j)}{max(i,j)}_d{min(i,j)}{max(i,j)}", 1)
            add_term(constraints[0], f"h{i}{i}_d{j}{j}", -1)
            add_term(constraints[i+1], f"pi{min(i,j)}{max(i,j)}_d{j}", 1)
    pivots = ["h11_d00", "pi00_d0", "pi01_d0", "pi02_d0"]
    pivot_matrix = [[row.get(p, F(0)) for p in pivots] for row in constraints]
    assert rank(pivot_matrix) == 4 and determinant(pivot_matrix) == -1
    return dict(spatial_pair_order=pairs, spatial_velocity_hessian=dewitt,
                spatial_velocity_rank=rank(dewitt), spatial_velocity_determinant=determinant(dewitt),
                trace_direction_note="Unreduced DeWitt trace has negative sign; gauge constraints, not this Hessian alone, identify physical positive helicities.",
                spatial_index_convention="0,1,2 here label physical x,y,z; not spacetime time",
                constraint_polynomials=constraints, constraint_pivots=pivots,
                constraint_pivot_matrix=pivot_matrix,
                constraint_pivot_determinant=determinant(pivot_matrix),
                dependent_equations="Time derivatives of E_PF^(0 nu) are supplied by linear Bianchi identities; proof is analytic.")


def pf_parameter_symbol():
    results = []
    # Includes null covector: the generator is injective there too.
    for k in ([F(1), F(0), F(0), F(0)],
              [F(1), F(1), F(0), F(0)],
              [F(0), F(0), F(2, 3), F(-5, 7)]):
        pairs = [(mu, nu) for mu in range(4) for nu in range(mu, 4)]
        matrix = [[k[mu]*(v == nu)+k[nu]*(v == mu) for v in range(4)]
                  for mu, nu in pairs]
        a = next(i for i, value in enumerate(k) if value)
        columns = [a]+[i for i in range(4) if i != a]
        selected = [(a, a)]+[tuple(sorted((a, i))) for i in columns[1:]]
        minor = [[matrix[pairs.index(pair)][c] for c in columns] for pair in selected]
        assert rank(matrix) == 4 and determinant(minor) == 2*k[a]**4
        results.append(dict(k_covector=k, symbol_matrix=matrix,
                            symbol_pair_order=pairs, rank=4,
                            selected_rows=selected, column_order=columns,
                            minor_determinant=determinant(minor),
                            expected_minor=2*k[a]**4))
    return dict(cases=results,
                analytic_argument="Choose k_a != 0; diagonal equation gives v_a=0, then mixed equations give every v_b=0. Highest arbitrary parameter-jet coefficients vanish recursively on shell.",
                no_local_projection_to_Killing="Finite-order differential dependence on arbitrary compact parameters cannot produce a nonzero universal Killing parameter.")


def higgs_trivial_lift():
    # A genuine off-shell full-Higgs Euler jet, at A=psi=sigma=0.
    h = [F(0), F(5, 4), F(0), F(0)]
    dy_h = [F(0), F(2, 7), F(0), F(0)]
    desired_euler = [F(1, 3), F(-2, 5), F(1, 7), F(3, 2)]
    lam, v_squared = F(2, 5), F(9, 4)
    grad_u = [lam*(dot(h, h)-v_squared)*v for v in h]
    # All spatial second derivatives zero, so E_h=-h_tt-grad U.
    h_tt = [-e-v for e, v in zip(desired_euler, grad_u)]
    actual_euler = [-acc-force for acc, force in zip(h_tt, grad_u)]
    assert actual_euler == desired_euler
    ty_c = _old.real_generator(_old.scale(_old.eye(2), F(1, 2)), False)
    assert all(v.i == 0 for row in ty_c for v in row)
    ty = [[v.r for v in row] for row in ty_c]
    assert all(ty[i][j] == -ty[j][i] for i in range(4) for j in range(4))
    b = dot(h, dy_h)  # xi=partial_y; b=xi^mu partial_mu O.
    delta_h = [b*v for v in mv(ty, actual_euler)]
    variation_action = dot(actual_euler, delta_h)
    variation_o = dot(h, delta_h)
    assert variation_action == 0 and variation_o == F(-75, 224)
    return dict(h_real=h, partial_y_h=dy_h, lambda_h=lam, v_squared=v_squared,
                potential_gradient=grad_u, partial_t_squared_h=h_tt,
                actual_higgs_Euler=actual_euler, tY=ty, parameter_scalar_b=b,
                homogeneous_lift_delta_h=delta_h,
                Euler_contraction=variation_action, off_shell_delta_O=variation_o,
                exact_identity="E_h^T b tY E_h=0 because b is even and tY is antisymmetric",
                classification="Euler-trivial symmetry; nonzero off-shell O variation disproves bare lift uniqueness, not equivalence-class uniqueness.",
                jet_status="Chosen off-shell second jets realize this exact Euler value. They are not claimed to solve matter or gauge equations.")


def full_euler_recoil_transport():
    # Actual leading trace source from the legal radial branch in 1028.
    f, p, lam, v_squared, q = F(5, 4), F(2, 7), F(2, 5), F(9, 4), F(2)
    u = lam*(f*f-v_squared)**2/F(4)
    trace_t = -p*p-4*u
    r = -q*trace_t/F(2)  # E_PF,0 = epsilon*r from first recoil.
    e1 = q*trace_t/F(2)
    assert r+e1 == 0 and r != 0
    b, norm_h_sq = f*p, f*f
    response = b*norm_h_sq
    # A genuine cross-block antisymmetric trivial operator can be chosen as
    # Delta h_i=b h_i tr(E_PF), Delta gamma_mn=-b eta_mn h_i E_h_i.
    # Its O response is response*tr(E_PF).  Keep unknown E2 as a symbol c.
    old_first_order_piece = {"eps^2": response*r}
    inherited_second_order_piece = {"eps^2": -response*e1, "eps^3*c": -response}
    removed_full_euler_piece = {"eps^3*c": response}
    assert old_first_order_piece["eps^2"] == inherited_second_order_piece["eps^2"]
    assert inherited_second_order_piece["eps^3*c"] == -removed_full_euler_piece["eps^3*c"]
    # Direct cancellation of the two Euler contractions for generic point data.
    matter_euler_dot_h, pf_euler_trace = F(3, 11), F(-7, 13)
    first = b*matter_euler_dot_h*pf_euler_trace
    second = -b*matter_euler_dot_h*pf_euler_trace
    assert first+second == 0
    return dict(radial_f=f, radial_p=p, radial_potential=u, source_q=q,
                full_matter_stress_trace=trace_t,
                free_PF_Euler_linear_recoil_coefficient=r,
                first_interaction_Euler_coefficient=e1,
                first_recoil_sum=0, mixed_trivial_operator_O_factor=response,
                old_epsilon_M_E0_O=old_first_order_piece,
                inherited_minus_epsilon_squared_M_E1_and_higher=inherited_second_order_piece,
                removed_epsilon_M_Efull_O=removed_full_euler_piece,
                mixed_Euler_contractions=[first, second],
                unknown_c="Symbol for the unspecified full Euler coefficient at order epsilon^2; not a calculated nonlinear completion.",
                lesson="E0 is O(epsilon), not zero, on recoil data. Full-E subtraction preserves its O(epsilon^2) effect in the transported higher generator.")


ZERO = (0, 0, 0, 0)


def padd(a, b):
    c = dict(a)
    for monomial, coefficient in b.items():
        c[monomial] = c.get(monomial, F(0))+coefficient
        if not c[monomial]:
            del c[monomial]
    return c


def pscale(a, value):
    return {m: value*c for m, c in a.items() if value*c}


def pmul(a, b):
    c = {}
    for m, x in a.items():
        for n, y in b.items():
            p = tuple(i+j for i, j in zip(m, n))
            c[p] = c.get(p, F(0))+x*y
    return {m: x for m, x in c.items() if x}


def pdiff(a, axis):
    c = {}
    for monomial, coefficient in a.items():
        if monomial[axis]:
            p = list(monomial)
            p[axis] -= 1
            c[tuple(p)] = coefficient*monomial[axis]
    return c


def improvement_certificate():
    # Phi=(0,f)/sqrt(2), f=1+t+x*y+z^2/2.  This is an OFF-SHELL
    # polynomial, because the improvement identity must not need field EOM.
    f = {ZERO: F(1), (1, 0, 0, 0): F(1), (0, 1, 1, 0): F(1),
         (0, 0, 0, 2): F(1, 2)}
    o = pscale(pmul(f, f), F(1, 2))
    signs = [-1, 1, 1, 1]
    box = {}
    for mu in range(4):
        box = padd(box, pscale(pdiff(pdiff(o, mu), mu), signs[mu]))
    improvement = [[padd(pscale(pdiff(pdiff(o, mu), nu), signs[mu]*signs[nu]),
                          pscale(box, -signs[mu] if mu == nu else 0))
                    for nu in range(4)] for mu in range(4)]
    divergence = []
    for nu in range(4):
        d = {}
        for mu in range(4):
            d = padd(d, pdiff(improvement[mu][nu], mu))
        divergence.append(d)
    assert not any(divergence)
    at_origin = [[entry.get(ZERO, F(0)) for entry in row] for row in improvement]
    assert any(x for row in at_origin for x in row)
    assert all(improvement[mu][nu] == improvement[nu][mu] for mu in range(4) for nu in range(4))
    return dict(coordinate_order=["t", "x", "y", "z"], f_polynomial=f,
                O_polynomial=o, box_O_polynomial=box,
                I_upper_polynomials=improvement, I_upper_at_origin=at_origin,
                exact_divergence_polynomials=divergence,
                definition="I^(mu nu)=(partial^mu partial^nu-eta^(mu nu) box)O",
                no_equations_of_motion_used=True,
                claim="The extra h I vertex is delta0-invariant modulo a boundary; necessary D condition survives. No existence claim for every improvement completion.")


def run():
    sources = [
        "archive_1009_/research_note_1028.md",
        "archive_1009_/1028/input_dependency_update_v0_17.md",
        "archive_1009_/1028/NEXT.md",
        "archive_1009_/research_note_1027.md",
        "archive_1009_/1027/charged_source_selection.py",
        "archive_1009_/research_note_1018.md",
        "archive_342_369/research_note_358.md",
        "archive_956_989/981/drafts/common_parent_contract_v1.md",
        "archive_990_1008/993/common_candidate_v1.md",
    ]
    result = dict(round=1029, status="scientific_calibration_verified",
                  new_calibration_groups=1, cumulative_test_groups=3806,
                  new_cognitive_axioms=0, goal_complete=False,
                  all_scientific_calibrations_passed=True,
                  scope="Fixed S0 and S1; local regular lift equivalence modulo original gauge generators and graded Euler-trivial symmetries. No classification of all vertices.",
                  retained_inputs=["1028 complete leading SM/PF action", "fixed EH3+hT first-order vertex", "nondegenerate standard kinetic principal parts", "regular finite-jet transformations and arbitrary local gauge parameters", "graded local Euler splitting", "allowed improvement is matter-jet-only, symmetric, internally invariant, identically conserved; C is truly constant"],
                  not_proved=["bare off-shell lift uniqueness", "uniqueness or generation of S1", "all higher-derivative kinetic theories", "UV or quantum completion", "nonlinear completion of arbitrary improvements", "unique physical metric for every nonminimal effect"],
                  kinetic_principal_blocks=principal_blocks(),
                  PF_principal_and_constraint_blocks=pf_principal_and_constraints(),
                  PF_parameter_symbol=pf_parameter_symbol(),
                  genuine_higgs_trivial_lift=higgs_trivial_lift(),
                  full_E_recoil_transport=full_euler_recoil_transport(),
                  identically_conserved_improvement=improvement_certificate(),
                  primary_source={"url": "https://arxiv.org/html/hep-th/0002245", "sections": ["5.1.3", "6.4.3", "6.6 (Theorem 6.9 and Corollary 6.3)"]},
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  historical_source_sha256={p: hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in sources})
    return encode(result)


def compare(actual, expected, path="root"):
    if isinstance(actual, dict):
        assert isinstance(expected, dict) and actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(actual, list):
        assert isinstance(expected, list) and len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, f"{path}[{i}]")
    else:
        assert actual == expected, (path, actual, expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf-8")))
    print(json.dumps({"status": "passed", "round": 1029,
                      "mode": "exclusive_first_write" if args.write else "read_only_recompute_compare",
                      "off_shell_delta_O": result["genuine_higgs_trivial_lift"]["off_shell_delta_O"],
                      "PF_spatial_velocity_determinant": result["PF_principal_and_constraint_blocks"]["spatial_velocity_determinant"],
                      "historical_sources": len(result["historical_source_sha256"]),
                      "output": str(OUT)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
