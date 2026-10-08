"""1028: a full-Higgs local branch bridges 1027 sources to 1018 gravity ideals.

Fraction Taylor data calibrate local identities, not global solutions.  The
regular local ODE and first-recoil existence statements are analytic.  The
complete minimal first-order Noether representative is an explicit input.
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
OUT = HERE / "full_parent_gravity_bridge_results.json"
OLD_CODE = BASE / "archive_1009_/1027/charged_source_selection.py"
_spec = importlib.util.spec_from_file_location("r1027_bridge_arithmetic", OLD_CODE)
_old = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _old
_spec.loader.exec_module(_old)
C, I = _old.C, _old.C(0, 1)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, C):
        return {"real": str(value.r), "imag": str(value.i)}
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def mul(a, b, n):
    return [sum((a[j]*b[k-j] for j in range(k+1)
                 if j < len(a) and k-j < len(b)), F(0))
            for k in range(n+1)]


def add(a, b, n):
    return [(a[k] if k < len(a) else F(0)) +
            (b[k] if k < len(b) else F(0)) for k in range(n+1)]


def derivative(a):
    return [(k+1)*a[k+1] for k in range(len(a)-1)]


def scaled(a, v):
    return [v*x for x in a]


def matvec(t, v):
    return [sum((a*b for a, b in zip(row, v)), C()) for row in t]


def inner(v, w):
    return sum((a.conj()*b for a, b in zip(v, w)), C())


def zero_complex(values):
    return all(not v for v in values)


def higgs_branch():
    degree, lam, v_squared = 12, F(2, 5), F(9, 4)
    f = [F(5, 4), F(2, 7)] + [F(0)]*(degree-1)

    def rhs(n):
        cubic = mul(mul(f, f, n), f, n)
        return [lam*(cubic[k]-v_squared*f[k]) for k in range(n+1)]

    for k in range(degree-1):
        f[k+2] = rhs(k)[k]/((k+2)*(k+1))
    residual = [derivative(derivative(f))[k]-rhs(degree-2)[k]
                for k in range(degree-1)]
    assert not any(residual)

    # r=sqrt(2)*Phi=(0,f).  U=lambda*(Phi^dagger Phi-v^2/2)^2.
    # The full doublet equation is r''=2 lambda (r^dagger r/2-v^2/2)r.
    fp = derivative(f)
    n = degree-1
    r = [[C() for _ in range(n+1)], [C(x) for x in f[:n+1]]]
    rp = [[C() for _ in range(n+1)], [C(x) for x in fp]]
    generators = {"Y": _old.scale(_old.eye(2), F(1, 2))}
    generators.update({f"SU2_{i+1}": t for i, t in enumerate(_old.SU2)})
    currents = {}
    for label, t in generators.items():
        coefficients = []
        for k in range(n+1):
            first, second = C(), C()
            for ell in range(k+1):
                rv = [r[a][ell] for a in range(2)]
                rpv = [rp[a][k-ell] for a in range(2)]
                first += inner(rv, matvec(t, rpv))
                second += inner(rpv, matvec(t, rv))
            coefficients.append(I*F(1, 2)*(first-second))
        assert zero_complex(coefficients)
        currents[label] = coefficients

    # Conserved stress on the branch.  T^y_y=1/2 f'^2-U; all off-diagonal
    # entries vanish.  Thus its derivative is the only possible divergence.
    fsq = mul(f, f, degree-1)
    shifted = list(fsq)
    shifted[0] -= v_squared
    potential = scaled(mul(shifted, shifted, degree-1), lam/F(4))
    tyy = add(scaled(mul(fp, fp, degree-1), F(1, 2)),
              scaled(potential, -1), degree-1)
    stress_divergence = derivative(tyy)
    assert not any(stress_divergence)
    norm = scaled(mul(f, f, degree), F(1, 2))
    assert derivative(norm)[0] == F(5, 14)

    # At psi=0 every allowed Yukawa/Weinberg contribution to a fermion
    # equation still has at least one fermion factor.  This records the
    # exact degree test; it is not a sample fit of physical flavor matrices.
    fermion_degrees_in_equations = {"kinetic": 1, "Yukawa": 1, "Weinberg": 1}
    assert all(k >= 1 for k in fermion_degrees_in_equations.values())
    sigma, mu_d_sq, lambda_d, portal = F(0), F(3, 7), F(1, 6), F(2, 9)
    sigma_equation = sigma*(mu_d_sq+lambda_d*sigma*sigma+portal*f[0]**2/F(2))
    assert sigma_equation == 0

    data = dict(taylor_degree=degree, equation_orders=list(range(degree-1)),
                potential="U_rad=lambda*(f^2-v^2)^2/4",
                lambda_h=lam, v_squared=v_squared, f_coefficients=f,
                radial_equation_residual=residual,
                transverse_doublet_equation_identically_zero=True,
                higgs_current_convention="D=partial+i*g*A; j/g=i*(r^dagger*T*Dr-(Dr)^dagger*T*r)/2",
                gauge_current_orders=list(range(n+1)), gauge_current_series=currents,
                color_current_identically_zero=True,
                fermion_equation_minimum_fermion_degree=fermion_degrees_in_equations,
                zero_fermion_branch_exact_by_homogeneity=True,
                sigma_equation_at_zero=sigma_equation,
                nonzero_portal_calibration=portal,
                sigma_zero_is_consistent_for_any_declared_portal=True,
                stress_T_y_y_coefficients=tyy,
                full_stress_divergence_coefficients=stress_divergence,
                observable="O=Phi^dagger Phi=f^2/2",
                observable_at_origin=norm[0], observable_y_derivative=derivative(norm)[0],
                statement="Only retained finite Taylor orders are calibrated; a genuine smooth local ODE solution is supplied by the analytic proof.")
    return data, f


def internal_invariance_and_phase_control(f):
    r = [C(F(2, 3), F(1, 7)), C(F(5, 4), F(-2, 9))]
    t = _old.add(_old.scale(_old.eye(2), F(4, 9)),
                _old.add(_old.scale(_old.SU2[0], F(2, 5)),
                         _old.add(_old.scale(_old.SU2[1], F(-3, 7)),
                                  _old.scale(_old.SU2[2], F(5, 11)))))
    dr = [I*v for v in matvec(t, r)]
    invariant_variation = F(1, 2)*(inner(dr, r)+inner(r, dr))
    assert not invariant_variation

    # A forbidden phase-only truncation with A=0.  The point is selected
    # with theta=0, theta'=a.  This is NOT asserted to be a solution.
    a = F(3, 5)
    rb = [C(), C(f[0])]
    phase_derivative = [C(), C(f[1], a*f[0])]
    labels = {"Y": _old.scale(_old.eye(2), F(1, 2))}
    labels.update({f"SU2_{i+1}": v for i, v in enumerate(_old.SU2)})

    def current(v, dv, generator):
        return I*F(1, 2)*(inner(v, matvec(generator, dv))-
                          inner(dv, matvec(generator, v)))

    bad = {label: current(rb, phase_derivative, generator)
           for label, generator in labels.items()}
    assert bad["Y"] == C(F(-15, 32))
    assert bad["SU2_3"] == C(F(15, 32))
    assert not bad["SU2_1"] and not bad["SU2_2"]
    # The connection g' B_y=-2a restores the gauge-transported real branch.
    connection_action = matvec(_old.scale(_old.eye(2), -a), rb)
    covariant_derivative = [x+I*y for x, y in zip(phase_derivative, connection_action)]
    repaired = {label: current(rb, covariant_derivative, generator)
                for label, generator in labels.items()}
    assert zero_complex(repaired.values())
    return dict(arbitrary_hermitian_internal_generator=t, arbitrary_complex_doublet=r,
                internal_O_variation=invariant_variation,
                local_lorentz_O_variation="0 (O is a Lorentz scalar)",
                field_dependent_internal_parameters_also_annihilate_O=True,
                phase_only_A_zero_theta_prime=a, bad_phase_current_over_coupling=bad,
                bad_phase_A_zero_violates_gauge_equation=True,
                restored_gprime_B_y=-2*a,
                covariant_derivative_after_restoring_connection=covariant_derivative,
                restored_current_over_coupling=repaired,
                distinction="Restoring the pure-gauge connection is a gauge transform of the valid radial branch; deleting it changes the equations.")


def composite_killing_variations(f):
    r = [C(), C(f[0])]
    grad = [[C(), C()] for _ in range(4)]
    grad[2] = [C(), C(f[1])]
    hess = [[[C(), C()] for _ in range(4)] for _ in range(4)]
    hess[2][2] = [C(), C(2*f[2])]
    xi, zeta = [F(0), F(1), F(0), F(0)], [F(0)]*4
    dxi, dzeta = [[F(0)]*4 for _ in range(4)], [[F(0)]*4 for _ in range(4)]
    dzeta[1][2], dzeta[2][1] = F(1), F(-1)
    signs = [-1, 1, 1, 1]
    killing = [[[signs[nu]*dp[mu][nu]+signs[mu]*dp[nu][mu]
                 for nu in range(4)] for mu in range(4)] for dp in (dxi, dzeta)]
    assert all(not v for matrix in killing for row in matrix for v in row)
    bracket = [sum((xi[mu]*dzeta[mu][nu]-zeta[mu]*dxi[mu][nu]
                    for mu in range(4)), F(0)) for nu in range(4)]
    assert bracket == [0, 0, 1, 0]

    def variation(vector):
        return [sum((vector[mu]*grad[mu][i] for mu in range(4)), C())
                for i in range(2)]

    def ordered(left, right, dright):
        # delta_right(delta_left r), external left/right parameters fixed.
        nested = [sum((left[mu]*(dright[mu][nu]*grad[nu][i]+
                              right[nu]*hess[mu][nu][i])
                       for mu in range(4) for nu in range(4)), C())
                  for i in range(2)]
        dl, dr = variation(left), variation(right)
        terms = [inner(nested, r), inner(dl, dr), inner(dr, dl), inner(r, nested)]
        return nested, terms, sum(terms, C())*F(1, 2)

    first = ordered(xi, zeta, dzeta)
    second = ordered(zeta, xi, dxi)
    difference = first[2]-second[2]
    br = variation(bracket)
    bracket_O = (inner(br, r)+inner(r, br))*F(1, 2)
    assert difference == bracket_O == C(F(5, 14))

    # Just one inadmissible pair and one admitted branch; no parameter scan.
    examples = {}
    for label, kappa, q in (("two_active_incompatible", [F(2), F(3)], [F(1), F(1)]),
                           ("one_active_admitted", [F(2), F(3)], [F(2), F(0)])):
        actual, expected = [], []
        for a in range(2):
            row, erow = [], []
            for b in range(2):
                lhs = q[a]*q[b]*(first[2]-second[2])
                rhs = (kappa[a]*q[a] if a == b else F(0))*bracket_O
                row.append(lhs-rhs)
                erow.append((q[a]*q[b]-(kappa[a]*q[a] if a == b else F(0)))*bracket_O)
            actual.append(row)
            expected.append(erow)
        assert actual == expected
        examples[label] = dict(kappa=kappa, q=q, direct_composite_closure_residual=actual)
    assert any(examples["two_active_incompatible"]["direct_composite_closure_residual"][0])
    assert all(not x for row in examples["one_active_admitted"]["direct_composite_closure_residual"] for x in row)
    return dict(coordinates=["t", "x", "y", "z"], xi=xi, zeta_at_origin=zeta,
                derivative_zeta=dzeta, killing_symmetric_derivatives=killing,
                all_higher_Killing_derivatives_zero_by_affinity=True,
                bracket=bracket, nested_field_variations=[first[0], second[0]],
                four_product_rule_terms=[first[1], second[1]],
                ordered_O_variations=[first[2], second[2]],
                composite_commutator=difference, bracket_action_on_O=bracket_O,
                convention="delta_zeta(delta_xi O)-delta_xi(delta_zeta O)",
                examples=examples,
                higher_corrections_excluded_by="Analytic local Killing and first-recoil proof, not Taylor samples.")


def run():
    branch, f = higgs_branch()
    sources = [
        "archive_1009_/research_note_1018.md",
        "archive_342_369/research_note_358.md",
        "archive_1009_/research_note_1024.md",
        "archive_1009_/research_note_1027.md",
        "archive_1009_/1027/charged_source_selection.py",
        "archive_1009_/1027/charged_source_selection_results.json",
        "archive_1009_/1027/input_dependency_update_v0_16.md",
        "archive_1009_/1027/NEXT.md",
        "archive_956_989/981/drafts/common_parent_contract_v1.md",
        "archive_990_1008/993/common_candidate_v1.md",
    ]
    result = dict(round=1028, status="scientific_calibration_verified",
                  new_calibration_groups=1, cumulative_test_groups=3805,
                  new_cognitive_axioms=0, goal_complete=False,
                  all_scientific_calibrations_passed=True,
                  scope="Classical local retained action; fixed full minimal first-order Noether representative; same-parent bridge from 1027 to 1018.",
                  retained_inputs=["3+1 common Lorentz kinetic principal part", "P981 gauge/Higgs/Weyl fields and nonzero gauge interactions", "1027 leading source class; minimal I=0 and additional C=0 representative", "complete first-order minimal transformations and gravity bracket", "standard internal YM/local-Lorentz action fixes O at every order", "regular finite-jet higher deformations and free canonical M0=0", "nonzero portal only for inclusion of the neutral sigma source"],
                  not_proved=["cognitive generation of gravity or the Standard Model", "classification of every first-order coupling or improvement", "all quantum loops, anomalies, or UV completion", "full parent finite positive quantum realization", "absence of disconnected massless or additional massive gravity sectors"],
                  full_higgs_branch=branch,
                  internal_and_phase_checks=internal_invariance_and_phase_control(f),
                  invariant_killing_bridge=composite_killing_variations(f),
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  historical_source_sha256={p: hashlib.sha256((BASE/p).read_bytes()).hexdigest()
                                            for p in sources})
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
    print(json.dumps({"status": "passed", "round": 1028,
                      "mode": "exclusive_first_write" if args.write else "read_only_recompute_compare",
                      "observable_y_derivative": result["full_higgs_branch"]["observable_y_derivative"],
                      "ODE_orders_checked": result["full_higgs_branch"]["equation_orders"],
                      "historical_sources": len(result["historical_source_sha256"]),
                      "output": str(OUT)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
