"""Unnumbered finite split-gate screen, not the complete square-root parent.
Default: read-only recomputation and comparison. --write: exclusive first save.
Fraction proves exact identities/lower bounds; NumPy checks a 4x4 Gaussian shear.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from math import exp, sin, isclose
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT_CODE = HERE.parent / "finite_grid_candidate_after1256" / "diagnostic.py"
EXPECTED_PARENT_SHA = "e9bbcda917d9b9cebcffbc7c474320c418f952e50fd2e48034c3da3f97a4a759"
ABS_TOL = 1e-14
REL_TOL = 1e-12


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def transpose(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Q(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def serial(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: serial(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(item) for item in value]
    return value


def calculate():
    parent_sha = digest(PARENT_CODE)
    assert parent_sha == EXPECTED_PARENT_SHA, "Frozen parent provenance changed"
    dx, v, a = Q(2), Q(8), Q(1)
    delta, eps, d = Q(1, 16), Q(1, 1024), Q(1, 4)
    sq, sp = Q(1, 16), Q(1, 1024)
    hcoef = v * a / (8 * dx**2)
    commutator_scale = eps / v
    g = 2 * hcoef * delta / v
    assert hcoef == Q(1, 4) and g == Q(1, 256)
    assert sq * sp == commutator_scale / 2 == Q(1, 16384)
    # Ordering (Qs, Qt, Ps, Pt); P are density momenta, not cell momenta.
    S = [[Q(1), Q(0), Q(0), Q(0)],
         [Q(0), Q(1), Q(0), Q(0)],
         [-g, g, Q(1), Q(0)],
         [g, -g, Q(0), Q(1)]]
    J = [[Q(0), Q(0), Q(1), Q(0)],
         [Q(0), Q(0), Q(0), Q(1)],
         [-Q(1), Q(0), Q(0), Q(0)],
         [Q(0), -Q(1), Q(0), Q(0)]]
    assert mul(mul(S, J), transpose(S)) == J
    C = [[Q(0) for _ in range(4)] for _ in range(4)]
    for i, value in enumerate((sq**2, sq**2, sp**2, sp**2)):
        C[i][i] = value
    Cout = mul(mul(S, C), transpose(S))
    pure = mul(mul(Cout, J), Cout)
    expected_pure = [[commutator_scale**2 / 4 * x for x in row] for row in J]
    assert pure == expected_pure
    variance = sp**2 + 2 * g**2 * sq**2
    assert Cout[3][3] == variance == Q(9, 8 * 1024**2)
    omega, offset = Q(1024), Q(1, 2)
    exponent = omega**2 * variance / 2
    assert exponent == Q(9, 16)
    phases = [-offset, omega * g * d - offset]
    assert phases == [Q(-1, 2), Q(1, 2)]
    # Taylor signs are justified in README; exact arithmetic certifies margins.
    exp_lower = 1 - exponent + exponent**2 / 2 - exponent**3 / 6
    sin_lower = offset - offset**3 / 6
    gap_lower = exp_lower * sin_lower
    assert exp_lower == Q(4637, 8192) > 0
    assert sin_lower == Q(23, 48) > 0
    assert gap_lower == Q(106651, 393216) > Q(1, 4)
    error_lower = gap_lower / 2
    assert error_lower > Q(1, 8)
    gap_float = exp(-float(exponent)) * sin(float(offset))
    assert gap_float > float(gap_lower)
    probs = [(1 - gap_float) / 2, (1 + gap_float) / 2]

    # Independent float matrix calculation, not copied rational output arrays.
    Sf = np.eye(4)
    Sf[2, 0], Sf[2, 1] = -float(g), float(g)
    Sf[3, 0], Sf[3, 1] = float(g), -float(g)
    Jf = np.block([[np.zeros((2, 2)), np.eye(2)],
                   [-np.eye(2), np.zeros((2, 2))]])
    Cf = np.diag([float(sq)**2, float(sq)**2, float(sp)**2, float(sp)**2])
    Cout_f = Sf @ Cf @ Sf.T
    sympl_error = float(np.max(np.abs(Sf @ Jf @ Sf.T - Jf)))
    cov_error = float(np.max(np.abs(Cout_f - np.array(Cout, dtype=float))))
    assert sympl_error <= ABS_TOL and cov_error <= ABS_TOL
    kinetic_in = v * sp**2
    assert kinetic_in == Q(1, 131072)
    branches = []
    for bit, mean in enumerate((Q(0), d)):
        mu = [[mean], [Q(0)], [Q(0)], [Q(0)]]
        mu_out = [row[0] for row in mul(S, mu)]
        assert mu_out[3] == g * mean
        h_in = hcoef * (2 * sq**2 + mean**2)
        kin_delta = v * g**2 * (2 * sq**2 + mean**2)
        # Compute both quadratic expectations directly from output covariance.
        h_out = hcoef * (Cout[0][0] + Cout[1][1] - 2*Cout[0][1]
                            + (mu_out[0] - mu_out[1])**2)
        kinetic_out = v * sum((Cout[j][j] + mu_out[j]**2 for j in (2, 3)), Q(0)) / 2
        assert h_in == h_out
        assert kinetic_out - kinetic_in == kin_delta
        assert h_in == Q((1, 9)[bit], 512)
        assert kin_delta == Q((1, 9)[bit], 1048576)
        muf = Sf @ np.array([float(mean), 0.0, 0.0, 0.0])
        pf = (1 + exp(-float(omega)**2 * Cout_f[3, 3] / 2)
              * sin(float(omega)*muf[3]-float(offset))) / 2
        assert isclose(pf, probs[bit], rel_tol=REL_TOL, abs_tol=ABS_TOL)
        branches.append({"input_bit": bit, "source_Q_mean": mean,
                         "receiver_P_mean": mu_out[3], "output_mean": mu_out,
                         "receiver_effect_probability_float": float(pf),
                         "gradient_energy_in": h_in, "gradient_energy_out": h_out,
                         "aux_kinetic_in": kinetic_in, "aux_kinetic_out": kinetic_out,
                         "aux_kinetic_increase": kin_delta})
    separation = min(2*dx, 5*dx-2*dx)
    assert separation == 4 and delta < separation
    return serial({
        "status": "PASS_WITH_SCOPE",
        "kind": "unnumbered_frozen_gradient_subgate_screen",
        "new_round": None,
        "adopted_axioms": 0,
        "check_groups": 3,
        "scope": {
            "exact_object": "two-mode frozen a=1 gradient subgate, density CCR",
            "not_complete_square_root_parent": True,
            "not_original_two_program_K_menu": True,
            "actual_preparation_or_readout_not_certified": True,
            "Gaussian_amplitude_tails_not_parent_positive_domain": True,
            "finite_moments_not_full_apparatus_resource_closure": True,
            "no_whole_Trotter_error_claim": True},
        "provenance": {"check_sha256": digest(Path(__file__)),
                       "readme_sha256": digest(HERE / "README.md"),
                       "parent_diagnostic_sha256": parent_sha,
                       "parent_diagnostic_not_executed": True},
        "parameters": {"dx": dx, "cell_weight": v, "scale_a": a,
                       "simulation_parameter_delta": delta, "epsilon": eps,
                       "density_CCR_scale": commutator_scale,
                       "sigma_Q": sq, "sigma_P": sp, "source_displacement": d,
                       "gradient_coefficient": hcoef, "shear_g": g,
                       "effect_frequency": omega, "effect_offset": offset},
        "gaussian": {"uncertainty_product": sq*sp,
                     "Heisenberg_lower_bound": commutator_scale/2,
                     "exact_symplectic_identity": True,
                     "exact_pure_covariance_identity": True,
                     "output_covariance": Cout, "receiver_P_variance": variance,
                     "float_symplectic_residual": sympl_error,
                     "float_covariance_residual": cov_error},
        "probability_certificate": {
            "attenuation_exponent": exponent, "phases": phases,
            "exp_strict_lower": exp_lower, "sin_strict_lower": sin_lower,
            "probability_gap_strict_lower": gap_lower,
            "gap_margin_over_quarter": gap_lower-Q(1,4),
            "one_input_error_strict_lower": error_lower,
            "error_margin_over_eighth": error_lower-Q(1,8),
            "probability_gap_float": gap_float,
            "one_input_error_half_gap_float": gap_float/2,
            "no_communication_receiver_common_probability": "q",
            "bound": "max(|q-p0|,|q-p1|)>=gap/2>1/8"},
        "branches": branches,
        "conditional_causal_interpretation": {
            "frozen_metric_scale": "a=1, c=1",
            "point_factor_shortest_separation": separation,
            "no_communication_if_actual_duration": "tau < 4",
            "simulation_delta_not_actual_tau": True,
            "finite_mode_extent_and_source_backreaction_require_new_budget": True},
        "comparison_tolerance": {"absolute": ABS_TOL, "relative": REL_TOL}})


def compare(actual, saved, path="root"):
    if isinstance(actual, dict):
        assert isinstance(saved, dict) and actual.keys() == saved.keys(), path
        for key in actual:
            compare(actual[key], saved[key], path + "." + key)
    elif isinstance(actual, list):
        assert isinstance(saved, list) and len(actual) == len(saved), path
        for i, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(actual, float):
        assert isinstance(saved, (int, float)) and isclose(actual, saved,
            rel_tol=REL_TOL, abs_tol=ABS_TOL), (path, actual, saved)
    else:
        assert actual == saved, (path, actual, saved)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusive first results save")
    args = parser.parse_args()
    path = HERE / "results.json"
    if args.write and path.exists():
        raise FileExistsError("Refusing to replace existing results.json")
    if not args.write and not path.exists():
        raise FileNotFoundError("No saved result; use --write only for first authorized save")
    result = calculate()
    if args.write:
        with path.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
        action = "CREATED_ONCE"
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        action = "READ_ONLY_MATCH"
    cert = result["probability_certificate"]
    print(json.dumps({"status": result["status"], "action": action,
                      "check_groups": result["check_groups"],
                      "gap": cert["probability_gap_float"],
                      "strict_error_lower": cert["one_input_error_strict_lower"],
                      "scope": "frozen subgate only; no original K-task failure"},
                      ensure_ascii=False))


if __name__ == "__main__":
    main()
