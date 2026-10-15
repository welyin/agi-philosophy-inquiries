"""1232 author certificate: exact rational guards plus finite diagnostics.

This file does not compute a Dirac spectral density or build an accelerator.
The atomic regression is a separate finite KMS measure, not Dirac data.
Default: read-only replay. --write-results: create a new receipt, never overwrite.
"""
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / "results.json"
SOURCES = (
    "archive_1219_/1232/proof.md",
    "archive_1219_/_admission/stage_scope_audit_after1230.md",
    "archive_1219_/_admission/after1230_finite_motion_thermal_admission.md",
    "archive_1219_/1222/proof.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1229/proof.md",
    "archive_1106_1119/1119/proof.md",
    "archive_1120_1158/1126/proof.md",
    "archive_301_341/research_note_312.md",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rat(value: F) -> dict:
    value = F(value)
    return {"exact": str(value), "display": float(value)}


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def response_interval(atoms: dict[F, F], energy: F, duration: F) -> tuple[F, F]:
    """A rigorous interval: retain the central atom and bound all sin^4 by one."""
    low = atoms.get(energy, F(0)) * duration**2
    high = low
    for nu, mass in atoms.items():
        offset = energy - nu
        if offset:
            high += mass * F(16) / (duration**2 * offset**4)
    return low, high


def finite_response_display(atoms: dict[F, F], energy: F, duration: F) -> float:
    ans = 0.0
    for nu, mass in atoms.items():
        z = float((energy - nu) * duration / 2)
        sinc = 1.0 if z == 0.0 else math.sin(z) / z
        ans += float(mass * duration**2) * sinc**4
    return ans


def calculate() -> dict:
    T, r, mass_core = F(128), F(1, 4), F(1, 64)
    coupling, op_error, rho = F(1, 65536), F(1, 10**10), F(1, 100)
    c0 = F(23, 24)**4
    A = c0 * T**2 * mass_core
    tail = F(16) / (r**4 * T**2)
    u = tail / A
    x = coupling * T
    remainder = x**3 / (6 * (1 - x*x/20))
    dyson = 2*x*remainder + remainder**2
    qmin = coupling**2 * A

    # Group 1: exact guards in the analytic proof.
    exp07_lower = sum((F(7, 10)**j / math.factorial(j) for j in range(5)), F(0))
    ensure(exp07_lower > 2, "finite positive exp series must prove log(2)<0.7")
    ensure(F(1) - F(1, 2)**2 / 6 == F(23, 24), "sinc core guard")
    r3_coefficient = F(1) / (6 * (1 - F(1, 20)))
    dyson_coefficient = 2*r3_coefficient + r3_coefficient**2
    ensure(dyson_coefficient < F(1, 2), "uniform d_D <= x^4/2 for x<=1")
    ensure(F(1)/T <= r and x <= 1, "core and series domains")
    # These polynomial integrals give int chi=T and int chi^2=2T/3.
    tent_L1 = 2*(T - T/2)
    tent_L2_sq = 2*(T - T + T/3)
    ensure(tent_L1 == T and tent_L2_sq == 2*T/3, "tent normalization")
    g1 = {
        "passed": True,
        "method": "exact rational inequalities; Fourier identity proved in proof.md",
        "tent_half_duration": rat(T),
        "tent_full_duration": rat(2*T),
        "tent_L1": rat(tent_L1),
        "tent_L2_squared": rat(tent_L2_sq),
        "sinc_core_factor": rat(c0),
        "exp_0p7_positive_partial_sum": rat(exp07_lower),
        "uniform_dyson_x4_coefficient": rat(dyson_coefficient),
        "coarse_x4_coefficient": rat(F(1, 2)),
    }

    # Group 2: finite KMS measure fixture, never identified with the Dirac field.
    atoms = {F(1): F(1, 32), F(-1): F(1, 16),
             F(2): F(1, 64), F(-2): F(1, 16)}
    total = sum(atoms.values(), F(0))
    ensure(total <= 1, "finite-measure total mass")
    bands = []
    for energy in (F(1), F(2)):
        target = F(2)**int(energy)
        ensure(atoms[-energy] == target*atoms[energy], "KMS reflection convention")
        core = sum((m for nu, m in atoms.items()
                    if abs(nu-energy) <= 1/T), F(0))
        ensure(core >= mass_core, "fixture core mass")
        plo, phi = response_interval(atoms, energy, T)
        mlo, mhi = response_interval(atoms, -energy, T)
        ratio_low, ratio_high = mlo/phi, mhi/plo
        ensure(ratio_low <= target <= ratio_high, "response ratio interval")
        # log y <= y-1 makes this an exact rational log-error upper bound.
        fixture_log_error = max(target/ratio_low-1, ratio_high/target-1, F(0))
        generic_spectral_error = F(7, 10)*r + u
        ensure(fixture_log_error < generic_spectral_error, "band envelope")
        pnum = finite_response_display(atoms, energy, T)
        mnum = finite_response_display(atoms, -energy, T)
        ensure(float(plo)-1e-10 <= pnum <= float(phi)+1e-10,
               "diagnostic forward response inside exact interval")
        ensure(float(mlo)-1e-10 <= mnum <= float(mhi)+1e-10,
               "diagnostic reverse response inside exact interval")
        bands.append({
            "energy": rat(energy), "target_ratio": rat(target),
            "core_mass": rat(core),
            "forward_response_interval": [rat(plo), rat(phi)],
            "reverse_response_interval": [rat(mlo), rat(mhi)],
            "ratio_interval": [rat(ratio_low), rat(ratio_high)],
            "log_error_rational_upper": rat(fixture_log_error),
            "response_display_only": {"forward": pnum, "reverse": mnum,
                                      "log_ratio": math.log(mnum/pnum)},
        })
    g2 = {
        "passed": True,
        "fixture_is_Dirac_spectrum": False,
        "fixture_type": "four-atom finite positive KMS measure at beta=log(2)",
        "atoms": [{"frequency": rat(nu), "mass": rat(m)}
                  for nu, m in sorted(atoms.items())],
        "total_mass": rat(total), "bands": bands,
        "exact_interval_method": "central atom exact; each other sin^4 bounded by one",
    }

    # Group 3: common, strictly positive finite certificate.
    total_error = dyson + op_error
    allowed = rho*qmin
    log_spectral_upper = F(7, 10)*r + u
    log_probability_upper = 2*rho/(1-rho)
    log_total_upper = log_spectral_upper + log_probability_upper
    ensure(total_error < allowed, "actual/full-Dyson joint budget")
    ensure(log_total_upper < F(1, 5), "two-gap log error < 1/5")
    ensure(qmin > 0 and op_error > 0 and coupling > 0, "nonzero finite window")
    g3 = {
        "passed": True,
        "conditional_input": "each actual positive-frequency core mass is at least 1/64",
        "actual_Dirac_core_mass_computed": False,
        "gaps": [1, 2], "beta": "log(2)", "a": "2*pi/log(2)",
        "T_half": rat(T), "r": rat(r), "m_star": rat(mass_core),
        "coupling": rat(coupling), "operator_norm_guard": 1,
        "rho": rat(rho), "op_error": rat(op_error),
        "core_response_lower_A": rat(A),
        "both_tail_integrals_upper": rat(tail), "tail_to_core_u": rat(u),
        "x": rat(x), "R3_rational_upper": rat(remainder),
        "dyson_probability_error_upper": rat(dyson),
        "leading_probability_lower": rat(qmin),
        "actual_probability_lower": rat(qmin-total_error),
        "joint_probability_error_upper": rat(total_error),
        "allowed_relative_error": rat(allowed),
        "joint_budget_strict_margin": rat(allowed-total_error),
        "spectral_log_error_upper": rat(log_spectral_upper),
        "probability_log_error_upper": rat(log_probability_upper),
        "total_log_error_upper": rat(log_total_upper),
        "strict_margin_to_one_fifth": rat(F(1, 5)-log_total_upper),
    }

    # Group 4: exactly soluble 2x2 blocks of the full two-qubit unitary.
    # This is a separate constant-switch diagnostic, not Dirac or the tent task.
    diagnostics = []
    for lam in (0.01, 0.05, 0.1):
        for duration in (0.25, 1.0):
            c = lam/2.0  # B=sigma_x/2
            resonant = math.sin(c*duration)**2
            off = c*c/(1+c*c)*math.sin(math.sqrt(1+c*c)*duration)**2
            p_up = (2*off+resonant)/3
            p_down = (2*resonant+off)/3
            q_up = c*c*(2*math.sin(duration)**2+duration**2)/3
            q_down = c*c*(2*duration**2+math.sin(duration)**2)/3
            xx = lam*duration  # uses the same conservative ||B||<=1 guard
            rr = xx**3/(6*(1-xx*xx/20))
            dd = 2*xx*rr+rr*rr
            err_up, err_down = abs(p_up-q_up), abs(p_down-q_down)
            ensure(max(err_up, err_down) <= dd+1e-16, "full-unitary Dyson diagnostic")
            diagnostics.append({
                "lambda": lam, "duration": duration,
                "exact_closed_form_up": p_up, "exact_closed_form_down": p_down,
                "leading_up": q_up, "leading_down": q_down,
                "max_probability_difference": max(err_up, err_down),
                "dyson_guard_display": dd,
            })
    g4 = {
        "passed": True, "certification": "analytic two-block formula; floats are diagnostics",
        "system": "two qubits, H0=diag(0,1)+diag(0,1), B=sigma_x/2",
        "field_initial_probabilities": ["2/3", "1/3"],
        "switch": "constant on one interval; not the Dirac/tent example",
        "cases": diagnostics,
    }

    # Group 5: deletions with exact rational witnesses where relevant.
    theta = F(1, 10)
    total_x = 2*theta
    error_real_lower = (theta-theta**3/6)**2
    wrong_odd_bound = total_x**3/(6*(1-total_x**2/20))
    ensure(error_real_lower > wrong_odd_bound, "diagonal-then-flip kills odd-only remainder")
    small = F(1, 10**8)
    ensure(abs(small-2*small) <= 2*small, "absolute error deletion")
    g5 = {
        "passed": True,
        "no_positive_core": {
            "measure": "delta_0", "KMS_for_any_beta": True,
            "Fminus_over_Fplus_for_real_tent_when_nonzero": 1,
            "positive_core_mass_at_E1_rquarter": 0,
        },
        "no_relative_probability_control": {
            "leading_plus": rat(small), "leading_minus": rat(2*small),
            "allowed_absolute_error": rat(2*small),
            "permitted_observed_plus": rat(small),
            "permitted_observed_minus": rat(small), "observed_ratio": 1,
        },
        "wrong_clock_scale": {
            "target_beta": "log(2)", "transport_rate_multiplier": 2,
            "actual_spectral_beta": "log(2)/2", "gap": 2,
            "target_ratio": 4, "actual_ideal_ratio": 2,
            "status": "analytical time-reparametrization deletion, not an actual device",
        },
        "diagonal_then_flip": {
            "unitary": "exp(-i theta X) exp(-i theta Z), input |0>",
            "theta": rat(theta),
            "flip_amplitude": "-sin(theta)^2-i sin(theta)cos(theta)",
            "first_order_flip_amplitude": "-i theta",
            "amplitude_error_rigorous_lower": rat(error_real_lower),
            "inapplicable_odd_only_remainder_upper": rat(wrong_odd_bound),
            "strict_violation_margin": rat(error_real_lower-wrong_odd_bound),
        },
    }

    source_hashes = {name: sha(ROOT/name) for name in SOURCES}
    groups = {
        "tent_and_series_guards": g1,
        "finite_measure_band_regression": g2,
        "joint_finite_certificate": g3,
        "full_finite_coupling_diagnostic": g4,
        "deletion_tests": g5,
    }
    ensure(len(groups) == 5 and all(g["passed"] for g in groups.values()), "all groups")
    return {
        "round": 1232, "status": "PASS", "all_checks_passed": True,
        "check_groups": len(groups), "groups": groups,
        "scope": {
            "full_probability_conditional_theorem": True,
            "Dirac_spectral_nonempty_existence": "analytic proof, not computed by this code",
            "fixed_table_is_conditional": True,
            "fixed_table_Dirac_core_mass_certified": False,
            "atomic_fixture_is_Dirac_data": False,
            "actual_motion_controller_constructed": False,
            "actual_operation_error_achieved": False,
            "global_vacuum_called_finite_temperature_bath": False,
            "TG2_or_stage_completion_signed": False,
        },
        "code_sha256": sha(Path(__file__)),
        "source_sha256": source_hashes,
    }


def main() -> None:
    args = sys.argv[1:]
    if args not in ([], ["--write-results"]):
        raise SystemExit("usage: check.py [--write-results]")
    result = calculate()
    if args:
        with RECEIPT.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        read_only = False
    else:
        saved = json.loads(RECEIPT.read_text(encoding="utf-8"))
        ensure(result == saved, "saved result differs from read-only replay")
        read_only = True
    print(json.dumps({
        "round": 1232, "status": "PASS", "all_checks_passed": True,
        "results_match": True, "check_groups": 5, "read_only": read_only
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
