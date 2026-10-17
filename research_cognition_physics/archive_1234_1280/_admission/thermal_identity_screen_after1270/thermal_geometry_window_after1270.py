"""Non-numbered finite geometry/window screen; no field, apparatus or source simulation.

Default: read-only exact replay. --write-results: exclusive first save.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve()
ROOT = next((p for p in HERE.parents if
             (p / "research_cognition_physics").is_dir() and
             (p / "scripts").is_dir()), None)
if ROOT is None:
    raise RuntimeError("Cannot locate repository root from script parents.")
OUT = HERE.with_suffix(".json")


def rat(x: F) -> dict:
    return {"numerator": str(x.numerator), "denominator": str(x.denominator)}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def calculate() -> dict:
    # e<3, 3<pi<22/7; exact rational upper bounds only.
    q, A = 1024, 2**1800
    zeta = F(1, 16)
    T, r_packet = F(q, A), zeta / A
    Hmax, amin = F(3, 10), F(1)
    period, half_time = F(1, 100), F(1, 8)
    exp_q1_upper = 3 ** (q + 1)
    # Bootstrap for actual constant-proper-acceleration curve on comparison FLRW.
    rapidity_error = Hmax * q * exp_q1_upper / A
    tube_half_time = (q + zeta) * exp_q1_upper / A
    tube_half_comoving_space = tube_half_time / amin
    assert rapidity_error < 1
    assert tube_half_time < half_time
    assert tube_half_comoving_space < period / 8
    assert 2 * A * r_packet < 1  # Flat reference base box strictly within wedge.
    assert tube_half_time < F(1, 2**160)

    # Full old sufficient thermal certificate is stronger than beta/T < eta.
    eta = F(1, 5)
    c0 = F(23, 24)**4
    b, mstar, rho = F(16), F(1, 64), F(1, 100)
    energy_times_T = F(128)
    assert b >= 1 and energy_times_T > b
    beta_r_upper = F(44, 7) * b / q
    u = F(16) / (c0 * mstar * b**4)
    spec_upper = beta_r_upper + u  # log(1+u) <= u
    log_transport_upper = 2 * rho / (1 - rho)
    log_total_upper = spec_upper + log_transport_upper
    assert log_total_upper < eta
    x = F(1, 512)  # |coupling| T, not the geometric A.
    q_leading_lower = x*x*c0*mstar
    d_dyson_upper = x**4/2
    op_error_input = F(1, 10**10)
    assert d_dyson_upper + op_error_input < rho*q_leading_lower
    true_probability_lower = (1-rho)*q_leading_lower
    assert true_probability_lower > F(1, 10**8)
    # beta E = pi/4 > 3/4; total error <1/5 leaves a nonzero thermal contrast.
    log_ratio_lower = F(3, 4) - log_total_upper
    assert log_ratio_lower > F(1, 2)

    # q=32 passes ONLY the necessary beta/T test but fails this full certificate.
    assert F(44, 7*32) < eta
    # If beta*r<eta with pi>3, b<eta*q/6=16/15.
    b_upper_q32 = eta*32/6
    assert b_upper_q32 == F(16, 15)
    # mstar<=1, c0<1; hence u>16/b_upper^4>12.
    u_lower_q32 = F(16)/b_upper_q32**4
    assert u_lower_q32 > 12
    # log(1+u)>eta: exp(eta)<1/(1-eta)=5/4<1+u.
    exp_eta_upper = 1/(1-eta)
    assert 1+u_lower_q32 > exp_eta_upper

    # Explicit scaling identities for arbitrary positive scale factor, sample rational
    # only a regression of the algebra in the memo, not proof for all factors.
    s = F(13, 7)
    beta_without_2pi = F(1, A)
    assert (s*beta_without_2pi)/(s*T) == beta_without_2pi/T
    assert (A/s)*(s*T) == A*T
    assert (s*tube_half_time)/(s*period) == tube_half_time/period

    source_files = {
        "1232_finite_thermal_proof": ROOT/"research_cognition_physics/archive_1219_1233/1232/proof.md",
        "1256_window_screen": ROOT/"research_cognition_physics/archive_1234_/_admission/thermal_window_screen_after1256.md",
        "1258_joint_screen": ROOT/"research_cognition_physics/archive_1234_/_admission/priority_and_joint_menu_after1258.md",
        "1270_continuation": ROOT/"research_cognition_physics/archive_1234_/1270/continuation_proof.md",
        "1270_output": ROOT/"research_cognition_physics/archive_1234_/1270/output_proof.md",
    }
    return {
        "status": "PASS_UNNUMBERED_GEOMETRY_AND_CONDITIONAL_BUDGET_ONLY",
        "new_scientific_rounds": 0,
        "source_bindings": {k: sha(v) for k, v in source_files.items()},
        "check_sha256": sha(HERE),
        "parameters": {
            "q_A_times_T": q, "A": "2^1800", "T": "2^-1790",
            "base_packet_halfwidth": "2^-1804",
            "period": rat(period), "H_upper": rat(Hmax), "a_lower": rat(amin),
            "b_frequency_halfwidth_times_T": rat(b),
            "energy_times_T": rat(energy_times_T),
            "mstar_conditional_input": rat(mstar),
            "coupling_times_T": rat(x), "op_error_conditional_input": rat(op_error_input),
        },
        "exact_geometry_bounds": {
            "rapidity_error_upper": rat(rapidity_error),
            "whole_transported_tube_half_time_upper": rat(tube_half_time),
            "whole_transported_tube_half_comoving_space_upper": rat(tube_half_comoving_space),
            "strong_display_upper": "2^-160",
            "flat_base_wedge_margin": rat(1-2*A*r_packet),
        },
        "exact_conditional_thermal_bounds": {
            "c0": rat(c0), "tail_ratio_upper": rat(u),
            "log_total_upper": rat(log_total_upper),
            "dyson_probability_error_upper": rat(d_dyson_upper),
            "leading_probability_lower": rat(q_leading_lower),
            "actual_probability_lower_if_inputs": rat(true_probability_lower),
            "log_ratio_lower_if_inputs": rat(log_ratio_lower),
        },
        "q32_certificate_failure": {
            "b_required_upper": rat(b_upper_q32),
            "tail_ratio_strict_lower": rat(u_lower_q32),
            "statement": "Fails old sufficient bound, not a physical thermometry impossibility."
        },
        "display_only": {
            "log_total_upper": float(log_total_upper),
            "positive_probability_lower_if_inputs": float(true_probability_lower),
            "log_ratio_lower_if_inputs": float(log_ratio_lower),
        },
        "unproved": [
            "Actual small-block reference supplies mstar and KMS comparison.",
            "Curved-reference covariance/CCR and clock/implementation match op_error.",
            "Acceleration, gap, coupling, preparation, output and backreaction have same-source budget.",
            "The selected high finite frequencies lie in an admitted effective material domain.",
            "No Wald entropy identity or original 1250 cap transplantation."
        ]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_results:
        with OUT.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
        mode = "EXCLUSIVE_FIRST_WRITE"
    else:
        assert json.loads(OUT.read_text(encoding="utf-8")) == result
        mode = "READ_ONLY_MATCH"
    print(json.dumps({"mode": mode, "status": result["status"],
                      "display_only": result["display_only"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()

