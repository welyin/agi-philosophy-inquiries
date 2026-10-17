"""1267 geometry: one exact-arithmetic author group; no ODE reintegration.

Default only recomputes and compares the saved JSON structure.
--write exclusively creates geometry_results.json for the first time.
All floating values are labelled display; authoritative bounds are Fraction.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "research_cognition_physics").is_dir())
ARCH = HERE.parent
OUT = HERE / "geometry_results.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(rel):
    return json.loads((ARCH / rel).read_text(encoding="utf-8"))

def interval(values):
    lo, hi = map(F, values)
    assert lo <= hi
    return lo, hi

def maxabs(values):
    return max(map(abs, interval(values)))

def inverse_interval(values):
    lo, hi = interval(values)
    assert lo > 0
    return 1 / hi, 1 / lo

def calculate():
    old47 = load("1247/results.json")
    old49 = load("1249/quantum_bounds_results.json")
    old52 = load("1252/results.json")
    family = old52["group_parameter_family"]
    assert family["alpha"] == 1 and family["parameter_cells"] == 1
    assert family["exact_closed_interval_union"] is True
    assert list(map(F, family["kappa_interval"])) == [F(3), F(3001, 1000)]
    cell = family["cells"][0]
    assert list(map(F, cell["kappa"])) == [F(3), F(3001, 1000)]
    assert cell["pulse_junctions_exact"] is True
    # Check saved upstream bindings, not merely a status label.
    for rel, digest in old52["inputs_sha256"].items():
        assert sha(ARCH / rel) == digest, ("1252 input changed", rel)
    assert sha(ARCH / "1252/check.py") == old52["author_code_sha256"]

    tube = cell["all_time_state_tube"]
    ulo, uhi = interval(tube["u"])
    assert ulo > F(1, 3) and uhi < F(3, 2)
    assert F(tube["H"][1]) < F(11, 10)
    for key in ("p", "v", "X", "V"):
        assert maxabs(tube[key]) < F(3, 2)
    for key in ("q", "r", "psi", "Y", "K"):
        assert maxabs(tube[key]) < 1
    assert maxabs(tube["z"]) + F(1, 5) < 1
    # Q core: exp(1/2)>1+1/2 and exp(3/2)>1+3/2+(3/2)^2/2>3.
    assert 1 + F(3, 2) + F(3, 2)**2 / 2 > 3
    piQ_bound = F(6, 3) * 27 * F(11, 10)
    pi_scalar_bound = 27 * F(3, 2)
    assert piQ_bound < 61 and pi_scalar_bound < 61
    positive = cell["positive_all_time"]
    thresholds = {"J": F(2, 25), "M_c": F(9, 10),
                  "F": F(3, 5), "M_r": F(1, 25)}
    positive_slacks = {}
    for key, threshold in thresholds.items():
        lo = F(positive[key][0])
        assert lo > threshold
        positive_slacks[key] = str(lo - threshold)

    ledger = old49["group_constructive_finite_resource_bound"]["ledger"]
    Hwide = int(ledger["H0"])
    assert Hwide == 2**36
    added_gravity = F(1, 12000) * 129**2 * 3**12
    h_upper = F(ledger["h_complex_upper"]) + added_gravity
    h_slack = Hwide - h_upper
    assert h_upper == F(79598612877987, 4000)
    assert h_slack == F(195279294066013, 4000) > 0
    source_base = 3**18 * (1 + Hwide) * (1 + 129**2)
    source_upper = 1024 * source_base
    source_slack = 2**96 - source_upper
    assert source_slack == 78774463429498191360375764992 > 0
    # The only new kappa coefficients in h, H and its first equations stay
    # within the explicit coefficient guards used in geometry_proof section 4.
    kmin, kmax = F(3), F(3001, 1000)
    assert kmax / 6 < 1 and (kmax - 3) / 12 == F(1, 12000)

    Lcut = int(ledger["Lcut_decimal"])
    A_source = 2**98 * (1 + Lcut)**2
    L_source = 4 * (Lcut + 14)
    assert A_source >= 27 and L_source >= Lcut + 1
    EH = int(ledger["Fourier_H_exponent"])
    E0 = int(ledger["Fourier_menu_base_exponent"])
    EF = int(ledger["Fourier_transport_exponent_decimal"])
    N = int(ledger["integer_n_power_exponent_decimal"])
    assert ledger["M_power_of_two_exponent"] == 3059
    assert EH == 3093 and E0 == 14008
    assert EF == 95928 + 100 * 2**3059
    assert N == 1 + max((EF + 44)//3, (EH + EF + 47)//6)
    linear_margin, quadratic_margin = 4*N - EF, 8*N - EH - EF
    assert linear_margin > 42 and quadratic_margin > 42
    # Do not expand 2^N or the physical resource s.
    eq = F(1, 2**40)
    assert F(1, 2**41) < eq

    eta, delta_ref = F(1, 10**8), F(1, 10**16)
    bin_width, tau, beta_sys = F(1, 10**10), F(5, 10**9), F(4, 10**9)
    mean_error_upper = tau + beta_sys + bin_width + eq
    assert mean_error_upper < eta
    repeats = 26 * A_source**2 * 200000000**2
    assert F(repeats) * tau**2 / (2*A_source**2) == 13
    four_task_failure_upper = 8 * F(3, 8)**13
    assert four_task_failure_upper < F(1, 100)
    outputs = 2 * A_source * 10**10 + 1

    intervals_a = []
    for branch in ("alpha0", "alpha1"):
        endpoint = old47["group_finite_memory"][branch]["checkpoints"]["2"]["state"]
        intervals_a.append(inverse_interval(endpoint["u"]))
    program_gap_lower = intervals_a[0][0] - intervals_a[1][1]
    program_gap_upper = intervals_a[0][1] - intervals_a[1][0]
    assert program_gap_lower > F(683, 10**7)
    program_after_error = F(683, 10**7) - 2*eta
    assert program_after_error > F(68, 10**6)

    final_state = cell["checkpoints"]["2"]["state"]
    final_sens = cell["checkpoints"]["2"]["sensitivity"]
    u2lo, u2hi = interval(final_state["u"])
    sulo, suhi = interval(final_sens["u"])
    assert u2lo > 0 and suhi < 0
    derivative_lower = -suhi / u2hi**2
    assert derivative_lower > F(1, 5)
    ma = F(1, 5)
    endpoint_gap = ma * (kmax - kmin)
    endpoint_gap_after_read = endpoint_gap - 2*eta
    assert endpoint_gap == F(1, 5000)
    assert endpoint_gap_after_read == F(9999, 50000000) > 0

    ref_lo, ref_hi = intervals_a[1]
    ref_mid = (ref_lo + ref_hi)/2
    ref_radius = (ref_hi - ref_lo)/2
    assert ref_radius < delta_ref
    param_bound = (2*eta + delta_ref) / ma
    field_bound = F(111, 100) * param_bound
    assert param_bound == F(200000001, 2000000000000000)
    assert field_bound == F(22200000111, 200000000000000000)

    # Analytic continuation guards for the CLASSICAL hold, not a claimed
    # new post-t=2 quantum tube. No requirement u>1/3 on this larger family.
    rho0 = F(111, 100)
    assert rho0 < F(9, 8)
    assert kmax*rho0/3 < F(11, 10)**2
    hold_duration = F(1, 4)
    hold_u_ratio = 1 - F(11, 10)*hold_duration
    assert hold_u_ratio == F(29, 40) > 0
    hold_F_lower = 1 - F(3, 2)/5
    hold_Mr_lower = F(1, 5) - F(3, 2)/10
    hold_dust_extra_lower = F(1, 50) - F(3, 2)/100
    assert hold_F_lower == F(7, 10)
    assert hold_Mr_lower == F(1, 20)
    assert hold_dust_extra_lower == F(1, 200)
    # Original h at t=0 with constraint-compatible pi_Q: grav=-rho0.
    material_h0 = F(1, 200) + F(3, 20) + F(71, 100) + F(1, 8) + F(1, 50)
    h0 = -rho0 + material_h0
    assert h0 == F(-1, 10)
    wrong_fixed_H0_constraint = rho0 - kmax*rho0/3
    assert wrong_fixed_H0_constraint == F(-37, 100000)
    base_energy_coefficient = F(111, 100) * F(1, 10**6)
    assert base_energy_coefficient == F(111, 10**8)

    dependencies = [
        "1245/proof.md", "1247/proof.md", "1247/check.py", "1247/results.json",
        "1248/proof.md",
        "1249/canonical.md", "1249/quantum_bridge.md", "1249/quantum_bounds.py",
        "1249/quantum_bounds_results.json",
        "1251/proof.md", "1251/check.py", "1251/results.json",
        "1252/proof.md", "1252/check.py", "1252/results.json",
        "1266/proof.md",
        "_admission/causal_block_after1265/proof.md",
        "_admission/causal_block_after1265/results.json",
        "_admission/finite_common_task_after1266/README.md",
    ]
    bindings = {rel: sha(ARCH / rel) for rel in dependencies}
    bindings["geometry_proof.md"] = sha(HERE / "geometry_proof.md")
    bindings["geometry_check.py"] = sha(Path(__file__))
    return {
        "round": 1267,
        "part": "fixed_geometry_uniform_quantum_parameter_bridge",
        "status": "PASS_AUTHOR_WITH_SCOPE",
        "author_groups": 1,
        "arithmetic": "exact integers and Fraction; no ODE rerun or parameter sampling",
        "same_object": {
            "canonical_pairs": 7,
            "kappa_interval": ["3", "3001/1000"],
            "parameter_program": "B=2",
            "read_time": "2",
            "gate_relation": "terminal prefix before the 1266 flag gate",
            "fixed_parameter_independent_menu": "a_c=chi*exp(Q), fixed binned marginal",
            "not_asserted": "the full (a,rho,P) binning is not parameter-independent",
            "spatial_period": "1/100",
            "base_coordinate_volume": "1/1000000",
            "resource_s": "1000*ceil(max(1,hbar))*2^(2*N)",
            "epsilon_upper": "2^(-4*N)",
            "native_block_and_systematic_matching": "explicit conditional C_block inputs",
        },
        "uniform_tube_and_bounds": {
            "closed_parameter_cells": 1,
            "tube_time": ["0", "2"],
            "saved_u_tube": tube["u"],
            "saved_H_tube": tube["H"],
            "pi_Q_absolute_upper": str(piQ_bound),
            "pi_scalar_absolute_upper": str(pi_scalar_bound),
            "positive_threshold_slacks": positive_slacks,
            "h_complex_upper": str(h_upper),
            "h_Hwide_slack": str(h_slack),
            "source_complex_upper": str(source_upper),
            "source_2power96_slack": str(source_slack),
            "Mflow_exponent": 3059,
            "EH": EH, "E0": E0, "EF_decimal": str(EF), "N_decimal": str(N),
            "linear_negative_exponent": str(linear_margin),
            "quadratic_negative_exponent": str(quadratic_margin),
            "uniform_menu_error_strict_upper": "2^-41",
            "conservative_menu_error": str(eq),
        },
        "finite_fixed_readout": {
            "bin_width": str(bin_width), "tau": str(tau), "beta_sys_input": str(beta_sys),
            "range_A_source_decimal": str(A_source),
            "repeats_per_task_decimal": str(repeats),
            "finite_output_count_decimal": str(outputs),
            "mean_error_upper": str(mean_error_upper),
            "eta": str(eta),
            "at_most_four_task_failure_upper": str(four_task_failure_upper),
            "coverage_scope": "a finite selected family of tasks, not simultaneous data for all kappa",
        },
        "geometry_and_inverse": {
            "B_program_a_intervals": [[str(x) for x in pair] for pair in intervals_a],
            "B_program_gap_exact_lower": str(program_gap_lower),
            "B_program_gap_exact_upper": str(program_gap_upper),
            "B_program_gap_certified_strict_lower": "683/10000000",
            "B_program_gap_after_read_strict_lower": str(program_after_error),
            "uniform_u2": final_state["u"],
            "uniform_du2_dkappa": final_sens["u"],
            "uniform_da_dkappa_exact_lower": str(derivative_lower),
            "uniform_da_dkappa_coarse_strict_lower": str(ma),
            "kappa_endpoint_gap_strict_lower": str(endpoint_gap),
            "kappa_endpoint_after_read_strict_lower": str(endpoint_gap_after_read),
            "reference_a_midpoint": str(ref_mid),
            "reference_a_exact_radius": str(ref_radius),
            "reference_prediction_error": str(delta_ref),
            "acceptance_rule": "abs(observed_a-reference_a_midpoint)<=eta",
            "parameter_error_bound": str(param_bound),
            "classical_orthonormal_field_residual_bound": str(field_bound),
        },
        "classical_hold_and_initial_constraint": {
            "h_initial_for_all_kappa": str(h0),
            "fixed_old_H0_bad_initial_constraint": str(wrong_fixed_H0_constraint),
            "classical_field_time_window": ["0", "9/4"],
            "source_component_bound": str(rho0),
            "hold_u_ratio_lower": str(hold_u_ratio),
            "hold_F_lower": str(hold_F_lower),
            "hold_Mr_lower": str(hold_Mr_lower),
            "hold_dust_extra_lower": str(hold_dust_extra_lower),
            "initial_material_inventory": "111*s/10^8",
            "physical_volume": "10^-6*s^3",
            "no_post_t2_quantum_extension_claim": True,
        },
        "display_only": {
            "B_gap_lower": float(program_gap_lower),
            "uniform_da_dkappa_lower": float(derivative_lower),
            "mean_error_upper": float(mean_error_upper),
            "four_task_failure_upper": float(four_task_failure_upper),
            "parameter_error_bound": float(param_bound),
            "classical_field_residual_bound": float(field_bound),
        },
        "scope": [
            "Known Gaussian preparation family; no all-input classical geometry claim.",
            "The actual fixed a marginal and systematic/source/clock matching are explicit inputs.",
            "No old local numerical K precision, no parameter-family flag-classification inheritance.",
            "The field residual is for the same declared classical comparison family through 9/4.",
            "No quantum Einstein operator equation, full GR, arbitrary-space-time inverse, or experiment.",
            "CT not adopted; finite resource sufficiency is not engineering feasibility.",
        ],
        "bindings_sha256": bindings,
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true",
                        help="First exclusive result creation; never overwrite.")
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("FIRST_EXCLUSIVE_WRITE PASS_AUTHOR_WITH_SCOPE; author_groups=1")
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        if saved != result:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS_AUTHOR_WITH_SCOPE; author_groups=1")
    print(json.dumps(result["display_only"], ensure_ascii=False))

if __name__ == "__main__":
    main()

