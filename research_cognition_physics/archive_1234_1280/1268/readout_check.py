"""1268 readout: one exact-arithmetic author group, no ODE reintegration.
Default recomputes and compares JSON. --write creates the first result only.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "research_cognition_physics").is_dir())
ARCH = HERE.parent
OUT = HERE / "readout_results.json"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load(p):
    return json.loads(p.read_text(encoding="utf-8"))

def calculate():
    geo = load(ARCH / "1267/geometry_results.json")
    for rel, digest in geo["bindings_sha256"].items():
        p = ARCH / ("1267/" + rel if rel.startswith("geometry_") else rel)
        assert sha(p) == digest, ("frozen upstream changed", rel)
    qb = load(ARCH / "1249/quantum_bounds_results.json")
    ledger = qb["group_constructive_finite_resource_bound"]["ledger"]
    G = geo["geometry_and_inverse"]
    U = geo["uniform_tube_and_bounds"]
    assert geo["author_groups"] == 1
    assert geo["same_object"]["kappa_interval"] == ["3", "3001/1000"]
    assert U["tube_time"] == ["0", "2"]
    assert F(U["saved_u_tube"][0]) > F(1, 3)
    assert F(U["saved_H_tube"][1]) < F(11, 10)
    # Actual frozen endpoints, including the whole closed parameter box.
    a_intervals = [list(map(F, v)) for v in G["B_program_a_intervals"]]
    ulo, uhi = map(F, G["uniform_u2"])
    a_intervals.append([1/uhi, 1/ulo])
    assert all(F(5, 2) < lo <= hi < F(29, 10) for lo, hi in a_intervals)
    assert F(G["B_program_gap_certified_strict_lower"]) == F(683, 10**7)
    assert F(G["uniform_da_dkappa_coarse_strict_lower"]) == F(1, 5)
    assert F(G["kappa_endpoint_gap_strict_lower"]) == F(1, 5000)

    # Same finite physical resource, never expand 2**N.
    M = 2**3059
    EH, E0 = 3093, 14008
    EF = 95928 + 100*M
    N = 1 + max((EF+44)//3, (EH+EF+47)//6)
    assert N == int(U["N_decimal"]) == int(ledger["integer_n_power_exponent_decimal"])
    B = 324 * (factorial(81)//factorial(40)**2) * 2**40
    assert B == int(ledger["B_decimal"])
    LF = 27*(1+B)
    assert B < 2**132 and LF < 2**136
    # Shared Gaussian residual and its Ramsey extension.
    exponents = {
        "free_Gaussian_norm": EH+15*M+16,
        "Ramsey_Duhamel_norm": EH+15*M+20,
        "full_Gaussian_comparison": EH+15*M+21,
        "bounded_old_source_norm": E0+EH+15*M+25,
    }
    assert all(2*N - x > 110 for x in exponents.values())
    assert 4*N-E0 > 110
    assert 4*N-(3*M+136+E0+4) > 110

    tau, alpha = F(1, 10**8), F(1, 16)
    # Exact normalization and maximum of 140u^3(1-u)^3.
    integral = sum(F(140*((-1)**j)*factorial(3), factorial(j)*factorial(3-j)*(4+j))
                   for j in range(4))
    assert integral == 1
    pulse_peak = F(140, 64)/tau
    assert pulse_peak == F(35, 16)/tau == 218750000
    assert F(5, 2)*(1-F(11, 10)*tau) > 2
    assert F(11, 10)*3 < 4
    cosine_lower = 1-F(3, 8)**2/2
    slope_lower = F(1, 32)
    assert cosine_lower/16 > slope_lower
    time_prob_error = tau/4

    # Full-space bounded source comparison completion, not a covariant claim.
    source_peak_coefficient = alpha*pulse_peak*max(162, 27*(1+B))
    assert source_peak_coefficient < 2**161
    # Clock compensation uses the inherited Weyl source completion,
    # not a false pointwise-sup bound for a phase-space Weyl symbol.
    Lcut = 56*B
    CF = 2**34*(5*256**13)*15**10*factorial(20)
    clock_Weyl_F0 = CF*81*(2*Lcut+2)**20
    assert clock_Weyl_F0 < 2**E0
    complete_source_peak = alpha*pulse_peak*(162+27*(1+B)+2**E0)
    assert complete_source_peak < 2**(E0+26)
    assert 4*N-(E0+26) > 110
    gate_peak = alpha*pulse_peak*27
    assert gate_peak < 2**29
    # Finite statistical and actual matching budgets.
    dR = F(1, 2**100)
    beta, stat = F(1, 10**7), F(1, 10**7)
    R = 5*10**14
    assert 2*R*stat*stat == 10
    four_failure = 8*F(3, 8)**10
    assert four_failure < F(1, 100)
    ep = beta+stat+time_prob_error+dR
    eta_a = F(1, 150000)
    assert 32*ep < eta_a
    B_ideal_gap = F(683, 32*10**7)
    K_ideal_gap = F(1, 160000)
    B_gap = B_ideal_gap-2*ep
    K_gap = K_ideal_gap-2*ep
    assert B_gap > F(17, 10**7)
    assert K_gap > F(58, 10**7)
    ref = F(1, 10**16)
    kappa_accept = 5*(2*eta_a+ref)
    field_prefix = F(111, 100)*kappa_accept
    assert kappa_accept == F(1, 15000)+F(5, 10**16)
    assert field_prefix < F(75, 10**6)

    paths = [
        "1247/proof.md", "1247/results.json",
        "1249/canonical.md", "1249/quantum_bridge.md",
        "1249/quantum_bounds.py", "1249/quantum_bounds_results.json",
        "1251/proof.md", "1252/proof.md", "1252/results.json",
        "1266/proof.md", "1267/geometry_proof.md",
        "1267/geometry_check.py", "1267/geometry_results.json",
        "_admission/causal_block_after1265/proof.md",
        "_admission/finite_semiclassical_budget_gate_after1256.md",
    ]
    bindings = {p: sha(ARCH/p) for p in paths}
    bindings.update({p: sha(HERE/p) for p in
                     ("admission.md", "readout_proof.md", "readout_check.py")})
    return {
        "round": 1268, "part": "finite_Ramsey_readout_on_same_quantum_prefix",
        "status": "PASS_AUTHOR_WITH_SCOPE", "author_groups": 1,
        "arithmetic": "exact integers and Fraction; no ODE or wavefunction simulation",
        "task": {
            "known_Gaussian_parameter_family": ["3", "3001/1000"],
            "family_program": "B=2",
            "extra_program_comparison": "kappa=3; B=1,2",
            "gate_window": [str(2-tau), "2"], "tau": str(tau),
            "alpha": str(alpha), "pulse_peak": str(pulse_peak),
            "coupling": "epsilon*alpha*g*F(Q)*sigma_z",
            "F": "xi_Q(Q)*exp(Q), Q-only multiplication",
            "Q_core": [-1, 2], "Q_support": [-2, 3],
            "same_resource": "s=1000*ceil(max(1,hbar))*2^(2*N); V0=10^-6",
            "epsilon_upper": "2^(-4*N); hbar fixed",
            "N_decimal_digits": len(str(N)),
            "readout_probability": "(1+sin(a/8))/2",
        },
        "analytic_checks": {
            "pulse_integral": str(integral),
            "F_sup_strict_upper": "27",
            "F_Lipschitz_bound": str(LF),
            "F_Lipschitz_power_bound": "2^136",
            "new_gate_intensity_peak_coefficient": str(gate_peak),
            "new_bounded_source_intensity_peak_coefficient": str(source_peak_coefficient),
            "single_Q_source_strict_bound": "epsilon*2^161",
            "clock_compensation_Weyl_F0_bound": "2^E0",
            "clock_Weyl_F0_explicit_integer_bits": clock_Weyl_F0.bit_length(),
            "new_source_strict_bound": "epsilon*2^(E0+26) < 2^-100",
            "Gaussian_and_Duhamel_margins_over_110_bits": True,
            "old_bounded_source_total_error_strict_upper": "2^-90",
            "classical_branch_core_and_source_margins_over_110_bits": True,
            "Ramsey_full_norm_error_strict_upper": "2^-100",
            "free_and_Ramsey_same_H": True,
            "a_short_window": ["2", "3"],
            "a_velocity_strict_upper": "4",
            "r_derivative_strict_lower": str(slope_lower),
            "cosine_rational_lower": str(cosine_lower),
        },
        "finite_prediction": {
            "new_actual_matching_input_per_task": str(beta),
            "statistics_halfwidth": str(stat),
            "independent_repeats_per_task": R,
            "at_most_four_tasks_failure_upper": str(four_failure),
            "time_averaging_probability_error": str(time_prob_error),
            "probability_to_a2_total_error": str(ep),
            "B_probability_gap_lower": str(B_gap),
            "B_probability_gap_coarse_strict_lower": "17/10000000",
            "kappa_endpoint_probability_gap_lower": str(K_gap),
            "kappa_endpoint_probability_gap_coarse_strict_lower": "58/10000000",
            "fixed_inverse_a_error": str(eta_a),
            "parameter_acceptance_bound": str(kappa_accept),
            "classical_common_pre_gate_field_residual_bound": str(field_prefix),
            "field_time_domain": ["0", str(2-tau)],
        },
        "scope": [
            "Known material Gaussian family; arbitrary independent pointer-reference extension.",
            "Specific probability formula requires initial pointer |+x>.",
            "Full mathematical instrument retains both outcomes and true conditional poststates.",
            "Actual Ramsey preparation, timing and readout mismatch is a new explicit input.",
            "No automatic transfer of old a-POVM systematic/source/causal permissions.",
            "Volumetric rod modulation is a candidate source dictionary, not a certified full material model.",
            "No spacetime local-field/collective-block identity from naming F(Q).",
            "Only declared bounded completed sources/tails; no uncompleted stress norm claim.",
            "New read history does not inherit the old post-t2 or 9/4 field bound.",
            "No CT adoption, full GR, experiment, or stage completion claim.",
        ],
        "display_only": {
            "B_probability_gap": float(B_gap),
            "kappa_probability_gap": float(K_gap),
            "a_inverse_error": float(eta_a),
            "field_prefix_bound": float(field_prefix),
            "four_task_failure": float(four_failure),
        },
        "bindings_sha256": bindings,
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("FIRST_EXCLUSIVE_WRITE PASS_AUTHOR_WITH_SCOPE; author_groups=1")
    else:
        saved = load(OUT)
        if saved != result:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS_AUTHOR_WITH_SCOPE; author_groups=1")
    print(json.dumps(result["display_only"], ensure_ascii=False))

if __name__ == "__main__":
    main()
