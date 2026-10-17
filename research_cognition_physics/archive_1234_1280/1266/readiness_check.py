"""1266 readiness: original Gaussian -> true spectral blank/K -> old flag gate.
Default reads and compares; --write or --write-once creates the first saved result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "research_cognition_physics/archive_1234_/1249/quantum_bounds_results.json").is_file())
BASE = ROOT / "research_cognition_physics/archive_1234_"
assert HERE == (BASE / "1266").resolve(), "Unexpected formal-round location"
OUT = HERE / "readiness_results.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def tr(a):
    return [list(x) for x in zip(*a)]


def calculate():
    qpath = BASE / "1249/quantum_bounds_results.json"
    cpath = BASE / "_admission/causal_block_after1265/results.json"
    old = json.loads(qpath.read_text(encoding="utf-8"))
    cap = json.loads(cpath.read_text(encoding="utf-8"))
    group = old["group_constructive_finite_resource_bound"]
    ledger = group["ledger"]
    N = int(ledger["integer_n_power_exponent_decimal"])
    EF = int(ledger["Fourier_transport_exponent_decimal"])
    EH = int(ledger["Fourier_H_exponent"])
    assert N == int(cap["finite_resource"]["N_decimal"])
    assert 4*N-EF > 42 and 8*N-EH-EF > 42
    assert 4*N > 101
    B = int(ledger["B_decimal"])
    Lcut = int(ledger["Lcut_decimal"])
    Lnew = Lcut + 128*B
    LK = Lcut + 4_000_000_000*B
    assert Lnew < LK
    eta = F(1, 10**9)
    centers = list(map(F, group["terminal_K"]["centers"]))
    assert centers[1] - centers[0] > 6*eta
    for interval, center in zip(group["terminal_K"]["old_K_intervals"], centers):
        assert max(abs(F(v)-center) for v in interval) < F(3, 10**15) < eta/2
    assert F(1,32)-F(1,64) == F(1,64)
    dB, dK = F(1,32), eta
    coeffB, coeffK = 1/(2*dB**2), 1/(2*dK**2)
    assert coeffB == 512
    assert coeffK == 500_000_000_000_000_000 < 2**59
    # epsilon < 2^-101 is a strict conservative consequence of 4N>101.
    # Each heat-tail contribution is strictly < 2^-42.
    assert coeffB < 2**59
    individual = F(3, 2**42)
    joint = 2*individual
    assert individual < F(1,2**40)
    assert joint == F(3,2**41) < F(1,2**39)
    assert F(1,2**20)**2 == F(1,2**40)
    # Exact full 3-flag rotations on each K core: blank index is 2.
    blank = [[0,0,0],[0,0,0],[0,0,1]]
    effects = []
    for j in (0,1):
        U = [[int(i==k) for k in range(3)] for i in range(3)]
        U[j][j] = U[2][2] = 0
        U[j][2], U[2][j] = 1, -1
        identity = [[int(i==k) for k in range(3)] for i in range(3)]
        assert mul(tr(U),U) == identity
        Q = [[int(i==j and k==j) for k in range(3)] for i in range(3)]
        E = mul(mul(tr(U),Q),U)
        assert E == blank
        effects.append(E)
    paths = [
        Path(__file__).resolve(), HERE / "readiness.md",
        qpath, cpath, BASE/"1249/quantum_bridge.md",
        BASE/"1249/canonical.md", BASE/"1259/proof.md",
    ]
    return {
        "round": 1266,
        "part": "readiness",
        "status": "PASS_AUTHOR_WITH_SCOPE",
        "author_groups": 1,
        "old_resource": {
            "N_decimal_digits": len(str(N)), "epsilon_upper": "2^(-4N)",
            "linear_margin_gt_42": True, "quadratic_margin_gt_42": True,
            "four_N_gt_101": True,
        },
        "new_AW_blank_menu": {
            "A": 1, "L_exact": str(Lnew), "L_le_old_LK": True,
            "core_halfwidth": "1/64", "support_halfwidth": "1/32",
            "blank_spectral_halfwidth": "1/16",
        },
        "true_spectral_K_menu": {
            "centers": list(map(str, centers)), "eta": str(eta),
            "AW_support_halfwidth": str(eta), "spectral_core_halfwidth": str(2*eta),
            "flag_angle_support_halfwidth": str(3*eta), "disjoint": True,
        },
        "rational_bounds": {
            "AW_transport_error_strict_upper": str(F(1,2**41)),
            "blank_spectral_tail_coefficient_times_epsilon": str(coeffB),
            "K_spectral_tail_coefficient_times_epsilon": str(coeffK),
            "each_spectral_failure_strict_upper": str(individual),
            "actual_flag_error_strict_upper": str(joint),
            "blank_normalized_comparison_half_trace_distance_strict_upper": str(F(1,2**20)),
        },
        "exact_flag_core_effects": effects,
        "scope": [
            "same original Gaussian and actual chi*h propagation",
            "fixed-time uniform bounds, not continuous no-exit monitoring",
            "whole native effective block, not per-label K",
            "no actual projection or postselection of blank",
            "probability applies to all compatible reference extensions",
            "hard blank comparison does not certify H1 or psi-squared domain",
            "no gate-time source, causal-access, or autonomous-clock certificate",
        ],
        "bindings": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", "--write-once", dest="write_once", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_once:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        receipt = "FIRST_WRITE_EXCLUSIVE"
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        if result != saved:
            raise SystemExit("READ_ONLY_MISMATCH")
        receipt = "READ_ONLY_MATCH"
    print(json.dumps({"status":result["status"],"receipt":receipt,
                      "bound":result["rational_bounds"]["actual_flag_error_strict_upper"],
                      "results_sha256":sha(OUT)},ensure_ascii=False))

if __name__ == "__main__":
    main()
