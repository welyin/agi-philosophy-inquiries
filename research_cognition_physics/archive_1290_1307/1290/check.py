"""Fixed-smearing heavy KG stress variance: exact guards plus diagnostics.

Default: read-only recomputation and comparison. --write: exclusive first result.
No historical checker is imported. Gaussian quadrature is NOT a certified finite
time-window calculation; the compact-window certificate is the Fraction block.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents
            if (p / "research_cognition_physics").is_dir()
            and (p / "scripts").is_dir())
HISTORY = [
    "research_cognition_physics/archive_629_652/research_note_630.md",
    "research_cognition_physics/archive_629_652/research_note_631.md",
    "research_cognition_physics/archive_894_922/research_note_899.md",
    "research_cognition_physics/archive_935_955/research_note_950.md",
    "research_cognition_physics/archive_956_989/research_note_984.md",
    "research_cognition_physics/archive_1234_1280/_admission/"
    "source_fluctuation_scope_after1272/README.md",
]


def frac(x: F | int) -> str:
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def pair_tensor(mass, p, ep, q, eq):
    """On-shell covectors, lower stress indices; signature (-,+,+,+)."""
    m, ep, eq = F(mass), F(ep), F(eq)
    p, q = tuple(map(F, p)), tuple(map(F, q))
    assert ep * ep - sum(x * x for x in p) == m * m
    assert eq * eq - sum(x * x for x in q) == m * m
    pc, qc = (-ep,) + p, (-eq,) + q
    eta = (-1, 1, 1, 1)
    dot = sum(F(eta[i]) * pc[i] * qc[i] for i in range(4))
    t = [[-(pc[i] * qc[j] + pc[j] * qc[i]) / 2
          + (F(eta[i]) * (dot - m * m) / 2 if i == j else 0)
          for j in range(4)] for i in range(4)]
    return t


def exact_checks():
    beta = sum(F((-1) ** j * math.comb(3, j), 4 + j)
               for j in range(4))
    assert beta == F(1, 140)
    prefactor_without_pi4 = F(9, 8) * beta
    assert prefactor_without_pi4 == F(9, 1120)

    cases = [
        (3, (4, 0, 0), 5, (0, 4, 0), 5),
        (1, (F(3, 4), 0, 0), F(5, 4), (0, 0, F(4, 3)), F(5, 3)),
        (2, (0, 0, F(3, 2)), F(5, 2), (0, 0, 0), 2),
        (1, (0, 0, 0), 1, (0, 0, 0), 1),
        (5, (12, 0, 0), 13, (0, 0, -12), 13),
    ]
    ward_rows = []
    max_ratio = F(0)
    for mass, p, ep, q, eq in cases:
        t = pair_tensor(mass, p, ep, q, eq)
        momentum = [F(ep) + F(eq)] + [F(p[i]) + F(q[i]) for i in range(3)]
        ward = [sum(momentum[i] * t[i][j] for i in range(4))
                for j in range(4)]
        assert ward == [0, 0, 0, 0]
        ratio = max(abs(x) / (F(ep) * F(eq)) for row in t for x in row)
        max_ratio = max(max_ratio, ratio)
        assert ratio <= 3
        opp = pair_tensor(mass, p, ep, tuple(-F(x) for x in p), ep)
        assert opp[0][0] == 0
        ward_rows.append({
            "mass": frac(mass), "p": [frac(x) for x in p],
            "q": [frac(x) for x in q],
            "ward": [frac(x) for x in ward], "opposite_t00": frac(opp[0][0]),
        })

    # The sixth derivative of six unit-mass boxes is a finite signed measure.
    width = F(1, 8)
    atoms = [(F(j - 3) * width,
              F((-1) ** j * math.comb(6, j)) / width ** 6)
             for j in range(7)]
    tv = sum(abs(c) for _, c in atoms)
    assert tv == 16 ** 6
    moment_checks = [sum(c * x ** k for x, c in atoms) for k in range(7)]
    assert moment_checks[:6] == [0] * 6
    assert moment_checks[6] == math.factorial(6)
    support_radius = 3 * width + F(1, 8)
    assert support_radius == F(1, 2)
    # pi > 3; strict < is analytic, stored rational is the looser upper bound.
    mass = 2 ** 20
    compact_bound = F(9, 1120 * 3 ** 4 * (2 * 6 - 8)) * tv ** 2 / (2 * mass) ** 4
    assert compact_bound == F(1, 40320 * 2 ** 36)
    assert compact_bound < F(1, 10 ** 15)
    exponents = [{"r": r, "mass_power": 8 - 2 * r,
                  "coefficient_without_pi4": frac(F(9, 1120 * (2 * r - 8)))}
                 for r in (5, 6, 8, 10)]

    # Rotational covariance of the symmetric rank-two tensor S.
    angular_rows = []
    for c in (F(-1), F(-1, 2), F(0), F(1, 2), F(1)):
        trace_s2_over_p2q2 = (1 + c * c) / 2
        trace_s_squared_over_p2q2 = c * c
        invariant_value = (trace_s2_over_p2q2
                           - trace_s_squared_over_p2q2 / 3) / 10
        explicit = F(1, 20) + c * c / 60
        assert invariant_value == explicit
        angular_rows.append({"cosine": frac(c), "factor": frac(explicit)})
    assert angular_rows[-1]["factor"] == "1/15"
    return {
        "beta_4_4": frac(beta),
        "radial_coefficient_without_pi4": frac(prefactor_without_pi4),
        "ward_exact_cases": ward_rows,
        "sample_max_vertex_to_energy_product": frac(max_ratio),
        "sample_bound_is_not_the_continuum_proof": True,
        "mass_power_family": exponents,
        "compact_example": {
            "time_support": ["-1/2", "1/2"],
            "uniform_kernel_full_width": frac(width),
            "convolution_count": 6,
            "mollifier_support_radius": "1/8",
            "sixth_derivative_L1_upper": frac(tv),
            "signed_measure_moments_0_through_6": [frac(x) for x in moment_checks],
            "mass": mass,
            "variance_upper_using_pi_gt_3": frac(compact_bound),
            "variance_target": "1/1000000000000000",
            "strict_margin": frac(F(1, 10 ** 15) - compact_bound),
        },
        "angular_tensor_invariant_checks": angular_rows,
    }


def gaussian_variance(mass: float, nodes: int, cutoff: float = 8.0):
    """Gauss-Legendre on a finite box: diagnostic, no quadrature certification."""
    z, weights = np.polynomial.legendre.leggauss(nodes)
    p = ((z + 1) * cutoff / 2)[:, None, None]
    q = ((z + 1) * cutoff / 2)[None, :, None]
    c = z[None, None, :]
    w = ((weights * cutoff / 2)[:, None, None]
         * (weights * cutoff / 2)[None, :, None]
         * weights[None, None, :])
    ep = np.sqrt(p * p + mass * mass)
    eq = np.sqrt(q * q + mass * mass)
    integrand = (p ** 4 * q ** 4 / (ep * eq)
                 * (1 / 20 + c ** 2 / 60)
                 * np.exp(-(ep + eq) ** 2 - (p * p + q * q + 2 * p * q * c)))
    return float(np.sum(w * integrand) / (16 * math.pi ** 4))


def diagnostics():
    rows = []
    for mass in (0.5, 1.0, 2.0, 4.0):
        coarse = gaussian_variance(mass, 32)
        fine = gaussian_variance(mass, 48)
        analytic_bound = (9 * math.factorial(6)
                          / (1120 * math.pi ** 4 * 4) * (2 * mass) ** (-4))
        assert coarse > 0 and fine > 0
        assert coarse < analytic_bound and fine < analytic_bound
        rows.append({
            "mass": mass,
            "nodes_each_axis": [32, 48],
            "radial_cutoff": 8.0,
            "gaussian_variance_coarse": coarse,
            "gaussian_variance_fine": fine,
            "relative_difference": abs(coarse - fine) / fine,
            "r6_analytic_bound_B6_squared_le_720": analytic_bound,
            "fine_to_bound_ratio": fine / analytic_bound,
        })
    assert all(rows[j + 1]["gaussian_variance_fine"]
               < rows[j]["gaussian_variance_fine"] for j in range(len(rows) - 1))
    return {
        "role": "Schwartz diagnostic only; not compact support or certified quadrature",
        "tau": 1.0, "spatial_L": 1.0, "tensor_component": "T12",
        "rows": rows,
    }


def calculate():
    sources = {"proof.md": digest(HERE / "proof.md"),
               "selection.md": digest(HERE / "selection.md"),
               "check.py": digest(Path(__file__).resolve())}
    sources.update({p: digest(ROOT / p) for p in HISTORY})
    return {
        "schema": "heavy-source-fixed-smearing-v1",
        "status": "PASS_CANDIDATE_SOURCE_BOUND_ONLY",
        "candidate_check_groups": 1,
        "formal_rounds_assigned_by_this_file": 0,
        "source_bindings": sources,
        "primary_source": {
            "url": "https://arxiv.org/abs/gr-qc/0001098",
            "read_scope": "II (7),(12)-(16); III.1 (21)-(29); (79)",
            "noise_normalization_not_imported": True,
        },
        "exact": exact_checks(),
        "diagnostic": diagnostics(),
        "scope": {
            "fixed_mass_independent_finite_spacetime_test_menu": True,
            "minimal_free_scalar_Minkowski_vacuum": True,
            "mass_family_normal_ordered_zero_mean_is_declared_renormalized_rule": True,
            "bare_counterterms_fixed_across_masses": False,
            "nonzero_portal_transport_proved": False,
            "curved_background_mean_source_or_geometry_proved": False,
            "joint_classical_distribution_for_noncommuting_tests_claimed": False,
            "actual_instrument_constructed": False,
            "Gaussian_is_a_strict_finite_window": False,
        },
    }


def compare(actual, saved, path="$"):
    if type(actual) is not type(saved):
        raise AssertionError(f"type mismatch at {path}")
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise AssertionError(f"key mismatch at {path}")
        for key in actual:
            compare(actual[key], saved[key], f"{path}.{key}")
    elif isinstance(actual, list):
        if len(actual) != len(saved):
            raise AssertionError(f"length mismatch at {path}")
        for j, (x, y) in enumerate(zip(actual, saved)):
            compare(x, y, f"{path}[{j}]")
    elif isinstance(actual, float):
        if not math.isclose(actual, saved, rel_tol=1e-10, abs_tol=1e-100):
            raise AssertionError(f"numeric mismatch at {path}: {actual} != {saved}")
    elif actual != saved:
        raise AssertionError(f"value mismatch at {path}: {actual} != {saved}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="Create results.json once; fail if it already exists.")
    args = parser.parse_args()
    result = calculate()
    destination = HERE / "results.json"
    if args.write:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True,
                      allow_nan=False)
            stream.write("\n")
        receipt = "CREATED_EXCLUSIVELY"
    else:
        compare(result, json.loads(destination.read_text(encoding="utf-8")))
        receipt = "READ_ONLY_MATCH"
    print(json.dumps({"status": result["status"], "receipt": receipt,
                      "compact_variance_upper":
                          result["exact"]["compact_example"]["variance_upper_using_pi_gt_3"],
                      "result_sha256": digest(destination)}, ensure_ascii=False))


if __name__ == "__main__":
    main()

