"""Round 1789. Identical WW partial wave, scalar spectral moments, finite windows."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


GF = 1.1663787e-5  # GeV^-2, conventional charged-current calibration
V2 = 1.0 / (math.sqrt(2.0) * GF)
RESULT = Path(__file__).with_name("independent_results.json")


def moment(spectrum):
    return sum(weight * mass**2 for weight, mass in spectrum)


def angular_projection(s, neutral, charged):
    """Numerical scalar Goldstone amplitude, identical initial and final pairs."""
    cosines, weights = np.polynomial.legendre.leggauss(160)
    t = -0.5 * s * (1.0 - cosines)
    u = -s - t
    amplitude = t + u
    for k2, mass in neutral:
        amplitude = amplitude - k2 * (t**2 / (t - mass**2) + u**2 / (u - mass**2))
    for z, mass in charged:
        amplitude = amplitude - z * s**2 / (s - mass**2)
    return float(np.dot(weights, amplitude) / (64 * math.pi * V2))


def finite_projection(s, neutral, charged):
    positive = 0.0
    for k2, mass in neutral:
        ratio = s / mass**2
        positive += 2 * k2 * mass**2 * (1 - math.log1p(ratio) / ratio)
    for z, mass in charged:
        positive += z * mass**2 / (1 - mass**2 / s)
    return -positive / (32 * math.pi * V2)


def window(budget, fraction, edge, neutral_moment_lower, charged_floor=0.0):
    total_lower = 2 * fraction
    mass_budget = budget - 2 * neutral_moment_lower
    if mass_budget < charged_floor**2 * total_lower:
        raise ValueError("the registered moments and lower thresholds are incompatible")
    weight = max(0.0, (edge**2 * total_lower - mass_budget) /
                 (edge**2 - charged_floor**2))
    return {"f": fraction, "edge_GeV": edge, "weight_at_masses_below_edge": weight,
            "fraction_of_all_doubly_charged_weight": weight / total_lower,
            "edge_above_which_some_weight_is_required_GeV": math.sqrt(mass_budget / total_lower)}


def run():
    # Deliberately use multiple poles and non-unit individual vertices.
    neutral = [(0.94, 125.0), (0.35, 650.0), (0.17, 950.0)]
    charged = [(0.10, 300.0), (0.19, 620.0), (0.17, 880.0)]
    assert abs(sum(a for a, _ in neutral) - sum(a for a, _ in charged) - 1) < 2e-15
    mass_sum = 2 * moment(neutral) + moment(charged)
    asymptotic = -mass_sum / (32 * math.pi * V2)
    projections = []
    largest_mass = max(m for _, m in neutral + charged)
    for q in (1.5, 3.0, 10.0, 30.0):
        s = q * largest_mass**2
        by_integration = angular_projection(s, neutral, charged)
        analytic = finite_projection(s, neutral, charged)
        assert abs(by_integration - analytic) < 1e-9
        ratio = analytic / asymptotic
        assert 1 - math.log1p(q) / q <= ratio <= q / (q - 1)
        projections.append({"s_over_largest_mass_squared": q, "analytic_a0": analytic,
                            "quadrature_a0": round(by_integration, 10), "finite_to_asymptotic": ratio})

    output = {}
    for name, factor in (("identical_charge_only", 16), ("also_opposite_charge", 8)):
        budget = factor * math.pi * V2
        windows = [window(budget, f, edge, 125.0**2)
                   for f, edge in ((0.10, 4000.0), (0.20, 3000.0), (0.25, 3000.0))]
        # Saturation as a two-atom measure when only the two moments are prescribed.
        floor = 200.0
        sharp = window(budget, 0.25, 3000.0, 125.0**2, floor)
        below = sharp["weight_at_masses_below_edge"]
        moments = 2 * 125.0**2 + moment([(below, floor), (0.50 - below, 3000.0)])
        assert abs(moments / budget - 1) < 1e-12
        output[name] = {"mass_budget_GeV2": budget, "windows": windows,
                        "sharp_two_atom_moment_example_with_200GeV_floor": sharp,
                        "if_no_weight_below_3TeV_then_f_upper": (budget - 2 * 125.0**2) / (2 * 3000.0**2),
                        "if_no_neutral_calibration_f_quarter_3TeV": window(budget, 0.25, 3000.0, 0)}

    return {"GF_GeV_minus2": GF, "v_GeV": math.sqrt(V2),
            "neutral_weight_mass_GeV": neutral, "charged_weight_mass_GeV": charged,
            "positive_mass_moment_GeV2": mass_sum, "asymptotic_same_sign_a0": asymptotic,
            "finite_energy_check": projections, "two_channel_conventions": output}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != json.loads(data):
            raise SystemExit("independent_results.json disagrees with the calculation")
    else:
        RESULT.write_text(data, encoding="utf-8")
    print("1789 WW identical partial wave", hashlib.sha256(data.encode()).hexdigest()[:12])
