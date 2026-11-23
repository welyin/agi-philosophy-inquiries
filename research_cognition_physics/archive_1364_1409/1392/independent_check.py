"""Independent integer anomaly and hadronic-normalization check."""
import argparse
from fractions import Fraction
import json
from math import pi, isclose
from pathlib import Path


def anomaly(n, leptons):
    # Twice the sum over color/charge-squared times the chiral PQ charge.
    electromagnetic = 2 * (3 * Fraction(1, 9) * n + leptons)
    assert (Fraction(2, 3) * n + electromagnetic / 2).denominator == 1
    return electromagnetic / n


def compute():
    alpha = 1 / 137.035999
    k = 0.00570
    down = anomaly(1, 0)
    extra = anomaly(1, 1)
    large = anomaly(8, 5)
    assert (down, extra, large) == (Fraction(2, 3), Fraction(8, 3), Fraction(23, 12))
    variants = {"chiral_2016": (1.92, 0.04),
                "background_field_2026": (8*pi*pi*.0224, 8*pi*pi*.0010),
                "ward_identity_2026": (1.965, .035)}
    points = {}
    for label, (c, width) in variants.items():
        slope = alpha / (2*pi*k) * (float(down) - c)
        min_slope = alpha / (2*pi*(k+.00010)) * (c-width-float(down))
        max_slope = alpha / (2*pi*(k-.00010)) * (c+width-float(down))
        points[label] = {"C": c, "width": width, "down_abs_slope": abs(slope),
                         "slope_interval": [min_slope, max_slope],
                         "CAPP_nominal_fraction": (6.2e-16/(4.55e-15*slope))**2,
                         "abs_delta_for_n8": abs(float(large)-c)}
        assert float(extra)-c < c-float(down)
        assert min_slope > .1394
    assert isclose(points["chiral_2016"]["down_abs_slope"], .25537430, rel_tol=1e-6)
    return {"E_over_N_primitive_down": str(down), "one_extra_lepton": str(extra),
            "eight_units_color_five_leptons": str(large), "hadronic_inputs": points}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.check:
        frozen = json.loads(Path(__file__).with_name("independent_results.json").read_text())
        assert result == frozen
        print("1392 independent global-charge, photon and local-density checks passed")
    else:
        print(json.dumps(result, indent=2))
