"""Independent four-isotope rate identity, using exact tabulated IBM-2 amplitudes."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

OUT = Path(__file__).with_name("independent_results.json")
NUCLEI = ("Ge76", "Te130", "Xe136", "Mo100")
LIGHT = list(map(F, ("4.68", "3.70", "3.05", "4.22")))
HEAVY = list(map(F, ("104", "91.8", "72.6", "164")))


def interpolation_weights(x, target):
    weights = []
    for j in range(len(x)):
        if j == target:
            continue
        weight = F(1)
        for k in range(len(x)):
            if k != target and k != j:
                weight *= (x[target] - x[k]) / (x[j] - x[k])
        weights.append(weight)
    return weights


def calculate():
    x = [light / heavy for light, heavy in zip(LIGHT, HEAVY)]
    held_out = {}
    for target, name in enumerate(NUCLEI):
        weights = interpolation_weights(x, target)
        others = [v for j, v in enumerate(x) if j != target]
        for degree in (0, 1, 2):
            assert sum(w * v**degree for w, v in zip(weights, others)) == x[target]**degree
        # Explicit complex relative phase: eta_light=1, eta_heavy=-1/30+i/100.
        def rate(v):
            return (v - F(1, 30))**2 + F(1, 100)**2
        assert sum(w * rate(v) for w,v in zip(weights, others)) == rate(x[target])
        held_out[name] = {
            "input_order": [n for n in NUCLEI if n != name],
            "weights_exact": [str(w) for w in weights],
            "max_absolute_noise_amplification": float(sum(map(abs, weights))),
        }
    # A genuine third amplitude with a distinct response can violate the identity.
    weights = interpolation_weights(x, 3)
    third_at_mo = (x[3] - F(1, 30) + F(1, 100))**2 + F(1, 100)**2
    assert third_at_mo != sum(w * rate(v) for w,v in zip(weights,x[:3]))
    ge, te, xe, mo = range(4)
    assert x[mo] < x[te] < x[xe] < x[ge]
    return {
        "reduced_rate_definition": "1/(T_i G_i gA^4 H_i^2)",
        "response_ratios": dict(zip(NUCLEI, map(str, x))),
        "four_isotope_holdouts_at_fixed_NME": held_out,
        "third_response_at_Mo_toy_signed_contrast": str(third_at_mo-rate(x[mo])),
        "realistic_data_used_as_positive_rates": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    options = parser.parse_args()
    result = calculate()
    if options.write_new:
        with OUT.open("x", encoding="utf-8") as file:
            json.dump(result, file, indent=2)
            file.write("\n")
    else:
        assert result == json.loads(OUT.read_text(encoding="utf-8"))
    print("Verified all four isotope contrasts and a third-response witness.")
