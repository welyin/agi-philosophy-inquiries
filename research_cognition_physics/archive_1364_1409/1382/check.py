"""Review the anomaly sum, the integer spectrum, and the coherent h -> photons holdout."""
import argparse
import json
import math
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def traces(z):
    """One colorless module: all hypercharges are expressed as left-handed Weyl charges."""
    y = F(z, 2)
    left = [y, y, -y - F(1, 2), -y + F(1, 2)]
    q2 = (y + F(1, 2)) ** 2 + (y - F(1, 2)) ** 2
    return y / 2, sum(a**3 for a in left), sum(left), sum(a*a for a in left), q2


def certify_integer(delta):
    # A trial score of 46 bounds each |z| <= 6 and D <= 4 (per term >= 10).
    charges = [i for i in range(-6, 7) if i % 2 and abs(i) >= 3]
    candidates = []
    for count in range(1, 5, 2):
        for group in combinations_with_replacement(charges, count):
            weight = sum(1 + z*z for z in group)
            if sum(group) == -delta and weight <= 46:
                candidates.append((weight, group))
    score, group = min(candidates)
    assert score == 46 and len(candidates) == 1
    return list(group)


def triangle(tau):
    if tau == 0:
        return 4/3
    return 2 * (tau + (tau - 1) * math.asin(math.sqrt(tau))**2) / tau**2


def certificate():
    spec = read_json(HERE / "inputs.json")
    old = read_json(HERE / spec["baseline_result_path"])
    inputs = read_json(HERE / spec["cms_input_path"])
    base = complex(*old["SM_gamma_amplitude"])
    ratios = {}
    for n in (1, 2):
        gamma, zz = inputs["cms"]["gammagamma"], inputs["cms"]["ZZ"]
        lower = (gamma[0]-n*gamma[1]) / (zz[0]+n*zz[2])
        upper = (gamma[0]+n*gamma[2]) / (zz[0]-n*zz[1])
        shift_max = -base.real + math.sqrt(upper*abs(base)**2 - base.imag**2)
        ratios[str(n)] = {"range": [lower, upper], "largest_real_shift": shift_max}
    cases = {}
    for delta, ordinary in ((1, F(5, 3)), (-1, F(1))):
        z = certify_integer(delta)
        mixed, cubic, grav, hyper, photon = [sum(row) for row in zip(*(traces(i) for i in z))]
        assert mixed + F(delta, 4) == cubic - F(3*delta, 4) == grav == 0
        assert photon == 23 and hyper == F(89, 2)
        shift = F(4, 3) * (ordinary + photon)
        cases[str(delta)] = {"z": z, "minimum_exotic_charge_squares": int(photon),
                             "extra_ordinary_charge_squares": str(ordinary),
                             "minimum_real_shift": float(shift),
                             "ratio_floor": abs(base + float(shift))**2 / abs(base)**2}

    frac = [F(s) for s in spec["fractional_reflection_example"]["z"]]
    assert sum(frac) == 1 and len(frac) % 2 == 1
    total_q2 = 1 + sum(traces(z)[-1] for z in frac)
    target = -2*base.real / float(total_q2)
    lo, hi = 0., 1.
    assert triangle(lo) < target < triangle(hi)
    for _ in range(64):
        mid = (lo+hi)/2
        if triangle(mid) < target:
            lo = mid
        else:
            hi = mid
    tau = (lo+hi)/2
    reflection = abs(base+float(total_q2)*triangle(tau))**2/abs(base)**2
    assert math.isclose(reflection, 1, rel_tol=2e-14)
    return {"baseline_complex": [base.real, base.imag], "cms_rectangles": ratios,
            "exo_beta_Y_min": str(F(2, 3)*F(89, 2)), "exo_beta_2_at_optimum": 1,
            "exo_heavy_amplitude_infimum": 92/3,
            "exo_only_floor_if_ordinary_omitted_conservatively": abs(base+92/3)**2/abs(base)**2,
            "two_sides_with_heavy_ordinary": cases,
            "fractional_competitor": {"total_additional_charge_squares": str(total_q2),
                "common_mass_GeV": inputs["higgs_mass_GeV"]/(2*math.sqrt(tau)),
                "ratio": reflection}}


def agrees(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(agrees(v, b[k]) for k, v in a.items())
    if isinstance(a, list):
        return isinstance(b, list) and len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-12) if isinstance(a, float) else a == b


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("--write", action="store_true")
    fresh = certificate()
    path = HERE / "results.json"
    if args.parse_args().write:
        with path.open("x", encoding="utf-8") as output:
            json.dump(fresh, output, indent=2, ensure_ascii=False)
            output.write("\n")
    else:
        assert agrees(fresh, read_json(path))
        print("1382: anomaly spectra, two unit mismatches and reflected competitor verified")
