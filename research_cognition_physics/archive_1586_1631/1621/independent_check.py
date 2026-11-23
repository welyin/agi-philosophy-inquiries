"""Independent tree portal, on-shell diphoton, and fixed-pole matching checks."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
MW, MZ, MT, MH, V = 80.3602, 91.1876, 172.5, 125.38, 246.22
CMS = (1.02, .11, .12)  # CMS 2602.18611 Table 5, diphoton/ZZ normalized


def triangle(x):
    return math.asin(math.sqrt(x)) ** 2


def scalar(mass):
    x = (MH / (2 * mass)) ** 2
    return (triangle(x) - x) / x ** 2


def z_scalar(mass):
    x, y = (MH / (2 * mass)) ** 2, (MZ / (2 * mass)) ** 2

    def g(z):
        return math.sqrt(1 / z - 1) * math.asin(math.sqrt(z))
    # 2 I_1(x,y), Arhrib et al. eqs. 3.13-15; same limit as A0.
    return (1 / (y - x) + (triangle(x) - triangle(y)) / (x - y) ** 2
            + 2 * y * (g(x) - g(y)) / (x - y) ** 2)


def top_w():
    w, t = (MH / (2 * MW)) ** 2, (MH / (2 * MT)) ** 2
    vector = -(2 * w ** 2 + 3 * w + 3 * (2 * w - 1) * triangle(w)) / w ** 2
    fermion = 2 * (t + (t - 1) * triangle(t)) / t ** 2
    return vector + 4 * fermion / 3


def increment(mass):
    return (MW / mass) ** 2 * scalar(mass)


def crossing(upper, norm, coupling=1, transfer=1):
    lower_mass, upper_mass = MH / 2, 1000.
    for _ in range(75):
        mid = (upper_mass + lower_mass) / 2
        rate = transfer * (1 + coupling * increment(mid) / norm) ** 2
        if rate > upper:
            lower_mass = mid
        else:
            upper_mass = mid
    return (lower_mass + upper_mass) / 2


def portal_example(single_mass):
    double_mass = math.hypot(single_mass, MW)
    a, b = scalar(single_mass) / single_mass ** 2, scalar(double_mass) / double_mass ** 2
    x0 = V ** 2 / 2
    c = MW ** 2 * a / (x0 * (a + 4 * b))
    # k_i = v g_hii/(2m_i²) ; both physical masses are held fixed.
    single = (-MW ** 2 + c * x0) / single_mass ** 2
    double = c * x0 / double_mass ** 2
    photon = single * scalar(single_mass) + 4 * double * scalar(double_mass)
    s2 = .231
    z = (single * (-.5 - s2) * z_scalar(single_mass)
         + 2 * double * (.5 - 2 * s2) * z_scalar(double_mass))
    common_mass_squared = double_mass ** 2 - c * x0
    assert common_mass_squared > 0 and single < 0 < double
    assert abs(photon) < 1e-13 and z > 0
    return dict(m_plus=single_mass, m_plus_plus=double_mass, portal=c,
                adjusted_hard_mass_squared=common_mass_squared,
                k_plus=single, k_plus_plus=double,
                photon_increment=photon,
                z_gamma_reduced_increment_at_s2_0231=z)


def calculate():
    inputs = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))
    couplings = inputs["central"]
    def inv(g):
        return 1 / couplings[g] ** 2
    d0 = inv("gY") - 3 * inv("g2")
    beta_gap = 41 / 6 + 3 * 19 / 6
    tree_cross = inputs["mu0"] * math.exp(8 * math.pi ** 2 * d0 / beta_gap)
    tree_gap = {}
    for mass in inputs["hard_doublet_examples_GeV"]:
        # Above a complete (2,3/2) threshold ΔbY−3Δb2=3/2−1/2=1.
        t = ((8 * math.pi ** 2 * d0 + math.log(mass / inputs["mu0"]))
             / (beta_gap + 1))
        tree_gap[str(mass)] = inputs["mu0"] * math.exp(t)
    assert 3550 < tree_cross < 3560 and all(m < tree_cross for m in tree_gap.values())
    nominal = abs(top_w())
    curves = {}
    for n in (1, 2):
        upper = CMS[0] + n * CMS[2]
        curves[str(n)] = {
            "CMS_upper": upper,
            "central_min_mplus": crossing(upper, nominal),
            "budget_min_mplus": crossing(upper, 7, .8, .95),
        }
    rates = {str(m): (1 + increment(m) / nominal) ** 2 for m in (70, 80, 100, 200, 300)}
    at_hard = {}
    for m0 in inputs["hard_doublet_examples_GeV"]:
        m = math.sqrt(m0 ** 2 - MW ** 2)
        at_hard[str(m0)] = dict(m_plus=m, branching_ratio=(1 + increment(m) / nominal) ** 2)
    return dict(masses_GeV=dict(h=MH, W=MW, Z=MZ, top=MT, v=V),
                CMS_ratio_centre_minus_plus=CMS,
                compact_scale_SM_GeV=tree_cross, compact_G2_split_1loop_GeV=tree_gap,
                reference_top_W_photon_amplitude=top_w(),
                sensitivity_GeV=curves, pole_rates= rates, at_hard_doublet=at_hard,
                fixed_mass_photon_blind=portal_example(70.))


def compare(stored, current):
    if isinstance(current, dict):
        return stored.keys() == current.keys() and all(compare(stored[k], v) for k, v in current.items())
    if isinstance(current, (tuple, list)):
        return len(stored) == len(current) and all(compare(a, b) for a, b in zip(stored, current))
    return math.isclose(stored, current, rel_tol=2e-12, abs_tol=1e-13)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    result = calculate()
    target = HERE / "independent_results.json"
    if args.check:
        assert compare(json.loads(target.read_text(encoding="utf-8-sig")), result)
    else:
        with target.open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, ensure_ascii=False)
            out.write("\n")
    print("Independent portal and diphoton checks passed.")
