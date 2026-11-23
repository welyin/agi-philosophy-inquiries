"""Round 1583: integer charges, mass volumes and conditional axion outputs."""
import argparse
from itertools import combinations, permutations, product
import json
from math import cos, exp, log, pi, sin, sqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
INPUT = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))


def min_order(left, right):
    """Least permutation valuation among all full column minors."""
    return min(sum(abs(left[i] + right[j]) for j, i in enumerate(rows))
               for rows in permutations(range(len(left)), len(right)))


def orders(p, s):
    a = [p, s, -s]
    b = [p, p, s, s, -s, -s]
    c = [p] * 3 + [s] * 3
    return [min_order(a, [-s] * 3), min_order(b, b), min_order(c, a)]


def crossing(g2, g3):
    g = INPUT["gauge"]
    a2, a3 = 4 * pi / g2**2, 4 * pi / g3**2
    step = 2*pi*(a2-a3)/(g["b2"]-g["b3"])
    return g["mz_gev"]*exp(step), a2-g["b2"]*step/(2*pi)


def volume_factor(cap):
    return (14/3*log(3**1.5*cap**3)
            + 10/3*log(720*cap**6)
            + 2*log(6*sqrt(20)*cap**3))


def thresholds():
    g = INPUT["gauge"]
    g2, w2 = g["g2"]
    g3, w3 = g["g3"]
    corners = [crossing(g2+i*g["widths"]*w2, g3+j*g["widths"]*w3)
               for i, j in product((-1, 1), repeat=2)]
    mlo, mhi = min(x[0] for x in corners), max(x[0] for x in corners)
    alo, ahi = min(x[1] for x in corners), max(x[1] for x in corners)
    mass, alpha = crossing(g2, g3)
    q = INPUT["qcd"]
    chi_hi = (q["chi_fourth_mev"] + 2*sum(q["chi_errors_mev"]))/1000
    ph = abs(8/3-q["photon_offset"]) + 2*q["photon_uncertainty"]
    rows = []
    omega = 188/3
    for cap, parent in product(g["coefficient_caps"], g["inverse_parent_floors"]):
        e = exp(-(2*pi*(ahi-3*parent)+volume_factor(cap))/omega)
        fa = sqrt(2)*mlo*e/3
        rows.append(dict(cap=cap, alpha15_inverse_floor=parent, epsilon_floor=e,
                         fa_floor_GeV=fa, qcd_mass_ceiling_eV=1.e9*chi_hi**2/fa,
                         photon_ceiling_inverse_GeV=ph/(2*pi*q["alpha_inverse"]*fa)))
    e = INPUT["sparse_epsilon"]
    return dict(central_M23_GeV=mass, central_inverse_alpha=alpha,
                M23_box_GeV=[mlo, mhi], inverse_alpha_box=[alo, ahi], scenarios=rows,
                sparse_inverse_parent=(alpha-omega*log(1/e)/(2*pi))/3,
                sparse_lightest_mass_GeV=mass*e**6)


def sparse_volume():
    e = INPUT["sparse_epsilon"]
    m15 = np.diag([e**2, 1., e**6])
    m24 = np.zeros((6, 6))
    for i in range(2):
        m24[i, i] = e**2
        m24[i+2, i+4] = m24[i+4, i+2] = 1.
    m10 = np.zeros((6, 3))
    m10[0, 0], m10[1, 1], m10[3, 2] = e**2, e**4, 1.
    vols = [float(np.prod(np.linalg.svd(m, compute_uv=False))) for m in (m15, m24, m10)]
    assert np.allclose(vols, [e**8, e**4, e**6], rtol=1.e-12, atol=0)
    return vols


def torque():
    q = INPUT["qcd"]
    theta = q["edm_ecm"] / ((q["edm_slope_efm"]-2*sum(q["edm_slope_errors_efm"]))*1.e-13)
    delta = pi*q["delta_floor_over_pi"]
    limit = 3*(1+q["torque_relative_allowance"])*theta/sin((delta-theta)/3)
    max_b_quarter_mev = (q["chi_fourth_mev"]+2*sum(q["chi_errors_mev"]))*limit**.25
    for x in np.linspace(-pi, pi, 67):
        for phase in [-2., -.1, .1, 2., pi]:
            expected = cos(((x+phase+pi) % (2*pi)-pi)/3)
            actual = max(cos((x+phase+2*k*pi)/3) for k in range(3))
            assert abs(actual-expected) < 1.e-14
    return dict(small_angle_limit=theta, B_over_chi_max_for_misalignment=limit,
                B_quarter_max_MeV=max_b_quarter_mev, aligned_B9chi_mass_ratio=sqrt(2))


def calculate():
    assert orders(1, 3) == orders(-1, -3) == [8, 4, 6]
    assert min(abs(2+9*k)+2*abs(3+13*k) for k in range(-10, 11) if k) == 27
    # Unit-instanton indices of the two 105's, conjugate 120 and three 15bars.
    assert 13 + 13 + 17 + 3 == 46 and 13*1 + 13*3 - 17*3 == 1
    # One-copy 5 has tr(Qem²)=4/3 and tr(Tcolor²)=1/2.
    assert abs(sum(v*v for v in [-1/3]*3+[1,0])/(.5)-8/3) < 1.e-14
    return dict(round=1583, positive_representative=[1,3], mass_powers_15_24_10=orders(1,3),
                weighted_threshold_power=188/3, other_representative_lower_bound=126,
                sparse_products=sparse_volume(), K15=1, Kcolor=3, E=8,
                primitive_qcd_vacua=3, instanton_zero_modes=46,
                matching=thresholds(), competition=torque())


def stable(obj):
    if isinstance(obj, dict):
        return {k: stable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [stable(v) for v in obj]
    return float(f"{obj:.12g}") if isinstance(obj, float) else obj


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = stable(calculate())
    target = ROOT / "results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
    elif json.loads(target.read_text(encoding="utf-8")) != result:
        raise SystemExit("1583 result mismatch")
    print("1583 verified")


if __name__ == "__main__":
    main()
