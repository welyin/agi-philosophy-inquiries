"""Independent profile of the CMS 2022 decay-strength pair (muon, tau)."""
import argparse
import json
from math import isclose, log1p, sqrt
from pathlib import Path

import numpy as np


OBS = np.array([1.21, 0.85])
MU_MINUS, MU_PLUS, TAU_SIGMA, RHO = 0.42, 0.45, 0.10, 0.03


def projected(q, sigma_mu, rho):
    """Whiten the measured vector and project onto the physical nonnegative ray."""
    cross = rho*sigma_mu*TAU_SIGMA
    covariance = np.array([[sigma_mu**2, cross], [cross, TAU_SIGMA**2]])
    factor = np.linalg.cholesky(covariance)
    y = np.linalg.solve(factor, OBS)
    direction = np.linalg.solve(factor, [q, 1])
    gain = max(0., float(direction @ y / (direction @ direction)))
    residual = y-gain*direction
    discrepancy = float(residual @ residual)
    # Independent scalar contrast. The unconstrained solution lies on the ray here.
    assert gain > 0
    variance = np.array([1, -q]) @ covariance @ np.array([1, -q])
    assert isclose(discrepancy, (OBS[0]-q*OBS[1])**2/variance, abs_tol=1e-13)
    return gain, discrepancy


def inverted_interval(target):
    """Find the two positive roots by bisection on the projected discrepancy."""
    center = OBS[0]/OBS[1]
    roots = []
    for edge in [0., 8.]:
        inner, outer = center, edge
        for _ in range(70):
            probe = 0.5*(inner+outer)
            if projected(probe, (MU_MINUS+MU_PLUS)/2, RHO)[1] <= target:
                inner = probe
            else:
                outer = probe
        roots.append(float(0.5*(inner+outer)))
    assert roots[0] < center < roots[1]
    for root in roots:
        assert isclose(projected(root, (MU_MINUS+MU_PLUS)/2, RHO)[1], target,
                       abs_tol=1e-12)
    return roots


def tail_variants(value):
    """Two monotone extensions with the same curvature and every contour up to four."""
    if value <= 4:
        return [value, value]
    excess = value - 4
    return [4+log1p(excess), value+excess**2]


def compute():
    # The Born width of a scalar to two identical-mass leptons is proportional
    # to m^2 times the cube of the magnitude of either decay 3-momentum.
    lepton_mev = np.array([105.6583755, 1776.91])
    momenta_mev = np.sqrt((125380/2)**2-lepton_mev**2)
    born_widths = lepton_mev**2 * momenta_mev**3
    sm_width_ratio = 0.0216/6.21   # CMS Table B.1 gives both in percent
    rows = []
    for q in (0.25, 1., 4.):
        center_gain, center_d = projected(q, 0.435, RHO)
        # Corner choices bound the displayed grid values; not a confidence band.
        sensitivity = [projected(q, width, coeff)[1]
                       for width in (MU_MINUS, MU_PLUS)
                       for coeff in (0.025, 0.035)]
        rows.append({"q": q, "common_gain": center_gain, "profile": center_d,
                     "grid_profile_extrema": [min(sensitivity), max(sensitivity)],
                     "tail_relabelings": tail_variants(center_d)})
    low, high = inverted_interval(3.841458820694124)
    rectangle = [(OBS[0]-MU_MINUS)/(OBS[1]+TAU_SIGMA),
                 (OBS[0]+MU_PLUS)/(OBS[1]-TAU_SIGMA)]
    assert rectangle[0] < 1 < rectangle[1]
    assert rows[1]["profile"] < 1 < 4 < min(rows[i]["grid_profile_extrema"][0]
                                             for i in (0, 2))
    return {
        "source": "https://arxiv.org/pdf/2207.00043",
        "figures": ["2b", "B.5b"],
        "observables_mu_tau": OBS.tolist(),
        "uncertainties_mu_minus_plus_tau": [MU_MINUS, MU_PLUS, TAU_SIGMA],
        "correlation_rounded": RHO,
        "independent_mass_ratio": float(lepton_mev[0]/lepton_mev[1]),
        "tree_width_ratio": float(born_widths[0]/born_widths[1]),
        "cms_table_sm_width_ratio": sm_width_ratio,
        "width_ratio_from_experimental_central_values": float(sm_width_ratio*OBS[0]/OBS[1]),
        "best_width_ratio": float(OBS[0]/OBS[1]),
        "best_coupling_modulus_ratio": sqrt(OBS[0]/OBS[1]),
        "symmetric_local_gaussian": rows,
        "fieller_1": inverted_interval(1.),
        "fieller_3.84146": [low, high],
        "sqrt_fieller_3.84146": [sqrt(low), sqrt(high)],
        "marginal_rectangle_ratio": rectangle,
        "bonferroni_bound_if_each_marginal_has_nominal_coverage": 2*0.682689492-1
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true", help="compare with frozen results (default)")
    args = parser.parse_args()
    file = Path(__file__).with_name("independent_results.json")
    # Rounding the output absorbs floating-point variation across LAPACK builds.
    def rounded(value):
        if isinstance(value, dict):
            return {k: rounded(v) for k, v in value.items()}
        if isinstance(value, list):
            return [rounded(v) for v in value]
        return round(value, 10) if isinstance(value, float) else value
    actual = rounded(compute())
    if args.write:
        with file.open("x", encoding="utf-8") as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert json.loads(file.read_text(encoding="utf-8")) == actual
    print("1366 independent profile verified.")


if __name__ == "__main__":
    main()
