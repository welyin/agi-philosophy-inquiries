"""Independent review: Weyl-gamma traces, Beta series, rational bounds.

Does not import the author's scientific functions. Default is read-only;
--write exclusively creates the separate reviewer result file.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
OUT = HERE / "review_independent_results.json"
I2, O2, I4 = np.eye(2), np.zeros((2, 2)), np.eye(4)
sigmas = [np.array([[0, 1], [1, 0]]),
          np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
gamma = [np.block([[O2, I2], [I2, O2]])]
gamma += [np.block([[O2, x], [-x, O2]]) for x in sigmas]
g5, eta = np.diag([-1, -1, 1, 1]), np.diag([1, -1, -1, -1])
zvertex = [g @ (-.04 * I4 + .5 * g5) for g in gamma]
me, mm, alpha = 3e-6, 6e-4, 1 / 137
e2 = 4 * math.pi * alpha
gz2 = e2 / (4 * .23 * .77)
mz2, mh2 = 2 * gz2, .72 ** 2
kappa = 8 * alpha / (3 * math.pi)


def slash(p):
    return sum(p[i] * eta[i, i] * gamma[i] for i in range(4))


def weight(s, c):
    root = math.sqrt(s)
    ee = (s + me * me - mm * mm) / (2 * root)
    em = (s + mm * mm - me * me) / (2 * root)
    k = math.sqrt(((s - me * me - mm * mm) ** 2 - 4 * me * me * mm * mm) / (4 * s))
    p1, p2 = np.array([ee, 0, 0, k]), np.array([em, 0, 0, -k])
    v3 = np.array([k * math.sqrt(1 - c * c), 0, k * c])
    p3, p4 = np.r_[ee, v3], np.r_[em, -v3]
    a, b = slash(p3) + me * I4, slash(p1) + me * I4
    cc, d = slash(p4) + mm * I4, slash(p2) + mm * I4
    q = p3 - p1
    t = q @ eta @ q
    tz = eta - np.outer(eta @ q, eta @ q) / mz2
    ge = np.array([[np.trace(a @ x @ b @ y) for y in gamma] for x in gamma])
    gm = np.array([[np.trace(cc @ x @ d @ y) for y in gamma] for x in gamma])
    ze = np.array([[np.trace(a @ x @ b @ y) for y in gamma] for x in zvertex])
    zm = np.array([[np.trace(cc @ x @ d @ y) for y in gamma] for x in zvertex])
    he = np.array([np.trace(a @ b @ y) for y in gamma])
    hm = np.array([np.trace(cc @ d @ y) for y in gamma])
    gg = np.einsum("ab,cd,ac,bd", eta, eta, ge, gm).real / 4
    zg = np.einsum("ab,cd,ac,bd", tz, eta, ze, zm).real / 4
    hg = np.einsum("ab,a,b", eta, he, hm).real / 4
    inter = (e2 / t) ** 2 * gg
    inter += (e2 / t) * (gz2 / (t - mz2)) * zg
    inter += (e2 / t) * (-me * mm / (2 * (t - mh2))) * hg
    return 2 * inter / (32 * math.pi * s), t


def integral():
    xs, ws = leggauss(20)
    xc, wc = leggauss(28)
    total = np.zeros(3)
    smallest = math.inf
    for x, w in zip(xs, ws):
        s = .0065 + .0025 * x
        for y, v in zip(xc, wc):
            c = -.1 + .5 * y
            ww, t = weight(s, c)
            a = -t
            pp, ss = 0., kappa / 3
            for j in range(1, 21):
                beta = math.factorial(j + 1) ** 2 / math.factorial(2 * j + 3)
                pp += kappa * (-1) ** (j + 1) * a ** j * beta / j
                ss += 2 * kappa * (-1) ** j * a ** j * beta
            total += (w / 2) * (v / 2) * ww * np.array([1, pp, ss])
            smallest = min(smallest, ww)
    return {"W": float(total[0]), "top_bin": float(total[1]),
            "source_bin": float(total[2]),
            "source_defect": float(total[0] * kappa / 3),
            "minimum_sampled_W": smallest,
            "algorithm": "Weyl gamma spin traces; 20 Beta moments; 20x28 bin quadrature"}


def rational_certificate():
    # Classical rational enclosure of pi; all subsequent decisions are exact.
    plo, phi = F(3141592653589793, 10**15), F(3141592653589794, 10**15)
    al, m1, m2 = F(1, 137), F(3, 10**6), F(6, 10**4)
    slo, shi, clo, chi = F(4, 1000), F(9, 1000), F(-6, 10), F(4, 10)
    sigma = m1*m1 + m2*m2
    k2lo = ((slo-sigma)**2 - 4*m1*m1*m2*m2)/(4*slo)
    qlo, qhi = 2*k2lo*(1-chi), shi*(1-clo)/2
    e2lo, e2hi = 4*plo*al, 4*phi*al
    mzlo = e2lo/(2*F(23,100)*F(77,100))
    glo = 2*e2lo**2*((slo-sigma)**2-2*qhi*sigma)/qhi**2
    ghi = 4*e2hi**2*shi**2/qlo**2
    scale = 10**20
    rootlo = F(math.isqrt(glo.numerator*scale**2//glo.denominator), scale)
    assert rootlo**2 <= glo
    zhi = F(1,2)*(16+12*qhi/mzlo)*(shi/4)*F(54,100)**2
    hhi = m1*m2*shi/(2*F(72,100)**2)
    etahi = 2*(zhi+hhi)/rootlo
    alo = (1-etahi)*glo/(16*phi*shi)*(chi-clo)
    ahi = (1+etahi)*ghi/(16*plo*slo)*(chi-clo)
    kaplo, kaphi = 8*al/(3*phi), 8*al/(3*plo)
    astar = qhi/F(9,10)**2
    b4 = F(math.factorial(4)**2, math.factorial(9))
    source_loss_lo = kaplo*alo/3
    source_error_hi = ahi*2*kaphi*astar**3*b4
    assert etahi < F(149,1000)
    assert alo > F(9767,10**6) and ahi < 11
    assert source_loss_lo > F(2,100000)
    assert source_error_hi < F(15,10**11)
    pairs = {"eta_upper": etahi, "A_lower": alo, "A_upper": ahi,
             "source_loss_lower": source_loss_lo, "source_error_upper": source_error_hi}
    return {"exact_inequality_checks_passed": True,
            "pi_enclosure": [str(plo), str(phi)],
            "decimal_diagnostics": {k: float(v) for k,v in pairs.items()},
            "exact_bounds": {k: str(v) for k,v in pairs.items()}}


def run():
    result = {"reviewer_calculation_only": True, "new_scientific_groups": 0,
              "trace_calculation": integral(), "rational_certificate": rational_certificate()}
    author = json.loads((HERE / "top_threshold_results.json").read_text(encoding="utf-8"))
    bin0 = author["finite_flux_bin"]
    row = bin0["kernel_bins_on_fixed_reference_external_basis"][1]
    refs = {"W": bin0["W_integral"], "top_bin": row["top_bin_coefficient"],
            "source_bin": row["radial_top_kernel_bin"]}
    differences = {k: result["trace_calculation"][k]-v for k,v in refs.items()}
    assert abs(differences["W"]) < 2e-13
    assert abs(differences["top_bin"]) < 2e-18
    assert abs(differences["source_bin"]) < 2e-16
    result["differences_from_author"] = differences
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        assert result == saved
    print(json.dumps({"all_passed": True, "trace": result["trace_calculation"],
                      "exact_rational_bounds": result["rational_certificate"]["decimal_diagnostics"]}))


if __name__ == "__main__":
    main()
