"""Independent directional derivatives of the gauge-invariant radial potential in round 1653."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly

ROOT = Path(__file__).resolve().parent


def experiment(v, f, s, heavy, light=125.1, shifts=(0., 0.)):
    c = sqrt(1-s*s)
    Q = np.array([[c, -s], [s, c]])
    mass = Q @ np.diag([light**2, heavy**2]) @ Q.T
    scales = np.diag([1/v, 1/f])
    potential = scales @ mass @ scales / 2
    lower = np.linalg.eigvalsh(potential)[0]
    alpha, beta = shifts
    stab = (alpha**2+beta**2)/(4*lower) + 1e-18

    def jet(direction):
        x, y = direction
        X = Poly([0, v*x, x*x/2])
        L = Poly([0, f*y, y*y/2])
        V = potential[0, 0]*X**2+2*potential[0, 1]*X*L+potential[1, 1]*L**2
        V += alpha*X**3+beta*X**2*L+stab*(X**2+L**2)**2
        return np.array([2*V.coef[2], 6*V.coef[3]])

    h, H = Q[:, 0], Q[:, 1]
    h2, hhh = jet(h)
    H2, HHH = jet(H)
    plus, minus = jet(h+H), jet(h-H)
    hH = (plus[0]-minus[0])/4
    Hhh = (plus[1]-minus[1]-2*HHH)/6
    np.testing.assert_allclose([h2, H2, hH], [light**2, heavy**2, 0], atol=1e-9)
    K = hhh*v/(3*light**2)
    A = Hhh*v/(-(heavy**2+2*light**2)*s)
    xi = s*v/f
    weak_shifts = np.array([6*alpha*v**3, 2*beta*v**2*f])
    transfer = np.array([[c**3, 3*c**2*s], [-s*c**2, c*(c*c-2*s*s)]])
    k0 = c**3+s**2*xi
    a0 = c*(c-xi)
    shift = transfer @ weak_shifts
    np.testing.assert_allclose(np.linalg.det(transfer), c**4, rtol=2e-12)
    np.testing.assert_allclose(hhh, k0*(3*light**2/v)+shift[0], rtol=3e-12)
    np.testing.assert_allclose(Hhh, -a0*(heavy**2+2*light**2)*s/v+shift[1], atol=3e-10)
    mz = 91.1876
    t = (mz/heavy)**2
    zphase = sqrt(1-4*t)*(1-4*t+12*t*t)
    hphase = sqrt(1-4*(light/heavy)**2)
    whh = Hhh**2*hphase/(32*np.pi*heavy)
    wzz = heavy**3*s*s*zphase/(32*np.pi*v**2)
    arate = sqrt(whh/wzz*zphase/hphase)/(1+2*(light/heavy)**2)
    np.testing.assert_allclose(arate, abs(A), atol=1e-12)
    return {"f": f, "s": s, "xi": xi, "c": c, "k_lambda": K, "A": A,
            "hhh_GeV": hhh, "Hhh_GeV": Hhh, "partial_width_ratio_hh_ZZ": whh/wzz,
            "normalized_quartic_eigenvalues": np.linalg.eigvalsh(potential).tolist(),
            "positive_radial_stabilizer_GeVminus4":stab}


def closed(x):
    if isinstance(x, dict):
        return {k: closed(v) for k,v in x.items()}
    if isinstance(x, list):
        return [closed(z) for z in x]
    return float(f"{x:.11g}") if isinstance(x, (float,np.floating)) else x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    opts = ap.parse_args()
    cfg = json.loads((ROOT.parent/"1653"/"inputs.json").read_bytes())
    v, f, s = cfg["v_GeV"], cfg["f_GeV"], cfg["sin_theta"]
    ref = experiment(v, f, s, 600)
    # Two positive breaking scales, two positive radial quartic forms, one identical heavy width.
    first = experiment(v, 246., .6, 350)
    second = experiment(v, 147.6, .6, 350)
    np.testing.assert_allclose(first["A"], -second["A"], atol=1e-13)
    np.testing.assert_allclose(first["partial_width_ratio_hh_ZZ"], second["partial_width_ratio_hh_ZZ"])
    assert abs(first["k_lambda"]-second["k_lambda"])>.14

    c = sqrt(1-s*s)
    T = np.array([[c**3,3*c*c*s],[-s*c*c,c*(c*c-2*s*s)]])
    # Move either cubic by +15 GeV while keeping both eigenmasses, mixing and the other cubic fixed.
    changed=[]
    for target in ([15.,0.],[0.,15.]):
        a,b = np.linalg.solve(T, target)
        alpha,beta = a/(6*v**3),b/(2*v*v*f)
        row = experiment(v,f,s,600,shifts=(alpha,beta))
        np.testing.assert_allclose([row["hhh_GeV"]-ref["hhh_GeV"],
                                    row["Hhh_GeV"]-ref["Hhh_GeV"]],target,atol=1e-9)
        changed.append({"alpha_beta_GeVminus2":[alpha,beta],"result":row})
    out=closed({"1653_parameters_at_nominal_600_GeV":ref, "equal_heavy_width_points":[first,second],
               "counterterms_independent":changed})
    path=ROOT/"independent_results.json"
    if opts.check:
        assert json.loads(path.read_bytes())==out
    else:
        path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print("1654 independent PASS: direct cubic jets, equal heavy widths, two stable deformations")

if __name__=="__main__":
    main()
