"""Check finite-mass scalar triangles, vacuum stability and a two-observable ambiguity."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def old(name):
    return json.loads((BASE.parent/"1373"/name).read_text(encoding="utf-8-sig"))


def form_factor(mass, higgs, neutral, order=45):
    """Positive Feynman-parameter representation, normalized to 1/3 at both zero momenta."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (1+nodes)/2, weights/2
    value = 0.
    for x, weight in zip(nodes, weights):
        y = (1-x)*nodes
        denom = 1-y*(higgs**2*x+neutral**2*(1-x-y))/mass**2
        value += weight*(1-x)*sum(weights*8*x*y/denom)
    return value


def analytic_ratio(m, h, z):
    a, b = 4*m*m/h**2, 4*m*m/z**2
    f = lambda k: math.asin(1/math.sqrt(k))**2
    g = lambda k: math.sqrt(k-1)*math.asin(1/math.sqrt(k))
    i1 = a*b/(2*(a-b)) + (a*b/(a-b))**2*(f(a)-f(b))/2 + a*a*b*(g(a)-g(b))/(a-b)**2
    return 2*i1/(a*(a*f(a)-1))


def traces(fields):
    q, t = [], []
    for spin, hypercharge in fields:
        m = np.arange(-spin, spin+1, 1.)
        q.extend(m+hypercharge)
        t.extend(m)
    q, t = np.array(q), np.array(t)
    return {"charges": sorted(q.tolist()), "Q2":sum(q*q), "QT":sum(q*t),
            "T2":sum(t*t), "Q4":sum(q**4), "Q3T":sum(q**3*t)}


def audit():
    data, known = old("inputs.json"), old("results.json")
    h, v = data["higgs_mass_GeV"], data["v_MZ_GeV"]
    w, z = known["input_W_Z_tree_GeV"]
    s2 = 1-w*w/(z*z)
    c2 = 1-s2
    m = 200.
    a0, z0 = form_factor(m,h,0), form_factor(m,h,z)
    ratio = z0/a0
    assert math.isclose(ratio, analytic_ratio(m,h,z), rel_tol=1e-12)
    for mass in [70., 100., 200., 550.]:
        numerical = form_factor(mass,h,z)/form_factor(mass,h,0)
        assert 1 < numerical < 1/(1-z*z/(4*mass**2))
        assert math.isclose(numerical, analytic_ratio(mass,h,z), rel_tol=1e-10)
    amplitude = complex(*known["SM_gamma_amplitude"])
    published = data["cms"]
    lambda_h = h*h/(2*v*v)
    mirror = {}
    for n in (1,2):
        gam, zz = published["gammagamma"], published["ZZ"]
        interval = [(gam[0]-n*gam[1])/(zz[0]+n*zz[2]),
                    (gam[0]+n*gam[2])/(zz[0]-n*zz[1])]
        min_kappa = (-amplitude.real+math.sqrt(interval[0]*abs(amplitude)**2-amplitude.imag**2))/(5*a0)
        lambda_min = (m/v)**4*(min_kappa-1)**2/lambda_h
        mirror[str(n)] = {"lower_gamma_ratio": interval[0],"mirror_kappa_min":min_kappa,
                          "necessary_charged_ray_quartic_min":lambda_min}
    max_kappa = 1+(v/m)**2*math.sqrt(lambda_h*4*np.pi)
    choices = [traces([(0,1),(0,2)]),traces([(.5,1.5)]),traces([(1,1)])]
    slopes = [(t["QT"]/t["Q2"]-s2)/math.sqrt(s2*c2)*ratio for t in choices]
    a, b = traces([(.5,.5),(0,2)]), traces([(.5,1.5),(0,0)])
    for moment in ["charges","Q2","QT","T2","Q4"]:
        assert a[moment] == b[moment]
    assert (a["Q3T"],b["Q3T"]) == (.5,3.5)
    energy = 500.
    f = energy**2/(energy**2-z*z)
    red = lambda third: -2+f*(third-2*s2)/c2
    production = (red(0)/red(.5))**2
    return {"mass_GeV": m, "weak_sine_squared_scheme1373":s2,
            "finite_scalar_covariant_hgg_and_hZgamma": [a0,z0],
            "Zgamma_to_diphoton_kinematic_transfer":ratio,
            "provable_upper_all_higgs_masses_below_threshold_at_200": 1/(1-z*z/(4*m*m)),
            "physical_tensor_NP_Zgamma_to_gammagamma_slopes": slopes,
            "mirror_scalar_vacuum_necessary":mirror,
            "kappa_ceiling_global_EW_and_lambda_dir_at_most_4pi":max_kappa,
            "identical_light_sums": {"doublet_half_plus_singlet_two": a,"doublet_three_halves_plus_neutral": b},
            "right_electron_Q2_pair_rate_ratio_at_500GeV_tree":production,
            "asymptotic_ratio":16/9}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
