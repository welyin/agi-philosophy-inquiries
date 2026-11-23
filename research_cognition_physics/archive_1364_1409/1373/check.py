"""Round 1373. Reproduce mass matching and the gamma-gamma/ZZ comparisons."""
import argparse
import json
from math import acosh, asin, pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def loop_functions(mh, mass):
    z = (mh/(2*mass))**2
    if z <= 1:
        f = asin(sqrt(z))**2
    else:
        f = -(2*acosh(sqrt(z))-1j*pi)**2/4
    scalar = (f-z)/z**2
    fermion = 2*(z+(z-1)*f)/z**2
    vector = -(2*z*z+3*z+3*(2*z-1)*f)/z**2
    return scalar, fermion, vector


def common_inputs(config):
    v = config["v_MZ_GeV"]
    gauge = read(config["gauge_and_heavy_yukawa"])
    prior = read(config["down_lepton"])
    light = dict(zip(prior["order"], np.array(prior["units"])*prior["points"]["MZ"]["central"]))
    up = read(config["up_charm"])["MZ"]
    up_y = [up["up"][0], up["charm"][0], gauge["yt"][0]]
    down_y = [light["yd"], light["ys"], gauge["yb"][0]]
    lepton_y = [light["ye"], light["ymu"], gauge["ytau"][0]]
    fermions = [(v*y/sqrt(2), 3, 2/3) for y in up_y]
    fermions += [(v*y/sqrt(2), 3, -1/3) for y in down_y]
    fermions += [(v*y/sqrt(2), 1, -1) for y in lepton_y]
    mw = gauge["g2"][0]*v/2
    mz = v*sqrt(gauge["g2"][0]**2 + 0.6*gauge["g1_GUT"][0]**2)/2
    return mw, mz, fermions


def profile(q, data, cov):
    direction = np.array([q, 1.0])
    best = direction @ np.linalg.solve(cov, data)/(direction @ np.linalg.solve(cov, direction))
    best = max(0., best)
    diff = data - best*direction
    return float(best), float(diff @ np.linalg.solve(cov, diff))


def evaluate():
    config = read("inputs.json")
    mh, v = config["higgs_mass_GeV"], config["v_MZ_GeV"]
    mw, mz, fermions = common_inputs(config)
    sm_trace = 6*mw**4 + 3*mz**4 - sum(4*nc*m**4 for m,nc,_ in fermions)
    missing = 8*pi**2*v**2*mh**2 - sm_trace
    assert sm_trace < 0 < missing
    scalar_mass = (missing/2)**0.25
    assert scalar_mass > mh/2
    standard = loop_functions(mh, mw)[2]
    standard += sum(nc*q**2*loop_functions(mh,m)[1] for m,nc,q in fermions)
    increment = loop_functions(mh, scalar_mass)[0]
    assert 1/3 < increment < 0.35 and standard.real < 0
    cms = config["cms"]
    g, gs, _ = cms["gammagamma"]
    z, zlo, zhi = cms["ZZ"]
    zs = (zlo+zhi)/2
    cov = np.array([[gs**2, cms["correlation"]*gs*zs],
                    [cms["correlation"]*gs*zs, zs**2]])
    data = np.array([g, z])
    hypotheses = {}
    for charge in config["scalar_charges"]:
        delta = charge**2*increment
        q = abs(standard+delta)**2/abs(standard)**2
        t, chisq = profile(q, data, cov)
        np.testing.assert_allclose(chisq, (g-q*z)**2/(cov[0,0]+q*q*cov[1,1]-2*q*cov[0,1]))
        reflection = -2*standard.real-delta
        np.testing.assert_allclose(abs(standard+delta),abs(standard+reflection))
        hypotheses[str(charge)] = {"width_ratio_to_SM":float(q), "profile_normalization":t,
                                   "local_gaussian_delta_chi_squared":chisq}
    # Four complex singlets of equal mass. Two charge assignments with equal sum Q^2.
    mass4 = (missing/8)**0.25
    same = []
    for charges in [[1,1,1,1],[2,0,0,0]]:
        same.append({"width_ratio":float(abs(standard+sum(q*q for q in charges)*
                   loop_functions(mh,mass4)[0])**2/abs(standard)**2),
                     "sum_Q4_over_M4_GeV^-4":sum(q**4 for q in charges)/mass4**4})
    np.testing.assert_allclose(same[0]["width_ratio"],same[1]["width_ratio"])
    np.testing.assert_allclose(same[1]["sum_Q4_over_M4_GeV^-4"]/same[0]["sum_Q4_over_M4_GeV^-4"], 4)
    return {"input_W_Z_tree_GeV":[mw,mz], "known_supertrace_GeV4":sm_trace,
            "required_positive_fourth_moment_GeV4":missing,"single_complex_mass_GeV":scalar_mass,
            "single_portal_lambda":2*scalar_mass**2/v**2,
            "SM_gamma_amplitude":[float(standard.real),float(standard.imag)],
            "A0":float(increment),"published_fit_central_ratio":g/z,
            "one_marginal_width_rectangle_ratio":[(g-gs)/(z+zhi),(g+gs)/(z-zlo)],
            "hypotheses":hypotheses,"reflected_charge_squared_at_Q0":-2*standard.real/increment,
            "four_complex_common_mass_GeV":mass4,"four_species_competitors":same}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a,list) and a and isinstance(a[0],dict):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            compare(x,y)
    else:
        np.testing.assert_allclose(a,b,rtol=2e-10,atol=1e-14)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    path = HERE/"results.json"
    values = evaluate()
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(values,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        compare(values,json.loads(path.read_text(encoding="utf-8")))
    print("1373 verified: fourth moment, finite mass loops, charge alternatives, joint fit and degeneracy.")
