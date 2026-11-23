"""Two physical weak currents of a common Pati-Salam neutral sector."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))


def from_asymmetry(asym):
    # Use the electron branch with |g_V/g_A| below one.
    return (1 - asym / (1 + math.sqrt(1 - asym**2))) / 4


def measured_box():
    anchor = json.loads((HERE / cfg["SLD_anchor"]).read_text(encoding="utf-8"))
    asym, error = anchor["SLD_Ae_1996_98"]
    widths = cfg["widths"]
    smin = from_asymmetry(asym + widths*error)
    smax = from_asymmetry(asym - widths*error)
    gv, ge = cfg["CHARM_II_gV"]
    return smin, smax, gv - widths*ge + .5 - 2*smax, gv + widths*ge + .5 - 2*smin


def budget_bounds(delta_s, delta_v):
    smin, _, _, upper = measured_box()
    upper += 2 * delta_s + delta_v
    r, z = cfg["r_upper"], cfg["z_upper"]
    alpha = (1 - 2*r*(1 - smin + delta_s)) * (1 - r/(1-z))
    beta = r*r/(1-z)
    assert alpha - upper*beta > 0
    dmax = upper / (alpha - upper*beta)
    v = 1 / math.sqrt(math.sqrt(2)*cfg["muon_GeVm2"])
    return {"difference_max": upper, "d_max": dmax, "sqrt_C_times_X_GeV_min": v/(2*math.sqrt(dmax))}


def physical_readouts(g4, gr, gl, breaking_norm, v):
    gb = math.sqrt(1.5)*g4
    total = gb*gb + gr*gr
    r = gr*gr/total
    gy = gr*gb/math.sqrt(total)
    gz = math.hypot(gy, gl)
    s = (gy/gz)**2
    m0sq = total*breaking_norm
    a = (gz*v)**2/4
    h = math.sqrt(total)*r/gz
    mass = np.array([[a, -a*h], [-a*h, m0sq + a*h*h]])
    eigen, vec = np.linalg.eigh(mass)
    # (T3_L, hypercharge, half B-L) for nuL, eL, eR.
    quantum = np.array([[.5,-.5,-.5], [-.5,-.5,-.5], [0.,-1.,-.5]])
    weak, y, charge = quantum.T
    currents = np.column_stack((gz*(weak-s*(weak+y)), math.sqrt(total)*(r*y-charge)))
    low = v*v * currents @ np.linalg.solve(mass, currents.T)
    # Low-energy Fermi coefficient of nuL with each electron chirality is gV +/- gA.
    gv = (low[0,1] + low[0,2])/2
    ga = (low[0,1] - low[0,2])/2
    pole = currents @ vec[:,0]
    ratio = (pole[1]+pole[2])/(pole[1]-pole[2])
    eff = (1-ratio)/4
    d, z = v*v/(4*breaking_norm), eigen[0]/m0sq
    diff = gv+.5-2*eff
    predicted = d*(1-2*r*(1-eff))*(1-z-r)/(1-z+r*r*d)
    np.testing.assert_allclose([diff,ga,gv],
        [predicted,-.5,-.5+2*s+d*(1-2*r*(1-s))],rtol=3e-11,atol=1e-12)
    return {"r":r,"d":d,"z":z,"pole_angle":eff,"gV":gv,"gA":ga,"difference":diff}


def evaluate():
    vals = {name:budget_bounds(*shifts) for name,shifts in cfg["transfer_sets"].items()}
    gauge = json.loads((HERE/cfg["SM_gauge_anchor"]).read_text(encoding="utf-8"))
    gy = math.sqrt(.6)*gauge["g1_GUT"][0]
    g4, gl = gauge["g3"][0], gauge["g2"][0]
    r = gy**2/(1.5*g4**2)
    gr = math.sqrt(gy**2/(1-r))
    v = 1/math.sqrt(math.sqrt(2)*cfg["muon_GeVm2"])
    examples = {str(w):physical_readouts(g4, gr, gl, w*w, v) for w in (270, 450, 3000)}
    assert examples["270"]["difference"] > vals["finite"]["difference_max"]
    assert examples["3000"]["difference"] < vals["ideal"]["difference_max"]
    # Without the independent color normalization a blind gV direction exists for every heavy norm.
    s = .231
    rb = 1/(2*(1-s))
    gb = gy/math.sqrt(rb)
    g4blind = gb/math.sqrt(1.5)
    grblind = gy/math.sqrt(1-rb)
    glblind = gy*math.sqrt((1-s)/s)
    blind = physical_readouts(g4blind,grblind,glblind,600**2,v)
    assert abs(blind["difference"]) < 1e-12 and blind["d"] > .03
    # Algebra also holds outside fit windows; independently varied gauge angles and mass ratios.
    for g4x in (.4,.9,1.25):
        for grx in (.27,.61,1.4):
            for wx in (380,780,5000):
                physical_readouts(g4x,grx,gl,wx*wx,v)
    return {"two_width_box":measured_box(),"bounds":vals,"MZ_tree_matching_examples":examples,
            "color_matching_released":blind,"gauge_reference_r":r}


def close(now, saved):
    if isinstance(now,dict):
        assert now.keys()==saved.keys()
        for key,value in now.items():
            close(value,saved[key])
    else:
        np.testing.assert_allclose(now,saved,rtol=2e-9,atol=2e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=evaluate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        close(out,json.loads(path.read_text(encoding="utf-8")))
    print("1544 pole and contact amplitudes, finite data bounds and blind direction verified")
