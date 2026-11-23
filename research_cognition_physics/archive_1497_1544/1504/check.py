"""Finite one-loop Higgs amplitudes and weak isospin at fixed charged spectrum."""

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("round1373", HERE.parent / "1373" / "check.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def load(p):
    return p.read_text(encoding="utf-8")


def fg(t):
    if t >= 1:
        angle = math.asin(1/math.sqrt(t))
        return angle**2, math.sqrt(t-1)*angle
    log_cut = 2*math.acosh(1/math.sqrt(t)) - 1j*math.pi
    return -log_cut**2/4, math.sqrt(1-t)*log_cut/2


def integrals(mh, mz, mass):
    a, b = (2*mass/mh)**2, (2*mass/mz)**2
    fa, ga = fg(a)
    fb, gb = fg(b)
    d = a-b
    i1 = a*b/(2*d)+(a*b/d)**2*(fa-fb)/2 + a*a*b*(ga-gb)/d**2
    i2 = -a*b*(fa-fb)/(2*d)
    return i1, i2


def group_traces(multiplets):
    states = [(m/2, (m+y)/2) for j,y in multiplets for m in range(-j,j+1,2)]
    assert sorted(q for _,q in states if q) == [1,2]
    return sum(q*q for _,q in states), sum(t*q for t,q in states)


def ratio(a, b, x):
    return abs(a+b*x)**2/abs(a)**2


def range_ratio(a, b, bounds):
    lo, hi = sorted(bounds)
    vertex = min(hi, max(lo,-a.real/b))
    values = [ratio(a,b,x) for x in (lo,hi,vertex)]
    return [min(values), max(values)]


def interval(measurement):
    value, minus, plus = measurement
    return [value-minus,value+plus]


def evaluate():
    config = json.loads(load(HERE/"inputs.json"))
    source = json.loads(load(HERE/config["baseline"]))
    mw, mz, fermions = prior.common_inputs(source)
    mh, v, mass = source["higgs_mass_GeV"],source["v_MZ_GeV"],config["mass_GeV"]
    cw = mw/mz
    sw = math.sqrt(1-cw*cw)
    tw2 = (sw/cw)**2
    i1w,i2w = integrals(mh,mz,mw)
    tau = (2*mw/mh)**2
    az = cw/sw*(4*(3-tw2)*i2w+((1+2/tau)*tw2-(5+2/tau))*i1w)
    ag = prior.loop_functions(mh,mw)[2]
    for mf,nc,q in fermions:
        i1,i2 = integrals(mh,mz,mf)
        t3 = .5 if q>0 else -.5
        az += 2*nc*q*(t3-2*q*sw**2)/(sw*cw)*(i1-i2)
        ag += nc*q*q*prior.loop_functions(mh,mf)[1]
    a0 = prior.loop_functions(mh,mass)[0]
    i1 = integrals(mh,mz,mass)[0]
    assert az.real>0 and ag.real<0 and a0>1/3 and i1>1/6
    datum = source["cms"]
    g, gl, gh = datum["gammagamma"]
    z, zl, zh = datum["ZZ"]
    measured = [(g-gl)/(z+zh), (g+gh)/(z-zl)]
    central = g/z
    bg = 5*a0

    def roots(r):
        absolute_real = math.sqrt(r*abs(ag)**2-ag.imag**2)
        return ((-ag.real-absolute_real)/bg,(-ag.real+absolute_real)/bg)

    roots_interval = list(zip(roots(measured[0]),roots(measured[1])))
    center_roots = roots(central)
    results = {}
    delta = config["relative_signal_transfer_budget"]
    atlas,cms = interval(config["atlas_zgamma"]),interval(config["cms_zgamma"])
    for name, multiplets in config["sets_twice_isospin_twice_Y"].items():
        qs,qt = group_traces(multiplets)
        bz = -2*i1*(qt-qs*sw**2)/(sw*cw)
        branches = {}
        for label, k, krange in zip(("near","reflected"),center_roots,roots_interval):
            support = range_ratio(az,bz,krange)
            mu = [(1-delta)*support[0],(1+delta)*support[1]]
            branches[label] = {
                "central_width_ratio":ratio(az,bz,k), "box_width_ratio":support,
                "box_mu_with_transfer":mu,
                "overlaps_ATLAS": bool(mu[0]<=atlas[1] and mu[1]>=atlas[0]),
                "overlaps_CMS":bool(mu[0]<=cms[1] and mu[1]>=cms[0])}
        results[name] = {"Q2":qs,"Q_T3":qt,"Z_amplitude_per_kappa":bz,"branches":branches}
    center_z = sorted(row["branches"]["near"]["central_width_ratio"] for row in results.values())
    closest = min(y-x for x,y in zip(center_z,center_z[1:]))
    lambda_h = .5*(mh/v)**2
    lam_bounds = [2*mass**2/v**2*k for k in roots_interval[0]]
    co_positive = max(x*x for x in lam_bounds)/(4*lambda_h)
    assert config["example_scalar_quartic"]>co_positive
    return {"SM_amplitudes_gamma_Z":[[ag.real,ag.imag],[az.real,az.imag]],
            "A0_I1_finite":[a0,i1], "CMS_calibrating_partial_ratio":[central,*measured],
            "kappa_central":center_roots, "kappa_box":[sorted(x) for x in roots_interval],
            "lambda_portal_central":[2*mass**2/v**2*k for k in center_roots],
            "weak_branch_sufficient_positive_quartic":co_positive,
            "ATLAS_and_CMS_working_intervals":[atlas,cms], "representations":results,
            "near_central_common_absolute_halfwidth_three_separated":closest/2,
            "near_central_halfwidth_extremes_separated":(center_z[-1]-center_z[0])/2}


def compare(new, old):
    if isinstance(new,dict):
        assert new.keys()==old.keys()
        for k in new:
            compare(new[k],old[k])
    elif isinstance(new, (list,tuple)):
        assert len(new)==len(old)
        for a,b in zip(new,old):
            compare(a,b)
    elif isinstance(new,bool):
        assert new is old
    else:
        np.testing.assert_allclose(new,old,rtol=3e-9,atol=1e-12)


if __name__=="__main__":
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args = p.parse_args()
    result = evaluate()
    target = HERE/"results.json"
    if args.write:
        serialized = json.dumps(result,indent=2,allow_nan=False)
        with target.open("x",encoding="utf-8") as stream:
            stream.write(serialized+"\n")
    else:
        compare(result,json.loads(load(target)))
    print("1504 verified", result["kappa_central"])
