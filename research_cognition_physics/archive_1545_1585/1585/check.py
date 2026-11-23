"""Inclusive t->ch and h->bs partial-width transfer at a single pole."""
import argparse
from itertools import product
import json
from math import cos, prod, radians, sqrt
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "1581"
tools = runpy.run_path(str(PRIOR / "check.py"))
read, box = tools["read"], tools["box"]


def tree_inputs():
    source = read(PRIOR / "inputs.json")
    y = read(PRIOR / source["ordinary_yukawa_MZ"])["yukawa"]
    tree = read(PRIOR / source["semileptonic"])
    cb = read(PRIOR / source["Belle_BGL"])["BelleII_Vcb_BGL"]
    gamma = read(PRIOR / source["gamma"])["chargedB"]["central_degrees"]
    a, b = box(tree["K_product"]), box(tree["K_form"])
    us = [a[0]/b[1], a[1]/b[0]]
    return y, tree, cb, gamma, us


def flavor_intervals(y, tree, cb, us):
    rho = box(tree["ub_over_cb"])
    cb = box(cb)
    V, trig = tools["abs_ckm"](us, cb, rho)
    s12, c12, _, c23, _, c13 = trig
    h_lo = c12[0]/c13[1] - s12[1]*c23[1]*rho[1]
    h_hi = c12[1]/c13[0] + s12[1]*c23[1]*rho[1]
    eta_u = V[1][0][1]*V[2][0][1]/(cb[0]*V[2][2][0])
    eta_d = us[1]*rho[1]/V[2][2][0]
    u, c, t = map(box, y["up"])
    d, s, b = map(box, y["down"])
    nlow = b[0]-s[1]-(d[1]+s[1])*eta_u
    nhigh = b[1]+s[1]+(d[1]+s[1])*eta_u
    dlow = (t[0]-c[1])*h_lo-(u[1]+c[1])*eta_d
    dhigh = (t[1]+c[1])*h_hi+(u[1]+c[1])*eta_d
    assert min(nlow, dlow) > 0
    return [nlow/dhigh, nhigh/dlow], dict(h=[h_lo,h_hi], eta_u=eta_u, eta_d=eta_d)


def phase_space(x,y):
    return sqrt((1-(x+y)**2)*(1-(x-y)**2))


def kin_exact(mt,mh,mc,mb,ms,utheta,dtheta):
    rc, rh, rb, rs = mc/mt, mh/mt, mb/mh, ms/mh
    t = phase_space(rh,rc)*(1+rc**2-rh**2+2*rc*utheta)
    h = phase_space(rb,rs)*(1-rb**2-rs**2-2*rb*rs*dtheta)
    return mt*t/(12*mh*h)


def kin_bound(top, mh, mc_hi, bs_hi):
    rt0, rt1 = mh[0]/top[1], mh[1]/top[0]
    c = mc_hi/top[0]
    tlo = sqrt((1-(rt1+c)**2)*(1-rt1**2))*((1-c)**2-rt1**2)
    thi = (1-rt0**2)*((1+c)**2-rt0**2)
    hlo = (1-(bs_hi/mh[0])**2)**2
    return [top[0]*tlo/(12*mh[1]), top[1]*thi/(12*mh[0]*hlo)]


def central_sign_examples(y,tree,cb,gamma,params,shape,kin):
    central_us = tree["K_product"][0]/tree["K_form"][0]
    table = []
    for offset in (0,180):
        v = np.array(tools["ckm"](central_us, cb[0], tree["ub_over_cb"][0], radians(gamma+offset)))
        y_u = np.array([i[0] for i in y["up"]])
        y_d = np.array([i[0] for i in y["down"]])
        ratios = []
        for signs in product((-1,1), repeat=6):
            su,sd = signs[:3],signs[3:]
            xu = v @ np.diag(y_d*sd) @ v.conj().T
            xd = v.conj().T @ np.diag(y_u*su) @ v
            amp_ratio = abs(xu[1,2]/xd[1,2])
            assert shape[0] <= amp_ratio <= shape[1]
            for phi in (0,.4,np.pi/2):
                k = kin_exact(params["kinematic_top_GeV"][1],params["h_mass_central_GeV"],
                              params["kinematic_charm_GeV"][1],4.2,.1,
                              su[1]*su[2]*cos(2*phi),sd[1]*sd[2]*cos(2*phi))
                assert kin[0] <= k <= kin[1]
                ratios.append(k*amp_ratio**2)
        table.append(dict(gamma_degree=gamma+offset, absolute_width_ratio=[min(ratios),max(ratios)]))
    return table


def calculate():
    p=read(HERE/"inputs.json")
    y,tree,cb,gamma,us=tree_inputs()
    amp, shape=flavor_intervals(y,tree,cb,us)
    mh=[p["h_mass_central_GeV"]+i*2*p["h_mass_uncertainty_GeV"] for i in (-1,1)]
    mt=p["kinematic_top_GeV"]
    out={}
    kt=[]
    for label,trange in [("fixed_nominal_top",[mt[1],mt[1]]),("outer_top_template_transfer",[mt[0],mt[2]])]:
        k=kin_bound(trange,mh,max(p["kinematic_charm_GeV"]),p["bottom_plus_strange_max_GeV"])
        kt=k
        tree_ratio=[a*b**2 for a,b in zip(k,amp)]
        rr=[a*b for a,b in zip(tree_ratio,p["ordinary_transfer_ratio"])]
        upper={}
        for name,w in p["h_width_upper_GeV"].items():
            base=rr[1]*w/p["top_width_direct_min_GeV"]
            response = sqrt(p["direct_FCNC_limit"]["ATLAS"]*p["top_width_direct_min_GeV"])-sqrt(rr[1]*w)
            upper[name]=dict(top_BR=base, relative_to_direct_ATLAS=p["direct_FCNC_limit"]["ATLAS"]/base,
                             additional_amplitude_to_ATLAS_sqrt_GeV=max(0,response))
        out[label]=dict(kinematic_coeff=k, born_width_ratio=tree_ratio, transferred_ratio=rr,
                        top_limits=upper)
    return dict(round=1585, amplitude_ratio=amp, flavor_controls=shape, scenarios=out,
                central_inputs_all_sign_examples=central_sign_examples(y,tree,cb,gamma,p,amp,kt))


def normalize(x):
    if isinstance(x,dict):
        return {k:normalize(v) for k,v in x.items()}
    if isinstance(x,list):
        return [normalize(v) for v in x]
    return float(f"{x:.12g}") if isinstance(x,(float,np.floating)) else x


def main():
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=normalize(calculate())
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,ensure_ascii=False)
            f.write("\n")
    elif read(path)!=result:
        raise SystemExit("1585 result mismatch")
    print("1585 verified")


if __name__=="__main__":
    main()
