"""Flavor-resolved, finite one-gluino contribution to the down-quark mass phase."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf8"))


def inputs_to_projectors(cfg):
    entry = (HERE / cfg["family"]).resolve()
    spec = importlib.util.spec_from_file_location("flavor1521", entry)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    frame, mother, right, _, left = mod.family(read(cfg["family_inputs"]))
    top = left.conj().T @ frame[2]
    sterile = right[0].conj()
    norm2 = np.vdot(sterile, sterile).real
    L = np.outer(top,top.conj())
    R = np.outer(sterile,sterile.conj())/norm2
    for p in (L, R):
        np.testing.assert_allclose(p@p,p,atol=1e-12)
    return L, R, norm2, np.linalg.svd(mother@right,compute_uv=False)[::-1]


def invariant(mass, left, right):
    d = np.diag(mass)
    return float(np.trace(np.linalg.solve(d,left @ d @ right)).imag)


def invariant_box(table, left, right, widths):
    low,high = table[:,0]-widths*table[:,1],table[:,0]+widths*table[:,1]
    lo,hi = 0.,0.
    for i,j in ((0,1),(0,2),(1,2)):
        ratios=np.array([low[j]/high[i],high[j]/low[i]])
        terms=(left[i,j]*right[j,i]).imag*(ratios-1/ratios)
        lo+=terms.min()
        hi+=terms.max()
    return [float(lo),float(hi)]


def root_of_increasing(func, target):
    lo,hi=0.,1.
    assert func(lo)<=target<func(hi)
    for _ in range(64):
        mid=(lo+hi)/2
        if func(mid)<target:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def compute():
    cfg=read("inputs.json")
    L,R,sterile_norm2,high=inputs_to_projectors(cfg)
    table=np.array(cfg["down_yukawa_dr_1TeV"])
    ds=table[:,0]
    inv=invariant(ds,L,R)
    box=invariant_box(table,L,R,cfg["box_half_widths"])
    assert box[0]*box[1]>0
    gauge=read(cfg["top_transport"])["sm_dr_at_one_TeV"]["g3"]
    pref=lambda g: 2*(g*g/(4*math.pi))/(3*math.pi)*abs(cfg["real_Aeff_over_M"])
    k=pref(gauge[0])*abs(inv)
    k_low=pref(gauge[0]-cfg["box_half_widths"]*gauge[1])*min(map(abs,box))
    pts,weights=np.polynomial.legendre.leggauss(48)
    z=(1+pts)/2
    w=weights/2
    def mixed(x,y):
        return float(x*y*(w @ (z*z*(1-z)/((1+x*z)*(1+y*z)))))
    def single(x,y):
        return float(w @ ((1-z)/((1+x*z)*(1+y*z))))
    for x in cfg["examine_splits"][1:]:
        np.testing.assert_allclose(single(x,x)-2*single(x,0)+single(0,0),mixed(x,x),
                                   rtol=2e-8,atol=1e-16)
        assert x*x/(12*(1+x)**2) <= mixed(x,x) <= x*x/12
    d=np.diag(ds)
    def full_match(x,y):
        Delta=np.zeros((3,3),complex)
        for a,p in [(0,np.eye(3)-L),(x,L)]:
            for b,q in [(0,np.eye(3)-R),(y,R)]:
                Delta += single(a,b)*(p @ d @ q)
        return np.trace(np.linalg.solve(d,Delta)).imag
    np.testing.assert_allclose(full_match(.27,.41),inv*mixed(.27,.41),atol=1e-12)
    np.testing.assert_allclose([full_match(.6,0),full_match(0,.6)],0,atol=1e-13)
    edm=read(cfg["edm"])["edm"]
    out={}
    slopes={"2023":edm["abs_slope_e_fm"],"2026":cfg["lattice_2026_abs_slope_e_fm"]}
    for name,(s,es,et) in slopes.items():
        lim=edm["limit_90_ecm"]/(edm["fm_in_cm"]*s)
        max_lim=edm["limit_90_ecm"]/(edm["fm_in_cm"]*(s-2*(es+et)))
        out[name]={"theta_cap_central":lim,"theta_cap_two_width_lattice_box":max_lim,
                   "equal_splits_cap_central":root_of_increasing(lambda t:k*mixed(t,t),lim),
                   "equal_splits_cap_generous_box":root_of_increasing(lambda t:k_low*mixed(t,t),max_lim),
                   "right_split_cap_if_left_point_one":root_of_increasing(lambda t:k*mixed(.1,t),lim)}
    return {"invariant_at_one_TeV":inv,"invariant_masses_two_width_outer":box,
            "invariant_if_original_high_scale_mass_ratios_used":invariant(high,L,R),
            "norm_squared_of_projected_parent_kernel":float(sterile_norm2),
            "theta_at_equal_positive_fractional_mass_squared_splits":{
                str(x):k*mixed(x,x) for x in cfg["examine_splits"]},
            "coefficient_of_small_xy":k/12,"edm_single_source":out}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=5e-8,atol=1e-15)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        compare(result,read("results.json"))
    print("Verified the finite gluino threshold, flavor invariant and neutron sensitivity.")
