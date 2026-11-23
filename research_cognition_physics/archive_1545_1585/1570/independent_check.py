"""Timelike spectral measure and all-open threshold packing, independent of the author."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


BASE=Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8-sig"))


def spectral_difference(m,v,q):
    # Optical integral over u>=4m^2; z=sqrt(1-4m^2/u).
    knots,weights=np.polynomial.legendre.leggauss(64)
    z=(knots+1)/2
    a=1-z*z
    mass2=4*m*m
    den=np.prod([mass2+n*q*q*a for n in (1,2,3)],axis=0)
    kernel=4*mass2*q**4*z*z*a/den
    return m**4/(8*math.pi**2*v*v)*np.dot(weights/2,kernel)


def integral_check(m,v,q):
    knots,weights=np.polynomial.legendre.leggauss(48)
    t=(knots+1)/2
    x=t*(1-t)*q*q/m**2
    return m**4/(8*math.pi**2*v*v)*np.dot(weights/2,np.log1p(x*x/((1+x)*(1+3*x))))


def compute():
    base=read("../1372/inputs.json")
    old=read("../1372/results.json")
    v,dv=base["matching"]["v_GeV"]
    h,dh=base["atlas"]["mass_GeV"]
    trace=old["SM_trace_over_v4"]["center"]
    lo,central,hi=[old["bosonic_fourth_root_with_f_equal_v_no_neutral"][k]["Higgs_matching_boundary_GeV"]
                   for k in ("box_low","center","box_high")]
    def deficit(mh):
        return 8*math.pi**2*v*v*mh*mh-trace*v**4
    def ratio(mh):
        return deficit(mh)/(mh/2)**4
    assert math.isclose(deficit(h),central**4,rel_tol=1e-14)
    y=ratio(h)
    rest=y-math.floor(y)
    gamma_floor=(h/2)**4/(8*math.pi*h*v*v)*rest*math.sqrt(1-math.sqrt(rest))
    beta=5e-7
    target_count=math.floor(y)
    target_ratio=target_count*(1-beta*beta)**2
    low,high=h-dh,h+dh
    assert ratio(low)>target_ratio>ratio(high)
    for j in range(80):
        mid=(low+high)/2
        if ratio(mid)>target_ratio:
            low=mid
        else:
            high=mid
    hmatch=(low+high)/2
    assert math.isclose(ratio(hmatch),target_ratio,abs_tol=5e-11)
    gamma=deficit(hmatch)*beta/(8*math.pi*hmatch*v*v)
    real_response={}
    for energy in (h+1,150.,200.):
        imag=central**4/(8*math.pi*v*v)*math.sqrt(1-(h/energy)**2)
        real_response[str(energy)]={"Im_Pi_lower_GeV2":imag,"virtual_width_lower_GeV":imag/energy}
    go=[]
    for step in (100.,250.,500.):
        row={}
        for key,m,vev,n in (("lower",lo,v+dv,1),("central_one",central,v,1),
                             ("central_four",central/math.sqrt(2),v,4)):
            val=n*spectral_difference(m,vev,step)
            assert math.isclose(val,n*integral_check(m,vev,step),rel_tol=6e-13)
            row[key]=val
        go.append(row)
    return {"central_X_over_threshold4":y,"minimum_real_species_at_central_all_open":math.ceil(y),
            "central_infimum_width_GeV":gamma_floor,
            "broad_count_interval":[(lo/((h+dh)/2))**4,(hi/((h-dh)/2))**4],
            "near_threshold_open_example":{"N_real":target_count,"h_GeV":hmatch,
                "m_GeV":hmatch/2*math.sqrt(1-beta*beta),"fixed_beta":beta,
                "tree_width_GeV":gamma,"tree_Br_for_visible_4_MeV":gamma/(.004+gamma)},
            "all_open_absorptive_if_exact_central_inputs":real_response,
            "independent_dispersive_GeV2":go}


def assert_equal(actual,expected):
    if isinstance(actual,dict):
        assert actual.keys()==expected.keys()
        for k in actual:
            assert_equal(actual[k],expected[k])
    elif isinstance(actual,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):
            assert_equal(a,b)
    else:
        assert math.isclose(actual,expected,rel_tol=1e-10,abs_tol=2e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=compute()
    file=BASE/"independent_results.json"
    if args.write_new:
        payload=json.dumps(data,indent=2,allow_nan=False)+"\n"
        with file.open("x",encoding="utf-8") as f:
            f.write(payload)
    else:
        assert_equal(data,read(file))
        print("1570 independent dispersive and timelike threshold checks passed.")
