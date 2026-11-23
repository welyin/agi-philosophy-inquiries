"""Inverse eta data, exterior-power mass bound and direct singular-value roots."""
import argparse
import json
from math import exp, hypot, log, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def mixing(values):
    s13,s23=values[1:3]
    s12=values[0]
    c12,c13,c23=np.sqrt(1-np.asarray([s12,s13,s23])**2)
    e=np.exp(1j*values[3])
    return np.array([[c12*c13,s12*c13,s13/e],
                     [-s12*c23-c12*s23*s13*e,c12*c23-s12*s23*s13*e,s23*c13],
                     [s12*s23-c12*c23*s13*e,-c12*s23-s12*c23*s13*e,c23*c13]])


def bisection(f,lo,hi):
    val_lo=f(lo)
    assert val_lo*f(hi)<0
    for _ in range(70):
        mid=(lo+hi)/2
        if f(mid)*val_lo>0:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def run():
    cfg=read("inputs.json")
    weak=read(cfg["low"])
    vu=np.array(weak["quark_sines_and_phase"])
    V=mixing(vu[:,0])
    up=np.array(weak["yukawa"]["up"])
    down=np.array(weak["yukawa"]["down"])
    def log_residual(ln_x):
        ys=np.array([exp(ln_x),*up[1:,0]])
        singular=np.linalg.svd(ys[:,None]*V*down[:,0],compute_uv=False)
        return float(3*log(singular[1])-sum(np.log(ys))-sum(np.log(down[:,0])))
    intervals=[(1e-7,1e-4),(1e-4,10),(10,1e8)]
    roots=[exp(bisection(log_residual,log(a),log(b))) for a,b in intervals]

    s0,estat,esys=cfg["ms_over_average"]
    s_lo,s_hi=s0-2*hypot(estat,esys),s0+2*hypot(estat,esys)
    q0,err=cfg["eta_q"]
    def r(q,s):
        return (4*q*q-s*s+1)/(4*q*q+s*s-1)
    r_lo,r_hi=r(q0-2*err,s_hi),r(q0+2*err,s_lo)
    c_lo=up[1,0]-2*up[1,1]
    c_hi=up[1,0]+2*up[1,1]
    t_lo=up[2,0]-2*up[2,1]
    strange_hi=down[1,0]+2*down[1,1]
    b_lo=down[2,0]-2*down[2,1]
    ud2_upper=1-(vu[0,0]-2*vu[0,1])**2
    tb2_lower=(1-(vu[1,0]+2*vu[1,1])**2)*(1-(vu[2,0]+2*vu[2,1])**2)
    def rho_cap(delta):
        rr=r_lo*(1-delta)
        rr_upper=r_hi*(1+delta)
        assert 0<rr<=rr_upper<1
        ds=2/(s_lo*(1-delta)*(1+rr))
        x_max=strange_hi*ds*rr_upper
        frob=ud2_upper+ds**2+(x_max/c_lo)**2+(x_max/t_lo)**2
        return c_hi**2*(s_hi*(1+delta))**2/(4*t_lo*b_lo)*(1+rr)**2/rr*(frob/tb2_lower)**1.5
    threshold=bisection(lambda delta:rho_cap(delta)-1,0,.12)
    assert rho_cap(.08)<1 and rho_cap(.09)>1
    return {"central_low_direct_svd_roots":roots, "r_QCD_interval_from_eta_and_BMW":[r_lo,r_hi],
            "rho_upper_if_two_ratio_transfers_have_2_percent_budget":rho_cap(.02),
            "rho_upper_if_two_ratio_transfers_have_8_percent_budget":rho_cap(.08),
            "symmetric_transfer_half_width_where_this_certificate_expires":threshold}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=run()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    else:
        expected=read("independent_results.json")
        assert out.keys()==expected.keys()
        for key,value in out.items():
            assert np.allclose(value,expected[key],rtol=2e-8,atol=1e-15)
        print("1440 direct singular roots, eta inverse and exterior-norm continuous box verified")
