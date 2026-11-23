"""Check the charged Hessian and its physical derivative in the explicit R1375 potential."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
def read(p):
    return json.loads((HERE/p).read_bytes())
ph=runpy.run_path(str(HERE/"../1373/check.py"))
ew=read("../1474/inputs.json")
v=ew["v_GeV"]
higgs=125.

def direct_check(m,y,heavy,a=0.):
    x=math.sqrt(1-y*y)
    vd,v7=x*v,y*v/4
    u=np.array([x,y/4,math.sqrt(15)*y/4])
    ell=np.array([-6*v7/vd,1.,0.])
    Jinv=np.diag([1/vd,1/v7,1.])
    K=heavy**2*np.eye(3)+(higgs**2-heavy**2)*np.outer(u,u)
    C=Jinv@(K-m*m*np.outer(ell,ell))@Jinv
    assert np.min(np.linalg.eigvalsh(C))>0
    def charged(h):
        d=(vd+x*h)/vd
        q=np.array([(vd+x*h)**2/2-vd**2/2,(v7+u[1]*h)**2/2-v7**2/2,u[2]*h])
        src=C@q
        scale=m*m*(1+2*a*u[2]*h/v+2*(a*u[2]*h/v)**2)
        p=math.sqrt(10)*v7/vd*d**5
        n=-math.sqrt(6)*v7/vd*d**5
        res=(v7+u[1]*h-v7*d**6)/math.sqrt(2)
        c=5*math.sqrt(2)*v7/vd**2*d**4
        block=scale*np.array([[p*p+n*n+2*res*c,-p,-n],[-p,1.,0.],[-n,0.,1.]])
        block+=np.diag([src[0],src[1],src[1]])
        return block,scale+src[1]
    d=.0003
    at,md=charged(0)
    vals,U=np.linalg.eigh(at)
    derivative=(charged(d)[0]-charged(-d)[0])/(2*d)
    analytic=(higgs**2+(5+2*a*u[2])*vals)/v
    np.testing.assert_allclose(vals,[0,m*m,(m/x)**2],atol=1e-8)
    np.testing.assert_allclose((U.T@derivative@U).diagonal(),analytic,atol=2e-6,rtol=2e-8)
    np.testing.assert_allclose((charged(d)[1]-charged(-d)[1])/(2*d),analytic[1],atol=2e-6)
    return [(1.,m),(1.,m/x),(2.,m),(3.,m),(4.,m),(5.,m)],u[2]


def calculate():
    mw,_,ff=ph["common_inputs"](ph["read"]("inputs.json"))
    loop=ph["loop_functions"]
    sm=loop(higgs,mw)[2]+sum(nc*q*q*loop(higgs,m)[1] for m,nc,q in ff)
    atlas=(.074+2*.012)/(.0860-2*.0010)
    out={}
    for m,y,heavy in [(200.,.44,475.),(300.,.32,1500.)]:
        poles,sh=direct_check(m,y,heavy)
        weights=np.array([q*q*loop(higgs,p)[0] for q,p in poles])
        k=np.array([(5+higgs**2/p**2)/2 for _,p in poles])
        extra=float(weights@k)
        gamma=abs(sm+extra)**2/abs(sm)**2
        cancel=-extra/(sh*sum(weights))
        direct_check(m,y,heavy,cancel)
        assert abs(np.sum(weights*(k+cancel*sh)))<1e-12
        c=m*m/(32*math.pi*v*v)
        u=ew["W_mass_GeV"]**2/m**2
        c*=math.sqrt(1-4*u)*(1-4*u+12*u*u)
        z=15*y*y
        max_ctau=1.973269804e-16/(c*z*m)
        out[str(int(m))]={"y":y,"pole_charges_and_masses_GeV":poles,
                 "amplitude_extra_real":extra,"diphoton_ZZ_double_ratio":gamma,
                 "minimum_abs_other_amplitude_for_two_widths":abs(sm+extra)-abs(sm)*math.sqrt(atlas),
                 "same_spectrum_lock_coefficient_a_for_SM_ratio":cancel,
                 "two_W_width_GeV":c*z*m,"max_proper_decay_length_m":max_ctau,
                 "B_WW_range_with_five_percent_width_cap":[c*z/.05,1.]}
        if m==300:
            out["300"]["DY_pair_times_WW_squared_fb_conditional_5pct"]=[5.6*(c*z/.05)**2,5.6]
    # Independently compare inclusive real Z production. The fit baseline omits Gamma_Z.
    sw=.2315
    gf=1.1663787e-5
    mz=91.1876
    budget=2.4955+2*.0023-(2.49423-2*.00069)
    zc={}
    for charges in ((2,),(2,3,4,5)):
        pref=gf*mz**3/(6*math.sqrt(2)*math.pi)*sum((q-2-sw*q)**2 for q in charges)
        zc[str(charges)]={"massless_width_GeV":pref,
                         "conditional_mass_floor_GeV":mz/2*math.sqrt(1-(budget/pref)**(2/3))}
    floor=56*5/6
    return {"SM_loop_amplitude":[sm.real,sm.imag],"ATLAS_upper_two_width_double_ratio":atlas,
            "conditional_all_closed_Higgs_channels_extra_amplitude_floor":floor,
            "conditional_all_M_above_h_over_2_double_ratio_floor":(floor-abs(sm))**2/abs(sm)**2,
            "examples":out,"Z_extra_width_budget_GeV":budget,"Z_windows":zc}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    result=calculate()
    path=HERE/"independent_results.json"
    payload=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    if parser.parse_args().check:
        assert json.loads(path.read_bytes())==json.loads(payload)
    else:
        path.write_text(payload,encoding="utf-8")
    print("1788 charged scalar eigenstate loop",hashlib.sha256(payload.encode()).hexdigest()[:12])

