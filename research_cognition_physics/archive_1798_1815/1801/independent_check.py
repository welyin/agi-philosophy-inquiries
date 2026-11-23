"""All-rephasing charged-mass tests of direct alignment in the 1799 effective model."""

import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def bisect(fun, low, high):
    assert fun(low)<0<fun(high)
    for _ in range(80):
        middle=(low+high)/2
        if fun(middle)<0:
            low=middle
        else:
            high=middle
    return (low+high)/2


def block_ceiling(a, ds, db, uc_over_t, tau_lo, tau_hi, tilt):
    light=a*(ds+db*tilt**2)+(tau_hi+a*db)*uc_over_t
    cross=a*tilt*(db+ds)
    if light>=tau_lo:
        return math.inf
    lower_tt=tau_lo-cross**2/(tau_lo-light)
    if lower_tt<=0:
        return math.inf
    return light+cross**2/lower_tt


def charged(data, width):
    u,d,e,v,_=data
    tilt=math.sqrt(1-abs(v[2,2])**2)
    if width:
        tilt=max(tilt,.05)
    ds,db=d[-2:]*(1+width)
    uc_over_t=(1+width)*u[1]/u[2]
    tau_lo,tau_hi=e[-1]*(1-width),e[-1]*(1+width)
    mu_lo=e[-2]*(1-width)
    def cap(a, slack=0.):
        return block_ceiling(a,ds,db,uc_over_t,tau_lo-slack,tau_hi+slack,tilt)
    edge=bisect(lambda a:cap(a)-mu_lo,1,3)
    shift=bisect(lambda t:cap(1,t)+t-mu_lo,0.,mu_lo)
    projective=shift/((np.linalg.norm(d)+np.linalg.norm(e))*(1+width))
    # Positive H lower bound uses the most favorable endpoints.
    h_floor=(3*d[-1]*(1-width)-tau_hi)/4
    return {"width":width,"upper_light_top_mixing":tilt,"observed_mu_lower":mu_lo,
            "direct_alignment_mu_upper":cap(1),"required_abs_a":edge,
            "required_abs_rH_over_rF":(edge-1)/(edge+3),
            "required_rH_H_over_top_in_operator_norm":(edge-1)*h_floor/(tau_hi+edge*db),
            "direct_E_minus_D_minus_cU_minimum_operator_norm":shift,
            "required_sine_of_U_F_complex_Frobenius_angle":projective}


def evaluate():
    old=runpy.run_path(str(HERE.parent/"1799/check.py"))
    cfg=read("../1799/inputs.json")
    data=old["backend"]()
    for vec,key in zip(data[:3],("up","down","lepton")):
        np.testing.assert_allclose(vec,cfg[key],rtol=3e-12)
    center=charged(data,0.)
    examples={}
    for group in ("mixed_type_I_II","fixed_delta_minus70_mass_readout"):
        for name,raw in cfg[group].items():
            x=np.array(raw)
            U,D,E,F,_,_,_,_=old["matrices"](x,data)
            a,b,rh,rf,Y,_=old["symmetric_sources"](x,data)
            norm=lambda z: np.linalg.norm(z,"fro")
            sine=norm(F-U*np.vdot(U,F)/norm(U)**2)/norm(F)
            direct=norm(Y-F*np.vdot(F,Y)/norm(F)**2)/norm(U)
            np.testing.assert_allclose(sine,direct,rtol=3e-13)
            t=abs(rh)*np.linalg.norm((3*D+E)/4,2)/np.linalg.norm(U,2)
            assert abs(a)>=center["required_abs_a"] and abs(rh/rf)>=center["required_abs_rH_over_rF"]
            assert t>=center["required_rH_H_over_top_in_operator_norm"] and sine>=center["required_sine_of_U_F_complex_Frobenius_angle"]
            examples[name]={"abs_a":abs(a),"abs_rH_over_rF":abs(rh/rf),
                            "rH_H_over_U_operator_norm":t,"sine_U_F_Frobenius":sine}
    spectrum=read("../1668/inputs.json")
    s,sm,_=spectrum["juno_NO_dm21"]
    a,_,ap=spectrum["juno_NO_dm31"]
    r=(s-2*sm)/(a+2*ap)
    u,d,e,v,t=data
    sr= np.linalg.cond(t)**2*u[1]/u[2]
    return {"central_10PeV":center,"separate_10pct_working_box":charged(data,.1),
            "five_existing_mixed_points":examples,"normal_order_neutral_shape":{
                    "condition_T":np.linalg.cond(t),
                    "pure_upper_overall_solar_fraction":sr**2,"registered_observed_lower":r,
                    "extra_light_mass_operator_norm_over_total_minimum":(math.sqrt(r)-sr)/(1+sr)}}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=clean(evaluate())
    raw=json.dumps(result,indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==result
    else:
        dest.write_text(raw,encoding="utf-8")
    print("1801 independent",hashlib.sha256(raw.encode()).hexdigest()[:12])
