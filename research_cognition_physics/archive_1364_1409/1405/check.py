"""Two physical clocks, untagged muon lifetime, normalization and complex amplitudes."""
import argparse
import json
from itertools import product
from math import acos, cos, pi, sqrt
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def clocks(g, delta):
    return 1/(g+delta/2), 1/(g-delta/2)


def lifetime(A, light, heavy):
    return ((1+A)*heavy**2+(1-A)*light**2)/((1+A)*heavy+(1-A)*light)


def inverse(t, light, heavy):
    x=(t-light)/(heavy-light)
    return (x/heavy-(1-x)/light)/(x/heavy+(1-x)/light)


def ranges(obs,k,light,heavy,g,dg,delta,d_delta):
    center,minus,sminus,plus,splus=obs
    lo,hi=center-k*(minus+sminus),center+k*(plus+splus)
    a,b=max(light,lo),min(heavy,hi)
    fixed=[inverse(a,light,heavy),inverse(b,light,heavy)] if a<=b else None
    all_clocks=[clocks(g+i*k*dg,delta+j*k*d_delta) for i,j in product((-1,1),repeat=2)]
    def allowed(ix):
        xs=[p[ix] for p in all_clocks]
        return max(lo,min(xs))<=min(hi,max(xs))
    return {"measured_life_window":[lo,hi],"A_window_at_central_clock":fixed,
            "pure_light_with_clock_outer_box":allowed(0),"pure_heavy_with_clock_outer_box":allowed(1)}


def effective_observables(f,A):
    phi=acos(A)/2
    p=sqrt(1-f)*np.exp(1j*phi)
    s=sqrt(f)*np.exp(1j*(phi+pi/2))
    z=p*p-s*s
    # Sum of helicities; the two tagged mixing components are the real/imaginary parts.
    return float(abs(p)**2+abs(s)**2),float(z.real),float(z.imag),float((p.conjugate()*s).real)


def run():
    inp=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    g,dg=inp["clock"]["Gamma"]
    delta,d_delta=inp["clock"]["DeltaGamma"]
    light,heavy=clocks(g,delta)
    y=delta/(2*g)
    dist={}
    for name,obs in inp["lifetime"].items():
        dist[name]={str(k):ranges(obs,k,light,heavy,g,dg,delta,d_delta) for k in inp["window_multipliers"]}
        dist[name]["central_inversion_without_physical_constraint"]=inverse(obs[0],light,heavy)
    assert dist["LHCb_R123"]["1"]["pure_light_with_clock_outer_box"]
    assert dist["LHCb_R123"]["1"]["pure_heavy_with_clock_outer_box"]
    for name in ("CMS_R2","LHCb_R3","LHCb_R123"):
        assert dist[name]["2"]["pure_light_with_clock_outer_box"] and dist[name]["2"]["pure_heavy_with_clock_outer_box"]
    rate=inp["branch_heavy_fit"]["LHCb_R3"][0]
    by_case={}
    grid=np.linspace(0,60.,120001)
    for A,multiplier in zip((1,0,-1),inp["LHCb_R3_branch_template_multipliers_for_plus_zero_minus"]):
        t=lifetime(A,light,heavy)
        density=(1+A)*np.exp(-grid/heavy)/2+(1-A)*np.exp(-grid/light)/2
        tnum=np.trapezoid(grid*density,grid)/np.trapezoid(density,grid)
        assert abs(t-tnum)<1e-7 and abs(inverse(t,light,heavy)-A)<1e-12
        measured=rate*multiplier
        by_case[str(A)]={"time_ps":t,"time_integrated_branch":measured,
                         "production_time_branch":measured*(1-y*y)/(1+A*y)}
    cm,sl=inp["CMS_single_exponential_rate_rescale_intercept_and_slope_ps_inv"]
    cms=inp["branch_heavy_fit"]["CMS_R2"][0]
    # Ratios cancel the rounding of the fitted correction at its nominal lifetime.
    c_ratio=(cm+sl*light)/(cm+sl*heavy)
    f_test=np.linspace(0,1,5)
    observables=[effective_observables(float(f),.6) for f in f_test]
    np.testing.assert_allclose(observables,np.tile([1.,.6,.8,0.],(5,1)),atol=5e-16)
    c=cos(.4)
    fwindow=[max(0.,(c-.6)/(1+c)),min(1.,(1-.6)/(1+c))]
    assert np.allclose([c-(1+c)*fwindow[0],1-(1+c)*fwindow[1]],[.6,.6])
    return {"HFLAV_light_heavy_ps":[light,heavy],"y":y,"HFLAV_pure_angular_light_heavy_ps":clocks(*inp["clock"]["angular_only"]),
            "time_windows":dist,"LHCb_Run3_three_dedicated_refits":by_case,
            "same_integrated_rate_pure_light_over_heavy_production_intensity":heavy/light,
            "LHCb_R3_same_events_pure_light_over_heavy_production_intensity":1.092*heavy/light,
            "CMS_pure_light_over_heavy_branch_from_lifetime_factor":c_ratio,
            "CMS_branch_at_pure_light_referenced_to_heavy":cms*c_ratio,
            "muon_helicity_sum_examples":{"example_target_A":.6,"scalar_fractions":f_test.tolist(),
                                         "N_A_ImZ_interference":observables,
                                         "scalar_fraction_with_phase_budget_0p4":fwindow}}


def verify(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            verify(x,y)
    elif a is None or isinstance(a,bool):
        assert a==b
    else:
        np.testing.assert_allclose(a,b,rtol=5e-10,atol=2e-13)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true")
    opts=p.parse_args()
    result=run()
    out=HERE/"results.json"
    if opts.write:
        out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    elif opts.check:
        verify(result,json.loads(out.read_text(encoding="utf-8")))
        print("1405 time-integrated rates, scalar ambiguity, refitted templates and physical bounds verified")
    else:
        print(json.dumps(result,indent=2))
