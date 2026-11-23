"""Joint radiative quadratic and finite stop-quartic budget on the frozen 1699 interval."""
from pathlib import Path
import argparse
import json
import math
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2


def read(path):
    return json.loads((HERE/path).read_bytes())


def tangent_state(runner, low, t):
    z = low.copy()
    z[:3]=1/np.sqrt(1/z[:3].real**2-np.array([0,2,3])/(48*math.pi**2))
    c=1/math.hypot(1,t)
    z[3]/=t*c
    z[4]/=c
    z[6:9]/=c
    z[9:]/=t*c
    return z


def couplings(z):
    return np.array([abs(z[3])**2,abs(z[4])**2,np.linalg.norm(z[6:9])**2,np.linalg.norm(z[9:])**2])


def integrate(runner,state,events,step):
    z=np.r_[state, np.zeros(4,complex)]
    def rhs(y,scale):
        return np.r_[runner["beta"](y[:21],scale,True),
                     np.array([6,6,2,2])*couplings(y[:21])/LOOP]
    for a,b in zip(events,events[1:]):
        chunks=max(1,math.ceil(math.log(b/a)/step))
        dl=math.log(b/a)/chunks
        at=math.sqrt(a*b)
        for _ in range(chunks):
            k1=rhs(z,at)
            k2=rhs(z+dl*k1/2,at)
            k3=rhs(z+dl*k2/2,at)
            k4=rhs(z+dl*k3,at)
            z+=dl*(k1+2*k2+2*k3+k4)/6
    return z[:21],z[21:].real


def compute():
    cfg=read("inputs.json")
    prev=runpy.run_path(str(HERE/cfg["previous_running"]))
    slopes=read(cfg["previous_interface"])["soft_interface"]
    q=read(cfg["quadratic_target"])["results_by_other_quartic_allowance"]
    soft,high=cfg["source_range_GeV"]
    budget=cfg["budget"]
    lo,hi=slopes["tan_at_max_stop"][1],slopes["tan_at_zero_stop"][1]
    z=prev["lower"]()
    events=sorted(set([soft,high]+[m for m,_ in prev["thresholds"] if soft<m<high]
                     +[m for m in prev["spec"]["neutral_scalar_poles_if_common_shift_GeV"] if soft<m<high]))
    gauge=tangent_state(prev,z,lo)[:3].real
    cascade=[gauge.copy()]
    for left,right in zip(events,events[1:]):
        coeff=np.array([33/5,1.,-3.])
        for m,delta in prev["thresholds"]:
            if m<math.sqrt(left*right):
                coeff+=delta
        gauge=1/np.sqrt(1/gauge**2-2*coeff*math.log(right/left)/LOOP)
        cascade.append(gauge.copy())
    caps=cfg["certified_coupling_box"]
    gl=np.array(caps["gauge_lo"])
    gh=np.array(caps["gauge_hi"])
    assert np.all(np.min(cascade,axis=0)>=gl) and np.all(np.max(cascade,axis=0)<=gh)
    top,bot,tau,ntr=(caps[key] for key in ["top_max","bottom_max","tau_max","neutrino_trace_max"])
    dt=math.log(high/soft)/LOOP
    ylo=couplings(tangent_state(prev,z,lo))
    yhi=couplings(tangent_state(prev,z,hi))
    assert ylo[0]<top**2 and yhi[1]<bot**2 and yhi[2]<tau**2
    # At any of the first three exits its beta points inward.
    top_face=6*top**2+bot**2+ntr-np.dot([13/15,3,16/3],gl**2)
    bottom_face=6*bot**2+top**2+tau**2-np.dot([7/15,3,16/3],gl**2)
    tau_face=4*tau**2+ntr+3*bot**2-np.dot([9/5,3,0],gl**2)
    ntrace_bound=ylo[3]*math.exp(2*dt*(4*ntr+tau**2+3*top**2))
    assert max(top_face,bottom_face,tau_face)<0 and ntrace_bound<ntr
    coeff=6*dt*top**2
    nu_allow=2*dt*ntr*budget["nu_Wilson_mass_and_trilinear_per_trace"]
    force=(budget["gaugino_abs_over_soft"]*np.dot([26/15,6,32/3],gh**2)
           +2*bot**2*budget["bottom_abs_A_over_soft"]
           +2*ntr*budget["nu_trilinear_Frobenius_over_Yukawa_Frobenius_and_soft"])
    gain=dt*force
    msum=budget["top_average_mass_sum_over_soft2"]
    finite=q["0"]["mechanisms"]["actual_neutral_soft_CW"]["max_abs_Hu_mass_squared_GeV2"]/soft**2
    gap=q["0"]["gap_GeV2"]/soft**2
    stopslope=3*slopes["yt"]**4/LOOP
    cot=budget["abs_mu_over_soft"]/lo

    def max_run(u):
        return coeff*(msum+u*u+u*gain+gain*gain/3)+nu_allow
    def requested(u,extra):
        mix=min(math.sqrt(6),u+gain+cot)
        shift=stopslope*mix**2*(1-mix**2/12)
        return gap*math.sqrt(max(0,(slopes["lambda"]-shift-extra)/slopes["D"]))-finite
    choices={}
    for key in ("0","0.01"):
        extra=float(key)
        lower=q[key]["required_diagonal_difference_GeV2"]/soft**2-finite
        rms2=(lower-nu_allow)/coeff-msum
        u0=0.
        u1=4.
        for _ in range(65):
            m=(u0+u1)/2
            if max_run(m)>=requested(m,extra):
                u1=m
            else:
                u0=m
        lack=requested(0,extra)-max_run(0)
        choices[key]=dict(
            min_full_diagonal_GeV2=q[key]["required_diagonal_difference_GeV2"],
            min_top_A_RMS_over_Ms=math.sqrt(rms2),
            lower_UV_At_over_Ms_using_quartic_same_A=u0,
            max_radiative_Delta_over_Ms2_at_zero_UV_At=max_run(0),
            min_quartic_Delta_after_finite_over_Ms2_at_zero_UV_At=requested(0,extra),
            unfilled_GeV2_for_zero_UV_At=lack*soft**2,
            minimum_weighted_mass_sum_over_Ms2_for_zero_UV_At=msum+lack/coeff,
            minimum_positive_S_log_average_over_Ms2_if_only_escape=lack/((6/5)*gh[0]**2*dt))
    examples={}
    err=0.
    for tval in (lo,3,hi):
        state=tangent_state(prev,z,tval)
        final,measures=integrate(prev,state,events,.05)
        _,small=integrate(prev,state,events,.025)
        err=max(err,float(max(abs(small-measures))))
        actual_request=gap*(tval*tval-1)/(tval*tval+1)-finite
        need_sq=(actual_request-budget["nu_Wilson_mass_and_trilinear_per_trace"]*small[3])/small[0]-msum
        examples[str(round(tval,6))]=dict(Y_top_bottom_UV=final[3:5].real.tolist(),
            log_exposures_top_bottom_tau_nu=small.tolist(),min_top_A_yukawa_weighted_RMS_over_Ms=math.sqrt(need_sq))
    return dict(tan_beta_range=[lo,hi],gauge_extrema_segmentwise=[np.min(cascade,axis=0).tolist(),np.max(cascade,axis=0).tolist()],
                Ynu_trace_uniform_upper_from_Gronwall=ntrace_bound,
                log_length_over_16pi2=dt,upper_log_top_exposure=coeff,upper_integrated_nu_soft_source=nu_allow,
                normalized_max_trilinear_drive=force,integrated_At_drive_upper_in_Ms=gain,
                quadratic_budget_with_same_trilinear=choices,
                running_exposure_controls=examples,max_control_step_refinement_absolute=err)


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=tidy(compute())
    if args.check:
        assert data==read("results.json")
        print("1700: coupled Yukawa and soft-mass budgets, UV trilinear floor and quartic verified")
    else:
        (HERE/"results.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
