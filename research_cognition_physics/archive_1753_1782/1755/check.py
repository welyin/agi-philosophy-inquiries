"""Project the same-charge quark masses onto the uncalibrated Cabibbo modulus."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


CFG = read("inputs.json")
TABLE = read(CFG["conditional_MILC_covariance"])["table_IV"]
MEANS = np.array(TABLE["means"][:3])
COV = np.array(TABLE["correlation"])[:3,:3]*np.outer(TABLE["sigmas"][:3], TABLE["sigmas"][:3])
T = math.log(CFG["color_scale_GeV"]/CFG["MZ_GeV"])/(16*math.pi**2)


def interval(measurement, radius):
    center, *uncertainties = measurement
    width = radius*math.hypot(*uncertainties)
    return np.array([center-width,center+width])


def transport(n):
    """Enclose the one-loop flavour flow for either CP sign without using us or ud."""
    ub = interval(CFG["inclusive_ub"], n)
    cb = interval(CFG["inclusive_cb"], n)
    yt, yc, yb, ys = 1., .005, .02, .00035
    w = math.hypot(ub[1],cb[1])*math.exp(1.52*(yt*yt+yc*yc+yb*yb+ys*ys)*T)
    assert w < .05
    up_drift = 1.52*(yb*yb*w+ys*ys)
    epsilon = 2/(35**2-1)
    residual = up_drift+1.5*epsilon*(yt*yt+yc*yc)*w
    lower_slope = 1.5*(.65**2*(1-w*w)-yc*yc)
    upper_slope = 1.5*yt*yt
    umin = ub[0]*math.exp(lower_slope*T)-residual*math.expm1(lower_slope*T)/lower_slope
    umax = ub[1]*math.exp(upper_slope*T)+residual*math.expm1(upper_slope*T)/upper_slope
    delta_us = T*(up_drift+1.52*(yc*yc+yt*yt*w*w/2+umax*(yt*yt*w+yc*yc)))
    logD = 1.5*(ys*ys+yc*yc+yt*yt*w*w)*T
    logB = 1.5*T*np.array([-yt*yt-yc*yc, yb*yb+yc*yc-.65**2*(1-2*w*w)])
    h = abs(read(CFG["electroweak_matching"])["thresholds"]["top_minus_light_one_loop_in_decoupling"])
    KD = np.array([1-h*w*w,1+h*w*w])*np.exp([-logD,logD])
    kblo,kbhi=read(CFG["electroweak_matching"])["thresholds"]["HZ_bs_print_rounding_range"]
    KB = np.array([kblo-2*h*w*w,kbhi])*np.exp(logB)
    return dict(weak_ub=ub, color_ub=[umin,umax], w_upper=w, delta_us=delta_us,
                unified_D_factors=KD, unified_B_factors=KB, log_running_D=logD, log_running_B=logB)


def mass_interval(source,n):
    """Two light-ratio measurements whose absolute scale was set by hadron masses."""
    if source=="QCDSF":
        D = 1/interval(CFG["QCDSF_s_d"],n)[::-1]
    else:
        S = interval(CFG["BMW_s_mud"],n)
        r = interval(CFG["BMW_u_d"],n)
        D = 2/(S*(1+r))[::-1]
    return D,CFG["strong_alternative_b_over_s_working_interval"]


def envelope(source,n,gate,sign):
    flow = transport(n)
    roots = [1.,0.]
    for u,kd,kb in itertools.product(flow["color_ub"],flow["unified_D_factors"],flow["unified_B_factors"]):
        q = u*u
        def residual(s,side):
            c = math.sqrt(1-s*s-q)
            f = math.sqrt((1+c)/2) if gate else 1.
            ad, ab = kd*(c if gate else c*c), sign*q*kb/f
            if source=="MILC":
                a=np.array([-ab,s*s,-ad])
                return a@MEANS + side*n*math.sqrt(a@COV@a)
            ds,bs = mass_interval(source,n)
            values=[s*s-ad*d-ab*b for d,b in itertools.product(ds,bs)]
            return max(values) if side==1 else min(values)
        for j,side in enumerate((1,-1)):
            lo,hi=.1,.35
            assert residual(lo,side)<0<residual(hi,side)
            for _ in range(54):
                v=(lo+hi)/2
                if residual(v,side)>0:
                    hi=v
                else:
                    lo=v
            roots[j] = (min if j==0 else max)(roots[j],(lo+hi)/2)
    low=roots+np.array([-1.,1.])*flow["delta_us"]
    ulo,uhi=flow["weak_ub"]
    ud=np.sqrt(1-low[::-1]**2-np.array([uhi,ulo])**2)
    return dict(color_us=roots, MZ_us=low, MZ_ud=ud, MZ_us_over_ud=low/ud[::-1])


def pulls(row):
    us,ub,_ = [row["low_moduli"][k] for k in ("Vus","Vub","Vcb")]
    ud=math.sqrt(1-us*us-ub*ub)
    def shift(x,m):
        return (x-m[0])/math.hypot(*m[1:])
    o=read("../1484/inputs.json")
    reactor,minus,plus=o["DayaBay_sin2_2theta13"]
    offset=row["angle_readings"][0]-reactor
    return dict(Kl3=shift(us,CFG["holdout_Kl3"]), axial=shift(us/ud,CFG["holdout_axial"]),
                nuclear=shift(ud,CFG["holdout_nuclear"]),reactor=offset/(plus if offset>0 else minus))


def calculate():
    flow={str(n):transport(n) for n in CFG["marginal_radii"]}
    bands={}
    for source in ("QCDSF","BMW","MILC"):
        bands[source]={}
        for name,gate,sign in (("signed_gate",True,1),("no_gate",False,1),("other_bottom_sign",True,-1)):
            bands[source][name]={str(n):envelope(source,n,gate,sign) for n in CFG["marginal_radii"]}
    shoot=runpy.run_path(str(HERE/"joint.py"))["shoot"]
    ub=CFG["inclusive_ub"]; cb=CFG["inclusive_cb"]
    targets=dict(QCDSF=1/CFG["QCDSF_s_d"][0],BMW=2/(CFG["BMW_s_mud"][0]*(1+CFG["BMW_u_d"][0])),
                 MILC=MEANS[2]/MEANS[1])
    rows={}
    for source,target in targets.items():
        trials=[("central",0,0)]
        if source=="QCDSF":
            trials.extend((f"ub{a:+d}_cb{b:+d}", a,b) for a,b in itertools.product((-1,1),repeat=2))
        cases={}
        for label,a,b in trials:
            row=shoot(HERE,target,ub[0]+a*math.hypot(*ub[1:]),cb[0]+b*math.hypot(*cb[1:]),.25)
            row["leave_out_pulls"]=pulls(row)
            cases[label]=row
        rows[source]=cases
    # Give an attainable clean BMW corner. This is a feasibility scan of two pre-existing
    # marginal widths, leaving direct first-row weak data unopened by the solver.
    Slo=interval(CFG["BMW_s_mud"],1)[0]
    rlo=interval(CFG["BMW_u_d"],1)[0]
    target=2/(Slo*(1+rlo))
    row=shoot(HERE,target,ub[0],cb[0],.25)
    row["leave_out_pulls"]=pulls(row)
    rows["BMW"]["both_lattice_inputs_minus_one_width"]=row
    fine=shoot(HERE,targets["QCDSF"],ub[0],cb[0],.125)
    diff={k:np.max(np.abs(np.asarray(fine[k])-np.asarray(rows["QCDSF"]["central"][k])))
          for k in ("trial_us","bottom_over_s","sun_calibrated_dm31_eV2","angle_readings","gamma_deg")}
    assert diff["trial_us"]<1e-7
    return dict(transport=flow,bands=bands,calibrating_strong_ratio=targets,full_two_threshold_shooting=rows,
                half_step_QCDSF_central=diff)


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(np.ndarray,list,tuple)): return [tidy(v) for v in x]
    if isinstance(x,(np.floating,float)): return float(f"{x:.10g}")
    if isinstance(x,np.bool_): return bool(x)
    return x


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--check",action="store_true")
    verify=args.parse_args().check
    data=tidy(calculate())
    if verify:
        assert data==read("results.json")
        print("1755 OK: independent mass-to-Cabibbo intervals, split matter solutions and reactor holdout")
    else:
        (HERE/"results.json").write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
