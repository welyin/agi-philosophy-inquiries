"""Abelian threshold needed by the negative radiative electron source."""
import argparse
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
MZ=91.1876
VEV=246.22


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def running(yukawa_inputs, flavor, old):
    gauge=[yukawa_inputs[k] for k in ("g1_GUT","g2","g3")]
    lower=[4*math.pi/(c+2*d)**2 for c,d in gauge]
    upper=[4*math.pi/(c-2*d)**2 for c,d in gauge]
    cs_cap=-old["recoil_ae_rectangles_one_two"]["Cs"][1][1]
    eta=old["scenarios"]["finite_transfers"]["nearest_kernel_cap_GeV_minus2"]
    kernel_floor=(cs_cap-flavor["stress_ae_additive"])/flavor["me_GeV"]**2-2*eta
    fermion_cap=math.sqrt(1.2*.4/kernel_floor)
    count=8
    rho=(upper[0]-lower[1])/(lower[1]-upper[2])
    target=109/15-23*rho/6
    doublet=(1+rho)/6-.1
    time_min=2*math.pi*(lower[1]-upper[2])/(23/6+count/6)
    loss=.8*(1-math.log(fermion_cap/MZ)/time_min)
    gap=target+loss-count*doublet
    assert gap>0
    limits={"negative_kernel_lower_abs":kernel_floor,"charged_fermion_cap_GeV":fermion_cap,
            "weak_color_time_min_for_up_to_eight":time_min,"strongest_gain_per_doublet":doublet,
            "required_base_slope_lower":target,"hypercharge_cost_lower":loss,
            "slope_shortfall_with_eight":gap,"inverse_coupling_shortfall_at_fixed_23":gap*time_min/(2*math.pi)}
    return np.array([4*math.pi/c**2 for c,d in gauge]),limits


def make_spectra():
    cfg=read(HERE.parent/"1369"/"inputs.json")
    flav=read(HERE.parent/"1591"/"inputs.json")
    old=read(HERE.parent/"1591"/"results.json")
    inv,limits=running(cfg,flav,old)
    # At the upper scale, both independent residuals (12,23) vanish.
    # Solve for t and the positive singlet's log mass relative to MZ directly.
    b=np.array([4.1,-19/6,-7])
    d=np.array([.1,1/6,0])
    q_e=np.array([.8,0,0])
    q_mu=np.array([.2,0,0])
    total=b+9*d+q_e+q_mu
    constants=2*math.pi*inv+9*d*math.log(100/MZ)+q_e*math.log(400/MZ)
    matrix=np.column_stack([np.diff(total),-np.diff(q_mu)])
    t, log_s = np.linalg.solve(matrix,np.diff(constants))
    mx=MZ*math.exp(t)
    s=MZ*math.exp(log_s)
    final=(constants-t*total+log_s*q_mu)/(2*math.pi)
    assert np.ptp(final)<1e-11 and final[0]>1 and 400<s<mx
    meet={"doublet_count":9,"active_lepton_sources":2,"inert_doublets":7,
          "scale_GeV":mx,"muon_hypercharged_scalar_threshold_GeV":s,
          "common_inverse_gauge_coupling":float(final[0])}
    return limits,meet,flav


def charged_kernels(fermion,low_sq,high_sq,mass):
    point,weight=np.polynomial.legendre.leggauss(72)
    t,w=(1+point)/2,weight/2
    low=(1-t)*fermion**2+t*low_sq
    high=(1-t)*fermion**2+t*high_sq
    difference=high_sq-low_sq
    integral=np.dot(w,np.log1p(t*difference/low))/difference
    kernel=-2*np.dot(w,t*(1-t)**2/(low*high))/integral
    return float(integral),float(kernel*mass**2)


def diagonal_channel(doublet,singlet,split,fermion,mass,positive):
    # For arbitrary Higgs radius r, determinant equals doublet^2 singlet^2.
    B=split*doublet*singlet
    portal=2*split**2*doublet**2/VEV**2
    D=doublet**2*(1+split**2)
    S=singlet**2
    high=(D+S+math.hypot(S-D,2*B))/2
    low=doublet**2*S/high
    if positive:
        # Electromagnetically neutral Dirac pole can be fixed to the heavy scalar pole;
        # evaluate the analytic x2=1 limit without subtracting equal doubles.
        fermion=math.sqrt(high)
        with localcontext() as context:
            context.prec=65
            x=Decimal(str(low))/Decimal(str(high))
            f=x*x.ln()/(x-1)
            dip=(3*x-1)/(1-x)**2+2*x**2*x.ln()/(1-x)**3+2*(1/(1-x)+x*x.ln()/(1-x)**2)
            I=float((1-f)/Decimal(str(high-low)))
            da=float((dip-Decimal(1)/3)/(1-f))*mass**2/high
    else:
        I,da=charged_kernels(fermion,low,high,mass)
    yprod=16*math.pi**2*mass/(fermion*B*I)
    assert 0<yprod<1 and low>0
    return {"unbroken_doublet_GeV":doublet,"unbroken_scalar_singlet_GeV":singlet,
            "Higgs_portal":portal,"soft_cubic_GeV":math.sqrt(2)*B/VEV,
            "mixed_scalar_poles_GeV":[math.sqrt(low),math.sqrt(high)],
            "Dirac_pole_GeV":fermion,"Yukawa_product_for_observed_mass":yprod,
            "leading_delta_a":da,"sin_2theta":2*B/(high-low)}


def compute():
    limits,meet,flav=make_spectra()
    e=diagonal_channel(100.,160.,.05,400.,flav["me_GeV"],False)
    mu=diagonal_channel(100.,meet["muon_hypercharged_scalar_threshold_GeV"],
                        .25,None,flav["mmu_GeV"],True)
    lower,upper=read(HERE.parent/"1591"/"results.json")["recoil_ae_rectangles_one_two"]["Cs"][1]
    assert lower<e["leading_delta_a"]<upper and 0<mu["leading_delta_a"]<164e-11
    return {"historical_Cs_two_active_sources":limits,"central_unification_witness":meet,
            "electron_source":e,"muon_source":mu,"flavor_changing_dipoles":0}


def close(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(close(i,j) for i,j in zip(a,b))
    return math.isclose(a,b,rel_tol=3e-11,abs_tol=1e-38)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=compute()
    dest=HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        assert close(result,read(dest))
    print(json.dumps(result))
