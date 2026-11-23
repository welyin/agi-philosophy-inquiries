"""Independent low-energy current, complex clock geometry and spin-trace checks."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE/name).read_bytes())


def ellipse_support(chi, k, mu, w, intervals=32768):
    """Bracket the farthest distance from the origin by a Lipschitz partition."""
    step = (math.pi/2)/intervals
    theta = np.arange(intervals+1)*step
    axis = 1/w-1
    vals = (1-k)*np.cos(theta)+k*np.sqrt(mu*(1+axis*np.sin(theta+chi)**2))
    lower = float(np.max(vals))
    slope = 1-k+k*math.sqrt(mu)*axis/2
    return lower, lower+slope*step/2+1e-13


def clock_and_flavor():
    d=load("../1406/inputs.json")
    lh=load("../1404/inputs.json")["LHCb_run3_2024_25"]["s_in_1e9"]
    cms=d["CMS"]
    kmin,kmax=d["sm_coefficient_ratio"]["range"]
    sm=d["sm_Bs"]
    b=d["BelleII"]
    records={}
    for n in (1,2):
        bel_lo=(b["mu"]-n*b["width"])/(1+n*(b["total_width_in_1e-6"]+b["tau_width_in_1e-6"])/b["reference_SD_in_1e-6"])
        experiments={
            "CMS":(cms["branching_in_1e-9"]+n*sum(cms["upper_widths"]),.75),
            "LHCb_R3":(lh[0]+n*(lh[2]+lh[4]),.78)}
        for lab,(mu_up,w) in experiments.items():
            mu=mu_up/(sm["branching_in_1e-9"]-n*sm["width"])
            # A right-handed Z companion enters the two hadronic currents with opposite sign.
            rhs=[]
            disk=[]
            aligned=[]
            for k in (kmin,kmax):
                radius=(1-k)+k*math.sqrt(mu/w)
                disk.append(radius**2)
                rhs.append(max(0.,(math.sqrt(bel_lo)-radius)/(2*k)))
                # On the aligned ellipse the squared norm is a concave quadratic in Re P.
                xmax=min(math.sqrt(mu),(1-k)*w/(k*(1-w)))
                aligned.append(((1-k)+k*xmax)**2+k*k*(mu-xmax*xmax)/w)
            rec={"Belle_ratio_lower":bel_lo,"mu_accepted_upper":mu,"unknown_clock_Rnu_upper":max(disk),
                 "aligned_clock_Rnu_upper":max(aligned),"extra_right_hat_min_unknown_clock":min(rhs)}
            if lab=="CMS" and n==2:
                for deg,side in [(81.1,"upper"),(81.25,"lower")]:
                    bounds=[ellipse_support(math.radians(deg),k,mu,w) for k in (kmin,kmax)]
                    idx=1 if side=="upper" else 0
                    edge=max(z[idx] for z in bounds)
                    rec[f"amplitude_{side}_at_{deg}_degrees"]=edge
                    assert (edge<math.sqrt(bel_lo)) if side=="upper" else (edge>math.sqrt(bel_lo))
            records[f"{lab}_{n}"]=rec
    return records


def absolute_amplitudes():
    rare=load("../1406/inputs.json")
    d=load("../1464/inputs.json")
    ext=load("../1404/inputs.json")
    hbar=load("../1407/inputs.json")["hbar_GeV_ps"]
    ref=d["HPQCD_tableI_eq15"]
    b=rare["BelleII"]
    Nref=ref["lambda_ref"]*ref["Xt_ref"]/(ref["alpha_inverse_ref"]*ref["sine_squared_ref"]*math.pi)
    ob=ext["const_GeV"]
    M=ob["M_Bs"]; mm=ob["m_mu"]; gf=ob["GF"]
    base=ob["GF"]**2*mm**2*M*math.sqrt(1-(2*mm/M)**2)/(16*math.pi*hbar)
    f,*ferr=d["fBs_GeV"]
    life,life_sd=ext["tau_ps"]["s_heavy"]
    nurows=load("../1464/results.json")["working_boxes"]
    c=.017
    C=math.hypot(.046,.005)
    k0,k1=rare["sm_coefficient_ratio"]["range"]
    r3=ext["LHCb_run3_2024_25"]["s_in_1e9"]
    cms=rare["CMS"]
    out={}
    for n in (1,2):
        nlo=Nref*math.sqrt((b["mu"]-n*b["width"])/(1+n*ref["hadronic_rate_width"]))
        nhi=Nref*math.sqrt((b["mu"]+n*b["width"])/(1-n*ref["hadronic_rate_width"]))
        f2=(f-n*sum(ferr))**2
        Hmin=base*(life-n*life_sd)*f2*(1-d["mu_rate_transfer_width"])
        for lab,btop,w in (("CMS",cms["branching_in_1e-9"]+n*sum(cms["upper_widths"]),.75),
                           ("LHCb_R3",r3[0]+n*(r3[2]+r3[4]),.78)):
            name=f"{lab}_{n}"
            pmax=math.sqrt(1e-9*btop/(Hmin*w))
            assert abs(nlo-nurows[name]["abs_N"][0])<1e-15 and abs(pmax-nurows[name]["abs_P_max"])<1e-15
            Lmin=(nlo-pmax)/(1-k0)
            Fmin=max(0.,(k0*nlo-pmax)/(1-k0))
            Fmax=(pmax+k1*nhi)/(1-k1)
            Emin=max(0.,(1-c*k0)*Lmin-c*(C+pmax))
            out[name]={"N_range":[nlo,nhi],"P_upper":pmax,"L_lower":Lmin,"F_range":[Fmin,Fmax],
                       "charged_C_min_if_epsilon_zero":(1/c-k0)*Lmin-pmax,
                       "minimum_additive_loop_remainder":Emin,
                       "max_c_allowed_for_incompatibility":Lmin/(C+pmax+k0*Lmin)}
    return out


def spin_trace():
    # Each direction uses a running-Yukawa benchmark identified with tree kinematic masses here.
    eye=np.eye(4,dtype=complex)
    zero=np.zeros((2,2))
    unit=np.eye(2)
    gamma0=np.diag([1.,1.,-1.,-1.])
    pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    gs=[np.block([[zero,s],[-s,zero]]) for s in pauli]
    g5=1j*gamma0@gs[0]@gs[1]@gs[2]
    PL=(eye-g5)/2; PR=(eye+g5)/2
    mw=125.1; mb=2.8; ms=.055
    def width_kernel(a,b):
        E=(mw*mw+a*a-b*b)/(2*mw)
        mom=math.sqrt(E*E-a*a)
        q1=E*gamma0-mom*gs[2]; q2=(mw-E)*gamma0+mom*gs[2]
        Y=b*PR+a*PL
        spin=np.trace((q1+a*eye)@Y@(q2-b*eye)@gamma0@Y.conj().T@gamma0)
        assert abs(spin.imag)<1e-8
        return mom*float(spin.real)
    ratio=(width_kernel(mb,ms)+width_kernel(ms,mb))/width_kernel(mb,mb)
    b,s=(mb/mw)**2,(ms/mw)**2
    lam=(1-b-s)**2-4*b*s
    analytic=math.sqrt(lam)*((b+s)*(1-b-s)-4*b*s)/(b*(1-4*b)**1.5)
    assert abs(analytic-ratio)<1e-12
    return {"benchmark_h_b_s_GeV":[mw,mb,ms],"hbs_pair_over_reference_hbb_per_F_squared":ratio,
            "scalar_over_axial_new_amplitude_bound":(load("../1404/inputs.json")["const_GeV"]["M_Bs"]/mw)**2}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result={"clock_and_flavor":clock_and_flavor(),"absolute":absolute_amplitudes(),"spin":spin_trace()}
    target=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(target.read_bytes())
    else:
        with target.open("x",encoding="utf-8") as file:
            file.write(json.dumps(result,indent=2)+"\n")
    print("1632 independent PASS: clock tilt certificate, flavor projections, charged-current and spin bounds")
