"""Eliminate a complex charged-Higgs coupling in K/pi and B/pi leptonic decays."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy


HERE = Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())

CFG = read("inputs.json")
W = read(CFG["weak_witness"])
ROW = read(CFG["flavour"])["exploratory_k3"][CFG["flavour_choice"]]
YUK = read(CFG["yukawas"])["yukawa"]
MASSES = list(CFG["meson_masses_GeV"][key] for key in ("pi","K","B"))


def belle_reference(h, lower=True):
    sign = -1 if lower else 1
    def endpoint(pair):
        return pair[0] + sign*h*pair[1]
    decay = CFG["decay_constants_MeV"]
    ratio = endpoint(decay["B"]) / (decay["pi"][0]-sign*h*decay["pi"][1])
    pi,K,B = MASSES
    mu,tau = CFG["mu_tau_masses_GeV"]
    def kin(meson,lepton):
        return meson*lepton**2*(1-(lepton/meson)**2)**2
    return (1e-6 * endpoint(CFG["B_lifetime_ps"]) * endpoint(CFG["pi_mu_width_per_microsec"])
            * (W["weak_Vub"]/W["weak_Vud"]*ratio)**2 * kin(B,tau)/kin(pi,mu))


def charged_vectors(q_ratios, propagator_mass=math.inf):
    u,d=[],[]
    for m,q in zip(MASSES,q_ratios):
        kin = m*m/(1-(m/propagator_mass)**2)
        u.append(kin*q/(1+q))
        d.append(kin/(1+q))
    return u,d


def eliminate(u,d):
    w=(d[1]-d[0])/(d[2]-d[0])
    C=(u[1]-u[0])-w*(u[2]-u[0])
    assert 0 < w < 1
    return w,C


def parameter_solution(u,d,r,T):
    w,C=eliminate(u,d)
    ahat=((1-w)+w*T-r)/C
    bhat=(1-T-(u[2]-u[0])*ahat)/(d[2]-d[0])
    n=1+u[0]*ahat+d[0]*bhat
    a,b=ahat/n,bhat/n
    A=[1-x*a-y*b for x,y in zip(u,d)]
    assert abs(A[1]/A[0]-r)<1e-12 and abs(A[2]/A[0]-T)<1e-12
    return {"a_GeV_inv2":a,"b_GeV_inv2":b,"a_over_Api_GeV_inv2":ahat,"A_pi":A[0]}


def envelope(ratios):
    # Enclose simultaneous +/-2% changes in each pure-QCD quark ratio, and
    # the entire single-pole propagator for M >= the named floor.
    frac=CFG["independent_quark_ratio_fraction"]
    floor=CFG["alternative_mass_floor_GeV"]
    low,high,dl,dh=[],[],[],[]
    for m,q0 in zip(MASSES,ratios):
        qlo,qhi=q0*(1-frac),q0*(1+frac)
        factor=1/(1-(m/floor)**2)
        low.append(m*m*qlo/(1+qlo))
        high.append(factor*m*m*qhi/(1+qhi))
        dl.append(m*m/(1+qhi))
        dh.append(factor*m*m/(1+qlo))
    wmin=(dl[1]-dh[0])/(dh[2]-dl[0])
    wmax=(dh[1]-dl[0])/(dl[2]-dh[0])
    possible=[k-(1-w)*p-w*B for k,p,B,w in itertools.product(
        (low[1],high[1]),(low[0],high[0]),(low[2],high[2]),(wmin,wmax))]
    return {"w_min":wmin,"w_max":wmax,"C_interval_GeV2":[min(possible),max(possible)],
            "C_abs_upper_GeV2":max(map(abs,possible)),"u_pi_upper":high[0],"d_pi_upper":dh[0]}


def main():
    q=[YUK["up"][0][0]/yd for yd in ROW["down_strange_bottom_MZ"]]
    u,d=charged_vectors(q)
    w,C=eliminate(u,d)
    env=envelope(q)
    weak=read("../1732/inputs.json")["exploratory_weak_offsets"][CFG["flavour_choice"]]
    ws,wb=ROW["down_strange_bottom_marginal_widths"][1:]
    neu=runpy.run_path(str(HERE/CFG["neutral_flow"]))
    pt=neu["endpoint"]("explicit_open",.25,explicit_offsets={**weak,"strange":ws,"bottom":wb})
    tan=pt["row"]["positive_zero_stop_tan_beta"]
    mH=CFG["charged_scalar_reference_GeV"]
    assert mH==read(CFG["prior_spectrum"])["ordinary_soft_scale_GeV"]
    uref,dref=charged_vectors(q,mH)
    A=[1+(x-y*tan*tan)/mH**2 for x,y in zip(uref,dref)]
    Bobs,*errors=CFG["b_branching"]
    Berr=math.hypot(*errors)
    cap=read(CFG["prior_yukawa_cap"])["perturbativity_Y_squared_upper"]
    cap_u=cap/(YUK["up"][2][0]*YUK["charged_leptons"][2][0])
    cap_d=cap/(ROW["down_strange_bottom_MZ"][2]*YUK["charged_leptons"][2][0])
    pi_cost=cap_u*env["u_pi_upper"]+cap_d*env["d_pi_upper"]
    floor=CFG["alternative_mass_floor_GeV"]
    abar_cap=cap_u/(floor**2-pi_cost)
    central_B=belle_reference(0)
    T0=math.sqrt(Bobs/central_B)
    cases={}
    for name in ("RM123_axial","RBC_axial"):
        R,sig=W["experimental_comparisons"][name]["central_ratio_and_marginal_error"]
        r0=R/W["weak_Vus_over_Vud"]
        down_min=1+(r0-1)/w
        signed=parameter_solution(u,d,r0,T0)
        sharp=(r0-(1-w+w*T0))/abs(C)
        assert abs(abs(signed["a_over_Api_GeV_inv2"])-sharp)<1e-9
        rows={}
        for h in (0,1,2):
            r=(R-h*sig)/W["weak_Vus_over_Vud"]
            # The extra finite-rate window stays in force even at h=0.
            H=CFG["B_over_pi_rate_matching"][0]
            T=math.sqrt((Bobs+h*Berr)/(H*belle_reference(h)))
            needed=max(0,(r-1-env["w_max"]*(T-1))/env["C_abs_upper_GeV2"])
            rows[str(h)]={"r_K_lower":r, "T_B_upper":T,"minimum_abs_ahat_uniform_GeV_inv2":needed,
                          "one_down_source_min_BR_B":H*belle_reference(h)*(1+(r-1)/env["w_max"])**2,
                          "excluded_with_common_cap_and_mass_floor":bool(needed>abar_cap)}
        cases[name]={"central_r_K":r0,"central_down_only_min_B_amplitude_over_pion":down_min,
                     "central_down_only_min_B_branching":central_B*down_min**2,
                     "central_general_equality_witness":signed,"outward_boxes":rows}
    return {"quark_ratios_u_over_d_s_b":q,"tan_beta_soft_from_joint_point":tan,"u_GeV2":u,"d_GeV2":d,
            "w":w,"C_GeV2":C,"uniform_bounds":env,"typeII_K_over_pi_rate_shift":(A[1]/A[0])**2-1,
            "typeII_B_over_pi_rate_shift":(A[2]/A[0])**2-1,"SM_B_branching_normalized_to_pi":central_B,
            "central_B_reduced_amplitude":T0, "bounded_alignment_abs_ahat_ceiling_GeV_inv2":abar_cap,
            "central_comparisons":cases}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    res=main()
    # Numerical portability for the underlying Yukawa shooting.
    def tidy(v):
        if isinstance(v,dict): return {k:tidy(x) for k,x in v.items()}
        if isinstance(v,list): return [tidy(x) for x in v]
        return float(f"{v:.11g}") if isinstance(v,float) else v
    res=tidy(res)
    if args.check:
        assert res==read("results.json")
        print("1735 verified charged scalar elimination, finite propagators, mass-ratio boxes and Belle comparison")
    else:
        (HERE/"results.json").write_text(json.dumps(res,indent=2,allow_nan=False)+"\n",encoding="utf8")
