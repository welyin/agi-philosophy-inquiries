"""Coherent three-meson test for the two common electroweak mass sources."""
import argparse
import hashlib
import itertools
import json
import math
import runpy
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())
CFG=read("inputs.json")
OLD=runpy.run_path(str(HERE/CFG["previous_vertices"]))
KA=runpy.run_path(str(HERE/CFG["kaon_matching"]))
KC=read(CFG["kaon_inputs"])


def pack(x):
    return [float(x.real),float(x.imag)]


def sources():
    old=runpy.run_path(str(OLD["HERE"]/OLD["CFG"]["charged_background"]))
    flow,rep=old["charged_background"](.1)
    unpack=old["unpack"]
    d=unpack(rep["down_texture_10PeV"])
    e=unpack(rep["charged_texture_10PeV"]).T
    zm,zt=e[0,1]/d[0,1],e[2,2]/d[2,2]
    elow=unpack(rep["charged_left_MZ"])
    y=flow["shoot"](rep["gamma_MZ_deg"],flow["angles"]("direct"),elow,.1)[0]
    yd,ye,yu=(y[i,:,:3] for i in (1,2,0))
    left=OLD["canonical_left"](d,yd,yu)
    dm,em=(flow["central"](k) for k in ("down","charged_leptons"))
    V=flow["low_ckm"](rep["gamma_MZ_deg"],flow["angles"]("direct"))
    result={}
    for n,z in {"canonical":1.,"first_edge":zm,"third_edge":zt,
                "A_vanishes_first_edge":-zm/3,"A_vanishes_third_edge":-zt/3}.items():
        b=(d-e/z)/4
        D=OLD["transported"](d,b,yd,left,dm)[0]/math.sqrt(2)
        E=OLD["transported"](e.T,(-3*z*b).T,ye,np.eye(3),em)[0]/math.sqrt(2)
        result[n]=(z,D,E)
    return result,V


def short(mh,alpha,had):
    cfg=dict(KC,mh_GeV=mh)
    wil,lr,ss=KA["contractions"](cfg,KA["scale_alpha"](cfg,alpha),had)
    running=(KA["scale_alpha"](cfg,alpha)/alpha)**(12/23)
    return float(lr),float(ss),running


def vertex_window():
    w=CFG["matrix_element_relative_half_window"]
    return itertools.product(CFG["neutral_poles_window_GeV"],CFG["QCD_alpha_s_MZ_interval"],
                             itertools.product([1-w,1+w],repeat=4))


def meson(j,D,V,eta):
    data=CFG["mesons"][j]
    i,k=data["indices"]
    phase=(V[0,i].conjugate()*V[0,k])
    phase=phase.conjugate()/abs(phase)
    a,b=D[i,k]*phase,D[k,i].conjugate()*phase
    if abs(a)*abs(b)==0:
        raise ValueError("Missing chiral partner")
    ab=abs(a*b)
    s=.5*(abs(a)**2+abs(b)**2)
    unit=(a*b)/ab
    had=np.array(data["matrix_elements"])
    ref=CFG["central_real_doublet"]
    lr,ss,qcd=short(ref["mh_GeV"],ref["alpha_s_MZ"],had)
    den=ref["mh_GeV"]**2
    # Return M12 rather than the conjugate appearing at the matching end.
    amp=-(qcd**2/den)*(lr*a*b*abs(eta)**2+.5*ss*((a*eta)**2+(b*eta.conjugate())**2)).conjugate()
    bound=(lr*ab-abs(ss)*s)*qcd**2*abs(eta)**2/den
    # The coefficient times common transfer^2 is enclosed independently of alpha for the running Yukawas.
    vals,real_vals=[],[]
    for m,al,moves in vertex_window():
        l,s1,_=short(m,al,had*np.array(moves))
        vals.append(l*ab-abs(s1)*s)
        real_vals.append(-l*(a*b).real-.5*abs(s1)*abs(a*a+b.conjugate()**2))
    qlo,qhi=sorted(short(125,al,had)[2] for al in CFG["QCD_alpha_s_MZ_interval"])
    scales=CFG["additional_common_quark_transfer"]
    conservative=min(vals)*((qlo*scales[0])**2 if min(vals)>=0 else (qhi*scales[1])**2)
    real_floor=min(real_vals)*((qlo*scales[0])**2 if min(real_vals)>=0 else (qhi*scales[1])**2)
    # An aligned two-doublet model has two perpendicular neutral fields with couplings exp(i phi),i exp(i phi).
    # The h pole has no FCNC here. Tune the two nonstandard masses and their mixing angle instead.
    phi=-.5*np.angle(a/b)
    def split(t):
        left,s1,_=short(250.,ref["alpha_s_MZ"],had)
        right,s2,_=short(250.*t,ref["alpha_s_MZ"],had)
        return ab*(left+right/t**2)+s*(s1-s2/t**2)
    cancel=None
    if split(1)*split(4)<0:
        lo,hi=1.,4.
        for _ in range(62):
            mid=(lo+hi)/2
            if split(mid)>0:
                lo=mid
            else:
                hi=mid
        residual=split(mid)/(ab*lr)
        assert abs(residual)<1e-11
        cancel={"masses_GeV":[250.,250.*mid],"source_axis_phase_deg":np.degrees(phi),
                "zero_residual_relative_to_LR":residual}
    # T=sum |eta|^2/m^2. Bound is a projection on the fixed LR phase, independent of every scalar phase.
    return {"right_left_MZ": [pack(a),pack(b)],"common_LR_direction_M12":pack(-unit.conjugate()),
            "rho_chiral_large_over_small":max(abs(a/b),abs(b/a)),
            "central_LR_S_LL_GeV3":[lr,ss],"single_h_M12_GeV":pack(amp),
            "single_h_abs_GeV":abs(amp),"central_unknown_h_phase_floor_GeV":max(0.,bound),
            "all_poles_projection_coeff_floor_GeV3":conservative,
            "all_poles_real_coeff_floor_GeV3":real_floor,
            "aligned_H_A_individual_cancellation":cancel,
            "all_poles_floor_from_h_alone_GeV":max(0,conservative)*abs(eta)**2/den}


def run():
    source,V=sources()
    obs=KA["oscillation"](KC)
    caps={"K":obs["total_abs_upper"]*1e-15}
    bd=read(CFG["B_clocks"])["data"]
    hbar=6.582119569e-13     # GeV ps
    # Bd,s quoted as dispersive norms only in the documented |Gamma12|<<|M12| approximation.
    caps.update({j:(bd[k][0]+2*bd[k][1])*hbar/2 for j,k in (("Bd","delta_md"),("Bs","delta_ms"))})
    examples={}
    eta=1/CFG["central_real_doublet"]["vB_over_v"]
    for n,(z,d,e) in source.items():
        rows={j:meson(j,d,V,eta) for j in CFG["mesons"]}
        for k,r in rows.items():
            r["single_h_abs_over_observed_cap"]=r["single_h_abs_GeV"]/caps[k]
            r["window_floor_over_observed_cap"]=r["all_poles_floor_from_h_alone_GeV"]/caps[k]
            r["minimum_external_amplitude_after_allowing_full_observed_cap_GeV"]=max(0,r["all_poles_floor_from_h_alone_GeV"]-caps[k])
        c=rows["K"]["all_poles_real_coeff_floor_GeV3"]
        assert c>0
        rbc=KC["RBC_preliminary_Delta_m_in_1e15GeV"]
        rest=(rbc[0]+2*sum(rbc[1:]))*.5e-15
        cap=caps["K"]+rest
        s=CFG["central_real_doublet"]["vB_over_v"]
        spec=2/(s*s*(1-s*s))
        hmass=CFG["central_real_doublet"]["mh_GeV"]
        T=cap/c
        emax=hmass*math.sqrt(T)
        props={"max_inverse_mass_moment_GeV_minus2":T,
               "min_extra_real_compensation_after_K_and_SM_GeV":max(0,c*(eta/hmass)**2-cap),
               "max_abs_eta_h_if_only_SM_baseline":emax,
               "min_abs_kappa_V_two_doublets_any_vevs":math.sqrt(max(0,1-emax**2/4)),
               "min_abs_kappa_V_two_doublets_registered_vev":math.sqrt(max(0,1-s**2*(1-s**2)*emax**2)),
               "two_doublets_all_vevs_min_largest_mass_GeV":math.sqrt(8/T),
               "at_registered_vev_formal_min_largest_mass_GeV":math.sqrt(spec/T),
               "all_registered_vev_coupled_masses_can_be_in_certified_window_by_this_bound":
               math.sqrt(spec/T)<=CFG["neutral_poles_window_GeV"][1]}
        examples[n]={"z":pack(z),"mesons":rows,"kaon_real_axis_independent_SM_reference":props}
    return {"observed_norm_working_caps_GeV":caps, "examples":examples}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=OLD["clean"](run())
    data=json.dumps(result,indent=2)+"\n"
    path=HERE/"results.json"
    if args.check:
        assert result==read("results.json")
    else:
        path.write_text(data,encoding="utf-8")
    print("1795",hashlib.sha256(data.encode()).hexdigest()[:12])
