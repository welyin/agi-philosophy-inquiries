"""Run the scalar Delta F=2 matching, and compare neutral kaons with muonium."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def read(n):
    return json.loads((HERE/n).read_bytes())
CFG=read("inputs.json")
M=runpy.run_path(str(HERE/CFG["matter"]))
EPS=runpy.run_path(str(HERE/CFG["weak_CP"]))


def pack(z):
    return [float(z.real),float(z.imag)]


def bounds():
    gamma=np.array([[-10.,1/6],[-40.,34/3]])
    mat=gamma.T+16*np.eye(2)
    low,high=sorted(np.linalg.eigvals(mat))
    gap=high-low
    assert 0<6<low<high and low*(high-6)>high*(low-6)
    w=CFG["tolerances"]["hadronic_fraction"]
    had=M["CFG"]["mesons"]["K"]["matrix_elements"]
    s_max=max(abs(had[2]),abs(had[3])/(6*gap))*(1+w)
    l_min=(1-w)*had[1]
    initial=np.array([1.,0.])
    for tau in (0.,.013,.1,.4,2.):
        x=(high-6)*math.exp(-low*tau)/gap+(6-low)*math.exp(-high*tau)/gap
        y=(math.exp(-high*tau)-math.exp(-low*tau))/(6*gap)
        modes,V=np.linalg.eig(mat)
        np.testing.assert_allclose(V @ (np.exp(-modes*tau)*np.linalg.solve(V,initial)),[x,y],atol=1e-14)
        assert 0<=x<=1+1e-14 and -1./(6*gap)<=y<=0
    return {"eigenvalues_RR":[low,high],"max_scalar_tensor_weight":1/(6*gap),
            "all_scale_S_upper_GeV3":s_max,"all_scale_L_lower_GeV3":l_min}


def clock(V):
    inputs={"Vus":abs(V[0,1]),"Vub":abs(V[0,2]),"Vcb":abs(V[1,2])}
    gamma=np.degrees(np.angle(-V[0,0]*V[0,2].conj()/(V[1,0]*V[1,2].conj())))
    w=EPS["amplitudes"](inputs,gamma)[0]
    return w@EPS["COEFFICIENTS"],math.sqrt(w@EPS["COV"]@w),w@EPS["ERRORS"]


def run():
    b=bounds()
    variants,V=M["sources"]()
    kc=M["KC"]
    yrun=M["short"](125.,0.1181,[0,1,0,0])[2]
    y_min=min(M["short"](125.,a,[0,1,0,0])[2] for a in M["CFG"]["QCD_alpha_s_MZ_interval"])*.95
    ref=M["KA"]["oscillation"](kc)
    ckm_eps,ckm_width,boxwidth=clock(V)
    obs_c=ref["total_abs_upper"]*1e-15
    phase=(V[0,1].conjugate()*V[0,0])
    s=CFG["vB_over_v"]
    c=math.sqrt(1-s*s)
    weak=runpy.run_path(str(HERE/CFG["old_thresholds"]))
    v=CFG["v_GeV"]
    m_e=M["OLD"]["CFG"]["lepton_masses_GeV"]
    read_ew=read("../1386/inputs.json")["yukawa"]["charged_leptons"]
    to_pole=np.sqrt(np.array(m_e)/(np.array(read_ew)[:,0]*v/math.sqrt(2)))
    result={}
    for label in ("A_vanishes_first_edge","A_vanishes_third_edge"):
        z,d,e=variants[label]
        e=to_pole[:,None]*e*to_pole[None,:]
        a,bb=d[1,0]*phase.conjugate()/abs(phase), d[0,1].conjugate()*phase.conjugate()/abs(phase)
        # K phase convention is inherited from the physical tree. Both muonium charges need no meson phase.
        ka=-2*.18*yrun**2*(a*bb).conjugate()/(s*c)**2
        g3=abs(2*e[1,0]*e[0,1].conjugate()/(s*c)**2)*math.sqrt(2)/8
        rho=abs(a/bb)
        coeff=y_min**2*(-b["all_scale_L_lower_GeV3"]*(a*bb).real-
                         b["all_scale_S_upper_GeV3"]*abs(a*a+bb.conjugate()**2)/2)
        assert coeff>0 and 0.95<rho<1.05
        rows=[]
        numerator=ka.imag*1e15*ref["central_eps_response_per_im"]
        extra=ckm_eps+2*boxwidth-(EPS["EPS_EXP"]-2*EPS["EPS_ERROR"])
        assert extra>0 and numerator<0
        for mass in CFG["witness_masses_GeV"]:
            amp=ka/mass**2
            rel=g3/(mass**2*CFG["GF_GeV_minus2"])
            eps_shift=amp.imag*1e15*ref["central_eps_response_per_im"]
            # A common normalization and CKM model are used; this checks compatibility rather than assigning a probability.
            required_msm=ref["Delta_m_central"]*.5e-15-amp.real
            pull=(ckm_eps+eps_shift-EPS["EPS_EXP"])/math.hypot(ckm_width,EPS["EPS_ERROR"])
            bottom=EPS["EPS_EXP"]-2*EPS["EPS_ERROR"]
            conversion=1e15*ref["central_eps_response_per_im"]
            need_im=(bottom-ckm_eps-2*boxwidth)/conversion
            direction=(math.asin(np.clip(need_im/abs(amp),-1,1))-np.angle(amp))
            rows.append({"mass_GeV":mass,"K_M12_GeV":pack(amp),
                         "K_fraction_of_experimental_abs_cap":abs(amp)/obs_c,
                         "SM_Re_M12_needed_GeV":required_msm,
                         "SM_epsilon_plus_new":ckm_eps+eps_shift,
                         "epsilon_difference_in_marginal_widths":pull,
                         "epsilon_K_at_common_plus_two_individual_theory_widths":ckm_eps+eps_shift+2*boxwidth,
                         "fits_two_individual_width_rectangles":bool(ckm_eps+eps_shift+2*boxwidth>=bottom),
                         "minimum_M12_rotation_deg_to_reach_two_width_boundary":max(0,np.degrees(direction)),
                         "muonium_abs_G3_over_GF":rel,"muonium_P_from_rounded_PSI_G3_translation":
                          CFG["PSI_0_1_T_conversion_limit"]*(rel/CFG["PSI_only_G3_limit_over_GF"])**2,
                         "delta_inverse_alpha_123_relative_to_old_30TeV_second_doublet":
                         (-np.array([.1,1/6,0])*math.log(weak["soft"]/mass)/(2*math.pi)).tolist()})
        result[label]={"z":pack(z),"K_uniform_LL_floor_GeV3":coeff,
                       "positive_floor_survives_additive_kernel_defect_below_GeV3_after_dividing_by_abs_ab":
                       coeff/(abs(a*bb)*y_min**2),
                       "K_constant_M2_times_M12_GeV3":pack(ka),
                       "inclusive_Mu_G3_over_K_absM12_GeV_minus3":g3/abs(ka),
                       "equal_mass_min_GeV_from_epsilon_with_two_marginal_box":math.sqrt(-numerator/extra),
                       "max_admissible_inverse_mass_T_if_only_SM_real_budget":(obs_c+5.8e-15)/coeff,
                       "cases":rows}
    lo,hi=list(result.values())
    ratio=hi["inclusive_Mu_G3_over_K_absM12_GeV_minus3"]/lo["inclusive_Mu_G3_over_K_absM12_GeV_minus3"]
    assert abs(ratio-(abs(complex(*hi["z"]))/abs(complex(*lo["z"])))**4)<1e-8
    return {"strong_LL":b,"original_SM_epsilon_and_width":[ckm_eps,ckm_width],
            "baseline_epsilon_distance_in_joint_quadrature":(EPS["EPS_EXP"]-ckm_eps)/math.hypot(ckm_width,EPS["EPS_ERROR"]),
            "epsilon_SM_max_plus_two_individual_widths":ckm_eps+2*boxwidth,
            "light_quark_mass_running_MZ_to_163GeV":yrun,
            "previous_colored_thresholds_below_Msoft":[[float(m),float(k[-1])] for m,k in weak["thresholds"]
                                                   if m<weak["soft"] and k[-1]!=0],
            "common_potential":{"Z1":(CFG["mh_GeV"]/v)**2,"Z2":0.5,
                                "Z3_Z4_Z5_Z6_Z7_Y3":0,"M2_geV2":"positive free invariant"},
            "ratio_of_Mu_amplitudes_at_equal_scalar_K":ratio,"ratio_of_Mu_probabilities_at_equal_scalar_K":ratio**2,
            "cases":result}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=M["OLD"]["clean"](run())
    data=json.dumps(result,indent=2)+"\n"
    if args.check:
        assert result==read("results.json")
    else:
        (HERE/"results.json").write_text(data,encoding="utf-8")
    print("1796",hashlib.sha256(data.encode()).hexdigest()[:12])
