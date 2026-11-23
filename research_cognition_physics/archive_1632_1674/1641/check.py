"""Tree-phase closure with four correlated, directly measured Bs time coefficients."""
import argparse
import cmath
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent

def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))

def degrees(z):
    return math.degrees(cmath.phase(z))

def circular(x):
    return (x+90)%180-90

def score_model(cfg):
    d=cfg["Bs_Run2"]; obs=np.array(d["central"])
    cov=sum(np.outer(d[k],d[k])*d[k+"_corr"] for k in ("stat","syst"))
    assert min(np.linalg.eigvalsh(cov))>0
    def evaluate(deg):
        t=math.radians(2*deg); c,s=math.cos(t),math.sin(t)
        A=np.array([[c,-1,s,0],[s,0,-c,-1]])
        b=A@obs
        return float(b@np.linalg.solve(A@cov@A.T,b))
    return evaluate

def profile_summary(score):
    grid=np.linspace(0,180,721)
    best=int(np.argmin([score(x) for x in grid]))
    assert 0<best<len(grid)-1
    a,b=grid[best-1],grid[best+1]
    for _ in range(75):
        l=(2*a+b)/3; r=(a+2*b)/3
        if score(l)<score(r): b=r
        else: a=l
    center=(a+b)/2; floor=score(center)
    def edge(delta,sign):
        inner=center
        for step in range(1,181):
            outer=center+sign*step/2
            if score(outer)>floor+delta:
                break
            inner=outer
        else:
            raise RuntimeError("not bracketed")
        for _ in range(55):
            mid=(outer+inner)/2
            if score(mid)>floor+delta: outer=mid
            else: inner=mid
        return (inner+outer)/2
    return {"minimum_deg":center,"minimum_score":floor,
            "profile_d1_deg":[edge(1,-1),edge(1,1)],"profile_d4_deg":[edge(4,-1),edge(4,1)]}

def ckm(s12,s23,s13,phi):
    c12,c23,c13=np.sqrt(1-np.array([s12,s23,s13])**2)
    e=cmath.exp(1j*phi)
    return np.array([[c12*c13,s12*c13,s13/e],[-s12*c23-c12*s23*s13*e,
        c12*c23-s12*s23*s13*e,s23*c13],[s12*s23-c12*c23*s13*e,
        -c12*s23-s12*c23*s13*e,c23*c13]])

def angles(K):
    ls=K[:,1].conj()*K[:,2]
    u,c=ls[:2]
    gamma_s=degrees(c/u)
    gamma_d=degrees(-K[0,0]*K[0,2].conjugate()/(K[1,0]*K[1,2].conjugate()))
    charm=degrees(-K[1,0]*K[0,1]/(K[1,1]*K[0,0]))
    assert abs(circular(gamma_s-gamma_d-charm))<2e-12
    # Insert D0 -> K+ K- as a single-tree phase reference for B+- -> DK+-.
    A_c=K[1,2]*K[0,1].conjugate(); A_u=K[0,2]*K[1,1].conjugate()
    d=K[1,1].conjugate()*K[0,1]
    assert abs(circular(gamma_s+degrees(A_u/A_c*d.conjugate()/d)))<2e-12
    return gamma_s,gamma_d,charm,ls

def main():
    cfg=read("inputs.json")
    states=read(cfg["witnesses"]); mix=read(cfg["mixing"]); traj=read(cfg["trajectory"])
    base=read(cfg["matrix_cfg"])["witness"]; pole=read(cfg["poles"])
    score=score_model(cfg); profile=profile_summary(score)
    U=ckm(*base["ordinary_rotation_s12_s23_s13_delta"])
    gs,gd,dc,ls=angles(U)
    # Uniform on every F and M in the same polar completion as 1639.
    M=cfg["mass_band_GeV"][0]; ms,mb=base["light_masses_GeV"][1:]
    rb=math.sqrt(pole["right_doublet_probability"])
    ellb=mb*rb/M; ells=ms/M
    h0=math.sqrt(1-2*pole["left_difference"]/pole["k"])
    disk=2*(2*ells*ellb+ells**2)/(math.sqrt(h0**2-ells**2)+h0)
    gamma_shift=math.degrees(sum(math.asin(disk/(h0*abs(z))) for z in ls[:2]))
    allF=[gs-gamma_shift,gs+gamma_shift]
    omega,down,up=cfg["Bs_Run1plus2_Omega_deg"]; chsig,cherr=cfg["charged_Run3_gamma_deg"]
    bands={}
    for n in (1,2):
        ph=list(map(math.degrees,traj["common_witness"][str(n)]["mix_phase_range_rad"]))
        direct=[omega-n*down,omega+n*up]
        gamma_box=[direct[0]-ph[1],direct[1]-ph[0]]
        charged=[chsig-n*cherr,chsig+n*cherr]
        joint=[max(gamma_box[0],charged[0]),min(gamma_box[1],charged[1])]
        omega_pred=[allF[0]+ph[0],allF[1]+ph[1]]
        bands[str(n)]={"measured_Omega_working":direct,"mixing_phase_working_deg":ph,
            "gamma_from_neutral_control":gamma_box,"with_charged_calibration":joint,
            "frozen_rotation_Omega":omega_pred,"fixed_rotation_gap":max(0.,joint[0]-allF[1],allF[0]-joint[1])}
    pts={}
    for k in ("common_witness","high_mass_witness"):
        w=states[k]; a=np.asarray(w["light_charged_K"]); K=a[...,0]+1j*a[...,1]
        g1,g2,g3,l=angles(K)
        # Arbitrary row/column rephasings must leave the three tree observables unchanged.
        U_test=np.diag(np.exp(1j*np.array([.3,-1.1,.4])))
        D_test=np.diag(np.exp(1j*np.array([.2,1.9,-.5])))
        test=angles(U_test@K@D_test)
        assert all(abs(circular(x-y))<5e-12 for x,y in zip(test[:3],(g1,g2,g3)))
        F=complex(*w["F_L_R_and_L_N_P"][1])
        nonclose=ms*mb/M**2*F if k=="common_witness" else ms*mb/cfg["mass_band_GeV"][1]**2*F
        assert abs(sum(l)-nonclose)<2e-13 and abs(g1-gs)<gamma_shift
        phis={}
        for tier,v in mix[k]["tiers"].items():
            pred=g1+math.degrees(v["phase"])
            phis[tier]={"Omega_phase_only_deg":pred,"Run2_profile_delta_chi2":score(pred)-profile["minimum_score"],
                        "corresponding_frequency_ps_inverse":v["frequency"]}
        pts[k]={"gamma_s":g1,"gamma_d":g2,"charm_quartet_angle_deg":g3,
                "strange_bottom_nonclosure_abs":abs(nonclose),"central_tier_phase_projections":phis}
    return {"Run2_four_raw_profile":profile,"polar_family":{"center_gamma_s":gs,
                "uniform_lambda_displacement":disk,"uniform_gamma_displacement_deg":gamma_shift,
                "allowed_gamma_fixed_rotation":allF},"workboxes":bands,"frozen_points":pts}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write-new",action="store_true"); p.add_argument("--check",action="store_true")
    opt=p.parse_args()
    result=json.loads(json.dumps(main()),parse_float=lambda x:float(f"{float(x):.10g}"))
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as fp:
            json.dump(result,fp,indent=2); fp.write("\n")
    else:
        assert result==read("results.json")
    print("1641 PASS: raw time observables, charm convention, and all-F tree-phase bound")
