"""One exotic doublet, the physical charged-current completion and recoil flavor."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

D=Path(__file__).resolve().parent
def read(name):
    return json.loads((D/name).read_bytes())

def enc(a):
    a=np.asarray(a)
    return np.stack((a.real,a.imag),axis=-1).tolist()

def polar_columns(col):
    proj=np.eye(len(col))-col@col.conj().T
    vals,vecs=np.linalg.eigh(proj[:3,:3])
    assert min(vals)>0
    return proj[:,:3]@(vecs@np.diag(vals**-.5)@vecs.conj().T)

def ckm(s12,s23,s13,delta):
    c12,c23,c13=np.sqrt(1-np.array([s12,s23,s13])**2)
    ph=np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13/ph],
        [-s12*c23-c12*s23*s13*ph,c12*c23-s12*s23*s13*ph,s23*c13],
        [s12*s23-c12*c23*s13*ph,-c12*s23-s12*c23*s13*ph,c23*c13]])

def width_factor(m,M,MW):
    x,w=(m/M)**2,(MW/M)**2
    if m+MW>=M:
        return 0.
    return math.sqrt((1-x-w)**2-4*x*w)*((1+x)*((1-x)**2+w*(1+x)-2*w*w)-12*w*x)

def build(cfg, high=False):
    w = dict(cfg["witness"])
    if high:
        w.update(cfg["high_mass"])
    M=w["M_D_GeV"]; m5=w["m5_physical_GeV"]
    m1,ms,mb=w["light_masses_GeV"]
    fr,ph=w["F_R_absolute_and_phase"]; F=fr*np.exp(1j*ph)
    p=read(cfg["Z_projection"])
    rb=math.sqrt(p["right_doublet_probability"]); rs=-F/rb
    sb=math.sqrt(2*(p["left_difference"]/p["k"]-(mb*rb/M)**2))
    available=1-rb**2-abs(rs)**2
    def trial(r5):
        r4=math.sqrt(available-r5**2)
        m4=math.sqrt((M*M-ms**2*abs(rs)**2-mb**2*rb**2-m5**2*r5**2)/r4**2)
        masses=np.asarray([m1,ms,mb,m4,m5])
        r=np.array([0,rs,rb,r4,r5])
        ell=masses*r/M
        bdot=ell[2].real*sb
        s4=-bdot*(m5*m5-mb*mb)/((m5*m5-m4*m4)*ell[3].real)
        s5= bdot*(m4*m4-mb*mb)/((m5*m5-m4*m4)*ell[4].real)
        return masses,r,ell,np.array([0,0,sb,s4,s5])
    lo,hi=1e-10,.01
    assert np.linalg.norm(trial(lo)[3])>1>np.linalg.norm(trial(hi)[3])
    for _ in range(80):
        mid=(lo+hi)/2
        if np.linalg.norm(trial(mid)[3])>1: lo=mid
        else: hi=mid
    m,r,l,s=trial((lo+hi)/2)
    t=m*s/np.linalg.norm(m*s)
    u=ckm(*w["ordinary_rotation_s12_s23_s13_delta"])
    left=np.column_stack((s,l)); right=np.column_stack((t,r))
    WL=np.column_stack((polar_columns(left)@u.conj().T,left))
    WR=np.column_stack((polar_columns(right),right))
    K=WL[:,:3].conj().T
    mass=WL.conj().T@np.diag(m)@WR
    np.testing.assert_allclose(WL.conj().T@WL,np.eye(5),atol=1e-12)
    np.testing.assert_allclose(WR.conj().T@WR,np.eye(5),atol=1e-12)
    assert max(abs(mass[:4,4]))<1e-9 and max(abs(mass[3,:3]))<1e-9 and abs(mass[4,3])<1e-9
    assert abs(mass[4,4]-M)<1e-9 and abs(np.linalg.norm(m*s)-mass[3,3])<1e-9
    T=np.outer(r,r.conjugate()); PL=np.outer(l,l.conjugate()); PS=np.outer(s,s.conjugate())
    left_F=(-PS-2*PL)[1,2]; right_F=-T[1,2]
    assert abs(F-right_F)<1e-13 and abs(left_F-2*ms*mb/M**2*F)<1e-13
    np.testing.assert_allclose(K.conj().T@K,np.eye(5)-PL-PS,atol=1e-12)
    Vprod=K[:,1].conjugate()*K[:,2]
    assert abs(sum(Vprod)-ms*mb/M**2*F)<1e-13
    z=read(cfg["charged_caps"]); ref=read("../1464/inputs.json")["HPQCD_tableI_eq15"]
    c=ref["Xt_ref"]/(math.pi*ref["sine_squared_ref"]*ref["alpha_inverse_ref"])
    k=read(cfg["rare_calibration"])["sm_coefficient_ratio"]["center"]
    L=-c*Vprod[-1]
    N=L+left_F+F; P=k*L+left_F-F
    wide=read(cfg["rare_amplitudes"]); Hnu=wide["Hnu"]; Hmu=wide["Hmu_heavy"]
    assert abs(K[0,2])<z["charged_current_caps"]["abs_Kub"] and abs(K[1,2])<z["charged_current_caps"]["abs_Kcb"]
    # The scalar from the same Higgs is spin-orthogonal and has nonnegative accepted yield.
    mBs=read("../1404/inputs.json")["const_GeV"]["M_Bs"]
    scalar_max=(mBs/w["h_mass_GeV"])**2*abs(F)
    obs={}
    for name,weight in [("CMS",wide["calculated_light_response"]),("LHCb_R3",wide["run3_central_light_response"])]:
        relative=P*np.conjugate(L)/abs(L)
        mu0=1e9*Hmu*(relative.real**2+weight*relative.imag**2)
        obs[name]=[mu0,mu0+1e9*Hmu*scalar_max**2]
    h=read(cfg["rare_calibration"])
    strength=Hnu*abs(N)**2/(1e-6*h["BelleII"]["reference_SD_in_1e-6"])
    cms=h["CMS"]; r3=read("../1404/inputs.json")["LHCb_run3_2024_25"]["s_in_1e9"]
    assert cms["branching_in_1e-9"]-sum(cms["lower_widths"])<min(obs["CMS"])<=max(obs["CMS"])<cms["branching_in_1e-9"]+sum(cms["upper_widths"])
    assert r3[0]-r3[1]-r3[3]<min(obs["LHCb_R3"])<=max(obs["LHCb_R3"])<r3[0]+r3[2]+r3[4]
    assert abs(strength-h["BelleII"]["mu"])<h["BelleII"]["width"]
    mw=w["W_mass_GeV"]
    weights=np.abs(r)**2*np.array([width_factor(mi,M,mw) for mi in m])
    g=w["weak_g"]
    widths=g**2*M**3/(64*math.pi*mw**2)*weights
    ratio=widths[1]/widths[2]
    hard=np.zeros_like(mass)
    hard[3,3]=mass[3,3]; hard[4,4]=M
    yukawa=np.sqrt(2)*(mass-hard)/w["Higgs_vev_GeV"]
    largest_yukawa=np.linalg.svd(yukawa,compute_uv=False)[0]
    assert largest_yukawa < np.sqrt(4*math.pi) and all(m[3:]>M)
    return {"right_D":enc(r),"left_D":enc(l),"left_S":enc(s),"right_S":enc(t),
            "light_charged_K":enc(K[:,:3]),"weak_basis_mass_GeV":enc(mass),"physical_masses_GeV":m.tolist(),
            "singlet_hard_mass_GeV":float(mass[3,3].real),"F_L_R_and_L_N_P":enc([left_F,right_F,L,N,P]),
            "abs_Vub_Vcb":[abs(K[0,2]),abs(K[1,2])],"Tbb":rb**2,"Tss":abs(rs)**2,
            "LEP_D_L_R":[p["k"]*(sb**2/2+abs(l[2])**2),p["k"]*rb**2/2],
            "kappa_b":1-sb**2-rb**2,"M_Y_GeV": M,
            "Y_to_Wq_widths_GeV":widths.tolist(),"bW_branching_tree":widths[2]/sum(widths),
            "width_over_mass_tree":sum(widths)/M,"largest_yukawa_singular_value":largest_yukawa,
            "Belle_short_distance_mu":strength,"fitted_muons_in_1e9":obs,"width_ratio_sW_to_bW":ratio}

def boxes(cfg):
    w=cfg["witness"]; M=w["M_D_GeV"]; ms,mb=w["light_masses_GeV"][1:]
    eta=ms*mb/M**2; shape=width_factor(ms,M,w["W_mass_GeV"])/width_factor(mb,M,w["W_mass_GeV"])
    pole=read(cfg["Z_pole"]); mean=np.array(pole["central"]); errs=np.array(pole["error"])
    cov=errs[:,None]*np.array(pole["correlation"])*errs
    a=np.array([0,0,-2/3,2]); den=np.array([-2,0,2,0])
    def rangeT(n):
        # Roots of the tangency equation for a ratio of linear functionals on an ellipsoid.
        x=float(den@mean); y=float(a@mean)
        assert x>n*np.sqrt(den@cov@den)
        return sorted(np.roots([x*x-n*n*(den@cov@den),-2*x*y+2*n*n*(a@cov@den),y*y-n*n*(a@cov@a)]).tolist())
    caps=read(cfg["charged_caps"]); c=caps["top_coefficient_cap"]; C=math.hypot(*caps["charged_current_caps"].values())
    k=read(cfg["rare_calibration"])["sm_coefficient_ratio"]["range"][1]
    res={}
    for name in ("LHCb_R3_1","LHCb_R3_2","CMS_1","CMS_2"):
        n=int(name[-1]); t0,t1=rangeT(n)
        row=read(cfg["rare_amplitudes"])["working_boxes"][name]
        pmax=row["abs_P_max"]; nl=row["abs_N"][0]
        fmin=max(0,(nl-c*C)/(1+(2+c)*eta))
        fmax=(pmax+k*c*C)/(1-(2+k*c)*eta)
        res[name]={"Tbb":[t0,t1],"F_abs_bounds":[fmin,fmax],
                   "width_ratio_bounds":[shape*(fmin/t1)**2,shape*(fmax/t0)**2]}
    return res

def main():
    cfg=read("inputs.json")
    return {"rate_envelopes":boxes(cfg),"common_witness":build(cfg),"high_mass_witness":build(cfg,high=True)}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-results",action="store_true",
                        help="regenerate this round's own numerical results")
    args=parser.parse_args()
    val=main()
    def stable(x):
        if isinstance(x,list): return [stable(y) for y in x]
        if isinstance(x,dict): return {k:stable(v) for k,v in x.items()}
        if isinstance(x,float): return 0.0 if abs(x)<1e-12 else float(f"{x:.12g}")
        return x
    val=stable(val)
    if args.check:
        assert val==read("results.json")
    else:
        mode="w" if args.write_results else "x"
        with (D/"results.json").open(mode,encoding="utf-8") as out: out.write(json.dumps(val,indent=2)+"\n")
    print("1635 PASS: five-state texture, shared FCNC, calibrated rates and exotic flavor bounds")
