"""Two physical radial insertions in coherent neutral-kaon mixing, in units 1e-15 GeV."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def polar_mass(cfg):
    obj=read(cfg["texture"])["reconstruction"]
    a,b,c,e,f=obj["normalized_real_edges_A_B_C_E_F"]
    r,v=obj["two_real_coefficients"]
    d=v*cmath.exp(1j*math.radians(obj["eta_example_deg"]))
    mass=np.array([[0,a,0],[b,c,r+d],[0,e,f]],complex)
    ul,s,vrh=np.linalg.svd(mass)
    ix=np.argsort(s)
    ul,vr,s=ul[:,ix],vrh.conj().T[:,ix],s[ix]
    assert np.max(abs(ul.conj().T@mass@vr-np.diag(s)))<1e-13
    unit=np.outer(ul[1,:].conj(),vr[2,:])
    lam=ul[0,1].conjugate()*ul[0,0]
    uv=lam.conjugate()/abs(lam)
    # Remove the phase of the actual s -> d tree current before squaring.
    R,L=unit[1,0]*uv,unit[0,1].conjugate()*uv
    J=(ul[0,0]*ul[1,1]*ul[0,1].conjugate()*ul[1,0].conjugate()).imag
    assert abs(J-obj["CP_quartet"])<1e-14
    table=read(cfg["radial"])["example_slopes"]
    assert table["phase_path_dim6"]==[0,0,0,0,0,0,2]
    assert table["real_path_dim6"]==[0,0,0,0,0,2,0]
    return r,d,R,L,J


def scale_alpha(cfg,alpha):
    return alpha/(1+alpha*23/(6*math.pi)*math.log(cfg["mu_GeV"]/cfg["MZ_GeV"]))


def contractions(cfg,alpha,had):
    a=alpha/(4*math.pi)
    log=2*math.log(cfg["mh_GeV"]/cfg["mu_GeV"])
    wil=[-1.5*a,1-a,1+a*(4.5-3*log),a*(1/8-log/12)]
    LR=wil[0]*had[0]+wil[1]*had[1]
    S=wil[2]*had[2]+wil[3]*had[3]
    return wil,LR,S


def pieces(w,gain,R,L,LR,S,mh):
    right=gain*R*w
    left=gain*L*w.conjugate()
    unit=-1e15/mh**2
    terms=unit*np.array([S*left**2/2,S*right**2/2,LR*left*right])
    return right,left,terms


def oscillation(cfg):
    ka=read(cfg["kaon"]); hb=ka["hbar_GeV_s"]
    dm,de=ka["ktev"]["delta_m_1e6_per_second"];ts,ds=ka["ktev"]["tau_s_ps"]
    tl,*dl=ka["kloe_tau_L_ns"]
    # M12 reported in 1e-15 GeV.
    d0=dm*1e21*hb
    short=1e27*hb/ts
    long=1e24*hb/tl
    dhi=(dm+2*de)*1e21*hb
    sh=1e27*hb/(ts-2*ds)
    lh=1e24*hb/(tl-2*sum(dl))
    # Central superweak phase and small-CP signed response at held measured dm.
    phi=math.atan2(2*d0,short-long)
    dmin=(dm-2*de)*1e21*hb
    sphi_min=math.sin(math.atan2(2*dmin,sh-1e24*hb/(tl+2*sum(dl))))
    eps,err=cfg["KLOE_epsilon"]
    return {"Delta_m_central":d0,"total_abs_upper":math.sqrt(dhi**2+sh*lh)/2,
            "phiSW_degrees":math.degrees(phi),"central_eps_response_per_im":math.sin(phi)/d0,
            "combined_im_xi_observed_abs_upper":(eps+2*err)*dhi/sphi_min}


def cp_parts(x):
    return [float(x.real),float(x.imag)]


def results():
    cfg=read("inputs.json"); r,d,R,L,J=polar_mass(cfg)
    spectrum=read(cfg["matching_inputs"])["yukawa"]["down"][-1]
    yb,dy=spectrum; mh=cfg["mh_GeV"]
    a0,amin,amax=cfg["alpha_s_MZ"]
    def run(a):
        return (scale_alpha(cfg,a)/a)**(12/23)
    had=cfg["table2_kaon_GeV3"]
    wil,lr,s=contractions(cfg,scale_alpha(cfg,a0),had)
    central_had=[lr,s]
    g=yb/math.sqrt(2)*run(a0)
    obs=oscillation(cfg)
    data={"phase_path":2*d,"real_path":2*r,"both_dim6":2*(d+r),"common_linear":0j}
    central={}
    for label,w in data.items():
        right,left,t=pieces(complex(w),g,R,L,lr,s,mh)
        star=sum(t)
        central[label]={"w":cp_parts(w),"physical_R_L":[cp_parts(right),cp_parts(left)],
                        "M12star_SLL_SRR_LR":[cp_parts(z) for z in t],"M12":cp_parts(star.conjugate()),
                        "M12_abs":abs(star),"delta_eps_signed_at_observed_dm":-star.imag*obs["central_eps_response_per_im"]}
    # Hadronic matrix elements at a common mu. Relax the coupling and Wilson running
    # independently. Each affine endpoint below contains the correlated theory choices.
    arange=[scale_alpha(cfg,x) for x in (amin,amax)]
    runs=[run(x) for x in (amin,amax)]
    vlo,vhi=cfg["relative_vertex_transfer"]
    gains=[(yb-2*dy)/math.sqrt(2)*min(runs)*vlo,(yb+2*dy)/math.sqrt(2)*max(runs)*vhi]
    width=cfg["relative_ME_halfwidth"]
    samples=[]
    for variation,alpha,gain in itertools.product(itertools.product([1-width,1+width],repeat=4),arange,gains):
        _,lr,s=contractions(cfg,alpha,np.array(had)*variation)
        def amplitude(w):
            return sum(pieces(complex(w),gain,R,L,lr,s,mh)[2]).conjugate()
        slopes={key:amplitude(w) for key,w in data.items()}
        # Re M12 defines a real symmetric two-dimensional quadratic form in (Re w, Im w).
        xx=amplitude(1).real; yy=amplitude(1j).real
        xy=(amplitude(1+1j).real-xx-yy)/2
        matrix=np.array([[xx,xy],[xy,yy]])
        eig=np.linalg.eigvalsh(matrix)[0]
        assert eig>0
        w0=np.array([2*r,0.]); direction=np.array([(d-r).real,(d-r).imag])
        minimum_x=-float(direction@matrix@w0)/float(direction@matrix@direction)
        minimum=amplitude(2*r+minimum_x*(d-r)).real
        assert abs(minimum-(w0@matrix@w0-(direction@matrix@w0)**2/(direction@matrix@direction)))<1e-12
        samples.append((slopes,float(eig),minimum,minimum_x))
    def minmax(seq):
        vals=list(seq)
        return [float(min(vals)),float(max(vals))]
    bound=obs["total_abs_upper"]
    rbc=cfg["RBC_preliminary_Delta_m_in_1e15GeV"]
    rbc_cap=(rbc[0]+2*sum(rbc[1:]))/2
    boxes={}
    for label in data:
        re=minmax(row[0][label].real for row in samples)
        im=minmax(row[0][label].imag for row in samples)
        abs_im_floor=max(0.,im[0],-im[1])
        boxes[label]={"Re_M12":re,"Im_M12":im,
                      "extra_Re_needed_if_others_centered_at_zero_at_least":max(0.,re[0]-bound),
                      "extra_Re_needed_even_with_abs_sm_rbc_cap_at_least":max(0.,re[0]-bound-rbc_cap),
                      "other_Im_plus_dm_xi_needed_at_least":max(0.,abs_im_floor-obs["combined_im_xi_observed_abs_upper"])}
    gap=min(row[0]["phase_path"].real-row[0]["real_path"].real for row in samples)
    # Stronger quantity if two theories are assigned separate short-distance nuisances.
    free_gap=boxes["phase_path"]["Re_M12"][0]-boxes["real_path"]["Re_M12"][1]
    eig=min(row[1] for row in samples); fixed=min(row[2] for row in samples)
    wmin=2*abs((r*d).imag)/abs(d-r)
    common={"positive_quadratic_lower_per_w2":eig,"min_abs_w_on_nsum2":wmin,
            "coarse_Re_bound_on_nsum2":eig*wmin**2,
            "sharp_box_Re_bound_on_nsum2":fixed,
            "minimum_real_part_minimizer_nD_range":minmax(row[3] for row in samples),
            "abs_other_Re_needed_all_nsum2":fixed-bound,"other_Re_needed_beyond_abs_sm_rbc_cap_all_nsum2":max(0.,fixed-bound-rbc_cap),
            "if_no_other_Re_max_abs_nsum":2*math.sqrt(bound/fixed),
            "if_abs_sm_rbc_cap_max_abs_nsum":2*math.sqrt((bound+rbc_cap)/fixed),
            "if_no_other_Re_max_abs_w":math.sqrt(bound/eig)}
    # Changing the hadronic input independently across the two hypothetical worlds
    # can erase the real-projection proof. The imaginary CP readout needs the
    # additional common DeltaS=1 decay-axis condition.
    cp_gap=boxes["real_path"]["Im_M12"][0]-boxes["phase_path"]["Im_M12"][1]
    assert gap>2*bound and cp_gap>2*obs["combined_im_xi_observed_abs_upper"]
    assert fixed>bound and boxes["common_linear"]["Re_M12"]==[0.,0.]
    return {"unit":"M12 in 1e-15 GeV; h couplings dimensionless","matching":{"R_shape":cp_parts(R),"L_shape":cp_parts(L),
            "J":J,"eta_QCD_vertex":run(a0),"alpha_s_mt":scale_alpha(cfg,a0),"Wilson_LR1_LR2_S1_S2":wil,
            "LR_S_GeV3": central_had,
            "physical_gain_center_bounds":[g,*gains]},"kaon_data":obs,"RBC_abs_Re_proxy_cap":rbc_cap,
            "central":central,"boxes":boxes,"separations":{"same_nuisance_min_delta_Re":gap,
            "separate_nuisance_min_delta_Re":free_gap,"identical_common_baseline_disk_diameter":2*bound,
            "additional_branch_dependent_min_difference_shared_inputs":max(0.,gap-2*bound),
            "separate_inputs_extra_difference_certified_by_Re_alone":max(0.,free_gap-2*bound),
            "separate_inputs_min_delta_Im":cp_gap,
            "separate_inputs_extra_Im_or_weak_decay_axis_if_both_fit":cp_gap-2*obs["combined_im_xi_observed_abs_upper"]},
            "continuous":common}


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true"); args=parser.parse_args()
    cooked=json.loads(json.dumps(results()),parse_float=lambda x:float(format(float(x),".10g")))
    if args.check:
        assert cooked==read("results.json")
    elif args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(cooked,f,indent=2); f.write("\n")
    else:
        print(json.dumps(cooked,indent=2))
    print("1649 PASS: coherent K amplitude, common-baseline separation and continuous insertion test")
