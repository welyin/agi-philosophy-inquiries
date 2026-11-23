"""One pair of vector doublets from Cabibbo data to K and B decays."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def rep(z):
    return [float(z.real),float(z.imag)]


def build(F,p,data,power,complete):
    x,y,r,z=(float(data[k]["mean"]) for k in ("x","y","r","z"))
    m_u,m_c,m_t=p["light_masses_up_GeV"]
    down=np.array(p["down_masses_GeV"])
    a=np.sqrt(p["u_right_components_squared"])
    phase=np.exp(1j*p["light_left_phase"])

    def trial(t):
        ld,ls=(x+t)/2,(y+r*t)/2
        d,s=(x-t)/(2*a[0]),(y-r*t)/(2*a[1])
        b=complex(F/s).conjugate()
        right=a[1]*b.conjugate()
        left=sqrt(z**2-abs(right)**2)/phase
        heavy=np.sqrt([1-d*d,1-s*s-abs(b)**2])
        B=np.array([[d,0],[0,s],[0,b],[heavy[0],0],[0,heavy[1]]],dtype=complex)
        lheavy=(m_u*a-np.array([down[0]*ld*d,down[1]*ls*s+down[2]*left*b]))/(down[3:]*heavy)
        return np.array([ld,ls,left,*lheavy]),B

    lo,hi=x+.0002,x+.003
    assert np.linalg.norm(trial(lo)[0])<1<np.linalg.norm(trial(hi)[0])
    for _ in range(55):
        mid=(lo+hi)/2
        if np.linalg.norm(trial(mid)[0])<1:
            lo=mid
        else:
            hi=mid
    first,B=trial((lo+hi)/2)
    X=down[:,None]*B
    C=np.eye(2)-np.outer(a,a)
    inv=power(C,-.5)
    heavy2,O=np.linalg.eigh(inv@(X.conj().T@X-m_u**2*np.outer(a,a))@inv)
    A=np.vstack([a,np.zeros((2,2)),O.conj().T@power(C,.5)])
    up=np.array([m_u,m_c,m_t,*np.sqrt(heavy2)])
    s13=abs(first[2])/np.linalg.norm(first[:3])
    c13=sqrt(1-s13**2)
    s12=abs(first[1])/np.linalg.norm(first[:2])
    c12=sqrt(1-s12**2)
    s23=p["light_left_s23"]
    c23=sqrt(1-s23*s23)
    ref=np.array([[-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13,0,0],
                  [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13,0,0]])
    q=np.linalg.qr(np.column_stack((X,first.conj())))[0]
    rows=ref@(np.eye(5)-q@q.conj().T)
    light=np.vstack([first,power(rows@rows.conj().T,-.5)@rows])
    n=np.linalg.svd(light,full_matrices=True)[2][3:]
    V=np.vstack([light,(up[3:,None]*A[3:])@np.linalg.solve(n@X,n)])
    R=A@B.conj().T
    f=B@B.conj().T
    err=np.max(abs((np.diag(up)@complete(A))[:,-2:]-(V@np.diag(down)@complete(B))[:,-2:]))
    np.testing.assert_allclose(V@V.conj().T,np.eye(5),atol=2e-10)
    np.testing.assert_allclose([f[0,1],f[0,2],f[1,2]],[0,0,F],atol=1e-14)
    assert err<1.e-8
    np.testing.assert_allclose([abs(V[0,0]+R[0,0]),abs(V[0,1]+R[0,1]),
      abs((V[0,1]-R[0,1])/(V[0,0]-R[0,0])), np.hypot(abs(V[0,2]),abs(R[0,2]))],[x,y,r,z],atol=1e-13)
    return V,R,f,float(err)


def run():
    c=read("inputs.json")
    p=read(c["construction"])
    rate=read(c["rare_B"])
    data=read("../1465/"+p["three_rates"])
    k=read(c["kaon"])
    h=read(c["strange_normalizations"])
    m=read(c["beauty_observations"])
    tau=read(c["charged_tau"])["charged_tau_control"]
    lib=runpy.run_path(str((HERE/c["shared_utilities"]).resolve()))
    power,complete=lib["hermitian_power"],lib["complete_basis"]
    target=h["nu_reference_amplitude"]*sqrt(rate["BelleII"]["mu"])
    g=k["Xtop"][0]/(pi*k["inverse_alpha_MZ"]*k["sW2"])
    kc=rate["sm_coefficient_ratio"]["center"]
    def left_B(V):
        return -V[2,1].conjugate()*V[2,2]*g
    def beta(V):
        return np.angle(-(V[1,0]*V[1,2].conjugate())/(V[2,0]*V[2,2].conjugate()))
    def solve_direction(sign, fixed_angle=None):
        # Match the observed K+ missing-energy strength; beta is the phase from this same V, not a fitted muon phase.
        F=0j
        for _ in range(35):
            V,R,f,err=build(F,p,data,power,complete)
            L=left_B(V)
            ang=sign*beta(V) if fixed_angle is None else fixed_angle
            mag=-abs(L)*np.cos(ang)+sqrt(target**2-abs(L)**2*np.sin(ang)**2)
            new=mag*L/abs(L)*np.exp(1j*ang)
            if abs(new-F)<1e-16:
                break
            F=new
        else:
            raise ArithmeticError("branch did not converge")
        assert abs(abs(L+F)-target)<1e-13
        return F,V,R,f,err

    results={}
    choices=[(str(s),s,None) for s in c["signed_orientations"]]
    choices.extend((f"free_{angle}",0,angle) for angle in c["retrospective_free_phase_rad"])
    for label,sign,fixed in choices:
        F,V,R,f,err=solve_direction(sign,fixed)
        L=left_B(V)
        P=kc*L-F
        angle=np.angle(P/L)
        def mu_rate(w):
            return h["Hmu_heavy"]*abs(P)**2*(1-(1-w)*np.sin(angle)**2)*1e9
        ls=V[2,1].conjugate()*V[2,0]
        charm=V[1,1].conjugate()*V[1,0]
        lam=k["lambda_reference"]
        SK=ls*k["Xtop"][0]/lam**5+charm.real*k["Pcharm"][0]/lam
        kaon=[k["kplus_e11"][0]*(1+k["delta_EM"])*abs(SK)**2,k["kL_e11"][0]*SK.imag**2]
        # All orientation changes preserve the two K projections. Bound all uncomputed common vector amplitudes.
        kp,ep=k["kplus_e11"]
        min_radius=sqrt(6/((kp-2*ep)*(1+k["delta_EM"])))
        max_radius=sqrt(13.4/((kp+2*ep)*(1+k["delta_EM"])))
        remaining_K=min(abs(SK)-min_radius,max_radius-abs(SK))-2*abs(ls)*k["Xtop"][1]/lam**5-2*abs(charm.real)*k["Pcharm"][1]/lam
        assert remaining_K>0 and kaon[1]<6
        left,right=V[0,2],R[0,2]
        zsq=abs(left)**2+abs(right)**2
        # Unknown LR integrated Gram overlap; expose tolerance for the inclusive extraction.
        zinfo=data["z"]
        dz=np.linalg.norm([float(s) for s in zinfo["sigma_components"]])
        z=sqrt(zsq)
        gram_bound=(zsq-(z-dz)**2)/(2*abs(left*right))
        axial=abs(left-right)**2/zsq
        vector=abs(left+right)**2/zsq
        mt=np.mean(tau["mass_tau_GeV_window"])
        mb=tau["mass_Bplus_GeV"][0]
        Btau=tau["tau_Bplus_ps"][0]/c["hbar_GeV_ps"]*m["const_GeV"]["GF"]**2*tau["f_Bplus_GeV"][0]**2*mt**2*mb*(1-(mt/mb)**2)**2/(8*pi)*abs(left-right)**2
        masses=m["const_GeV"]
        Lt_d=-g*V[2,0].conjugate()*V[2,2]
        fbs=read("../1464/inputs.json")["fBs_GeV"][0]
        ratio=m["tau_ps"]["d"][0]/m["tau_ps"]["s_heavy"][0]*(c["fBd_GeV"]/fbs)**2*sqrt((masses["M_Bd"]**2-4*masses["m_mu"]**2)/(masses["M_Bs"]**2-4*masses["m_mu"]**2))
        b_d=1e10*h["Hmu_heavy"]*ratio*abs(kc*Lt_d)**2
        cms=mu_rate(h["calculated_light_response"])
        lhcb=mu_rate(h["run3_central_light_response"])
        # Any additional P amplitude shifts sqrt(the fitted rate) by at most sqrt(Hmu)*|delta P|.
        def safe_margin(pred,window):
            margin=min(sqrt(pred)-sqrt(window[0]),sqrt(window[1])-sqrt(pred))
            return margin/sqrt(1e9*h["Hmu_heavy"])
        mu_margin=min(safe_margin(cms,[3.83-.65,3.83+.71]),safe_margin(lhcb,[3.31-.43,3.31+.46]))
        # Exact phase interval at THIS left modulus, independently of how the UV model produces it.
        # psi=arg(N/L). CMS and LHCb use different acceptance functions.
        def angle_for_rate(value,w):
            def rate_at(u):
                return h["Hmu_heavy"]*1e9*(((1+kc)*abs(L)-target*u)**2+w*target**2*(1-u*u))
            assert (1+kc)*abs(L)>(1-w)*target and rate_at(1)<value<rate_at(-1)
            lo,hi=-1.,1.
            for _ in range(62):
                u=(lo+hi)/2
                if rate_at(u)>value:
                    lo=u
                else:
                    hi=u
            return np.arccos((lo+hi)/2)
        sources=[(h["calculated_light_response"],[3.83-.65,3.83+.71]),
                 (h["run3_central_light_response"],[3.31-.43,3.31+.46])]
        psi_low=max(angle_for_rate(bounds[0],w) for w,bounds in sources)
        psi_high=min(angle_for_rate(bounds[1],w) for w,bounds in sources)
        assert psi_low<psi_high
        belle_error=min(target-h["nu_reference_amplitude"]*sqrt(4.6-1.3),
                        h["nu_reference_amplitude"]*sqrt(4.6+1.3)-target)
        results[label]={"beta_rad":beta(V), "actual_phase_abs_arg_N_over_L":abs(np.angle((L+F)/L)),
            "conditional_phase_abs_arg_N_over_L_common_interval":[psi_low,psi_high],
            "Belle_one_width_additional_vector_amplitude_margin":belle_error,"Fsb":rep(F),"Fbb":float(f[2,2].real),
            "Fss":float(f[1,1].real),"right_ud_us_ub":[rep(z) for z in R[0,:3]],"left_ub":rep(left),
            "left_B_amplitude":rep(L),"P_mu":rep(P),"common_columns_error_GeV":err,
            "delta_Qp_Qcs":[2*.0022-f[0,0].real,188*.0022-211*f[0,0].real],
            "K_plus_and_long_times_1e11":kaon,"extra_K_amplitude_sufficient_two_width":remaining_K,
            "Belle_neutrino_strength":abs(L+F)**2/h["nu_reference_amplitude"]**2,
            "fitted_Bs_times_1e9_CMS_Run3":[cms,lhcb],
            "signed_minimum_margin_in_mu_amplitude_for_both_one_widths":mu_margin,
            "Bd_mumu_times_1e10":b_d,"charged_tau_times_1e4":Btau*1e4,
            "tau_cascade_times_1e6":Btau*c["tau_to_Knu_branch"]*1e6,
            "cascade_shift_in_Belle_strength_if_equal_efficiency":(Btau*c["tau_to_Knu_branch"]*1e6-rate["BelleII"]["tau_background_in_1e-6"])/rate["BelleII"]["reference_SD_in_1e-6"],
            "charged_vector_axial_over_inclusive":[vector,axial],"inclusive_gram_overlap_upper_assuming_equal_diagonals":gram_bound,
            "exact_min_max_abs_P_over_arbitrary_phase_at_this_L_and_target":[abs((1+kc)*abs(L)-target),(1+kc)*abs(L)+target]}
    return json.loads(json.dumps(results),parse_float=lambda x:float(f"{float(x):.10g}"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=run()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as fp:
            json.dump(data,fp,indent=2)
            fp.write("\n")
    else:
        assert data==read("results.json")
        print("1467: one five-dimensional K/B reconstruction and independent flavour rates checked")
