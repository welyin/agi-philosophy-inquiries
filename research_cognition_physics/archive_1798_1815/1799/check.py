"""Actual charged-spectrum embedding with two symmetric matter tensors."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
UP = 0
DOWN = 1
LEPTON = 2


def backend():
    gen=runpy.run_path(str(HERE.parent/"1715/check.py"))
    fl,rep=gen["charged_background"](.1)
    transfer,summary=gen["congruence"](fl,rep,.1)
    old_e=gen["unpack"](rep["charged_left_10PeV"])
    t=old_e.conj().T@transfer
    t=t*(np.diag(t).conj()/abs(np.diag(t)))[None,:]
    return (np.array(rep["Y_up_at_10PeV"]),np.array(rep["Y_down_at_10PeV"]),
            np.array(rep["Y_charged_at_10PeV"]),
            fl["low_ckm"](rep["gamma_10PeV_deg"],rep["CKM_angles_at_10PeV_rad"]),t)


def phase_fix(u):
    # Physical reference convention for matching the additional L-conserving ports.
    a=u[np.argmax(abs(u),axis=0),np.arange(3)]
    return u*(a.conj()/abs(a))[None,:]


def matrices(x, data):
    u,d,e,v,t=data
    tagged=np.exp(1j*np.r_[x[:2],0.])[:,None]*v*np.exp(1j*np.r_[x[2:4],0.])[None,:]
    D=(tagged*d)@tagged.T
    U=np.diag(u)
    w=complex(x[4],x[5])
    E0=D/d[-1]+w*U/u[-1]
    le,ss,_=np.linalg.svd(E0)
    E=np.exp(1j*x[6])*e[-1]*E0/ss[0]
    Ue=phase_fix(le[:,::-1])
    F=(D-E)/4
    C0=Ue.conj().T@F@Ue.conj()
    C=t.conj().T@C0@t.conj()
    return U,D,E,F,C,Ue,ss[::-1]/ss[0]*e[-1],tagged


def neutral(C, dm21=7.388e-5):
    L,ss,_=np.linalg.svd(C)
    s=ss[::-1]
    V=L[:,::-1]
    f=abs(V)**2
    g=f[0,2]
    ratio=(s[1]**2-s[0]**2)/(s[2]**2-s[0]**2)
    s12,s23=f[0,1]/(1-g),f[1,2]/(1-g)
    j=float((V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag)
    den=math.sqrt(max(1e-30,s12*(1-s12)*s23*(1-s23)*g))
    sind=j/((1-g)*den)
    cosd=(f[1,0]-s12*(1-s23)-(1-s12)*s23*g)/(2*den)
    unit=math.sqrt(dm21/(s[1]**2-s[0]**2))
    return np.array([ratio,4*g*(1-g),s12,s23]),dict(
        delta_deg=math.degrees(math.atan2(sind,cosd)),J_CP=j,
        masses_eV=(s*unit).tolist(),sum_eV=float(s.sum()*unit),
        mbb_eV=float(abs(C[0,0])*unit),mbeta_eV=float(np.linalg.norm(C[0])*unit),
        delta_atm_eV2=float(dm21/ratio),common_light_scale_eV=float(unit))


def numbers(x,data):
    U,D,E,F,C,Ue,s,tagged=matrices(x,data)
    nu,additional=neutral(C)
    flat=np.r_[np.log(s[:2]/data[2][:2]),np.log(nu[0]),nu[1:]]
    return flat,additional


def symmetric_sources(x,data,require_charged=True):
    U,D,E,F,C,Ue,s,tagged=matrices(x,data)
    w=complex(x[4],x[5])
    a=np.exp(1j*x[6])*data[2][-1]/(data[1][-1]*np.linalg.norm(D/data[1][-1]+w*U/data[0][-1],2))
    b=a*w*data[1][-1]/data[0][-1]
    rH,rF=(1-a)/b,(-3-a)/b
    H=(3*D+E)/4
    Y=U-4*rF*F
    I=Y@np.linalg.solve(F,Y.T)
    np.testing.assert_allclose(a*D+b*U,E,atol=3e-14)
    np.testing.assert_allclose(rH*H+rF*F,U,atol=5e-14)
    np.testing.assert_allclose(rH*H-3*rF*F,Y,atol=5e-14)
    np.testing.assert_allclose(Y,Y.T,atol=1e-14)
    if require_charged:
        np.testing.assert_allclose(np.linalg.svd(E,compute_uv=False)[::-1],data[2],rtol=2e-8,atol=0)
    np.testing.assert_allclose(abs(tagged),abs(data[3]),atol=1e-14)
    return a,b,rH,rF,Y,I


def mixed(x,data,require_charged=True):
    U,D,E,F,C,Ue,s,tagged=matrices(x,data)
    a,b,rh,rf,Y,I=symmetric_sources(x,data,require_charged)
    z=complex(x[7],x[8])
    K=F/np.linalg.norm(F)+z*I/np.linalg.norm(I)
    L=Ue@data[4]
    m=L.conj().T@K@L.conj()
    obs,more=neutral(m)
    pack=lambda x:[float(x.real),float(x.imag)]
    return dict(rH=pack(rh),rF=pack(rf),a=pack(a),b=pack(b),relative_normalized_complex_I=pack(z),
                dirac_neutrino_singular=np.linalg.svd(Y,compute_uv=False).tolist(),
                norm_F=float(np.linalg.norm(F)),norm_I=float(np.linalg.norm(I)),
                relative_residual_norm_at_high=float(np.linalg.norm(K)/(1+abs(z))),
                charged_relative_max=float(max(abs(s/data[2]-1))),
                low_ratio_A13_solar_atmospheric=obs.tolist(),low=more)


def pure_cover(x,data,r_interval,ncells=32768):
    U,D,E,F,C,Ue,s,tagged=matrices(x,data)
    L=Ue@data[4]
    A=L.conj().T@D@L.conj()/4
    B=L.conj().T@E@L.conj()/4
    nb=np.linalg.norm(B,2)
    angle=np.arange(ncells)*2*np.pi/ncells
    def spectrum(ph):
        return np.linalg.svd(A-np.exp(1j*ph)[...,None,None]*B)
    u,d,_=spectrum(angle)
    d=d[:,::-1]
    # Nearest center within pi/ncells; Weyl for singular values.
    eps=math.pi/ncells*nb+1e-14*(np.linalg.norm(A,2)+nb)
    lo=np.maximum(d-eps,0)**2
    hi=(d+eps)**2
    rlo=np.maximum(lo[:,1]-hi[:,0],0)/np.maximum(hi[:,2]-lo[:,0],1e-100)
    rhi=(hi[:,1]-lo[:,0])/np.maximum(lo[:,2]-hi[:,0],1e-100)
    mask=(rlo<r_interval[1])&(rhi>r_interval[0])
    assert np.any(mask)
    # Project the exact top eigenvector of G onto its central orthogonal complement.
    err=eps*(2*d[:,2]+eps)
    gap=d[:,2]**2-d[:,1]**2
    eta=np.minimum(1,err/np.maximum(gap-err,1e-100))
    pe,pm=abs(u[:,0,0])**2,abs(u[:,1,0])**2
    atm_upper=np.minimum(1,(pm+eta)/np.maximum(1-pe-eta,1e-100))
    atm_lower=np.maximum(0,(pm-eta)/np.maximum(1-pe+eta,1e-100))
    reactor_low=np.maximum(0,pe-eta)
    reactor_high=np.minimum(1,pe+eta)
    assert np.max(reactor_high[mask])<.5
    ar=[float(min((4*reactor_low*(1-reactor_low))[mask])),
        float(max((4*reactor_high*(1-reactor_high))[mask]))]
    atm=[float(min(atm_lower[mask])),float(max(atm_upper[mask]))]
    assert atm[1]<.49
    center=.00007388/.002509
    def ratio(t):
        q=np.linalg.svd(A-np.exp(1j*t)*B,compute_uv=False)[::-1]**2
        return (q[1]-q[0])/(q[2]-q[0])-center
    ratios=(d[:,1]**2-d[:,0]**2)/(d[:,2]**2-d[:,0]**2)-center
    roots=[]
    for j in np.flatnonzero(ratios*np.roll(ratios,-1)<0):
        low=angle[j]; high=low+2*np.pi/ncells
        left=ratio(low)
        for _ in range(48):
            mid=(high+low)/2
            if left*ratio(mid)>0:
                low=mid
            else:
                high=mid
        xx=x.copy();xx[6]=mid
        q,more=neutral(matrices(xx,data)[4])
        assert abs(q[0]-center)<1e-12
        roots.append(dict(chi_rad=mid,three_angles=q[1:].tolist(),delta_deg=more["delta_deg"]))
    return dict(relative_charged_spectral_residual=float(max(abs(s/data[2]-1))),
                chi_cells=ncells,retained_cells=int(mask.sum()),A13_enclosure=ar,
                sin2theta23_enclosure=atm,central_frequency_crossings=roots)


def calculate():
    cfg=json.loads((HERE/"inputs.json").read_bytes())
    data=backend()
    for i,k in enumerate(["up","down","lepton"]):
        np.testing.assert_allclose(data[i],cfg[k],rtol=3e-12)
    old=json.loads((HERE/ cfg['frequencies']).read_bytes())
    l0,lm,lp=old['juno_NO_dm21'];h0,hm,hp=old['juno_NO_dm31']
    rr=((l0-2*lm)/(h0+2*hp),(l0+2*lp)/(h0-2*hm))
    wins={}
    for label,x in cfg["mixed_type_I_II"].items():
        r=mixed(np.array(x),data)
        obs=r["low_ratio_A13_solar_atmospheric"]
        np.testing.assert_allclose(obs,[l0/h0,.0855,.309,.55],atol=3e-9)
        assert abs(r["low"]["delta_deg"]+int(label.replace("delta_minus","")))<1e-6
        # Constructive local recalibration: fix two family phase labels; fit the two
        # independent charged ratios, the mass ratio, three angles, and signed delta.
        fixed_delta=r["low"]["delta_deg"]
        chosen=np.array([1,3,4,5,6,7,8])
        def residual(z):
            values=mixed(z,data,require_charged=False)
            s=matrices(z,data)[6]
            q,ar,sol,atm=values["low_ratio_A13_solar_atmospheric"]
            return np.r_[np.log(s[:2]/data[2][:2]),np.log(q/(l0/h0)),ar-.0855,
                         sol-.309,atm-.55,np.radians(values["low"]["delta_deg"]-fixed_delta)]
        basis=np.eye(9)[chosen]
        x1=np.array(x)+np.sum(basis,axis=0)*1e-6
        for _ in range(4):
            h=2e-7
            j=np.column_stack([(residual(x1+h*d)-residual(x1-h*d))/(2*h) for d in basis])
            x1[chosen]-=np.linalg.solve(j,residual(x1))
        fit_residual=np.linalg.norm(residual(x1))
        assert fit_residual<1e-8
        r["local_seven_equation_Newton"]={"smallest_singular_value":np.linalg.svd(j,compute_uv=False)[-1],
                                          "residual_norm":fit_residual,"largest_recovered_coordinate_difference":max(abs(x1-x))}
        wins[label]=r
    counter={}
    for label,x in cfg["fixed_delta_minus70_mass_readout"].items():
        row=mixed(np.asarray(x),data)
        np.testing.assert_allclose(row["low_ratio_A13_solar_atmospheric"],[l0/h0,.0855,.309,.55],atol=2e-9)
        assert abs(row["low"]["delta_deg"]+70)<1e-6
        counter[label]=row
    assert counter["larger_mbb"]["low"]["mbb_eV"]/counter["smaller_mbb"]["low"]["mbb_eV"]>1.15
    pencils={}
    for label,x in cfg["pencils"].items():
        symmetric_sources(np.array(x),data)
        pencils[label]=pure_cover(np.array(x),data,rr)
    return dict(normal_ordering_two_width_frequency_ratio=rr,
                transfer_distance_after_low_row_rephase=np.linalg.norm(data[4]-np.eye(3),2),
                exact_charged_finite_precision_witnesses_and_all_chi_pureII=pencils,
                same_two_tensor_mixed_type_I_and_II=wins,
                at_identical_signed_CP_and_identical_charged_and_oscillation_inputs=counter)


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true");args=ap.parse_args()
    out=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    target=HERE/"results.json"
    if args.check:
        assert json.loads(target.read_bytes())==json.loads(out)
    else:
        target.write_text(out,encoding="utf-8")
    print("1799 OK",hashlib.sha256(out.encode()).hexdigest()[:12])

