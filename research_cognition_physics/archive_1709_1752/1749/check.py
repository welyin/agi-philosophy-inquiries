"""Charged-current dipole, including both real partners of each neutral Dirac pair."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())
CFG=read("inputs.json")
PORT=runpy.run_path(str(HERE/CFG["port"]))
def svd(y):
    left,s,rh=np.linalg.svd(y)
    return left[:,::-1],s[::-1],rh.conj().T[:,::-1]


def functions(x):
    g=np.zeros_like(x)
    h=np.full_like(x,2.)
    use=x>1e-22
    z=x[use]
    g[use]=-z*(2*z*z+5*z-1)/(4*(1-z)**3)-1.5*z**3*np.log(z)/(1-z)**4
    h[use]=(z**3-12*z*z+15*z-4+6*z*z*np.log(z))/(2*(z-1)**3)
    g[~use]=x[~use]/4
    return g,h


def common_matrices(b,a,port,F,derived):
    u=CFG["v_GeV"]/math.sqrt(2)
    mass=b["Mweak"]
    y=u*a*b["yW"][2,:,:3]
    n=u*b["N"]
    cg=np.array(port["physical_mass_rows_over_M_real"])+1j*np.array(port["physical_mass_rows_over_M_imag"])
    Mch=np.column_stack([np.vstack([y,np.zeros((2,3))]),mass*cg.conj().T])
    rc,ec,sc=svd(Mch)
    paired=np.array(b["pt"]["flow"]["fields"]["neutral_poles_GeV"]).reshape(2,2)
    peaks=paired.mean(axis=1)
    assert max(np.ptp(paired,axis=1))<.001
    # In the original gauge basis light and vector-doublet fields have exactly the
    # same charged left current. The two other left states are neutral weak singlets.
    MN=np.column_stack([np.vstack([mass*cg.conj().T,np.zeros((2,2))]),
                        np.vstack([n,np.zeros((2,2)),np.diag(peaks)])])
    pt=b["pt"]; T=np.linalg.inv(b["S"]@pt["low_E"])
    C=F.T@T.T@pt["C"]@T@F
    vec=F.T@T.T@pt["b"]
    K=C+b["nu"]["FROZEN"][0]*np.outer(vec,vec)
    solar_GeV=derived["neutrinos"]["scale_eV"]*1e-9
    P=(F@n).conj()/peaks
    pplus=np.linalg.pinv(P)
    mu=pplus@(solar_GeV*K)@pplus.T
    assert np.linalg.norm(P@mu@P.T-solar_GeV*K)<1e-9*np.linalg.norm(solar_GeV*K)
    return rc,ec,sc,peaks,MN,mu


def dipole(rc,ec,sc,peaks,MN,mu,weak_doublet_mass):
    # Matrix in left-handed Weyl notation. This 11x11 permits LNV without
    # independently choosing a phase for either member of a quasi-Dirac cluster.
    M=np.zeros((11,11),complex)
    M[:7,7:]=MN.conj()
    M[7:,:7]=MN.conj().T
    M[5:7,5:7]=mu
    # SVD supplies both chiralities with identical unitary conventions even on a
    # degenerate singular subspace. No arbitrary Takagi column rephasing is used.
    ol,d,orr=svd(M)
    left=rc.conj().T@orr[:5]
    right=sc[3:,:].conj().T@ol[7:9]
    np.testing.assert_allclose(left@left.conj().T,np.eye(5),atol=3e-12)
    np.testing.assert_allclose((left*d)@right.conj().T,ec[:,None]*(sc[3:].conj().T@sc[3:]),atol=2e-10)
    g,h=functions((d/CFG["weak_mass_GeV"])**2)
    pref=3/(CFG["inverse_alpha"]*2*np.pi)
    light=sum(abs(left[:2,:3].T)**2)
    denom=np.prod(light)
    pieces=[]
    for i,j in [(0,1),(1,0)]:
        ll=ec[j]/ec[1]*left[i]*left[j].conj()*g
        rr=ec[i]/ec[1]*right[i]*right[j].conj()*(g-5/6)
        lr=left[i]*right[j].conj()*d*h/ec[1]
        pieces.append(np.array([ll,rr,lr]))
    pieces=np.array(pieces)
    both=np.sum(pieces,axis=1)
    al=np.sum(pieces[:,0,:],axis=1)
    amp=np.sum(both,axis=1)
    labels=np.where(d<1.,0,1+np.argmin(abs(d[:,None]-np.array([*peaks,weak_doublet_mass])[None,:]),axis=1))
    assert list(np.bincount(labels))==[3,2,2,4]
    # Scalar products of two physical cluster amplitudes, invariant under charged-field phases.
    sub=np.array([np.sum(pieces[:,0,labels==k],axis=1) for k in range(1,4)])
    phases=np.vdot(sub[0],sub[1])/(np.linalg.norm(sub[0])*np.linalg.norm(sub[1]))
    return dict(BR_W_all_chiralities=pref*np.linalg.norm(amp)**2/denom,
                BR_W_left_only=pref*np.linalg.norm(al)**2/denom,
                cluster_mass_ranges_GeV=[[min(d[labels==k]),max(d[labels==k])] for k in range(1,4)],
                BR_left_one_cluster_at_a_time=pref*np.linalg.norm(sub,axis=1)**2/denom,
                singlet_cluster_interference_cosine=phases.real,
                physical_light_weights=light,LL_cluster_amplitude_moduli=np.linalg.norm(sub,axis=1),
                relative_other_chiralities=np.linalg.norm(amp-al)/np.linalg.norm(al),
                right_light_weight_max=np.linalg.norm(right[:2,:3],2)**2)


def pole_ports(b,F,a):
    y=b["yW"].copy()
    y[2]=a*F@y[2]
    y[3]=F@y[3]
    src={}
    flow=b["pt"]["flow"]
    poles=flow["fields"]["neutral_poles_GeV"]
    model=b["model"]
    for scale,susy,dt,g in reversed(b["segments"]):
        y,_=model["advance"](y,np.eye(3,dtype=complex),-dt,g[::-1],flow)
        if scale==b["Mcolor"]:
            l,_,r=svd(y[1,:,:3])
            c=abs(svd(y[0,:,:3])[0][:,0].conj()@l[:,0])
            y[1,:,:3]=y[1,:,:3]@((r*[c,1.,math.sqrt((1+c)/2)])@r.conj().T)
        for i,pole in enumerate(poles):
            if scale==pole:
                if i in (0,2):
                    src[i]=y[3,:,i].copy()
                y[3,:,i]=0.
    l,_,_=svd(y[2,:,:3])
    return l,np.column_stack([src[0],src[2]])


def leading(z,peaks):
    u=CFG["v_GeV"]/math.sqrt(2)
    mix=z*u/peaks
    g,_=functions((peaks/CFG["weak_mass_GeV"])**2)
    amp=(mix[0]*mix[1].conj())@g
    return 3/(CFG["inverse_alpha"]*2*np.pi)*abs(amp)**2


def scan(step):
    b=PORT["OLD"]["mother"](step)
    out=[]
    for j in range(2):
        for signs in PORT["CFG"]["family_signs"]:
            X,F,port,_=PORT["terms"](b,j,signs)
            a,derived=PORT["daughter"](b,F)
            rc,ec,sc,peaks,MN,mu=common_matrices(b,a,port,F,derived)
            base=dipole(rc,ec,sc,peaks,MN,0*mu,b["Mweak"])
            broken=dipole(rc,ec,sc,peaks,MN,mu,b["Mweak"])
            low,src=pole_ports(b,F,a)
            layer=svd(F@b["yW"][2,:,:3])[0]
            pred=leading(layer.conj().T@F@b["N"],peaks)
            mixed=leading(low.conj().T@src,peaks)
            base.update(root=j+1,signs=signs,physical_tau_calibration=a,
                example_terminal_LNV_spectral_norm_GeV=np.linalg.norm(mu,2),
                BR_with_finite_example_LNV=broken["BR_W_all_chiralities"],
                relative_LNV_numerical_shift=(broken["BR_W_all_chiralities"]/base["BR_W_all_chiralities"]-1),
                left_leading_common_scale_BR=pred,left_leading_with_each_column_at_own_pole_and_MZ_frame=mixed,
                pole_over_common_leading_rate=mixed/pred,
                ratio_to_MEGII_90=base["BR_W_all_chiralities"]/CFG["MEGII_BR90"])
            out.append(base)
    return out


def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):return float(f"{x:.11g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    coarse,fine=[scan(s) for s in CFG["steps"]]
    err=max(abs(c["BR_W_all_chiralities"]/f["BR_W_all_chiralities"]-1) for c,f in zip(coarse,fine))
    assert err<2e-5
    assert all(max(abs(x["relative_LNV_numerical_shift"]),x["right_light_weight_max"])<1e-6 for x in fine)
    slope=fine[1]["BR_W_all_chiralities"]
    result=clean(dict(cases=fine,weak_pole_calculation="Single-scale finite W and Goldstone matching; LL, RR and internal chirality flip on the physical neutral states.",
                      accepted_1748_branch_index=1,coarse_fine_max_rel_BR=err,
                      free_same_convention_dipole_radius_in_units_of_conditional_W=math.sqrt(CFG["MEGII_BR90"]/slope)))
    if args.check:
        assert result==read("results.json")
        print("1749 verified weak doublets, neutral pairs, GIM cancellation and transported pole columns")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")
