"""Physical five charged doublet components, their right Z/H couplings and the mixed W loop."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
MODEL=runpy.run_path(str(HERE/"../1748/check.py"))
VEV=246.22
ALPHA=1/137.036
MW=80.379
MZ=91.1876
WIDTH_Z=2.4952
SW2=.23122


def svd_up(matrix):
    l,d,rh=np.linalg.svd(matrix)
    return l[:,::-1],d[::-1],rh.conj().T[:,::-1]


def dipole(b,F,B,Le,Re,lepton,angle_rotation=None):
    poles=np.array(b["pt"]["flow"]["fields"]["neutral_poles_GeV"])
    peaks=(poles[::2]+poles[1::2])/2
    dn=VEV*b["N"]/math.sqrt(2)
    D=b["Mweak"]*np.eye(2) if angle_rotation is None else b["Mweak"]*angle_rotation
    N=np.block([[np.zeros((3,2)),F@dn],[D,B@dn],[np.zeros((2,2)),np.diag(peaks)]])
    Ln,heavy,Vh=np.linalg.svd(N,full_matrices=True)
    Rn=Vh.conj().T
    masses=np.r_[heavy,np.zeros(3)]
    # Left-handed charged current sees all five weak doublets, right-handed
    # current sees only the two vectorlike neutral components.
    left=Le.conj().T@Ln[:5]
    right=np.column_stack([Re[3:,:].conj().T@Rn[:2,:],np.zeros((5,3))])
    rotated=Re[3:,:].conj().T@Re[3:,:]
    np.testing.assert_allclose((left*masses)@right.conj().T,lepton[:,None]*rotated,atol=1e-10)
    norm_light=np.sum(np.abs(left[:2,4:])**2,axis=1)
    denominator=np.prod(norm_light)
    x=(masses[:4]/MW)**2
    # Hue et al. B1-B5, Q_F=0: use the unsubtracted functions for the right current.
    wmix=np.r_[(x**3-12*x**2+15*x-4+6*x**2*np.log(x))/(4*(x-1)**3),np.ones(3)]
    wext=np.r_[(-4*x**4+49*x**3-78*x**2+43*x-10-18*x**3*np.log(x))/(24*(x-1)**4),
               np.full(3,-10/24)]
    parts=[]
    for i,j in ((0,1),(1,0)):
        ll=lepton[j]*np.dot(left[i]*left[j].conj(),wext-wext[-1])/lepton[1]
        rr=lepton[i]*np.dot(right[i]*right[j].conj(),wext)/lepton[1]
        lr=np.dot(left[i]*right[j].conj()*masses,wmix)/lepton[1]
        parts.append([ll,rr,lr])
    parts=np.array(parts)
    B0=6*ALPHA/math.pi/denominator
    channels=[float(B0*np.sum(abs(parts[:,:j].sum(axis=1))**2)) for j in (1,2,3)]
    return dict(neutral_Dirac_masses_GeV=heavy,light_norm_e_mu=norm_light,ordinary_width_factor=denominator,
                BR_W_LL_RR_mixed_cumulative=channels,
                abs_terms_rows_e_mu__mu_e_cols_LL_RR_LR=np.abs(parts),
                norm_LR_over_LL=np.linalg.norm(parts[:,2])/np.linalg.norm(parts[:,0]),
                mass_sum_heavy_cross_abs_GeV=np.abs(((left*masses)@right.conj().T)[:2,:2]))


def neutral_boson_bound(m,K):
    """Termwise absolute orbital bounds (same normalization as the W form factor)."""
    H=m[:,None]*(np.eye(5)-K)/VEV
    gweak=2*MW/VEV
    gz=2*MZ/VEV
    sectors={}
    for name,pole,L,R in (("Z",MZ,gz*(SW2-.5)*np.eye(5),gz*(SW2*np.eye(5)-K/2)),
                          ("h",125.25,H.conj(),H.T)):
        x=(m/pole)**2
        if name=="Z":
            # The F_Z and G_Z of Kawamura et al. (3.6,3.7); allowance
            # 32 covers the normalizations of both orbital terms.
            direct=(5*x**4-14*x**3+39*x**2-38*x-18*x**2*np.log(x)+8)/(12*(1-x)**4)
            mixed=(x**3+3*x-6*x*np.log(x)-4)/(2*(1-x)**3)
            cap=np.full(5,32.)
            fv=(x**3-12*x*x+15*x-4+6*x*x*np.log(x))/(4*(x-1)**3)
            gv=(x*x-5*x+4+3*x*np.log(x))/(2*(x-1)**2)
            ft=(-4*x**4+49*x**3-78*x*x+43*x-10-18*x**3*np.log(x))/(24*(x-1)**4)
            gt=-3*(x**3-6*x*x+7*x-2+2*x*x*np.log(x))/(x-1)**3
            assert np.all(np.maximum(abs(fv-gv),abs(ft-gt))<cap)
        else:
            # B2 of Hue et al., Q_F=-1, neutral real scalar.
            mixed=(x*x-1-2*x*np.log(x))/(4*(x-1)**3)-(x-1-np.log(x))/(2*(x-1)**2)
            direct=(2*x**3+3*x*x-6*x+1-6*x*x*np.log(x))/(24*(x-1)**4)-(x*x-1-2*x*np.log(x))/(8*(x-1)**3)
            cap=32/np.maximum(1,x)
        assert np.all(2*np.maximum(abs(mixed),abs(direct))<cap)
        bounds=[]
        for i,j in ((0,1),(1,0)):
            # Rows are external flavour, columns are internal charged eigenstate.
            products=np.abs(L[i]*R[j])*m+abs(L[i]*L[j])*m[j]+abs(R[i]*R[j])*m[i]
            bounds.append(2*(MW/pole)**2*(cap@products)/(gweak**2*m[1]))
        sectors[name]=np.linalg.norm(bounds)
    return sectors


def compute(step):
    b=MODEL["OLD"]["mother"](step)
    X,F,port,_=MODEL["terms"](b,0,(1,-1,1))
    a,row=MODEL["daughter"](b,F)
    Y=a*b["yW"][2,:,:3]
    M=b["Mweak"]
    r=np.array(port["physical_mass_rows_over_M_real"])+1j*np.array(port["physical_mass_rows_over_M_imag"])
    B=r[:,:3]
    np.testing.assert_allclose(B.conj().T@B,np.eye(3)-F@F,atol=4e-15)
    hats=VEV*Y/math.sqrt(2)
    A=F@hats
    # Unitary embedding of the 3 light and 2 heavy left components into
    # gauge coordinates. This independently reconstructs the full pre-rotation mass.
    light=np.vstack([np.eye(3),-X])@F
    change=np.column_stack([light,r.conj().T])
    np.testing.assert_allclose(change.conj().T@change,np.eye(5),atol=4e-15)
    mass=np.block([[A,np.zeros((3,2))],[B@hats,M*np.eye(2)]])
    raw=np.column_stack([np.vstack([hats,np.zeros((2,3))]),M*r.conj().T])
    np.testing.assert_allclose(change@mass,raw,atol=4e-11)
    Le,m,Re=svd_up(mass)
    P=np.diag([0,0,0,1,1])
    K=Re[3:,:].conj().T@Re[3:,:]
    exact_h=Le.conj().T@(mass@(np.eye(5)-P))@Re/VEV
    np.testing.assert_allclose(exact_h,m[:,None]*(np.eye(5)-K)/VEV,atol=3e-14)
    # Solve algebraically for the small heavy entries of each right singular
    # vector; their dot product avoids roundoff from subtracting large matrices.
    right=Re[:3,:3]
    heavy=-M*(B@hats@right)/(M*M-m[:3]**2)
    Ks=heavy.conj().T@heavy
    np.testing.assert_allclose(K[:3,:3],Ks,rtol=1e-6,atol=1e-22)
    left,d,right0=svd_up(A)
    G=left.conj().T@np.linalg.solve(F,np.linalg.solve(F,left))-np.eye(3)
    lead=G*np.outer(d,d)/M**2
    np.testing.assert_allclose(abs(Ks),abs(lead),rtol=1e-6,atol=1e-22)
    # MZ convention of the tree-level observables and weak-pole Yukawa convention are explicit.
    normalization=dipole(b,F,B,Le,Re,m)
    coupling=abs(Ks[0,1])/2
    GF=1/(math.sqrt(2)*VEV**2)
    zBR=2*GF*MZ**3/(3*math.sqrt(2)*math.pi*WIDTH_Z)*coupling**2
    geL=-.5+SW2
    geR=SW2-Ks[0,0].real/2
    width=4*(2*geR**2+geL**2)*coupling**2/normalization["ordinary_width_factor"]
    # Direct rotation of the full v-derivative subtracts O(M) contributions;
    # use the exact mass-projector identity for tiny light off-diagonals.
    real_h=abs(m[:3,None]*(np.eye(3)-Ks)/VEV)
    factors=neutral_boson_bound(m,K)
    ll_norm=np.linalg.norm(normalization["abs_terms_rows_e_mu__mu_e_cols_LL_RR_LR"][:,0])
    zhover={name:val/ll_norm for name,val in factors.items()}
    dipole3factor=ALPHA/(3*math.pi)*(2*math.log(105.6583755/.51099895)-11/4)
    # Joint change of the heavy-left labels also rotates mass rows and both
    # Yukawa bottom blocks, leaving the Z projection and W dipole unchanged.
    u=np.array([[1,1j],[1j,1]],complex)/np.sqrt(2)
    rot=np.eye(5,dtype=complex)
    rot[3:,3:]=u
    Le2,m2,Re2=svd_up(rot@mass)
    alt=dipole(b,F,u@B,Le2,Re2,m2,u)
    np.testing.assert_allclose(alt["BR_W_LL_RR_mixed_cumulative"],normalization["BR_W_LL_RR_mixed_cumulative"],rtol=2e-6)
    np.testing.assert_allclose(abs(Re2[3:,:3].conj().T@Re2[3:,:3]),abs(Ks),rtol=2e-6,atol=1e-22)
    return dict(weak_pole_GeV=M,charged_scale=a,sshape_central_solar=[row["neutrinos"]["dm31_eV2"],row["neutrinos"]["dm32_eV2"]],
                prebreak_leading_charged_masses_GeV=d,broken_charged_masses_GeV=m,
                right_doublet_overlap_abs=np.abs(Ks),neutral_coupling_abs_e_mu=coupling,
                gram_G_abs=np.abs(G),leading_vs_exact_e_mu_rel=abs(abs(lead[0,1]/Ks[0,1])-1),
                higgs_physical_vertex_abs=real_h,BR_Z_to_e_mu_both_signs=zBR,
                effective_Z_only_ratio_mu_to_3e=width,Z_to_e_mu_ATLAS95=2.62e-7,mu3e_SINDRUM90=1e-12,
                neutral_boson_dipole_over_left_W_norm_upper=zhover,
                three_body_interconversion_if_only_W_dipole=dipole3factor*normalization["BR_W_LL_RR_mixed_cumulative"][-1],
                maximal_common_W_amplitude_change_from_RR_and_LR=(np.linalg.norm(normalization["abs_terms_rows_e_mu__mu_e_cols_LL_RR_LR"][:,1],2)
                    +np.linalg.norm(normalization["abs_terms_rows_e_mu__mu_e_cols_LL_RR_LR"][:,2],2))
                    /np.linalg.norm(normalization["abs_terms_rows_e_mu__mu_e_cols_LL_RR_LR"][:,0],2),
                W_chiral_check=normalization)


def clean(o):
    if isinstance(o,dict):return {k:clean(v) for k,v in o.items()}
    if isinstance(o,np.ndarray):return clean(o.tolist())
    if isinstance(o,list):return [clean(v) for v in o]
    if isinstance(o,(float,np.floating)):return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    cfg=p.parse_args()
    fine,coarse=compute(.125),compute(.25)
    ratio=abs(fine["BR_Z_to_e_mu_both_signs"]/coarse["BR_Z_to_e_mu_both_signs"]-1)
    assert ratio<1e-5
    result=clean(dict(fine=fine,step_halving_relative_change_Z=ratio,
                conventions=dict(v_GeV=VEV,alpha=ALPHA,MW=MW,MZ=MZ,widthZ=WIDTH_Z,sin2_weak=SW2,
                        sources=["https://arxiv.org/abs/2204.10783","https://doi.org/10.1016/0550-3213(88)90462-2",
                        "https://doi.org/10.1140/epjc/s10052-024-13624-6"])))
    path=HERE/"independent_results.json"
    if cfg.check:
        archived=json.loads(path.read_bytes())
        assert archived==result
        print("1749 vectorlike right Z, aligned Higgs and both W chiralities verified")
    else:
        with path.open("w",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")

