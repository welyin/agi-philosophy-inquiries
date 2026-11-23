"""Hermitian surgery confined to the two physical Dirac columns."""
import json
import math
import runpy
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())
CFG=read("inputs.json")
MODEL=runpy.run_path(str(HERE/CFG["model"]))
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
I=np.eye(2)

def setup(step=.125):
    b=MODEL["mother"](step)
    U=np.linalg.qr(b["N"],mode="complete")[0]
    y=U.conj().T@b["yW"][2,:,:3]
    l,e,_=np.linalg.svd(y)
    weights=[(l*w)@l.conj().T for w in (e**2,e**-2)]
    weights=[w/w.trace().real for w in weights]
    pt=b["pt"]
    inv=np.linalg.inv(b["S"]@pt["low_E"])
    mass=inv.T@(pt["C"]+b["nu"]["FROZEN"][0]*np.outer(pt["b"],pt["b"]))@inv
    K=U.T@mass@U
    assert np.linalg.norm(K[:,-1])<1e-12*np.linalg.norm(K)
    return b,U,y,e,weights,K

def orientations(weights, spectrum, rel=(1.,1.)):
    f1,f2,f0=spectrum
    g=(f1*f2*f0)**(1/3)
    A=[]; rhs=[]
    for p,desired_frac,w in zip((2,-2),rel,weights):
        mean=(f1**p+f2**p)/2
        diff=(f1**p-f2**p)/2
        axis=np.einsum('aij,ji->a',PAULI,w[:2,:2]).real
        A.append(axis)
        rhs.append((g**p*desired_frac-mean*np.trace(w[:2,:2]).real-f0**p*w[2,2].real)/diff)
    A=np.array(A)
    p=np.linalg.lstsq(A,rhs,rcond=None)[0]
    rad=1-p@p
    if rad < 0:
        return [],p,rad,A
    direction=np.cross(A[0],A[1])
    direction/=np.linalg.norm(direction)
    return [p+sign*np.sqrt(rad)*direction for sign in (-1,1)],p,rad,A

def gate(axis,spectrum):
    f1,f2,f0=spectrum
    a=np.zeros((3,3),complex)
    a[:2,:2]=(f1+f2)/2*I+(f1-f2)/2*np.einsum("a,aij->ij",axis,PAULI)
    a[2,2]=f0
    return a


def charged_relative(e,x,y):
    d=np.array([1.,y,x])/(x*y)**(1/3)
    return [sum((e*d)**p)/sum(e**p) for p in (2,-2)]


def trace_box(e,weights,spec,lower,upper):
    """Enclose the smallest possible squared Bloch norm for every charged-ratio cell."""
    intervals=[np.linspace(x,y,25) for x,y in zip(lower,upper)]
    mid=list(spec)
    rods=orientations(weights,mid)[3]
    pinv=np.linalg.pinv(rods)
    assert np.linalg.svd(rods,compute_uv=False)[-1]>.25
    g=np.prod(mid)**(1/3)
    best=np.inf
    err=1e-10
    for ix in range(24):
        for iy in range(24):
            lx,ux=intervals[0][ix:ix+2]
            ly,uy=intervals[1][iy:iy+2]
            # Desired normalized physical masses from tau to e, after factoring g.
            lo=np.array([1,ly,lx])/(ux*uy)**(1/3)
            hi=np.array([1,uy,ux])/(lx*ly)**(1/3)
            bounds=[]
            for power,w in zip((2,-2),weights):
                dl,du=(lo,hi) if power>0 else (hi,lo)
                tar=np.array([sum((e*z)**power)/sum(e**power) for z in (dl,du)])*g**power
                base=(mid[0]**power+mid[1]**power)/2*np.trace(w[:2,:2]).real+mid[2]**power*w[2,2].real
                half=(mid[0]**power-mid[1]**power)/2
                bounds.append(np.sort((tar-base)/half))
            bounds=np.array(bounds)
            point=pinv@bounds.mean(axis=1)
            size=np.linalg.norm(pinv,axis=0)@((bounds[:,1]-bounds[:,0])/2)
            best=min(best,max(0.,np.linalg.norm(point)-size-err)**2)
    return best


def view(b,U,spec,axis,threshold_charged):
    a=gate(axis,spec)
    F=U@a@U.conj().T
    row=MODEL["daughter"](b,F)
    o=row["neutrinos"]
    _,v=np.linalg.eigh(a[:2,:2])
    line=U[:,:2]@v[:,0]
    norms=np.linalg.norm(b["N"],axis=0)
    physical=b["N"]/norms
    coeff=np.linalg.lstsq(physical,line,rcond=None)[0]
    coeff*=np.conj(coeff[0])/abs(coeff[0])
    assert abs(np.linalg.norm(physical@coeff)-1)<1e-12
    # One row of a heavy-light mixing X. A second, unmixed heavy doublet is allowed.
    return dict(spectrum=spec,axis=axis,charged_target_gains=threshold_charged,
                line_in_normalized_dirac_sources_real=coeff.real,line_in_normalized_dirac_sources_imag=coeff.imag,
                relative_source_phase_deg=np.degrees(np.angle(coeff[1])),
                charged_MZ_after_tau=row["MZ_relative_e_mu_after_tau_calibration"],
                dm31_eV2=o["dm31_eV2"],dm32_eV2=o["dm32_eV2"],
                angles=o["angle_readings"],angle_in_windows=o["angles_in_box"],
                delta_deg=o["delta_deg"],FC90=o["SK_T2K_FC68_90"][1],
                quark_four_ratio_radius=row["absolute_four_ratio_radius"],
                charged_source_residual=row["matching_residual_neutral_columns"])


def solve(step):
    b,U,y,e,w,K=setup(step)
    corr=read(CFG["wide_charged_box"])["arbitrary_positive_weak_rows_tagged_GJ"]
    lower,upper=np.array(corr["allowed_charged_ratios_at_MW"])+np.array([[-1e-10],[1e-10]])
    length,active=corr["integration_lengths"]
    # Eigenvector velocity for the electron. The spectral gap is bounded by the recorded
    # low e/mu ratio enlarged by the full anisotropic flow; use a much wider 0.006.
    cap_e,cap_n=corr["allowed_lepton_norms"]
    etime=1.5*(cap_e**2*length+cap_n**2*active)
    nd=1.5*cap_e**2*length+.5*cap_n**2*active
    ratio_e_mu=b["model"]["JOINT"]["LEPTONS"][0]/b["model"]["JOINT"]["LEPTONS"][1,0]
    assert (ratio_e_mu[0]+2*ratio_e_mu[1])*math.exp(etime)*1.002<.006
    delta=etime*(1+.006**2)/(1-.006**2)+nd
    mixing=b["nu"]["ANG"]
    wemin=(1-mixing[1,1])*(1+math.sqrt(1-mixing[0,1]))/2
    # The small-reactor-angle branch is part of the registered normal-ordering data.
    # A reactor disappearance measurement alone is invariant under theta13 -> pi/2-theta13.
    source_weights=abs(np.linalg.svd(y)[0][-1,::-1])**2
    ei=e[::-1]
    weighted=sum(source_weights*(ei[0]/ei)**2)
    c,h=b["cW"],b["fW"]
    tests={}
    for tag,spec in dict(null_unaffected=(c,h,1.),null_half=(c,1.,h),null_full=(h,1.,c)).items():
        ns,p,rad,A=orientations(w,spec)
        g=np.prod(spec)**(1/3)
        we_upper=weighted*(g/spec[2])**2*upper[0]**(4/3)/lower[1]**(2/3)
        low=trace_box(e,w,spec,lower,upper)
        tests[tag]=dict(spectrum=spec,exact_center_axis_norm_squared=1-rad,
                        all_charged_windows_axis_norm_squared_lower=low,
                        electron_zero_weight_Mw_upper=we_upper,electron_zero_weight_MZ_upper=min(1,we_upper+delta),
                        exact_charged_center_branches=[view(b,U,spec,n,[1,1]) for n in ns])
    assert tests["null_unaffected"]["all_charged_windows_axis_norm_squared_lower"]>1
    assert tests["null_full"]["all_charged_windows_axis_norm_squared_lower"]>1
    assert tests["null_half"]["electron_zero_weight_MZ_upper"]<wemin
    arbitrary_floor=((wemin-delta)/weighted)**1.5*lower[1]/upper[0]**2
    dy=1e-4
    charged=(1+dy,1-dy)
    cases={}
    for tag,second in (("two_passes_on_one_line",1.),("both_lines_active",.99995)):
        spec=(h*h,second,1.)
        axes,_,_,_=orientations(w,spec,charged_relative(e,*charged))
        assert axes
        rows=[view(b,U,spec,n,charged) for n in axes]
        # Explicit candidate points: solar central; both independent atmospheric
        # frequencies, charged ratios and the separate three-angle test.
        table=b["nu"]["FREQ"]
        for out in rows:
            ok=[]
            for field,src in (("dm31_eV2","juno_NO_dm31"),("dm32_eV2","accelerator_NO_dm32")):
                ctr,down,up=table[src]
                ok.append(ctr-2*down<out[field]<ctr+2*up)
            out["two_width_atmospheric_windows"]=ok
            assert max(abs(out["charged_MZ_after_tau"]))<.00025
        assert any(all(r["two_width_atmospheric_windows"]) for r in rows)
        assert all(not r["angle_in_windows"][1] for r in rows)
        cases[tag]=rows
    return dict(poles_GeV=[b["Mweak"],b["Mcolor"]],charged_ratio_transport_limits=[lower,upper],
                zero_line_Mw_charged_weights=source_weights,charged_inverse_moment_electron_units=weighted,
                nr_geometric_leak=np.linalg.norm(K[:,-1])/np.linalg.norm(K),
                charged_null_probability_drift=delta,measured_electron_weight_lower=wemin,
                arbitrary_plane_determinant_lower_if_null_unaffected=arbitrary_floor,
                three_original_strengths=tests,plane_only_new_eigenvalue_witnesses=cases)


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):
        return [tidy(z) for z in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.12g}")
    if isinstance(x,np.bool_):
        return bool(x)
    return x


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    co,fi=[solve(step) for step in CFG["steps"]]
    label="two_passes_on_one_line"
    discrepancy=max(abs(a["dm31_eV2"]-b["dm31_eV2"]) for a,b in zip(
        co["plane_only_new_eigenvalue_witnesses"][label],fi["plane_only_new_eigenvalue_witnesses"][label]))
    assert discrepancy<1e-8
    fine=tidy(dict(**fi,step_difference_dm31=discrepancy))
    if args.check:
        assert fine==read("results.json")
        print("1747 verified source-plane trace enclosure, protected zero line and coherent two-source witnesses")
    else:
        (HERE/"results.json").write_text(json.dumps(fine,indent=2,allow_nan=False)+"\n",encoding="utf8")
