"""Reconstruct the LNV spurion on the two existing Dirac ports at a stated common matching scale."""
import argparse
import json
import math
import runpy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial import Polynomial


HERE=Path(__file__).resolve().parent
def read(path):return json.loads((HERE/path).read_bytes())
CFG=read("inputs.json")
WEAK=read(CFG["strong_and_B"])
PROFILE=read(CFG["alt_frozen"])["combined_central_solution"]
OBS=runpy.run_path(str(HERE/"../1731/check.py"))["observables"]


def endpoint(ub,step):
    model=runpy.run_path(str(HERE/CFG["pipeline"]))
    found={}
    def load_module(path):
        mod=runpy.run_path(str(path))
        if Path(path).parent.name!="1748":
            return mod
        daughter=mod["daughter"]
        def with_ports(bg,F):
            # Save the actual last downward step and final charged phase frame of the tau shooting.
            matter=bg["model"]
            advance=matter["advance"]
            def traced(y,T,dt,g,flow):
                out=advance(y,T,dt,g,flow)
                if dt<0:
                    found["last"]=out
                return out
            matter["advance"]=traced
            try:
                gain,row=daughter(bg,F)
            finally:
                matter["advance"]=advance
            y,T=found.pop("last")
            pt=bg["pt"]
            flow=pt["flow"]
            final_basis=flow["extraction"](y)[3]
            phases=np.diag(final_basis.conj().T@pt["low_E"])
            final_basis=final_basis@np.diag(phases/abs(phases))
            transfer=F@np.linalg.solve(T,final_basis)
            carrier=np.linalg.inv(bg["S"]@pt["low_E"])
            C=carrier.T@pt["C"]@carrier
            b=carrier.T@pt["b"]
            M=np.array(flow["fields"]["neutral_poles_GeV"])[[0,2]]
            Y=bg["N"].conj()
            B=CFG["matter_vev_GeV"]*Y/M
            L=transfer.T@B
            C0=transfer.T@C@transfer
            b0=transfer.T@b
            np.testing.assert_allclose(OBS(C0,b0,bg["nu"]["FROZEN"][0])["angle_readings"],
                                       row["neutrinos"]["angle_readings"],rtol=2.e-10)
            found.update(C=C,b=b,Y=Y,M=M,B=B,L=L,C0=C0,b0=b0,
                         weak_scale=bg["Mweak"],R=bg["nu"]["FROZEN"][0],tau_gain=gain)
            return gain,row
        mod["daughter"]=with_ports
        return mod
    model["point"].__globals__["runpy"]=SimpleNamespace(run_path=load_module)
    row=model["shoot"](HERE,1/WEAK["QCDSF_s_d"][0],ub,WEAK["inclusive_cb"][0],step)
    return row,found


def pack(z):
    a=np.asarray(z)
    return dict(real=a.real,imag=a.imag)


def diagonal_scan(L):
    """For universal |mu_11|=|mu_22| let even an arbitrary relative phase vary."""
    x,y=L.T
    A,B=np.outer(x,x),np.outer(y,y)
    t0=np.sum(abs(A)**2)+np.sum(abs(B)**2)
    cross=abs(np.vdot(A,B))
    area=np.linalg.det(L.conj().T@L).real
    extrema=(np.array([t0-2*cross,t0+2*cross])/area)
    assert min(extrema)>2
    R=(extrema+np.sqrt(extrema**2-4))/2
    # Check limiting phases in the original three-dimensional basis.
    ang=-np.angle(np.vdot(A,B))
    for a,exp in zip([ang+np.pi,ang],R):
        sv=np.linalg.svd(A+np.exp(1j*a)*B,compute_uv=False)
        np.testing.assert_allclose(sv[0]/sv[1],exp,rtol=5e-11)
    return R


def pole_splits(B,C,b,C0,b0,rho):
    dual=np.linalg.pinv(B)
    Q=dual@C@dual.T
    q=dual@b
    np.testing.assert_allclose(B@Q@B.T,C,rtol=1e-11,atol=1e-13)
    np.testing.assert_allclose(B@q,b,rtol=1e-11,atol=1e-13)
    arm=np.outer(q,q)
    ab=Q[0,1]; cd=arm[0,1]
    best=-float((ab*cd.conjugate()).real/abs(cd)**2)
    val=Q+rho*arm
    scale=OBS(C0,b0,rho)["scale_eV"]
    ratio=abs(val[0,1])/math.sqrt(abs(val[0,0]*val[1,1]))
    return dict(base=Q,arm=arm,source=q,best=best,resid=abs(ab+best*cd),
                fraction=abs(ab+best*cd)/abs(ab),complex_root=-ab/cd,
                mu_eV=scale*val,scale_eV=scale,cross_to_geometric_diagonal=ratio)


def diagonal_clock_reactor(L, ratio, a):
    """All complex diagonal ratios at fixed nonzero mass ratio and reactor reading, up to global scale/phase."""
    # Normalize columns only to improve the conditioning; recover the actual heavy coefficients below.
    norms=np.linalg.norm(L,axis=0)
    U=L/norms
    A=np.outer(U[:,0],U[:,0]); B=np.outer(U[:,1],U[:,1])
    E=np.array([1.,0.,0.]); n=np.linalg.svd(L.T)[2].conj().T[:,-1]
    pe=1-abs(n[0])**2
    h=np.vdot(A,B); e=np.vdot(A@E,B@E)
    det=np.linalg.det(U.conj().T@U).real
    small=2*a/(4*(1+math.sqrt(1-a)))
    out=[]
    for p in (small,1-small):
        if p>pe:continue
        # For t=r exp(i phi), the trace and <e|H|e> equations are linear in Re(t),Im(t).
        K=np.array([[2*h.real,-2*h.imag],[2*e.real,-2*e.imag]])
        assert abs(np.linalg.det(K))>1.e-5
        target_trace=det*(ratio+1/ratio)
        target_e=det*((ratio-1/ratio)*p+pe/ratio)
        rhs=np.array([[-1,target_trace,-1],[-abs(U[0,0])**2,target_e,-abs(U[0,1])**2]])
        real,imag=[Polynomial(row) for row in np.linalg.solve(K,rhs)]
        eq=real**2+imag**2-Polynomial([0,0,1])
        roots=eq.roots()
        realcount=runpy.run_path(str(HERE/"../1757/independent_check.py"))["sturm_real_count"](eq.coef)
        assert sum(abs(v.imag)<1.e-8 for v in roots)==realcount
        for r in roots:
            if abs(r.imag)>1.e-8 or r.real<=0:continue
            t=complex(real(r.real),imag(r.real))
            mat=A+t*B
            obs=OBS(mat,np.zeros(3),0.)
            np.testing.assert_allclose(obs["angle_readings"][0],a,atol=1e-10)
            sv=np.linalg.svd(mat,compute_uv=False)
            np.testing.assert_allclose(sv[0]/sv[1],ratio,rtol=1e-10)
            out.append(dict(physical_mu22_over_mu11=pack(t*(norms[0]/norms[1])**2),
                            angles=obs["angle_readings"],angles_in_box=obs["angles_in_box"],
                            m_bb_meV_KamLAND_scale=obs["mbb_meV"]*math.sqrt(.0000754/.00007388),
                            delta_deg=obs["delta_deg"],SK_T2K_68_90=obs["SK_T2K_FC68_90"]))
    return out


def window_test(L,Y,M):
    x,y=L.T
    a,b=np.outer(x,x),np.outer(y,y)
    n=np.linalg.svd(L.T)[2].conj().T[:,-1]
    p0=1-abs(n[0])**2
    area=np.linalg.det(L.conj().T@L).real
    lo,hi=math.sqrt(1+.00219/.0000773),math.sqrt(1+.00263/.0000736)
    pmax=(1-math.sqrt(1-.0899))/2
    assert p0<1-pmax
    def necessary(A,B,upper):
        return np.sort(np.roots([np.vdot(B,B).real,-(area*upper+2*abs(np.vdot(A,B))),np.vdot(A,A).real]))
    electron=necessary(a[:,0],b[:,0],(hi-1/hi)*pmax+p0/lo)
    trace=necessary(a,b,hi+1/hi)
    e=np.linalg.norm(a[:,0])**2+np.linalg.norm(b[:,0])**2
    de=2*abs(np.vdot(a[:,0],b[:,0]))
    eunit=np.array([((e-de)/area-p0/lo)/(hi-1/hi),((e+de)/area-p0/hi)/(lo-1/lo)])
    assert 0<eunit[0]<eunit[1]<.5
    ranges=np.array([[.0803,.0899],[.2908,.3164],[.49,.63]])
    def arc(k):
        bins=4096
        half=np.pi/bins
        eps=2*k*np.linalg.norm(b,2)*math.sin(half/2)*(1+1e-12)
        t=(np.arange(bins)+.5)*(2*half)-np.pi
        _,s,vh=np.linalg.svd(a+k*np.exp(1j*t[:,None,None])*b)
        gap=s[:,0]**2-s[:,1]**2
        dh=(2*s[:,0]+eps)*eps
        err=np.ones(bins)
        stable=gap>dh
        err[stable]=np.minimum(1,dh[stable]/(gap[stable]-dh[stable])+1e-11)
        # All evaluated two-source models have nondegenerate nonzero eigenvalues.
        assert min(s[:,1])>eps
        lok=(s[:,0]-eps)/(s[:,1]+eps)-1e-11
        hik=(s[:,0]+eps)/(s[:,1]-eps)+1e-11
        p=abs(vh[:,0,:])**2
        pmin=np.maximum(0,p[:,0]-err)
        pmax=np.minimum(p0,p[:,0]+err)
        # p0<1/2, so the reactor map is increasing throughout the physical plane.
        assert p0<.5
        low=np.column_stack([4*pmin*(1-pmin),(p0-pmax)/(1-pmax),np.maximum(0,p[:,1]-err)/(1-pmin)])-1e-11
        high=np.column_stack([4*pmax*(1-pmax),(p0-pmin)/(1-pmin),np.minimum(1,p[:,1]+err)/(1-pmax)])+1e-11
        freq=(hik>=lo)&(lok<=hi)
        passed=(low<=ranges[:,1])&(high>=ranges[:,0])
        central=np.column_stack([4*p[:,0]*(1-p[:,0]),(p0-p[:,0])/(1-p[:,0]),p[:,1]/(1-p[:,0])])
        ok=(s[:,0]/s[:,1]>=lo)&(s[:,0]/s[:,1]<=hi)&np.all((central>=ranges[:,0])&(central<=ranges[:,1]),axis=1)
        centers=np.flatnonzero(ok)
        closure=None
        if len(centers):
            i=centers[0]
            ref=OBS(a+k*np.exp(1j*t[i])*b,np.zeros(3),0.)
            closure=dict(relative_phase_degrees=np.degrees(t[i]),
                         rephasing_invariant_degrees=np.degrees(np.angle(np.exp(1j*t[i])*np.vdot(*Y.T)**2)),
                         angles=central[i],mass_ratio=s[i,0]/s[i,1],delta_deg=ref["delta_deg"],
                         SK_T2K_68_90=ref["SK_T2K_FC68_90"])
        return dict(k=k,cells=bins,possible_clock_and_reactor=int(sum(freq&passed[:,0])),
                    possible_clock_and_all_angles=int(sum(freq&np.all(passed,axis=1))),
                    witness=closure)
    yn=np.linalg.norm(Y,axis=0)
    candidates=dict(equal=1.,common_fractional_splitting=M[1]/M[0],inverse_mass=M[0]/M[1],
                    two_inverse_mass_insertions=(M[0]/M[1])**2,equal_LNV_times_Yukawa_weight=(yn[0]/yn[1])**2,
                    equal_LNV_times_Higgs_weight=(yn[0]*M[0]/(yn[1]*M[1]))**2,
                    equal_individual_light_mass_norm=(np.linalg.norm(x)/np.linalg.norm(y))**2)
    return dict(trace_bound=trace,electron_bound=electron,combined_bound=[max(trace[0],electron[0]),min(trace[1],electron[1])],
                p0=p0,universal_electron=eunit,universal_reactor_lower=4*eunit[0]*(1-eunit[0]),
                tied_models={name:arc(k) for name,k in candidates.items()})


def normalized_source(d):
    B,L=d["B"],d["L"]
    g=B.conj().T@B
    eig,vec=np.linalg.eigh(g)
    inv=(vec/np.sqrt(eig))@vec.conj().T
    orthogonal=B@inv
    carried=L@inv
    np.testing.assert_allclose(orthogonal.conj().T@orthogonal,np.eye(2),atol=2.e-14)
    s=np.linalg.svd(carried,compute_uv=False)
    distortion=(s[0]/s[1])**2
    lower=math.sqrt(1+.00219/.0000773)/distortion
    # Explicit saturating symmetric unitary; no orientation scan is needed.
    _,_,vh=np.linalg.svd(carried,full_matrices=False)
    polar=vh.conj().T@vh.conj()
    np.testing.assert_allclose(polar@polar.conj().T,np.eye(2),atol=1.e-14)
    sv=np.linalg.svd(carried@polar@carried.T,compute_uv=False)
    np.testing.assert_allclose(sv[0]/sv[1],distortion,rtol=3.e-13)
    phase=np.vdot(*d["Y"].T)
    rephase=np.diag([1.,phase.conjugate()/abs(phase)])
    Yn=d["Y"]@rephase
    G=Yn.conj().T@Yn
    np.testing.assert_allclose(G.imag,0.,atol=1e-15)
    G=G.real
    W=L@rephase
    examples={}
    for label,J in [("common_identity",np.eye(2)),("shared_inverse_Y_Gram",np.linalg.inv(G)),
                    ("real_diagonal_inverse_Y_norms",np.diag(1/np.diag(G)))]:
        M=W@J@W.T
        o=OBS(M,np.zeros(3),0.)
        ss=np.linalg.svd(M,compute_uv=False)
        examples[label]=dict(mass_ratio=ss[0]/ss[1],angles=o["angle_readings"],delta=o["delta_deg"])
    return dict(common_scale_B_singular_values=np.sqrt(eig)[::-1],low_transport_of_unit_plane_singular_values=s,
                all_symmetric_unitary_sources_mass_ratio_upper=distortion,
                dimensionless_source_condition_number_minimum=lower,explicit_CP_frame_tests=examples)


def evaluate(u,step):
    row,d=endpoint(u,step)
    Q=pole_splits(d["B"],d["C"],d["b"],d["C0"],d["b0"],d["R"])
    rk=math.sqrt(1+.0025/.0000754)
    obs=diagonal_clock_reactor(d["L"],rk,.0851)
    return dict(us=row["low_moduli"]["Vus"],ub=u,M_dirac_GeV=d["M"],weak_match_GeV=d["weak_scale"],
                yukawa_norms=np.linalg.norm(d["Y"],axis=0),polar_yukawa_overlap=abs(np.vdot(*d["Y"].T))/np.prod(np.linalg.norm(d["Y"],axis=0)),
                B_condition=np.linalg.cond(d["B"]),low_port_norms=np.linalg.norm(d["L"],axis=0),
                frozen_rho=d["R"],eV_per_shape_unit=Q["scale_eV"],mu_eV=pack(Q["mu_eV"]),
                base_heavy=pack(Q["base"]),rank_one_heavy=pack(Q["arm"]),
                real_rho_best_offdiagonal=Q["best"],offdiag_min_fraction_of_base_all_real_rho=Q["fraction"],
                complex_rho_for_exactly_diagonal=pack(Q["complex_root"]),cross_to_geometric_diagonal=Q["cross_to_geometric_diagonal"],
                equal_abs_majorana_all_phase_mass_ratio=diagonal_scan(d["L"]),
                constraints=window_test(d["L"],d["Y"],d["M"]),
                normalized_two_source=normalized_source(d),
                arbitrary_complex_diagonal_roots_for_two_clocks_and_reactor=obs)


def tidy(obj):
    if isinstance(obj,dict):return {k:tidy(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):return [tidy(v) for v in obj]
    if isinstance(obj,np.ndarray):return tidy(obj.tolist())
    if isinstance(obj,(bool,np.bool_)):return bool(obj)
    if isinstance(obj,(float,np.floating)):return float(f"{obj:.10g}")
    return obj


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    ub=WEAK["inclusive_ub"][0]
    primary=evaluate(ub,.25)
    compare=evaluate(PROFILE["Vub"],.25)
    fine=evaluate(ub,.125)
    results=tidy(dict(primary_inclusive=primary,profiled_comparison=compare,
                      step_halving_mu_eV=np.max(abs(np.array(fine["mu_eV"]["real"])-primary["mu_eV"]["real"])),
                      independent_frequency_mass_ratio_one_width=[math.sqrt(1+.00219/.0000773),math.sqrt(1+.00263/.0000736)]))
    dest=HERE/"results.json"
    if args.check:
        assert results==json.loads(dest.read_bytes())
        print("1758 OK: normalized two-port inversion, all real rho, universal diagonal phase circle, diagonal mass-and-reactor roots.")
    else:
        dest.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf8")

