"""Independent spin traces, charged-current bound and SVD of the reported weak-basis witness."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def mixing_bound():
    coeff=read("../1464/inputs.json")["HPQCD_tableI_eq15"]
    p_belle=read("../1406/inputs.json")["BelleII"]
    cal=read("../1404/inputs.json")
    source=read("../1464/inputs.json")
    f,*df=source["fBs_GeV"]
    th,dth=cal["tau_ps"]["s_heavy"]
    const=cal["const_GeV"]
    mass,mu,GF=const["M_Bs"],const["m_mu"],const["GF"]
    hbar=read("../1407/inputs.json")["hbar_GeV_ps"]
    r=cal["LHCb_run3_2024_25"]["s_in_1e9"]
    hmu=(GF*mu)**2*mass*math.sqrt(1-(2*mu/mass)**2)/(16*math.pi*hbar)
    hmu*=(f-2*sum(df))**2*(th-2*dth)*(1-source["mu_rate_transfer_width"])
    p=math.sqrt((r[0]+2*(r[2]+r[4]))*1e-9/(hmu*source["run3_working_light_response"]))
    ref=coeff["lambda_ref"]*coeff["Xt_ref"]/(math.pi*coeff["alpha_inverse_ref"]*coeff["sine_squared_ref"])
    n=ref*math.sqrt((p_belle["mu"]-2*p_belle["width"])/(1+2*coeff["hadronic_rate_width"]))
    c,k=.017,.60
    C=math.hypot(.005,.046)
    a_max=.2*5/1000**2
    h=(1-k)/(1-c*k)
    # Replace (1-a) by 1 for a uniform conservative positive lower endpoint.
    any_left=max(0.,n-p-h*c*(C+p))/(2+h*c)
    aligned=max(0.,n-c*C)/(1+(2+c)*a_max)
    upper=(p+.70*c*C)/(1-(2+.70*c)*a_max)
    return n,p,any_left,aligned,upper


def spin_weights(M,masses,mw):
    I=np.eye(2,dtype=complex)
    zero=np.zeros((2,2),complex)
    pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    gamma=[np.block([[I,zero],[zero,-I]])]+[np.block([[zero,a],[-a,zero]]) for a in pauli]
    g5=1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]
    one=np.eye(4)
    pl,pr=(one-g5)/2,(one+g5)/2
    def trace_weight(x,y):
        lam=(1-x-y)**2-4*x*y
        d=math.sqrt(lam)/2
        qslash=gamma[0]*((1+x-y)/2)-d*gamma[3]
        phys_k=[(1+y-x)/2,0,0,d]
        pol=np.outer(phys_k,phys_k)/y-np.diag([1,-1,-1,-1])
        vertex=[g@(math.sqrt(x)*pl+pr) for g in gamma]
        pair=0j
        for i,j in np.ndindex(4,4):
            adj=gamma[0]@vertex[j].conj().T@gamma[0]
            pair+=pol[i,j]*np.trace((qslash+math.sqrt(x)*one)@vertex[i]@(gamma[0]+one)@adj)/2
        assert abs(pair.imag)<1.e-9
        return pair.real*2*y*math.sqrt(lam)
    for m in (1000,1600,1701,5400):
        y=(80.38/m)**2
        for f in (.1,4.8,5):
            x=(f/m)**2
            expected=math.sqrt((1-x-y)**2-4*x*y)*((1+x)*((1-x)**2+(1+x)*y-2*y*y)-12*x*y)
            assert math.isclose(trace_weight(x,y),expected,rel_tol=2e-14)
    return np.array([trace_weight((mq/M)**2,(mw/M)**2) if mq+mw<M else 0.
                     for mq in masses])


def uniform_spin_range():
    # A deliberately coarse interval valid for M>=1 TeV, 0<=m_q<=5 GeV.
    x,y=(5/1000)**2,(80.38/1000)**2
    floor=math.sqrt((1-x-y)**2-4*x*y)*((1-x)**2-2*y*y-12*x*y)
    ceiling=(1+x)*(1+(1+x)*y)
    return floor/ceiling,ceiling/floor


def read_complex(a):
    a=np.asarray(a,dtype=float)
    return a[...,0]+1j*a[...,1]


def mass_and_currents(w, cfg):
    """SVD begins with the weak-basis mass matrix only; phase fixing is conventional."""
    mass=read_complex(w["weak_basis_mass_GeV"])
    u,d,vh=np.linalg.svd(mass)
    ix=np.argsort(d)
    u,d,v=u[:,ix],d[ix],vh.conj().T[:,ix]
    np.testing.assert_allclose(d,w["physical_masses_GeV"],rtol=1e-10,atol=1e-10)
    hard=mass[4,4].real
    assert hard>0 and np.max(abs(mass[:4,4]))<1e-9 and abs(mass[4,3])<1e-9
    assert np.max(abs(mass[3,:3]))<1e-9
    rec_r=read_complex(w["right_D"])
    # An independent mass-eigenstate phase per column cancels out of all rates.
    refk=read_complex(w["light_charged_K"])
    phase=np.ones(5,dtype=complex)
    phase[0]=refk[0,0]/u[0,0]
    phase[1:]=rec_r[1:].conj()/v[-1,1:]
    phase/=abs(phase)
    u*=phase
    v*=phase
    np.testing.assert_allclose(u[:3,:3],refk,atol=3e-10)
    np.testing.assert_allclose(u[-1,:],d*v[-1,:]/hard,atol=2e-11)
    t=abs(v[-1,2])**2
    assert abs(t-w["Tbb"])<2e-10
    ps=u[-2,1].conj()*u[-2,2]
    assert abs(ps)<2e-11
    # Weak T3 operators in the Q1,Q2,Q3,S,D basis, normalized by -2.
    gl=u.conj().T@np.diag([1.,1.,1.,0.,-1.])@u
    gr=v.conj().T@np.diag([0.,0.,0.,0.,-1.])@v
    f,fl=gr[1,2],gl[1,2]
    a=float(d[1]*d[2]/hard**2)
    g=(u[:3].conj().T@u[:3])[1,2]
    for x,y in [(g,fl-a*f),(fl,2*a*f-ps)]:
        assert abs(x-y)<2e-13
    theta=read("../1464/inputs.json")["HPQCD_tableI_eq15"]
    c=theta["Xt_ref"]/(math.pi*theta["alpha_inverse_ref"]*theta["sine_squared_ref"])
    k=read("../1406/inputs.json")["sm_coefficient_ratio"]["center"]
    L=-c*u[2,1].conj()*u[2,2]
    N=L+fl+f
    P=k*L+fl-f
    np.testing.assert_allclose([fl,f,L,N,P],read_complex(w["F_L_R_and_L_N_P"]),atol=3e-11)
    wide=read("../1464/results.json")
    belle=read("../1406/inputs.json")["BelleII"]
    neutral=abs(N)**2*wide["Hnu"]/(belle["reference_SD_in_1e-6"]*1e-6)
    assert abs(neutral-w["Belle_short_distance_mu"])<1e-7
    z=P*L.conjugate()/abs(L)
    for name,kind in [("CMS","calculated_light_response"),("LHCb_R3","run3_central_light_response")]:
        signal=1e9*wide["Hmu_heavy"]*(z.real**2+wide[kind]*z.imag**2)
        assert abs(signal-w["fitted_muons_in_1e9"][name][0])<1e-7
    unit=cfg["witness"]
    mw=unit["W_mass_GeV"]
    weights=spin_weights(hard,d,mw)*np.abs(v[-1,:])**2
    widths=weights*(unit["weak_g"]**2*hard**3/(64*math.pi*mw**2))
    np.testing.assert_allclose(widths,w["Y_to_Wq_widths_GeV"],rtol=2.e-9,atol=1.e-12)
    ratio=weights[1]/weights[2]
    assert abs(ratio-w["width_ratio_sW_to_bW"])<1e-11
    higgs=mass.copy()
    higgs[3,3]=higgs[4,4]=0.
    yuk=math.sqrt(2*np.linalg.eigvalsh(higgs.conj().T@higgs)[-1])/unit["Higgs_vev_GeV"]
    assert abs(yuk-w["largest_yukawa_singular_value"])<1e-9
    return {"partner_mass_GeV":round(hard,8),"recovered_Tbb":round(t,10),
            "recovered_F_abs":round(abs(f),10),"width_ratio":round(ratio,12),
            "total_width_GeV":round(sum(widths),8),"max_yukawa":round(yuk,9)}


def run():
    cfg=read("inputs.json")
    results=read("results.json")
    n,p,arbitrary,aligned,upper=mixing_bound()
    box=results["rate_envelopes"]["LHCb_R3_2"]
    assert aligned<=box["F_abs_bounds"][0] and upper>=box["F_abs_bounds"][1]
    spin_floor,spin_ceil=uniform_spin_range()
    t_cap=.07
    assert box["Tbb"][1]<t_cap
    # Positive contraction on the light subspace, completable by heavy columns
    # to two orthonormal right rows. Row 1 has zero sW while the summed FCNC is nonzero.
    rows=np.array([[0.,.1],[.02,.15]])
    rest=np.eye(2)-rows@rows.T
    rows=np.column_stack((rows,np.linalg.cholesky(rest)))
    assert np.allclose(rows@rows.T,np.eye(2))
    assert abs((rows.T@rows)[0,1])>0 and rows[0,0]==0.
    witnesses={k:mass_and_currents(v,cfg) for k,v in results.items()
               if isinstance(v,dict) and "weak_basis_mass_GeV" in v}
    assert set(witnesses)=={"common_witness","high_mass_witness"}
    return {"two_width_abs_N_lower":round(n,12),"two_width_abs_P_upper":round(p,12),
            "uniform_F_right_lower_arbitrary_left":round(arbitrary,12),
            "uniform_F_right_bounds_aligned":[round(aligned,12),round(upper,12)],
            "spin_ratio_safe_interval":[round(spin_floor,8),round(spin_ceil,8)],
            "width_ratio_lower_arbitrary_left":float(f"{spin_floor*(arbitrary/t_cap)**2:.9g}"),
            "direct_SVD_checks":witnesses,"two_doublet_rank_counterexample":True}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    options=parser.parse_args()
    result=run()
    if options.check:
        assert result==read("independent_results.json")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as file:
            file.write(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("1635 independent SVD, weak generators, Dirac trace and uniform bounds PASS")
