"""Test measured charged-lepton axes on the two polar-normalized physical ports."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
PREV=runpy.run_path(str(HERE/CFG["old"]))
CHARGED=runpy.run_path(str(HERE/CFG["charged"]))["LEPTONS"]
OBS=PREV["PARENT"]["OBS"]
ANG=PREV["ANG"]
DATA=json.loads((HERE/CFG["observed"]).read_bytes())
KL,SK=DATA["independent_KamLAND"],DATA["independent_SK"]
R0=math.sqrt(1+SK[0]/KL[0])
RLO=math.sqrt(1+(SK[0]-SK[1])/(KL[0]+KL[2]))
RHI=math.sqrt(1+(SK[0]+SK[2])/(KL[0]-KL[1]))


def frame(A,weights,low=False):
    B=A
    if low:
        s,v=np.linalg.eigh(A.conj().T@A)
        B=A@((v/np.sqrt(s))@v.conj().T)
    g=B.conj().T@(weights[:,None]*B)
    eigen,O=np.linalg.eigh(g)
    return eigen,O


def upper(A):
    p=abs(A)**2
    p/=p.sum(axis=0)
    zbox=(1-np.sqrt(1-ANG[0]))/2
    plane=A@np.linalg.pinv(A)
    assert plane[0,0].real<1-zbox[1]
    out=[]
    for e,mu,tau in p.T:
        t=np.clip(mu/(mu+tau),*ANG[2])
        k=math.sqrt(mu*t)+math.sqrt(tau*(1-t))
        z=np.clip(e/(e+k*k),*zbox)
        out.append((math.sqrt(e*z)+k*math.sqrt(1-z))**2)
    return np.array(out)


def axis_table(A):
    measured=CHARGED[:,0]/CHARGED[-1,0]
    weights={"charged_yukawa_squared":measured**2,"charged_yukawa_magnitude":measured,
             "electron":np.array([1.,0.,0.]),"muon":np.array([0.,1.,0.]),
             "tau":np.array([0.,0.,1.]),"equal_response":np.ones(3)}
    result={}
    for label,w in weights.items():
        for low in [False,True]:
            vals,O=frame(A,w,low)
            name=label+("_low_polar" if low else "")
            if vals[1]-vals[0]<1e-13:
                result[name]=dict(eigen=vals,axis_undetermined=True)
                continue
            out=upper(A@O)
            result[name]=dict(eigen=vals,upper=out,
                              rejects_every_filter=bool(max(out)<RLO/(1+RLO)),
                              flavor_probabilities=abs(A@O)**2/np.linalg.norm(A@O,axis=0)**2)
    gap=np.linalg.norm(A[2,:])**2
    target=RLO/(1+RLO)
    margin=math.acos(math.sqrt(max(result["tau"]["upper"])))-math.acos(math.sqrt(target))
    d=math.sin(margin)
    cap=d*gap/(np.linalg.norm(A,2)**2*(np.linalg.cond(A)+d))
    return dict(tests=result,tau_axis_margin_degrees=math.degrees(margin),any_subleading_weight_below=cap,
                measured_mu_to_tau_squared=measured[1]**2,minimum_required_overlap=target)


def fixed_filter(M,gap,n):
    return ((M[1]**2+gap**2)/(M[0]**2+gap**2))**n


def circle(A,O,r,R):
    X=A@O@np.diag([r,1.])
    U,s,vh=np.linalg.svd(X,full_matrices=False)
    cap=(s[0]/s[1])**2
    if cap<R:return None
    c2=(R+1/R-2)/(cap+1/cap-2)
    det=np.prod(s)**2
    p=(s[0]**2*(s[1]**2+(s[0]**2-s[1]**2)*c2)-det/R)/(det*(R-1/R))
    v=math.sqrt(p)*U[:,0]; w=math.sqrt(1-p)*U[:,1]
    base=abs(v)**2+abs(w)**2
    harm=v.conj()*w
    ezero=1-np.sum(abs(U[0])**2)
    def equation(c,d):
        if abs(d)<1e-15:return []
        rhs=-c/(2*abs(d))
        if not -1<=rhs<=1:return []
        return [(k*math.acos(rhs)-np.angle(d))%(2*math.pi) for k in [-1,1]]
    def roots_atm(t):return equation(base[1]+t*base[0]-t,harm[1]+t*harm[0])
    def direct(phi,save=False):
        P=np.diag([1.,np.exp(1j*(phi-math.pi/2))])
        J0=np.array([[math.sqrt(c2),1j*math.sqrt(1-c2)],[1j*math.sqrt(1-c2),math.sqrt(c2)]])
        J=vh.conj().T@P@J0@P@vh.conj()
        assert np.linalg.norm(J@J.conj().T-np.eye(2))<5e-12
        # Propagating-frame basis. This must be O† J_source O*.
        M=X@J@X.T
        obs=OBS(M,np.zeros(3),0.)
        heavy,sv,_=np.linalg.svd(M)
        np.testing.assert_allclose(sv[0]/sv[1],R,rtol=2e-12)
        np.testing.assert_allclose(abs(heavy[:,0])**2,base+2*np.real(harm*np.exp(1j*phi)),atol=4.e-12)
        o=dict(angles=obs["angle_readings"],angle_windows=obs["angles_in_box"],delta=obs["delta_deg"],
               FC_68_90=obs["SK_T2K_FC68_90"],mbb_meV_KamLAND=obs["mbb_meV"]*math.sqrt(KL[0]/.00007388))
        if save:
            matrix=O@J@O.T
            o.update(J_real=matrix.real,J_imag=matrix.imag)
        return o
    extrema=[]
    # Derivative of the ratio (m0+2Re m1 z)/(1−e0−2Re e1 z).
    # For an exact attained atmospheric reading, discriminant is zero.
    bmu,be=base[1],base[0]-1
    hmu,he=harm[1],harm[0]
    for t in sorted(np.roots([be**2-4*abs(he)**2,2*bmu*be-8*(hmu*he.conjugate()).real,
                              bmu**2-4*abs(hmu)**2])):
        extrema.append(float(t.real))
    z=np.array([base[0]-2*abs(harm[0]),base[0]+2*abs(harm[0])])
    assert z[1]<.5
    return dict(z=z,sol=np.sort(1-ezero/(1-z)),atmosphere=np.array(extrema),
                roots_atm=roots_atm,direct=direct,reactor=4*z*(1-z),R=R)


def continuous_R(A,O,r):
    X=A@O@np.diag([r,1.])
    U,s,_=np.linalg.svd(X,full_matrices=False)
    C=(s[0]/s[1])**2
    hi=min(RHI,C)
    if hi<RLO:return None
    # The heaviest singular-frame weight is R/(R+1) + (R-1)/(R+1)/(C-1).
    def w(R):return (C*R-1)/((R+1)*(C-1))
    pol=np.array([w(RLO),w(hi)])
    row=abs(U[0])
    e1=(row[0]*np.sqrt(pol)-row[1]*np.sqrt(1-pol))**2
    e2=(row[0]*np.sqrt(pol)+row[1]*np.sqrt(1-pol))**2
    if pol[0]<=row[1]**2/(row@row)<=pol[1]: e1=np.append(e1,0.)
    if pol[0]<=row[0]**2/(row@row)<=pol[1]: e2=np.append(e2,row@row)
    z=np.array([min(e1),max(e2)])
    assert z[1]<.5
    # Invert the same 2D projection: for fixed electron fraction q in its plane,
    # vary the relative phase to bound the leading left-singular-frame probability.
    pe=row@row
    f=row[0]**2/pe
    measured_q=(1-np.sqrt(1-ANG[0]))/(2*pe)
    small=(np.sqrt(f*measured_q)-np.sqrt((1-f)*(1-measured_q)))**2
    large=(np.sqrt(f*measured_q)+np.sqrt((1-f)*(1-measured_q)))**2
    if measured_q[0]<=1-f<=measured_q[1]:small=np.append(small,0.)
    if measured_q[0]<=f<=measured_q[1]:large=np.append(large,1.)
    pband=np.array([min(small),max(large)])
    rband=(pband*(C-1)+1)/(C-pband*(C-1))
    check=circle(A,O,r,R0) if R0<=hi else None
    if check is not None:assert z[0]-1e-12<=check["z"][0]<=check["z"][1]<=z[1]+1e-12
    return dict(eligible_R=[RLO,hi],fixed_all_source_all_R_reactor=4*z*(1-z),
                R_allowed_by_reactor_only=[max(1.,rband[0]),min(C,rband[1])],
                intersection_reactor_and_clocks=[max(RLO,rband[0]),min(RHI,C,rband[1])])


def evaluate(u,step):
    _,material=PREV["PARENT"]["endpoint"](u,step)
    A=PREV["frame"](material)
    table=axis_table(A)
    masses=material["M"]
    options={}
    _,e=frame(A,np.array([1.,0.,0.]))
    _,equal=frame(A,np.ones(3))
    for role,O in [("electron_null_strong",e),("electron_visible_strong",e[:,::-1]),
                   ("equal_weak_response_strongest_eigenvalue",equal[:,::-1])]:
        for n,tag,g in [(2,"M1",masses[0]),(3,"M1",masses[0]),(3,"M2",masses[1]),
                         (4,"M2",masses[1]),(4,"Mweak",material["weak_scale"])]:
            r=fixed_filter(masses,g,n)
            key=f"{role}_{n}_{tag}"
            c=circle(A,O,r,R0)
            opts=dict(mediator_GeV=g,steps=n,ratio=r,mass_ratio_cap=PREV["hierarchy_cap"](A@O,r))
            if c is not None:
                opts.update(reactor_at_R0=c["reactor"],solar_at_R0=c["sol"],atm_at_R0=c["atmosphere"])
            options[key]=opts
    r=fixed_filter(masses,masses[1],4)
    c=circle(A,e,r,R0)
    direct=c["direct"]
    eq=c["roots_atm"](.5)
    examples=[direct(v,True) for v in eq]
    for ex in examples:np.testing.assert_allclose(ex["angles"][2],.5,atol=2e-12)
    weak_atm=c["roots_atm"](.57)
    return dict(ub=u,charged_yukawas_MZ=CHARGED[:,0],actual_Dirac_poles=masses,
                required_R=[RLO,RHI],central_R=R0,axes=table,independent_mass_options=options,
                electron_null_four_upper=dict(all_R=continuous_R(A,e,r),at_R_central={k:c[k] for k in
                ("z","sol","reactor","atmosphere")},atmospheric_central_057_roots=len(weak_atm),
                mu_tau_equal_modulus_all_roots=examples))


def tidy(v):
    if isinstance(v,dict):return {k:tidy(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [tidy(x) for x in v]
    if isinstance(v,np.ndarray):return tidy(v.tolist())
    if isinstance(v,(float,np.floating)):return float(f"{v:.10g}")
    if isinstance(v,np.bool_):return bool(v)
    return v


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    u1=PREV["PARENT"]["WEAK"]["inclusive_ub"][0]
    u2=DATA["combined_central_solution"]["Vub"]
    primary=evaluate(u1,.25)
    other=evaluate(u2,.25)
    fine=evaluate(u1,.125)
    diff=np.max(abs(fine["axes"]["tests"]["charged_yukawa_squared"]["upper"]-primary["axes"]["tests"]["charged_yukawa_squared"]["upper"]))
    result=tidy(dict(inclusive=primary,profiled=other,step_halving_overlap=diff))
    file=HERE/"results.json"
    if args.check:
        assert result==json.loads(file.read_bytes())
        print("1762 OK: all-strength charged Yukawa exclusion, flavor-tagged electron filter, finite mass and both CP roots.")
    else:
        file.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

