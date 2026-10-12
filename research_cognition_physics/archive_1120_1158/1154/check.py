"""1154: cumulative force / spectral-integral tool calibrations, not a field simulation."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent

def opnorm(x):
    return float(np.linalg.norm(x, 2))

def exp_h(h, t=1.0):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T

def rational(x):
    return dict(numerator=x.numerator, denominator=x.denominator, float=float(x))

def parameter_recipe():
    # Exact g=eta=1 constants, independent of any uncomputed actual material spectrum.
    D=F(15); R0=F(1,60); MF=2*R0*D
    k2=8*D*D+128*D
    k4=8*MF/R0**4
    kt=8*(10+F(5,2))*D*D
    assert (MF,k2,k4,kt)==(F(1,2),F(3720),F(51840000),F(22500))
    epsilon=F(1,1000)
    rows=[]
    for g,C,nu,eta in [(1., .001, .03, 1.), (.2, .5, .0001, 5.),
                       (.01, 1., .2, 32.), (1., 1e-16, .001, 1.)]:
        Df=10/math.sqrt(g)+5/g
        R=1/(4*Df)
        old=min(R/1e4, math.sqrt(C*R**4/(1e6*g)),
                R*math.sqrt(nu/(1e6*g)))
        M=2*g*R*eta*Df
        K2=8*Df**2+128*eta*Df/g
        K4=8*M/R**4
        Kt=8*eta*(10*math.sqrt(g)+2.5)*Df**2
        lam=min(old, math.sqrt(float(epsilon)/(3*K2)),
                math.sqrt(float(epsilon)*C/(3*K4)),
                float(epsilon)*C/(3*Kt))
        budgets=[K2*lam*lam, K4*lam*lam/C, Kt*lam/C]
        assert 0 < lam <= old and lam <= R/1e4
        assert max(budgets) <= float(epsilon/3)*(1+2e-15)
        assert sum(budgets) <= float(epsilon)*(1+2e-15)
        T=math.pi/(C*lam*lam)
        assert math.isfinite(T)
        rows.append(dict(g=g,C=C,nu=nu,eta=eta,
                         D=Df,R0=R,MF=M,K2=K2,K4=K4,Kt=Kt,
                         previous_coupling=old,new_coupling=lam,
                         new_duration=T,
                         dressing_budget=budgets[0],
                         fourth_order_budget=budgets[1],
                         tail_budget=budgets[2],
                         total_budget=sum(budgets),
                         old_linear_tail_budget=Kt*old/C,
                         new_to_old_coupling=lam/old))
    assert rows[-1]["old_linear_tail_budget"] > 1e5
    return dict(name="one_shot_parameter_recipe_and_budget",passed=True,
                exact_g_eta_one=dict(MF=rational(MF),K2=rational(k2),
                                    K4=rational(k4),Kt=rational(kt)),
                epsilon=rational(epsilon), each_budget_cap=rational(epsilon/3),
                sample_inputs_not_actual_spectral_constants=rows,
                old_lambda_over_C_not_uniformly_controlled=True,
                changed_coupling_and_duration_explicit=True)

def radial_profile(n, radius=.5):
    nodes,weights=np.polynomial.legendre.leggauss(n)
    r=(nodes+1)*radius/2
    wr=weights*radius/2
    q=np.exp(-1/(1-(r/radius)**2))
    q/=4*math.pi*np.dot(wr,r*r*q)
    return r,wr,q

def compact_yukawa():
    def evaluate(n,d):
        r,wr,q=radial_profile(n)
        u,wu=np.polynomial.legendre.leggauss(n)
        # First angular integration is exact, as in 1149:
        # integral_-1^1 e^{-sqrt(s^2+r^2-2sru)}/sqrt(...) du
        # = (e^{-|s-r|}-e^{-(s+r)})/(s*r).
        s=np.sqrt(d*d+r[:,None]**2+2*d*r[:,None]*u[None,:])
        sr=s[:,:,None]; ra=r[None,None,:]
        angular=(np.exp(-np.abs(sr-ra))-np.exp(-(sr+ra)))/(sr*ra)
        phiA=.5*np.einsum("bur,r->bu",angular,wr*r*r*q)
        direct=2*math.pi*np.dot(wr*r*r*q,phiA@wu)
        Q=4*math.pi*np.dot(wr,r*q*np.sinh(r))
        exact=Q*Q*math.exp(-d)/(4*math.pi*d)
        return direct,exact,Q
    d=3.
    direct,exact,Q=evaluate(80,d)
    coarse,_,_=evaluate(48,d)
    step=2e-4
    derivative=(evaluate(80,d+step)[0]-evaluate(80,d-step)[0])/(2*step)
    exact_derivative=-(1+1/d)*exact
    impulse=math.pi*derivative/direct
    expected=-4*math.pi/3
    assert abs(direct-exact) < 1e-12
    assert abs(coarse-direct) < 1e-10
    assert abs(derivative-exact_derivative) < 1e-10
    assert abs(impulse-expected) < 1e-7
    return dict(name="compact_radial_yukawa_and_force_ratio",passed=True,
                calibration_profile_radius=.5,center_separation=d,
                radial_charge=Q,direct_integral=direct,exterior_formula=exact,
                quadrature_difference=abs(coarse-direct),
                derivative_step=step,finite_difference_derivative=derivative,
                analytic_derivative=exact_derivative,
                impulse_ratio=impulse,expected_impulse_ratio=expected,
                same_actual_material_profile_not_claimed=True,
                exact_angular_formula_reused_from_1149=True)

def spectral_integral(e,F_eig,T):
    # Stable divided difference for integral exp(i (E_i-E_j)t) dt.
    half=.5*T*(e[:,None]-e[None,:])
    kernel=T*np.exp(1j*half)*np.sinc(half/math.pi)
    return F_eig*kernel

def force_calibration():
    I4=np.eye(4)
    H0=np.diag(np.repeat([0.,1.,2.],4))
    W0=np.eye(12,4)
    P0=W0@W0.T
    parity=np.diag(np.repeat([1.,-1.,1.],4))
    X=np.array([[0.,1.],[1.,0.]])
    Z=np.diag([1.,-1.])
    rot=exp_h(.31*np.kron(X,X)+.23*np.kron(Z,X)+.19*np.kron(X,Z))
    B=rot@np.diag([.2,.4,.6,.8])
    d=3.
    # Amplitude calibration Q=1; this is not the actual 1153 profile.
    C=math.exp(-d)/(4*math.pi*d)
    Cp=-(1+1/d)*C
    a=math.sqrt(C)
    ap=Cp/(2*a)
    V=np.zeros((12,12),complex); Vp=V.copy()
    V[:4,4:8]=a*I4; V[4:8,:4]=a*I4
    V[4:8,8:]=B; V[8:,4:8]=B.conj().T
    Vp[:4,4:8]=ap*I4; Vp[4:8,:4]=ap*I4
    assert opnorm(parity@V@parity+V)<1e-14
    assert opnorm(V)<2.5 and opnorm(Vp)<2.5
    target=(math.pi*Cp/C)*I4
    K4=C*C*I4-(C/2)*(B@B.conj().T)
    assert opnorm(K4-np.diag(np.diag(K4)))>1e-5
    rows=[]
    R0=1/60
    for fraction in [.125, .2, .25]:
        lam=fraction*R0
        assert lam<=R0/4
        H=H0+lam*V
        force=-lam*Vp
        T=math.pi/(C*lam*lam)
        e,v=np.linalg.eigh(H)
        F_eig=v.conj().T@force@v
        integrated=spectral_integral(e,F_eig,T)
        compress=v.conj().T@W0
        full=compress.conj().T@integrated@compress
        pieces=[]
        for first,last in [(slice(0,4),slice(0,4)),
                           (slice(0,4),slice(4,12)),
                           (slice(4,12),slice(0,4)),
                           (slice(4,12),slice(4,12))]:
            pieces.append(compress[first].conj().T@
                          integrated[first,last]@compress[last])
        pp,pq,qp,qq=pieces
        cross=pq+qp
        assert opnorm(full-(pp+cross+qq))<1e-12
        assert opnorm(full-full.conj().T)<1e-12
        P=v[:,:4]@v[:,:4].conj().T
        delta=P-P0
        ev,basis=np.linalg.eigh(np.eye(12)-delta@delta)
        invroot=(basis*(1/np.sqrt(ev)))@basis.conj().T
        U=(P@P0+(np.eye(12)-P)@(np.eye(12)-P0))@invroot
        f_low=(U.conj().T@force@U)[:4,:4]
        p=opnorm(compress[4:])
        gap=float(e[4]-e[3])
        pp_bound=T*opnorm(f_low-lam*lam*Cp*I4)+abs(math.pi*Cp/C)*p*p
        cross_bound=4*p*opnorm(F_eig[:4,4:])/gap
        qq_bound=T*p*p*opnorm(F_eig[4:,4:])
        error=opnorm(full-target)
        bound=pp_bound+cross_bound+qq_bound
        theorem_bound=3720*lam*lam+51840000*lam*lam/C+22500*lam/C
        assert error<=theorem_bound
        assert opnorm(pp-target)<=pp_bound+2e-11
        assert opnorm(cross)<=cross_bound+1e-12
        assert opnorm(qq)<=qq_bound+1e-12
        assert error<=bound+2e-11
        difference=opnorm(full-pp)
        assert difference>1e-10
        assert opnorm(compress[4:])>1e-5
        ref=opnorm(np.kron(full-target,np.eye(3)))
        assert abs(ref-error)<1e-12
        # Exact integral consistency across two windows, without resetting any component.
        integral_half=spectral_integral(e,F_eig,T/2)
        transport=np.exp(1j*e*T/2)
        stitched=integral_half+transport[:,None]*integral_half*transport.conj()[None,:]
        assert opnorm(compress.conj().T@(stitched-integrated)@compress)<1e-10
        rows.append(dict(coupling=lam,lambda_fraction_of_R0=fraction,duration=T,gap=gap,
                         conservative_force_theorem_bound=theorem_bound,
                         bare_tail_norm=p,complete_impulse_eigenvalues=np.linalg.eigvalsh(full).tolist(),
                         full_impulse_error=error,full_spectral_bound=bound,
                         pp_error=opnorm(pp-target),pp_bound=pp_bound,
                         cross_norm=opnorm(cross),cross_endpoint_bound=cross_bound,
                         qq_norm=opnorm(qq),qq_tail_bound=qq_bound,
                         complete_minus_P_only=difference,
                         passive_reference_norm_difference=abs(ref-error),
                         two_window_integral_residual=opnorm(compress.conj().T@(stitched-integrated)@compress)))
    return dict(name="complete_rank_four_force_integral_and_PQ_terms",passed=True,
                hilbert_dimension=12,low_cluster_rank=4,
                bare_input_includes_tail=True,KG_parity_proxy=True,
                nontrivial_fourth_order_spin_matrix_offdiag_norm=opnorm(K4-np.diag(np.diag(K4))),
                Q_calibration=1.,C=C,C_derivative=Cp,
                leading_impulse_scalar=math.pi*Cp/C,
                samples=rows,
                diagnostic_couplings_not_actual_lambda_star=True,
                diagnostic_sharp_bound_uses_exact_spectral_decomposition=True,
                conservative_theorem_bound_also_checked_within_R0_over_4=True,
                diagnostic_is_not_physical_epsilon_certificate=True,
                field_model_simulation=False,field_reset_performed=False)

def compare(a,b,trail="root"):
    if isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=2e-8,abs_tol=1e-10),(trail,a,b)
    elif isinstance(a,dict):
        assert a.keys()==b.keys(),trail
        for k in a: compare(a[k],b[k],trail+"."+k)
    elif isinstance(a,list):
        assert len(a)==len(b),trail
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,trail+"."+str(i))
    else:assert a==b,(trail,a,b)

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    sources=HERE/"source_hashes.json"
    if sources.exists():
        for name,digest in json.loads(sources.read_text(encoding="utf-8")).items():
            assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    result=dict(round=1154,all_passed=True,numerical_check_groups=3,
                groups=[parameter_recipe(),compact_yukawa(),force_calibration()],
                interpretation="Analytic-tool calibrations; no complete field simulation or actual spectral parameter estimate.",
                physical_new_experiments=0,new_adopted_axioms=0,full_goal_complete=False)
    dest=HERE/"results.json"
    if args.write_results and not dest.exists():
        dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        assert dest.exists(),"Run --write-results once."
        compare(result,json.loads(dest.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="PASS",round=1154,groups=3,
                          source_hashes_checked=sources.exists(),
                          force_ratio=result["groups"][1]["impulse_ratio"],
                          max_matrix_error=max(x["full_impulse_error"] for x in result["groups"][2]["samples"]))))
if __name__=="__main__":main()
