"""871: full conditional history weight of the original870 finite preparation.
Exact finite CAR calculations retain exterior modes and all pulse ordering.
No infinite-dimensional Fredholm determinant or finite-coupling QG is claimed.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'870'))
import receiver_four_point_backreaction as old
TARGET=HERE/'receiver_history_influence_results.json'
I=old.I;X=old.X;Y=old.Y

def gamma(U):
    n=len(U);basis=[[i for i in range(n) if z&(1<<(n-1-i))] for z in range(2**n)]
    out=np.zeros((2**n,2**n),complex)
    for i,a in enumerate(basis):
        for j,b in enumerate(basis):
            if len(a)==len(b):out[i,j]=np.linalg.det(U[np.ix_(a,b)]) if a else 1.
    return out

def exp_derivative(h,dh,t):
    ev,v=np.linalg.eigh(h);ex=np.exp(-1j*t*ev)
    divided=np.zeros((len(ev),len(ev)),complex)
    for i in range(len(ev)):
        for j in range(len(ev)):
            divided[i,j]=(ex[i]-ex[j])/(ev[i]-ev[j]) if abs(ev[i]-ev[j])>1e-12 else -1j*t*ex[i]
    return (v*ex)@v.conj().T, v@(divided*(v.conj().T@dh@v))@v.conj().T

def cofactor_derivative(A,dA):
    return sum((-1)**(i+j)*np.linalg.det(np.delete(np.delete(A,i,axis=0),j,axis=1))*dA[i,j]
               for i in range(len(A)) for j in range(len(A)))

def weights(r,m):return np.array([(1-r)/2,(1-r)/2,(r+m)/2,(r-m)/2])

def branch_covariances(Cperp):
    n=len(Cperp)+2;out=[]
    for Cq in (np.zeros((2,2)),I,(I+X)/2,(I-X)/2):
        C=np.zeros((n,n),complex);C[:2,:2]=Cq;C[2:,2:]=Cperp;out.append(C)
    return out

def four_det(U,Cperp,r,m,dU=None):
    ans=0j;der=0j
    for w,C in zip(weights(r,m),branch_covariances(Cperp)):
        M=np.eye(len(U))-C+C@U
        ans+=w*np.linalg.det(M)
        if dU is not None:der+=w*cofactor_derivative(M,C@dU)
    return ans,der

def schur(U,C,r,m):
    D=np.eye(len(C))-C+C@U[2:,2:]
    A=U[:2,:2]-U[:2,2:]@np.linalg.solve(D,C@U[2:,:2])
    z=np.linalg.det(D)
    out=z*((1-r)*(1+np.linalg.det(A))/2+r*np.trace(A)/2+m*np.trace(X@A)/2)
    return out,z,A

def finite_case(n):
    rng=np.random.default_rng(871+n)
    def herm():
        raw=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        return (raw+raw.conj().T)/(2*np.sqrt(n))
    H1,H2,V,W=[herm() for _ in range(4)]
    def evolution(s):
        h1=H1+s*V+.3*s*s*W;h2=H2+.4*s*W+.2*s*s*V
        u1,d1=exp_derivative(h1,V+.6*s*W,.31)
        u2,d2=exp_derivative(h2,.4*W+.4*s*V,.47)
        return u2@u1,d2@u1+u2@d1
    raw=rng.normal(size=(n-2,n-2))+1j*rng.normal(size=(n-2,n-2))
    rot=np.linalg.qr(raw)[0];occupations=np.linspace(.17,.74,n-2)
    Cperp=(rot*occupations)@rot.conj().T
    probs=np.array([np.prod([occupations[i] if z&(1<<(n-3-i)) else 1-occupations[i] for i in range(n-2)]) for z in range(2**(n-2))])
    grot=gamma(rot);rho_ext=(grot*probs)@grot.conj().T
    m=.5;rs=(.5,.625,1.)
    modes=old.annihilators(n)
    def second(A):
        return sum(A[i,j]*(modes[i].conj().T@modes[j]) for i in range(n) for j in range(n))
    um,_=evolution(-.21)
    errors=[];source_errors=[];wrong=[];cov_errors=[];delta_errors=[]
    for a in (-.12,.19,.43):
        up,dup=evolution(a);U=um.conj().T@up;dU=um.conj().T@dup;g=gamma(U)
        dg=g@second(U.conj().T@dU)
        vals=[]
        for r in rs:
            rho=np.kron(old.receiver_state(r,m),rho_ext)
            expected=np.trace(rho@g);expected_d=np.trace(rho@dg)
            got,der=four_det(U,Cperp,r,m,dU);sc,z,A=schur(U,Cperp,r,m)
            errors.extend([abs(got-expected),abs(sc-expected)])
            source_errors.append(abs(der-expected_d))
            assert min(weights(r,m))>=0 and abs(sum(weights(r,m))-1)<1e-15
            Cmean=sum(w*C for w,C in zip(weights(r,m),branch_covariances(Cperp)))
            Ctarget=np.zeros((n,n),complex);Ctarget[:2,:2]=(I+m*X)/2;Ctarget[2:,2:]=Cperp
            cov_errors.append(np.max(abs(Cmean-Ctarget)))
            # Illegally omitting exterior-mode excursions changes the weight.
            Aq=U[:2,:2]
            bare=z*((1-r)*(1+np.linalg.det(Aq))/2+r*np.trace(Aq)/2+m*np.trace(X@Aq)/2)
            wrong.append(abs(bare-got))
            vals.append(got)
        _,z,A=schur(U,Cperp,rs[0],m)
        delta_errors.append(abs(vals[2]-vals[0]+(rs[2]-rs[0])*z*np.linalg.det(I-A)/2))
    # First source at ANY equal pair of histories is insensitive to r.
    equal_source_errors=[]
    for s in (-.18,.07,.31):
        u,du=evolution(s);tangent=u.conj().T@du
        ders=[four_det(np.eye(n,dtype=complex),Cperp,r,m,tangent)[1] for r in rs]
        equal_source_errors.append(abs(ders[2]-ders[0]))
    # Second unequal-history coefficient knows the missing fourth moment.
    ub,dub=evolution(.07);tangent=ub.conj().T@dub
    coefficient=-(rs[2]-rs[0])*np.linalg.det(tangent[:2,:2])/2
    tail=[]
    for eps in (.02,.01,.005):
        out=[]
        for e in (-eps,eps):
            u,_=evolution(.07+e);U=ub.conj().T@u
            out.append(four_det(U,Cperp,rs[2],m)[0]-four_det(U,Cperp,rs[0],m)[0])
        estimate=sum(out)/(2*eps*eps)
        tail.append(dict(step=eps,error=float(abs(estimate-coefficient))))
    assert max(errors)<3e-14 and max(source_errors)<3e-14 and max(delta_errors)<3e-14
    assert max(cov_errors)<1e-14 and max(equal_source_errors)<2e-14 and max(wrong)>1e-5
    assert tail[-1]['error']<tail[0]['error']/8
    return dict(CAR_modes=n,Fock_dimension=2**n,
        full_Fock_vs_four_determinants_and_Schur_error=float(max(errors)),
        full_source_insertion_vs_cofactor_derivative_error=float(max(source_errors)),
        complete_initial_covariance_error=float(max(cov_errors)),
        explicit_nonGaussian_history_difference_error=float(max(delta_errors)),
        same_history_first_source_difference=float(max(equal_source_errors)),
        omitted_exterior_excursion_weight_error=float(max(wrong)),
        second_history_coefficient=[float(coefficient.real),float(coefficient.imag)],
        second_coefficient_convergence=tail)

def singular_case():
    C=np.diag([.5,.3]);U=np.eye(4,dtype=complex)
    U[:2,:2]=old.previous.prior.herm_fun(.4*Y,lambda x:np.exp(-1j*x))
    U[2,2]=-1
    dU=U@np.diag([0,0,-1j,0])
    modes=old.annihilators(4)
    h=U.conj().T@dU
    dg=gamma(U)@sum(h[i,j]*modes[i].conj().T@modes[j] for i in range(4) for j in range(4))
    ext=np.diag([.5*.7,.5*.3,.5*.7,.5*.3])
    rho=np.kron(old.receiver_state(.8,.5),ext)
    val,der=four_det(U,C,.8,.5,dU)
    err=abs(der-np.trace(rho@dg))
    assert abs(val)<1e-15 and err<1e-14 and abs(der)>.1
    return dict(weight_abs=float(abs(val)),source_abs=float(abs(der)),
        source_error=float(err),division_by_Schur_or_weight_used=False)

def run():
    cases=[finite_case(n) for n in (4,6)]
    singular=singular_case()
    # Actual870 current coefficient with original863 source tangent.
    _,_,_,d=old.previous.original_menu()
    m=.5;r0=.5;r1=1.;x=-old.previous.GAINS[2]*float(d[2])/(2*m)
    coefficient=-(r1-r0)*x*x/2
    oldvar=json.loads((HERE.parent/'870/receiver_four_point_backreaction_results.json').read_text('utf-8'))['original_companion_lambda2_hbar_variance_difference_coefficient']
    assert abs(-2*coefficient-oldvar)<1e-22
    return dict(round=871,date='2026-10-06',formal_reports=871,
        cumulative_numbered_groups=3656,fresh_numbered_groups=1,all_checks_passed=True,
        exact_positive_quasifree_input_branches=4,finite_history_cases=cases,
        zero_weight_source_check=singular,
        original_companion_characteristic_quadratic_coefficient=coefficient,
        original870_variance_recovered=-2*coefficient,
        continuum_preparation_identity='exact normal-state convex decomposition; not Gaussian time evolution',
        determinant_scope='finite CAR regulator and conditional classical bilinear receiver histories',
        same_fixed_formal_S_sources_records_and_renormalization=True,
        arbitrary_classical_history_quadratic_mean_independent_of_r=True,
        continuous_all_order_Fredholm_determinant_proved=False,
        original_four_leg_total_source_evaluated=False,
        total_mean_geometry_or_finite_coupling_proved=False,
        full_goal_completed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
