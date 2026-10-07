"""Working875: conserved-current transport for a declared finite clock constraint.
Finite history and its source derivative are transported together. No original
field quantum constraint or complete graph-continuum matching is inferred.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'nonstationary_clock_current_probe_results.json'
sys.path.insert(0,str(HERE.parent/'874'))
from constraint_clock_kernel_probe import lift,annihilation
I=np.eye(4,dtype=complex);D=np.diag([.025,-.02]);O=np.array([[0,.04],[.04,0]])
BASE=np.diag([.08,.12]);SRC=lift(np.diag([1.,0.]))
def hermfun(A,f):
    vals,Q=np.linalg.eigh((A+A.conj().T)/2)
    return (Q*f(vals))@Q.conj().T
def sylvester(T,B):
    vals,Q=np.linalg.eigh(T);bt=Q.conj().T@B@Q
    return Q@(bt/(vals[:,None]+vals[None,:]))@Q.conj().T
def profile(t):
    x=t/2
    if x<=0 or x>=1:return 0.,0.
    d=x*(1-x);g=np.exp(4-1/d);gp=g*(1-2*x)/(2*d*d)
    return g,gp
def coeff(t,s=0):
    g,gp=profile(t)
    A=I-lift(BASE+g*D+g*g*O)-s*g*SRC
    Ap=-lift(gp*D+2*g*gp*O)-s*gp*SRC
    As=-g*SRC
    return A,Ap,As
def state_parts(t,y):
    u,w,us,ws,V=y
    inv=np.linalg.inv(u);A,_,_=coeff(t)
    P=-1j*w@inv;Pp=1j*(A-P@P);G=P+P.conj().T
    vals=np.linalg.eigvalsh(G);assert vals.min()>0
    T=hermfun(G,np.sqrt);Ti=np.linalg.inv(T)
    Tp=sylvester(T,Pp+Pp.conj().T)
    H=I-T@P@Ti+1j*Tp@Ti
    U=np.exp(-1j*t)*T@u
    Ps=-1j*(ws@inv-w@inv@us@inv);Gs=Ps+Ps.conj().T
    Ts=sylvester(T,Gs)
    dU=np.exp(-1j*t)*(Ts@u+T@us)
    wrong=np.exp(-1j*t)*T@us
    return dict(U=U,dU=dU,wrong=wrong,H=H,
        current=-1j*(u.conj().T@w-w.conj().T@u),
        metric_min=float(vals.min()),
        missing_connection_nonhermiticity=float(np.linalg.norm((-T@P@Ti)-(-T@P@Ti).conj().T)),
        naive_U=np.exp(-1j*t)*V)
def integrate(steps,s=0,variations=True):
    A0,_,_=coeff(0,s);E0=hermfun(A0,np.sqrt)
    u0=hermfun(2*E0,lambda x:1/np.sqrt(x));w0=1j*E0@u0
    y=np.stack([u0,w0,np.zeros((4,4),complex),np.zeros((4,4),complex),I.copy()])
    dt=2/steps
    def rhs(t,z):
        u,w,us,ws,V=z;A,_,As=coeff(t,s);E=hermfun(A,np.sqrt)
        return np.stack([w,-A@u,ws,-A@us-As@u if variations else np.zeros_like(u),1j*E@V])
    samples=[]
    for j in range(steps):
        t=j*dt
        k1=rhs(t,y);k2=rhs(t+dt/2,y+dt*k1/2);k3=rhs(t+dt/2,y+dt*k2/2);k4=rhs(t+dt,y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        if s==0 and j+1 in (steps//4,steps//2,3*steps//4,steps):
            r=state_parts((j+1)*dt,y)
            samples.append(dict(time=(j+1)*dt,
                current_error=float(np.linalg.norm(r['current']-I)),
                metric_min=r['metric_min'],
                hermitian_generator_error=float(np.linalg.norm(r['H']-r['H'].conj().T)),
                missing_connection_nonhermiticity=r['missing_connection_nonhermiticity']))
    if s==0:final=state_parts(2,y)
    else:
        u,w=y[:2];P=-1j*w@np.linalg.inv(u);T=hermfun(P+P.conj().T,np.sqrt)
        final=dict(U=np.exp(-2j)*T@u)
    return final,samples
def run():
    low,_=integrate(256);high,samples=integrate(512)
    numerical_convergence=float(np.linalg.norm(high['U']-low['U']))
    assert numerical_convergence<2e-10
    assert max(x['current_error'] for x in samples)<1e-10
    assert max(x['hermitian_generator_error'] for x in samples)<1e-11
    assert max(x['missing_connection_nonhermiticity'] for x in samples)>1e-3
    U,dU=high['U'],high['dU']
    tangent=float(np.linalg.norm(dU.conj().T@U+U.conj().T@dU))
    wrong_tangent=float(np.linalg.norm(high['wrong'].conj().T@U+U.conj().T@high['wrong']))
    assert tangent<1e-9 and wrong_tangent>1e-3
    eps=1e-4
    plus,_=integrate(512,eps,False);minus,_=integrate(512,-eps,False)
    derivative_error=float(np.linalg.norm((plus['U']-minus['U'])/(2*eps)-dU))
    assert derivative_error<1e-8
    rho=np.zeros((4,4),complex);rho[0,0]=rho[3,3]=.1
    rho[1:3,1:3]=np.array([[.8,.5],[.5,.8]])/2
    c=annihilation(0);P=c.conj().T@c
    probability=float(np.trace(P@U@rho@U.conj().T).real)
    dp=float(2*np.trace(P@dU@rho@U.conj().T).real)
    pplus=float(np.trace(P@plus['U']@rho@plus['U'].conj().T).real)
    pminus=float(np.trace(P@minus['U']@rho@minus['U'].conj().T).real)
    assert abs((pplus-pminus)/(2*eps)-dp)<1e-8
    A,Ap,_=coeff(.5);E=hermfun(A,np.sqrt);Ep=sylvester(E,Ap)
    instant_residual=float(np.linalg.norm(Ep))
    naive_difference=float(np.linalg.norm(U-high['naive_U']))
    assert instant_residual>1e-3 and naive_difference>1e-3
    old=json.loads((HERE.parent/'861/magnetic_reduced_hamiltonian_results.json').read_text('utf-8'))['original_859_background_clock_reduction'][1]
    return dict(working_round=875,formal_reports=874,cumulative_numbered_groups=3659,
        fresh_numbered_groups=0,all_probe_checks_passed=True,
        dimensionless_time_scaling_hstar=old['actual_clock_speed']/(2*old['a']),
        smooth_history_is_declared_calibration=True,transport_checks=samples,
        RK4_256_512_process_difference=numerical_convergence,
        final_unitarity_error=float(np.linalg.norm(U.conj().T@U-I)),
        source_unitarity_tangent_error=tangent,
        omitted_metric_source_tangent_defect=wrong_tangent,
        source_variational_finite_difference_error=derivative_error,
        record_probability=probability,record_source_derivative=dp,
        instantaneous_root_constraint_residual_at_time_half=instant_residual,
        transported_vs_instantaneous_process_difference=naive_difference,
        full_original_quantum_constraint_ordering_fixed=False,
        formal_current_construction_is_full_Hadamard_or_QG_proof=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
