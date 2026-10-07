"""Working875: noncommuting clock coefficients and complete endpoint/record transport.
A declared finite constraint model; not the original field quantum ordering.
Constant mass-source labels change both the equation and its initial current data.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from nonstationary_clock_current_probe import hermfun,sylvester,profile
sys.path.insert(0,str(HERE.parent/'874'))
from constraint_clock_kernel_probe import lift,annihilation
TARGET=HERE/'clock_boundary_record_transport_results.json'
I=np.eye(4,dtype=complex)
X=lift(np.array([[0.,1.],[1.,0.]]))
Y=lift(np.array([[0.,-1j],[1j,0.]]))
N1=lift(np.diag([1.,0.]));N2=lift(np.diag([0.,1.]))
BASE=lift(np.diag([.08,.12]))
D=lift(np.diag([.025,-.02]));O=.04*X
SOURCES=[N1,N2]
def coeff(t,s):
    g,gp=profile(t)
    alpha=I+.08*g*X
    beta=2*I+.1*g*Y
    zeta=BASE+g*D+g*g*O+s[0]*N1+s[1]*N2
    return alpha,beta,zeta
def initial(s,transport=True):
    _,_,z=coeff(0,s)
    R=hermfun(I-z,np.sqrt);T=hermfun(2*R,np.sqrt);u=np.linalg.inv(T)
    m=R@u
    data=[u,m]
    for source in SOURCES:
        Rs=sylvester(R,-source)
        Ts=sylvester(T,2*Rs)
        us=-u@Ts@u if transport else np.zeros_like(u)
        ms=Rs@u+R@us if transport else np.zeros_like(u)
        data.extend([us,ms])
    return np.stack(data)
def parts(t,y,s):
    u,m=y[:2];ui=np.linalg.inv(u)
    a,b,z=coeff(t,s);ai=np.linalg.inv(a)
    R=m@ui;G=R+R.conj().T
    assert np.linalg.eigvalsh(G).min()>0
    T=hermfun(G,np.sqrt);Ti=np.linalg.inv(T)
    Rp=-1j*((R+b/2)@ai@(R-b/2)+z)
    Tp=sylvester(T,Rp+Rp.conj().T)
    H=-T@ai@(R-b/2)@Ti+1j*Tp@Ti
    U=T@u
    dU=[];wrong=[]
    for j in range(2):
        us,ms=y[2+2*j:4+2*j]
        Rs=ms@ui-R@us@ui
        Ts=sylvester(T,Rs+Rs.conj().T)
        dU.append(Ts@u+T@us)
        wrong.append(T@us)
    return dict(U=U,dU=dU,wrong=wrong,u=u.copy(),m=m.copy(),T=T,R=R,H=H,
        current=u.conj().T@m+m.conj().T@u,
        metric_min=float(np.linalg.eigvalsh(G).min()),
        coefficient_commutator=float(np.linalg.norm(a@b-b@a)))
def integrate(steps=512,s=(0.,0.),transport=True):
    y=initial(s,transport);dt=2/steps;snapshots=[]
    def rhs(t,y):
        a,b,z=coeff(t,s);ai=np.linalg.inv(a);rows=[]
        for j in range(3):
            u,m=y[2*j:2*j+2]
            du=ai@(m-b@u/2)
            dm=-b@du/2-z@u
            if j:dm-=SOURCES[j-1]@y[0]
            rows.extend([1j*du,1j*dm])
        return np.stack(rows)
    for k in range(steps):
        t=k*dt
        k1=rhs(t,y);k2=rhs(t+dt/2,y+dt*k1/2)
        k3=rhs(t+dt/2,y+dt*k2/2);k4=rhs(t+dt,y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        if k+1 in (steps//2,steps):snapshots.append(parts((k+1)*dt,y,s))
    return snapshots
def menu(data):
    first,last=data;U1,U2=first['U'],last['U'];U21=U2@U1.conj().T
    rho=np.zeros((4,4),complex);rho[0,0]=rho[3,3]=.1
    rho[1:3,1:3]=np.array([[.8,.5],[.5,.8]])/2
    c=(annihilation(0)+annihilation(1))/np.sqrt(2)
    P=c.conj().T@c;Q=N1
    probs=[];derivatives=[[],[]]
    for proj in (P,I-P):
        for final in (Q,I-Q):
            K=final@U21@proj@U1
            probs.append(float(np.trace(K@rho@K.conj().T).real))
            for j in range(2):
                d1,d2=first['dU'][j],last['dU'][j]
                d21=d2@U1.conj().T+U2@d1.conj().T
                dK=final@(d21@proj@U1+U21@proj@d1)
                derivatives[j].append(float(2*np.trace(dK@rho@K.conj().T).real))
    # Correct lift of a middle-record projection to both Cauchy components.
    t=first['T'];r=first['R'];u=first['u'];m=first['m']
    up=np.linalg.solve(t,P@t@u);mp=r@up
    actual=up.conj().T@mp+mp.conj().T@up
    target=U1.conj().T@P@U1
    # Applying the same matrix to both raw components is a different instrument.
    bad=u.conj().T@P@m+m.conj().T@P@u
    return dict(probabilities=np.array(probs),derivatives=np.array(derivatives),
        lifted_record_error=float(np.linalg.norm(actual-target)),
        untransported_record_minimum_current=float(np.linalg.eigvalsh((bad+bad.conj().T)/2).min()),
        untransported_record_error=float(np.linalg.norm(bad-target)))
def run():
    lo=integrate(256);hi=integrate(512);record=menu(hi)
    convergence=float(np.linalg.norm(hi[-1]['U']-lo[-1]['U']))
    assert convergence<2e-10
    eps=1e-4;process_errors=[];record_errors=[]
    for j in range(2):
        s=np.zeros(2);s[j]=eps
        plus=integrate(512,tuple(s));minus=integrate(512,tuple(-s))
        numeric=(plus[-1]['U']-minus[-1]['U'])/(2*eps)
        process_errors.append(float(np.linalg.norm(numeric-hi[-1]['dU'][j])))
        record_errors.append(float(np.max(np.abs((menu(plus)['probabilities']-menu(minus)['probabilities'])/(2*eps)-record['derivatives'][j]))))
    assert max(process_errors)<1e-8 and max(record_errors)<1e-8
    wrong=integrate(512,transport=False)
    missing_initial=max(float(np.linalg.norm(wrong[-1]['dU'][j]-hi[-1]['dU'][j])) for j in range(2))
    metric_omission=max(float(np.linalg.norm(d.conj().T@hi[-1]['U']+hi[-1]['U'].conj().T@d)) for d in hi[-1]['wrong'])
    assert missing_initial>1e-3 and metric_omission>1e-3
    assert record['lifted_record_error']<1e-11
    assert record['untransported_record_minimum_current']< -1e-6
    assert abs(record['probabilities'].sum()-1)<1e-10 and record['probabilities'].min()>0
    assert np.max(np.abs(record['derivatives'].sum(axis=1)))<1e-9
    current=max(float(np.linalg.norm(d['current']-I)) for d in hi)
    hermitian=max(float(np.linalg.norm(d['H']-d['H'].conj().T)) for d in hi)
    unitary=max(float(np.linalg.norm(d['U'].conj().T@d['U']-I)) for d in hi)
    assert max(current,hermitian,unitary)<1e-10
    return dict(working_round=875,formal_reports=874,cumulative_numbered_groups=3659,
        fresh_numbered_groups=0,all_checks_passed=True,
        declared_finite_noncommuting_coefficient_model=True,
        constant_mass_source_initial_data_transported=True,
        current_error=current,hermitian_generator_error=hermitian,unitarity_error=unitary,
        coefficient_commutator_norms=[d['coefficient_commutator'] for d in hi],
        RK4_256_512_process_difference=convergence,
        constant_mass_source_process_derivative_errors=process_errors,
        constant_mass_source_record_derivative_errors=record_errors,
        omitted_initial_transport_process_derivative_error=missing_initial,
        omitted_metric_transport_unitarity_tangent_defect=metric_omission,
        joint_record_probabilities=record['probabilities'].tolist(),
        joint_record_mass_derivatives=record['derivatives'].tolist(),
        correctly_lifted_record_current_error=record['lifted_record_error'],
        raw_component_projection_minimum_current=record['untransported_record_minimum_current'],
        raw_component_projection_current_error=record['untransported_record_error'],
        original_full_field_quantum_ordering_matched=False,
        continuum_or_scale_limit_proved=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
