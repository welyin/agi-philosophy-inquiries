"""581 entry only: local EOM decomposition of the same scalar heat coefficient."""
from pathlib import Path
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import joint_scalar_effective_geometry as parent

HERE=Path(__file__).resolve().parent
TARGET=HERE/'operator_reduction_entry_results.json'
original=parent.entry.original


def metric_factorization_check():
    rng=np.random.default_rng(581)
    rows=[]
    for j in range(4):
        X=rng.normal(size=(4,5))*.3
        A=rng.normal(size=(5,5))*.2;A=(A+A.T)/2
        B=X@X.T;M=X.T@X;S=float(np.trace(B));Q=float(np.trace(B@B))
        U=.13+.07*j
        Ric=rng.normal(size=(4,4))*.3;Ric=(Ric+Ric.T)/2
        if j==3:Ric=B+U*np.eye(4)
        R=float(np.trace(Ric));ric2=float(np.trace(Ric@Ric))
        riem2=1.37+2*ric2-R*R/3
        euler=riem2-4*ric2+R*R
        full=parent.entry.invariant_coefficient(X,A,R)+5*((riem2-ric2)/180+R*R/72)
        reduced=(euler/36-7*S*S/216+11*Q/108+U*S/18+U*U
            -2*U*np.trace(A)/3+np.trace(A@A)/2-np.trace(A@M)/6)
        E=Ric-B-U*np.eye(4);e=float(np.trace(E))
        D=(Ric+B+U*np.eye(4))/12+np.eye(4)*((R+S+4*U)/24-np.trace(A)/6-S/9)
        H=-(E-e*np.eye(4)/2)/2
        T=-2*D+np.trace(D)*np.eye(4)
        factored=float(np.sum(E*D));via_euler=float(np.sum(H*T))
        assert abs(full-reduced-factored)<1e-14
        assert abs(factored-via_euler)<1e-14
        if j==3:assert abs(full-reduced)<1e-14
        rows.append(dict(case=j,full=float(full),reduced=float(reduced),
                         E_contraction=factored,metric_euler_contraction=via_euler,
                         residual=float(full-reduced-factored),leading_metric_equation=j==3))
    return rows


def original_potential_chain_check():
    x=np.linspace(0,2*np.pi,4096,endpoint=False)
    q=.73+.17*np.sin(x);qp=.17*np.cos(x);qpp=-.17*np.sin(x)
    direction=np.array([0.,np.sqrt(.6),0.,0.,np.sqrt(.4)])
    def phi(t):return np.sqrt(6*original.M)*np.tanh(t[...,None]/np.sqrt(6))*direction
    def U(t):return original.node_potential(phi(t))
    dq=1e-4
    first=(U(q+dq)-U(q-dq))/(2*dq)
    second=(U(q+dq)-2*U(q)+U(q-dq))/(dq*dq)
    left=-2*np.pi*np.mean(second*qp*qp)/6
    right=2*np.pi*np.mean(first*qpp)/6
    assert abs(left-right)<2e-8
    # Direct complete-target gradient, not a one-field replacement of U.
    position=phi(q);gradient=np.zeros_like(position)
    for a in range(5):
        step=np.zeros(5);step[a]=2e-5
        gradient[:,a]=(original.node_potential(position+step)-original.node_potential(position-step))/(4e-5)
    inverse=original.inverse(position)
    gradnorm=np.einsum('...a,...ab,...b->...',gradient,inverse,gradient)
    tangent=np.sqrt(original.M)/np.cosh(q[...,None]/np.sqrt(6))**2*direction
    grad_tension=np.sum(gradient*tangent,axis=1)*qpp
    grad_times_Ephi=gradnorm-grad_tension
    via_euler=2*np.pi*np.mean((gradnorm-grad_times_Ephi)/6)
    assert abs(via_euler-left)<2e-8
    return dict(original_U_negative_Hessian_gradient_integral=float(left),
                integration_by_parts=float(right),full_target_Ephi_identity=float(via_euler),
                error=float(abs(via_euler-left)),geodesic_profile_not_an_on_shell_solution=True)


def formal_metric_shift_check():
    # Pure RS block: polynomial coefficient movement at first order in lambda.
    alpha=1/27;beta=1/54;c=-1/9
    shift=-2*c
    rng=np.random.default_rng(5812)
    errors=[]
    for _ in range(6):
        S,Q,U,R=rng.normal(size=4)
        old=alpha*S*S+beta*Q+c*R*S
        delta=shift*(R*S/2-S*S/2-2*U*S)
        new=(alpha+c)*S*S+beta*Q+4*c*U*S
        errors.append(float(abs(old+delta-new)))
    assert max(errors)<1e-15
    return dict(old_RS_coefficient=c,inverse_metric_shift_coefficient=shift,
                shifted_S_squared_coefficient=alpha+c,shifted_U_S_coefficient=4*c,
                max_polynomial_error=max(errors),
                exact_variable_change_and_second_order_remainder_not_yet_tested=True)


def run():
    return dict(status='581 candidate entry only',checks_passed=3,
        metric_factorization=metric_factorization_check(),
        scalar_chain=original_potential_chain_check(),RS_block=formal_metric_shift_check(),
        scope='four-dimensional dynamic Euclidean leading Einstein-scalar action, kappa=1, same K and U; scalar-loop local coefficient, zero gauge curvature; formal first-order reduction, no full independent operator basis or quantum-frame equivalence')


if __name__=='__main__':
    result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if '--write' in sys.argv:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
