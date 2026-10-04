"""580 entry: scalar-loop heat-kernel structures in the same Einstein-frame H5 target.
This is a continuum background-field input, not a renormalized limit of the graph.
"""
from pathlib import Path
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'scalar_effective_action_entry_results.json'
CURVATURE=-1/6
_,vacuum,_=original.lattice.scalar.parameters()
CENTER=np.array([0.,np.sqrt(vacuum[0]),0.,0.,np.sqrt(vacuum[1])])


def trace_coefficient(X,A=None,spacetime_R=0.):
    dim=X.shape[1]
    if A is None:A=np.zeros((dim,dim))
    M=X.T@X;S=np.trace(M)
    P=A-CURVATURE*(S*np.eye(dim)-M)
    omega_trace=0.
    for x in X:
        for y in X:
            omega=CURVATURE*(np.outer(x,y)-np.outer(y,x))
            omega_trace+=np.trace(omega@omega)
    # Pure background Riemann/Ricci/R^2 terms are excluded on both sides.
    return float(np.trace(P@P)/2+omega_trace/12-spacetime_R*np.trace(P)/6)


def invariant_coefficient(X,A=None,spacetime_R=0.):
    dim=X.shape[1]
    if A is None:A=np.zeros((dim,dim))
    M=X.T@X;S=np.trace(M);Q=np.trace(M@M);k=CURVATURE
    return float(np.trace(A@A)/2-k*(S*np.trace(A)-np.trace(A@M))
        +k*k*((3*dim-7)*S*S/6+2*Q/3)
        -spacetime_R*np.trace(A)/6+k*(dim-1)*spacetime_R*S/6)


def potential_hessian(step):
    E=np.eye(5)*step;out=np.zeros((5,5))
    for a in range(5):
        for b in range(5):
            out[a,b]=(original.node_potential(CENTER+E[a]+E[b])
                -original.node_potential(CENTER+E[a]-E[b])
                -original.node_potential(CENTER-E[a]+E[b])
                +original.node_potential(CENTER-E[a]-E[b]))/(4*step*step)
    eigen,vectors=np.linalg.eigh(original.inverse(CENTER))
    root=(vectors*np.sqrt(eigen))@vectors.T
    return root@out@root


def run():
    rng=np.random.default_rng(580)
    X=rng.normal(size=(4,5))*.23
    A=potential_hessian(1e-4);A2=potential_hessian(5e-5)
    eig=np.linalg.eigvalsh(A2)
    assert np.max(abs(eig[:3]))<2e-9
    assert np.max(abs(eig[3:]-np.array([.04683851043525822,.13805672648598155])))<2e-8
    assert np.max(abs(A2-A))<3e-8
    trace=trace_coefficient(X,A2,.31);invariant=invariant_coefficient(X,A2,.31)
    assert abs(trace-invariant)<1e-15
    rot=np.linalg.qr(rng.normal(size=(5,5)))[0]
    spacetime_rot=np.linalg.qr(rng.normal(size=(4,4)))[0]
    covariant=trace_coefficient(spacetime_rot@X@rot,rot.T@A2@rot,.31)
    assert abs(covariant-trace)<1e-14
    S=float(np.sum(X*X));Q=float(np.trace((X.T@X)@(X.T@X)))
    four_gradient=S*S/27+Q/54
    direct_four=trace_coefficient(X)
    assert abs(direct_four-four_gradient)<1e-15 and direct_four>0
    rank_one=np.zeros((4,5));rank_one[0,0]=.7
    assert abs(trace_coefficient(rank_one)-.7**4/18)<1e-15
    # Scaling X by t isolates the quartic derivative contribution independently of A.
    z=np.array([0.,1.,2.])
    values=np.array([trace_coefficient(np.sqrt(v)*X,A2,.31) for v in z])
    quartic_fit=(values[2]-2*values[1]+values[0])/2
    shifted_values=np.array([trace_coefficient(np.sqrt(v)*X,A2+.73*np.eye(5),.31) for v in z])
    shifted_fit=(shifted_values[2]-2*shifted_values[1]+shifted_values[0])/2
    assert abs(quartic_fit-four_gradient)<1e-14 and abs(shifted_fit-four_gradient)<1e-14
    # R-X mixed finite difference removes potential-only and pure-curvature pieces.
    mixed=(trace_coefficient(X,A2,1)-trace_coefficient(X,A2,0)
           -trace_coefficient(np.zeros_like(X),A2,1)+trace_coefficient(np.zeros_like(X),A2,0))
    assert abs(mixed+S/9)<1e-14
    direction=np.array([.13,-.07,.04,.06,.11])
    deltaF=-float(CENTER@direction)/3;F=float(original.F(CENTER))
    eps=1e-5
    jordan_fixed=np.eye(4)/F
    changed_E=(original.F(CENTER+eps*direction)-original.F(CENTER-eps*direction))/(2*eps)*jordan_fixed
    expected_E=deltaF/F*np.eye(4)
    changed_J=(1/original.F(CENTER+eps*direction)-1/original.F(CENTER-eps*direction))/(2*eps)*np.eye(4)
    expected_J=-deltaF/F**2*np.eye(4)
    assert np.max(abs(changed_E-expected_E))<1e-10
    assert np.max(abs(changed_J-expected_J))<1e-10 and np.linalg.norm(changed_E)>1e-4
    return dict(status='580 candidate entry only',checks_passed=4,
        inherited_potential=dict(vacuum=CENTER.tolist(),full_five_field_hessian_eigenvalues=eig.tolist(),
            finite_step_change=float(np.max(abs(A2-A))),angular_modes_not_claimed_physical_particles=True),
        covariant_trace=dict(direct=trace,invariant=invariant,error=abs(trace-invariant),
                             basis_change_error=abs(trace-covariant)),
        derivative_structures=dict(S=S,Q=Q,four_gradient_coefficient=four_gradient,
            rank_one_coefficient=trace_coefficient(rank_one),scaled_gradient_fit=float(quartic_fit),
            constant_mass_shift_fit=float(shifted_fit),R_gradient_cross=mixed,expected_R_gradient_cross=-S/9,
            heat_kernel_a4_density_without_overall_factor_not_Gamma_pole_or_finite_matching=True),
        frozen_frame_check=dict(deltaF=deltaF,fixed_J_delta_E=expected_E[0,0],
             fixed_E_delta_J=expected_J[0,0],error=max(float(np.max(abs(changed_E-expected_E))),
                                                   float(np.max(abs(changed_J-expected_J))))),
        scope='only scalar determinant on fixed Euclidean Einstein background, zero background gauge curvature; target H5 and original U; no full gauge/fermion/gravity loops, on-shell operator-basis classification, graph-continuum renormalization or physical frame inequivalence')


if __name__=='__main__':
    result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if '--write' in sys.argv:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
