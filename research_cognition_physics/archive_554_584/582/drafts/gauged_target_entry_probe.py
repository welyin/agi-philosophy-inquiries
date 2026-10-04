"""582 candidate: original U(1) charge, H5 Killing connection and scalar heat kernel."""
from pathlib import Path
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'gauged_target_entry_results.json'


def run():
    phi=np.array([.11,.68,-.07,.13,.52])
    K=original.metric(phi);inv=original.inverse(phi)
    eigen,vectors=np.linalg.eigh(inv)
    frame=(vectors*np.sqrt(eigen))@vectors.T
    invframe=np.linalg.inv(frame)
    generator=np.zeros((5,5))
    generator[0,2]=generator[1,3]=-3
    generator[2,0]=generator[3,1]=3
    killing=generator@phi
    J=generator+np.outer(killing,phi)/(6*original.F(phi))
    orth=invframe@J@frame
    assert np.max(abs(K@J+J.T@K))<1e-14
    assert np.max(abs(orth+orth.T))<1e-14
    # A finite original U(1) transformation preserves K and U.
    theta=.217;rot=np.eye(5)
    for a,b in ((0,2),(1,3)):
        rot[a,a]=rot[b,b]=np.cos(3*theta)
        rot[a,b]=-np.sin(3*theta);rot[b,a]=np.sin(3*theta)
    moved=rot@phi
    assert np.max(abs(rot.T@original.metric(moved)@rot-K))<1e-14
    assert abs(original.node_potential(moved)-original.node_potential(phi))<1e-15
    rng=np.random.default_rng(582)
    X=rng.normal(size=(4,5))*.3
    f=np.zeros((4,4));f[0,1]=.29;f[1,0]=-.29;f[2,3]=-.17;f[3,2]=.17
    zero=plus=minus=quadratic=cross=0.
    k=-1/6
    for mu in range(4):
        for nu in range(4):
            omega=k*(np.outer(X[mu],X[nu])-np.outer(X[nu],X[mu]))
            addition=f[mu,nu]*orth
            zero+=np.trace(omega@omega)/12
            plus+=np.trace((omega+addition)@(omega+addition))/12
            minus+=np.trace((omega-addition)@(omega-addition))/12
            quadratic+=np.trace(addition@addition)/12
            cross+=np.trace(omega@addition)/6
    assert abs(plus-zero-quadratic-cross)<1e-14
    assert abs(minus-zero-quadratic+cross)<1e-14
    assert abs(cross)>1e-4 and abs(plus-minus-2*cross)<1e-14
    # Euclidean gauge stress shape: normalization belongs to the common gauge action.
    tau=f@f.T-np.eye(4)*np.sum(f*f)/4
    reversed_tau=(-f)@(-f).T-np.eye(4)*np.sum((-f)*(-f))/4
    assert np.max(abs(tau-reversed_tau))==0 and abs(np.trace(tau))<1e-15
    return dict(status='582 candidate entry only',checks_passed=2,
        killing_connection=dict(K_skew_error=float(np.max(abs(K@J+J.T@K))),
                                trace_J_squared=float(np.trace(orth@orth)),original_Higgs_charge=3),
        scalar_connection_coefficient=dict(zero=float(zero),positive_field_strength=float(plus),
            reversed_field_strength=float(minus),pure_gauge_quadratic=float(quadratic),mixed=float(cross)),
        classical_gauge_stress=dict(sign_reversal_difference=float(np.max(abs(tau-reversed_tau))),
                                   trace=float(np.trace(tau))),
        scope='local jets, original Abelian generator in full H5 target, scalar determinant connection term only; no full on-shell Gauss/Einstein source, full gauge loop or independent operator classification')


if __name__=='__main__':
    result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if '--write' in sys.argv:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
