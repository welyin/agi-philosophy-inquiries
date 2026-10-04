"""582: original gauge action, scalar-loop connection and common metric reduction."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_scalar_effective_geometry as scalar

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauged_effective_geometry_results.json'
spec=importlib.util.spec_from_file_location('gauged_entry',HERE/'round582_drafts/gauged_target_entry_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def real_generator(z):
    out=np.zeros((5,5))
    out[:4,:4]=np.block([[z.real,-z.imag],[z.imag,z.real]])
    return out


GENERATORS=np.stack([real_generator(1j*t) for t in original.lattice.T]+[real_generator(3j*np.eye(2))])
INDEX=np.array([1.,1.,1.,36.])


def target_connection(phi):
    K=original.metric(phi)
    eigen,vectors=np.linalg.eigh(original.inverse(phi))
    frame=(vectors*np.sqrt(eigen))@vectors.T
    inverse_frame=np.linalg.inv(frame)
    killing=np.einsum('Aab,b->Aa',GENERATORS,phi)
    J=GENERATORS+killing[:,:,None]*phi[None,None,:]/(6*original.F(phi))
    orth=inverse_frame@J@frame
    return K,frame,inverse_frame,killing,J,orth


def connection_correction(X,strength,J):
    zero=full=pure=mixed=0.
    for mu in range(4):
        for nu in range(4):
            omega=-(np.outer(X[mu],X[nu])-np.outer(X[nu],X[mu]))/6
            G=np.einsum('A,Aab->ab',strength[:,mu,nu],J)
            zero+=np.trace(omega@omega)/12
            full+=np.trace((omega+G)@(omega+G))/12
            pure+=np.trace(G@G)/12
            mixed+=np.trace(omega@G)/6
    assert abs(full-zero-pure-mixed)<1e-13
    return dict(correction=float(full-zero),pure=float(pure),mixed=float(mixed))


def mass_connection_gram_check():
    _,vacuum,_=original.lattice.scalar.parameters()
    points=[np.zeros(5),np.array([.11,.68,-.07,.13,.52]),
            np.array([0.,np.sqrt(vacuum[0]),0.,0.,np.sqrt(vacuum[1])])]
    rows=[]
    for phi in points:
        K,_,_,killing,_,J=target_connection(phi)
        mass=killing@K@killing.T
        gram=-np.einsum('Aab,Bba->AB',J,J)
        expected=np.diag(INDEX)+mass/3
        assert np.max(abs(gram-expected))<3e-14
        assert np.max(abs(J+J.transpose(0,2,1)))<3e-14
        rows.append(dict(phi=phi.tolist(),trace_gram_error=float(np.max(abs(gram-expected))),
                         gram_eigenvalues=np.linalg.eigvalsh(gram).tolist()))
    phi=points[-1];K,_,_,killing,_,_=target_connection(phi)
    mass=killing@K@killing.T
    _,gw2,gy2=original.lattice.PAR['g2']
    scale=np.diag([np.sqrt(gw2)]*3+[np.sqrt(gy2)/6])
    physical=scale@mass@scale
    masses=np.linalg.eigvalsh(physical)
    F=original.F(phi);h2=vacuum[0]
    expected=np.sort([0.,gw2*h2/(4*F),gw2*h2/(4*F),(gw2+gy2)*h2/(4*F)])
    assert np.max(abs(masses-expected))<1e-15
    return dict(rows=rows,canonical_vacuum_mass_squared=masses.tolist(),
        expected_from_same_Higgs_and_F=expected.tolist(),
        old_connection_normalization_preserved=True)


def gauge_covariance_check():
    phi=np.array([.11,.68,-.07,.13,.52])
    _,frame,_,_,_,J=target_connection(phi)
    rng=np.random.default_rng(58222)
    X=rng.normal(size=(4,5))*.2
    f=rng.normal(size=(4,4,4))*.13;f=f-f.transpose(0,2,1)
    # A finite weak rotation, implemented through the original complex generator.
    angles=np.array([.21,-.17,.09]);H=np.einsum('a,aij->ij',angles,original.lattice.T)
    eig,vec=np.linalg.eigh(H);unitary=(vec*np.exp(1j*eig))@vec.conj().T
    R=real_generator(unitary);R[4,4]=1
    rotated_phi=R@phi
    _,newframe,newinverse,_,_,newJ=target_connection(rotated_phi)
    coordinate=X@frame.T
    newX=(coordinate@R.T)@newinverse.T
    newf=np.zeros_like(f)
    for mu in range(4):
        for nu in range(4):
            G=np.einsum('A,Aab->ab',f[:,mu,nu],GENERATORS)
            rotated=R@G@R.T
            newf[:,mu,nu]=-np.einsum('Aab,ba->A',GENERATORS,rotated)/INDEX
    before=connection_correction(X,f,J);after=connection_correction(newX,newf,newJ)
    err=max(abs(before[k]-after[k]) for k in before)
    assert err<1e-13
    assert np.max(abs(newframe.T@original.metric(rotated_phi)@newframe-np.eye(5)))<1e-14
    return dict(before=before,after=after,max_gauge_error=float(err),
                transformation_in_original_SU2_not_arbitrary_target_rotation=True)


def gauge_stress_reduction_check():
    rng=np.random.default_rng(58233)
    phi=np.array([.11,.68,-.07,.13,.52]);_,_,_,_,_,J=target_connection(phi)
    U=float(original.node_potential(phi))
    weights=np.array([original.lattice.PAR['K'][1]]*3+[original.lattice.PAR['K'][2]])
    rows=[]
    for case in range(3):
        X=rng.normal(size=(4,5))*.2;B=X@X.T;S=float(np.trace(B));Q=float(np.trace(B@B));M=X.T@X
        A=rng.normal(size=(5,5))*.11;A=(A+A.T)/2
        f=rng.normal(size=(4,4,4))*.011;f=f-f.transpose(0,2,1)
        tau=sum(weights[a]*(f[a]@f[a].T-np.eye(4)*np.sum(f[a]*f[a])/4) for a in range(4))
        assert abs(np.trace(tau))<1e-13
        correction=connection_correction(X,f,J)['correction']
        Ric=rng.normal(size=(4,4))*.1;Ric=(Ric+Ric.T)/2
        if case==2:Ric=B+U*np.eye(4)+tau
        R=float(np.trace(Ric));ric2=float(np.trace(Ric@Ric));riem2=1.1+2*ric2-R*R/3
        euler=riem2-4*ric2+R*R
        full=scalar.entry.invariant_coefficient(X,A,R)+correction+5*((riem2-ric2)/180+R*R/72)
        oldreduced=(euler/36-7*S*S/216+11*Q/108+U*S/18+U*U
                    -2*U*np.trace(A)/3+np.trace(A@A)/2-np.trace(A@M)/6)
        source_terms=float(np.sum(B*tau)/6+np.sum(tau*tau)/12)
        reduced=oldreduced+correction+source_terms
        E=Ric-B-U*np.eye(4)-tau;e=float(np.trace(E))
        D=(Ric+B+U*np.eye(4)+tau)/12+np.eye(4)*((R+S+4*U)/24-np.trace(A)/6-S/9)
        T=-2*D+np.trace(D)*np.eye(4);H=-(E-e*np.eye(4)/2)/2
        error=float(full-reduced-np.sum(H*T))
        assert abs(error)<2e-13
        if case==2:
            assert abs(full-reduced)<2e-13
            assert abs(full-(oldreduced+correction))>1e-4
        rows.append(dict(case=case,metric_identity_error=error,source_cross_terms=source_terms,
                         norm_gauge_stress=float(np.linalg.norm(tau)),on_leading_metric_relation=case==2))
    return dict(gauge_weights=weights.tolist(),rows=rows,
                scalar_Euler_elimination_uses_gauge_invariance_of_original_U=True)


def local_connection_realization_check():
    phi0=np.array([.11,.68,-.07,.13,.52]);_,frame,inverse_frame,_,_,J=target_connection(phi0)
    rng=np.random.default_rng(58244);X=rng.normal(size=(4,5))*.17;V=X@frame.T
    f=np.zeros((4,4));f[0,1]=.23;f[1,0]=-.23
    T=GENERATORS[3]
    def connection(mu,point):
        phi=phi0+point@V;value=original.F(phi)
        a=-f@point/2
        nabla=T+np.outer(T@phi,phi)/(6*value)
        return (np.outer(V[mu],phi)+np.dot(V[mu],phi)*np.eye(5))/(6*value)+a[mu]*nabla
    zero=np.zeros(4);G=[connection(mu,zero) for mu in range(4)]
    rows=[]
    for step in (2e-4,1e-4,5e-5):
        worst=0.
        for mu in range(4):
            for nu in range(mu+1,4):
                em=np.eye(4)[mu]*step;en=np.eye(4)[nu]*step
                omega=(connection(nu,em)-connection(nu,-em)-connection(mu,en)+connection(mu,-en))/(2*step)+G[mu]@G[nu]-G[nu]@G[mu]
                orth=inverse_frame@omega@frame
                expected=-(np.outer(X[mu],X[nu])-np.outer(X[nu],X[mu]))/6+f[mu,nu]*J[3]
                worst=max(worst,float(np.max(abs(orth-expected))))
        assert worst<1e-8
        rows.append(dict(step=step,max_connection_commutator_error=worst))
    assert rows[-1]['max_connection_commutator_error']<rows[0]['max_connection_commutator_error']/4
    return dict(rows=rows,explicit_connection='a_mu(x)=-f_mu_nu*x_nu/2; phi(x)=phi0+V_mu*x_mu',
                only_a_local_off_shell_realization=True)


def run():
    initial=entry.run()
    assert initial==json.loads((HERE/'round582_drafts/gauged_target_entry_results.json').read_text('utf8'))
    evidence=dict(original_Abelian_connection=initial['killing_connection'],
        same_stress_sign_witness=dict(coefficient=initial['scalar_connection_coefficient'],stress=initial['classical_gauge_stress']),
        mass_connection_gram=mass_connection_gram_check(),gauge_covariance=gauge_covariance_check(),
        full_gauge_stress_reduction=gauge_stress_reduction_check(),
        local_connection_realization=local_connection_realization_check())
    deps=('joint_curved_quantum_source.py','joint_scalar_effective_geometry.py','research_note_569.md',
          'research_note_572.md','research_note_581.md','research_round_581_checks.json',
          'round582_drafts/gauged_target_entry_probe.py','round582_drafts/gauged_target_entry_results.json')
    return dict(round=582,tests_run=len(evidence),failures=0,errors=0,checks=list(evidence),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='same five-component H5 matter with original SU2 and charge-3 integer U1 connection; fixed-background scalar-loop local coefficient and first-order reduction using leading Einstein-scalar-YM action; explicit off-shell local jets and same-coupling gauge stress, not complete gauge/fermion/gravity loops, an independent operator basis, on-shell quantum source or graph-continuum matching')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
