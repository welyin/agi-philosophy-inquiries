"""821: original homogeneous completion difference and common quantum shifts.
Original background/constraints are used in the first check; finite oscillators
only calibrate the joint preparation identities, not an original W evaluation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'820'))
import color_all_mode_certificate as fixed
import controlled_shift_probe as early
previous=fixed.previous;geo=previous.geo
TARGET=HERE/'homogeneous_preparation_bridge_results.json'

def original_homogeneous_difference():
    d=previous.setup();q=d['q'];psi=d['psi'];tensor=d['tensor']
    A,E,F,C=fixed.probe.data()
    pm=previous.old.momentum(q,d['dp'],d['dE'],d['dE0'])
    bc,bw,b0=geo.old.PAR['b']
    pe=(psi**-12*np.sum(previous.old.kinverse(q,q['p'])*d['dp'],axis=-1)
        +2*psi**-8*(bw*np.sum(q['f']['E']*d['dE'],axis=(-1,-2))
                      +b0*np.sum(q['f']['E0']*d['dE0'],axis=-1)))
    ec_new=fixed.solve_color(-d['sigma_c'],-d['J']-pm,-(d['rho']+pe)*psi**8/(2*bc))
    dtensor,_=geo.solve_momentum(dict(q,mom=d['dmom']))
    dA=2*np.sum(tensor*dtensor,axis=(-1,-2))+2*np.sum(previous.old.kinverse(q,q['p'])*d['dp'],axis=-1)
    dY=(2*bw*np.sum(q['f']['E']*d['dE'],axis=(-1,-2))
        +2*b0*np.sum(q['f']['E0']*d['dE0'],axis=-1)+2*bc*np.sum(E*d['dEc'],axis=(-1,-2)))
    oldA=np.sum(tensor*tensor,axis=(-1,-2))+q['pKp']
    pot=5*q['C']*psi**4-q['B']+7*oldA*psi**-8+6*q['Y']*psi**-4
    op=lambda v:-8*geo.laplace(v)+pot*v
    rhs=dA*psi**-7+2*dY*psi**-3+2*d['rho']*psi**5
    k2=np.sum(geo.waves(q['N'])**2,axis=-1)
    v,_=geo.cg(op,rhs,lambda x:geo.ifft(geo.fft(x)/(8*k2+pot.mean())))
    hec=ec_new-d['dEc'];hpsi=-v;htensor=-dtensor
    dAh=2*np.sum(tensor*htensor,axis=(-1,-2))
    dYh=2*bc*np.sum(E*hec,axis=(-1,-2))
    errors=dict(color_Gauss=float(np.max(abs(sum(previous.D(hec[...,i,:],i,C) for i in range(3))))),
        momentum=float(np.max(abs(sum(geo.derivative(htensor[..., :,i],i) for i in range(3))
                                 +np.einsum('...ja,ija->...i',hec,F)))),
        Hamiltonian=float(np.max(abs(op(hpsi)-dAh*psi**-7-2*dYh*psi**-3))))
    assert max(errors.values())<2e-10
    assert np.max(abs(hpsi))>1e-3
    return dict(N=q['N'],homogeneous_constraint_residuals=errors,
        conformal_h_max=float(np.max(abs(hpsi))),color_momentum_h_max=float(np.max(abs(hec))),
        all_material_configuration_differences_zero=True,
        same_declared754_source_in_both_completions=True,
        actual819_logical_source_not_numerically_computed=True)

def joint_shift_and_ordering():
    N=14
    a=np.zeros((N,N),complex)
    for n in range(1,N):a[n-1,n]=np.sqrt(n)
    q=(a+a.conj().T)/np.sqrt(2);p=(a-a.conj().T)/(1j*np.sqrt(2))
    I,X,Y,Z=early.I,early.X,early.Y,early.Z
    vac=np.zeros((N,N));vac[0,0]=1
    L=[p,.7*q,.3*p+.2*q]
    V=sum((np.kron(s,l) for s,l in zip((X,Y,Z),L)),np.zeros((2*N,2*N),complex))
    C=np.array([[np.trace(vac@lr@ls) for ls in L] for lr in L])
    eig=np.linalg.eigvalsh(C);assert eig[0]>-1e-14
    assert np.max(abs(C.imag))>.1
    # Direct coefficient of U rho U* is checked against the ordered Gram,
    # including both positive and anticommutator terms.
    sigma=[X,Y,Z];max_error=0.;wrong_error=0.;mean_error=0.
    for rho in ((I+X)/2,(I+Z)/2,(I+.2*X-.4*Y+.1*Z)/2):
        state=np.kron(rho,vac)
        full=V@state@V-.5*(V@V@state+state@V@V)
        reduced=early.reduce_record(full,N)
        def dissipator(ordered):
            out=np.zeros((2,2),complex)
            for r in range(3):
                for s in range(3):
                    term=sigma[r]@rho@sigma[s]-.5*(sigma[s]@sigma[r]@rho+rho@sigma[s]@sigma[r])
                    out+=ordered[s,r]*term
            return out
        max_error=max(max_error,float(np.linalg.norm(reduced-dissipator(C))))
        wrong_error=max(wrong_error,float(np.linalg.norm(reduced-dissipator(C.T))))
        for observable,target in ((q,X+.3*Z),(p,-.7*Y-.2*Z)):
            b=np.kron(I,observable)
            exact=np.trace(state@(1j*(V@b-b@V))).real
            predicted=np.trace(rho@target).real
            mean_error=max(mean_error,abs(exact-predicted))
    assert max_error<1e-12 and mean_error<1e-12 and wrong_error>.1
    # Original process plus controlled preparation: keep shared covariance.
    B=[.4*p-.1*q,.2*p+.5*q,-.3*q]
    gram=lambda A:np.array([[np.trace(vac@r@s) for s in A] for r in A])
    total=gram([b+l for b,l in zip(B,L)])
    omitted=gram(B)+C
    cross_error=float(np.linalg.norm(total-omitted))
    assert cross_error>.1
    return dict(noncommuting_three_axis_preparation=True,
        source_mean_coefficient_error=mean_error,ordered_logical_coefficient_error=max_error,
        error_if_order_reversed=wrong_error,
        error_if_original_process_preparation_cross_covariance_omitted=cross_error,
        compensation_noise_eigenvalues=eig.tolist(),compensation_noise_trace=float(np.trace(C).real),
        original_physical_W_evaluated=False)

def run():
    return dict(round=821,all_checks_passed=True,
        original_homogeneous_difference=original_homogeneous_difference(),
        joint_preparation=joint_shift_and_ordering(),
        construction_scope='Original global free physical representation and its local formal interacting image.',
        original_generator_coefficients_computed=False,
        autonomous_local_preparation_proven=False,finite_coupling_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
