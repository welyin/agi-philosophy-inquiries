"""652: original smooth material references, common energy and form compression.

The full-Gauss self-adjointness, energy and Ritz statements are analytic.
Finite matrices reuse the exact 625 radial diagnostic, not a full graph/CAR
spectrum. The matrix reference tests only the normal derivative (W_f=0).
No detector design, continuum limit or quantum Einstein solution is tested.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_source_preserving_compression as old625

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_quantum_reference_forms_results.json'
M=original.M
HBAR=old625.HBAR
VOLUME=old625.VOLUME


def smooth_data(phi,which):
    inv=original.inverse(phi)
    f=original.F(phi)
    grad=np.zeros(5)
    if which=='T':
        grad[:4]=phi[:4]
        lap=f*(4-np.sum(phi[:4]**2)/(2*M))
    else:
        grad[4]=1
        lap=-f*phi[4]/(3*M)
    return inv@grad, float(grad@inv@grad), float(lap)


def calculus_check():
    rng=np.random.default_rng(652)
    points=[np.zeros(5),np.array([0.,0.,0.,0.,.7])]
    for frac in (.15,.5,.85,.96):
        for _ in range(6):
            v=rng.normal(size=5);points.append(v/np.linalg.norm(v)*np.sqrt(6*M)*frac)
    rows=[]
    for which in ('T','s'):
        bound_g=8*M*M/9 if which=='T' else M
        bound_d=4*M if which=='T' else np.sqrt(6*M)/3
        errors=[]
        for step in (4e-5,2e-5,1e-5):
            error=0.
            for phi in points:
                x,g,d=smooth_data(phi,which)
                # Independent Cartesian divergence with original full measure.
                measure=np.sqrt(M)/original.F(phi)**3
                div=0.
                for j in range(5):
                    delta=np.eye(5)[j]*step
                    pp=phi+delta;pm=phi-delta
                    vp=smooth_data(pp,which)[0];vm=smooth_data(pm,which)[0]
                    div+=(np.sqrt(M)*vp[j]/original.F(pp)**3
                         -np.sqrt(M)*vm[j]/original.F(pm)**3)/(2*step*measure)
                error=max(error,abs(div-d)/(1+abs(d)))
                assert g<=bound_g+1e-12 and abs(d)<=bound_d+1e-12
            errors.append(float(error))
        assert errors[-1]<1e-7 and errors[-1]<errors[0]/10
        rows.append(dict(reference=which,gradient_squared_bound=bound_g,
                         laplacian_bound=bound_d,divergence_errors=errors))
    # Old classical chart is unchanged up to a regular coordinate map at h>0.
    p=np.array([.7,.54,-.012,.031])  # h,s,A_h,B_s: algebraic chart test only.
    def chart(z):return np.array([z[0]**2/2,z[1],z[0]**2*z[2],z[3]])
    eps=1e-5
    jac=np.column_stack([(chart(p+eps*e)-chart(p-eps*e))/(2*eps) for e in np.eye(4)])
    determinant=float(np.linalg.det(jac))
    assert abs(determinant-p[0]**3)<1e-10
    return dict(points=len(points),rows=rows,chart_determinant=determinant,
                expected_determinant=float(p[0]**3),smooth_at_zero=True)


def diagnostic():
    # Same 64-dimensional quadratic-form diagnostic as 625, with its basis
    # retained so original T=h^2/2 and s multipliers can be represented.
    n=8
    hs=np.linspace(.15,1.55,n+2)[1:-1]
    ss=np.linspace(-1.1,1.1,n+2)[1:-1]
    h,s=np.meshgrid(hs,ss,indexing='ij');h=h.ravel();s=s.ravel()
    f=M-(h*h+s*s)/6
    measure=np.sqrt(M)*h**3/f**3
    rinv=1/np.sqrt(measure)
    def derivative(step):return (np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/(2*step)
    dh=np.kron(derivative(hs[1]-hs[0]),np.eye(n))*rinv[None,:]
    ds=np.kron(np.eye(n),derivative(ss[1]-ss[0]))*rinv[None,:]
    ghh=f*(1-h*h/(6*M));gss=f*(1-s*s/(6*M));ghs=-f*h*s/(6*M)
    kinetic=HBAR**2/(2*VOLUME)*(dh.T@((measure*ghh)[:,None]*dh)
        +ds.T@((measure*gss)[:,None]*ds)+dh.T@((measure*ghs)[:,None]*ds)
        +ds.T@((measure*ghs)[:,None]*dh))
    phi=np.zeros((n*n,5));phi[:,1]=h;phi[:,4]=s
    potential=np.diag(VOLUME*original.node_potential(phi))
    e,u=np.linalg.eigh(kinetic+potential)
    refs=[u.T@(a[:,None]*u) for a in (h*h/2,s)]
    velocities=[1j/HBAR*(e[:,None]-e[None,:])*a for a in refs]
    probabilities=np.exp(-old625.BETA*(e-e[0]));probabilities/=probabilities.sum()
    inherited=old625.diagnostic()
    error=float(np.max(abs(e-inherited['E'])))
    assert error<1e-12
    return dict(E=e,refs=refs,V=velocities,prob=probabilities,dimension=n*n,
                kinetic=kinetic,potential=potential,old625_spectrum_error=error)


def compression_check(d):
    rows=[]
    for label,v in zip(('T','s'),d['V']):
        assert np.max(abs(v-v.conj().T))<1e-12
        squared=v.conj().T@v
        for rank in (6,12,24,40):
            small=v[:rank,:rank]
            correct=squared[:rank,:rank]
            wrong=small.conj().T@small
            tail=v[rank:,:rank].conj().T@v[rank:,:rank]
            error=float(np.max(abs(correct-wrong-tail)))
            minimum=float(np.linalg.eigvalsh(tail).min())
            assert error<1e-11 and minimum>-1e-11
            step=1e-5
            # w(u)=exp(6u)w: original normal reference Q_f(u)=-exp(-12u)V_f^2.
            fd=(-np.exp(-12*step)*correct+np.exp(12*step)*correct)/(2*step)
            derivative_error=float(np.linalg.norm(fd-12*correct)/(1+np.linalg.norm(12*correct)))
            assert derivative_error<3e-9
            rows.append(dict(reference=label,rank=rank,tail_norm=float(np.linalg.norm(tail,2)),
                tail_minimum_eigenvalue=minimum,tail_identity_error=error,
                source_derivative_relative_error=derivative_error))
    assert all(row['tail_norm']>.01 for row in rows)
    return dict(old625_spectrum_error=d['old625_spectrum_error'],rows=rows,
                omitted_tail_changes_reference_and_source=True)


def common_state_check(d):
    n=d['dimension'];prob=d['prob'];energy=d['E']
    rfull=[v.conj().T@v for v in d['V']]
    kfull=[np.linalg.inv(np.eye(n)+r) for r in rfull]
    phase=np.exp(-1j*.17*energy/HBAR)
    def word(k):return (phase.conj()[:,None]*k[0]*phase[None,:])@k[1]
    target=complex(prob@np.diag(word(kfull)))
    means=[float(np.real(prob@np.diag(r))) for r in rfull]
    rows=[]
    rhs=np.exp(-np.arange(n)/9)*(1+1j*np.sin(np.arange(n)))
    rhs/=np.linalg.norm(rhs)
    solutions=[k@rhs for k in kfull]
    for rank in (6,12,24,40,56,64):
        q=prob[:rank].copy();q[0]+=prob[rank:].sum() # same 625 state map
        ks=[];residuals=[];ritz_errors=[]
        for r,u in zip(rfull,solutions):
            kr=np.linalg.inv(np.eye(rank)+r[:rank,:rank])
            lifted=np.zeros((n,n),complex);lifted[:rank,:rank]=kr;ks.append(lifted)
            ur=lifted@rhs
            residuals.append(float(np.linalg.norm(((np.eye(n)+r)@ur-rhs)[:rank])))
            delta=u-ur
            ritz_errors.append(float(np.real(delta.conj()@(np.eye(n)+r)@delta)))
            assert np.linalg.eigvalsh(kr).min()>0 and np.linalg.eigvalsh(kr).max()<=1+1e-11
        value=complex(q@np.diag(word(ks))[:rank])
        mean_r=[float(np.real(q@np.diag(r)[:rank])) for r in rfull]
        rows.append(dict(rank=rank,resolvent_vector_errors=[float(np.linalg.norm(k@rhs-u)) for k,u in zip(ks,solutions)],
            ritz_energy_errors=ritz_errors,variational_residuals=residuals,
            ordered_word=[value.real,value.imag],ordered_word_error=float(abs(value-target)),
            velocity_square_means=mean_r,velocity_square_mean_errors=[abs(a-b) for a,b in zip(mean_r,means)]))
    for j in (0,1):
        assert all(b['ritz_energy_errors'][j]<=a['ritz_energy_errors'][j]+1e-10 for a,b in zip(rows,rows[1:]))
    assert max(rows[-1]['resolvent_vector_errors'])<1e-11 and rows[-1]['ordered_word_error']<1e-11
    assert rows[0]['ordered_word_error']>1e-6
    return dict(full_velocity_square_means=means,full_ordered_word=[target.real,target.imag],rows=rows,
        maximum_variational_residual=max(max(r['variational_residuals']) for r in rows),
        scope='Fixed 64-dimensional original radial diagnostic; no infinite-model numerical error bound.')


def run():
    calculus=calculus_check();d=diagnostic()
    compression=compression_check(d);state=common_state_check(d)
    names=('research_note_574.md','research_note_575.md','research_note_588.md','research_note_598.md',
           'research_note_603.md','research_note_617.md','research_note_625.md','research_note_649.md',
           'research_note_651.md','joint_curved_quantum_source.py','joint_source_preserving_compression.py')
    return dict(date='2026-10-02',round=652,tests_run=3,failures=0,errors=0,
        calculus=calculus,compression=compression,common_state=state,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
        scope='Analytic original finite-graph Gauss reference forms, shared energy and region transport; numerical full-target calculus and old 64D normal-sector diagnostic. No joint sharp coordinates, quantum GR or continuum completion.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
