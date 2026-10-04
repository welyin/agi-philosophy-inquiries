"""609: functorial Fock completion of the actual 606/607 flat Wilson fibre.

The four-mode mass is a declared normal/pair diagnostic, not the SM mass map.
The analytic construction treats general finite-rank smooth projectors.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_chiral_fibre_source as prior
from joint_projection_local_composition import unitary, norm, comm

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_fock_covariant_completion_results.json'
C=prior.annihilators(4)
def dg(b):
    return sum(b[i,j]*C[i].conj().T@C[j] for i in range(4) for j in range(4))
def mass(b,delta):
    pair=.5*sum(delta[i,j]*C[i].conj().T@C[j].conj().T for i in range(4) for j in range(4))
    return dg(b)+pair+pair.conj().T
def exterior(U):
    # Columns and rows use the occupation-bit order of annihilators().
    m=U.shape[0];out=np.zeros((2**m,2**m),complex)
    for r in range(2**m):
        rr=[i for i in range(m) if r>>i&1]
        for c in range(2**m):
            cc=[i for i in range(m) if c>>i&1]
            if len(rr)==len(cc):
                out[r,c]=np.linalg.det(U[np.ix_(rr,cc)]) if rr else 1
    return out
def data():
    e,v=np.linalg.eigh(-prior.G5)
    W0=np.column_stack([v[:,e>0],v[:,e<0]])
    K=W0.conj().T@(1j*prior.GAMMA[0]/2)@W0
    p=np.diag([1,1,0,0]).astype(complex);dp=comm(K,p)
    J=np.eye(16,dtype=complex)[:,:4];P=J@J.conj().T
    nq=dg(np.eye(4)-p)
    return W0,K,p,dp,J,P,nq

def connection_check():
    W0,K,p,dp,J,P,nq=data();KF=dg(K)
    C1=comm(p,dp);CF=dg(C1);dPF=comm(KF,P)
    globalC=comm(P,dPF);globalV=dPF@dPF;localV=dg(dp@dp)
    assert norm(C1+K)<1e-12
    assert norm(dPF+comm(CF,P))<1e-12
    assert norm(comm(KF,nq)+comm(CF,nq))<1e-12
    phi=J.conj().T@localV@J
    assert norm(phi-J.conj().T@KF.conj().T@KF@J)<1e-12
    assert norm((globalC-CF)@J)<1e-12
    assert norm(J.conj().T@(globalV-localV)@J)<1e-12
    assert norm(globalC-CF)>.1 and norm(globalV-localV)>.1
    eigen=np.linalg.eigvalsh(phi)
    assert np.max(abs(eigen-np.array([0,.25,.25,.5])))<1e-12
    rng=np.random.default_rng(609)
    f=rng.normal(size=4)+1j*rng.normal(size=4)
    df=rng.normal(size=4)+1j*rng.normal(size=4)
    psi=J@f;dpsi=KF@J@f+J@df
    old=float(np.linalg.norm(dpsi)**2)
    new=float(np.linalg.norm(dpsi+CF@psi)**2+np.vdot(psi,localV@psi).real)
    assert abs(old-new)<1e-12
    F=np.kron(np.ones((16,1))/4,np.eye(4));actual_error=0.
    for theta in (-.2,0,.3):
        full=prior.projector(theta,1)[0]
        U=np.cos(theta/2)*np.eye(4)+1j*np.sin(theta/2)*prior.GAMMA[0]
        expected=U@W0@p@W0.conj().T@U.conj().T
        actual_error=max(actual_error,norm(F.conj().T@full@F-expected))
        assert norm(exterior(U).conj().T@exterior(U)-np.eye(16))<1e-12
    assert actual_error<1e-12
    return dict(actual_Wilson_projector_error=actual_error,
        allowed_form_error=abs(old-new),allowed_geometric_eigenvalues=eigen.tolist(),
        full_connection_difference=norm(globalC-CF),
        full_positive_term_difference=norm(globalV-localV),
        covariant_constraint_error=norm(dPF+comm(CF,P)),
        full_positive_term_minimum=float(np.linalg.eigvalsh(localV)[0]))

def masses():
    b=np.array([[.2,.13+.07j,.19,0],[.13-.07j,-.1,.11j,.17],
                [.19,-.11j,.12,.08],[0,.17,.08,-.16]],complex)
    delta=np.zeros((4,4),complex)
    for i,j,z in ((0,1,.21+.08j),(0,2,.15),(1,3,.11j),(2,3,.27)):
        delta[i,j]=z;delta[j,i]=-z
    return b,delta
def mass_check():
    _,_,p,_,J,P,nq=data();q=np.eye(4)-p
    b,delta=masses();B=mass(b,delta)
    numbers=np.real(np.diag(nq));mask=abs(numbers[:,None]-numbers[None,:])<.1
    averaged=B*mask
    exact=mass(p@b@p+q@b@q,p@delta@p.T)
    wrong=mass(p@b@p+q@b@q,p@delta@p.T+q@delta@q.T)
    err=norm(averaged-exact)
    assert err<1e-12 and norm(comm(averaged,nq))<1e-12
    assert norm(J.conj().T@(averaged-B)@J)<1e-12
    leak=norm((np.eye(16)-P)@wrong@J);assert leak>.2
    assert norm((np.eye(16)-P)@averaged@J)<1e-12
    assert norm(averaged)<=norm(B)+1e-12
    # Explicit composition of two independent four-mode even mass systems.
    NB=np.kron(nq,np.eye(16))+np.kron(np.eye(16),nq)
    nval=np.real(np.diag(NB));keep=abs(nval[:,None]-nval[None,:])<.1
    total=np.kron(B,np.eye(16))+np.kron(np.eye(16),B)
    predicted=np.kron(averaged,np.eye(16))+np.kron(np.eye(16),averaged)
    composition=float(np.max(abs(total*keep-predicted)));assert composition<1e-12
    return dict(normal_pair_formula_error=err,original_mass_norm=norm(B),
        completed_mass_norm=norm(averaged),wrong_qq_pair_leakage=leak,
        completed_sector_leakage=norm((np.eye(16)-P)@averaged@J),
        independent_mass_composition_error=composition,
        mass_spectrum_change=float(np.max(abs(np.linalg.eigvalsh(averaged)-np.linalg.eigvalsh(B)))),
        diagnostic_mass_is_not_original_32_mode_SM_embedding=True)

def grid(N=7,gamma=.12):
    W0,K,p,dp,J,P,nq=data();kappa=np.exp(-2*gamma)
    step=2/(N+1);theta=-1+step*np.arange(1,N+1)
    T=(2*np.eye(N)-np.eye(N,k=1)-np.eye(N,k=-1))/step**2
    b,delta=masses();q=np.eye(4)-p
    B=mass(p@b@p+q@b@q,p@delta@p.T)
    kinetic=np.kron(T,np.eye(16))+np.kron(np.eye(N),dg(dp@dp))
    Hf=kappa*kinetic+np.kron(np.eye(N),B)
    W=np.zeros_like(Hf)
    for i,t in enumerate(theta):
        # This is the full exterior frame, not just the allowed columns.
        U=np.cos(t/2)*np.eye(4)+2*np.sin(t/2)*K
        W[16*i:16*i+16,16*i:16*i+16]=exterior(U)
    Jg=W@np.kron(np.eye(N),J)
    H=W@Hf@W.conj().T;G=W@(-2*kappa*kinetic)@W.conj().T
    return H,G,Jg,theta
def history_source_check():
    H,G,J,theta=grid();h=J.conj().T@H@J;g=J.conj().T@G@J
    intertwining=norm(H@J-J@h);assert intertwining<1e-11
    assert norm(G@J-J@g)<1e-11
    rng=np.random.default_rng(1609);v=rng.normal(size=(len(h),2))+1j*rng.normal(size=(len(h),2))
    v/=np.linalg.norm(v)
    instruments=[]
    for function in (np.cos,np.sin):
        eff=.5+.2*function(theta)
        instruments.append([np.kron(np.diag(np.sqrt(eff)),np.eye(16)),
                            np.kron(np.diag(np.sqrt(1-eff)),np.eye(16))])
    history_error=0.;state_error=0.;source_error=0.;total=0.
    for a,b in itertools.product(range(2),repeat=2):
        L1=instruments[0][a];L2=instruments[a][b]
        K=unitary(H,.05)@L2@unitary(H,.08+.02*a)@L1@unitary(H,.11)
        k=unitary(h,.05)@(J.conj().T@L2@J)@unitary(h,.08+.02*a)@(J.conj().T@L1@J)@unitary(h,.11)
        history_error=max(history_error,norm(K@J-J@k))
        f=K@J@v;s=J@k@v
        diff=np.outer(f.ravel(),f.ravel().conj())-np.outer(s.ravel(),s.ravel().conj())
        state_error+=float(np.sum(abs(np.linalg.eigvalsh(diff))))
        source_error=max(source_error,abs(float(np.trace(f.conj().T@G@f).real-np.trace((k@v).conj().T@g@(k@v)).real)))
        total+=float(np.linalg.norm(f)**2)
    assert history_error<1e-11 and state_error<1e-11 and source_error<1e-9 and abs(total-1)<1e-11
    eps=2e-6;fd=(grid(gamma=.12+eps)[0]-grid(gamma=.12-eps)[0])/(2*eps)
    derivative=norm(fd-G);assert derivative<1e-7
    L=instruments[1];l=[J.conj().T@x@J for x in L]
    injection=sum(x@H@x for x in L)-H
    small=sum(x@h@x for x in l)-h
    injerror=norm(J.conj().T@injection@J-small);assert injerror<1e-11
    return dict(grid_points=7,ambient_dimension=len(H),allowed_dimension=len(h),
        generator_intertwining_error=intertwining,operator_history_error=history_error,
        complete_record_trace_norm=state_error,total_probability=total,
        branch_source_error=source_error,geometric_derivative_error=derivative,
        energy_injection_mapping_error=injerror,
        mass_independent_of_gamma_not_included_in_minus_two_source=True,
        discretizes_new_covariant_form_not_old_grid_pinching=True)
def run():
    result=dict(connection=connection_check(),mass=mass_check(),history_and_sources=history_source_check())
    deps=('research_note_598.md','research_note_603.md','research_note_606.md',
          'research_note_607.md','research_note_608.md','joint_chiral_fibre_source.py',
          'joint_projection_local_composition.py')
    return dict(round=609,tests_run=3,failures=0,errors=0,**result,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(conditional_Fock_form_completion=True,independent_composition_preserved=True,
            allowed_original_form_and_histories_preserved=True,
            complementary_pair_creation_must_be_removed=True,
            fixed_graph_Gauss_heat_proof_conditional=True,
            no_full_SM_mass_embedding_GW_measure_or_unbounded_LR=True,
            no_dimension_continuum_or_GR_completion=True))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=609,tests=3,all_passed=True,connection=result['connection'],
                         mass=result['mass'],histories=result['history_and_sources'])))
