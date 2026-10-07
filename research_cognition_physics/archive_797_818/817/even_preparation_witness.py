"""817: parity-preserving finite preparation witnesses.
The finite CAR identities are checked independently. Original neutral mass
coefficients calibrate the source change; this does not simulate the original
constrained spacetime or an autonomous preparation device.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'815'))
import record_retention_criterion as previous
TARGET=HERE/'even_preparation_witness_results.json'

def data():
    a,xi,Q,cov=previous.system(5);eye=np.eye(32,dtype=complex)
    gamma=[]
    for c in a:gamma += [c+c.conj().T,1j*(c-c.conj().T)]
    parity=np.eye(32,dtype=complex)
    for c in a:parity=parity@(eye-2*c.conj().T@c)
    units=[];rotations=[];subsets=[]
    for mask in range(256):
        indices=[i for i in range(8) if mask>>i&1]
        if len(indices)%2:indices.append(8)
        U=eye.copy();R=np.eye(10,dtype=complex)
        for k in indices:
            U=U@gamma[k]
            e=np.zeros(10,complex);mode=k//2
            if k%2:e[mode]=1j/np.sqrt(2);e[mode+5]=-1j/np.sqrt(2)
            else:e[mode]=e[mode+5]=1/np.sqrt(2)
            R=R@(np.eye(10)-2*np.outer(e,e.conj()))
        assert np.linalg.norm(U.conj().T@U-eye)<1e-12
        assert np.linalg.norm(U@parity-parity@U)<1e-12
        units.append(U);rotations.append(R);subsets.append(mask.bit_count()%2==0)
    return a,xi,Q,cov,gamma,parity,units,rotations,subsets

def trace_mean(O):return np.trace(O)/len(O)

def retention(d):
    a,xi,Q,cov,gamma,parity,units,rotations,even_subsets=d
    eye=np.eye(32,dtype=complex)
    H=sum((a[i+2].conj().T@a[i]+a[i].conj().T@a[i+2] for i in (0,1)),np.zeros((32,32),complex))
    HI=.3j*(a[2].conj().T@a[0]-a[0].conj().T@a[2])
    HI+=.2*(a[3].conj().T@a[0]+a[0].conj().T@a[3])
    K=H+1j*HI  # W_B(0,a)=(1,i) for finite X,Y boson diagnostic.
    basis=[1j*gamma[i]@gamma[j] for i in range(4) for j in range(i+1,4)]
    basis.append(gamma[0]@gamma[1]@gamma[2]@gamma[3])
    response=[1j*(A@K-K@A) for A in basis]
    G=np.array([[trace_mean(X.conj().T@Y).real for Y in response] for X in response])
    eig=np.linalg.eigvalsh(G)
    assert eig.min()>1.
    vacuum=np.eye(32)[:,0];rho=np.zeros((32,32),complex);naive=rho.copy()
    for U,keep in zip(units,even_subsets):
        v=U@vacuum
        rho+=np.outer(v,v.conj())/len(units)
        if keep:naive+=np.outer(v,v.conj())/sum(even_subsets)
    # Even preparations preserve TOTAL parity, while both local parities occur.
    local_parity=eye.copy()
    for c in a[:4]:local_parity=local_parity@(eye-2*c.conj().T@c)
    bad=(eye-local_parity)/2
    assert abs(np.trace(naive@bad))<1e-12
    assert abs(np.trace(rho@bad)-.5)<1e-12
    assert abs(np.trace(rho@parity)-1)<1e-12
    # Check all 128 even monomials, not just the selected response.
    twirl_error=0.
    for r in (0,2,4,6,8):
        for indices in itertools.combinations(range(8),r):
            O=eye.copy()
            for i in indices:O=O@gamma[i]
            twirl_error=max(twirl_error,float(abs(np.trace(rho@O)-trace_mean(O))))
    assert twirl_error<1e-12
    Z=eye-2*a[0].conj().T@a[0]
    X=1j*(Z@K-K@Z);square=X.conj().T@X
    reference=float(np.vdot(vacuum,square@vacuum).real)
    mixed=float(np.trace(rho@square).real)
    witness_values=[float(np.vdot(U@vacuum,square@(U@vacuum)).real) for U in units]
    assert reference<1e-12 and mixed>1.
    assert abs(mixed-trace_mean(square))<1e-12
    Xb=np.array([[0,1],[1,0]],complex);Yb=np.array([[0,-1j],[1j,0]],complex)
    D=np.kron(Xb,1j*(Z@H-H@Z))+np.kron(Yb,1j*(Z@HI-HI@Z))
    boson=np.diag([1.,0.]);actual=float(np.trace(np.kron(boson,rho)@D@D).real)
    assert abs(actual-mixed)<1e-12
    return dict(preparation_branches=256,all_branches_total_parity_even=True,
        even_operator_trace_identity_max_error=twirl_error,
        local_odd_sector_weight_naive_even_twirl=float(np.trace(naive@bad).real),
        local_odd_sector_weight_with_spectator=float(np.trace(rho@bad).real),
        traceless_internal_even_dimension=7,
        real_response_Gram_eigenvalues=eig.tolist(),
        reference_Z_change_square=reference,prepared_Z_change_square=mixed,
        maximum_branch_Z_change_square=max(witness_values),
        boson_contraction_square_error=abs(actual-mixed),
        original_continuous_response_numerically_evaluated=False)

def source_check(d):
    a,xi,Q,cov,gamma,parity,units,rotations,even_subsets=d
    ids=[24,25,30,31,56,57,62,63];dest=[0,1,2,3,5,6,7,8]
    phi=np.array([0.,.6987151138647895,0.,0.,.545863690226873])
    def embed(M):
        k=np.zeros((10,10),complex);k[np.ix_(dest,dest)]=M[np.ix_(ids,ids)];return k
    mass=embed(previous.old.vertex.mass(phi));H=Q(mass)
    # Fix the spectator vacuum and choose the lowest even vector. This is
    # explicitly a finite diagnostic preparation, not the original continuum P.
    eye=np.eye(32);number=a[4].conj().T@a[4]
    allowed=np.where((np.diag(parity).real>.5)&(np.diag(number).real<.5))[0]
    es,vs=np.linalg.eigh(H[np.ix_(allowed,allowed)])
    v=np.zeros(32,complex);v[allowed]=vs[:,0];P=cov(v)
    post=np.zeros((32,32),complex);P_rot=np.zeros((10,10),complex)
    covariance_error=0.
    for i,(U,R) in enumerate(zip(units,rotations)):
        w=U@v;post+=np.outer(w,w.conj())/256;P_rot+=R@P@R.conj().T/256
        if i in (0,1,2,3,7,31,255):
            covariance_error=max(covariance_error,float(np.max(abs(cov(w)-R@P@R.conj().T))))
    P_post=np.array([[np.trace(post@x@y.conj().T) for y in xi] for x in xi])
    covariance_error=max(covariance_error,float(np.max(abs(P_post-P_rot))))
    assert covariance_error<1e-12
    delta=P_post-P
    C=np.block([[np.zeros((5,5)),np.eye(5)],[np.eye(5),np.zeros((5,5))]])
    selfdual=float(np.max(abs(C@P_post.conj()@C-(np.eye(10)-P_post))))
    assert selfdual<1e-12 and np.linalg.eigvalsh(P_post).min()>-1e-12
    changes=[]
    for index in (1,4):
        dm=embed(previous.old.vertex.dmass(phi,np.eye(5)[index]))
        action_source=-Q(dm)
        actual=float((np.trace(post@action_source)-np.vdot(v,action_source@v)).real)
        expected=float((.5*np.trace(delta@dm)).real)
        assert abs(actual-expected)<1e-12
        changes.append(dict(scalar_index=index,action_mass_source_change=actual,
            covariance_contraction_change=expected))
    energy=float((np.trace(post@H)-np.vdot(v,H@v)).real)
    assert energy>1e-3
    return dict(original_neutral_mass_coefficients_used=True,
        finite_diagnostic_preparation_explicitly_changed=True,
        branch_and_mixture_covariance_max_error=covariance_error,
        covariance_self_duality_error=selfdual,
        covariance_difference_rank=int(np.linalg.matrix_rank(delta,tol=1e-11)),
        original_auxiliary_H_energy_change=energy,scalar_source_changes=changes,
        all_local_spacetime_sources_computed=False,
        preparation_device_or_compatible_geometry_simulated=False)

def run():
    d=data()
    return dict(round=817,all_checks_passed=True,retention=retention(d),source=source_check(d),
        scope='Parity-preserving finite preparation ensemble; original continuum connection is analytic and does not supply an autonomous device.',
        original_reference_state_replaced_silently=False,
        information_loss_or_recovery_impossibility_proven=False,
        finite_coupling_retention_time_computed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

