"""926: distinguish rotation, role exchange, quadratic dispersion and task-metric additivity.
Small matrices verify explicit classifications and counterexamples. Physical coefficients remain inputs.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.]).astype(complex)
PAULI=(X,Y,Z)
def op(a):return float(np.linalg.norm(a,2))
def com(a,b):return a@b-b@a
def anti(a,b):return a@b+b@a
def direct(a,b):
    return np.block([[a,np.zeros((len(a),len(b)),complex)],[np.zeros((len(b),len(a)),complex),b]])
def hermbasis(d):
    out=[]
    for i in range(d):
        a=np.zeros((d,d),complex);a[i,i]=1;out.append(a)
    for i in range(d):
        for j in range(i+1,d):
            a=np.zeros((d,d),complex);a[i,j]=a[j,i]=1;out.append(a)
            a=np.zeros((d,d),complex);a[i,j]=1j;a[j,i]=-1j;out.append(a)
    return out

def linear_dimension(basis,maps):
    cols=[]
    for b in basis:
        col=np.concatenate([np.asarray(fn(b)).reshape(-1) for fn in maps])
        cols.append(np.r_[col.real,col.imag])
    s=np.linalg.svd(np.column_stack(cols),compute_uv=False)
    return len(basis)-int(np.count_nonzero(s>1e-10))

def directions(m,n):
    alpha=[direct(np.kron(s,np.eye(m)),-np.kron(s,np.eye(n))) for s in PAULI]
    spin=[direct(np.kron(s,np.eye(m)),np.kron(s,np.eye(n)))/2 for s in PAULI]
    return alpha,spin

def mass(k):
    m,n=k.shape;b=np.kron(I,k)
    return np.block([[np.zeros((2*m,2*m),complex),b],[b.conj().T,np.zeros((2*n,2*n),complex)]])

def stabilizer_dimension(k):
    m,n=k.shape
    basis=[(a,np.zeros((n,n),complex)) for a in hermbasis(m)]+[(np.zeros((m,m),complex),b) for b in hermbasis(n)]
    cols=[]
    for a,b in basis:
        c=(a@k-k@b).reshape(-1);cols.append(np.r_[c.real,c.imag])
    return len(basis)-int(np.count_nonzero(np.linalg.svd(np.column_stack(cols),compute_uv=False)>1e-10))

def expi(a):
    e,u=np.linalg.eigh(a);return (u*np.exp(1j*e))@u.conj().T

def mean(psi,a):return float(np.vdot(psi,a@psi).real)
def cov(psi,a,b):return mean(psi,anti(a,b)/2)-mean(psi,a)*mean(psi,b)
def var(psi,a):return mean(psi,a@a)-mean(psi,a)**2

def run():
    rng=np.random.default_rng(926)
    classification=[]
    for m,n in ((1,1),(2,2),(3,2)):
        alpha,spin=directions(m,n);basis=hermbasis(2*(m+n))
        r_maps=[lambda a,s=s:com(a,s) for s in spin]
        q_maps=[lambda a,s=s:anti(a,s) for s in alpha]
        dr=linear_dimension(basis,r_maps);dq=linear_dimension(basis,r_maps+q_maps)
        assert dr==(m+n)**2 and dq==2*m*n
        row=dict(positive_role_multiplicity=m,negative_role_multiplicity=n,rotation_dimension=dr,quadratic_separation_dimension=dq)
        if m==n:
            parity=np.kron(X,np.eye(2*m));dp=linear_dimension(basis,r_maps+[lambda a:com(a,parity)])
            dqp=linear_dimension(basis,r_maps+q_maps+[lambda a:com(a,parity)])
            assert dp==2*m*m and dqp==m*m
            row.update(rotation_and_role_exchange_dimension=dp,with_quadratic_separation_and_role_exchange_dimension=dqp)
        classification.append(row)
    # Rotational and parity invariant, both roles mix, but energy square does not separate.
    a=.4;b=.7;p=.6
    alpha,spin=directions(2,2)
    diagonal=np.kron(np.eye(4),a*Z)
    mix=np.kron(np.kron(X,I),b*X)
    M=diagonal+mix;P=np.kron(X,np.eye(4));H=p*alpha[2]+M
    err_rot=max(op(com(M,s)) for s in spin);err_parity=op(com(M,P))
    square_defect=op(H@H-p*p*np.eye(8)-M@M)
    ehi=np.sqrt((p+a)**2+b*b);elo=np.sqrt((p-a)**2+b*b)
    expected=np.array([-ehi,-ehi,-elo,-elo,elo,elo,ehi,ehi])
    spectrum_error=float(np.max(np.abs(np.linalg.eigvalsh(H)-expected)))
    assert err_rot==err_parity==0 and abs(square_defect-2*a*p)<1e-14 and spectrum_error<1e-14
    assert op(P@H@P-(-p*alpha[2]+M))<1e-14
    assert op(M@M-(a*a+b*b)*np.eye(8))<1e-14
    # One-dimensional internal roles are the positive special case after one common energy shift.
    aa,ss=directions(1,1);m0=.7;center=.3
    minM=center*np.eye(4)+np.kron(X,m0*I);hc=p*aa[2]+minM-center*np.eye(4)
    minimal_error=op(hc@hc-(p*p+m0*m0)*np.eye(4));assert minimal_error<1e-14
    # General rectangular mixing: one K controls gap, unmatched role count and stabilizer.
    cases=[]
    for name,k,want_stab in (
        ('balanced_distinct',np.diag([.4,.9]).astype(complex),2),
        ('balanced_degenerate',.7*np.eye(2,dtype=complex),4),
        ('unbalanced_full_rank',np.array([[.4,0],[0,.9],[0,0]],complex),3),
        ('balanced_rank_one',np.diag([.4,0]).astype(complex),3)):
        m,n=k.shape;A,_=directions(m,n);mm=mass(k);dim=len(mm)
        rank=int(np.linalg.matrix_rank(k,tol=1e-12));nullity=int(np.count_nonzero(abs(np.linalg.eigvalsh(mm))<1e-12))
        assert nullity==2*(m+n-2*rank)
        pp=np.array([.17,-.29,.31]);hh=sum(pp[i]*A[i] for i in range(3))+mm
        error=op(hh@hh-np.dot(pp,pp)*np.eye(dim)-mm@mm)
        stab=stabilizer_dimension(k);assert error<1e-14 and stab==want_stab
        u=expi(np.diag(np.arange(m)*.17));v=expi(np.ones((n,n))*.11)
        transform=direct(np.kron(I,u),np.kron(I,v));kk=u@k@v.conj().T
        covariance_error=op(mass(kk)-transform@mm@transform.conj().T)
        assert covariance_error<1e-14
        cases.append(dict(name=name,rank=rank,mass_nullity=nullity,unpaired_role_index=m-n,
            nonzero_singular_values=[float(x) for x in np.linalg.svd(k,compute_uv=False) if x>1e-12],
            internal_stabilizer_lie_dimension=stab,square_identity_error=error,comparison_covariance_error=covariance_error))
    # Pure-state distinguishability metric: variance additivity is not uncentred second-moment additivity.
    A=np.kron(Z,Z);B=np.kron(X,I)
    theta=np.pi/4
    psi=np.kron(np.array([np.cos(theta/2),np.sin(theta/2)]),np.array([1.,0.])).astype(complex)
    witness=dict(anticommutator_norm=op(anti(A,B)),mean_A=mean(psi,A),mean_B=mean(psi,B),
        covariance=cov(psi,A,B),variance_A=var(psi,A),variance_B=var(psi,B),variance_sum=var(psi,A+B),
        uncentred_second_moment_additivity_error=abs(mean(psi,(A+B)@(A+B))-mean(psi,A@A)-mean(psi,B@B)))
    assert abs(witness['covariance']+.5)<1e-14 and abs(witness['variance_sum'])<1e-14
    assert witness['anticommutator_norm']==0 and witness['uncentred_second_moment_additivity_error']<1e-14
    # Numerical calibration of the general all-pure-state theorem; proof is in the note.
    pure_constraints=[]
    for _ in range(48):
        x=rng.normal(size=4)+1j*rng.normal(size=4);x/=np.linalg.norm(x)
        pure_constraints.append(lambda b,x=x:np.array([cov(x,A,b)]))
    var_dimension=linear_dimension(hermbasis(4),pure_constraints)
    assert var_dimension==1
    # A maximally mixed-only test is too weak even for the rotation/parity counterexample.
    tr_cross=[float(abs(np.trace(anti(ai,M)/2)/8-np.trace(ai)*np.trace(M)/64)) for ai in alpha]
    assert max(tr_cross)<1e-14
    return dict(round=926,date='2026-10-06',all_scientific_checks_passed=True,
        finite_operator_classification=classification,
        rotational_role_exchange_counterexample=dict(a=a,mixing=b,momentum=p,rotation_error=err_rot,role_exchange_error=err_parity,
            rest_square=(a*a+b*b),square_separation_defect=square_defect,positive_energies=[elo,ehi],
            analytic_spectrum_error=spectrum_error,mean_reference_cross_terms=tr_cross),
        minimal_balanced_role_case_square_error=minimal_error,mixing_gap_and_stabilizer_cases=cases,
        distinguishability_additivity_counterexample=witness,all_pure_variance_additivity_solution_dimension=var_dimension,
        square_separation_is_additional_input=True,variance_cost_not_identified_with_physical_energy_cost=True,
        classical_dirac_classification_not_claimed_new=True,full_quantum_unknown_state_contract_kept=True,
        all_matter_and_gauge_representations_selected=False,physical_mass_values_selected=False,
        Lorentz_covariance_derived_from_cognition=False,GR_generated=False,whole_stage_completed=False,full_goal_completed=False,
        source_hashes={str(Path('926')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();data=run()
    target=HERE/'direction_update_selection_results.json'
    if a.write:
        assert not target.exists();target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in data.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
