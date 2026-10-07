"""814: separate record change from reading of an existing record.
Exact-size matrix checks and original Nambu principal symbols only.
No original continuum covariance or autonomous write/hold protocol is simulated.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
TARGET=HERE/'write_read_split_results.json'
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
I=np.eye(2,dtype=complex)
comm=lambda a,b:a@b-b@a
anti=lambda a,b:a@b+b@a

def algebra_case(label,J,rho_f,B,phi,rho_b):
    p=(I+Z)/2;V=np.kron(B,J);pt=np.kron(I,p);phit=np.kron(phi,I)
    state=np.kron(rho_b,rho_f);exp=lambda a:np.trace(state@a)
    mean=float(np.trace(rho_f@p).real)
    Dp=1j*comm(pt,V);Dphi=1j*comm(phit,V)
    write=exp(anti(phit,Dp)/2).real
    read=exp(anti(Dphi,pt-mean*np.eye(4))/2).real
    cov=np.trace(rho_f@p@J)-np.trace(rho_f@p)*np.trace(rho_f@J)
    w=np.trace(rho_b@phi@B)
    predicted_write=-2*w.real*cov.imag
    predicted_read=-2*w.imag*cov.real
    together=1j*exp(comm(phit@pt,V))-mean*exp(Dphi)
    square=float(exp(Dp@Dp).real)
    variance=float(exp(phit@phit).real)
    ordered_write=exp(phit@Dp);ordered_read=exp(Dphi@(pt-mean*np.eye(4)))
    errors=dict(write=float(abs(write-predicted_write)),read=float(abs(read-predicted_read)),
        total=float(abs(together-write-read)),
        imaginary_cancellation=float(abs((ordered_write+ordered_read).imag)),
        positive_Cauchy_bound=float(max(0,write**2-variance*square)))
    assert max(errors.values())<1e-12
    eps=.01;ev,vec=np.linalg.eigh(V);S=(vec*np.exp(1j*eps*ev))@vec.conj().T
    pout=S.conj().T@pt@S;q=np.eye(4)-pt
    disagreement=pt@(np.eye(4)-pout)@pt+q@pout@q
    identity_error=float(np.max(abs(disagreement-(pout-pt)@(pout-pt))))
    assert identity_error<1e-12
    return dict(label=label,write=float(write),read=float(read),total=float(together.real),
        ordered_write_imag=float(ordered_write.imag),ordered_read_imag=float(ordered_read.imag),
        record_change_operator_norm=float(np.linalg.norm(Dp,2)),record_change_square=square,
        symmetrized_witness_lower_bound=float(write**2/variance),
        errors=errors,formal_disagreement_identity_error=identity_error)

def principal_joint_witness():
    idx=[30,31];gam=[g[np.ix_(idx,idx)] for g in original.GAMMA]
    gr=lambda n:sum((g*x for g,x in zip(gam,n)),np.zeros((2,2),complex))
    r=np.array([[.4,.3,-.2],[.3,-.1,.25],[-.2,.25,-.3]])
    s=np.array([[.2,-.15,.4],[-.15,.05,.1],[.4,.1,-.25]])
    directions=[np.array(n,float)/np.linalg.norm(n) for n in
        ([1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],[1,1,0],[-1,-1,0])]
    spinors=[np.array([1,0],complex),np.array([0,1],complex),
        np.array([1,1],complex)/np.sqrt(2),np.array([1,1j],complex)/np.sqrt(2),
        np.array([1,2+1j],complex)/np.sqrt(6)]
    witnesses=[]
    for n in directions:
        pp=(I+gr(n))/2;pm=I-pp
        ar=pp@gr(r@n)@pm/2;ak=pp@gr((r+1j*s)@n)@pm/2
        lr=(ar-ar.conj().T)/1j;lk=(ak-ak.conj().T)/1j
        for f in spinors:
            a=float(np.vdot(f,lr@f).real);b=float(np.vdot(f,lk@f).real)
            if abs(a)>.01 and abs(b)>.01:
                witnesses.append(dict(direction=n.tolist(),write_symbol=a,total_symbol=b))
                break
    assert witnesses
    # Deliberate cancellation in one orientation: A(r+is,+z)=0,
    # while the opposite orientation retains -2 A(r,+z)^dagger.
    n=np.array([0.,0.,1.]);pp=(I+gr(n))/2;pm=I-pp
    r2=np.zeros((3,3));r2[0,2]=r2[2,0]=1
    target=pp@gr(r2@n)@pm
    # Solve only for the real y coefficient that cancels this block.
    y=np.zeros((3,3));y[1,2]=y[2,1]=1
    unit=pp@gr(y@n)@pm
    coeff=float((np.vdot(1j*unit,-target)/np.vdot(unit,unit)).real)
    k=r2+1j*coeff*y
    plus=pp@gr(k@n)@pm
    minus=pm@gr(k@(-n))@pp
    error=float(np.max(abs(minus+2*target.conj().T)))
    assert np.linalg.norm(plus)<1e-12 and np.linalg.norm(minus)>.1 and error<1e-12
    return dict(original_particle_sterile_indices=idx,joint_witnesses=witnesses,
        cancelled_orientation_norm=float(np.linalg.norm(plus)),
        opposite_orientation_norm=float(np.linalg.norm(minus)),opposite_identity_error=error,
        original_W_B_or_continuum_signal_computed=False)

def run():
    vacuum=(I+Z)/2
    cases=[
        algebra_case('record_changes',Y,(I+.4*X+.2*Z)/2,X,X,vacuum),
        algebra_case('existing_record_read_without_change',(I+Z)/2,(I+.2*Z)/2,Y,X,vacuum),
        algebra_case('both_contributions',.3*Z+.7*Y,(I+.4*X+.2*Z)/2,X+.5*Y,X,vacuum)]
    assert cases[0]['write']!=0 and cases[0]['read']==0
    assert cases[1]['record_change_operator_norm']==0 and abs(cases[1]['total'])>.1
    assert abs(cases[2]['write'])>.1 and abs(cases[2]['read'])>.1
    return dict(round=814,all_checks_passed=True,cases=cases,
        original_principal_symbol=principal_joint_witness(),
        scope='Direct finite operator expansion and original sterile principal-symbol compatibility.',
        original_continuum_numerical_signal=False,arbitrary_fixed_process_nonzero_proven=False,
        autonomous_writing_or_stable_storage_proven=False,finite_coupling_convergence_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

