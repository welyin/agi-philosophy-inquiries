"""742 entry: actual original record exceeds deterministic quadratic dilation.

Reconstruct the partner from the original full pure covariance, then verify
its actual fourth moment independently in the two-mode Fock representation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_retarded_reference_response as response
TARGET=HERE/'gaussian_record_entry_results.json'

def covariance(P):
    n=len(P)//2;I=np.eye(n)
    T=np.block([[I,I],[-1j*I,1j*I]])
    gamma=1j*(np.eye(2*n)-T@P@T.conj().T)
    assert np.max(abs(gamma.imag))<2e-12
    return gamma.real

def pf4(A):
    return A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]

def destroy(m):
    A=np.zeros((4,4),complex)
    for j in range(4):
        bits=[(j>>1)&1,j&1]
        if bits[m]:A[j^(1<<(1-m)),j]=(-1)**sum(bits[:m])
    return A

def run():
    P,_,_=response.flow(0.)
    G=covariance(P);n=len(P)//2
    purity=float(np.max(abs(P@P-P)));assert purity<2e-12
    f=np.zeros(n);f[30]=f[62]=1/np.sqrt(2)
    u=np.r_[f,np.zeros(n)];v=np.r_[np.zeros(n),f]
    outside=np.eye(2*n)-np.outer(u,u)-np.outer(v,v)
    a=float(u@G@v);p=(1+a)/2;kappa=np.sqrt(1-a*a)
    b1=u@G@outside/kappa;b2=v@G@outside/kappa
    W=np.array([u,v,b1,b2]);small=W@G@W.T
    expected=np.array([[0,a,kappa,0],[-a,0,0,kappa],[-kappa,0,0,-a],[0,-kappa,a,0]])
    frame_error=float(np.max(abs(W@W.T-np.eye(4))))
    covariance_error=float(np.max(abs(small-expected)))
    assert max(frame_error,covariance_error)<3e-12
    post=P+response.family.record_change(P)
    bar=W@covariance(post)@W.T
    reflection=np.diag([-1,-1,1,1])
    assert np.max(abs(bar-(small+reflection@small@reflection)/2))<2e-12
    actual_fourth=-.5*(pf4(small)+pf4(reflection@small@reflection))
    wick_fourth=-pf4(bar)
    defect=actual_fourth-wick_fourth
    assert abs(defect-(1-a*a))<4e-12
    # Independent exact four-dimensional CAR calculation for the actual pair.
    c,d=destroy(0),destroy(1)
    majorana=[c+c.conj().T,-1j*(c-c.conj().T),d+d.conj().T,-1j*(d-d.conj().T)]
    H=-.25j*sum(small[i,j]*majorana[i]@majorana[j] for i in range(4) for j in range(4))
    ev,V=np.linalg.eigh(H);psi=V[:,0];rho=np.outer(psi,psi.conj())
    occ=c.conj().T@c;other=d.conj().T@d;I=np.eye(4)
    out=occ@rho@occ+(I-occ)@rho@(I-occ)
    rho_cov=np.array([[float((.5j*np.trace(rho@(majorana[i]@majorana[j]-majorana[j]@majorana[i]))).real)
                      for j in range(4)] for i in range(4)])
    moment=float(np.trace(out@majorana[0]@majorana[1]@majorana[2]@majorana[3]).real)
    n1=float(np.trace(out@occ).real);n2=float(np.trace(out@other).real)
    n12=float(np.trace(out@occ@other).real)
    independent_error=max(float(np.max(abs(rho_cov-small))),abs(moment-actual_fourth),abs(n1-p),abs(n2-(1-p)))
    assert independent_error<4e-12 and abs(n12)<3e-12
    gaussian=np.diag([(1-p)*p,(1-p)**2,p*p,p*(1-p)])
    trace_distance=.5*float(np.sum(abs(np.linalg.eigvalsh(out-gaussian))))
    assert abs(trace_distance-2*p*(1-p))<4e-12
    deps=('research_note_633.md','research_note_634.md','research_note_716.md','research_note_730.md',
          'research_note_735.md','research_note_741.md','joint_retarded_reference_response.py')
    return dict(entry_round=742,latest_completed_round=741,formal_test_count_unchanged=3421,
        original_Nambu_dimension=len(P),original_pure_covariance_error=purity,
        actual_record_occupation_probability=p,majorana_single_mode_parameter=a,
        actual_partner_frame_error=frame_error,actual_pair_covariance_error=covariance_error,
        actual_post_record_fourth_moment=actual_fourth,
        Wick_fourth_moment_from_actual_two_point=wick_fourth,fourth_moment_Wick_defect=defect,
        occupation_connected_defect=n12-n1*n2,
        same_covariance_Gaussian_trace_distance=trace_distance,
        rigorous_lower_bound_to_any_even_Gaussian_state=(1-a*a)/14,
        independent_two_mode_Fock_error=independent_error,
        restricted_implementation='Fixed quadratic parity-even unitary, initially independent even quasifree probe, then discard probe; no external random choice, nonGaussian resource, measurement primitive or feedback.',
        no_claim_that_all_FLO_including_primitive_measurements_is_excluded=True,
        no_claim_that_original_quantized_boson_fermion_interactions_are_quadratic=True,
        nonGaussian_post_state_does_not_invalidate_bilinear_source_in_round741=True,
        dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
