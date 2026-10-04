"""719 entry: reuse633 witness as an accuracy condition, not a new no-go."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
sys.path.insert(0,str(ARCHIVE/'round718_drafts'))
import joint_smooth_mode_contract as shared
from reflected_history_entry import quadratic_matrix
TARGET=HERE/'causal_accuracy_entry_results.json'


def fixture(lam):
    # The same volume-normalised original sterile mode, restricted only for the
    # two-mode measurement algebra. No bosonic or field dynamics are simulated.
    w=np.array([.75,.25])*np.exp(np.array([1.1,-.7])*lam)
    f,_=shared.profile(w,np.array([[1.,0.],[1.,0.]]))
    f=f[[30,62]]
    p=float(abs(f[0])**2);q=1-p
    n=quadratic_matrix(np.outer(f,f.conj()),np.zeros((2,2)))
    R=np.eye(4)-2*n
    nA=np.diag([0.,1.,0.,1.]);nB=np.diag([0.,0.,1.,1.])
    # These are the already established633 continuous admissible witnesses
    # restricted to their two CAR modes.
    plus=np.array([0.,1.,1.,0.])/np.sqrt(2)
    minus=np.array([0.,1.,-1.,0.])/np.sqrt(2)
    rho=[np.outer(x,x.conj()) for x in (plus,minus)]
    out=[.5*(x+R@x@R) for x in rho]
    probs=[float(np.trace(nB@x).real) for x in out]
    contrast=abs(probs[0]-probs[1])
    analytic=2*np.sqrt(p*q)*abs(2*p-1)
    pp=1.8*p*q
    derivative=(4*np.sqrt(p*q)-(2*p-1)**2/np.sqrt(p*q))*pp
    V=np.diag(np.exp(1j*np.diag(.31*nA-.79*nB)))
    Rp=V@R@V.conj().T;Bp=V@nB@V.conj().T
    transported=[]
    for x in rho:
        xp=V@x@V.conj().T
        transported.append(float(np.trace(Bp@(.5*(xp+Rp@xp@Rp))).real))
    error=float(max(abs(np.array(probs)-transported)))
    assert abs(contrast-analytic)<1e-13 and error<1e-13
    # Identity is causal; it saturates just this binary-witness lower bound,
    # not the distance between the entire two channels.
    identity_errors=[abs(z-.5) for z in probs]
    assert abs(max(identity_errors)-contrast/2)<1e-13
    return dict(parameter=lam,mode_weight_A=p,original_B_probabilities=probs,
                signalling_witness=contrast,analytic_witness=float(analytic),
                minimum_uniform_probability_or_trace_distance_error=contrast/2,
                witness_derivative=float(derivative),
                representation_transport_error=error,
                identity_binary_witness_errors=identity_errors)


def run():
    rows=[fixture(x) for x in (-.1,0.,.15)]
    eps=1e-5
    fd=(fixture(eps)['signalling_witness']-fixture(-eps)['signalling_witness'])/(2*eps)
    assert abs(fd-rows[1]['witness_derivative'])<1e-9
    names=('research_note_633.md','research_note_667.md','research_note_718.md',
           'joint_smooth_mode_contract.py','round718_drafts/reflected_history_entry.py')
    return dict(entry_round=719,new_formal_round=False,rows=rows,
                central_source_derivative=fd,source_derivative_error=abs(fd-rows[1]['witness_derivative']),
                conditional_bound_for_candidate_signal_s='epsilon >= max(0, witness-s)/2',
                exact_old_633_witness_reused=True,no_new_lightcone_or_optimal_device_claim=True,
                not_a_new_instantaneous_signalling_counterexample=True,
                continuous_witness_regularity_inherited_from_633=True,
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
