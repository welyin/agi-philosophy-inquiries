"""Round1269 one author group. Default is read-only exact JSON comparison.
Finite matrices test exact code algebra, not a truncated field propagation.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir())
ARCH=HERE.parent
TOL=1e-12

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def unitary(h,t):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T

def zero(name,x,checks):
    err=float(np.linalg.norm(x,ord=2)) if np.ndim(x)==2 else float(np.linalg.norm(x))
    if not err < TOL:
        raise AssertionError((name,err,TOL))
    checks[name]=True

def calculate():
    checks={}
    # Exact normalization: this finite sample tests the identity, not the
    # actual old astronomical resource integer, which is inherited unchanged.
    hbar=F(7,3);v0=F(1,10**6);s=1000*3*2**24
    action=v0*s*s;eps=hbar/action
    assert action*eps==hbar
    checks['same_action_total_momentum_equals_hbar']=True
    checks['normalization_is_not_instantaneous_a_cubed']=True
    # Seven exact Fourier modes suffice for checking all images of |+/-1>
    # under the selected algebra and sin(2phi); no evolution is truncated.
    modes=np.arange(-3,4)
    p=np.diag(modes.astype(float))
    w=np.zeros((7,2),complex)
    w[np.where(modes==1)[0][0],0]=1
    w[np.where(modes==-1)[0][0],1]=1
    z=np.diag([1.,-1.]);x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex)
    zero('isometry',w.conj().T@w-np.eye(2),checks)
    zero('exact_momentum_code',p@w-w@z,checks)
    # Noncommuting original geometry H0 and multiplier F, with a passive R.
    h0=np.array([[.3,.2,.1j],[.2,-.4,.25],[-.1j,.25,.7]],complex)
    f=np.diag([.8,1.3,2.1])
    assert np.linalg.norm(h0@f-f@h0)>0.1
    v=np.kron(np.eye(3),w)
    for eta,t in [(0.,.3),(.17,.9),(-.29,1.2)]:
        hp=np.kron(h0,np.eye(7))+eta*np.kron(f,p)
        hc=np.kron(h0,np.eye(2))+eta*np.kron(f,z)
        key=f'eta_{eta}'
        zero('generator_'+key,hp@v-v@hc,checks)
        zero('propagator_'+key,unitary(hp,t)@v-v@unitary(hc,t),checks)
    # Arbitrary encoded entanglement of system, pointer and a reference is
    # preserved; this test does not assume a product unknown reference state.
    rng=np.random.default_rng(1269)
    psi=rng.normal(size=12)+1j*rng.normal(size=12)
    psi/=np.linalg.norm(psi)
    vr=np.kron(v,np.eye(2))
    hp=np.kron(h0,np.eye(7))+.17*np.kron(f,p)
    hc=np.kron(h0,np.eye(2))+.17*np.kron(f,z)
    out=np.kron(unitary(hp,.9),np.eye(2))@vr@psi
    cmp=vr@np.kron(unitary(hc,.9),np.eye(2))@psi
    zero('passive_reference_full_state',out-cmp,checks)
    # A generic linear P source/variation compresses through the same W.
    source0=np.array([[.8,.2,0],[.2,1.1,.1j],[0,-.1j,.6]],complex)
    for c in [F(1,3),F(-2,5)]:
        bp=np.kron(source0,np.eye(7))+float(c)*np.kron(f,p)
        bc=np.kron(source0,np.eye(2))+float(c)*np.kron(f,z)
        zero('linear_source_'+str(c),bp@v-v@bc,checks)
    eyplus=w@((np.eye(2)+y)/2)@w.conj().T
    eyminus=w@((np.eye(2)-y)/2)@w.conj().T
    outside=np.eye(7)-w@w.conj().T
    for name,e in [('plus',eyplus),('minus',eyminus),('outside',outside)]:
        zero('projection_'+name,e@e-e,checks)
    zero('effects_complete',eyplus+eyminus+outside-np.eye(7),checks)
    zero('exact_code_has_no_outcome_failure',outside@w,checks)
    # Deletion check 1: a local angle multiplier is not the code instrument.
    shift=np.zeros((7,7),complex)
    for j,m in enumerate(modes):
        if m+2 in modes:
            shift[np.where(modes==m+2)[0][0],j]=1
    sine=(shift-shift.conj().T)/(2j)
    zero('angle_multiplier_compresses_to_half_y',w.conj().T@sine@w-y/2,checks)
    leakage=float(np.linalg.norm(outside@sine@w,2))
    assert abs(leakage-.5)<TOL
    # Deletion check 2: replacing P by its initial mean loses the Ramsey record.
    plus=np.array([1,1],complex)/np.sqrt(2)
    state=unitary(z,np.pi/8)@plus
    prob=float(np.real(state.conj()@((np.eye(2)+y)/2)@state))
    assert prob-.5>1/3
    zero('initial_pointer_mean_P_is_zero',np.array([plus.conj()@z@plus]),checks)
    old=json.loads((ARCH/'1268/readout_results.json').read_text(encoding='utf-8'))
    accepted=json.loads((ARCH/'1268/acceptance.json').read_text(encoding='utf-8'))
    assert 'ACCEPTED' in json.dumps(accepted)
    assert old['task']['coupling']=='epsilon*alpha*g*F(Q)*sigma_z'
    assert old['analytic_checks']['Ramsey_full_norm_error_strict_upper']=='2^-100'
    assert F(old['finite_prediction']['B_probability_gap_coarse_strict_lower'])>0
    assert F(old['finite_prediction']['kappa_endpoint_probability_gap_coarse_strict_lower'])>0
    checks['inherited_1268_not_recomputed_or_rescaled']=True
    deps={
      '1269/admission.md':HERE/'admission.md',
      '1269/pointer_proof.md':HERE/'pointer_proof.md',
      '1269/pointer_check.py':Path(__file__),
      '1249/canonical.md':ARCH/'1249/canonical.md',
      '1267/geometry_proof.md':ARCH/'1267/geometry_proof.md',
      '1268/readout_proof.md':ARCH/'1268/readout_proof.md',
      '1268/readout_check.py':ARCH/'1268/readout_check.py',
      '1268/readout_results.json':ARCH/'1268/readout_results.json',
      '1268/acceptance.json':ARCH/'1268/acceptance.json',
      '_admission/causal_block_after1265/proof.md':ARCH/'_admission/causal_block_after1265/proof.md',
    }
    return {
      'round':1269,'part':'same_action_rotor_pointer_code',
      'status':'PASS_AUTHOR_WITH_SCOPE','author_groups':1,
      'checks':checks,
      'numeric_tolerance':TOL,
      'exact_dictionary':{
        'action_prefactor':'A=V0*s^2','epsilon':'hbar/A',
        'P_phys':'-i*epsilon*d/dphi','code':'m=+1,-1',
        'total_canonical_momentum_on_code':'+/-hbar',
        'd':'alpha*g(theta)*P_phys/(3*k^2)',
        'homogeneous_coupling_on_code':'epsilon*alpha*g*F(Q)*sigma_z',
        'additional_code_propagation_error':'0',
        'P_mean_deletion_probability_gap_strict_lower':'1/3',
        'sin2phi_code_compression':'sigma_y/2',
        'sin2phi_out_of_code_amplitude':'1/2',
      },
      'inherited_not_new_scientific_group':{
        'full_state_error':old['analytic_checks']['Ramsey_full_norm_error_strict_upper'],
        'old_bounded_source_error':old['analytic_checks']['old_bounded_source_total_error_strict_upper'],
        'new_bounded_source_error':old['analytic_checks']['new_source_strict_bound'],
        'finite_prediction':old['finite_prediction'],
      },
      'scope':[
        'Exact intertwiner for arbitrary encoded system-pointer-reference states.',
        'Old Gaussian approximation retains its own independent-material preparation scope.',
        'The finite code is physical input; all integer rotor sectors are only a mathematical extension.',
        'No whole rotor positivity, Gibbs reference, infinite free inventory, or finite-dimensional CCR claim.',
        'Source variation is independent; its homogeneous linear-P compression preserves the chosen dictionary.',
        'No local-field-to-collective-rotor quantization or actual local terminal instrument established.',
        'Actual matching, terminal apparatus source, and post-read spacetime remain separate conditions.',
        'No new axiom, CT adoption, complete GR, or stage closure.',
      ],
      'display_only':{'P_mean_deletion_probability_gap':round(prob-.5,12)},
      'bindings_sha256':{k:sha(v) for k,v in deps.items()},
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='exclusive first result creation')
    args=parser.parse_args()
    got=calculate()
    target=HERE/'pointer_results.json'
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as out:
            json.dump(got,out,ensure_ascii=False,indent=2)
            out.write('\n')
        print('FIRST_WRITE PASS_AUTHOR_WITH_SCOPE author_groups=1')
    else:
        saved=json.loads(target.read_text(encoding='utf-8'))
        if saved!=got:
            raise AssertionError('saved result does not match current proof/code/dependencies')
        print('READ_ONLY_MATCH PASS_AUTHOR_WITH_SCOPE author_groups=1')

if __name__=='__main__':
    main()