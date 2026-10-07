"""939: operational access to the actual 853 receiver flavour record.

Symmetric measurement principles are inherited literature. New application:
the 853 action, its original source distribution and its late readout contract.
The finite diagnostic is not a continuum or autonomous-apparatus simulation.
"""
from pathlib import Path
from fractions import Fraction
import argparse, hashlib, itertools, json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=HERE.parent/'853/native_receiver_content_results.json'
TARGET=HERE/'receiver_relational_access_results.json'
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)

def norm(a): return float(np.linalg.norm(a,2))
def real(a):
    v=complex(a)
    assert abs(v.imag)<2e-12
    return float(v.real)
def hermitian_function(a,f):
    e,v=np.linalg.eigh(a)
    return (v*f(e))@v.conj().T
def twirl(a,q):
    # Exact spectral pinching, not a sampled group average.
    e,v=np.linalg.eigh(q)
    b=v.conj().T@a@v
    b[np.abs(e[:,None]-e[None,:])>1e-11]=0
    return v@b@v.conj().T
def annihilator(i,n=4):
    a=np.zeros((2**n,2**n),complex)
    for k in range(2**n):
        if k>>i&1:
            a[k^(1<<i),k]=(-1)**((k&((1<<i)-1)).bit_count())
    return a

def run():
    old=json.loads(OLD.read_text('utf-8-sig'))
    assert old['round']==853 and old['all_checks_passed']
    c=[annihilator(i) for i in range(4)]
    n=[a.conj().T@a for a in c]
    def bilinear(offset,m):
        return sum(m[a,b]*c[offset+a].conj().T@c[offset+b]
                   for a in range(2) for b in range(2))
    xs,ys,zs=(bilinear(0,m) for m in (X,Y,Z))
    xa,ya,za=(bilinear(2,m) for m in (X,Y,Z))
    q=(ys+ya)/2
    bcar=(zs@xa-xs@za)/2
    old_z=2*n[0]-np.eye(16)
    full_car=dict(
        flavour_twirl_original_record_error=norm(twirl(old_z,q)-(n[0]+n[1]-np.eye(16))),
        original_record_charge_commutator=norm(q@old_z-old_z@q),
        relational_record_charge_commutator=norm(q@bcar-bcar@q),
        relational_record_norm=norm(bcar))
    assert full_car['flavour_twirl_original_record_error']<2e-14
    assert full_car['original_record_charge_commutator']>.9
    assert full_car['relational_record_charge_commutator']<2e-14
    assert abs(full_car['relational_record_norm']-1)<2e-14

    # The selected one-particle sector of each pair of CAR modes.
    q2=(np.kron(Y,I)+np.kron(I,Y))/2
    b=(np.kron(Z,X)-np.kron(X,Z))/2
    rho_x=(I+X)/2
    rho_pair=twirl(np.kron(rho_x,rho_x),q2)
    rho_independent=np.eye(4)/4
    assert norm(twirl(rho_x,Y/2)-I/2)<2e-14
    assert norm(rho_pair@q2-q2@rho_pair)<2e-14
    assert norm(b@q2-q2@b)<2e-14
    # A jointly invariant state still carries relative coherence.
    rt=rho_pair.reshape(2,2,2,2)
    assert norm(np.trace(rt,axis1=1,axis2=3)-I/2)<2e-14
    assert norm(np.trace(rt,axis1=0,axis2=2)-I/2)<2e-14
    assert abs(real(np.trace(rho_pair@np.kron(X,X)))-.5)<2e-14
    assert abs(real(np.trace(rho_pair@np.kron(Z,Z)))-.5)<2e-14

    eplus=(np.eye(4)+b)/2
    eminus=(np.eye(4)-b)/2
    a=hermitian_function(eplus,lambda e:np.sqrt(np.maximum(e,0)))
    d=hermitian_function(eminus,lambda e:np.sqrt(np.maximum(e,0)))
    # Pointer is neutral. This is a finite allowed measurement extension,
    # not a claim that the old continuum action generates this unitary.
    uread=np.block([[a,-d],[d,a]])
    qread=np.kron(I,q2)
    dilation=dict(
        unitary_error=norm(uread.conj().T@uread-np.eye(8)),
        total_flavour_conservation_error=norm(uread@qread-qread@uread),
        first_effect_error=norm(a.conj().T@a-eplus),
        generator_norm_bound=float(np.pi/2))
    assert max(v for k,v in dilation.items() if k.endswith('error'))<3e-14
    scan=[]
    for theta in (-.7,-.2,0.,.2,.7):
        u=np.cos(theta/2)*I+1j*np.sin(theta/2)*Y
        us=np.kron(u,I)
        rho=us@rho_pair@us.conj().T
        ideal=real(np.trace(u@rho_x@u.conj().T@Z))
        rel=real(np.trace(rho@b))
        lost=real(np.trace(us@rho_independent@us.conj().T@b))
        assert abs(rel-ideal/2)<2e-14 and abs(lost)<2e-14
        # End-point probabilities from the actual neutral-pointer dilation.
        rho8=np.kron(np.diag([1.,0.]),rho)
        after=uread@rho8@uread.conj().T
        assert abs(real(np.trace(after[:4,:4]))-(1+rel)/2)<2e-14
        scan.append(dict(theta=theta,ideal_record=ideal,
                         relational_record=rel,independent_twirl_record=lost))

    # Reuse the exact 853 source distribution and physical diagnostic inputs.
    fs=lambda key:[float(Fraction(x)) for x in old[key]]
    lam,eta,amp,nu=(fs(k) for k in ('lambdas','etas','source_amplitudes','probe_variances'))
    signs=list(itertools.product((-1,1),repeat=3))
    inputs={key:[float(Fraction(x)) for x in probs]
            for key,probs in old['original_code_joint_probabilities'].items()}
    hbar=1/3
    # Gaussian integration computed through independent 4x4 density matrices.
    nodes,weights=np.polynomial.hermite.hermgauss(32)
    weights=weights/np.sqrt(np.pi)
    conditional={}
    max_kernel_error=0.
    for i in range(3):
        for sign in (-1,1):
            rho=np.zeros((4,4),complex)
            for x,w in zip(nodes,weights):
                p=np.sqrt(2*nu[i])*x-lam[i]*np.sqrt(hbar)*amp[i]*sign
                angle=eta[i]*np.sqrt(hbar)*p
                u=np.kron(np.cos(angle/2)*I+1j*np.sin(angle/2)*Y,I)
                rho+=w*u@rho_pair@u.conj().T
            val=real(np.trace(rho@b))
            exact=-.5*sign*np.exp(-eta[i]**2*hbar*nu[i]/2)*np.sin(eta[i]*lam[i]*hbar*amp[i])
            max_kernel_error=max(max_kernel_error,abs(val-exact))
            conditional[i,sign]=[(1+z*val)/2 for z in (-1,1)]
    assert max_kernel_error<3e-14
    distributions={}
    for r,source_weights in inputs.items():
        probabilities=[]
        for z in signs:
            value=sum(pa*np.prod([conditional[i,aa[i]][(z[i]+1)//2]
                                  for i in range(3)])
                      for aa,pa in zip(signs,source_weights))
            probabilities.append(float(value))
        assert min(probabilities)>0 and abs(sum(probabilities)-1)<2e-14
        distributions[r]=probabilities
    delta=np.array(distributions['3/5'])-np.array(distributions['0'])
    parity=float(sum(np.prod(z)*v for z,v in zip(signs,delta)))
    old_delta=old['preparation_rows'][0]['parity_content_difference']
    assert abs(parity-old_delta/8)<2e-14 and abs(parity)>1e-7
    leading=Fraction(old['exact_hbar_cubed_leading_coefficient'])/8
    covariance_loss=norm(rho_pair-rho_independent)
    sourcefiles=[OLD,HERE.parent/'853/native_receiver_content.py',
                 HERE.parent/'research_note_853.md',HERE.parent/'research_note_934.md']
    return dict(round=939,date='2026-10-07',all_scientific_checks_passed=True,
        full_CAR=full_car,neutral_pointer_dilation=dilation,
        relational_state_vs_independent_twirl_norm=covariance_loss,
        single_receiver_scan=scan,max_original_gaussian_kernel_error=max_kernel_error,
        original_source_probabilities=old['original_code_joint_probabilities'],
        actual_finite_diagnostic_distributions=distributions,
        original_parity_content_difference=old_delta,
        relational_parity_content_difference=parity,
        maximum_joint_probability_difference=float(max(abs(delta))),
        exact_visibility_ratio='1/8',exact_hbar_cubed_leading_coefficient=str(leading),
        shared_reference_state_is_globally_invariant=True,
        lost_relative_correlations_are_not_replaced_by_local_marginals=True,
        old_853_formal_writing_theorem_preserved=True,
        extra_readout_unitary_is_an_explicit_interface_choice=True,
        full_original_action_generates_readout_proved=False,
        continuum_finite_coupling_or_full_backreaction_computed=False,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcefiles})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        saved=json.loads(TARGET.read_text('utf-8'))
        # Numeric comparisons are explicit; source hashes/claims exact.
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a:compare(a[k],b[k])
            elif isinstance(a,list):
                assert len(a)==len(b)
                for x,y in zip(a,b):compare(x,y)
            elif isinstance(a,float):assert abs(a-b)<5e-13
            else:assert a==b,(a,b)
        compare(out,saved)
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','single_receiver_scan','actual_finite_diagnostic_distributions','original_source_probabilities')},ensure_ascii=False,indent=2))
