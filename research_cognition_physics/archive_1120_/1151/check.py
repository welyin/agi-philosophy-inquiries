"""1151: bounded one-call transport tester; diagnostics are not a field simulation."""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
def pack(q): return dict(fraction=str(q), value=float(q))
def exact_budget():
    eps=Q(4,39999)+Q(2,10**6)
    ratio=Q(3,160000000)/(1-Q(1,75000))
    gap=Q(1,2)-eps
    robust=gap-2*Q(1,1000)
    assert ratio<Q(188,10**10)
    assert eps<Q(102003,10**9)
    assert gap>Q(499897,10**6) and robust>Q(497897,10**6)
    return dict(name="original_parameter_and_time_window_bounds",passed=True,
      no_boost_full_state_error=pack(eps),clock_time_over_original_T0_upper=pack(ratio),
      ideal_tester_probability_gap_lower=pack(gap),
      gap_with_each_test_error_1e_minus3=pack(robust),
      original_coupling_mass_inertia_and_packet_unchanged=True)

def expi(H,t):
    e,v=np.linalg.eigh(H)
    return (v*np.exp(1j*t*e))@v.conj().T

def car_tester():
    # Two flavors x both spins. Other physical modes are not being simulated.
    dim=16
    def annih(j):
        a=np.zeros((dim,dim),complex)
        for n in range(dim):
            if n>>j&1:
                a[n^(1<<j),n]=(-1)**((n&((1<<j)-1)).bit_count())
        return a
    b=[annih(j) for j in range(4)]
    number=[x.conj().T@x for x in b]
    G=sum((b[s].conj().T@b[2+s]+b[2+s].conj().T@b[s] for s in range(2)),
          np.zeros((dim,dim),complex))
    U=expi(G,-math.pi/4)
    W=np.zeros((dim,4),complex); W[[1,2,4,8],np.arange(4)]=1
    X=np.array([[0,1],[1,0]],complex)
    uc=expi(np.kron(X,np.eye(2)),-math.pi/4)
    omega=.7; H=omega*(number[2]+number[3])
    t=math.pi/omega
    T=expi(H,t)
    C=U.conj().T@T@U@T.conj().T
    compressed=W.conj().T@C@W
    real=(compressed+compressed.conj().T)/2
    expected=1j*np.kron(X,np.eye(2))
    residual=float(np.linalg.norm(compressed-expected,2))
    assert residual<2e-14 and np.linalg.norm(real,2)<2e-14
    assert np.linalg.norm(U@W-W@uc,2)<2e-14
    # Unknown flavor/spin entangled with a three-dimensional passive reference.
    rng=np.random.default_rng(1151)
    q=rng.normal(size=(4,3))+1j*rng.normal(size=(4,3))
    q/=np.linalg.norm(q); psi=W@q
    # CT(-t), unconditional U, CT(t): exactly one operation call.
    path0=U@psi; path1=T@U@T.conj().T@psi
    out=np.stack([path0,path1])/math.sqrt(2)
    plus=(out[0]+out[1])/math.sqrt(2)
    probability=float(np.linalg.norm(plus)**2)
    formula=(1+np.vdot(psi,C@psi).real)/2
    assert abs(probability-.5)<2e-14 and abs(probability-formula)<2e-14
    full_real=(C+C.conj().T)/2
    assert abs(np.linalg.norm(full_real,2)-1)<2e-14
    return dict(name="full_spin_CAR_compression_and_single_call_tester",passed=True,
      full_Fock_dimension=16,code_dimension=4,reference_dimension=3,
      code_loop_iX_residual=residual,
      code_Hermitian_part_norm=float(np.linalg.norm(real,2)),
      one_call_probability=probability,probability_formula_residual=float(abs(probability-formula)),
      full_Fock_Hermitian_part_norm=float(np.linalg.norm(full_real,2)),
      statement_restricted_to_declared_one_particle_code=True,
      controlled_candidate_or_inverse_required=False,physical_sea_simulation=False)

def continuous_boost():
    # A normal 1+1 dimensional two-mass representation, used only as an
    # independent change-of-variables and one-call probability diagnostic.
    gamma=1.25; sh=.75
    masses=np.array([50.,50.3]); omega=masses[1]-masses[0]
    sigma=.3; t=math.pi/omega; s=t/sh
    X=np.array([[0,1],[1,0]],complex)
    U=(np.eye(2)-1j*X)/math.sqrt(2)
    rng=np.random.default_rng(1152)
    q=rng.normal(size=(2,4))+1j*rng.normal(size=(2,4)); q/=np.linalg.norm(q)
    norm=(2*math.pi*sigma*sigma)**(-.25)
    f=lambda p: norm*np.exp(-p*p/(4*sigma*sigma))
    trap=getattr(np,"trapezoid",np.trapz)
    rows=[]
    for count in (15001,30001):
        bound=9*sigma
        lower=gamma*(-bound)+sh*math.sqrt(masses[0]**2+bound**2)
        upper=gamma*bound+sh*math.sqrt(masses[1]**2+bound**2)
        k=np.linspace(lower,upper,count)
        path0=np.zeros((count,2,4),complex)
        path1=np.zeros_like(path0); boosted=np.zeros_like(path0)
        for a,ma in enumerate(masses):
            Ea=np.sqrt(ma*ma+k*k)
            p=gamma*k-sh*Ea
            Epa=np.sqrt(ma*ma+p*p)
            jac=np.sqrt(Epa/Ea)
            boosted[:,a,:]=jac[:,None]*f(p)[:,None]*q[a]
            path0[:,a,:]=jac[:,None]*f(p)[:,None]*(U@q)[a]
            for b,mb in enumerate(masses):
                kb=gamma*p+sh*np.sqrt(mb*mb+p*p)
                path1[:,a,:]+=jac[:,None]*f(p)[:,None]*U[a,b]*q[b]*np.exp(1j*s*(k-kb))[:,None]
        overlap=trap(np.sum(path0.conj()*path1,axis=(1,2)),k)
        prob_lab=float((1+overlap.real)/2)
        norms=[float(trap(np.sum(abs(v)**2,axis=(1,2)),k)) for v in (boosted,path0,path1)]
        # Independently integrate the rest-frame clock commutator.
        p=np.linspace(-bound,bound,count)
        E=np.sqrt(p[:,None]**2+masses[None,:]**2)
        clock0=(U@q)[None,:,:]*np.exp(-1j*t*E)[:,:,None]*f(p)[:,None,None]
        before=q[None,:,:]*np.exp(-1j*t*E)[:,:,None]
        clock1=np.einsum('ab,nbj->naj',U,before)*f(p)[:,None,None]
        rest=trap(np.sum(clock0.conj()*clock1,axis=(1,2)),p)
        prob_rest=float((1+rest.real)/2)
        # Bare U commutes pointwise with every K phase.
        bare0=np.einsum('ab,nbj->naj',U,boosted)
        bare1=np.exp(1j*s*k)[:,None,None]*np.einsum(
            'ab,nbj->naj',U,np.exp(-1j*s*k)[:,None,None]*boosted)
        bare_prob=float(trap(np.sum(abs((bare0+bare1)/2)**2,axis=(1,2)),k))
        assert max(abs(x-1) for x in norms)<2e-12
        assert abs(prob_lab-prob_rest)<2e-11
        assert abs(bare_prob-1)<2e-12 and bare_prob-prob_lab>.499
        rows.append(dict(points=count,norm_residual=max(abs(x-1) for x in norms),
          transported_probability=prob_lab,rest_clock_probability=prob_rest,
          Lorentz_pullback_probability_residual=abs(prob_lab-prob_rest),
          bare_probability=bare_prob,probability_gap=bare_prob-prob_lab))
    assert abs(rows[0]['transported_probability']-rows[1]['transported_probability'])<1e-11
    return dict(name="normal_continuous_boost_and_lab_tester_diagnostic",passed=True,
      masses=masses.tolist(),momentum_standard_deviation=sigma,
      finite_translation=s,comparison_time=t,rows=rows,
      physical_1148_parameters=False,full_interacting_field_or_domain_simulation=False)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    sources=json.loads((HERE/'source_hashes.json').read_text(encoding='utf-8'))
    for name,digest in sources.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    result=dict(round=1151,all_passed=True,numerical_check_groups=3,
      groups=[exact_budget(),car_tester(),continuous_boost()],frozen_source_count=len(sources),
      new_adopted_axioms=0,scientific_calibration_increment=0,
      new_controlled_translation_menu_is_input=True,
      prior_narrow_endpoint_contract_refuted=False,goal_complete=False)
    dest=HERE/'results.json'
    if args.write_results:
        if dest.exists():
            assert json.loads(dest.read_text(encoding='utf-8'))==result,'Refuse different result overwrite'
        else:
            dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert json.loads(dest.read_text(encoding='utf-8'))==result,'Saved result mismatch'
    print(json.dumps(dict(status='PASS',groups=3,
      probability_gap_lower=result['groups'][0]['ideal_tester_probability_gap_lower']['value'],
      diagnostic_pullback_residual=result['groups'][2]['rows'][-1]['Lorentz_pullback_probability_residual'])))
if __name__=='__main__': main()
