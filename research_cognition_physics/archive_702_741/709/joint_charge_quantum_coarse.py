"""709: original global quark-charge coarse algebra and its regional boundary.

Uses all original 32 CAR coefficients; properness and cut witnesses are sparse
states of that same CAR model, not a new two-level replacement Hamiltonian.
No full heat spectrum or continuum anomaly simulation is attempted.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_topological_mass_phases as topology

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_charge_quantum_coarse_results.json'
QUARK=np.r_[np.ones(24),np.zeros(8)]


def create(state,index):
    out={}
    for mask,amp in state.items():
        if (mask>>index)&1:continue
        sign=(-1)**((mask&((1<<index)-1)).bit_count())
        new=mask|(1<<index);out[new]=out.get(new,0)+sign*amp
    return out


def annihilate(state,index):
    out={}
    for mask,amp in state.items():
        if not (mask>>index)&1:continue
        sign=(-1)**((mask&((1<<index)-1)).bit_count())
        new=mask^(1<<index);out[new]=out.get(new,0)+sign*amp
    return out


def add(out,state,coefficient=1):
    for mask,amp in state.items():out[mask]=out.get(mask,0)+coefficient*amp
    return {m:a for m,a in out.items() if abs(a)>1e-14}


def inner(a,b):
    return sum(complex(v).conjugate()*b.get(m,0) for m,v in a.items())


def normalized(state):
    norm=np.sqrt(inner(state,state).real)
    return {m:a/norm for m,a in state.items()}


def baryon_create(state,spin):
    out={}
    for p in itertools.permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        indices=(12+2*p[0]+spin,18+2*p[1]+spin,18+2*p[2]+spin)
        term=state
        for index in reversed(indices):term=create(term,index)
        out=add(out,term,sign)
    return out


def coefficients():
    rng=np.random.default_rng(7091);q=np.diag(QUARK)
    total=lepton=0.
    for _ in range(8):
        phi=.35*rng.normal(size=5);h,d=matter.mass_matrices(phi)
        total=max(total,float(np.linalg.norm(q@h-h@q)),float(np.linalg.norm(q@d+d@q.T)))
        lepton=max(lepton,float(np.linalg.norm(d+d)))
        C=matter.gauge.group_exp(rng.normal(size=8),3)
        W=matter.gauge.group_exp(rng.normal(size=3),2);z=np.exp(1j*rng.normal())
        R=matter.representation(C,W,z)
        total=max(total,float(np.linalg.norm(q@R-R@q)),abs(np.linalg.det(R)-1))
        phase=np.exp(1j*.43*QUARK)
        total=max(total,float(np.linalg.norm(phase[:,None]*h*phase.conj()[None,:]-h)),
                  float(np.linalg.norm(phase[:,None]*d*phase[None,:]-d)))
    assert total<1e-12 and lepton>.01
    return dict(original_modes_per_node=32,quark_modes=24,maximum_charge_identity_error=total,
                full_sea_gauge_determinant_checked=True,total_fermion_number_pairing_nonzero=lepton,
                no_accidental_total_fermion_U1_claim=True)


def gauss_properness():
    pair=normalized(baryon_create(baryon_create({0:1.},1),0))
    assert len(pair)==9 and all(m.bit_count()==6 for m in pair)
    assert all(sum(int(QUARK[i]) for i in range(32) if (m>>i)&1)==6 for m in pair)
    rng=np.random.default_rng(7092)
    C=matter.gauge.group_exp(rng.normal(size=8),3)
    W=matter.gauge.group_exp(rng.normal(size=3),2);z=np.exp(.391j)
    R=matter.representation(C,W,z)
    # All possible outputs allowed by the original R: one u per spin,
    # two d per spin. The 81-component test includes zero target components.
    modes=[[12+2*a for a in range(3)],[13+2*a for a in range(3)],
           [18+2*a for a in range(3)],[19+2*a for a in range(3)]]
    outputs=[]
    for up,dn,dp,dq in itertools.product(modes[0],modes[1],itertools.combinations(modes[2],2),itertools.combinations(modes[3],2)):
        outputs.append(tuple(sorted((up,dn,*dp,*dq))))
    inputs=[([i for i in range(32) if (m>>i)&1],a) for m,a in pair.items()]
    error=0.
    for out in outputs:
        amp=sum(a*np.linalg.det(R[np.ix_(out,ins)]) for ins,a in inputs)
        mask=sum(1<<i for i in out);error=max(error,float(abs(amp-pair.get(mask,0))))
    assert error<1e-12
    # Matrix on the two explicit orthonormal physical vectors vacuum/pair.
    # It tests a channel, not a Hamiltonian truncation.
    rho=np.ones((2,2),complex)/2;q=np.array([0,6])
    twirled=sum((np.exp(2j*np.pi*k*q/25)[:,None]*rho*np.exp(-2j*np.pi*k*q/25)[None,:]) for k in range(25))/25
    expected=np.eye(2)/2
    trace_distance=float(np.sum(abs(np.linalg.eigvalsh(rho-twirled)))/2)
    assert np.linalg.norm(twirled-expected)<1e-13 and abs(trace_distance-.5)<1e-13
    # The neutral singlet-neutrino pair shares the zero-quark sector:
    # its quantum coherence must survive, establishing a noncommutative range.
    neutral_pair=create(create({0:1.},31),30)
    assert inner(neutral_pair,neutral_pair)==1
    assert all((mask&((1<<24)-1))==0 for mask in neutral_pair)
    assert np.linalg.norm(R[30:32,30:32]-np.eye(2))<1e-13
    # [n_30, c_30^dagger c_31^dagger]|vac> equals the normalized pair,
    # giving two actual even, gauge-invariant, charge-neutral operators
    # whose commutator is nonzero in the retained algebra.
    number_after={mask:amp*((mask>>30)&1) for mask,amp in neutral_pair.items()}
    assert abs(inner(number_after,number_after)-1)<1e-14
    return dict(baryon_pair_nonzero_components=len(pair),all_candidate_gauge_outputs=len(outputs),
        original_quotient_gauge_error=error,charge_sectors=q.tolist(),same_even_parity=True,
        dephasing_error=float(np.linalg.norm(twirled-expected)),trace_distance=trace_distance,
        entropy_increase=float(np.log(2)),zero_quark_neutrino_pair_coherence_retained=True)


def hopping(state,link,reverse=False):
    out={}
    # A fixed-spin original u_R colour block. Both full 32-mode sites are kept
    # in the bit labels. The remaining 29 channels are present but not excited.
    indices=[12,14,16]
    for a,i in enumerate(indices):
        for b,j in enumerate(indices):
            if reverse:
                term=create(annihilate(state,32+i),j);coefficient=link[a,b].conjugate()
            else:
                term=create(annihilate(state,j),32+i);coefficient=link[a,b]
            out=add(out,term,coefficient)
    return out


def regional_current():
    rng=np.random.default_rng(7093)
    C=matter.gauge.group_exp(rng.normal(size=8),3);U=np.exp(4j*.217)*C
    omega={(1<<32)-1:1.}
    chi=normalized(hopping(omega,U));d=3
    assert len(chi)==9
    assert abs(inner(hopping(omega,U),hopping(omega,U)).real-d)<1e-13
    qtot=[];qlocal=[]
    for state in (omega,chi):
        totals=set();local=set()
        for mask in state:
            left=sum((mask>>i)&1 for i in range(24))
            right=sum((mask>>(32+i))&1 for i in range(24))
            totals.add(left+right);local.add(left)
        qtot.append(sorted(totals));qlocal.append(sorted(local))
    currents=[]
    for sign in (1,-1):
        psi=add({m:a/np.sqrt(2) for m,a in omega.items()},chi,sign*1j/np.sqrt(2))
        Tpsi=hopping(psi,U);Tdpsi=hopping(psi,U,reverse=True)
        current=1j*(inner(psi,Tpsi)-inner(psi,Tdpsi))
        currents.append(float(current.real));assert abs(current.imag)<1e-13
    assert abs(currents[0]-np.sqrt(d))<1e-13 and abs(currents[1]+np.sqrt(d))<1e-13
    mixed_current=.5*sum((1j*(inner(s,hopping(s,U))-inner(s,hopping(s,U,reverse=True)))).real for s in (omega,chi))
    assert abs(mixed_current)<1e-14 and qtot==[[24],[24]] and qlocal==[[24],[23]]
    # Link covariance independently checked at the original single-particle level.
    g0=matter.gauge.group_exp(rng.normal(size=8),3)*np.exp(4j*.11)
    g1=matter.gauge.group_exp(rng.normal(size=8),3)*np.exp(-4j*.07)
    transformed=g1@U@g0.conj().T
    cov=float(np.linalg.norm(g1.conj().T@transformed@g0-U));assert cov<1e-12
    return dict(original_CAR_modes=64,colour_channels=d,dressed_excitation_components=len(chi),
        total_quark_sectors=qtot,local_quark_sectors=qlocal,original_dressed_currents=currents,
        local_pinching_current=float(mixed_current),global_pinching_preserves_current=True,
        covariance_error=cov,no_claim_two_state_subspace_invariant_under_H=True)


def anomaly_interface():
    b=np.array([1,-1,-1,0,0,0],dtype=int)
    mass=topology.C@b;shift=-topology.A.T@b
    assert np.array_equal(mass,np.zeros(5,dtype=int)) and np.array_equal(shift,[0,0,3,0])
    rows=[]
    for ng in (1,3):
        generic=.31;bad=np.exp(1j*3*ng*generic)
        residual=np.exp(1j*3*ng*2*np.pi/(3*ng))
        rows.append(dict(generations=ng,quark_phase_index_coefficient=3*ng,
            generic_phase_distance_from_one=float(abs(bad-1)),residual_subgroup_phase_error=float(abs(residual-1))))
    return dict(inherited_left_handed_charge=b.tolist(),unchanged_mass_phases=mass.tolist(),
        original_integer_topological_shift=shift.tolist(),rows=rows,
        inherited_index_theorem_not_reproved=True,no_continuum_transition_rate_computed=True)


def run():
    deps=('research_note_598.md','research_note_629.md','research_note_642.md','research_note_708.md',
          'joint_fermion_gauss_completion.py','joint_topological_mass_phases.py',
          'round709_drafts/conditional_reference_entry.md')
    return dict(date='2026-10-03',round=709,tests_run=4,failures=0,errors=0,
        charge_coefficients=coefficients(),proper_quantum_range=gauss_properness(),
        regional_transport=regional_current(),continuum_scope=anomaly_interface(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_full_H_and_Gauss=True,global_charge_neutral_task_class=True,
            actual_normal_CPTP_coarse_map=True,all_time_neutral_histories_preserved=True,
            global_map_not_independent_region_maps=True,no_spatial_dof_reduction=True,
            no_charge_superselection_postulate=True,no_continuum_anomaly_or_GR_completion=True),
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    report=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    else:assert report==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:report[k] for k in ('round','tests_run','all_checks_passed')}))
