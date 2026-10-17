"""1270 finite future-boundary output: one author group.
Default read-only; --write creates result exclusively. No old ODE rerun.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir())
ARCH=HERE.parent
TOL=2e-12

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def unitary(h,t=1.):
    a,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*a))@u.conj().T

def td(a,b): return float(np.sum(np.abs(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2)))/2)
def dens(v): return np.outer(v,v.conj())
def kron(*items):
    out=np.array([[1]],complex)
    for v in items: out=np.kron(out,v)
    return out

def partial(rho,dims,keep):
    a=rho.reshape(tuple(dims)+tuple(dims))
    live=list(dims)
    for i in sorted(set(range(len(dims)))-set(keep),reverse=True):
        a=np.trace(a,axis1=i,axis2=i+len(live))
        live.pop(i)
    size=int(np.prod(live))
    return a.reshape(size,size)

def close(name,a,b,checks):
    err=float(np.linalg.norm(a-b))
    assert err<TOL,(name,err)
    checks[name]=True

def calculate():
    checks={}
    i=np.eye(2);x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex);z=np.diag([1.,-1.])
    p0=(i+z)/2;p1=(i-z)/2;plus=(i+y)/2;minus=(i-y)/2
    # Internal order C,F,E. Both finite generators, including noise source.
    go=np.pi/2*kron(minus,i-x,i)
    gr=np.pi/2*kron(i,i-x,p1)
    uo=unitary(go);ur=unitary(gr)
    ideal=kron(plus,i,i)+kron(minus,x,i)
    flip=kron(i,i,p0)+kron(i,x,p1)
    close('exact_sigma_y_flag_gate',uo,ideal,checks)
    close('retained_noise_control_gate',ur,flip,checks)
    close('two_internal_generators_commute',go@gr,gr@go,checks)
    assert abs(np.linalg.norm(go,2)-np.pi)<TOL
    assert abs(np.linalg.norm(gr,2)-np.pi)<TOL
    # Full S,C,F,E,R state, with arbitrary reference correlations.
    rng=np.random.default_rng(1270)
    a=rng.normal(size=(32,32))+1j*rng.normal(size=(32,32))
    rho=a@a.conj().T;rho/=np.trace(rho)
    full=kron(i,ur@uo,i)
    actual=full@rho@full.conj().T
    dims=[2,2,2,2,2]
    close('old_material_marginal_unchanged',partial(actual,dims,[0]),partial(rho,dims,[0]),checks)
    close('passive_reference_unchanged',partial(actual,dims,[4]),partial(rho,dims,[4]),checks)
    projs=[kron(i,i,p,i,i) for p in [p0,p1]]
    after=sum(p@actual@p for p in projs)
    close('full_instrument_complete',sum(p@p for p in projs),np.eye(32),checks)
    close('nonselective_source_state_unchanged',partial(after,dims,[0]),partial(rho,dims,[0]),checks)
    # Source cancellation with a noncommuting geometry coefficient and Jbase.
    hb=np.array([[.8,.2j],[-.2j,1.3]],complex)
    ac=np.array([[1.1,.4],[.4,.7]],complex)
    assert np.linalg.norm(ac@hb-hb@ac)>.01
    jbase=kron(-hb,np.eye(8));aa=kron(ac,np.eye(8))
    vv=kron(i,.003*(go+gr));jout=jbase-vv
    tbase=(aa@jbase+jbase@aa)/2
    tout=(aa@jout+jout@aa)/2+(aa@vv+vv@aa)/2
    close('same_linear_source_completion_cancels',tout,tbase,checks)
    # A finite noisy diagnostic checks cq and normalized conditional bounds.
    # eta is an arithmetic test perturbation, not the physical d_R.
    psi_s=np.array([np.sqrt(.7),1j*np.sqrt(.3)])
    psi_c=unitary(z,.17)@np.array([1,1],complex)/np.sqrt(2)
    b=.0001;r=.0002;xi=.0003
    cmp=kron(dens(psi_s),dens(psi_c),np.diag([1-b,b]),np.diag([1-r,r]),p0)
    ent=unitary(kron(x,z,i,i,i),.001)
    real=ent@cmp@ent.conj().T
    din=td(real,cmp)
    noisy=kron(i,ur@unitary(go,1+xi),i)
    real=noisy@real@noisy.conj().T;cmp=noisy@cmp@noisy.conj().T
    qr=[p@real@p for p in projs];qc=[p@cmp@p for p in projs]
    cq_distance=sum(td(a,b) for a,b in zip(qr,qc))
    assert cq_distance<=din+TOL
    for index,(branch_real,branch_cmp) in enumerate(zip(qr,qc)):
        pa=float(np.trace(branch_real).real);pb=float(np.trace(branch_cmp).real)
        assert pa>.25 and pb>.25
        condition_distance=td(branch_real/pa,branch_cmp/pb)
        assert condition_distance<=2*din/(pb-din)+TOL
        assert condition_distance<=8*din+TOL
        old_s=partial(branch_cmp/pb,dims,[0])
        close('comparison_condition_material_'+str(index),old_s,dens(psi_s),checks)
    checks['full_cq_and_conditional_reference_bound_diagnostic']=True
    # Probability error is different from full purified-noise state distance.
    for angle in [.125,.17,.1875]:
        state=unitary(z,angle)@np.array([1,1],complex)/np.sqrt(2)
        inp=kron(dens(state),np.diag([1-b,b]),np.diag([1-r,r]))
        out=ur@unitary(go,1+xi)@inp@unitary(go,1+xi).conj().T@ur.conj().T
        p=float(np.trace(kron(i,p0,i)@out).real)
        expected=float(np.real(state.conj()@plus@state))
        assert abs(p-expected)<=b+np.pi*abs(xi)+r+TOL
    checks['explicit_output_error_model']=True
    # Authoritative exact budgets.
    delta=F(1,10**9);pi_upper=F(22,7)
    beta=2*delta+pi_upper*delta
    assert beta==F(36,7*10**9)<F(6,10**9)<F(1,10**7)
    dR=F(1,2**100)
    assert F(5,16)-beta-dR>F(1,4)
    travel=F(21,400);nominal=F(33,16);jitter=F(1,256)
    duration=F(1,64);deadline=F(17,8)
    early=nominal-jitter-F(2)-travel
    late=deadline-(nominal+jitter+2*duration)
    assert early==F(39,6400)>0 and late==F(7,256)>0
    pulse_peak=F(35,16)/duration
    assert pulse_peak==140
    gate_peak=pulse_peak*pi_upper*(1+delta)
    assert gate_peak==440*(1+delta)
    # Same finite integers as 1268, with a wider common conditional-source bound.
    M=2**3059;EH=3093;E0=14008;EF=95928+100*M
    N=1+max((EF+44)//3,(EH+EF+47)//6)
    conditional_exp=E0+EH+15*M+27
    source_peak_exp=E0+60
    assert 2*N-conditional_exp>82
    assert 4*N-source_peak_exp>82
    assert 4*N-E0-12>100
    old=json.loads((ARCH/'1268/readout_results.json').read_text(encoding='utf-8'))
    assert F(old['finite_prediction']['new_actual_matching_input_per_task'])>beta
    assert old['task']['coupling']=='epsilon*alpha*g*F(Q)*sigma_z'
    deps={
      '1270/admission.md':HERE/'admission.md',
      '1270/output_proof.md':HERE/'output_proof.md',
      '1270/output_check.py':Path(__file__),
      '1259/proof.md':ARCH/'1259/proof.md',
      '1268/readout_proof.md':ARCH/'1268/readout_proof.md',
      '1268/readout_results.json':ARCH/'1268/readout_results.json',
      '1269/material_proof.md':ARCH/'1269/material_proof.md',
      '1269/pointer_proof.md':ARCH/'1269/pointer_proof.md',
      '1269/acceptance.json':ARCH/'1269/acceptance.json',
      '_admission/causal_block_after1265/proof.md':ARCH/'_admission/causal_block_after1265/proof.md',
      '_admission/positive_model_alternative_admission_after1252.md':ARCH/'_admission/positive_model_alternative_admission_after1252.md',
    }
    return {
      'round':1270,'part':'finite_future_block_output_and_conditional_source',
      'status':'PASS_AUTHOR_WITH_SCOPE','author_groups':1,'tolerance':TOL,
      'matrix_checks':checks,
      'rational_authoritative_bounds':{
        'output_probability_error_strict_upper':str(beta),
        'inherited_matching_allowance':str(F(1,10**7)),
        'minimum_comparison_outcome_probability':str(F(5,16)-beta),
        'minimum_actual_outcome_probability_strict_lower':'1/4',
        'conditional_state_multiplier_of_dR':'8',
        'conditional_bounded_source_multiplier':'16*norm(B)*dR',
        'pulse_duration':str(duration),'pulse_peak':str(pulse_peak),
        'clock_offset_max':str(jitter),'processing_start':str(nominal),
        'arrival_before_processing_margin':str(early),
        'output_deadline':str(deadline),'final_time_margin':str(late),
        'normalized_gate_mass_peak_coefficient':str(gate_peak),
        'physical_gate_total_energy_peak':'440*(1+10^-9)*hbar/s',
        'same_source_completed_cancellation':'exact with common linear completion',
        'common_source_exponent_margin_over_82_bits':True,
        'common_source_error_strict_upper':'2^-80',
      },
      'finite_resources':{
        'pointer_code_dimension':2,'fresh_flag_dimension':2,'noise_bit_dimension':2,
        'same_positive_baseline_extension':True,
        'fresh_or_retained_inventory_per_repeat':True,
        'no_free_reset':True,'noise_and_purification_retained':True,
      },
      'scope':[
        'New finite effective material primitive on the same whole-block future boundary.',
        'No finer in-window spatial ports or microscopic cell implementation inferred.',
        'Gate and bit-flip sources use the same action and same clock/source completion.',
        'Probability beta is not a full purified-environment distance or source norm error.',
        'The same actual noisy map transports d_R, so no beta times large source norm occurs.',
        'Conditional uniform bounds apply to the declared nonrare Ramsey preparations.',
        'Arbitrary encoded pointer-reference inputs retain full cq bounds; rare conditioning is separate.',
        'Known Gaussian material preparation remains independent of the added finite inventory.',
        'No new axiom, CT adoption, full GR, or experimental implementation claim.',
      ],
      'bindings_sha256':{k:sha(v) for k,v in deps.items()},
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true')
    a=p.parse_args();got=calculate();out=HERE/'output_results.json'
    if a.write:
        with out.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(got,f,ensure_ascii=False,indent=2);f.write('\n')
        print('FIRST_WRITE PASS_AUTHOR_WITH_SCOPE author_groups=1')
    else:
        assert json.loads(out.read_text(encoding='utf-8'))==got,'saved output does not match current bindings'
        print('READ_ONLY_MATCH PASS_AUTHOR_WITH_SCOPE author_groups=1')

if __name__=='__main__': main()