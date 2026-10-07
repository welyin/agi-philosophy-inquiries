"""925: finite autonomous record transfer and the same internal resource ledger.
The conserved ledger and additive bare energies are explicit model inputs.
No continuum limit, Gauss physical-state restriction, or physical gauge group is assumed.
"""
from pathlib import Path
import argparse, hashlib, json
from fractions import Fraction
import numpy as np
HERE=Path(__file__).resolve().parent

def op(a):return float(np.linalg.norm(a,2))
def comm(a,b):return a@b-b@a
def unitary(h,t):
    e,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*e))@v.conj().T

def rational_rank(rows):
    a=[[Fraction(x) for x in row] for row in rows];rank=0
    for col in range(len(a[0])):
        pivot=next((j for j in range(rank,len(a)) if a[j][col]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank]
        lead=a[rank][col];a[rank]=[x/lead for x in a[rank]]
        for j in range(rank+1,len(a)):
            f=a[j][col]
            if f:a[j]=[x-f*y for x,y in zip(a[j],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank

def classification(nlevels,r):
    rows=[]
    for q in range(r):
        for n in range(nlevels-q):
            row=[0]*(nlevels+r);row[n+q]+=1;row[n]-=1;row[nlevels+q]=1;rows.append(row)
    rank=rational_rank(rows)
    assert rank==nlevels+r-2
    intercept=[1]*nlevels+[0]*r
    slope=list(range(nlevels))+[-q for q in range(r)]
    assert all(sum(x*y for x,y in zip(row,basis))==0 for row in rows for basis in (intercept,slope))
    return dict(reference_levels=nlevels,payload_dimension=r,linear_constraints=len(rows),exact_rational_rank=rank,
                free_parameters=2,affine_spectrum_and_additive_endpoint_gap_span_nullspace=True)

def model(nlevels=10,r=3,g=.4,a=.7,kappa=0):
    d=2*r*nlevels
    idx=lambda loc,q,n:(loc*r+q)*nlevels+n
    ea=np.array([1.1,1.8,2.6])[:r];eb=ea-a*np.arange(r)
    f=2+a*np.arange(nlevels)+kappa*np.arange(nlevels)**2
    nref=np.diag(np.tile(np.arange(nlevels),2*r)).astype(complex)
    qa=np.zeros((d,d),complex);qb=qa.copy();href=qa.copy();hm=qa.copy();j=qa.copy();pairs=[]
    for q in range(r):
        for n in range(nlevels):
            ia=idx(0,q,n);ib=idx(1,q,n)
            qa[ia,ia]=q;qb[ib,ib]=q
            href[ia,ia]=href[ib,ib]=f[n]
            hm[ia,ia]=ea[q];hm[ib,ib]=eb[q]
            if n+q<nlevels:
                dst=idx(1,q,n+q);j[dst,ia]=1
                pairs.append((q,n,ia,dst,float(eb[q]+f[n+q]-ea[q]-f[n])))
    h0=hm+href;hint=g*(j+j.conj().T);h=h0+hint
    return locals()

def run():
    classified=[classification(n,r) for n,r in ((5,2),(6,3),(10,3),(8,4))]
    z=model();h=z['h'];idx=z['idx'];d=z['d'];r=z['r'];nlevels=z['nlevels'];g=z['g'];a=z['a']
    t=np.pi/(2*g);u=unitary(h,t)
    conservation={
        'left_ledger':op(comm(h,z['qa']+z['nref'])),
        'right_ledger':op(comm(h,z['qb']-z['nref'])),
        'bare_energy':op(comm(h,z['h0']))}
    assert max(conservation.values())<1e-13
    current_error=op(1j*comm(h,z['nref'])-1j*comm(h,z['qb']))
    assert current_error<1e-14
    initial_cols=[];target_cols=[]
    for q,n,src,dst,delta in z['pairs']:
        e=np.zeros(d,complex);e[src]=1;initial_cols.append(e)
        f=np.zeros(d,complex);f[dst]=-1j*np.exp(-1j*t*(z['ea'][q]+z['f'][n]));target_cols.append(f)
    initial=np.column_stack(initial_cols);target=np.column_stack(target_cols)
    isometry_error=op(u@initial-target)
    assert isometry_error<3e-14 and np.min(np.linalg.eigvalsh(h))>0
    # Unknown qutrit with an old reference. The resource starts at n=3, within all task margins.
    psi=np.zeros((d,r),complex)
    for q in range(r):psi[idx(0,q,3),q]=1/np.sqrt(r)
    out=u@psi
    # Trace location and resource; undo only the known endpoint energy phase on the payload.
    amp=out.reshape(2,r,nlevels,r)
    amp=amp*np.exp(1j*t*z['ea'])[None,:,None,None]
    receiver_old=np.einsum('pqnr,psnt->qrst',amp,amp.conj()).reshape(r*r,r*r)
    bell=np.eye(r).reshape(-1)/np.sqrt(r)
    fidelity=float(np.vdot(bell,receiver_old@bell).real)
    td=float(np.linalg.norm(receiver_old-np.outer(bell,bell.conj()),ord='nuc')/2)
    # Relational decoding on the output code restores all amplitudes, not just a chosen input.
    recovered=target.conj().T@u@initial
    recovery_error=op(recovered-np.eye(len(z['pairs'])))
    assert abs(fidelity-1/r)<1e-13 and abs(td-(1-1/r))<1e-13 and recovery_error<3e-14
    expectation=lambda state,obs:float(np.einsum('ij,ik,kj->',state.conj(),obs,state).real)
    flow=dict(reference_energy_change=expectation(out,z['href'])-expectation(psi,z['href']),
              endpoint_energy_change=expectation(out,z['hm'])-expectation(psi,z['hm']),
              interaction_energy_change=expectation(out,z['hint'])-expectation(psi,z['hint']))
    assert abs(flow['reference_energy_change']-a)<2e-13 and abs(flow['endpoint_energy_change']+a)<2e-13
    assert abs(sum(flow.values()))<2e-13
    # Same architecture and same transfer term: a nonlinear resource spectrum breaks universal resonance.
    bad=model(kappa=.04);ub=unitary(bad['h'],t)
    diagnostic=[]
    for q,n,src,dst,delta in bad['pairs']:
        omega=np.sqrt(g*g+(delta/2)**2);cap=g*g/(g*g+(delta/2)**2)
        analytic=cap*np.sin(omega*t)**2;actual=float(abs(ub[dst,src])**2)
        assert abs(actual-analytic)<2e-14
        diagnostic.append(dict(q=q,n=n,detuning=delta,maximum_transfer_probability=cap,
                               probability_at_original_time=actual,formula_error=abs(actual-analytic)))
    witness=next(x for x in diagnostic if x['q']==2 and x['n']==5)
    assert abs(witness['maximum_transfer_probability']-25/61)<1e-13
    # Fixed comparison parameter with no resource shift: charges fail if q is moved.
    bare_j=np.zeros((d,d),complex)
    for q in range(r):
        for n in range(nlevels):bare_j[idx(1,q,n),idx(0,q,n)]=1
    bare_h=bad['h0']+g*(bare_j+bare_j.conj().T)
    bare_ledger_error=op(comm(bare_h,z['qa']+z['nref']))
    assert abs(bare_ledger_error-.8)<1e-14
    return dict(round=925,date='2026-10-06',all_scientific_checks_passed=True,
        exact_classification=classified,model_dimension=d,positive_total_H_minimum=float(np.min(np.linalg.eigvalsh(h))),
        transfer_time=t,conservation_residuals=conservation,same_transfer_current_identity_error=current_error,
        all_input_isometry_error=isometry_error,relational_recovery_error=recovery_error,
        bare_receiver_entanglement_fidelity=fidelity,bare_receiver_trace_distance_from_ideal=td,
        resource_energy_flow=flow,nonlinear_spectrum_counterexample=witness,
        nonlinear_probability_formula_max_error=max(x['formula_error'] for x in diagnostic),
        fixed_dictionary_without_shift_ledger_violation=bare_ledger_error,
        same_autonomous_process_couples_record_and_resource=True,
        additive_bare_energy_and_shift_architecture_are_inputs=True,
        total_energy_conservation_alone_forces_resonance=False,
        autonomous_does_not_mean_apparatus_preparation_derived=True,
        gauge_group_or_angular_reference_derived=False,universal_mass_law_derived=False,
        selected_three_dimensional_spacetime=False,GR_generated=False,
        continuum_or_unlimited_reference_required=False,
        no_battery_size_optimization_performed=True,whole_stage_completed=False,full_goal_completed=False,
        source_hashes={str(Path('925')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();data=run()
    target=HERE/'record_resource_exchange_results.json'
    if a.write:
        assert not target.exists();target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in data.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
