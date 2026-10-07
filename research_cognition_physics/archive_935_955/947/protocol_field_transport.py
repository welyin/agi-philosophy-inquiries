"""947: actual two-event protocol in the common portal/Newton process.

Reuses the frozen 929 protocol and 946 field calculation.  No claim that the
terminal logical PVM permissions are already manufactured spatial instruments.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'protocol_field_transport_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def norm(a):return float(np.linalg.norm(a,2))
def read(p):return json.loads(p.read_text('utf-8'))

def run():
    old=load('protocol929',STAGE/'929/joint_protocol_selection.py')
    field=read(STAGE/'946/portal_common_process_results.json')
    assert field['all_scientific_checks_passed']
    theta=np.pi/2;scale=np.pi/30;t=30.;hmax=np.pi
    n=30;couplings=np.array([.5*np.sqrt((j+1)*(n-j)) for j in range(n)])
    hc=scale*(15*np.eye(31)+np.diag(couplings,1)+np.diag(couplings,-1))
    vals,vec=np.linalg.eigh(hc)
    end=vec@(np.exp(-1j*t*vals)*vec[0].conj())
    clock_error=float(np.linalg.norm(end-np.eye(31)[:,-1]))
    assert clock_error<1e-12
    gates,hist=old.history(theta,(2,3))
    sign=(1-2*((old.IND>>1)&1)).astype(float)
    def ZE(s):return sign.reshape((-1,)+(1,)*(s.ndim-1))*s
    def B_at(s,ell):
        # After the first weak event, every remaining gate commutes with Z_E1.
        if ell>=2:return ZE(s)
        relevant=gates[ell:2]
        out=s
        for g in relevant:out=old.act(out,g)
        out=ZE(out)
        for g in reversed(relevant):out=old.act(out,g,inverse=True)
        return out
    def full_B(s):return np.asarray([B_at(s[l],l) for l in range(31)])
    def full_h(s):
        out=scale*15*s.copy()
        for j,g in enumerate(gates):
            out[j+1]+=scale*couplings[j]*old.act(s[j],g)
            out[j]+=scale*couplings[j]*old.act(s[j+1],g,inverse=True)
        return out
    rng=np.random.default_rng(947)
    s=rng.normal(size=(31,old.DIM,2))+1j*rng.normal(size=(31,old.DIM,2))
    s/=np.linalg.norm(s.reshape(-1,2),axis=0)[None,None,:]
    b_square_error=norm((full_B(full_B(s))-s).reshape(-1,2))
    commuting_error=norm((full_B(full_h(s))-full_h(full_B(s))).reshape(-1,2))
    naive=ZE(np.moveaxis(s,0,1))
    naive=np.moveaxis(naive,1,0)
    naive_after=ZE(np.moveaxis(full_h(s),0,1))
    naive_after=np.moveaxis(naive_after,1,0)
    naive_commutator=norm((naive_after-full_h(naive)).reshape(-1,2))
    assert b_square_error<1e-12 and commuting_error<1e-12 and naive_commutator>.01
    # Intertwining tested on all columns of an unknown qubit, not one state.
    initial=hist[0];evolved=B_at(initial,0)
    export_error=0.
    for l in range(31):
        if l:evolved=old.act(evolved,gates[l-1])
        export_error=max(export_error,norm(evolved-B_at(hist[l],l)))
    assert export_error<1e-12
    # Every original declared fault is preserved. This is reuse, not a new code.
    inherited_error=0.
    for f1 in range(4):
        for f2 in range(4):
            _,h=old.history(theta,(f1,f2))
            inherited_error=max(inherited_error,norm(h[-1]-old.expected(theta,(f1,f2))))
    assert inherited_error<1e-12

    # Reuse the old field tables with an EXPLICIT transport error: those tables
    # contain h_old in the Newton mass. The reference below uses the base mass
    # only; the restored actual h correction is separately bounded.
    table=np.array(field['common_process']['joint_probability_tables'])
    assert table.shape==(2,8) and np.min(table)>0
    params=field['parameters']['inherited_gravity_and_moving_parameters']
    G=params['G'];mass=params['base_mass'];mr=params['receiver_mass']
    core=params['soft_core'];delta=params['position_std'];rmin=min(params['radii'])
    u_max=1/math.sqrt(rmin*rmin+core*core)
    table_transport=t*G*mr*30*u_max
    internal_redshift=t*G*mr*hmax*u_max
    # Complete CP instrument on Q plus any passive reference, jointly with
    # (x,y) and the three commuting physical effects z.
    ks=old.kraus(theta)
    bell=np.array([1.,0,0,1.],complex)/math.sqrt(2)
    choi_blocks=[];native_blocks=[];completeness=np.zeros((2,2),complex)
    instrument_error=0.;min_branch_probability=1.
    for x in (0,1):
        for y in (0,1):
            k=ks[y]@ks[x];completeness+=k.conj().T@k
            min_branch_probability=min(min_branch_probability,float(np.linalg.eigvalsh(k.conj().T@k).min()))
            v=np.kron(k,np.eye(2))@bell;block=np.outer(v,v.conj())
            added=[table[x,z]*block for z in range(8)]
            instrument_error=max(instrument_error,norm(sum(added)-block))
            choi_blocks.extend(added);native_blocks.append(block)
    assert norm(completeness-np.eye(2))<1e-12
    assert instrument_error<1e-12 and min_branch_probability>=.04-1e-13
    assert abs(sum(np.trace(b).real for b in choi_blocks)-1)<1e-12
    # Choi equality covers all unknown data and references. Testing particular
    # input states below is only to provide concrete signal numbers.
    probe_inputs=[old.ry(-theta)[:,q] for q in (0,1)]
    input_rows=[]
    for psi in probe_inputs:
        event=np.array([[np.linalg.norm(ks[y]@ks[x]@psi)**2 for y in (0,1)] for x in (0,1)])
        weight_x=event.sum(axis=1)
        joint=np.einsum('xy,xz->xyz',event,table)
        observed=joint.sum(axis=(0,1))
        input_rows.append(dict(event_probabilities=event.tolist(),
            first_event_probabilities=weight_x.tolist(),joint_physical_outcomes=observed.tolist(),
            record=float(observed[[2,3,6,7]].sum()),higgs=float(observed[[1,3,5,7]].sum()),
            path=float(observed[4:].sum())))
    imbalance=input_rows[0]['first_event_probabilities'][0]-input_rows[1]['first_event_probabilities'][0]
    assert abs(imbalance-.6)<1e-13
    record_gap=input_rows[0]['record']-input_rows[1]['record']
    higgs_gap=input_rows[0]['higgs']-input_rows[1]['higgs']
    assert abs(record_gap-.6*field['common_process']['record_probability_contrast'])<1e-12
    assert abs(higgs_gap-.6*field['common_process']['higgs_probability_contrast'])<1e-12
    # A common finite readout window, bounded using the same autonomous H and
    # initial vacuum, rather than holding the classical phase fixed by hand.
    half_window=1e-5
    scalar_min=math.sqrt(min(field['parameters']['mass_squared']))
    I0=1/(8*np.pi**1.5)
    initial_kinetic_norm=2*math.sqrt(15)/(8*mass*delta*delta)
    interaction_vacuum_norm=math.sqrt(I0/(2*scalar_min)) # |gS|+|gR|=1
    newton_norm=G*mr*(mass+hmax)/core
    initial_H_norm=hmax+initial_kinetic_norm+interaction_vacuum_norm+newton_norm
    time_error=half_window*initial_H_norm
    movement=field['finite_moving_control']['state_error_GN'] # safe for h<=pi<30
    total=movement+internal_redshift+time_error
    numerical_table_budget=table_transport+2e-12 # loose diagnostic display budget, not a new physical hypothesis
    higgs_lower=.6*field['finite_moving_control']['higgs_static_contrast_lower']-2*total
    record_lower=.6*field['finite_moving_control']['record_static_contrast_lower']-2*total
    gravity_static=(field['finite_moving_control']['gravity_moving_contrast_lower']
        +field['finite_moving_control']['state_error_GN']+field['finite_moving_control']['state_error_G0'])
    gravity_lower=gravity_static-2*total
    assert min(higgs_lower,record_lower,gravity_lower)>0
    conditional_bound=2*total/(.04-total)
    assert total<.000635 and conditional_bound<.033
    files=[Path(__file__),STAGE/'929/joint_protocol_selection.py',
        STAGE/'929/joint_protocol_selection_results.json',STAGE/'946/portal_common_process_results.json']
    return dict(round=947,date='2026-10-07',all_scientific_checks_passed=True,
        same_physical_protocol=dict(theta=float(theta),duration=t,internal_clock_scale=float(scale),
            internal_mass_energy_range=[0.,float(hmax)],internal_fault_sectors=16,repeats=2,
            original_work_and_fault_and_clock_dimension=4063232,
            autonomous_clock_finish_error=clock_error,event_coupling_square_error=b_square_error,
            event_coupling_h_commutator_error=commuting_error,
            naive_undressed_event_h_commutator_norm=naive_commutator,
            complete_event_intertwining_error=export_error,
            original_fault_protocol_transport_error=inherited_error,
            coupling_changes_only_clock_slots_0_and_1=True),
        joint_instrument=dict(choi_branch_recovery_error=instrument_error,
            effect_completeness_error=norm(completeness-np.eye(2)),
            minimum_event_effect_eigenvalue=min_branch_probability,
            unknown_input_and_passive_reference_identity=True,physical_outcome_condition_depends_only_on_first_event=True,
            input_examples=input_rows,record_probability_contrast=record_gap,higgs_probability_contrast=higgs_gap,
            recycled_field_table_transport_bound=table_transport,
            displayed_probability_budget_vs_base_mass_reference=numerical_table_budget),
        finite_joint_control=dict(inherited_two_field_moving_bound=movement,
            actual_internal_energy_redshift_bound=internal_redshift,
            readout_window_half_width=half_window,initial_H_norm_upper=initial_H_norm,
            autonomous_time_window_error=time_error,full_cq_instrument_distance_bound=total,
            all_unknown_input_and_declared_faults_covered=True,
            postevent_normalized_private_state_distance_upper=conditional_bound,
            record_probability_contrast_lower=record_lower,higgs_probability_contrast_lower=higgs_lower,
            gravity_probability_contrast_lower=gravity_lower),
        scope=dict(original_929_finite_six_protocol_mathematical_contract_jointly_transported=True,
            event_sensitive_material_coupling_and_timescale_are_inputs=True,
            terminal_logical_access_permissions_are_inherited_inputs=True,
            field_readout_of_actual_first_event_is_nonzero=True,
            autonomous_construction_of_all_access_instruments_certified=False,
            full_SM_Einstein_parent_matching_certified=False,
            spatially_separated_organizations_or_all_attacks_certified=False,
            arbitrary_long_repetition_claimed=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-10,(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
