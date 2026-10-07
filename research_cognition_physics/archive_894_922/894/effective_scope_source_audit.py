"""894: finite source witness and effective-model scope audit.
The smooth prescribed-history comparison is NOT an autonomous quantum-coordinate
comparison. No continuum or physical minimum length is asserted by this test.
"""
from pathlib import Path
from itertools import combinations
import json, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import actual_source_transport_probe as working
TARGET=HERE/'effective_scope_source_audit_results.json'


def positive_readout(U,V):
    # The original filled negative-energy Slater state; evaluate positive
    # transition amplitudes, never subtract nearly equal vacuum traces.
    W=V.conj().T@U@V[:,:2]
    n=float(np.linalg.norm(W[2:,:],'fro')**2)
    event=0.; fock_n=0.; norm=0.
    for ix in combinations(range(4),2):
        weight=float(abs(np.linalg.det(W[list(ix),:]))**2)
        norm+=weight
        positives=sum(i>=2 for i in ix)
        if positives:event+=weight
        fock_n+=positives*weight
    assert abs(norm-1)<2e-11
    assert abs(fock_n-n)<2e-12
    return dict(pair_event_probability=event,positive_particle_number=n,
                bounded_effect_probability=n/2,slater_norm_error=abs(norm-1),
                fock_vs_one_body_number_error=abs(fock_n-n))


def run():
    old=working.old;prior=working.prior
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.prior.prior.prior.last.em_charges(m,matter)
    ids=[26,27,28,29];a0=-1/12
    full=old.H(m,q,17,[0,0,0],a0,kind='continuum')
    h0=full[np.ix_(ids,ids)]
    G=(m.GAMMA[0]@np.diag(q))[np.ix_(ids,ids)]
    ev,V=np.linalg.eigh(h0);E=float(ev[-1])
    P=V[:,:2]@V[:,:2].conj().T
    coupling=float(np.linalg.norm(V[:,2:].conj().T@G@V[:,:2],'fro')**2)
    assert coupling>0 and 0<E<np.pi/2
    coeffs=[]
    for order in (80,160,320):
        xs,ws=np.polynomial.legendre.leggauss(order)
        ft=sum(w/2*working.pulse((x+1)/2)[0]*np.exp(1j*E*(x+1)) for x,w in zip(xs,ws))
        coeffs.append(float(coupling*abs(ft)**2))
    coefficient=coeffs[-1]
    # f >= exp(-1/3) on [1/4,3/4]; midpoint-rotated cosine >= cos(E)>0.
    lower_coefficient=float(coupling*(.5*np.exp(-1/3)*np.cos(E))**2)
    identity_errors=[]; finite_cutoff_errors=[]
    for a in (a0,a0+.001,a0+.004):
        h=h0+(a-a0)*G;en=np.sqrt(np.trace(h@h).real/4)
        enp=np.trace(h@G).real/(4*en)
        S=h/en;Sp=G/en-h*enp/en**2
        Pa=-Sp/2;Proj=(np.eye(4)-S)/2
        A=Pa@Proj-Proj@Pa
        identity_errors.append(float(np.max(abs(A@Proj-Proj@A-Pa))))
        for N in (17,33):
            finite=old.H(m,q,N,[0,0,0],a,kind='buffer')
            finite_cutoff_errors.append(float(np.max(abs(finite[np.ix_(ids,ids)]-h))))
    assert max(identity_errors)<2e-13
    assert max(finite_cutoff_errors)<2e-13
    rows=[]
    for eps in (.004,.002,.001):
        for steps in (1024,2048):
            row=dict(amplitude=eps,time_steps=steps)
            for transport,label in ((False,'original'),(True,'horizontal')):
                U=working.propagate(h0,G,eps,steps,transport)
                row[label]=positive_readout(U,V)
            # Dyson bound uses integral|f| <= 1 and ||G||=6. Numerical
            # evaluation of an analytic bound is not an interval certificate.
            b=abs(eps)*float(np.linalg.norm(G,2))
            r=float(np.expm1(b)-b)
            n_lower=max(0.,abs(eps)*np.sqrt(lower_coefficient)-np.sqrt(2)*r)**2
            row['analytic_event_lower_bound']=float(n_lower/2)
            row['original_event_over_amplitude_squared']=row['original']['pair_event_probability']/eps**2
            row['original_energy_deposited']=2*E*row['original']['positive_particle_number']
            assert row['original']['pair_event_probability']>row['analytic_event_lower_bound']>0
            assert row['horizontal']['pair_event_probability']<1e-14
            rows.append(row)
    assert abs(rows[-1]['original_event_over_amplitude_squared']-coefficient)<.002*coefficient
    step_error=max(abs(rows[j]['original']['pair_event_probability']-rows[j+1]['original']['pair_event_probability']) for j in (0,2,4))
    return dict(round=894,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3679,
        argument_scope='Audit finite-domain common-model acceptance, with one original finite-mode prescribed-EM-history witness. No universal autonomous no-go.',
        fixed_physical_protocol=dict(pulse_duration=1.,base_alpha=a0,momentum=[0,0,0],particle_indices=ids,gap=2*E),
        perturbative_pair_event_coefficient=coefficient,coefficient_quadrature_refinement=coeffs,
        strictly_positive_analytic_lower_coefficient=lower_coefficient,
        projector_identity_max_error=max(identity_errors),
        cutoff_17_33_same_retained_generator_error=max(finite_cutoff_errors),
        actual_finite_history_rows=rows,largest_time_step_refinement_change=step_error,
        audit_classification={
            'not_stage_prerequisites':['single_global_normal_folium_882','uniform_bare_coordinate_budget_888','exact_full_continuous_Gibbs_limit_889_893','arbitrary_point_or_zero_gap_region_resolution'],
            'must_bound_if_in_declared_protocol_domain':['preparation_and_record_probabilities','first_second_source_jets_and_contacts','Ward_constraints_and_backreaction','holonomy_dependent_vacuum_sources_886','cut_band_leakage_and_source_variance_884','horizontal_same_source_response_894'],
            'reuse_without_reproving':['875_876_fixed_graph_local_capacity_C2','881_finite_correlated_mode_source_approximation','882_finite_physical_CAR_C2','735_joint_one_loop_Ward','851_compatible_initial_mean']},
        microscopic_continuity_assumed=False,physical_minimum_scale_assumed=False,
        finite_coupling_remainder_for_full_E_proved=False,
        full_common_model_completed=False,visual_checks_performed=False)

if __name__=='__main__':
    result=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
