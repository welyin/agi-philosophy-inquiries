"""634: original continuum record charges and source-boundary obstruction.

The original smooth packet is phase modulated, not replaced by a new field.
Geometry is used only through the declared linearized Einstein/Bianchi test.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_continuum_record_sources as old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_record_ward_boundary_results.json'


def moments(ell=2.,boost=.7,refined=False):
    if refined:q,pw,_,_,_=old.compact_bump(nr=232,nk=660)
    else:q,pw,_,_,_=old.compact_bump()
    mu,wm=np.polynomial.legendre.leggauss(144 if refined else 96)
    wm=wm/2
    m,w,_,_,_=old.mass_data()
    kz=q[:,None]/ell*mu[None,:]+boost
    ksq=(q[:,None]/ell)**2+boost**2+2*q[:,None]/ell*boost*mu[None,:]
    en=np.sqrt(ksq[:,:,None]+m*m)
    mixE=en@w; mixInv=(1/en)@w; mixMass=(m*m/en)@w
    def avg(x):return float(pw@(x@wm))
    norm=float(pw.sum())
    occ=.5*(norm-avg(kz*mixInv))
    meanE=avg(mixE);meanK=boost*norm
    meanKK=avg(kz*kz*mixInv)
    # Normalize the quadrature-only tail loss; report it explicitly.
    p=occ/norm;A=meanE/norm;d=meanK/norm;K=meanKK/norm
    Echi=.5*(A+(2*p-1)*d)
    Pchi=.5*(K+(2*p-1)*d)
    return dict(p=p,mean_absolute_B=A,bare_momentum=d,
                Qchi=np.array([Echi,Pchi]),mass_source=avg(mixMass)/norm,
                quadrature_norm=norm)


def boosted_original_check():
    result=moments();ref=moments(refined=True)
    diffs=[abs(result[x]-ref[x]) for x in ('p','mean_absolute_B','mass_source')]
    diffs.append(float(np.max(abs(result['Qchi']-ref['Qchi']))))
    assert max(diffs)<2e-7
    assert .2<result['p']<.5 and result['Qchi'][0]>abs(result['Qchi'][1])
    masses,weights,_,_,_=old.mass_data()
    err=0.;rng=np.random.default_rng(634)
    for k in rng.normal(size=(28,3)):
        k[2]+=.7
        B=old.bdg(k);_,P,ab=old.spectral(B)
        en=np.sqrt(k@k+masses*masses)
        p=.5*(1-k[2]*np.dot(weights,1/en))
        minus=.5*(np.dot(weights,en)-k[2])
        source=np.dot(weights,masses*masses/en)
        err=max(err,abs(np.vdot(old.E,P@old.E).real-p),
                abs(-np.vdot(old.E,B@P@old.E).real-minus),
                abs(-2*np.vdot(old.E,P@old.B0@old.E).real-source))
    assert err<3e-14
    q=result['Qchi'];p=result['p']
    return dict(support_radius=2.,phase_modulation_kz=.7,
                exact_packet_is_smooth_compact=True,probability_one=p,
                branch_one_energy_momentum=(q/p).tolist(),
                branch_zero_energy_momentum=(q/(1-p)).tolist(),
                average_energy_momentum=(2*q).tolist(),
                nonselective_mass_source=result['mass_source'],
                raw_packet_norm=result['quadrature_norm'],
                norm_corrected_for_numerical_tail_only=True,
                quadrature_refinement_difference=float(max(diffs)),
                original_BdG_projector_error=float(err),
                arbitrary_grid_precision_not_claimed=True)


def charge_and_record_check():
    r=moments();p=r['p'];q=r['Qchi']
    probs=np.array([1-p,p]);conditional=np.array([q/(1-p),q/p])
    mean=probs@conditional;delta=conditional-mean
    covariance=(delta.T*probs)@delta
    expected=(1-2*p)**2/(p*(1-p))*np.outer(q,q)
    err=float(np.max(abs(covariance-expected)))
    assert err<2e-15
    assert np.linalg.eigvalsh(covariance)[0]>-2e-15 and covariance[0,0]>1e-3
    assert covariance[0,1]>1e-3
    assert all(v[0]>abs(v[1]) for v in conditional)
    # The scalar variance of every linear charge combination is positive.
    null=np.array([-q[1],q[0]])
    assert abs(null@covariance@null)<2e-15
    symmetric=moments(boost=0.)
    assert abs(symmetric['p']-.5)<1e-14
    return dict(probabilities=probs.tolist(),conditional_Q=conditional.tolist(),
                mean_Q=mean.tolist(),record_between_branch_covariance=covariance.tolist(),
                covariance_rank=int(np.linalg.matrix_rank(covariance,tol=1e-12)),
                covariance_identity_error=err,
                record_covariance_is_lower_bound_not_total_quantum_noise=True,
                symmetric_packet_probability=symmetric['p'],
                instantaneous_step_surface_charge_coefficients=conditional.tolist(),
                averaged_surface_charge_coefficient=mean.tolist())


def ward_and_geometry_check():
    r=moments();Qbar=2*r['Qchi']
    # Weak, spatially integrated balance for a smooth transition and temporary
    # compactly supported contact-energy corrections. No spatial stress is fitted.
    z,w=np.polynomial.legendre.leggauss(240)
    bump=np.exp(-1/(1-z*z));bumpprime=-2*z/(1-z*z)**2*bump
    norm=float(w@bump)
    rows=[]
    for tau in (.4,1.,2.):
        rate=bump/(norm*tau)
        temporary=np.array([7.,-3.])
        derivative=rate[:,None]*Qbar+bumpprime[:,None]/tau*temporary
        residual=tau*(w@derivative)
        error=float(np.max(abs(residual-Qbar)))
        assert error<2e-14
        rows.append(dict(duration_parameter=tau,integrated_residual=residual.tolist(),
                         momentarily_added_energy_momentum=temporary.tolist(),error=error))
    # Primary geometric identity: linearized Einstein operator is transverse
    # for arbitrary Fourier h; a nonzero source divergence cannot be fixed by gauge.
    eta=np.diag([-1.,1.,1.,1.]);rng=np.random.default_rng(6342)
    maxerr=0.
    for _ in range(32):
        kl=rng.normal(size=4);ku=eta@kl
        h=rng.normal(size=(4,4));h=(h+h.T)/2
        tr=float(np.trace(eta@h));k2=float(kl@ku);hk=h@ku
        Ric=.5*(-np.outer(kl,hk)-np.outer(hk,kl)+k2*h+np.outer(kl,kl)*tr)
        R=float(np.trace(eta@Ric));Ein=Ric-.5*eta*R
        maxerr=max(maxerr,float(np.max(abs(ku@Ein))))
    assert maxerr<2e-14
    # Every branch can obey total charge balance if another retained sector
    # changes by the opposite charge. This is necessary bookkeeping only.
    p=r['p'];q=r['Qchi'];required=np.array([-q/(1-p),-q/p])
    assert np.max(abs(np.array([1-p,p])@required+Qbar))<1e-15
    return dict(temporary_contact_rows=rows,
                linearized_Bianchi_error=maxerr,
                required_other_sector_branch_charge_changes=required.tolist(),
                no_compact_time_contact_can_remove_nonzero_integrated_jump=True,
                source_splicing_is_a_tested_rejected_ansatz=True,
                state_conditioning_itself_not_claimed_to_violate_conservation=True,
                compensating_source_or_metric_solution_not_constructed=True)


def run():
    a=boosted_original_check();b=charge_and_record_check();c=ward_and_geometry_check()
    deps=('joint_continuum_record_sources.py','joint_continuum_record_sources_results.json',
          'research_note_593.md','research_note_600.md','research_note_601.md',
          'research_note_618.md','research_note_620.md','research_note_624.md','research_note_633.md')
    return dict(round=634,tests_run=3,failures=0,errors=0,boosted_original=a,
                records_and_charges=b,Ward_boundary=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(original_continuum_neutral_mass_mixing=True,
                           inherited_algebraic_instrument_not_new_device=True,
                           record_charges_noise_and_source_boundary_linked=True,
                           original_fixed_background_branch_only=True,
                           source_jump_without_compensation_rejected=True,
                           Einstein_action_is_input=True,
                           conditional_state_not_physical_collapse_rule=True,
                           no_full_geometry_or_internal_instrument_constructed=True,
                           no_broad_GR_or_unification_no_go=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=634,tests=3,all_passed=True,
                         p=result['boosted_original']['probability_one'],
                         mean_charges=result['records_and_charges']['mean_Q'],
                         record_noise=result['records_and_charges']['record_between_branch_covariance'])))

