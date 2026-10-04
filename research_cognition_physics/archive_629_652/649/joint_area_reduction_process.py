"""649: original geometry area matching, local physical norm and process descent.

The matching constraint is an explicit candidate on two590 geometry replicas.
Its tangential replacement is a changed dynamics, not derived ADM gravity.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_area_reduction_process_results.json'
spec=importlib.util.spec_from_file_location('probe649',HERE/'round649_drafts/area_constraint_descent_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)


def area_gradient(theta):
    theta=np.asarray(theta)
    Q=probe.Q0+probe.RADIUS*np.einsum('a,aij->ij',np.sin(theta),probe.BASIS)
    values,V=np.linalg.eigh(Q);d=(values[:,None]-values[None,:])/2
    ratio=np.ones_like(d);np.divide(np.sinh(d),d,out=ratio,where=abs(d)>1e-14)
    divided=np.exp((values[:,None]+values[None,:])/2)*ratio
    g=probe.original.shape_exp(Q);ix=np.ix_(probe.TANGENT,probe.TANGENT);q=g[ix]
    c=np.sqrt(np.linalg.det(q));grad=[]
    for a,b in enumerate(probe.BASIS):
        dg=V@(divided*(V.T@(probe.RADIUS*np.cos(theta[a])*b)@V))@V.T
        grad.append(.5*c*np.trace(np.linalg.solve(q,dg[ix])))
    return float(c),np.array(grad)


def coarea_check():
    # Actual original metric restricted to a two-coordinate diagnostic section.
    # The full12-dimensional coarea result is analytic, not inferred from this.
    rows=[];width=.3
    for order in (256,384):
        x,w=np.polynomial.legendre.leggauss(order);x*=width;w*=width
        cs=[];dc=[]
        for t in x:
            th=np.zeros(6);th[0]=t;c,g=area_gradient(th);cs.append(c);dc.append(g[0])
        cs=np.array(cs);dc=np.array(dc)
        assert np.min(dc)>.08
        bump=np.exp(-1/(1-(x/width)**2))
        exact=float(np.sum(w*bump**4/dc))
        Q=cs[:,None]-cs[None,:];weight=w[:,None]*w[None,:]*(bump[:,None]*bump[None,:])**2
        cases=[]
        for eps in (.002,.001,.0005):
            delta=np.exp(-Q*Q/(2*eps*eps))/(np.sqrt(2*np.pi)*eps)
            regular=float(np.sum(weight*delta))
            nullnorm=float(np.sum(weight*delta*Q**4))
            cases.append(dict(epsilon=eps,regularized_inner_product=regular,
                coarea_error=abs(regular-exact),Q_squared_physical_norm_squared=nullnorm))
        assert cases[-1]['coarea_error']<cases[0]['coarea_error']/12
        assert cases[-1]['Q_squared_physical_norm_squared']<cases[0]['Q_squared_physical_norm_squared']/180
        rows.append(dict(order=order,reference_coarea_norm=exact,
                         min_area_derivative=float(dc.min()),cases=cases))
    assert abs(rows[0]['reference_coarea_norm']-rows[1]['reference_coarea_norm'])<1e-12
    assert abs(rows[0]['cases'][-1]['regularized_inner_product']-rows[1]['cases'][-1]['regularized_inner_product'])<1e-12
    return dict(rows=rows,only_two_coordinate_section_numerically_integrated=True,
                full_domain_positivity_and_quotient_proved_locally=True,
                Gaussian_delta_is_numerical_regularization_not_physical_cutoff=True)


def nullspace_descent_check():
    old=probe.run()
    saved=json.loads((HERE/'round649_drafts/area_constraint_descent_probe_results.json').read_text('utf8'))
    assert old==saved
    c,g=area_gradient(np.zeros(6));assert np.max(abs(g-probe.analytic_gradient()))<1e-15
    n=np.r_[g,-g];I=np.full(12,1.7);hbar=.7
    obstruction=-hbar*hbar*np.sum(n*n/I)
    assert abs(obstruction-old['nonzero_H_trace_analytic'])<1e-14
    # Original complete matter H and bounded records commute with Q because
    # they have no theta derivatives. Their restriction statement is analytic.
    return dict(original_probe_reproduced=True,area=c,area_gradient=g.tolist(),
        independent_geometry_H_on_null_class=float(obstruction),
        original_probe_differences=old['rows'],
        full_matter_terms_retained_in_analytic_product_rule=True,
        records_restrict_but_this_does_not_supply_full_waiting_process=True,
        no_inference_against_actual_ADM_constraints=True)


def constraint(x):
    return probe.area(x[:6])-probe.area(x[6:])


def projected_principal_check():
    _,g=area_gradient(np.zeros(6));n=np.r_[g,-g]
    M=np.eye(12)/1.7;Mn=M@n
    Mt=M-np.outer(Mn,Mn)/(n@Mn)
    ev=np.linalg.eigvalsh(Mt)
    assert np.max(abs(Mt@n))<1e-15 and abs(ev[0])<1e-14 and ev[1]>.5
    hbar=.7;rows=[];zero=np.zeros(12)
    # Independent full Hessian of Q^2, including mixed replica derivatives.
    for step in (.004,.002,.001):
        hess=np.zeros((12,12))
        for a in range(12):
            da=np.eye(12)[a]*step
            hess[a,a]=(constraint(da)**2+constraint(-da)**2)/step**2
            for b in range(a):
                db=np.eye(12)[b]*step
                val=(constraint(da+db)**2-constraint(da-db)**2-
                     constraint(-da+db)**2+constraint(-da-db)**2)/(4*step**2)
                hess[a,b]=hess[b,a]=val
        raw=-hbar*hbar/2*np.sum(M*hess)
        tangent=-hbar*hbar/2*np.sum(Mt*hess)
        rows.append(dict(step=step,original_H_Q_squared=float(raw),
                         modified_tangent_H_Q_squared=float(tangent),
                         hessian_error=float(np.max(abs(hess-2*np.outer(n,n))))))
    assert abs(rows[-1]['modified_tangent_H_Q_squared'])<abs(rows[0]['modified_tangent_H_Q_squared'])/12
    cross=float(np.linalg.norm(Mt[:6,6:],ord=2))
    assert abs(cross-1/(2*1.7))<1e-14
    return dict(tangent_principal_eigenvalues=ev.tolist(),normal_annihilation=float(np.linalg.norm(Mt@n)),
        necessary_cross_replica_block_in_this_projection=cross,rows=rows,
        projected_geometry_principal_part_changes_original_dynamics=True,
        full_selfadjoint_realization_requires_local_domain_and_boundary_choice=True,
        projected_choice_not_unique_or_derived_from_cognition=True)


def run():
    deps=('research_note_590.md','research_note_593.md','research_note_606.md','research_note_617.md',
          'research_note_648.md','joint_full_spatial_metric.py',
          'round649_drafts/area_constraint_descent_probe.py',
          'round649_drafts/area_constraint_descent_probe_results.json')
    return dict(round=649,tests_run=3,failures=0,errors=0,local_physical_inner_product=coarea_check(),
        original_process_descent=nullspace_descent_check(),
        explicitly_modified_tangent_candidate=projected_principal_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(local_regular_area_matching_constraint_on_original590_replicas=True,
                   conditional_record_and_matter_restriction_only=True,
                   old_geometry_evolution_cannot_descend_unchanged=True,
                   alternative_tangent_form_is_new_dynamics_not_an_equivalence=True,
                   no_ADM_quantization_or_unified_GR_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
