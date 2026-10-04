"""552: one algebraic coefficient controls curved-vacuum response and radial inertia.

RG endpoint decimal proxies are exact models, not rigorous ODE enclosures.
No curved-space beta functions or quantum effective vacuum are assumed solved.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_vacuum_hierarchy_matching as prior
import joint_reference_gravity_constraints as frame

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_split_mass_geometry_results.json'


def inverse_product(L,b):
    a,c,d=L
    det=a*d-c*c
    assert a>0 and det>0
    return [(d*b[0]-c*b[1])/det,(a*b[1]-c*b[0])/det]


def exact_background(L,C,M02,R):
    u=inverse_product(L,C)
    v=inverse_product(L,[Q(1),Q(1)])
    A=sum(v)
    F0=M02-sum(u)/6
    x=[u[i]-R*v[i]/6 for i in range(2)]
    F=M02-sum(x)/6
    energy=sum(C[i]*u[i] for i in range(2))
    V0=(F0*R+energy)/4
    a,p,b=L
    V=V0-(C[0]*x[0]+C[1]*x[1])/2+(a*x[0]**2+2*p*x[0]*x[1]+b*x[1]**2)/4
    assert a*x[0]+p*x[1]==C[0]-R/6
    assert p*x[0]+b*x[1]==C[1]-R/6
    assert F*R==4*V and F==F0+R*A/36
    assert M02*R==4*V0-sum(C[i]*x[i] for i in range(2))
    out=dict(L=L,C=C,M02=M02,R=R,u=u,v=v,A=A,F0=F0,x=x,F=F,V0=V0,V=V)
    if F!=0:
        factor=1-R*A/(36*F)
        z=R/(36*F)
        assert factor==F0/F
        assert (a-z)*(b-z)-(p-z)**2==(a*b-p*p)*factor
        out['factor']=factor
    return out


def scalar_frame(b):
    assert min(b['x'])>0 and b['F']>0
    p=dict(lh=b['L'][0],p=b['L'][1],ls=b['L'][2],a=b['C'][0],b=b['C'][1],
           M02=b['M02'],V0=b['V0'])
    data=dict(h2=b['x'][0],s2=b['x'][1],R=b['R'])
    return p,prior.stability_data(p,data)


def string_model(b):
    return {key:([str(x) for x in value] if isinstance(value,list) else str(value))
            for key,value in b.items()}


def running_inputs():
    saved=json.loads((HERE/'joint_singlet_common_mass_rg_results.json').read_text('utf8'))
    q=lambda x:Q(str(x))
    boundary=saved['common_boundary']
    C0=Q(1,4)
    M02=C0*(3+q(boundary['r']))/q(boundary['T'])  # n_g=3; declared fixed-endpoint input
    outputs=[]
    for index in (2,4,5):
        row=saved['examples'][index]
        s=row['state']
        L=[q(s[k]) for k in ('lambda_H','p','lambda_s')]
        d=[q(s[k]) for k in ('x','y')]
        w=inverse_product(L,d)
        assert min(w)>0
        kc=6*M02/(C0*sum(w))
        outputs.append(dict(source_index=index,q0=row['q0'],mu_GeV=row['mu_GeV'],
                            C0=C0,M02=M02,L=L,d=d,w=w,kappa_critical=kc,
                            source_step_doubling_error=row['step_doubling_max_error']))
    return outputs


def run():
    checks=[]
    base=prior.parameters()
    L=[base[k] for k in ('lh','p','ls')]
    M02=base['M02']
    admissible=[]
    rejected=0
    for C,R in itertools.product(([Q(3,10),Q(1,5)],[Q(9,5),Q(17,10)],[Q(3,2),Q(3,2)]),
                                 (Q(-1,10),Q(0),Q(1,10),Q(1),Q(3),Q(6),Q(8))):
        b=exact_background(L,C,M02,R)
        if min(b['x'])>0 and b['F']>0:
            if b['F0']!=0:
                assert R==(4*b['V0']-sum(C[i]*b['u'][i] for i in range(2)))/b['F0']
            if b['F0']<=0: assert R>0
            admissible.append(string_model(b))
        else: rejected+=1
    assert admissible and rejected
    checks.append('unequal_mass_full_equations_admissible_domain_and_curvature_classification')

    examples=[exact_background(L,[Q(3,10),Q(1,5)],M02,Q(1,10)),
              exact_background(L,[Q(9,5),Q(17,10)],M02,Q(3)),
              exact_background(L,[Q(3,2),Q(3,2)],M02,Q(3))]
    derivative_error=0.
    algebra_error=0.
    inertia=[]
    for b in examples:
        p,(phi,d,hu,hr,masses)=scalar_frame(b)
        algebra_error=max(algebra_error,float(np.max(abs(hu-hr)))/(1+float(np.max(abs(hu)))))
        fd=np.zeros((2,2))
        for j in range(2):
            delta=np.eye(2)[j]*1e-6
            fd[:,j]=(frame.fields(phi+delta,p)['dU']-frame.fields(phi-delta,p)['dU'])/2e-6
        derivative_error=max(derivative_error,float(np.max(abs(fd-hu)))/(1+float(np.max(abs(hu)))))
        if b['F0']>0: assert masses[0]>0
        elif b['F0']<0: assert masses[0]<0<masses[1] and b['R']>0
        else: assert abs(masses[0])<1e-9 and masses[1]>0
        inertia.append(dict(exact=string_model(b),canonical_radial_mass_squared=masses.tolist()))
    assert derivative_error<1e-7 and algebra_error<1e-12
    checks.append('independent_Einstein_Hessian_and_all_three_radial_inertia_cases')

    degenerate=[]
    for R in (Q(1),Q(3),Q(6),Q(8)):
        b=exact_background(L,[Q(3,2),Q(3,2)],M02,R)
        assert min(b['x'])>0 and b['F']>0 and b['F0']==0
        assert b['V0']==3 and b['V']/b['F']**2==9/b['A']
        p,(phi,d,hu,hr,masses)=scalar_frame(b)
        tangent=-np.array(list(map(float,b['v'])))/(12*phi)
        assert np.max(abs(hu@tangent))<1e-9
        degenerate.append(dict(R=str(R),F=str(b['F']),U=str(b['V']/b['F']**2)))
    # The complete-square inequality is checked away from the stationary line too.
    cauchy=[]
    b=examples[2]
    for x in ([Q(1),Q(1)],[Q(1),Q(2)],[Q(4),Q(4,3)]):
        y=[x[i]-b['u'][i] for i in range(2)]
        F=M02-sum(x)/6
        assert F>0
        qform=L[0]*y[0]**2+2*L[1]*y[0]*y[1]+L[2]*y[1]**2
        gap=b['A']*qform-sum(y)**2
        assert gap>=0
        U=qform/(4*F*F)
        assert U>=9/b['A']
        cauchy.append(dict(x=list(map(str,x)),positive_Cauchy_gap=str(gap),U=str(U)))
    checks.append('degenerate_branch_is_an_exact_flat_valley_on_its_feasible_domain')

    response=[]
    delta=Q(1,1000)
    for b in examples[:2]:
        shifted=exact_background(L,b['C'],M02,b['R']+4*delta/b['F0'])
        assert shifted['V0']-b['V0']==delta
        assert min(shifted['x'])>0 and shifted['F']>0
        response.append(dict(F0=str(b['F0']),dR_dV0=str(4/b['F0'])))
    flat=examples[2]
    # If F0=0, any unmatched V0 leaves a nonzero residual independent of R.
    for R in (Q(-2),Q(0),Q(2)):
        residual=flat['F0']*R-(4*(flat['V0']+delta)-sum(flat['C'][i]*flat['u'][i] for i in range(2)))
        assert residual==-4*delta
    checks.append('same_coefficient_controls_susceptibility_and_singular_unmatched_no_solution')

    running=[]
    endpoints=running_inputs()
    for endpoint in endpoints:
        for kappa in (Q(1),Q(1,100),endpoint['kappa_critical']/2):
            C=[kappa*endpoint['C0']*x for x in endpoint['d']]
            b=exact_background(endpoint['L'],C,endpoint['M02'],Q(0))
            assert min(b['x'])>0 and b['F0']==b['F']>0 and b['V']==0
            assert b['x']==[kappa*endpoint['C0']*x for x in endpoint['w']]
            p,(phi,d,hu,hr,masses)=scalar_frame(b)
            assert masses[0]>0
            assert b['x'][1]/b['x'][0]==endpoint['w'][1]/endpoint['w'][0]
            running.append(dict(source_index=endpoint['source_index'],q0=endpoint['q0'],
                mu_GeV=endpoint['mu_GeV'],kappa=str(kappa),kappa_critical=float(endpoint['kappa_critical']),
                F=float(b['F']),F_over_h2=float(b['F']/b['x'][0]),s2_over_h2=float(b['x'][1]/b['x'][0]),
                canonical_radial_mass_squared=masses.tolist(),
                exact_proxy_equations_verified=True,source_step_doubling_error=endpoint['source_step_doubling_error']))
    checks.append('frozen_running_endpoint_proxies_satisfy_joint_flat_geometry_and_radial_conditions')

    target=Q(10**12)
    ceilings=[e['M02']/(e['C0']*(target*e['w'][0]+sum(e['w'])/6)) for e in endpoints]
    common_kappa=min(ceilings)/2
    simultaneous=[]
    for e in endpoints:
        b=exact_background(e['L'],[common_kappa*e['C0']*d for d in e['d']],e['M02'],Q(0))
        assert b['F']>0 and b['F']/b['x'][0]>=target
        assert common_kappa<e['kappa_critical']
        # With kappa above critical, flat F fails; curved positive-F solutions, if any, have F0<0.
        bad=exact_background(e['L'],[2*e['kappa_critical']*e['C0']*d for d in e['d']],e['M02'],Q(0))
        assert bad['F0']<0 and bad['F']<0
        simultaneous.append(dict(source_index=e['source_index'],F_over_h2=float(b['F']/b['x'][0]),
                                 kappa_critical=float(e['kappa_critical'])))
    checks.append('one_common_matching_scale_works_for_all_listed_endpoints_and_has_explicit_target_bound')

    # A positive sum and F0 do not replace componentwise physical condensate positivity.
    phase_L=[Q(2),Q(1),Q(1)]
    w=inverse_product(phase_L,[Q(1),Q(3)])
    assert w==[Q(-2),Q(5)] and sum(w)>0
    wrong=exact_background(phase_L,[Q(1,10),Q(3,10)],Q(1),Q(0))
    assert wrong['F0']>0 and min(wrong['x'])<0
    # New theorem specializes to the frozen 551 branch, without modifying that branch.
    for t in (Q(1,10),Q(1),Q(3)):
        p,old=prior.background(base,t)
        b=exact_background(L,[base['C0'],base['C0']],M02,old['R'])
        assert b['x']==[old['h2'],old['s2']] and b['V0']==old['V0']
        assert b['factor']==Q(2*base['ng']-1)/(2*base['ng']-t)
    checks.append('missing_phase_condition_counterexample_and_exact_prior_branch_reduction')

    deps=('research_note_544.md','joint_singlet_common_mass_rg.py','joint_singlet_common_mass_rg_results.json',
          'research_note_551.md','joint_vacuum_hierarchy_matching.py','joint_vacuum_hierarchy_matching_results.json',
          'joint_reference_gravity_constraints.py','research_round_551_checks.json')
    return dict(round=552,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_admissible_examples=admissible,rejected_algebraic_examples=rejected,
        inertia_examples=inertia,maximum_Hessian_relative_difference=derivative_error,
        maximum_analytic_Hessian_relative_residual=algebra_error,
        exact_degenerate_valley=degenerate,off_valley_Cauchy_certificates=cauchy,
        vacuum_response_examples=response,running_endpoint_examples=running,
        declared_endpoint_gravity_coefficient=str(endpoints[0]['M02']),
        common_target_ratio=str(target),common_kappa=str(common_kappa),common_kappa_approx=float(common_kappa),
        simultaneous_endpoint_certificate=simultaneous,phase_counterexample=string_model(wrong),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(fixed_four_dimensional_classical_two_derivative_action=True,
            L_positive_definite_and_all_condensates_nonzero=True,common_conformal_curvature_coupling_input=True,
            curvature_uniqueness_and_radial_inertia_use_same_coefficient=True,
            flat_valley_only_on_nonempty_admissible_domain=True,radial_linear_sector_not_full_stability=True,
            running_endpoint_decimals_are_exact_proxy_models_not_ODE_enclosures=True,
            mass_homogeneity_inherited_not_new=True,finite_endpoint_list_not_entire_RG_trajectory=True,
            fixed_endpoint_gravity_and_vacuum_matching_are_declared_inputs=True,
            quantum_effective_vacuum_and_curved_space_running_not_completed=True,
            small_common_mass_not_derived=True,completed_unification=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else: assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','maximum_Hessian_relative_difference','common_kappa_approx')}))
