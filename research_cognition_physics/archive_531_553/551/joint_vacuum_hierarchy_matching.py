"""Round 551: constant-background matching and the matter/gravity hierarchy.

Analytic proofs are in research_note_551.md. Finite grids check, not prove,
the universal statements. Only the declared two-derivative radial sector.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_reference_gravity_constraints as previous

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_vacuum_hierarchy_matching_results.json'


def parameters(q=Q(1,4), S=Q(3,28), r=Q(7,3), ng=3, C0=Q(1,4)):
    assert q > S >= 0 and r > 0 and ng >= 1 and int(ng) == ng and C0 > 0
    T = 3*q+r*S
    A = (3+r)/T
    return dict(q=q,S=S,r=r,ng=ng,C0=C0,T=T,A=A,
                M02=C0*A*ng/3,lh=(3*q*q+r*S*S)/T,p=S,ls=T/r,
                a=C0,b=C0,V0=Q(0),gw2=3*T/(ng*(3+r)),gc2=3*T/(4*ng))


def background(pars, t, kappa=Q(1)):
    """Effective scalar coefficient is t*C0; bare common mass is kappa*C0."""
    p = dict(pars)
    assert 0 < t < 2*p['ng'] and kappa > 0
    c0, A, ng = p['C0'], p['A'], p['ng']
    cm = kappa*c0
    R = 6*c0*(kappa-t)
    h2 = t*c0/p['q']
    s2 = t*c0*p['r']/p['T']*(1-p['S']/p['q'])
    rho2 = h2+s2
    F = p['M02']-rho2/6
    V0 = (p['M02']*R+cm*rho2)/4
    p.update(a=cm,b=cm,V0=V0)
    V = V0-cm*rho2/2+(p['lh']*h2*h2+2*p['p']*h2*s2+p['ls']*s2*s2)/4
    assert rho2 == t*c0*A and F == c0*A*(2*ng-t)/6
    return p, dict(t=t,kappa=kappa,R=R,h2=h2,s2=s2,rho2=rho2,F=F,V=V,V0=V0)


def verify_exact_equations(p, b):
    h2,s2,R,F,V = (b[k] for k in ('h2','s2','R','F','V'))
    assert -p['a']+p['lh']*h2+p['p']*s2+R/6 == 0
    assert -p['b']+p['ls']*s2+p['p']*h2+R/6 == 0
    # For constant scalars on an Einstein metric, this IS the full metric equation.
    assert F*R == 4*V
    assert p['M02']*R == 4*p['V0']-p['a']*(h2+s2)


def stability_data(p, b):
    phi = np.sqrt([float(b['h2']),float(b['s2'])])
    d = previous.fields(phi,p)
    h,s = phi
    F,R = d['F'],float(b['R'])
    L = np.array([[float(p['lh']),float(p['p'])],[float(p['p']),float(p['ls'])]])
    hV = np.array([[-float(p['a'])+3*L[0,0]*h*h+L[0,1]*s*s,2*L[0,1]*h*s],
                   [2*L[0,1]*h*s,-float(p['b'])+3*L[1,1]*s*s+L[0,1]*h*h]])
    # Differentiate U=V/F^2 BEFORE applying stationarity, independently of rank-one formula.
    df,ddf = -phi/3,-np.eye(2)/3
    hu_full = (hV/F**2-2*(np.outer(d['dV'],df)+np.outer(df,d['dV']))/F**3
               -2*d['V']*ddf/F**3+6*d['V']*np.outer(df,df)/F**4)
    hu_rankone = 2*np.diag(phi)@(L-R*np.ones((2,2))/(36*F))@np.diag(phi)/F**2
    chol = np.linalg.cholesky(d['G'])
    inverse = np.linalg.inv(chol)
    canonical = inverse@hu_full@inverse.T
    masses = np.linalg.eigvalsh(canonical)
    return phi,d,hu_full,hu_rankone,masses


def run():
    checks = []
    base = parameters()
    assert base['M02'] == previous.bare_parameters()['M02'] == Q(4,3)
    rows = []
    for t in (Q(11,10),Q(1),Q(9,10),Q(1,10),Q(1,100)):
        p,b = background(base,t)
        verify_exact_equations(p,b)
        assert b['V0'] == p['C0']**2*p['A']*(2*p['ng']-(2*p['ng']-1)*t)/4
        reconstructed = (2*p['ng']-4*b['V0']/(p['C0']**2*p['A']))/(2*p['ng']-1)
        assert reconstructed == t
        rows.append({k:str(b[k]) for k in ('t','R','h2','s2','F','V0')}
                    |dict(R_over_mu2=str(4*b['R']/p['C0']),F_over_h2=str(b['F']/b['h2'])))
    checks.append('exact_full_constant_scalar_and_Einstein_equations_and_inverse_vacuum_matching')

    certificate_count = 0
    for q,alpha,r,ng in itertools.product((Q(1,8),Q(1,2)),(Q(0),Q(1,3),Q(9,10)),
                                          (Q(1,2),Q(7,3),Q(7)),(1,3,5)):
        p0 = parameters(q,q*alpha,r,ng)
        determinant = p0['lh']*p0['ls']-p0['p']**2
        assert determinant == 3*q*q/r > 0
        inverse_sum = (p0['lh']+p0['ls']-2*p0['p'])/determinant
        assert inverse_sum == p0['A']
        assert 3+r == 4*p0['gc2']/p0['gw2'] and p0['M02']*p0['gw2'] == p0['C0']
        for t in (Q(1,10),Q(1),Q(2*ng)-Q(1,10)):
            p,b = background(p0,t)
            verify_exact_equations(p,b)
            factor = 1-b['R']*p['A']/(36*b['F'])
            assert factor == Q(2*ng-1)/(2*ng-t) > 0
            z = b['R']/(36*b['F'])
            det_updated = (p['lh']-z)*(p['ls']-z)-(p['p']-z)**2
            assert p['lh']-z > 0 and det_updated == determinant*factor > 0
            certificate_count += 1
    checks.append('exact_rank_one_positive_Hessian_certificate_and_same_scale_gauge_identities')

    derivative_error = 0.
    matrix_error = 0.
    stationary_error = 0.
    stable = []
    for t in (Q(1,100),Q(1,5),Q(1),Q(2),Q(5)):
        p,b = background(base,t)
        phi,d,hu,hformula,masses = stability_data(p,b)
        matrix_error = max(matrix_error,float(np.max(abs(hu-hformula)))/(1+float(np.max(abs(hu)))))
        stationary_error = max(stationary_error,float(np.max(abs(d['dU']))))
        fd = np.zeros((2,2))
        for j in range(2):
            step = np.eye(2)[j]*1e-5
            fd[:,j] = (previous.fields(phi+step,p)['dU']-previous.fields(phi-step,p)['dU'])/2e-5
        derivative_error = max(derivative_error,float(np.max(abs(fd-hu)))/(1+float(np.max(abs(hu)))))
        assert np.linalg.eigvalsh(d['G'])[0] > 0 and masses[0] > 0
        stable.append(dict(t=str(t),R=str(b['R']),canonical_radial_masses_squared=masses.tolist()))
    assert matrix_error < 1e-12 and stationary_error < 1e-10 and derivative_error < 1e-7
    checks.append('independent_potential_Hessian_differences_and_canonical_radial_mass_positivity')

    weak_certificates = 0
    bound_examples = []
    for eta,ng,r,alpha in itertools.product((Q(0),Q(1),Q(12)),(1,3),(Q(1),Q(7,3)),
                                           (Q(0),Q(1,2),Q(99,100))):
        p0 = parameters(Q(1,4),Q(1,4)*alpha,r,ng)
        upper = (3+r)*(2*ng/(1-eta/24)-1)/18
        for sign in (-1,0,1):
            t = 1+sign*eta/24
            p,b = background(p0,t)
            ratio = b['F']/b['h2']
            assert abs(4*b['R']/p['C0']) <= eta and ratio <= upper
            if alpha == 0 and sign == -1:
                assert ratio == upper  # sharp in declared parameter family
            weak_certificates += 1
        if ng == 3 and r == Q(7,3) and alpha == Q(1,2):
            bound_examples.append(dict(eta=str(eta),uniform_F_over_h2_upper_bound=str(upper)))
    checks.append('weak_curvature_hierarchy_bound_with_exact_saturating_examples')

    critical = []
    for power in (1,2,4,6):
        p,b = background(base,Q(1,10**power))
        ratio = b['F']/b['h2']
        normalized_R = 4*b['R']/p['C0']
        assert normalized_R == 24*(1-b['t'])
        assert ratio == p['q']*p['A']*(2*p['ng']/b['t']-1)/6
        critical.append(dict(t=str(b['t']),F_over_h2=str(ratio),R_over_mu2=str(normalized_R)))
    assert Q(critical[-1]['F_over_h2']) > 10**6 and Q(critical[-1]['R_over_mu2']) > 23
    checks.append('large_hierarchy_vacuum_only_escape_leaves_weak_curvature_regime')

    common_mass = []
    mass_error = 0.
    for kappa in (Q(1),Q(1,100),Q(1,10000),Q(1,1000000)):
        p,b = background(base,kappa,kappa)
        verify_exact_equations(p,b)
        assert b['R'] == b['V'] == 0
        assert b['V0'] == (kappa*p['C0'])**2*p['A']/4
        assert b['F']/b['h2'] == p['q']*p['A']*(2*p['ng']/kappa-1)/6
        # Opening the gravity/matter hierarchy leaves the old matter-only ratio unchanged.
        assert b['s2']/b['h2'] == p['r']*(p['q']-p['S'])/p['T']
        phi,d,hu,hformula,masses = stability_data(p,b)
        mass_error = max(mass_error,float(np.max(abs(hu-hformula))))
        assert masses[0] > 0 and d['F'] > 0
        common_mass.append(dict(kappa=str(kappa),F_over_h2=str(b['F']/b['h2']),
                                R=str(b['R']),V0=str(b['V0']),s2_over_h2=str(b['s2']/b['h2']),
                                canonical_radial_masses_squared=masses.tolist()))
    assert mass_error < 1e-12
    checks.append('one_extra_common_mass_matching_allows_flat_large_hierarchy_with_radial_stability')

    reference_examples = []
    for kappa in (Q(1),Q(1,100),Q(1,10000)):
        p,b = background(base,kappa,kappa)
        phi = np.sqrt([float(b['h2']),float(b['s2'])])
        d,grad,hess,je,jj,psi_zz = previous.origin_jets(phi,.125,p)
        expected = -4*(.125)**6*d['F']**2
        assert abs(np.linalg.det(jj)/expected-1) < 1e-12
        wave = np.einsum('mn,imn->i',previous.ETA,hess)
        wave += np.einsum('ijk,jk->i',d['connection'],grad@previous.ETA@grad.T)-d['Gi']@d['dU']
        assert np.max(abs(wave)) < 1e-12
        reference_examples.append(dict(kappa=str(kappa),field_values=phi.tolist(),
                                       nonuniform_reference_determinant_J=float(np.linalg.det(jj)),
                                       inherited_local_development_not_constant_vacuum=True))
    checks.append('inherited_nonuniform_reference_can_use_matched_field_values_but_constant_vacuum_has_no_reference_rank')

    # Scope failures: no hidden negative fields, no F=0 continuation, and no V-minimum shortcut.
    endpoint_data = []
    for t in (Q(0),Q(6)):
        h2 = t*base['C0']/base['q']
        F = base['C0']*base['A']*(6-t)/6
        endpoint_data.append(dict(t=str(t),h2=str(h2),F=str(F)))
    assert endpoint_data[0]['h2'] == '0' and endpoint_data[1]['F'] == '0'
    p,b = background(base,Q(1))
    bare = previous.bare_parameters()
    # Keep flat-potential minimum fields but restore bare V0: metric demands R=33/10,
    # whereas scalar equations with nonzero fields demand R=0.
    wrong_V = b['V']+bare['V0']-p['V0']
    metric_R = 4*wrong_V/b['F']
    assert metric_R == Q(33,10) and metric_R/6 != 0
    checks.append('endpoint_and_unmatched_flat_minimum_counterchecks_prevent_false_full_solution')

    dependencies = ('research_note_531.md','research_note_543.md','research_note_544.md',
                    'research_note_549.md','research_note_550.md','joint_reference_gravity_constraints.py',
                    'joint_reference_gravity_constraints_results.json','research_round_550_checks.json')
    return dict(round=551,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_constant_backgrounds=rows,rank_one_certificates=certificate_count,
        radial_stability_examples=stable,maximum_Hessian_relative_difference=derivative_error,
        maximum_analytic_Hessian_relative_residual=matrix_error,maximum_stationarity_residual=stationary_error,
        weak_curvature_certificates=weak_certificates,weak_curvature_bounds=bound_examples,
        near_critical_vacuum_only_sequence=critical,independent_common_mass_sequence=common_mass,
        inherited_reference_examples=reference_examples,endpoint_checks=endpoint_data,
        inconsistent_bare_flat_minimum_metric_R=str(metric_R),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in dependencies},
        scope=dict(classical_four_dimensional_two_derivative_model=True,
            constant_scalar_Einstein_branch_only_for_hierarchy_theorem=True,
            fixed_r_and_generation_number_for_uniform_bound=True,
            positive_F_and_two_nonzero_condensates=True,independent_V0_is_new_matching_input=True,
            additional_common_mass_matching_is_new_input=True,
            scalar_radial_linear_stability_only=True,constant_vacuum_is_not_full_rank_reference=True,
            nonuniform_reference_existence_inherited_from_549=True,
            F_over_h2_is_action_coefficient_ratio_not_measured_Newton_over_pole_mass=True,
            matter_internal_mass_ratios_not_released=True,quantum_or_radiative_stability_not_proved=True,
            higher_derivative_and_full_spectral_validity_not_proved=True,
            uniqueness_or_minimal_input_count_not_proved=True,completed_unification=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','maximum_Hessian_relative_difference')}))
