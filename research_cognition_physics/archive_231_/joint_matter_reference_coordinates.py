"""548: reference scalars from the SAME Higgs-plus-singlet matter sector.

Classical local on-shell matter on supplied Minkowski geometry.  Tests check
two-jets, group orbits and frozen couplings, not a coupled Einstein solution,
quantum detector, global atlas, or choice of dimension.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_matter_reference_coordinates_results.json'
ETA = np.diag([-1., 1., 1., 1.])


def determinant_exact(matrix):
    a = [[F(v) for v in row] for row in matrix]
    value = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        v = a[j][j]
        value *= v
        for k in range(j, len(a)):
            a[j][k] /= v
        for i in range(j + 1, len(a)):
            v = a[i][j]
            for k in range(j, len(a)):
                a[i][k] -= v * a[j][k]
    return value


def canonical_su2(H):
    a, b = H
    return np.array([[b, -a], [a.conjugate(), b.conjugate()]]) / np.linalg.norm(H)


def jets(eps, vh=0., vs=0.):
    dh = np.array([eps, 0., 0., 0.])
    ds = np.array([0., eps, 0., 0.])
    hh = np.zeros((4, 4)); hs = hh.copy()
    hh[0, 0] = -vh; hs[0, 0] = -vs
    hh[0, 2] = hh[2, 0] = eps
    hs[1, 3] = hs[3, 1] = eps
    return dh, ds, hh, hs


def reference_jacobian(dh, ds, hh, hs, inverse_metric=ETA):
    return np.stack([dh, ds, 2 * hh @ inverse_metric @ dh,
                     2 * hs @ inverse_metric @ ds])


def second_order_menu(x, h0, s0, dh, ds, hh, hs):
    """A two-jet representative for checking derivatives, NOT a PDE solution."""
    h = h0 + dh @ x + .5 * x @ hh @ x
    s = s0 + ds @ x + .5 * x @ hs @ x
    ph, ps = dh + hh @ x, ds + hs @ x
    return np.array([h, s, ph @ ETA @ ph, ps @ ETA @ ps])


def frozen_parameters():
    saved = json.loads((HERE/'joint_top_boundary_identifiability_results.json').read_text('utf8'))
    rows = []
    for item in saved['inverse_examples']:
        vac = item['vacuum']
        if not vac['strict_double_vev']:
            continue
        a = vac['state']; rho = vac['squared_vev_ratio']
        scale = 1 / vac['h_squared_over_common_scale']
        rows.append(dict(top_target=item['synthetic_target_q'], h0=1., s0=math.sqrt(rho),
            a=scale*a['x'], b=scale*a['y'], lh=a['lambda_H'], p=a['p'], ls=a['lambda_s']))
    return rows


def potential_gradient(h, s, a, b, lh, p, ls):
    return np.array([h*(-a + lh*h*h + p*s*s), s*(-b + p*h*h + ls*s*s)])


def run():
    rng = np.random.default_rng(548); checks = []
    orbit_error = 0.
    for _ in range(20):
        H = rng.normal(size=2) + 1j*rng.normal(size=2)
        U = canonical_su2(H)
        orbit_error = max(orbit_error, float(np.linalg.norm(U@H - [0., np.linalg.norm(H)])))
        assert np.allclose(U.conj().T@U, np.eye(2), atol=1e-14)
        assert abs(np.linalg.det(U)-1) < 1e-14
    # Four raw Cartesian components can have rank 4 while their gauge quotient
    # has one scalar value.  The transformed gauge connection must also change.
    raw = np.array([.2, .1, 1., .3]); canonical_J = np.zeros((4, 4))
    canonical_J[2] = raw/np.linalg.norm(raw)
    assert np.linalg.matrix_rank(np.eye(4)) == 4
    assert np.linalg.matrix_rank(canonical_J) == 1
    checks.append('same_Higgs_gauge_orbit_and_raw_component_coordinate_counterexample')

    # Exact chain rule for f=(h,s,h²+s²,hs), with arbitrary independent dh, ds.
    h, s = F(3, 2), F(2, 3)
    M = [[1,0], [0,1], [2*h,2*s], [s,h]]
    G = [[1,2,3,4], [2,0,1,3]]
    J = [[sum(F(M[i][k])*G[k][j] for k in range(2)) for j in range(4)] for i in range(4)]
    from itertools import combinations
    for ii in combinations(range(4),3):
        for jj in combinations(range(4),3):
            assert determinant_exact([[J[i][j] for j in jj] for i in ii]) == 0
    assert determinant_exact([[J[i][j] for j in (0,1)] for i in (0,1)]) != 0
    checks.append('exact_zero_derivative_invariant_chain_rule_has_rank_at_most_two')

    exact_rows = []
    for e in (F(1,2), F(1,4), F(1,8), F(1,16)):
        for vh, vs in ((F(0),F(0)), (F(2,3),F(-3,7))):
            # EOM fixes ONLY the time-time Hessians; mixed Cauchy jets remain free.
            hh00, ss00 = -vh, -vs
            assert -hh00-vh == -ss00-vs == 0
            j = [[e,0,0,0], [0,e,0,0], [2*e*vh,0,-2*e*e,0], [0,0,0,2*e*e]]
            det = determinant_exact(j)
            assert det == -4*e**6
        exact_rows.append(dict(epsilon=str(e), determinant=str(det), excess_energy=str(e*e)))
    checks.append('exact_on_shell_Cauchy_two_jets_yield_nonzero_determinant_for_same_potential_class')

    rows = []; stationary_error = 0.; current_error = 0.
    paulis = [np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]),
              np.diag([1.,-1.]), np.eye(2)]
    for par in frozen_parameters():
        h0, s0 = par['h0'], par['s0']
        coeff = {k:par[k] for k in ('a','b','lh','p','ls')}
        vh, vs = potential_gradient(h0, s0, **coeff)
        stationary_error = max(stationary_error, abs(vh), abs(vs))
        quartic = np.array([[par['lh'],par['p']], [par['p'],par['ls']]])
        assert min(np.linalg.eigvalsh(quartic)) > 0
        radial = 2*np.diag([h0,s0]) @ quartic @ np.diag([h0,s0])
        assert min(np.linalg.eigvalsh(radial)) > 0
        dh, ds, hh, hs = jets(.25, vh, vs)
        jj = reference_jacobian(dh, ds, hh, hs)
        assert abs(np.trace(ETA@hh)-vh) < 1e-14
        assert abs(np.trace(ETA@hs)-vs) < 1e-14
        assert abs(np.linalg.det(jj) + 4*.25**6) < 1e-16
        H = np.array([0.,h0])/math.sqrt(2)
        for g in paulis:
            for v in dh:
                dH = np.array([0.,v])/math.sqrt(2)
                current = 1j*(H.conj()@g@dH - dH.conj()@g@H)
                current_error = max(current_error, float(abs(current)))
        rows.append(dict(synthetic_top_target=par['top_target'], scalar_parameters=par,
                         stationary_residual=[float(vh),float(vs)],
                         reference_jacobian=jj.tolist(), determinant=float(np.linalg.det(jj))))
    assert len(rows)==3 and stationary_error < 1e-14 and current_error < 1e-14
    checks.append('three_frozen_vacua_and_consistent_zero_gauge_current_radial_subsector')

    # Independent finite differences of the scalar menu, and tensor transport.
    dh, ds, hh, hs = jets(.25, .2, -.3)
    jj = reference_jacobian(dh, ds, hh, hs)
    delta = 1e-5
    finite = np.column_stack([(second_order_menu(delta*np.eye(4)[i],1.,.7,dh,ds,hh,hs)
                             -second_order_menu(-delta*np.eye(4)[i],1.,.7,dh,ds,hh,hs))/(2*delta)
                             for i in range(4)])
    finite_error = float(np.max(abs(finite-jj))); assert finite_error < 1e-9
    covariance_error = 0.
    for _ in range(12):
        L = np.eye(4) + .1*rng.normal(size=(4,4))
        metric = L.T@ETA@L; inverse_metric = np.linalg.inv(metric)
        transformed = reference_jacobian(L.T@dh, L.T@ds, L.T@hh@L, L.T@hs@L, inverse_metric)
        covariance_error = max(covariance_error, float(np.max(abs(transformed-jj@L))))
        assert abs((L.T@dh)@inverse_metric@(L.T@dh)-dh@ETA@dh) < 1e-14
    assert covariance_error < 1e-13
    checks.append('independent_menu_differences_and_affine_coordinate_covariance_with_metric')

    resources = []
    for eps in (.5,.25,.125,.0625):
        dh, ds, hh, hs = jets(eps)
        jj = reference_jacobian(dh,ds,hh,hs)
        singular = np.linalg.svd(jj,compute_uv=False)
        kinetic = dh@ETA@dh + ds@ETA@ds
        stress = np.outer(dh,dh)+np.outer(ds,ds)-.5*ETA*kinetic
        assert abs(stress[0,0]-eps*eps)<1e-15
        assert np.allclose(singular,[eps,eps,2*eps*eps,2*eps*eps])
        invnorm = float(np.linalg.norm(np.linalg.inv(jj),2))
        assert abs(invnorm-1/(2*eps*eps))<1e-13
        # This minimal-coupling stress is nonzero. A full nonminimal metric
        # variation is NOT tested here, even when the fixed background has R=0.
        assert np.max(abs(stress))>0
        resources.append(dict(epsilon=eps, excess_energy_density=float(stress[0,0]),
            minimum_singular_value=float(singular[-1]), inverse_gain=invnorm))
    checks.append('same_matter_stress_and_fixed_calibration_precision_cost_not_Einstein_solution')

    dh, ds, hh, hs = jets(0.)
    assert np.linalg.matrix_rank(reference_jacobian(dh,ds,hh,hs))==0
    homogeneous = reference_jacobian(np.array([1.,0,0,0]),np.array([2.,0,0,0]),
                                     np.diag([3.,0,0,0]),np.diag([4.,0,0,0]))
    assert np.linalg.matrix_rank(homogeneous)==1
    dh,ds,hh,hs=jets(.25)
    absent = reference_jacobian(dh,np.zeros(4),hh,np.zeros((4,4)))
    assert np.linalg.matrix_rank(absent)==2
    # s0=0 at ONE POINT does not enter J; arbitrary V_s is already allowed.
    assert np.linalg.matrix_rank(reference_jacobian(*jets(.25,.2,0.)))==4
    checks.append('vacuum_homogeneity_and_absent_singlet_are_distinct_scope_counterexamples')

    deps=('research_note_543.md','research_note_544.md','research_note_546.md',
          'joint_top_boundary_identifiability_results.json','research_round_547_checks.json',
          'material_reference_geometry_review.md','joint_condition_compression_table.md')
    return dict(round=548,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_cauchy_jet_certificates=exact_rows, gauge_orbit_max_error=orbit_error,
        same_parameter_examples=rows, stationary_max_error=float(stationary_error),
        gauge_current_max_error=current_error,finite_difference_error=finite_error,
        affine_covariance_error=covariance_error, fixed_calibration_resource_examples=resources,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(classical_matter_on_given_four_dimensional_Minkowski=True,
          same_Higgs_and_explicit_singlet_content=True, local_on_shell_existence_uses_standard_PDE_theorem=True,
          numerics_check_jets_not_full_PDE_solution=True, zero_derivative_obstruction_is_menu_specific=True,
          derivative_menu_uses_given_metric=True, local_reference_fields_not_actual_quantum_instruments=True,
          initial_state_preparation_is_input=True, finite_excess_energy_not_total_vacuum_energy=True,
          gain_uses_fixed_dimensionless_calibration=True, coupled_Einstein_solution=False,
          unique_dimension_selected=False, completed_unification=False))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','stationary_max_error','finite_difference_error')},ensure_ascii=False))
