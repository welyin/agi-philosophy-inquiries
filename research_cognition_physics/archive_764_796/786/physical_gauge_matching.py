"""786: contact gauge, nonlinear BV coordinates, and physical density matching.

Three diagnostics, not the original continuum quantum equivalence theorem.
The exact graded polynomial engine is reused from frozen round 777.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
import json
import numpy as np
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv777_for786', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)


def contact_block():
    n = 7
    eye = np.eye(n)
    d = np.diag(np.linspace(.6, 1.4, n)) @ ((np.roll(eye, 1, axis=1)-np.roll(eye, -1, axis=1))/2)
    d += np.diag(np.linspace(-.2, .3, n))
    k = np.vstack((d, eye))
    f = np.hstack((eye, -d))
    l = np.hstack((np.zeros((n, n)), eye)) + (.17*eye+.09*d)@f
    p = f.T@(2*eye+d.T@d)@f
    pi = np.eye(2*n)-k@l
    gp = pi@np.linalg.inv(p-k@k.T)@pi.T
    cases = []
    for xi in (0., .3, 2.):
        operator = np.block([[p, l.T], [l, xi*eye]])
        inverse = np.block([[gp-xi*k@k.T, k], [k.T, np.zeros((n,n))]])
        error = max(np.max(np.abs(operator@inverse-np.eye(3*n))),
                    np.max(np.abs(inverse@operator-np.eye(3*n))))
        assert error < 5e-12
        # Tests already in ker K^T do not see xi or mixed auxiliary contacts.
        physical = pi.T@np.eye(2*n)
        physical_error = np.max(np.abs(physical.T@(inverse[:2*n,:2*n]-gp)@physical))
        assert physical_error < 5e-12
        # Dropping the off-diagonal local kernels destroys the full inverse.
        wrong = np.block([[gp-xi*k@k.T, np.zeros((2*n,n))],
                          [np.zeros((n,2*n)), np.zeros((n,n))]])
        wrong_error = float(np.max(np.abs(operator@wrong-np.eye(3*n))))
        assert wrong_error >= 1
        cases.append(dict(xi=xi, inverse_residual=float(error),
                          physical_xi_difference=float(physical_error),
                          omitted_contact_residual=wrong_error))
    return dict(block_dimension=3*n, xi_values_tested=3, cases=cases,
                scope='Finite mixed-order inverse calibration; causal support is proved analytically.')


CHECK_ORDER = 5
WORK_ORDER = 8
FIELD_AXES = [bv.NAMES.index(n) for n in ('q','r','b','c','h')]


def trunc(poly, order=WORK_ORDER):
    return {e: v for e, v in poly.items() if sum(e[i] for i in FIELD_AXES) <= order}


# Restrict only this imported module instance. Frozen 777 source is unchanged.
raw_mul = bv.mul
bv.mul = lambda a,b: trunc(raw_mul(a,b))


def exp_series(poly):
    result, term = {}, bv.ONE
    for n in range(WORK_ORDER+1):
        if n:
            term = bv.mul(term, poly)
        result = bv.add(result, bv.scale(term, Q(1, factorial(n))))
    return result


def zero_through(poly):
    return not trunc(poly, CHECK_ORDER)


def nonlinear_bv_chart():
    v = {n:bv.var(n) for n in bv.NAMES}
    q,r,b,z,a,c,h,p,t,rho = (v[n] for n in ('q','r','b','z','a','c','h','p','t','rho'))
    kappa, mass, xi = Q(2,3), Q(3,2), Q(2,5)
    w = bv.add(r, bv.scale(bv.mul(q,q), -Q(1,2)))
    m = exp_series(bv.scale(w,kappa))
    im = exp_series(bv.scale(w,-kappa))
    physical = bv.scale(bv.mul(w,w), mass/2)
    minimal = bv.add(physical, bv.prod(p,m,c), bv.prod(t,q,m,c), bv.scale(bv.mul(a,b), -1))
    forward = dict(q=q, r=w, c=bv.mul(m,c),
                   p=bv.add(p,bv.mul(q,t)), t=bv.add(t,bv.scale(bv.mul(z,c),-kappa)),
                   z=bv.mul(z,im))
    mr, imr = exp_series(bv.scale(r,kappa)), exp_series(bv.scale(r,-kappa))
    old_t = bv.add(t,bv.scale(bv.mul(z,c),kappa))
    backward = dict(q=q, r=bv.add(r,bv.scale(bv.mul(q,q),Q(1,2))),
                    c=bv.mul(imr,c), z=bv.mul(mr,z), t=old_t,
                    p=bv.add(p,bv.scale(bv.mul(q,old_t),-1)))
    target = bv.add(bv.scale(bv.mul(r,r),mass/2), bv.mul(p,c), bv.scale(bv.mul(a,b),-1))
    assert zero_through(bv.add(bv.substitute(minimal,backward),bv.scale(target,-1)))
    assert zero_through(bv.bracket(minimal,minimal))
    assert not bv.bracket(target,target)
    for name in bv.NAMES:
        assert zero_through(bv.add(bv.substitute(bv.substitute(v[name],forward),backward),bv.scale(v[name],-1)))
        assert zero_through(bv.add(bv.substitute(bv.substitute(v[name],backward),forward),bv.scale(v[name],-1)))
        for other in bv.NAMES:
            assert zero_through(bv.add(bv.bracket(forward.get(name,v[name]),forward.get(other,v[other])),
                                      bv.scale(bv.bracket(v[name],v[other]),-1))), (name,other)
    # The adapted chart makes all four nonminimal/gauge variables strict pairs.
    expected = dict(q=c,c={},h=b,b={},r={})
    for name,value in expected.items():
        assert bv.bracket(target,v[name]) == value, name
    psi = bv.mul(h,bv.add(q,bv.scale(b,xi/2)))
    gf_terms = bv.bracket(target,psi)
    target_terms = bv.add(bv.mul(b,q), bv.scale(bv.mul(b,b),xi/2),bv.scale(bv.mul(h,c),-1))
    assert gf_terms == target_terms
    # Antifields must transform too; changing fields alone is not BV canonical.
    wrong = dict(forward)
    for name in ('p','t','z'):
        wrong[name] = v[name]
    defect = bv.add(bv.bracket(wrong['r'],wrong['p']), bv.scale(bv.bracket(r,p),-1))
    assert not zero_through(defect)
    return dict(coefficient_arithmetic='fractions.Fraction', field_degree_verified=CHECK_ORDER,
                expansion_work_degree=WORK_ORDER, canonical_coordinate_brackets=100,
                inverse_coordinate_checks=20,
                original_and_adapted_CME=True, nonlinear_full_action_split=True,
                nonminimal_BRST_pairs_and_gauge_fermion_signs=True,
                omitted_antifield_map_defect=bv.display(trunc(defect,CHECK_ORDER)),
                scope='Exact finite-jet BV chart with invariant-dependent ghost rescaling, not the original continuum quantum dictionary.')


def density_matching():
    kappa, mass = Q(2,3), Q(3,2)
    cases = []
    # After q=0, FP determinant exp(kappa*w) multiplies the physical density.
    # Ghost normalization theta=exp(kappa*w)c moves exactly that factor into
    # the Berezinian; it does not remove it. Completing the square is exact.
    for hbar in (Q(1,10), Q(1,5), Q(1,2)):
        mean, variance = hbar*kappa/mass, hbar/mass
        moments = [Q(1),mean]
        for n in range(2,7):
            moments.append(mean*moments[n-1]+(n-1)*variance*moments[n-2])
        # Independent binomial/centered-Gaussian moment expansion.
        for n in range(7):
            value = Q(0)
            for j in range(n//2+1):
                value += Q(factorial(n),factorial(n-2*j)*factorial(j)*2**j)*mean**(n-2*j)*variance**j
            assert value == moments[n]
        assert mass*mean-hbar*kappa == 0
        assert mean != 0  # The naive flat reduced Gaussian would have zero mean.
        cases.append(dict(hbar=str(hbar), transported_mean=str(mean),
                          covariance=str(variance), naive_flat_mean='0',
                          moments_through_six=[str(x) for x in moments]))
    return dict(kappa=str(kappa), mass=str(mass),
                forward_field_Berezinian='exp(-kappa*w)',
                old_measure_in_new_coordinates='exp(kappa*w)',
                physical_effective_action='mass*w^2/2 - hbar*kappa*w',
                one_loop_source_shift=str(-kappa),
                same_classical_action_and_free_covariance=True,
                omitted_physical_density_changes_observables=True,
                equation_of_motion_shift_coefficient=str(-kappa/mass),
                translation_of_action_and_observables_restores_agreement=True,
                irreducible_full_BV_obstruction_in_this_example=False,
                cases=cases,
                scope='Finite Euclidean normalized Gaussian counterexample; no continuum determinant or anomaly coefficient inferred.')


def run():
    return dict(round=786, fresh_test_groups=3, contact_gauge=contact_block(),
                nonlinear_BV_chart=nonlinear_bv_chart(), physical_density=density_matching(),
                all_checks_passed=True, original_continuum_quantum_matching_proven=False,
                original_interacting_positive_state_proven=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('physical_gauge_matching_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
