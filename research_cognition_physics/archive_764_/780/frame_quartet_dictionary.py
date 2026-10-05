"""780: local metric/frame cocycle, original spinor covariance, graded quartet.

Finite matrices and exact Grassmann polynomials are diagnostics. They do not
compute a continuum anomaly or justify a nonlinear quantum field redefinition.
The frozen 777 algebra supplies signs, not the action tested here.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'bv777', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
add, mul, scale, diff = bv.add, bv.mul, bv.scale, bv.diff
ETA = np.diag([-1., 1., 1., 1.])
I4 = np.eye(4)


def exp_derivative(x, dx=None):
    """Matrix exponential and Frechet derivative via convergent power series."""
    dx = np.zeros_like(x) if dx is None else dx
    term, dterm = np.eye(len(x), dtype=x.dtype), np.zeros_like(x)
    value, dvalue = term.copy(), dterm.copy()
    for n in range(1, 100):
        term, dterm = term@x/n, (dterm@x+term@dx)/n
        value, dvalue = value+term, dvalue+dterm
        if max(np.linalg.norm(term), np.linalg.norm(dterm)) < 1e-17:
            return value, dvalue
    raise AssertionError('matrix series did not converge')


def sqrt_near_one(a):
    x = a-I4
    assert np.linalg.norm(x, 2) < .8
    value, power, coefficient = I4.copy(), I4.copy(), 1.
    for n in range(1, 200):
        power = power@x
        coefficient *= (1.5-n)/n
        term = coefficient*power
        value += term
        if np.linalg.norm(term) < 1e-17:
            return value
    raise AssertionError('square-root series did not converge')


def section(g):
    return sqrt_near_one(ETA@g)


def compensator(j, g):
    return section(g)@j@np.linalg.inv(section(j.T@g@j))


def lorentz_generator(rng, size=.025):
    a = rng.normal(size=(4, 4))*size
    return ETA@(a-a.T)


def metric_frame_cocycle():
    rng = np.random.default_rng(780)
    errors = dict(section_metric=0., lorentz_factor=0., reconstruction=0.,
                  finite_cocycle=0., sylvester_derivative=0.)
    frozen_errors = []
    for _ in range(24):
        h = rng.normal(size=(4, 4))*.012
        g = ETA+h+h.T
        e = section(g)
        lam = exp_derivative(lorentz_generator(rng))[0]
        full = lam@e
        recovered = full@np.linalg.inv(e)
        j1, j2 = I4+rng.normal(size=(4, 4))*.015, I4+rng.normal(size=(4, 4))*.015
        g1 = j1.T@g@j1
        l1, l2 = compensator(j1, g), compensator(j2, g1)
        values = dict(section_metric=np.linalg.norm(e.T@ETA@e-g),
                      lorentz_factor=np.linalg.norm(recovered.T@ETA@recovered-ETA),
                      reconstruction=np.linalg.norm(recovered@e-full),
                      finite_cocycle=np.linalg.norm(compensator(j1@j2, g)-l1@l2))
        # dB solves B*dB+dB*B=eta*dg. Compare to a centered finite difference.
        dg = rng.normal(size=(4, 4))*.1
        dg = dg+dg.T
        operator = np.kron(I4, e)+np.kron(e.T, I4)
        db = np.linalg.solve(operator, (ETA@dg).reshape(-1, order='F')).reshape((4, 4), order='F')
        eps = 1e-5
        numeric = (section(g+eps*dg)-section(g-eps*dg))/(2*eps)
        values['sylvester_derivative'] = np.linalg.norm(db-numeric)
        for key, value in values.items():
            errors[key] = max(errors[key], float(value))
        frozen_errors.append(float(np.linalg.norm(compensator(j1@j2, g)-l1@compensator(j2, g))))
    assert max(errors.values()) < 2e-9, errors
    a, b = np.zeros((4, 4)), np.zeros((4, 4))
    a[1, 1], a[2, 2] = 1, -1
    b[1, 2] = b[2, 1] = 1
    projection = lambda m: (m-ETA@m.T@ETA)/2
    commutator = a@b-b@a
    assert not np.any(projection(a)) and not np.any(projection(b))
    assert np.array_equal(projection(commutator), commutator)
    assert min(frozen_errors) > 1e-6
    return dict(samples=24, maximum_residuals=errors,
                freeze_metric_in_finite_cocycle_minimum_defect=min(frozen_errors),
                frozen_infinitesimal_cocycle_defect=commutator.tolist(),
                result_scope='Near-background local section; no global square-root or topology assertion.')


def gamma_matrices():
    zero, eye = np.zeros((2, 2)), np.eye(2)
    pauli = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    return [1j*np.block([[eye, zero], [zero, -eye]])]+[
        1j*np.block([[zero, s], [-s, zero]]) for s in pauli]


def spinor_connection():
    rng = np.random.default_rng(781)
    gamma = gamma_matrices()
    gamma5 = 1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]
    charge = 1j*gamma[2]@gamma[0]
    def rho(lam):
        lower = ETA@lam
        return sum((lower[a, b]*gamma[a]@gamma[b]/4 for a in range(4) for b in range(4)), np.zeros((4, 4), complex))
    errors = dict(clifford=0., spin_lift=0., chirality=0., majorana=0.,
                  connection_lift=0., dirac_covariance=0.)
    for a in range(4):
        for b in range(4):
            errors['clifford'] = max(errors['clifford'], float(np.linalg.norm(gamma[a]@gamma[b]+gamma[b]@gamma[a]-2*ETA[a,b]*I4)))
    omitted = []
    for _ in range(16):
        h = rng.normal(size=(4, 4))*.01
        e = section(ETA+h+h.T)
        generator = lorentz_generator(rng, .07)
        lam, _ = exp_derivative(generator)
        spin, _ = exp_derivative(rho(generator))
        inv_spin = np.linalg.inv(spin)
        spin_err = max(np.linalg.norm(inv_spin@gamma[a]@spin-sum(lam[a,b]*gamma[b] for b in range(4))) for a in range(4))
        errors['spin_lift'] = max(errors['spin_lift'], float(spin_err))
        errors['chirality'] = max(errors['chirality'], float(np.linalg.norm(spin@gamma5-gamma5@spin)))
        errors['majorana'] = max(errors['majorana'], float(np.linalg.norm(spin.T@charge@spin-charge)))
        psi = rng.normal(size=4)+1j*rng.normal(size=4)
        dpsi = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        old_cov, new_cov, wrong_cov = [], [], []
        for mu in range(4):
            velocity = lorentz_generator(rng, .11)
            _, dlam = exp_derivative(generator, velocity)
            _, dspin = exp_derivative(rho(generator), rho(velocity))
            omega = lorentz_generator(rng, .09)
            big_omega = rho(omega)
            transformed = spin@big_omega@inv_spin-dspin@inv_spin
            lorentz_transformed = lam@omega@np.linalg.inv(lam)-dlam@np.linalg.inv(lam)
            errors['connection_lift'] = max(errors['connection_lift'], float(np.linalg.norm(transformed-rho(lorentz_transformed))))
            old_cov.append(dpsi[mu]+big_omega@psi)
            new_cov.append(dspin@psi+spin@dpsi[mu]+transformed@spin@psi)
            wrong_cov.append(dspin@psi+spin@dpsi[mu]+spin@big_omega@psi)
        inv_e, inv_new = np.linalg.inv(e), np.linalg.inv(lam@e)
        contract = lambda inverse, cov: sum(inverse[mu,a]*(gamma[a]@cov[mu]) for mu in range(4) for a in range(4))
        old_dirac = contract(inv_e, old_cov)
        errors['dirac_covariance'] = max(errors['dirac_covariance'], float(np.linalg.norm(contract(inv_new,new_cov)-spin@old_dirac)))
        omitted.append(float(np.linalg.norm(contract(inv_new,wrong_cov)-spin@old_dirac)))
    assert max(errors.values()) < 1e-11, errors
    assert min(omitted) > 1e-3
    return dict(samples=16, maximum_residuals=errors,
                omit_inhomogeneous_connection_minimum_defect=min(omitted),
                same_chiral_and_majorana_bilinears=True,
                independent_torsion_or_new_fermions_added=False)


def quartet_bv():
    q,r,c,h,p,a = [bv.var(x) for x in ('q','r','c','h','p','a')]
    action = add(mul(p,c), scale(mul(a,r),-1), mul(r,q), scale(mul(h,c),-1))
    s = lambda f: bv.bracket(action, f)
    contraction = lambda f: add(diff(diff(f,'q'),'r'), scale(diff(diff(f,'c'),'h'),-1))
    pairs = (('q','p'),('r','t'),('c','z'),('h','a'))
    laplacian = lambda f: add(*(scale(diff(diff(f,anti),field),-1 if field in bv.ODD else 1) for field,anti in pairs))
    assert not bv.bracket(action,action)
    expected = dict(q=c,r={},c={},h=r,p=r,t=add(q,scale(a,-1)),z=add(p,scale(h,-1)),a=c)
    for name,result in expected.items():
        assert s(bv.var(name)) == result, (name,bv.display(s(bv.var(name))))
    count = 0
    names = tuple(expected)
    for degree in range(5):
        for term in combinations_with_replacement(names, degree):
            f = bv.prod(*(bv.var(name) for name in term))
            if not f:
                continue
            assert not s(s(f))
            assert not laplacian(laplacian(f))
            assert not add(s(laplacian(f)),laplacian(s(f)))
            assert not add(s(contraction(f)),scale(contraction(s(f)),-1),laplacian(f))
            count += 1
    return dict(supermonomials=count, classical_master_equation_exact=True,
                brst_on_generators={name:bv.display(poly) for name,poly in expected.items()},
                graded_commutator='[s,C]=-Delta',
                bosonic_hessian=[[0,1],[1,0]],
                propagating_degrees_of_freedom_added=0)


def smul(left,right,order):
    return [add(*(mul(left[i],right[k-i]) for i in range(k+1))) for k in range(order+1)]


def compose(poly, mapping, order):
    out = [{} for _ in range(order+1)]
    for key,value in poly.items():
        term = [scale(bv.ONE,value)]+[{} for _ in range(order)]
        for name,power in zip(bv.NAMES,key):
            factor = mapping.get(name, [bv.var(name)]+[{} for _ in range(order)])
            for _ in range(power):
                term = smul(term,factor,order)
        out = [add(x,y) for x,y in zip(out,term)]
    return out


def graded_stationary(vertex,order):
    fields = ('q','r','c','h')
    u = {name:[{} for _ in range(order+1)] for name in fields}
    # Left derivatives: S0=q*r+c*h, so S0_c=h and S0_h=-c.
    inverse = dict(q=('r',-1),r=('q',-1),c=('h',1),h=('c',-1))
    for _ in range(order):
        mapping = {name:[bv.var(name)]+value[1:] for name,value in u.items()}
        fresh = {}
        for name,(source,sign) in inverse.items():
            gradient = compose(diff(vertex,source),mapping,order)
            fresh[name] = [{}]+[scale(f,sign) for f in gradient[:-1]]
        u = fresh
    mapping = {name:[bv.var(name)]+value[1:] for name,value in u.items()}
    vr = [{}]+compose(vertex,mapping,order)[:-1]
    stationary_cost = [add(a,b) for a,b in zip(smul(u['q'],u['r'],order),smul(u['c'],u['h'],order))]
    result = [add(a,b) for a,b in zip(vr,stationary_cost)]
    return u,result,mapping,vr


def graded_contact_envelope():
    q,r,c,h,p,a = [bv.var(x) for x in ('q','r','c','h','p','a')]
    vertex = add(scale(bv.prod(q,q,q),Q(1,6)),bv.prod(q,c,h),mul(a,q),
                 bv.prod(q,p,c),scale(bv.prod(q,r,r),Q(1,2)))
    order = 4
    u,reduced,mapping,without_cost = graded_stationary(vertex,order)
    checks = 0
    for name in ('q','r','c','h','p','a'):
        expected = [{}]+compose(diff(vertex,name),mapping,order)[:-1]
        for degree in range(order+1):
            assert diff(reduced[degree],name) == expected[degree], (name,degree)
            checks += 1
    stationary = dict(q=('r',1),r=('q',1),c=('h',1),h=('c',-1))
    for name,(partner,sign) in stationary.items():
        grad = [{}]+compose(diff(vertex,name),mapping,order)[:-1]
        for degree in range(order+1):
            assert not add(scale(u[partner][degree],sign),grad[degree])
            checks += 1
    expected = [{}]+compose(diff(vertex,'a'),mapping,order)[:-1]
    defect = add(diff(without_cost[2],'a'),scale(expected[2],-1))
    assert defect
    # On a pure spectator interaction there is no quartet tree to contract.
    pure = add(bv.prod(a,a,a),bv.prod(p,bv.var('t')))
    _,same,_,_ = graded_stationary(pure,order)
    assert same[1] == pure and not any(same[2:])
    return dict(order=order, exact_envelope_and_stationarity_checks=checks,
                odd_and_antifield_mixed_vertices=True,
                stationary_cost_omission_defect=bv.display(defect),
                second_order_contact=bv.display(reduced[2]),
                physical_spectator_action_unchanged=True,
                continuum_prescription_proven_by_finite_test=False)


def run():
    return dict(round=780,groups={
        'local_metric_frame_cocycle':metric_frame_cocycle(),
        'same_spinor_and_connection':spinor_connection(),
        'adapted_quartet_bv':quartet_bv(),
        'graded_tree_contact':graded_contact_envelope()},
        all_checks_passed=True, original_full_N1_N2_proven=False,
        nonlinear_original_chart_quantum_transport_proven=False,
        new_physical_structure_derived_from_cognition=False,
        original_loop_anomaly_computed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = run()
    destination = HERE/'frame_quartet_dictionary_results.json'
    if args.write:
        destination.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if args.check:
        assert result == json.loads(destination.read_text(encoding='utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
