"""708: one original block fluctuation controls reading cost and mass contrasts.

The full-H double commutator is analytic on the compact smooth Gauss core.
Numerics use the inherited H5 metric and all 32 original CAR modes per node;
they do not replace the quantum model by the sampled point configurations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_geodesic_spatial_block as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_block_fluctuation_contract_results.json'
R = old.R
HBAR = .7
WEIGHT = .8


def embedding_from_spatial(y):
    return np.r_[np.sqrt(R * R + y @ y), y]


def points_from_spatial(ys):
    return np.array([old.klein(embedding_from_spatial(y)) for y in ys])


def singlets(points):
    return np.array([old.hyper(p)[5] for p in points])


def kraus_derivative_sum(z):
    return np.cos(z)**2 / (4 * (4 - np.sin(z)**2))


def real_uniform_basis(n):
    # Deterministic Householder reflection; first row is the uniform vector.
    v = np.eye(n)[0] - np.ones(n) / np.sqrt(n)
    return np.eye(n) - 2 * np.outer(v, v) / (v @ v)


def full_mass_transform(points):
    n = len(points)
    h = np.zeros((32*n, 32*n), complex)
    d = np.zeros_like(h)
    for i, p in enumerate(points):
        hi, di = old.matter.mass_matrices(p)
        sl = slice(32*i, 32*i+32)
        h[sl, sl] = hi
        d[sl, sl] = di
    U = np.kron(real_uniform_basis(n), np.eye(32))
    return U @ h @ U.T, U @ d @ U.T


def mass_coefficients():
    out = []
    for a in range(5):
        y = np.eye(5)[a]
        out.append(old.matter.mass_matrices(old.klein(embedding_from_spatial(y))))
    return out


def geometry_and_budget():
    rng = np.random.default_rng(70831)
    rows = []
    max_gradient_error = max_budget_error = 0.
    for n in (3, 4, 7):
        points = .29 * rng.normal(size=(n, 5))
        z = singlets(points)
        Z, P2 = float(z.sum()), float(z @ z)
        finite_gamma = sum(old.grad(lambda x: old.hyper(x)[5], p)
                           @ np.linalg.inv(old.original.metric(p))
                           @ old.grad(lambda x: old.hyper(x)[5], p) for p in points)
        analytic_gamma = n + P2/R**2
        ge = abs(finite_gamma-analytic_gamma)
        # Differentiate the actual Kraus multipliers in each original Klein chart.
        actual_cost = 0.
        for i, p in enumerate(points):
            other = Z-z[i]
            inv = np.linalg.inv(old.original.metric(p))
            for sign in (-1, 1):
                grad = old.grad(lambda x: np.sqrt(.5+sign*.25*np.sin(other+old.hyper(x)[5])), p)
                actual_cost += HBAR**2/(2*WEIGHT) * (grad @ inv @ grad)
        predicted = HBAR**2/(2*WEIGHT)*kraus_derivative_sum(Z)*analytic_gamma
        be = abs(actual_cost-predicted)
        max_gradient_error = max(max_gradient_error, ge)
        max_budget_error = max(max_budget_error, be)
        rows.append(dict(nodes=n, Z=Z, P2=P2, gamma=analytic_gamma,
                         actual_kraus_cost=float(actual_cost), predicted_cost=predicted,
                         gradient_error=float(ge), budget_error=float(be)))
    assert max_gradient_error < 2e-8 and max_budget_error < 2e-9
    return dict(rows=rows, maximum_gradient_error=max_gradient_error,
                maximum_budget_error=max_budget_error)


def mass_and_common_witness():
    coeff = mass_coefficients()
    Ds2 = float(np.linalg.norm(coeff[4][1])**2)
    gram = np.array([[np.real(np.vdot(h, g)+np.vdot(d, e))
                      for g, e in coeff] for h, d in coeff])
    rng = np.random.default_rng(70832)
    points = .34*rng.normal(size=(4, 5))
    ys = np.array([old.hyper(p)[1:] for p in points])
    cov = ys.T@ys/4 - np.outer(ys.mean(0), ys.mean(0))
    h, d = full_mass_transform(points)
    actual = float(np.linalg.norm(h[:32, 32:])**2 + np.linalg.norm(d[:32, 32:])**2)
    predicted = float(np.sum(cov*gram))
    z = ys[:, 4]
    actual_majorana = float(np.linalg.norm(d[:32, 32:])**2)
    predicted_majorana = float(Ds2*(z@z/4-z.mean()**2))
    assert abs(actual-predicted)<2e-12 and abs(actual_majorana-predicted_majorana)<2e-13
    # Same entire aggregate hyperboloid vector, different internal singlet spread.
    a = .65
    length = R*np.sinh(a)
    witnesses = []
    aggregate = []
    for axis in (0, 4):
        ys = np.zeros((4, 5)); ys[:, axis] = length*np.array([1., -1., 1., -1.])
        ps = points_from_spatial(ys)
        X = np.array([old.hyper(p) for p in ps]); aggregate.append(X.sum(0))
        z = X[:, 5]; Z = float(z.sum()); P2 = float(z@z)
        hh, dd = full_mass_transform(ps)
        mass_leak = float(np.linalg.norm(dd[:32, 32:])**2)
        gamma = 4+P2/R**2
        cost = HBAR**2/(2*WEIGHT)*kraus_derivative_sum(Z)*gamma
        joint_cost = HBAR**2/(2*WEIGHT)*kraus_derivative_sum(Z)*(4+(Z*Z/4+4*mass_leak/Ds2)/R**2)
        assert abs(joint_cost-cost)<1e-14
        witnesses.append(dict(axis='Higgs' if axis==0 else 'singlet', count=4,
            amplitude=float(np.sqrt(-aggregate[-1]@old.ETA@aggregate[-1])/R),
            Z=Z, P2=P2, gamma=gamma, cost=cost, majorana_contrast_squared=mass_leak,
            average_mass_block_norm=float(np.linalg.norm(hh[:32,:32])+np.linalg.norm(dd[:32,:32])),
            joint_cost_error=abs(joint_cost-cost)))
    error=float(np.linalg.norm(aggregate[0]-aggregate[1]))
    assert error<1e-13
    assert witnesses[1]['cost']>witnesses[0]['cost']+.01
    assert witnesses[1]['majorana_contrast_squared']>.1
    return dict(original_modes=128, coefficient_gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
        full_Dirac_Majorana_contrast=actual, covariance_prediction=predicted,
        covariance_error=abs(actual-predicted), majorana_error=abs(actual_majorana-predicted_majorana),
        singlet_mass_coefficient_squared=Ds2, witnesses=witnesses,
        equal_aggregate_error=error, exact_model_parameter_status='inherited diagnostic Y, not empirical fit')


def moment_hierarchy():
    # Independent geodesic finite differences of the original 5D target Laplacian.
    p = np.array([.11, -.07, .12, .04, .37]); X = old.hyper(p); E = old.frame(p)
    z = X[5]; step = 2e-4
    laplace_rows=[]
    for k in (1, 2, 3, 4):
        actual=0.
        for a in range(5):
            plus=np.cosh(step/R)*X+R*np.sinh(step/R)*E[:,a]
            minus=np.cosh(step/R)*X-R*np.sinh(step/R)*E[:,a]
            actual += (plus[5]**k+minus[5]**k-2*z**k)/step**2
        predicted=(k*(k-1)*z**(k-2) if k>=2 else 0.)+k*(k+4)*z**k/R**2
        laplace_rows.append(dict(k=k, geodesic_difference=actual, analytic=predicted, error=abs(actual-predicted)))
    assert max(row['error'] for row in laplace_rows)<2e-6
    # Same count, A, m, P1, P2; P4 differs while each node has the same X0.
    length=R*np.sinh(.7); b=length*.55
    rows=[]; vectors=[]
    for heights in (np.array([b,0.,0.,0.]),np.array([b,b,0.,0.])/np.sqrt(2)):
        ys=[]
        for q in heights:
            y=np.array([np.sqrt(length*length-q*q),0.,0.,0.,q])
            ys.extend((y,-y))
        ps=points_from_spatial(np.array(ys)); XX=np.array([old.hyper(p) for p in ps])
        vectors.append(XX.sum(0)); zz=XX[:,5]
        P2=float(np.sum(zz**2)); P4=float(np.sum(zz**4))
        finite=sum(old.grad(lambda x: old.hyper(x)[5]**2,p)
                   @np.linalg.inv(old.original.metric(p))
                   @old.grad(lambda x: old.hyper(x)[5]**2,p) for p in ps)
        expected=4*(P2+P4/R**2)
        assert abs(finite-expected)<3e-8
        rows.append(dict(nodes=8,Z=float(zz.sum()),P2=P2,P4=P4,
            gamma_P2=expected,independent_gradient_error=float(abs(finite-expected))))
    assert np.linalg.norm(vectors[0]-vectors[1])<1e-13
    assert abs(rows[0]['P2']-rows[1]['P2'])<1e-13
    assert abs(rows[0]['gamma_P2']-rows[1]['gamma_P2'])>.1
    return dict(laplacian=laplace_rows, same_first_two_moments=rows,
                aggregate_error=float(np.linalg.norm(vectors[0]-vectors[1])),
                no_claim_all_finite_state_compressions_fail=True)


def run():
    deps=('research_note_574.md','research_note_598.md','research_note_642.md','research_note_707.md',
          'joint_geodesic_spatial_block.py','joint_fermion_gauss_completion.py',
          'round708_drafts/nested_mass_entry.md','round708_drafts/nested_mass_entry_results.json')
    return dict(date='2026-10-03',round=708,tests_run=3,failures=0,errors=0,
        geometry_and_budget=geometry_and_budget(),mass_and_common_witness=mass_and_common_witness(),
        moment_hierarchy=moment_hierarchy(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_H5_and_all_CAR=True,fixed_graph_and_equal_positive_node_weights=True,
            full_H_singlet_double_commutator_analytic=True,mass_statements_are_original_onsite_blocks=True,
            target_information_contract_is_additional=True,physical_spatial_dimension_not_inferred=True,
            no_full_Gibbs_spectrum_or_autonomous_apparatus=True,no_universal_cognition_no_go=True),
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    report=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    else:
        assert report==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:report[k] for k in ('round','tests_run','all_checks_passed')}))
