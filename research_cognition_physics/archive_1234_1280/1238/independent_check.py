"""Independent 1238 audit: exact old calibration and direct Hilbert stress.

The audited author's functions are never imported.  Only the frozen round 802
representation matrices/physical mass and the saved 900/907 inputs are reused.
Local off-shell jets test a coefficient, not an apparatus or continuum solution.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
TARGET = HERE / 'independent_results.json'
AUTHOR = ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']
HISTORY = [
    'archive_797_818/802/original_bff_vertex.py',
    'archive_894_922/907/material_reference_jets_results.json',
    'archive_894_922/900/current_reference_initial_budget_results.json',
    'archive_894_922/research_note_920.md',
]


def require(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def solve_fraction(matrix, rhs):
    """Independent exact elimination, solving J.T w=e_z directly."""
    rows = [[Fraction(x) for x in row] + [Fraction(b)]
            for row, b in zip(matrix, rhs)]
    n = len(rows)
    for col in range(n):
        pivot = next(r for r in range(col, n) if rows[r][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        p = rows[col][col]
        rows[col] = [x / p for x in rows[col]]
        for r in range(n):
            if r != col:
                p = rows[r][col]
                rows[r] = [x - p*y for x, y in zip(rows[r], rows[col])]
    return [row[-1] for row in rows]


def load_fixture():
    data = json.loads((BASE / HISTORY[1]).read_text('utf-8'), parse_float=Fraction)
    ini = json.loads((BASE / HISTORY[2]).read_text('utf-8'), parse_float=Fraction)
    row = data['rows'][0]
    initial = next(x for x in ini['rows'] if x['N'] == 24)
    j = [[Fraction(x) for x in r] for r in row['data']['old']['Jacobian']]
    psi = initial['psi']
    w = solve_fraction(list(zip(*j)), [0, 0, 0, 1])
    velocity = [Fraction(0), Fraction(0), -j[0][2]**2 / psi**8,
                -j[1][2]**2 / psi**8]
    zeta = sum(x*y for x, y in zip(w, velocity))
    c = 2*psi**4
    return row, psi, j, w, velocity, zeta, c


def calibration_check():
    row, psi, j, w, velocity, zeta, c = load_fixture()
    require(zeta != 0, 'fixed calibration must see metric velocity')
    coefficients = []
    for tau in [Fraction(-3, 100), Fraction(0), Fraction(2, 100)]:
        q = [[x + (tau*velocity[i] if mu == 0 else 0)
              for mu, x in enumerate(r)] for i, r in enumerate(j)]
        projected = [sum(w[i]*q[i][mu] for i in range(4)) for mu in range(4)]
        require(projected == [tau*zeta, 0, 0, 1], 'exact old calibration')
        coefficients.append(float(c*projected[0]**2))
    jf, pf = np.array(j, float), float(psi)
    gi = np.diag([-1., pf**-4, pf**-4, pf**-4])
    v = np.diag([0., 0., 1., 0.])
    # Differentiating the inverse metric in X3,X4 gives -g^-1 V g^-1.
    direct = np.array([0., 0., -jf[0]@gi@v@gi@jf[0],
                       -jf[1]@gi@v@gi@jf[1]])
    chain_error = float(np.max(abs(direct - np.array(velocity, float))))
    require(chain_error < 1e-18, 'inverse metric chain rule')
    require(c*zeta*zeta > 0, 'quadratic coefficient positive')
    return dict(status='PASS', old_N=24, saved_data_only_not_continuum_certificate=True,
                coefficient_C00_tau_squared=float(c*zeta*zeta),
                exact_coefficient_nonzero=True,
                inverse_metric_chain_error=chain_error,
                tiny_saved_spatial_h_derivatives_retained=True,
                fixed_B_no_refit=True, exact_C00_samples=coefficients)


def representation():
    sys.path.insert(0, str(ROOT / 'scripts'))
    from run_research_gr_response_active import GRResponseLayout, GRResponseRuntime
    with GRResponseRuntime(GRResponseLayout()).installed():
        spec = importlib.util.spec_from_file_location('independent_legacy802', BASE / HISTORY[0])
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def direct_stress_check():
    old = representation()
    row, p, _, _, _, a, c = load_fixture()
    psi, zeta, c = map(float, (p, a, c))
    frame = np.diag([1., psi**2, psi**2, psi**2])
    eta = np.diag([-1., 1., 1., 1.])
    metric = frame.T@eta@frame
    invmetric, invframe = np.linalg.inv(metric), np.linalg.inv(frame)
    volume = float(np.linalg.det(frame))
    alpha = np.array(old.ALPHA)
    gamma = np.array(old.GAMMA)
    clifford_error = max(float(np.max(abs(gamma[i]@gamma[j] + gamma[j]@gamma[i]
                                        - 2*eta[i,j]*np.eye(64))))
                          for i in range(4) for j in range(4))
    require(clifford_error < 1e-12, 'original Clifford representation')
    C = np.einsum('ma,aij->mij', invframe, alpha)
    Clower = np.einsum('mn,nij->mij', metric, C)
    values, vectors = np.linalg.eigh(alpha[1])
    u = vectors[:, np.argmax(values)]
    cx = float((u.conj()@C[1]@u).real)
    require(abs(cx-psi**-2) < 1e-13, 'original spatial eigenspinor')
    phi = np.array([0., float(row['data']['old']['value'][0]), 0., 0.,
                    float(row['data']['old']['value'][1])])
    mass = old.mass(phi)
    # Original representation but a different independent nonzero local jet.
    rng = np.random.default_rng(381238)
    internal = old.gauge_connection(rng)
    dg0 = rng.normal(size=(4, 4, 4))*.07
    dg0 += dg0.transpose(0, 2, 1)
    dpsi0 = (rng.normal(size=(4, 64))+1j*rng.normal(size=(4, 64)))*.031
    factor = float(Fraction(2, 7)*Fraction(3, 11))

    def spin_connection(dg):
        # Independent explicit Levi-Civita formula and tetrad compatibility.
        gamma_coord = np.empty((4, 4, 4))
        for mu in range(4):
            for r in range(4):
                for nu in range(4):
                    gamma_coord[mu,r,nu] = sum(invmetric[r,s]*(dg[mu,s,nu]
                        + dg[nu,s,mu] - dg[s,mu,nu]) for s in range(4))/2
        result = []
        metricity = 0.
        for mu in range(4):
            dframe = frame@invmetric@dg[mu]/2
            omega = frame@gamma_coord[mu]@invframe-dframe@invframe
            lower = eta@omega
            metricity = max(metricity, float(np.max(abs(lower+lower.T))))
            result.append(sum((lower[a,b]*gamma[a]@gamma[b]/4
                               for a in range(4) for b in range(4)),
                              np.zeros((64,64),complex)))
        require(metricity < 2e-13, 'Lorentz spin connection')
        return np.array(result)

    def evaluate(tau, k, keep_density=True, source_factor=factor):
        dg = dg0.copy()
        dg[0,2,2] += tau
        spin = spin_connection(dg)
        deriv = dpsi0.copy()
        deriv[1] += 1j*k*u
        D = deriv + np.einsum('mij,j->mi', spin+internal, u)
        ell = -sum(float(np.imag(u.conj()@C[mu]@D[mu])) for mu in range(4))
        ell -= float((u.conj()@mass@u).real)
        # Off-shell Hilbert stress of 1/2 times the Hermitian Dirac density.
        stress = np.empty((4,4))
        for mu in range(4):
            for nu in range(4):
                stress[mu,nu] = np.imag(u.conj()@Clower[mu]@D[nu]
                                    + u.conj()@Clower[nu]@D[mu])/4
                if keep_density:
                    stress[mu,nu] += metric[mu,nu]*ell/2
        stress_up = invmetric@stress@invmetric
        z = np.array([tau*zeta, 0., 0., 1.])
        Cref = c*np.outer(z,z)
        result = -source_factor*volume*np.sum(Cref*stress_up)/2
        return float(result), float(np.linalg.norm(spin)), stress[0,0]

    expected = -factor*volume*c*zeta*zeta*cx/4
    errors, rows, density_deleted = [], [], []
    for tau_step in [.011, .027, .049]:
        slopes, wrong = [], []
        for tau in [-tau_step, 0., tau_step]:
            slopes.append(evaluate(tau, 1.7)[0]-evaluate(tau, .7)[0])
            wrong.append(evaluate(tau, 1.7, False)[0]-evaluate(tau, .7, False)[0])
        coeff = (slopes[0]+slopes[2]-2*slopes[1])/(2*tau_step*tau_step)
        deleted = (wrong[0]+wrong[2]-2*wrong[1])/(2*tau_step*tau_step)
        errors.append(abs(coeff-expected)); density_deleted.append(deleted)
        rows.append(dict(step=tau_step, coefficient=coeff, error=abs(coeff-expected)))
    require(max(errors) < 2e-10, 'direct full-stress mixed coefficient')
    require(max(abs(x) for x in density_deleted) < 2e-10, 'density term essential')
    require(evaluate(.03, 1.7, source_factor=0)[0] == 0, 'delete read coupling')
    # The slope of the canonical momentum response is the derivative of this
    # quadratic coefficient: its k-dependent part is 2*expected*tau.
    t0, h = .021, .0001
    def kslope(t):
        return evaluate(t, 1.7)[0]-evaluate(t, .7)[0]
    legendre = (kslope(t0+h)-kslope(t0-h))/(2*h)
    legendre_error = abs(legendre-2*expected*t0)
    require(legendre_error < 2e-10, 'nonzero original metric momentum change')
    require(evaluate(0., 1.7)[1] > 0 and np.linalg.norm(mass) > 0
            and np.linalg.norm(internal) > 0, 'full connection and mass retained')
    T00slope = evaluate(0., 1.7)[2]-evaluate(0., .7)[2]
    require(abs(T00slope-cx/2) < 2e-13, 'spatial derivative survives off-shell T00')
    return dict(status='PASS', original_Nambu_dimension=64, nambu_half_weight=.5,
                independent_method='direct complete off-shell Hilbert stress, no author import',
                nonzero_mass_LC_and_internal_connection_retained=True,
                original_clifford_residual=clifford_error,
                full_stress_rows=rows, expected_coefficient=expected,
                omitted_density_coefficients=density_deleted,
                spatial_derivative_T00_slope=T00slope,
                metric_momentum_k_slope=legendre, metric_momentum_k_slope_error=legendre_error,
                original_background_jet_reconstructed=False,
                offshell_coefficient_not_actual_readout=True,
                canonical_redefinition_no_go=False)


def digest(paths):
    return {str(p.relative_to(BASE)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def run(bind=True):
    watched = [BASE/p for p in HISTORY]
    if bind:
        watched += [HERE/p for p in AUTHOR]
    before = digest(watched)
    result = dict(round=1238, status='PASS', independent_scientific_groups=2,
                  fixed_original_calibration=calibration_check(),
                  full_stress_and_original_Legendre=direct_stress_check(),
                  direct_CT_dictionary_inheritance_refuted=True,
                  GR_refuted=False, CT_role_adopted=False,
                  actual_total_quantum_process_certified=False,
                  source_sha256=before)
    require(before == digest(watched), 'read-only historical and author audit')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-once', action='store_true')
    parser.add_argument('--diagnostic', action='store_true')
    args = parser.parse_args()
    require(not (args.write_once and args.diagnostic), 'choose one mode')
    output = run(bind=not args.diagnostic)
    if args.diagnostic:
        print(json.dumps(output, ensure_ascii=False, indent=2))
    elif args.write_once:
        with TARGET.open('x',encoding='utf-8') as stream:
            json.dump(output,stream,ensure_ascii=False,indent=2); stream.write('\n')
        print('WROTE_ONCE: 1238 independent_results.json; two independent groups PASS')
    else:
        require(output == json.loads(TARGET.read_text('utf-8')), 'READ_ONLY_MATCH')
        print('READ_ONLY_MATCH: 1238 independent; two groups PASS')
