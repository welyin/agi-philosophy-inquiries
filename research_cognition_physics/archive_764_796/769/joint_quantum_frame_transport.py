"""769: original nonlinear frame, same linear BRST object, and mean/source jets.

Finite algebra calibrations only. No continuum loop source, renormalized
Jacobian, or actual interacting state is computed by these checks.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_brst_relative_source as old
import joint_covariant_gauge_complex as cov
import joint_scalar_propagation_matching as scalar

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_quantum_frame_transport_results.json'
PHI = np.array([.13, .64, -.11, .08, .37])


def maximum(x):
    return float(np.max(np.abs(x)))


def F(phi):
    return 2. - phi @ phi / 6.


def frame(q):
    return np.r_[F(q[10:]) * q[:10], q[10:]]


def frame_jets(q, v):
    g, phi, h, chi = q[:10], q[10:], v[:10], v[10:]
    t = np.eye(15)
    t[:10, :10] *= F(phi)
    t[:10, 10:] = np.outer(g, -phi/3.)
    dt = np.zeros((15, 15))
    dt[:10, :10] = np.eye(10) * (-phi @ chi/3.)
    dt[:10, 10:] = np.outer(h, -phi/3.) + np.outer(g, -chi/3.)
    ddt = np.zeros((15, 15))
    ddt[:10, :10] = -np.eye(10) * (chi @ chi/3.)
    ddt[:10, 10:] = -2*np.outer(h, chi)/3.
    k = np.r_[(-phi @ chi/3.)*h - (chi @ chi/6.)*g, np.zeros(5)]
    return t, dt, ddt, k


def full_linear_transport():
    # Use the original metric packing and original 111-field principal blocks.
    q = np.r_[cov.pack(cov.ETA/F(PHI)), PHI]
    t15 = frame_jets(q, np.zeros(15))[0]
    t = np.eye(111)
    indices = np.r_[np.arange(10), np.arange(58, 63)]
    t[np.ix_(indices, indices)] = t15
    ti = np.linalg.inv(t)
    _, _, _, _, _, _, le, ge, _, be, _ = old.setup([1.3, .2, -.3, .4])
    # Density Hessians and their dual Green matrices: no hidden pairing flip.
    he = be @ le
    hj = t.T @ he @ t
    grj = ti @ np.linalg.inv(he) @ ti.T
    gj = ti @ ge @ t
    bj = t.T @ be @ t
    lj = np.linalg.solve(bj, hj)
    gjstar = np.linalg.solve(bj, gj.conj().T @ bj)
    errors = dict(
        inverse=maximum(hj @ grj - np.eye(111)),
        BRST_nilpotency=maximum(gj @ gj),
        BRST_adjoint_identity=maximum(gjstar @ lj-lj @ gj),
        Hessian_self_adjoint=maximum(hj-hj.conj().T),
    )
    # A positive physical scalar polarization at one null principal symbol.
    _, _, _, _, _, _, ln, _, _, bn, _ = old.setup([1., 1., 0., 0.])
    ue = np.zeros(111)
    ue[58:63] = [.27, -.14, .31, .18, -.09]
    uj = ti @ ue
    we = np.outer(ue, ue)
    wj = ti @ we @ ti.T
    test_e = np.linspace(-.12, .24, 111)
    test_j = t.T @ test_e
    errors.update(
        scalar_principal_solution=maximum((bn @ ln) @ we),
        pulled_principal_solution=maximum((t.T @ bn @ ln @ t) @ wj),
        same_test_covariance=abs(float(test_j @ wj @ test_j-test_e @ we @ test_e)),
        physical_scalar_pairing=abs(float((uj @ bj @ uj-ue @ be @ ue).real)),
    )
    assert max(errors.values()) < 2e-13
    assert float((ue @ be @ ue).real) > 0
    return dict(errors=errors, F=float(F(PHI)),
                metric_scalar_mixing_norm=maximum(t[:10, 58:63]),
                scalar_physical_pairing=float((ue @ be @ ue).real),
                scope='Full 111-field principal algebra plus one physical scalar polarization; no actual spacetime covariance or loop integral is sampled.')


def mean_and_record():
    q = np.r_[cov.pack(cov.ETA/F(PHI)), PHI]
    h = cov.pack(np.array([[.08, .02, 0., 0.], [.02, -.03, .01, 0.],
                          [0., .01, .06, -.01], [0., 0., -.01, .04]]))
    chi = np.array([.27, -.14, .31, .18, -.09])
    v = np.r_[h, chi]
    t, _, _, k = frame_jets(q, v)
    q0, u = frame(q), t @ v
    rows = []
    # Symmetric two-atom calibration: exact mean of a cubic frame map.
    # It is not a substitute for the inherited continuum Hadamard state.
    for a in [.2, .1, .05]:
        plus, minus = frame(q+a*v), frame(q-a*v)
        mean = (plus+minus)/2
        error = maximum((mean-q0)/a**2-k)
        record = .5*(plus[0]**2+minus[0]**2)
        predicted_coefficient = 2*q0[0]*k[0]+u[0]**2
        remainder = (record-q0[0]**2)/a**2-predicted_coefficient
        rows.append(dict(amplitude=a, mean_contact_error=error,
                         record_scaled_remainder=float(remainder),
                         minimum_F=min(float(F(q[10:]+a*chi)), float(F(q[10:]-a*chi)))))
        assert error < 3e-13
        assert np.linalg.eigvalsh(cov.unpack(q[:10]+a*h))[0] < 0
    assert maximum(k[:10]) > 1e-3
    assert abs(rows[-1]['record_scaled_remainder']) < abs(rows[0]['record_scaled_remainder'])/12
    # Inverse-frame curvature cancels the forward mean contact (second-order chain rule).
    def inverse_tangent(y):
        phi = y[10:]
        ti = np.eye(15); ti[:10, :10] /= F(phi)
        ti[:10, 10:] = np.outer(y[:10], phi/(3*F(phi)**2))
        return ti
    ti = inverse_tangent(q0)
    # For f^{-1}: q_g=Q_g/F, scalar variables unchanged.
    fp = -PHI @ chi/3.
    inv_k_g = q0[:10]*(fp*fp/F(PHI)**3+(chi@chi)/(6*F(PHI)**2))-u[:10]*fp/F(PHI)**2
    inverse_contact_error = maximum(ti@k+np.r_[inv_k_g, np.zeros(5)])
    assert inverse_contact_error < 1e-14
    return dict(rows=rows, metric_mean_contact=k[:10].tolist(),
                omitted_record_mean_term=float(2*q0[0]*k[0]),
                record_covariance_term=float(u[0]**2),
                inverse_contact_error=inverse_contact_error,
                scope='Original metric/5-scalar nonlinear map and a composite-coordinate test; the latter is not asserted to be a diffeomorphism-invariant physical record.')


class Jet:
    """Exact truncated bivariate Taylor algebra, x-degree <=1, y-degree <=2."""
    def __init__(self, value=0., coeff=None):
        self.c = np.zeros((2, 3)) if coeff is None else np.array(coeff, float)
        if coeff is None:
            self.c[0, 0] = value

    @staticmethod
    def lift(v):
        return v if isinstance(v, Jet) else Jet(v)

    def __add__(self, other):
        return Jet(coeff=self.c+self.lift(other).c)

    __radd__ = __add__

    def __neg__(self):
        return Jet(coeff=-self.c)

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        out = np.zeros((2, 3))
        for i, j, k, l in itertools.product(range(2), range(3), range(2), range(3)):
            if i+k < 2 and j+l < 3:
                out[i+k, j+l] += self.c[i, j]*other.c[k, l]
        return Jet(coeff=out)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        base = self.c[0, 0]
        assert base > 0
        z = (self-base)*(1/base)
        power, result, factor = Jet(1.), Jet(1.), 1.
        for n in range(1, 4):
            power = power*z
            factor *= (exponent-n+1)/n
            result += factor*power
        return base**exponent*result


def density_jet(q, test, v, einstein):
    coords = []
    for qi, ai, vi in zip(q, test, v):
        c = np.zeros((2, 3)); c[0, 0] = qi; c[1, 0] = ai; c[0, 1] = vi
        coords.append(Jet(coeff=c))
    metric = [[sum(cov.T[a, i, j]*coords[a] for a in range(10)) for j in range(4)] for i in range(4)]
    det = Jet()
    for perm in itertools.permutations(range(4)):
        inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
        term = Jet((-1.)**inversions)
        for i in range(4):
            term = term*metric[i][perm[i]]
        det += term
    phi = coords[10:]
    matrix, vacuum, _ = scalar.parameters()
    delta = [sum(x*x for x in phi[:4])-vacuum[0], phi[4]*phi[4]-vacuum[1]]
    potential = sum(.25*matrix[i, j]*delta[i]*delta[j] for i in range(2) for j in range(2))
    density = -((-det)**.5)*potential
    if einstein:
        density = density*(2.-sum(x*x for x in phi)*(1/6))**-2
    return density.c


def derivative_vectors(q, v, einstein):
    jets = np.array([density_jet(q, e, v, einstein) for e in np.eye(15)])
    return jets[:, 1, 0], jets[:, 1, 1], jets[:, 1, 2]


def original_density_source_chain():
    # All 15 metric/scalar directions are differentiated independently in the
    # original Jordan and Einstein potential densities. This point is off shell.
    q = np.r_[cov.pack(cov.ETA/F(PHI)), PHI]
    v = np.r_[np.linspace(-.04, .06, 10), [.27, -.14, .31, .18, -.09]]
    t, dt, ddt, k = frame_jets(q, v)
    qe, u = frame(q), t@v
    ej, hjv, jj = derivative_vectors(q, v, False)
    ee, heu, je = derivative_vectors(qe, u, True)
    _, hek, _ = derivative_vectors(qe, k, True)
    contacts = [t.T@hek, dt.T@heu, .5*ddt.T@ee]
    corrected = t.T@je+sum(contacts)
    errors = dict(
        first_variation=maximum(ej-t.T@ee),
        second_variation=maximum(hjv-t.T@heu-dt.T@ee),
        third_source_chain=maximum(jj-corrected),
    )
    naive = maximum(jj-t.T@je)
    assert max(errors.values()) < 5e-15 and naive > 1e-6
    return dict(errors=errors, naive_source_pullback_error=naive,
                mean_contact_norm=maximum(contacts[0]),
                linear_equation_contact_norm=maximum(contacts[1]),
                background_equation_contact_norm=maximum(contacts[2]),
                scope='Exact bivariate Taylor coefficients of original potential densities, at an off-shell point; full kinetic/ghost loop sources and physical wave solutions are not numerically computed.')


def run():
    result = dict(round=769, tests_run=3, failures=0, errors=0,
                  linear_BRST_frame=full_linear_transport(),
                  nonlinear_mean_and_record=mean_and_record(),
                  source_chain=original_density_source_chain())
    deps = ('research_note_583.md', 'research_note_600.md', 'research_note_619.md',
            'research_note_664.md', 'research_note_767.md', 'research_note_768.md',
            'joint_brst_relative_source.py', 'joint_covariant_gauge_complex.py',
            'joint_scalar_propagation_matching.py',
            'round769_drafts/research_note_769_working.md')
    result['dependency_hashes'] = {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps}
    result['scope'] = ('Same on-shell free physical state and BRST extension transport through the original F-positive frame map. '
                       'A specified smooth physical covariance difference induces a necessary nonlinear mean contact; '
                       'the complete relative source and formal response transform together. Absolute quantum Ward, '
                       'renormalized measures and all-sector source existence are not proved.')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--write-results', action='store_true')
    args = p.parse_args(); result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('round', 'tests_run', 'failures', 'errors')}))
