"""Round 399: metric-free, finite-time certificates for a given receiver.

Known Duhamel/contractivity tools, applied to a reconstructed finite generator.
No exact support threshold, input spatial distance, or dimension inference.
Bounds cover deterministic operations confined to the region complement.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('operational_region_certificate_audit_results.json')
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
PS = (I, X, Y, Z)


def kron_all(xs):
    out = np.ones((1, 1), complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def word(s):
    return kron_all([dict(I=I, X=X, Y=Y, Z=Z)[c] for c in s])


def opnorm(a):
    return float(np.linalg.norm(a, ord=2))


def comm(a, b):
    return a@b-b@a


def unitary(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def expectation(h, keep):
    """Normalized conditional expectation, obtained by explicit Pauli twirls."""
    n = int(round(np.log2(len(h))))
    out = h.copy()
    for j in set(range(n))-set(keep):
        out = sum(p@out@p for p in
                  (kron_all([q if k == j else I for k in range(n)]) for q in PS))/4
    return out


def parts(h, region):
    n = int(round(np.log2(len(h))))
    k = expectation(h, region)
    outside = expectation(h, set(range(n))-set(region))
    scalar = np.trace(h)/len(h)*np.eye(len(h))
    return k, h-k-outside+scalar


def certificate(h, region, obs, t, delta=0., steps=128):
    k, boundary = parts(h, region)
    vals, vecs = np.linalg.eigh(k)
    rotated = vecs.conj().T@obs@vecs
    def residual(s):
        phase = np.exp(1j*s*vals)
        evolved = vecs@(phase[:, None]*rotated*phase.conj()[None, :])@vecs.conj().T
        return opnorm(comm(boundary, evolved))
    dt = t/steps
    midpoint = dt*sum(residual((j+.5)*dt) for j in range(steps))
    lipschitz = 4*opnorm(boundary)*opnorm(k)*opnorm(obs)
    remainder = lipschitz*t*t/(4*steps)
    return dict(midpoint_integral=float(midpoint), integration_remainder=float(remainder),
                calibration_allowance=float(2*delta*t*opnorm(obs)),
                upper_bound=float(midpoint+remainder+2*delta*t*opnorm(obs)),
                residual_at_zero=residual(0.), lipschitz_bound=float(lipschitz))


def local_kraus(matrices, site, n=3):
    return [kron_all([matrix if j == site else I for j in range(n)]) for matrix in matrices]


def operations(h, t, mode):
    if mode == 'none':
        return [[unitary(h, t)]]
    rotation = local_kraus([unitary(.37*X+.59*Y, 1.7)], 2)
    flip = local_kraus([X], 2)
    damp = local_kraus([np.diag([1, np.sqrt(.27)]),
                       np.array([[0, np.sqrt(.73)], [0, 0]])], 2)
    reset = local_kraus([np.array([[1, 0], [0, 0]]),
                        np.array([[0, 1], [0, 0]])], 2)
    middle = damp if mode == 'dissipative' else flip
    last = reset if mode == 'dissipative' else rotation
    return [[unitary(h, .21*t)], rotation, [unitary(h, .32*t)],
            middle, [unitary(h, .19*t)], last, [unitary(h, .28*t)]]


def adjoint(sequence, obs):
    out = obs.copy()
    for kraus in reversed(sequence):
        out = sum(k.conj().T@out@k for k in kraus)
    return out


def apply_reference(sequence, rho, refdim):
    for kraus in sequence:
        lifted = [np.kron(k, np.eye(refdim)) for k in kraus]
        rho = sum(k@rho@k.conj().T for k in lifted)
    # Three system qubits, receiver is 0, retain the entire untouched reference.
    return np.einsum('aerbes->arbs', rho.reshape(2, 4, refdim, 2, 4, refdim)).reshape(2*refdim, 2*refdim)


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)


def receiver_axes(n=3):
    return [kron_all([p]+[I]*(n-1)) for p in (X, Y, Z)]


@lru_cache(None)
def report():
    estimate = .7*word('XXI')+.5*word('IZZ')+.2*word('ZII')+.13*word('IIX')
    true = estimate+.01*word('XIZ')
    axes = receiver_axes()
    region = (0, 1)
    cases = []
    for t in (.15, .3, .6):
        certs = [certificate(estimate, region, p, t, .01) for p in axes]
        k, boundary = parts(estimate, region)
        ideal = [unitary(k, t)]
        local_ideal = [adjoint([ideal], p) for p in axes]
        errors = {}
        for mode in ('none', 'unitary', 'dissipative'):
            sequence = operations(true, t, mode)
            errors[mode] = [opnorm(adjoint(sequence, p)-q) for p, q in zip(axes, local_ideal)]
        # Choi input plus independent complex entangled input; not a diamond optimizer.
        maximally_entangled = np.eye(8).reshape(-1)/np.sqrt(8)
        rng = np.random.default_rng(399)
        random_state = rng.normal(size=64)+1j*rng.normal(size=64)
        random_state /= np.linalg.norm(random_state)
        ref_errors = {}
        for label, state in [('choi', maximally_entangled), ('random_complex', random_state)]:
            rho = np.outer(state, state.conj())
            target = apply_reference([ideal], rho, 8)
            ref_errors[label] = {mode: trace_distance(apply_reference(operations(true, t, mode), rho, 8), target)
                                 for mode in ('none', 'unitary', 'dissipative')}
        cases.append(dict(time=t, axis_certificates=certs, actual_axis_errors=errors,
                          complete_channel_half_diamond_upper_bound=sum(c['upper_bound'] for c in certs)/2,
                          reference_output_trace_distances=ref_errors))
    # Independently compare certified midpoint integration with Gauss quadrature.
    k, boundary = parts(estimate, region)
    nodes, weights = np.polynomial.legendre.leggauss(96)
    quadratures = []
    for p in axes:
        t = .6
        integral = t/2*sum(w*opnorm(comm(boundary, unitary(k, -t*(x+1)/2)@p@unitary(k, t*(x+1)/2)))
                             for x, w in zip(nodes, weights))
        coarse = certificate(estimate, region, p, t, steps=64)
        fine = certificate(estimate, region, p, t, steps=256)
        quadratures.append(dict(gauss_integral=float(integral), coarse=coarse, fine=fine))
    # A missed nonzero term: no exact support gap is used by the certificate.
    decoupled = .4*word('ZII')+.3*word('IXI')+.7*word('IIZ')
    weak_true = decoupled+.02*word('XIZ')
    weak_t = .8
    weak = [certificate(decoupled, (0,), p, weak_t, .02) for p in axes]
    wk, _ = parts(decoupled, (0,))
    weak_actual = [opnorm(adjoint(operations(weak_true, weak_t, 'dissipative'), p)
                          -adjoint([[unitary(wk, weak_t)]], p)) for p in axes]
    # Conditional expectation identity, independent support expansion and covariance.
    twirl_residual = opnorm(expectation(estimate, region)-(.7*word('XXI')+.2*word('ZII')))
    w = kron_all([unitary(X+.3*Z, .21), unitary(Y, -.37), unitary(Z, .81)])
    covariance = opnorm(expectation(w@estimate@w.conj().T, region)-w@expectation(estimate, region)@w.conj().T)
    cov_cert = certificate(w@estimate@w.conj().T, region, w@axes[1]@w.conj().T, .3, .01)
    baseline = certificate(estimate, region, axes[1], .3, .01)
    # Envelopes need not be optimal sets. Enumerate all receiver-containing subsets.
    regions = []
    for flags in itertools.product((False, True), repeat=2):
        selected = (0,)+tuple(i+1 for i, include in enumerate(flags) if include)
        c = [certificate(estimate, selected, p, .3, .01) for p in axes]
        regions.append(dict(region=list(selected), complete_upper_bound=sum(v['upper_bound'] for v in c)/2))
    return dict(round=399, scope='Finite-time operational regions from a reconstructed generator, using known Duhamel contraction tools; no input metric or exact support gap. Arbitrary deterministic operations confined to the complement, full unknown inputs and untouched reference.',
                calibrated_operator_error=opnorm(true-estimate), cases=cases, quadrature_audit=quadratures,
                weak_missing_term=dict(delta=.02, time=weak_t, certificates=weak, actual_errors=weak_actual,
                                       zero_error_assumption_would_be_false=True),
                conditional_expectation_error=twirl_residual, local_unitary_covariance_error=covariance,
                certificate_covariance_error=abs(cov_cert['upper_bound']-baseline['upper_bound']),
                region_enumeration=regions,
                assumed_spatial_distance=False, exact_support_gap_required=False,
                all_background_states_covered_by_proof=True, arbitrary_outside_control_covered_by_proof=True,
                postselected_conditional_state_guarantee=False,
                complete_channel_bound_proved_not_numerically_optimized=True,
                finite_geometry_or_dimension_generated=False,
                floating_point_roundoff_interval_certified=False)


class Audit(unittest.TestCase):
    def test_conditional_expectation_and_calibrated_operator_error(self):
        self.assertLess(report()['conditional_expectation_error'], 1e-13)
        self.assertAlmostEqual(report()['calibrated_operator_error'], .01)

    def test_time_integral_remainder_and_independent_quadrature(self):
        for case in report()['quadrature_audit']:
            for grid in ('coarse', 'fine'):
                item = case[grid]
                self.assertLessEqual(abs(item['midpoint_integral']-case['gauss_integral']),
                                     item['integration_remainder']+1e-12)
            self.assertAlmostEqual(case['coarse']['integration_remainder'], 4*case['fine']['integration_remainder'])

    def test_natural_evolution_and_repeated_outside_unitaries(self):
        for case in report()['cases']:
            for mode in ('none', 'unitary'):
                for actual, certificate_ in zip(case['actual_axis_errors'][mode], case['axis_certificates']):
                    self.assertLessEqual(actual, certificate_['upper_bound']+1e-12)

    def test_nonunitary_outside_channels_do_not_invalidate_bound(self):
        for case in report()['cases']:
            for actual, certificate_ in zip(case['actual_axis_errors']['dissipative'], case['axis_certificates']):
                self.assertLessEqual(actual, certificate_['upper_bound']+1e-12)
        for channel in operations(np.zeros((8, 8)), .3, 'dissipative'):
            self.assertLess(opnorm(sum(k.conj().T@k for k in channel)-np.eye(8)), 1e-12)

    def test_unknown_reference_outputs_under_complete_channel_bound(self):
        for case in report()['cases']:
            for row in case['reference_output_trace_distances'].values():
                for value in row.values():
                    self.assertLessEqual(value, case['complete_channel_half_diamond_upper_bound']+1e-12)

    def test_weak_omitted_link_requires_error_allowance_not_exact_graph(self):
        case = report()['weak_missing_term']
        self.assertGreater(max(case['actual_errors']), .01)
        for actual, certificate_ in zip(case['actual_errors'], case['certificates']):
            self.assertAlmostEqual(certificate_['midpoint_integral'], 0.)
            self.assertLessEqual(actual, certificate_['upper_bound']+1e-12)

    def test_local_basis_covariance(self):
        self.assertLess(report()['local_unitary_covariance_error'], 1e-12)
        self.assertLess(report()['certificate_covariance_error'], 1e-12)

    def test_region_and_horizon_certificates(self):
        rows = report()['region_enumeration']
        self.assertEqual(len(rows), 4)
        full = next(row for row in rows if len(row['region']) == 3)
        self.assertAlmostEqual(full['complete_upper_bound'], 3*.01*.3)
        self.assertLess(next(row for row in rows if row['region'] == [0, 1])['complete_upper_bound'], .15)
        for a, b in zip(report()['cases'], report()['cases'][1:]):
            self.assertLess(a['complete_channel_half_diamond_upper_bound'], b['complete_channel_half_diamond_upper_bound'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checks.testsRun, failures=len(checks.failures), errors=len(checks.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
