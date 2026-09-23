"""Round 400: compose near-unitary local region certificates.

Uses the frozen round-399 residual integrator on new four-qubit instances.
No old experiment is rerun. Multiplicative-domain/Stinespring tools are known.
"""
import argparse
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from operational_region_certificate_audit import (
    I, X, Y, Z, PS, word, kron_all, unitary, opnorm, comm, parts, certificate, adjoint)

TARGET = Path(__file__).with_name('joint_region_composition_audit_results.json')


def on_site(p, j, n=4):
    return kron_all([p if k == j else I for k in range(n)])


def sequences(h, t, mode):
    if mode == 'none':
        return [[unitary(h, t)]]
    rot = [on_site(unitary(.4*X+.7*Z, 1.2), 3)]
    damp = [on_site(np.diag([1., np.sqrt(.3)]), 3),
            on_site(np.array([[0, np.sqrt(.7)], [0, 0.]]), 3)]
    reset = [on_site(np.array([[1., 0], [0, 0]]), 3),
             on_site(np.array([[0., 1], [0, 0]]), 3)]
    middle = damp if mode == 'dissipative' else [on_site(X, 3)]
    last = reset if mode == 'dissipative' else rot
    return [[unitary(h, .2*t)], rot, [unitary(h, .35*t)],
            middle, [unitary(h, .15*t)], last, [unitary(h, .3*t)]]


def td(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def reference_output(sequence, state, refdim=3):
    rho = np.outer(state, state.conj())
    for kraus in sequence:
        lifted = [np.kron(k, np.eye(refdim)) for k in kraus]
        rho = sum(k@rho@k.conj().T for k in lifted)
    return np.einsum('aerbes->arbs', rho.reshape(4, 4, refdim, 4, 4, refdim)).reshape(4*refdim, 4*refdim)


def word_error(eps):
    ds = [np.sqrt(min(1., 2*e)) for e in eps]
    return float(sum(eps)+sum(a*b for a, b in itertools.combinations(ds, 2)))


def choi_from_adjoint(images):
    return sum(np.kron(np.kron(PS[a], PS[b]), images[(a, b)].T)
               for a, b in itertools.product(range(4), repeat=2))/(4*16)


def general_stinespring_audit():
    rng = np.random.default_rng(400)
    raw = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    h = (raw+raw.conj().T)/2
    u = unitary(h, .43)
    probability = .07
    v = np.vstack([np.sqrt(1-probability)*np.eye(8), np.sqrt(probability)*u])
    projection = v@v.conj().T
    embed = lambda a: np.kron(np.eye(2), a)
    cp = lambda a: v.conj().T@embed(a)@v
    generators = [word('XII'), word('IYI'), word('IIZ')]
    eps = [opnorm(cp(a)-a) for a in generators]
    defects = [opnorm(np.eye(8)-cp(a)@cp(a)) for a in generators]
    leaks = [opnorm((np.eye(16)-projection)@embed(a)@v) for a in generators]
    product = generators[0]@generators[1]@generators[2]
    image_product = cp(generators[0])@cp(generators[1])@cp(generators[2])
    return dict(isometry_error=opnorm(v.conj().T@v-np.eye(8)),
                local_errors=eps, schwarz_defects=defects, stinespring_leakages=leaks,
                product_leakage=opnorm((np.eye(16)-projection)@embed(product)@v),
                product_leakage_bound=sum(leaks),
                product_multiplicativity_error=opnorm(cp(product)-image_product),
                multiplicativity_bound=sum(a*b for a, b in itertools.combinations(leaks, 2)),
                full_word_error=opnorm(cp(product)-product), word_bound=word_error(eps))


@lru_cache(None)
def report():
    estimate = (.4*word('XIXI')+.3*word('IZZI')+.2*word('IIYY')
                +.17*word('ZIII')+.11*word('IXII')+.09*word('IIIZ'))
    true = estimate+.002*word('ZIIX')+.001*word('IYIZ')
    delta = .003
    regions = [(0, 2), (1, 2)]
    cases = []
    for t in (.03, .06, .12):
        target, eps = {}, {}
        for i in range(2):
            k, _ = parts(estimate, regions[i])
            for a in range(1, 4):
                p = on_site(PS[a], i)
                target[(i, a)] = unitary(k, -t)@p@unitary(k, t)
                eps[(i, a)] = certificate(estimate, regions[i], p, t, delta)['upper_bound']
        target_words, bounds = {}, {}
        for a, b in itertools.product(range(4), repeat=2):
            active = [(i, axis) for i, axis in enumerate((a, b)) if axis]
            value = np.eye(16, dtype=complex)
            for key in active:
                value = value@target[key]
            target_words[(a, b)] = value
            bounds[(a, b)] = word_error([eps[key] for key in active])
        images = {}
        rows = {}
        max_multiplicativity = 0.
        for mode in ('none', 'unitary', 'dissipative'):
            sequence = sequences(true, t, mode)
            ims = {key: adjoint(sequence, kron_all([PS[key[0]], PS[key[1]], I, I])) for key in bounds}
            images[mode] = ims
            errors = {str(key): opnorm(ims[key]-target_words[key]) for key in bounds if key != (0, 0)}
            defects = {str(key): opnorm(ims[key]-ims[(key[0], 0)]@ims[(0, key[1])])
                       for key in bounds if key[0] and key[1]}
            max_multiplicativity = max(max_multiplicativity, max(defects.values()))
            rows[mode] = dict(word_errors=errors, multiplicativity_errors=defects)
        choi = {mode: choi_from_adjoint(ims) for mode, ims in images.items()}
        rng = np.random.default_rng(401)
        state = rng.normal(size=48)+1j*rng.normal(size=48)
        state /= np.linalg.norm(state)
        ref = {mode: reference_output(sequences(true, t, mode), state) for mode in images}
        pair_distances = []
        for p, q in itertools.combinations(images, 2):
            pair_distances.append(dict(protocols=[p, q], choi_trace_distance=td(choi[p], choi[q]),
                                       random_reference_trace_distance=td(ref[p], ref[q]),
                                       maximum_joint_pauli_difference=max(opnorm(images[p][key]-images[q][key]) for key in bounds)))
        commutators = [opnorm(comm(target[(0, a)], target[(1, b)])) for a, b in itertools.product(range(1, 4), repeat=2)]
        comm_bounds = [2*bounds[(a, b)] for a, b in itertools.product(range(1, 4), repeat=2)]
        cases.append(dict(time=t, local_errors={str(key): e for key, e in eps.items()},
                          word_bounds={str(key): e for key, e in bounds.items() if key != (0, 0)},
                          protocol_audits=rows, pair_reference_distances=pair_distances,
                          joint_half_diamond_upper_bound=min(1., sum(bounds.values())),
                          uncapped_joint_bound=sum(bounds.values()),
                          target_commutators=commutators, target_commutator_bounds=comm_bounds,
                          maximum_actual_multiplicativity_defect=max_multiplicativity,
                          choi_trace_errors={mode: float(abs(np.trace(c)-1)) for mode, c in choi.items()},
                          choi_min_eigenvalues={mode: float(np.min(np.linalg.eigvalsh(c))) for mode, c in choi.items()}))
    # Exact common near-unitary reference: an entangling receiver operation, no outside coupling.
    exact_u = unitary(.31*word('XXII')+.2*word('ZYII'), .7)
    exact_defect = 0.
    for a, b in itertools.product(range(1, 4), repeat=2):
        pa, pb = on_site(PS[a], 0), on_site(PS[b], 1)
        transformed = lambda p: exact_u.conj().T@p@exact_u
        exact_defect = max(exact_defect, opnorm(transformed(pa@pb)-transformed(pa)@transformed(pb)))
    scale = []
    for m in (1, 2, 3, 4, 8):
        epsilon = 1e-4
        enumeration = sum(math.comb(m, k)*3**k*k*k*epsilon for k in range(1, m+1))
        closed = (9*m*m+3*m)*4**(m-2)*epsilon
        scale.append(dict(receivers=m, uniform_local_error=epsilon, pauli_word_bound_sum=enumeration,
                          closed_formula=closed, nontrivial=closed < 1.))
    return dict(round=400, scope='Near-unitary local operator certificates imply joint-receiver control bounds via known multiplicative-domain tools; no independent-noise assumption or preassigned spatial metric. Full channels include untouched references.',
                reused_frozen_integrator_from_round=399, calibration_error=opnorm(true-estimate),
                calibration_budget=delta, receiver_regions=[list(s) for s in regions], union_region=[0, 1, 2],
                cases=cases, general_stinespring=general_stinespring_audit(),
                exact_entangling_composition_defect=exact_defect, uniform_error_scaling=scale,
                merely_matching_marginals_sufficient=False, near_unitary_reference_required=True,
                actual_common_channel_required=True, independent_receiver_noise_assumed=False,
                ideal_local_references_assumed_to_commute=False,
                uniform_infinite_receiver_bound_proved=False, spatial_dimension_generated=False)


class Audit(unittest.TestCase):
    def test_local_certificates_on_new_overlapping_regions(self):
        self.assertLessEqual(report()['calibration_error'], report()['calibration_budget']+1e-13)
        for case in report()['cases']:
            for data in case['protocol_audits'].values():
                for key, value in data['word_errors'].items():
                    self.assertLessEqual(value, case['word_bounds'][key]+1e-12)

    def test_stinespring_leakage_and_product_estimate(self):
        s = report()['general_stinespring']
        self.assertLess(s['isometry_error'], 1e-12)
        for e, defect, leak in zip(s['local_errors'], s['schwarz_defects'], s['stinespring_leakages']):
            self.assertLessEqual(defect, min(1., 2*e)+1e-12)
            self.assertAlmostEqual(leak*leak, defect)
        self.assertLessEqual(s['product_leakage'], s['product_leakage_bound']+1e-12)
        self.assertLessEqual(s['product_multiplicativity_error'], s['multiplicativity_bound']+1e-12)
        self.assertLessEqual(s['full_word_error'], s['word_bound']+1e-12)

    def test_nonunitary_joint_multiplicativity_defect(self):
        for case in report()['cases']:
            for data in case['protocol_audits'].values():
                for a, b in itertools.product(range(1, 4), repeat=2):
                    eps0 = case['local_errors'][str((0, a))]
                    eps1 = case['local_errors'][str((1, b))]
                    bound = np.sqrt(min(1., 2*eps0)*min(1., 2*eps1))
                    self.assertLessEqual(data['multiplicativity_errors'][str((a, b))], bound+1e-12)
        self.assertGreater(report()['cases'][-1]['maximum_actual_multiplicativity_defect'], 1e-10)

    def test_incompatible_local_references_only_approximately_commute(self):
        for case in report()['cases']:
            for value, bound in zip(case['target_commutators'], case['target_commutator_bounds']):
                self.assertLessEqual(value, bound+1e-12)
        self.assertGreater(max(report()['cases'][-1]['target_commutators']), .001)

    def test_joint_channel_choi_and_reference_witnesses(self):
        for case in report()['cases']:
            for value in case['choi_trace_errors'].values():
                self.assertLess(value, 1e-12)
            for value in case['choi_min_eigenvalues'].values():
                self.assertGreater(value, -1e-12)
            for pair in case['pair_reference_distances']:
                self.assertLessEqual(pair['choi_trace_distance'], case['joint_half_diamond_upper_bound']+1e-12)
                self.assertLessEqual(pair['random_reference_trace_distance'], case['joint_half_diamond_upper_bound']+1e-12)

    def test_exact_entangling_case_has_zero_product_defect(self):
        self.assertLess(report()['exact_entangling_composition_defect'], 1e-12)
        self.assertEqual(word_error([0., 0., 0.]), 0.)

    def test_uniform_error_scaling_is_not_size_independent(self):
        for row in report()['uniform_error_scaling']:
            self.assertAlmostEqual(row['pauli_word_bound_sum'], row['closed_formula'])
        self.assertFalse(report()['uniform_error_scaling'][-1]['nontrivial'])

    def test_union_certificates_are_nontrivial_at_tested_windows(self):
        for case in report()['cases']:
            self.assertLess(case['uncapped_joint_bound'], 1.)
        self.assertFalse(report()['merely_matching_marginals_sufficient'])
        self.assertFalse(report()['independent_receiver_noise_assumed'])


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
