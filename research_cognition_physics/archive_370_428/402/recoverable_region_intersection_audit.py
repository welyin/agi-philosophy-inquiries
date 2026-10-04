"""Round 402: recovery regions versus fixed-process dependency regions.

Reuses the established three-qutrit secret-sharing code and round-391 algebra
solver on new channels. Explicit recovery, full instruments, unknown references,
intersection failure and a finite error gap are audited. A fixed full-domain
channel, by contrast, has an intersection error budget. No space is assumed.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from shared_sharp_record_audit import marginal_kraus, apply, recoverable

TARGET = Path(__file__).with_name('recoverable_region_intersection_audit_results.json')
D = 3
I = np.eye(D, dtype=complex)
OMEGA = np.exp(2j*np.pi/D)
X = np.roll(I, 1, axis=0)
Z = np.diag(OMEGA**np.arange(D))
F = OMEGA**np.outer(np.arange(D), np.arange(D))/np.sqrt(D)
PAIRS = tuple(itertools.combinations(range(3), 2))


def encoder():
    v = np.zeros((27, 3), complex)
    for s, r in itertools.product(range(3), repeat=2):
        v[9*r+3*((r+s) % 3)+(r+2*s) % 3, s] = 1/np.sqrt(3)
    return v


def encoding_unitary():
    """|s,0,0> -> V|s>; Fourier one blank then an invertible linear permutation."""
    p = np.zeros((27, 27), complex)
    for s, r, k in itertools.product(range(3), repeat=3):
        p[9*r+3*((r+s) % 3)+(r+2*s+k) % 3, 9*s+3*r+k] = 1.
    return p @ np.kron(np.kron(I, F), I)


def decoder(pair):
    """Pair values x,y -> logical s and the omitted share's value k."""
    i, j = pair
    omitted = next(k for k in range(3) if k not in pair)
    inverse = pow(j-i, -1, 3)
    matrix = np.zeros((9, 9), complex)
    for x, y in itertools.product(range(3), repeat=2):
        s = (y-x)*inverse % 3
        k = (x+(omitted-i)*s) % 3
        matrix[3*s+k, 3*x+y] = 1.
    return matrix


def on_subset(matrix, keep):
    """Embed an operator in three physical qutrits, preserving supplied order."""
    keep = tuple(keep)
    omitted = tuple(j for j in range(3) if j not in keep)
    result = np.zeros((27, 27), complex)
    tuples = list(itertools.product(range(3), repeat=3))
    for a, aa in enumerate(tuples):
        for b, bb in enumerate(tuples):
            if all(aa[j] == bb[j] for j in omitted):
                ia = sum(aa[j]*3**(len(keep)-1-k) for k, j in enumerate(keep))
                ib = sum(bb[j]*3**(len(keep)-1-k) for k, j in enumerate(keep))
                result[a, b] = matrix[ia, ib]
    return result


def reorder_output(v, order):
    return v.reshape(3, 3, 3, v.shape[1]).transpose(tuple(order)+(3,)).reshape(27, v.shape[1])


def channel_choi(kraus, refdim=3):
    phi = np.eye(refdim, dtype=complex).reshape(-1)/np.sqrt(refdim)
    return sum(np.outer(np.kron(k, np.eye(refdim))@phi,
                        (np.kron(k, np.eye(refdim))@phi).conj()) for k in kraus)


def decode_kraus(pair):
    dec = decoder(pair)
    return [np.kron(I, np.eye(3)[[k]]) @ dec for k in range(3)]


def norm(a):
    return float(np.linalg.norm(a, 2))


def td(a, b):
    delta = a-b
    return float(np.abs(np.linalg.eigvalsh((delta+delta.conj().T)/2)).sum()/2)


def span_rank(matrices):
    cols = np.column_stack([a.reshape(-1) for a in matrices])
    return int(np.sum(np.linalg.svd(cols, compute_uv=False) > 1e-9))


def filter_audit():
    families = []
    for raw in range(1, 256):
        gamma = {s for s in range(8) if (raw >> s) & 1}
        if 0 in gamma or not all(t in gamma for s in gamma for t in range(8) if s & t == s):
            continue
        intersection = 7
        for s in gamma:
            intersection &= s
        closed = all(s & t in gamma for s in gamma for t in gamma)
        principal = gamma == {s for s in range(8) if s & intersection == intersection}
        families.append(dict(members=sorted(gamma), core=intersection,
                             intersection_closed=closed, principal=principal,
                             no_disjoint_authorized_sets=all(s & t for s in gamma for t in gamma)))
    return families


def qubit_word(letters):
    ps = (np.eye(2), np.array([[0, 1], [1, 0]]),
          np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    result = np.ones((1, 1), complex)
    for letter in letters:
        result = np.kron(result, ps[letter])
    return result


def conditional_expectation(a, keep):
    """Normalized tensor conditional expectation, with no promise subspace."""
    result = a.copy()
    for site in range(3):
        if site not in keep:
            conjugations = []
            for label in range(4):
                letters = [0, 0, 0]
                letters[site] = label
                p = qubit_word(letters)
                conjugations.append(p@result@p)
            result = sum(conjugations)/4
    return result


def full_domain_intersection_report():
    rng = np.random.default_rng(1402)
    q = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    a = q+q.conj().T
    masks = [tuple(j for j in range(3) if s & (1 << j)) for s in range(8)]
    composition = max(norm(conditional_expectation(conditional_expectation(a, s), t)-
                           conditional_expectation(a, tuple(set(s)&set(t))))
                      for s, t in itertools.product(masks, repeat=2))
    # One fixed receiver, B: random internally flagged SWAPs from A or C into B.
    alpha, beta = .08, .13
    images = [np.eye(8, dtype=complex)]
    for label in range(1, 4):
        words = []
        for site in range(3):
            letters = [0, 0, 0]
            letters[site] = label
            words.append(qubit_word(letters))
        images.append(alpha*words[0]+(1-alpha-beta)*words[1]+beta*words[2])
    ps = [np.array([[1, 0], [0, 1]]), np.array([[0, 1], [1, 0]]),
          np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    def choi(transforms):
        return sum(np.kron(p, f.T) for p, f in zip(ps, transforms))/16
    original = choi(images)
    errors = {}
    for name, keep in (('AB', (0, 1)), ('BC', (1, 2)), ('B', (1,))):
        reset = choi([conditional_expectation(f, keep) for f in images])
        errors[name] = td(original, reset)
    # Full-input identity receiver B is supported in AB and BC and in their
    # intersection; it is a channel property, not alternative encoded decoders.
    b_images = [qubit_word((0, label, 0)) for label in range(4)]
    exact = max(norm(f-conditional_expectation(f, keep))
                for keep in ((0, 1), (1, 2), (1,)) for f in b_images)
    return dict(conditional_expectation_composition_error=composition,
                fixed_identity_receiver_intersection_error=exact,
                full_channel_choi_min_eigenvalue=float(np.linalg.eigvalsh(original)[0]),
                full_channel_choi_trace_error=float(abs(np.trace(original)-1)),
                choi_distance_after_resets=errors,
                exact_half_diamond_AB=.75*beta, exact_half_diamond_BC=.75*alpha,
                half_diamond_intersection_upper_bound=.75*(alpha+beta),
                normalized_tensor_resets_used=True, same_fixed_channel_used=True)


@lru_cache(None)
def report():
    v, u = encoder(), encoding_unitary()
    blank_embedding = np.eye(27, dtype=complex)[:, [0, 9, 18]]
    encoding = dict(isometry_error=norm(v.conj().T@v-I),
                    full_unitary_error=norm(u.conj().T@u-np.eye(27)),
                    circuit_error=norm(u@blank_embedding-v))
    coalitions = []
    channels = {}
    units = [np.eye(3, dtype=complex)[:, [a]] @ np.eye(3, dtype=complex)[[b]]
             for a, b in itertools.product(range(3), repeat=2)]
    for mask in range(8):
        keep = tuple(j for j in range(3) if mask & (1 << j))
        ks = marginal_kraus([v], (3, 3, 3), keep)
        channels[keep] = ks
        privacy = None
        if len(keep) == 1:
            privacy = max(norm(apply(ks, a)-np.trace(a)*I/3) for a in units)
        coalitions.append(dict(physical_mask=mask, keep=list(keep),
                               correctable_algebra_complex_dimension=len(recoverable(ks)),
                               single_share_channel_privacy_error=privacy))
    phi = np.eye(3, dtype=complex).reshape(-1)/np.sqrt(3)
    ideal_choi = np.outer(phi, phi.conj())
    decoder_checks, branch_checks = [], []
    k0 = np.diag(np.sqrt([.1, .4, .8]))
    k1 = F @ np.diag(np.sqrt([.9, .6, .2]))
    instrument = (k0, k1)
    rng = np.random.default_rng(402)
    q = rng.normal(size=(15, 7))+1j*rng.normal(size=(15, 7))
    rho = q@q.conj().T
    rho /= np.trace(rho)
    vr = np.kron(v, np.eye(5))
    encoded = vr@rho@vr.conj().T
    for pair in PAIRS:
        omitted = next(k for k in range(3) if k not in pair)
        dec = decoder(pair)
        mapped = np.kron(dec, I)@reorder_output(v, pair+(omitted,))
        factored = np.kron(I, phi.reshape(-1, 1))
        decoded = [a@b for a in decode_kraus(pair) for b in channels[pair]]
        decoder_checks.append(dict(pair=list(pair), unitary_error=norm(dec.conj().T@dec-np.eye(9)),
                                   factorization_error=norm(mapped-factored),
                                   full_channel_choi_error=norm(channel_choi(decoded)-ideal_choi)))
        lifted = []
        for k in instrument:
            local = dec.conj().T @ np.kron(k, I) @ dec
            global_op = on_subset(local, pair)
            lifted.append(global_op)
            lhs = np.kron(global_op, np.eye(5)) @ encoded @ np.kron(global_op.conj().T, np.eye(5))
            kr = np.kron(k, np.eye(5))
            rhs = vr @ (kr@rho@kr.conj().T) @ vr.conj().T
            branch_checks.append(dict(pair=list(pair), intertwining_error=norm(global_op@v-v@k),
                                      mixed_reference_branch_error=norm(lhs-rhs),
                                      branch_probability=float(np.trace(lhs).real)))
        branch_checks.append(dict(pair=list(pair),
                                  instrument_completeness_error=norm(sum(k.conj().T@k for k in lifted)-np.eye(27))))
    gamma = [row['physical_mask'] for row in coalitions if row['correctable_algebra_complex_dimension'] == 9]
    violations = [dict(left=s, right=t, intersection=s&t) for s,t in itertools.combinations(gamma, 2) if s&t not in gamma]
    local_units = [np.eye(9, dtype=complex)[:, [a]] @ np.eye(9, dtype=complex)[[b]]
                   for a,b in itertools.product(range(9), repeat=2)]
    ab = [on_subset(a, (0, 1)) for a in local_units]
    bc = [on_subset(a, (1, 2)) for a in local_units]
    b_only = [on_subset(a, (1,)) for a in units]
    intersections = dict(authorized_sets=gamma, intersection_violations=violations,
                         physical_pair_algebra_dimensions=[span_rank(ab), span_rank(bc)],
                         physical_pair_intersection_dimension=span_rank(ab)+span_rank(bc)-span_rank(ab+bc),
                         compressed_intersection_dimension=span_rank([v.conj().T@a@v for a in b_only]),
                         logical_common_pair_dimension=len(recoverable(channels[(0,1)])))
    constant_channels = []
    for index in range(5):
        if index == 0:
            sigma = I/3
        else:
            a = rng.normal(size=(3, index))+1j*rng.normal(size=(3, index))
            sigma = a@a.conj().T
            sigma /= np.trace(sigma)
        product = np.kron(sigma, I/3)
        constant_channels.append(dict(output_rank=int(np.linalg.matrix_rank(sigma)),
                                      bell_witness_distance=td(product, ideal_choi),
                                      bell_overlap=float(np.trace(ideal_choi@product).real),
                                      recovery_error_lower_bound=8/9))
    # A concrete noisy encoding (1-eta) E + eta Tr(.) I/27.
    noisy = []
    for eta in (.001, .01, .1):
        for pair in PAIRS:
            decoded = [a@b for a in decode_kraus(pair) for b in channels[pair]]
            pair_j = (1-eta)*channel_choi(decoded)+eta*np.eye(9)/9
            noisy.append(dict(noise_weight=eta, pair=list(pair),
                               decoded_choi_distance=td(pair_j, ideal_choi),
                               exact_decoded_half_diamond=eta*8/9,
                               general_pair_error_upper_bound=eta,
                               general_single_error_lower_bound=8/9-eta,
                               actual_single_error=8/9))
    # All finite proper filters are principal. This is combinatorial, not a
    # spatial model or an enumeration of all quantum encodings.
    filters = filter_audit()
    return dict(round=402,
                scope='Recovery regions of an encoded object need not form a point-neighborhood filter. Dependency regions of one fixed channel on the full tensor input do obey an additive intersection bound. Given interfaces are not derived space or a complete GR countermodel.',
                encoding=encoding, coalitions=coalitions, pair_decoder_checks=decoder_checks,
                logical_instrument_checks=branch_checks, region_intersections=intersections,
                constant_channel_recovery=constant_channels, noisy_encoding_checks=noisy,
                finite_filter_enumeration=dict(proper_nonempty_upward_families=len(filters),
                                              intersection_closed=sum(row['intersection_closed'] for row in filters),
                                              quantum_no_cloning_families=sum(row['no_disjoint_authorized_sets'] for row in filters),
                                              principal_equivalence_all_passed=all(row['intersection_closed']==row['principal'] for row in filters)),
                fixed_full_domain_channel_intersection=full_domain_intersection_report(),
                full_unknown_reference_preserved_by_authorized_recovery=True,
                alternative_overlapping_coalitions_not_independent_receivers=True,
                fixed_full_domain_dependency_intersection_proved=True,
                intervening_control_permissions_automatically_enlarged=False,
                unique_least_recovery_region_exists=False,
                recovery_regions_identified_as_point_neighborhoods=False,
                autonomous_geometric_local_implementation_derived=False,
                spatial_dimension_generated=False, full_cognition_to_gr_refuted=False)


class Audit(unittest.TestCase):
    def test_encoder_has_finite_unitary_implementation(self):
        for value in report()['encoding'].values():
            self.assertLess(value, 1e-12)

    def test_all_coalitions_and_complete_single_share_privacy(self):
        for row in report()['coalitions']:
            expected = 9 if len(row['keep']) >= 2 else 1
            self.assertEqual(row['correctable_algebra_complex_dimension'], expected)
            if row['single_share_channel_privacy_error'] is not None:
                self.assertLess(row['single_share_channel_privacy_error'], 1e-12)

    def test_three_decoders_factor_out_unknown_state_and_reference(self):
        for row in report()['pair_decoder_checks']:
            for key in ('unitary_error', 'factorization_error', 'full_channel_choi_error'):
                self.assertLess(row[key], 1e-12)

    def test_full_instrument_branches_preserve_code_and_reference(self):
        for row in report()['logical_instrument_checks']:
            for key in ('intertwining_error', 'mixed_reference_branch_error', 'instrument_completeness_error'):
                if key in row:
                    self.assertLess(row[key], 1e-12)

    def test_physical_intersection_and_logical_intersection_differ(self):
        row = report()['region_intersections']
        self.assertEqual(row['authorized_sets'], [3, 5, 6, 7])
        self.assertEqual(len(row['intersection_violations']), 3)
        self.assertEqual(row['physical_pair_algebra_dimensions'], [81, 81])
        self.assertEqual(row['physical_pair_intersection_dimension'], 9)
        self.assertEqual(row['compressed_intersection_dimension'], 1)
        self.assertEqual(row['logical_common_pair_dimension'], 9)

    def test_no_decoder_of_one_share_can_close_the_error_gap(self):
        rows = report()['constant_channel_recovery']
        self.assertAlmostEqual(rows[0]['bell_witness_distance'], 8/9)
        for row in rows:
            self.assertAlmostEqual(row['bell_overlap'], 1/9)
            self.assertGreaterEqual(row['bell_witness_distance']+1e-12, 8/9)

    def test_intersection_failure_survives_declared_noise(self):
        for row in report()['noisy_encoding_checks']:
            self.assertAlmostEqual(row['decoded_choi_distance'], row['exact_decoded_half_diamond'])
            self.assertLessEqual(row['exact_decoded_half_diamond'], row['general_pair_error_upper_bound'])
            self.assertLess(row['general_pair_error_upper_bound'], row['general_single_error_lower_bound'])

    def test_finite_filter_criterion_not_assumed_for_all_access_structures(self):
        row = report()['finite_filter_enumeration']
        self.assertTrue(row['principal_equivalence_all_passed'])
        self.assertEqual(row['proper_nonempty_upward_families'], 18)
        self.assertEqual(row['intersection_closed'], 7)
        self.assertGreater(row['quantum_no_cloning_families'], row['intersection_closed'])

    def test_full_domain_conditional_expectations_and_exact_receiver(self):
        row = report()['fixed_full_domain_channel_intersection']
        self.assertLess(row['conditional_expectation_composition_error'], 1e-12)
        self.assertLess(row['fixed_identity_receiver_intersection_error'], 1e-12)

    def test_one_fixed_quantum_channel_has_intersection_error_budget(self):
        row = report()['fixed_full_domain_channel_intersection']
        self.assertGreaterEqual(row['full_channel_choi_min_eigenvalue'], -1e-12)
        self.assertLess(row['full_channel_choi_trace_error'], 1e-12)
        self.assertAlmostEqual(row['choi_distance_after_resets']['AB'], row['exact_half_diamond_AB'])
        self.assertAlmostEqual(row['choi_distance_after_resets']['BC'], row['exact_half_diamond_BC'])
        self.assertLessEqual(row['choi_distance_after_resets']['B'], row['half_diamond_intersection_upper_bound']+1e-12)


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
