"""Round 451: consistent pair comparison frames for the existing record terms.
Exact integer certificates; no spatial dimension or full SoCA implementation.
"""
import argparse
import io
import itertools as it
import json
from pathlib import Path
import platform
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "global_comparison_frames_audit_results.json"
OBS = {}


def sylvester(n):
    assert n > 0 and n & (n-1) == 0
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    return h


def paley12():
    q = 11
    squares = {x*x % q for x in range(1, q)}
    chi = lambda x: 0 if x % q == 0 else (1 if x % q in squares else -1)
    c = np.array([[chi(i-j) for j in range(q)] for i in range(q)], dtype=np.int64)
    h = np.ones((q+1, q+1), dtype=np.int64)
    h[1:, 1:] = -c-np.eye(q, dtype=np.int64)
    return c, h


def frame_numerators(h):
    """U_i = A_i/N. H first row is all +1."""
    return [(h*row[None, :])@h.T for row in h]


def pair_numerator(h, i, j):
    """F_ij = B_ij/N, including i=j."""
    return (h*(h[i]*h[j])[None, :])@h.T


def with_blank(a, denominator):
    result = np.zeros((len(a)+1, len(a)+1), dtype=np.int64)
    result[:-1, :-1] = a
    result[-1, -1] = denominator
    return result


def xor_frame(m, i):
    a = np.zeros((m, m), dtype=np.int64)
    for x in range(m):
        a[x ^ i, x] = 1
    return a


def permutation_families(n):
    """Exhaustive involutive dictionaries with U_i(0)=i and all U_i commuting."""
    identity = tuple(range(n))
    candidates = []
    for i in range(1, n):
        candidates.append([p for p in it.permutations(range(n))
                           if p[0] == i and all(p[p[x]] == x for x in range(n))])
    count = 0

    def visit(chosen):
        nonlocal count
        if len(chosen) == n:
            count += 1
            return
        for p in candidates[len(chosen)-1]:
            if all(all(p[q[x]] == q[p[x]] for x in range(n)) for q in chosen):
                visit(chosen+[p])
    visit([identity])
    return count, [len(v) for v in candidates]


class Audit(unittest.TestCase):
    def test_01_classical_dictionary_exhaustion(self):
        rows = []
        for n, expected in ((2, 1), (3, 0), (4, 1), (6, 0)):
            count, sizes = permutation_families(n)
            self.assertEqual(count, expected)
            rows.append(dict(subjects=n, candidate_counts=sizes, compatible_families=count))
        OBS["permutation_branch"] = dict(exhaustive_small_cases=rows,
            arbitrary_N_power_of_two_requirement_is_analytic=True,
            binary_address_length_is_not_spatial_dimension=True)

    def test_02_hadamard_and_all_pair_composition(self):
        rows = []
        for n in (2, 4, 8, 12):
            h = paley12()[1] if n == 12 else sylvester(n)
            eye = np.eye(n, dtype=np.int64)
            np.testing.assert_array_equal(h@h.T, n*eye)
            a = frame_numerators(h)
            np.testing.assert_array_equal(a[0], n*eye)
            for i in range(n):
                np.testing.assert_array_equal(a[i], a[i].T)
                np.testing.assert_array_equal(a[i]@a[i], n*n*eye)
                np.testing.assert_array_equal(a[i][:, 0], n*eye[:, i])
                for j in range(n):
                    b = pair_numerator(h, i, j)
                    np.testing.assert_array_equal(a[j]@a[i], n*b)
                    np.testing.assert_array_equal(a[i]@a[j], n*b)
                    np.testing.assert_array_equal(b@b, n*n*eye)
                    np.testing.assert_array_equal(b[:, j], n*eye[:, i])
                    for k in range(n):
                        np.testing.assert_array_equal(pair_numerator(h, j, k)@b,
                                                      n*pair_numerator(h, i, k))
            # H diag(w) H^T = I has w_alpha=1/N; check without floats.
            np.testing.assert_array_equal(np.sum(h*h, axis=0), n*np.ones(n, dtype=int))
            rows.append(dict(order=n, ordered_pair_identities=n*n, ordered_triple_identities=n**3,
                             uniform_joint_eigenbasis_weight=f"1/{n}"))
        OBS["coherent_branch"] = dict(exact_integer_certificates=rows,
            necessity_uses_square_N_by_N_matrix=True, existence_iff_real_Hadamard_proved=True)

    def test_03_full_local_record_and_packet_tensor_identities(self):
        rows = []
        for n in (4, 12):
            h = paley12()[1] if n == 12 else sylvester(n)
            d = n+1
            local_swap = np.arange(d*d).reshape(d, d).T.reshape(-1)
            eye = np.eye(d*d, dtype=np.int64)
            for i, j in it.combinations(range(n), 2):
                b = with_blank(pair_numerator(h, i, j), n)
                f = np.kron(b, b)  # numerator of hat F, denominator n^2
                np.testing.assert_array_equal(f@f, n**4*eye)
                np.testing.assert_array_equal(f[local_swap, :], f[:, local_swap])
                # These two entire columns certify F-hat R(j) F-hat*=R(i).
                for col, target in ((j*d+n, i*d+n), (n*d+j, n*d+i)):
                    np.testing.assert_array_equal(f[:, col], n*n*eye[:, target])
                # Active marker projection has an arbitrary second register.
                select_j = np.arange(j*d, (j+1)*d)
                select_i = np.arange(i*d, (i+1)*d)
                transformed = f[:, select_j]@f[:, select_j].T
                expected = np.zeros_like(transformed)
                expected[select_i, select_i] = n**4
                np.testing.assert_array_equal(transformed, expected)
            rows.append(dict(subjects=n, local_raw_dimension=d*d,
                             all_unordered_pairs=n*(n-1)//2,
                             tensor_factor_certificates_exact=True))
        OBS["raw_operator_scope"] = dict(rows=rows,
            packet_commutator_reduced_by_exact_tensor_identity=True,
            full_unknown_states_and_references_covered=True,
            huge_pair_matrix_or_random_state_sampling_used=False)

    def test_04_twelve_subject_coherent_counterexample_to_permutation_restriction(self):
        c, h = paley12()
        np.testing.assert_array_equal(c.T, -c)
        np.testing.assert_array_equal(c@c.T, 11*np.eye(11, dtype=int)-np.ones((11, 11), dtype=int))
        witness = None
        for i, a in enumerate(frame_numerators(h)):
            for j in range(12):
                support = np.flatnonzero(a[:, j])
                if len(support) > 1:
                    witness = dict(frame=i, input_marker=j, denominator=12,
                                   output_numerators=a[:, j].tolist(), support=len(support))
                    break
            if witness:
                break
        self.assertIsNotNone(witness)
        self.assertNotEqual(12 & (12-1), 0)
        OBS["paley_twelve"] = dict(character_modulus=11, hadamard=h.tolist(),
            coherent_column_witness=witness, permutation_dictionary_possible=False,
            coherent_dictionary_possible=True, all_multiples_of_four_existence_claimed=False)

    def test_05_six_subject_padding_and_unknown_information(self):
        n, m = 6, 8
        frames = [xor_frame(m, i) for i in range(n)]
        for i, j, k in it.product(range(n), repeat=3):
            b = frames[j]@frames[i].T
            np.testing.assert_array_equal(b@b, np.eye(m, dtype=int))
            np.testing.assert_array_equal(b[:, j], np.eye(m, dtype=int)[:, i])
            np.testing.assert_array_equal((frames[k]@frames[j].T)@b, frames[k]@frames[i].T)
        e = np.zeros((m+1, n+1), dtype=int)
        e[:n, :n] = np.eye(n, dtype=int)
        e[m, n] = 1  # old blank -> new blank, not reserve marker 6
        w = np.kron(e, e)
        np.testing.assert_array_equal(w.T@w, np.eye((n+1)**2, dtype=int))
        wr = np.kron(w, np.eye(3, dtype=int))
        np.testing.assert_array_equal(wr.T@wr, np.eye(3*(n+1)**2, dtype=int))
        # An old marker can leave the old alphabet under a legitimate comparison.
        self.assertEqual(int(np.argmax(frames[3][:, 5])), 6)
        OBS["padding"] = dict(subjects=n, old_marker_dimension=n, padded_marker_dimension=m,
            old_raw_subject_dimension=49, new_raw_subject_dimension=81, extra_marker_levels=2,
            all_ordered_triple_identities=216, old_unknown_input_embedding_isometric=True,
            reference_dimension_checked=3, arbitrary_reference_extension_is_analytic=True,
            old_subspace_closed_under_all_new_comparisons=False,
            witness=dict(frame=3, marker=5, becomes_reserve_marker=6),
            autonomous_joining_or_free_preparation_proved=False)

    def test_06_pair_contract_is_not_full_network_symmetry(self):
        n = 4
        words = list(it.product(range(n), repeat=n))
        def energy(w):
            return n-2*sum(w[i] == j and w[j] == i for i,j in it.combinations(range(n), 2))
        witness = None
        for w in words:
            v = list(w)
            v[0], v[1] = w[1] ^ 1, w[0] ^ 1
            if energy(w) != energy(v):
                witness = dict(markers_before=list(w), markers_after=v,
                               h_R_before=energy(w), h_R_after=energy(v))
                break
        self.assertIsNotNone(witness)
        OBS["spectator_boundary"] = dict(subjects=4, all_D_registers_blank=True,
            exact_basis_counterexample=witness,
            pairwise_compatible_frames_imply_full_network_Tij_symmetry=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=451, baseline_round=450, date="2026-09-24",
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(consistent_pair_frame_contract_classified=True,
            permutation_and_coherent_dictionaries_distinguished=True,
            minimal_marker_dimension_requirement_explicit=True,
            full_local_record_and_exchange_tensor_identities_verified=True,
            explicit_paley_twelve_and_padded_six_certificates=True,
            pair_contract_not_full_network_symmetry_verified=True,
            soca_functional_interface_motivation_distinguished_from_full_state_identification=True,
            comparison_frames_dynamically_generated=False,
            new_record_potentials_derived_from_429=False,
            full_soca_architecture_implemented=False,
            all_Hadamard_orders_classified=False,
            exact_global_identification_adopted_as_cognitive_axiom=False,
            full_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))

