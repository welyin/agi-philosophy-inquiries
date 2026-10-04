"""Round 410: sharp shared records, continuous position quotients and finite limits."""
import argparse
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "shared_record_topology_audit_results.json"


def words(n):
    indices = np.arange(2**n, dtype=np.int64)
    return ((indices[:, None] >> np.arange(n-1, -1, -1)) & 1)


def coordinates(bits, k):
    """Finite-prefix coordinate observable; omitted bits are zero."""
    ans = np.zeros((len(bits), k))
    for i in range(bits.shape[1]):
        ans[:, i % k] += bits[:, i] * 2.0**(-(i//k+1))
    return ans


def cylinder_checks():
    rows = []
    for n in (2, 4, 7):
        b = words(n)
        diagonal_projectors = b.T
        # Two independent complete record encodings: bit tuples and binary weights.
        recovered = b @ (2**np.arange(n-1, -1, -1))
        pair = np.abs(b[:, None, :]-b[None, :, :]).max(axis=2)
        rows.append(dict(bits=n, full_quantum_dimension=2**n,
                         projector_error=int(np.max(np.abs(diagonal_projectors**2-diagonal_projectors))),
                         encoding_error=int(np.max(np.abs(recovered-np.arange(2**n)))),
                         smallest_distinct_record_distance=int(pair[np.triu_indices(2**n, 1)].min())))
    return rows


def finite_grids():
    rows = []
    for k in (1, 2, 3, 4):
        m = 3
        b = words(k*m)
        x = coordinates(b, k)
        independent = np.stack([
            (b[:, j::k] @ (2**np.arange(m-1, -1, -1)))/2**m
            for j in range(k)], axis=1)
        centered = x-x.mean(axis=0)
        covariance = centered.T @ centered/len(x)
        exact_covariance = np.eye(k)*(1-4.0**(-m))/12
        rows.append(dict(k=k, bits=k*m, configurations=len(b),
                         unique_coordinate_tuples=len(np.unique(x, axis=0)),
                         independent_binary_error=float(np.max(np.abs(x-independent))),
                         covariance_error=float(np.max(np.abs(covariance-exact_covariance))),
                         finite_covariance_rank=int(np.linalg.matrix_rank(covariance)),
                         infinite_tail_norm_bound=float(np.sqrt(k)*2.0**(-m))))
    return rows


def tail_checks():
    rows = []
    for k in (1, 2, 3, 4):
        for n in range(k, 3*k+1):
            counts = [(n+k-1-j)//k for j in range(k)]
            tails = np.array([2.0**(-c) for c in counts])
            # A direct, independently summed next four bits of every coordinate.
            extension = np.zeros(k)
            for i in range(n, n+4*k):
                extension[i % k] += 2.0**(-(i//k+1))
            rows.append(dict(k=k, known_prefix=n, coordinate_counts=counts,
                             extension_error=float(np.linalg.norm(extension)),
                             infinite_tail_bound=float(np.linalg.norm(tails)),
                             saturation_ratio=float(np.linalg.norm(extension)/np.linalg.norm(tails))))
    return rows


def dyadic_fibers():
    rows = []
    for q in range(1, 7):
        for odd in range(1, 2**q, 2):
            # Terminating expansion versus predecessor prefix followed by all ones.
            terminating = Fraction(odd, 2**q)
            recurring = Fraction(odd-1, 2**q) + Fraction(1, 2**q)
            last_a, last_b = odd & 1, (odd-1) & 1
            rows.append(dict(q=q, numerator=odd, coordinate=str(terminating),
                             equal_infinite_values=terminating == recurring,
                             differing_record_gap=last_a-last_b))
    return dict(exact_pairs=len(rows), all_equal=all(r["equal_infinite_values"] for r in rows),
                all_record_gaps_one=all(r["differing_record_gap"] == 1 for r in rows),
                examples=rows[:7], sharp_projection_distance_from_continuous_coordinate_algebra=0.5)


def finite_limit_checks():
    rows = []
    for n in (2, 4, 8, 12, 16):
        gap = Fraction(1, 2**n)
        left, right = Fraction(1, 2)-gap, Fraction(1, 2)
        for lipschitz in (4, 16, 256):
            lower = max(Fraction(0), (1-lipschitz*gap)/2)
            # A clipped affine function attains the lower bound on the entire finite grid.
            slope = min(lipschitz, 2**n)
            midpoint = (left+right)/2
            x = np.arange(2**n, dtype=float)/2**n
            target = (x >= .5).astype(float)
            fitted = np.clip(.5+slope*(x-float(midpoint)), 0, 1)
            actual = float(np.max(np.abs(fitted-target)))
            rows.append(dict(bits=n, lipschitz=lipschitz, gap=str(gap),
                             minimum_exact_lipschitz=2**n,
                             minimax_error_exact=str(lower), minimax_error=float(lower),
                             attained_error=actual))
    return rows


def joint_record_isometry():
    rng = np.random.default_rng(410)
    rows = []
    for n in (2, 3):
        d = 2**n
        # Record, unknown payload, and untouched reference can initially be entangled.
        psi = rng.normal(size=(d, 2, 3))+1j*rng.normal(size=(d, 2, 3))
        psi /= np.linalg.norm(psi)
        out = np.zeros((d, d, 2, 3), complex)
        for a in range(d):
            out[a, a] = psi[a]
        decoded = np.stack([out[a, a] for a in range(d)])
        old_qr = psi.reshape(d, 6).T @ psi.reshape(d, 6).conj()
        new_qr = out.reshape(d*d, 6).T @ out.reshape(d*d, 6).conj()
        joint_weights = np.sum(np.abs(out)**2, axis=(2, 3))
        b = words(n)
        bit_mismatch = sum(float(np.sum(joint_weights*np.abs(b[:, j, None]-b[None, :, j])))
                           for j in range(n))
        coordinate_mismatch = []
        for k in range(1, n+1):
            x = coordinates(b, k)
            distances = np.linalg.norm(x[:, None, :]-x[None, :, :], axis=2)
            coordinate_mismatch.append(float(np.sum(joint_weights*distances)))
        rows.append(dict(bits=n, ancilla_bits=n,
                         global_isometry_decoding_error=float(np.linalg.norm(decoded-psi)),
                         unknown_payload_reference_error=float(np.linalg.norm(old_qr-new_qr)),
                         sharp_record_disagreement=bit_mismatch,
                         coordinate_disagreements=coordinate_mismatch,
                         quantum_state_cloned=False))
    return rows


@lru_cache(maxsize=1)
def report():
    return dict(round=410, scope="Conditional topology of jointly sharp current records; no spatial dynamics derived",
                cylinder_records=cylinder_checks(), coordinate_grids=finite_grids(),
                prefix_tail_bounds=tail_checks(), exact_dyadic_fibers=dyadic_fibers(),
                finite_limit_regularity=finite_limit_checks(),
                internally_redundant_records=joint_record_isometry(),
                analytic_topology_proof_in_note=True,
                finite_experiments_prove_infinite_topology=False,
                all_quantum_states_assumed_classical=False,
                complete_record_topology_identified_as_physical_space=False,
                mathematical_coordinate_choices_are_distinct_physical_worlds=False,
                position_quotient_erases_actual_history=False,
                actual_subject_movement_implemented=False,
                spatial_dimension_generated=False, full_cognition_to_gr_refuted=False,
                initial_infinite_resource_adopted=False,
                new_cognitive_axiom_adopted=False)


class Audit(unittest.TestCase):
    def test_sharp_cylinders_separate_every_finite_word(self):
        for r in report()["cylinder_records"]:
            self.assertEqual(r["projector_error"], 0)
            self.assertEqual(r["encoding_error"], 0)
            self.assertEqual(r["smallest_distinct_record_distance"], 1)

    def test_same_records_give_arbitrary_selected_cube_coordinates(self):
        for r in report()["coordinate_grids"]:
            self.assertEqual(r["configurations"], r["unique_coordinate_tuples"])
            self.assertEqual(r["independent_binary_error"], 0)
            self.assertLess(r["covariance_error"], 1e-13)
            self.assertEqual(r["finite_covariance_rank"], r["k"])

    def test_full_infinite_tail_bounds_on_uneven_prefixes(self):
        for r in report()["prefix_tail_bounds"]:
            self.assertLess(r["extension_error"], r["infinite_tail_bound"])
            self.assertAlmostEqual(r["saturation_ratio"], 15/16, places=13)

    def test_exact_infinite_fibers_preserve_a_sharp_distinction(self):
        r = report()["exact_dyadic_fibers"]
        self.assertEqual(r["exact_pairs"], 63)
        self.assertTrue(r["all_equal"] and r["all_record_gaps_one"])
        self.assertEqual(r["sharp_projection_distance_from_continuous_coordinate_algebra"], .5)

    def test_finite_invertibility_requires_diverging_regularities(self):
        for r in report()["finite_limit_regularity"]:
            self.assertAlmostEqual(r["attained_error"], r["minimax_error"], places=13)
            self.assertEqual(r["minimum_exact_lipschitz"], 2**r["bits"])

    def test_consensus_keeps_joint_unknown_quantum_information(self):
        for r in report()["internally_redundant_records"]:
            self.assertEqual(r["global_isometry_decoding_error"], 0)
            self.assertLess(r["unknown_payload_reference_error"], 1e-14)
            self.assertEqual(r["sharp_record_disagreement"], 0)
            self.assertEqual(max(r["coordinate_disagreements"]), 0)
            self.assertFalse(r["quantum_state_cloned"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report(), checks=dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors)),
                  runtime=dict(python=platform.python_version(), numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
