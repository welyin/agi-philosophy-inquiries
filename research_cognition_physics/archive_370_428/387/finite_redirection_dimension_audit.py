"""Round 387: a finite actual redirection group gives a dimension lower bound.

The old tetrahedral matrices are reused. The new claim concerns their action
on a full position sphere, not the dimension of a finite program register.
"""
import unittest

import numpy as np

from direction_contrast_tomography_audit import tetrahedral_group, PAULI, rho
from growing_stream_audit import main


I2 = np.eye(2, dtype=complex)
RA = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], int)
RB = np.diag([1, -1, -1])
UA = (I2 - 1j * sum(PAULI)) / 2
UB = -1j * PAULI[0]


def key(matrix):
    return tuple(np.asarray(matrix, int).ravel())


GROUP = tetrahedral_group()
INDEX = {key(g): i for i, g in enumerate(GROUP)}
IDENTITY = INDEX[key(np.eye(3, dtype=int))]
IA, IB = INDEX[key(RA)], INDEX[key(RB)]
TABLE = np.array([[INDEX[key(g @ h)] for h in GROUP] for g in GROUP])


def normal_forms():
    result = {}
    v4 = [np.eye(3, dtype=int), RB, RA @ RB @ RA.T,
          RA @ RA @ RB @ RA.T @ RA.T]
    lifts = [I2, UB, UA @ UB @ UA.conj().T,
             UA @ UA @ UB @ UA.conj().T @ UA.conj().T]
    for v, w in zip(v4, lifts):
        for j in range(3):
            g = v @ np.linalg.matrix_power(RA, j)
            result[INDEX[key(g)]] = (j, w @ np.linalg.matrix_power(UA, j))
    return result


FORMS = normal_forms()
QUOTIENT = np.array([FORMS[g][0] for g in range(12)])
LIFTS = [FORMS[g][1] for g in range(12)]


def permutation(h):
    p = np.zeros((12, 12))
    for g in range(12):
        p[TABLE[h, g], g] = 1
    return p


def position_probability(n, visibility=0.6, probe=np.array([0., 0., 1.])):
    return float((1 + visibility * np.dot(n, probe)) / 2)


def position_gap_lower(observed_first, observed_second, error):
    return float(abs(observed_first - observed_second) - 2 * error)


def circle_action(g, theta):
    return theta + 2 * np.pi * QUOTIENT[g] / 3


def chord(theta, phi):
    return float(abs(np.exp(1j * theta) - np.exp(1j * phi)))


def cyclic_presentations(modulus):
    return [(a, b) for a in range(modulus) for b in range(modulus)
            if 3 * a % modulus == 0 and 2 * b % modulus == 0
            and 3 * (a + b) % modulus == 0]


def controlled(unitary):
    return np.block([[I2, np.zeros((2, 2))], [np.zeros((2, 2)), unitary]])


def report():
    n = np.array([0., 0., 1.])
    p0, p1 = position_probability(n), position_probability(RB @ n)
    negatives = 0
    for g in range(12):
        for h in range(12):
            overlap = np.trace(LIFTS[TABLE[g, h]].conj().T @ LIFTS[g] @ LIFTS[h]) / 2
            negatives += int(overlap.real < 0)
    return {
        "round": 387,
        "scope": {
            "lower_bound_theorem": "A full L homeomorphic to S^(n-1), n>=1, carrying continuous endpoint maps a,b with a^3=b^2=(ab)^3=id and b not identity has n>=3. Finite order already makes the maps homeomorphisms. No qubit covariance, contrast tomography, minimality, or continuous SO3 control is needed.",
            "contract_change": "This is an alternative sufficient contract, not automatically weaker than round 384: it adds exact global finite-word relations on actual position identities. Dense approximate rotations need not supply them.",
            "upper_bound": "Dimension exactly three still needs the separately certified inverse separation of round 383 or 386 on the same full shell.",
            "position_witness": "A single probability difference can witness b nonidentity only for the same specified position readout under matched internal/controller preparation. Reading a program label is not a position witness.",
            "finite_program": "The known round-383 conjugation construction can use 12 finite labels and a two-dimensional payload. This is only recalled analytically, not rerun or counted as a new program construction. It does not establish a continuous position-shell action, a cost in physical time, or an autonomous natural Hamiltonian.",
            "historical_reuse": "The 12 integer tetrahedral matrices come from round 380. Their old closure and four-point-orbit results are not new results of this round.",
            "not_claimed": ["Exact global endpoint identities follow from finitely many noisy samples", "Internal qubit gate availability implies actual spatial redirection", "An SU2 channel identity is automatically an identity of a phase-sensitive source", "The finite group chooses three-dimensional geometry without an upper bound", "All extra contracts have been derived from FUCP"],
        },
        "finite_presentation": {"relations": "a^3=b^2=(ab)^3=e",
                                "normal_form": "v*a^j, v in {e,b,aba^-1,a^2ba^-2}, j=0,1,2",
                                "normal_form_count": len(FORMS),
                                "tested_cyclic_targets": list(range(2, 16)),
                                "all_tested_cyclic_images_kill_b": all(b == 0 for m in range(2, 16) for _, b in cyclic_presentations(m))},
        "single_position_probe": {"visibility": 0.6, "port": n.tolist(),
                                  "before": p0, "after_b": p1,
                                  "probability_error_per_reading": 0.01,
                                  "guaranteed_nonzero_gap_from_exact_centers": position_gap_lower(p0, p1, 0.01)},
        "circle_counterexamples": {"A4_to_C3_a_motion_chord": chord(circle_action(IA, 0), 0),
                                   "A4_to_C3_b_motion_chord": chord(circle_action(IB, 0), 0),
                                   "dihedral_a3_b2_hold_but_ab3_endpoint_defect": chord(2 * np.pi / 3, 0)},
        "sample_relations_are_not_global_relations": {
            "circle_samples": [float(np.pi / 3), float(4 * np.pi / 3)],
            "all_three_relations_hold_at_samples": True,
            "b_readable_probability_gap_at_first_sample": float(0.6 * np.sin(np.pi / 3)),
            "mixed_relation_defect_at_zero": chord(2 * np.pi / 3, 0)},
        "approximate_relations_need_a_relative_margin": [
            {"b_rotation_angle": t, "a3_defect": 0., "b2_defect": chord(2 * t, 0),
             "ab3_defect": chord(3 * t, 0), "b_motion": chord(t, 0),
             "one_probe_probability_gap": float(0.3 * np.sin(t))}
            for t in (0.01, 0.001)],
        "qubit_lift_audit": {"A_cubed": "-I", "B_squared": "-I", "AB_cubed": "+I",
                             "negative_lift_product_signs_in_144_pairs": negatives,
                             "controlled_B_squared_changes_plus_to_minus": True},
        "program_only_false_witness": {"real_position_action": "identity for every group label",
                                       "memory_readout_probability_change": 1.,
                                       "position_readout_probability_change": 0.},
        "primary_source": "https://www.math.uchicago.edu/~farb/papers/mcgbook.pdf",
    }


class Checks(unittest.TestCase):
    def test_01_presentation_normal_forms_and_conjugate_involutions(self):
        a, b = RA, RB
        c, d = a @ b @ a.T, a @ a @ b @ a.T @ a.T
        np.testing.assert_array_equal(np.linalg.matrix_power(a, 3), np.eye(3, dtype=int))
        np.testing.assert_array_equal(b @ b, np.eye(3, dtype=int))
        np.testing.assert_array_equal(np.linalg.matrix_power(a @ b, 3), np.eye(3, dtype=int))
        np.testing.assert_array_equal(b @ c, d)
        np.testing.assert_array_equal(c @ b, d)
        self.assertEqual(len(FORMS), 12)

    def test_02_any_tested_cyclic_quotient_kills_the_second_generator(self):
        for modulus in range(2, 16):
            solutions = cyclic_presentations(modulus)
            self.assertTrue(solutions)
            self.assertTrue(all(b == 0 for _, b in solutions))
        self.assertEqual(cyclic_presentations(2), [(0, 0)])
        self.assertIn((1, 0), cyclic_presentations(3))

    def test_03_nontrivial_circle_action_of_the_whole_group_is_not_enough(self):
        for g in range(12):
            for h in range(12):
                self.assertAlmostEqual(chord(circle_action(g, circle_action(h, 0.37)),
                                             circle_action(TABLE[g, h], 0.37)), 0, places=14)
        self.assertGreater(chord(circle_action(IA, 0.37), 0.37), 1.)
        self.assertEqual(circle_action(IB, 0.37), 0.37)

    def test_04_omitting_the_mixed_word_allows_a_visible_reflection_on_the_circle(self):
        a = lambda t: t + 2 * np.pi / 3
        b = lambda t: -t
        self.assertAlmostEqual(chord(a(a(a(0.))), 0.), 0., places=14)
        self.assertEqual(b(b(0.5)), 0.5)
        self.assertGreater(chord(b(0.5), 0.5), 0.9)
        word = lambda t: a(b(t))
        self.assertGreater(chord(word(word(word(0.))), 0.), 1.7)

    def test_05_one_fixed_position_probe_certifies_a_nonidentity_involution(self):
        n = np.array([0., 0., 1.])
        p0, p1 = position_probability(n), position_probability(RB @ n)
        self.assertAlmostEqual(p0, 0.8)
        self.assertAlmostEqual(p1, 0.2)
        self.assertAlmostEqual(position_gap_lower(p0, p1, 0.01), 0.58)
        for e0 in (-0.01, 0.01):
            for e1 in (-0.01, 0.01):
                self.assertLessEqual(position_gap_lower(p0 + e0, p1 + e1, 0.01), p0 - p1 + 1e-15)

    def test_06_qubit_channel_relations_do_not_erase_controlled_lift_phases(self):
        np.testing.assert_allclose(np.linalg.matrix_power(UA, 3), -I2, atol=2e-15)
        np.testing.assert_allclose(UB @ UB, -I2, atol=2e-15)
        np.testing.assert_allclose(np.linalg.matrix_power(UA @ UB, 3), I2, atol=2e-15)
        r = rho([0.2, -0.3, 0.4])
        np.testing.assert_allclose((UB @ UB) @ r @ (UB @ UB).conj().T, r)
        plus = np.array([1., 1.]) / np.sqrt(2)
        minus = np.array([1., -1.]) / np.sqrt(2)
        start = np.kron(plus, np.array([1., 0.]))
        output = controlled(UB) @ controlled(UB) @ start
        self.assertAlmostEqual(abs(np.vdot(np.kron(minus, [1., 0.]), output)), 1.)

    def test_07_exact_word_returns_at_sample_ports_are_not_exact_global_relations(self):
        a = lambda t: t + 2 * np.pi / 3
        b = lambda t: -t
        word = lambda t: a(b(t))
        for t in (np.pi / 3, 4 * np.pi / 3):
            self.assertAlmostEqual(chord(a(a(a(t))), t), 0., places=14)
            self.assertAlmostEqual(chord(b(b(t)), t), 0., places=14)
            self.assertAlmostEqual(chord(word(word(word(t))), t), 0., places=14)
            self.assertGreater(chord(b(t), t), 1.7)
        self.assertGreater(chord(word(word(word(0.))), 0.), 1.7)

    def test_08_approximate_global_relations_and_a_nonzero_b_still_allow_a_circle(self):
        for t in (0.01, 0.001):
            # a is identity, b is an everywhere small rotation.
            self.assertLessEqual(chord(2 * t, 0), 2 * t)
            self.assertLessEqual(chord(3 * t, 0), 3 * t)
            self.assertGreater(chord(t, 0), 0.)
            self.assertGreater(0.3 * np.sin(t), 0.)

    def test_09_exact_memory_group_and_a_readable_label_are_not_a_position_action(self):
        e = np.eye(12)[IDENTITY]
        final = permutation(IB) @ e
        self.assertAlmostEqual(abs(np.vdot(e, final)) ** 2, 0.)
        n = np.array([0., 0., 1.])
        # In this countermodel all physical endpoint maps are identity.
        self.assertEqual(position_probability(n), position_probability(n.copy()))


if __name__ == '__main__':
    main(__name__, 'finite_redirection_dimension_audit', report)
