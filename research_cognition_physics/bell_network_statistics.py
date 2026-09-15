"""Round 66: one fixed two-source Bell network and its complete statistics.

Matrix probabilities and source tensor independence are declared assumptions.
The score uses all Bob outcomes, without conditional renormalization.
"""

import argparse
import json
import math
import unittest
from itertools import combinations, product
from pathlib import Path

import numpy as np

from bipartite_composition import interaction
from independent_source_alignment import keep_systems
from quantum_interface_audit import IDENTITY, PAULI_X, PAULI_Y, PAULI_Z, ALPHA, effective_readout


AXES = (PAULI_X, PAULI_Y, PAULI_Z)
PAIRS = tuple(combinations(range(3), 2))
SETTINGS = tuple((i, j, sign) for i, j in PAIRS for sign in (1, -1))
BELL_SIGNS = ((1, -1, 1), (1, 1, -1), (-1, -1, -1), (-1, 1, 1))


def bell_basis():
    return interaction(math.pi/2, "YX")


def bell_projectors():
    basis = bell_basis()
    return tuple(np.outer(basis[:, b], basis[:, b].conj()) for b in range(4))


def independent_bell_sources(visibility_left=1., visibility_right=1.):
    phi = bell_projectors()[0]
    return np.kron(visibility_left*phi+(1-visibility_left)*np.eye(4)/4,
                   visibility_right*phi+(1-visibility_right)*np.eye(4)/4)


def bob_effects(pointer_visibility=1.):
    if not 0 <= pointer_visibility <= 1:
        raise ValueError("Use pointer visibility in [0,1].")
    basis = bell_basis()
    effects = []
    for b in range(4):
        bits = (b//2, b%2)
        effects.append(basis @ np.kron(*[(IDENTITY+(-1)**bit*pointer_visibility*PAULI_Z)/2 for bit in bits]) @ basis.conj().T)
    return tuple(effects)


def endpoint_observables(visibility=1., imaginary_bias=1.):
    alice = tuple(visibility*(imaginary_bias if i == 1 else 1)*axis for i, axis in enumerate(AXES))
    charlie = tuple((alice[i]+sign*alice[j])/math.sqrt(2) for i, j, sign in SETTINGS)
    return alice, charlie


def network_table(alice=None, charlie=None, bob=None, source=None):
    if alice is None:
        alice = endpoint_observables()[0]
    if charlie is None:
        charlie = endpoint_observables()[1]
    if bob is None:
        bob = bob_effects()
    if source is None:
        source = independent_bell_sources()
    return {(x, z, a, b, c): float(np.trace(source @ np.kron(np.kron((IDENTITY+a*alice[x])/2, bob[b]), (IDENTITY+c*charlie[z])/2)).real)
            for x, z, a, b, c in product(range(3), range(6), (-1, 1), range(4), (-1, 1))}


def network_score(records):
    def correlator(x, z, b):
        return sum(a*c*records[x, z, a, b, c] for a, c in product((-1, 1), repeat=2))
    return sum(BELL_SIGNS[b][i]*(correlator(i, 2*p, b)+correlator(i, 2*p+1, b))
               +BELL_SIGNS[b][j]*(correlator(j, 2*p, b)-correlator(j, 2*p+1, b))
               for b in range(4) for p, (i, j) in enumerate(PAIRS))


def diagonal_score(visibility_a, visibility_c, bob_visibility, bias_a=1., bias_c=1.):
    return 2*math.sqrt(2)*visibility_a*visibility_c*(bob_visibility*(1+bias_a*bias_c)+bob_visibility**2)


class BellNetworkStatisticsTests(unittest.TestCase):
    def test_old_real_gate_gives_the_four_bell_sign_patterns(self):
        for projector, signs in zip(bell_projectors(), BELL_SIGNS):
            self.assertEqual(math.prod(signs), -1)
            for axis, sign in zip(AXES, signs):
                self.assertAlmostEqual(np.trace(projector @ np.kron(axis, axis)).real, sign, places=14)

    def test_two_sources_are_independent_and_outer_marginal_is_maximally_mixed(self):
        source = independent_bell_sources()
        np.testing.assert_allclose(keep_systems(source, (0, 1), 4), bell_projectors()[0], atol=3e-16)
        np.testing.assert_allclose(keep_systems(source, (0, 3), 4), np.eye(4)/4, atol=2e-16)

    def test_complete_probability_table_matches_the_analytic_swapping_formula(self):
        records = network_table()
        for (x, z, a, b, c), probability in records.items():
            i, j, sign = SETTINGS[z]
            correlation = BELL_SIGNS[b][x]*((x == i)+sign*(x == j))/math.sqrt(2)
            self.assertAlmostEqual(probability, (1+a*c*correlation)/16, places=14)

    def test_probabilities_are_positive_normalized_and_keep_all_bob_outputs(self):
        records = network_table()
        self.assertGreaterEqual(min(records.values()), 0.)
        for x, z in product(range(3), range(6)):
            self.assertAlmostEqual(sum(records[x, z, a, b, c] for a, b, c in product((-1, 1), range(4), (-1, 1))), 1., places=14)
            for b in range(4):
                self.assertAlmostEqual(sum(records[x, z, a, b, c] for a, c in product((-1, 1), repeat=2)), .25, places=14)

    def test_score_is_six_sqrt_two_without_conditional_rescaling(self):
        self.assertAlmostEqual(network_score(network_table()), 6*math.sqrt(2), places=13)

    def test_noisy_bell_decoder_is_the_actual_old_inverse_gate_plus_z_effects(self):
        basis = bell_basis()
        for b, effect in enumerate(bob_effects(ALPHA)):
            direct = basis @ np.kron(effective_readout(0, 1-b//2, 1), effective_readout(0, 1-b%2, 1)) @ basis.conj().T
            np.testing.assert_allclose(effect, direct, atol=3e-16)
        np.testing.assert_allclose(sum(bob_effects(.7)), np.eye(4), atol=3e-16)

    def test_anisotropic_endpoint_noise_and_two_bit_bob_noise_match_the_score_formula(self):
        for a, c, b, nu_a, nu_c in ((.7, .8, .6, .3, .4), (ALPHA, ALPHA, ALPHA, 1., 1.), (1., 1., 1., .2, .8)):
            alice = endpoint_observables(a, nu_a)[0]
            charlie = endpoint_observables(c, nu_c)[1]
            actual = network_score(network_table(alice, charlie, bob_effects(b)))
            self.assertAlmostEqual(actual, diagonal_score(a, c, b, nu_a, nu_c), places=13)

    def test_depolarized_independent_sources_multiply_the_score_by_both_visibilities(self):
        actual = network_score(network_table(source=independent_bell_sources(.8, .7)))
        self.assertAlmostEqual(actual, .56*6*math.sqrt(2), places=13)

    def test_ignored_bob_output_gives_uniform_outer_records_for_all_settings(self):
        records = network_table()
        for x, z, a, c in product(range(3), range(6), (-1, 1), (-1, 1)):
            self.assertAlmostEqual(sum(records[x, z, a, b, c] for b in range(4)), .25, places=14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BellNetworkStatisticsTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    records = network_table()
    report = {"round": 66, "source_count": 2, "settings": "Alice X,Y,Z; Charlie six pairwise sum/difference directions", "bob_signs": BELL_SIGNS,
              "probability_entries": len(records), "ideal_score": network_score(records),
              "one_raw_read_each_score": diagonal_score(ALPHA, ALPHA, ALPHA),
              "all_bob_outputs_retained": True, "extra_shared_target_source": False,
              "records": [{"x":x,"z":z,"a":a,"b":b,"c":c,"probability":p} for (x,z,a,b,c),p in records.items()],
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("bell_network_statistics_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
