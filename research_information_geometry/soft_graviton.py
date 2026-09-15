"""Check the conditional leading-soft spin-two universality argument."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SIGNATURE = np.array([1.0, -1.0, -1.0, -1.0])


@dataclass(frozen=True)
class HardEvent:
    momenta: np.ndarray
    species: tuple
    signs: np.ndarray


def minkowski_dot(first, second):
    return np.sum(np.asarray(first) * SIGNATURE * np.asarray(second), axis=-1)


def elastic_event(first_species="A", second_species="B", first_mass=0.6,
                  second_mass=1.1, momentum=0.8, angle=1.1, azimuth=0.4):
    parameters = (first_mass, second_mass, momentum, angle, azimuth)
    if not all(math.isfinite(value) for value in parameters):
        raise ValueError("Scattering parameters must be finite.")
    if min(first_mass, second_mass, momentum) <= 0:
        raise ValueError("This benchmark uses positive masses and momentum.")
    incoming = np.array([0.0, 0.0, momentum])
    outgoing = momentum * np.array([
        math.sin(angle) * math.cos(azimuth),
        math.sin(angle) * math.sin(azimuth),
        math.cos(angle),
    ])
    first_energy = math.hypot(first_mass, momentum)
    second_energy = math.hypot(second_mass, momentum)
    momenta = np.array([
        np.concatenate(([first_energy], incoming)),
        np.concatenate(([second_energy], -incoming)),
        np.concatenate(([first_energy], outgoing)),
        np.concatenate(([second_energy], -outgoing)),
    ])
    return HardEvent(momenta, (first_species, second_species, first_species, second_species),
                     np.array([-1.0, -1.0, 1.0, 1.0]))


def null_momentum(energy=0.02, direction=(0.3, -0.4, 1.0)):
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("A finite three-dimensional direction is required.")
    if not math.isfinite(energy) or energy <= 0 or np.linalg.norm(direction) == 0:
        raise ValueError("The soft energy and direction norm must be positive.")
    return energy * np.concatenate(([1.0], direction / np.linalg.norm(direction)))


def transverse_basis(soft_momentum):
    direction = np.asarray(soft_momentum)[1:] / soft_momentum[0]
    if not math.isclose(float(np.linalg.norm(direction)), 1.0, abs_tol=1e-12):
        raise ValueError("The soft four-momentum must be null.")
    reference_axis = np.eye(3)[int(np.argmin(abs(direction)))]
    first = reference_axis - np.dot(reference_axis, direction) * direction
    first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    return np.concatenate(([0.0], first)), np.concatenate(([0.0], second))


def pure_gauge_polarization(soft_momentum, gauge_vector):
    soft_covector = SIGNATURE * soft_momentum
    gauge_covector = SIGNATURE * gauge_vector
    return np.outer(soft_covector, gauge_covector) + np.outer(gauge_covector, soft_covector)


def plus_polarization(soft_momentum):
    first, second = transverse_basis(soft_momentum)
    first_covector = SIGNATURE * first
    second_covector = SIGNATURE * second
    return (np.outer(first_covector, first_covector)
            - np.outer(second_covector, second_covector)) / math.sqrt(2)


def leg_couplings(event, couplings):
    return np.array([couplings[species] for species in event.species], dtype=float)


def soft_tensor(event, couplings, soft_momentum):
    denominators = minkowski_dot(event.momenta, soft_momentum)
    if np.any(denominators <= 0):
        raise ValueError("Future-directed massive legs require positive soft denominators.")
    weights = event.signs * leg_couplings(event, couplings) / denominators
    return np.einsum("l,li,lj->ij", weights, event.momenta, event.momenta)


def soft_factor(event, couplings, soft_momentum, polarization):
    return float(np.sum(soft_tensor(event, couplings, soft_momentum) * polarization))


def ward_residual(event, couplings):
    weights = event.signs * leg_couplings(event, couplings)
    return np.sum(weights[:, None] * event.momenta, axis=0)


def photon_soft_current(event, charges, soft_momentum):
    weights = event.signs * leg_couplings(event, charges) / minkowski_dot(event.momenta, soft_momentum)
    return np.sum(weights[:, None] * event.momenta, axis=0)


def lorentz_boost(velocity):
    velocity = np.asarray(velocity, dtype=float)
    if velocity.shape != (3,) or not np.all(np.isfinite(velocity)):
        raise ValueError("A finite three-dimensional boost velocity is required.")
    speed_squared = float(np.dot(velocity, velocity))
    if speed_squared >= 1:
        raise ValueError("The boost velocity must be subluminal.")
    transformation = np.eye(4)
    if speed_squared == 0:
        return transformation
    gamma = 1 / math.sqrt(1 - speed_squared)
    transformation[0, 0] = gamma
    transformation[0, 1:] = -gamma * velocity
    transformation[1:, 0] = -gamma * velocity
    transformation[1:, 1:] += (gamma - 1) * np.outer(velocity, velocity) / speed_squared
    return transformation


def coupling_design(events, species_order):
    column_indices = {species: index for index, species in enumerate(species_order)}
    blocks = []
    for event in events:
        block = np.zeros((4, len(species_order)))
        for momentum, species, sign in zip(event.momenta, event.species, event.signs):
            block[:, column_indices[species]] += sign * momentum
        blocks.append(block)
    return np.vstack(blocks) if blocks else np.zeros((0, len(species_order)))


def design_summary(events, species_order):
    design = coupling_design(events, species_order)
    singular_values = np.linalg.svd(design, compute_uv=False)
    threshold = 1e-12 * max(float(singular_values.max(initial=0)), 1.0)
    rank = int(np.sum(singular_values > threshold))
    return {
        "species": list(species_order),
        "assumed_nonzero_hard_channels": [list(event.species[:2]) for event in events],
        "rank": rank,
        "nullity": len(species_order) - rank,
        "singular_values": singular_values.tolist(),
        "universal_direction_residual_max": float(np.max(abs(design @ np.ones(len(species_order))), initial=0)),
    }


def interaction_events():
    return (
        elastic_event("A", "B", 0.6, 1.1, momentum=0.8, angle=1.1, azimuth=0.4),
        elastic_event("B", "C", 1.1, 1.7, momentum=0.6, angle=0.9, azimuth=-0.2),
        elastic_event("C", "D", 1.7, 2.2, momentum=1.2, angle=0.7, azimuth=0.8),
    )


def interaction_diagnostics():
    events = interaction_events()
    species_order = ("A", "B", "C", "D")
    disconnected = (events[0], events[2])
    separate_couplings = np.array([1.0, 1.0, 1.4, 1.4])
    return {
        "single_pair": design_summary((events[0],), ("A", "B")),
        "connected": design_summary(events, species_order),
        "disconnected": design_summary(disconnected, species_order),
        "forward_only": design_summary((elastic_event(angle=0),), species_order),
        "separate_sector_couplings": separate_couplings.tolist(),
        "disconnected_residual_norm": float(np.linalg.norm(coupling_design(disconnected, species_order) @ separate_couplings)),
        "connected_residual_norm": float(np.linalg.norm(coupling_design(events, species_order) @ separate_couplings)),
        "scope": "Ranks describe this supplied list of assumed nonzero hard channels. No hard amplitudes or complete interacting theory are constructed. Omitted channels, including graviton-mediated interactions, can connect the apparent sectors.",
    }


def spin_one_diagnostics():
    event = elastic_event()
    soft = null_momentum()
    charges = {"A": 1.0, "B": -2.0}
    current = photon_soft_current(event, charges, soft)
    return {
        "charges": charges,
        "signed_charge_sum": float(np.dot(event.signs, leg_couplings(event, charges))),
        "photon_ward_contraction": float(minkowski_dot(soft, current)),
        "interpretation": "Spin-one gauge decoupling requires charge conservation. The elastic channel conserves each species count and therefore allows unequal charges; it does not require a common charge.",
    }


def lorentz_diagnostics():
    event = elastic_event()
    soft = null_momentum()
    gauge_vector, _ = transverse_basis(soft)
    couplings = {"A": 1.0, "B": 1.25}
    velocity = np.array([0.17, -0.13, 0.21])
    transformation = lorentz_boost(velocity)
    boosted_event = HardEvent(event.momenta @ transformation.T, event.species, event.signs)
    boosted_soft = transformation @ soft
    boosted_gauge = transformation @ gauge_vector
    original = soft_factor(event, couplings, soft, pure_gauge_polarization(soft, gauge_vector))
    boosted = soft_factor(boosted_event, couplings, boosted_soft,
                          pure_gauge_polarization(boosted_soft, boosted_gauge))
    residual_error = ward_residual(boosted_event, couplings) - transformation @ ward_residual(event, couplings)
    return {
        "boost_velocity": velocity.tolist(),
        "max_ward_vector_covariance_error": float(np.max(abs(residual_error))),
        "original_gauge_variation": original,
        "boosted_gauge_variation": boosted,
        "gauge_scalar_difference": boosted - original,
        "interpretation": "Component norms are coordinate diagnostics, not Lorentz scalars. This checks the weighted four-vector and its gauge contraction under one explicit boost.",
    }


def regular_term_diagnostics():
    event = elastic_event()
    couplings = {"A": 1.0, "B": 1.25}
    reference_energy = 0.1
    reference_soft = null_momentum(reference_energy)
    gauge_vector, _ = transverse_basis(reference_soft)
    reference_shift = pure_gauge_polarization(reference_soft, gauge_vector)
    reference_variation = soft_factor(event, couplings, reference_soft, reference_shift)
    contact_tensor = -reference_variation * reference_shift / np.sum(reference_shift**2)
    rows = []
    for energy in (0.1, 0.03, 0.01, 0.003, 0.001):
        soft = null_momentum(energy)
        shift = pure_gauge_polarization(soft, gauge_vector)
        leading = soft_factor(event, couplings, soft, shift)
        regular = float(np.sum(contact_tensor * shift))
        rows.append({
            "soft_energy": energy,
            "leading_gauge_variation": leading,
            "regular_gauge_variation": regular,
            "combined_gauge_variation": leading + regular,
        })
    return {
        "reference_soft_energy": reference_energy,
        "fixed_tensor_max_absolute_component": float(np.max(abs(contact_tensor))),
        "rows": rows,
        "scope": "A deliberately tuned bounded tensor cancels one gauge scalar at one reference energy only; it is not a covariant scattering completion. For any tensor bounded as q tends to zero, contraction with q_mu xi_nu is O(q), unlike the nonuniversal O(1) leading residue.",
    }


def gauge_diagnostics():
    event = elastic_event()
    soft = null_momentum()
    gauge_vector, _ = transverse_basis(soft)
    polarization = plus_polarization(soft)
    gauge_shift = pure_gauge_polarization(soft, gauge_vector)
    rows = []
    for name, couplings in (
        ("universal", {"A": 1.0, "B": 1.0}),
        ("nonuniversal", {"A": 1.0, "B": 1.25}),
    ):
        original = soft_factor(event, couplings, soft, polarization)
        shifted = soft_factor(event, couplings, soft, polarization + gauge_shift)
        residual = ward_residual(event, couplings)
        rows.append({
            "case": name,
            "couplings": couplings,
            "ward_residual": residual.tolist(),
            "coordinate_euclidean_residual_norm": float(np.linalg.norm(residual)),
            "original_soft_factor": original,
            "gauge_shifted_soft_factor": shifted,
            "gauge_variation": shifted - original,
            "predicted_gauge_variation": float(2 * minkowski_dot(residual, gauge_vector)),
        })
    return {
        "hard_momenta": event.momenta.tolist(),
        "species": list(event.species),
        "incoming_outgoing_signs": event.signs.tolist(),
        "soft_momentum": soft.tolist(),
        "gauge_vector": gauge_vector.tolist(),
        "momentum_conservation_residual": np.sum(event.signs[:, None] * event.momenta, axis=0).tolist(),
        "rows": rows,
        "scope": "Leading soft factors with the nonzero hard amplitude divided out; hard momentum conservation is imposed at zero soft momentum. No finite-recoil scattering amplitude is computed.",
    }


class SoftGravitonTests(unittest.TestCase):
    def setUp(self):
        self.event = elastic_event()
        self.soft = null_momentum()
        self.gauge_vector, _ = transverse_basis(self.soft)

    def test_hard_legs_are_on_shell_and_conserve_momentum(self):
        np.testing.assert_allclose(minkowski_dot(self.event.momenta, self.event.momenta),
                                   [0.36, 1.21, 0.36, 1.21], atol=1e-14)
        np.testing.assert_allclose(np.sum(self.event.signs[:, None] * self.event.momenta, axis=0),
                                   0, atol=1e-14)

    def test_soft_momentum_and_allowed_gauge_vector(self):
        self.assertAlmostEqual(float(minkowski_dot(self.soft, self.soft)), 0)
        self.assertAlmostEqual(float(minkowski_dot(self.soft, self.gauge_vector)), 0)
        polarization = plus_polarization(self.soft)
        np.testing.assert_allclose(self.soft @ polarization, 0, atol=1e-14)
        self.assertAlmostEqual(float(np.sum(SIGNATURE * np.diag(polarization))), 0)

    def test_soft_tensor_is_symmetric(self):
        tensor = soft_tensor(self.event, {"A": 1.0, "B": 1.25}, self.soft)
        np.testing.assert_allclose(tensor, tensor.T, atol=1e-14)

    def test_ward_contraction_is_the_weighted_momentum_sum(self):
        couplings = {"A": 1.0, "B": 1.25}
        contraction = (SIGNATURE * self.soft) @ soft_tensor(self.event, couplings, self.soft)
        np.testing.assert_allclose(contraction, ward_residual(self.event, couplings), atol=1e-14)

    def test_universal_coupling_decouples_pure_gauge_polarization(self):
        factor = soft_factor(self.event, {"A": 1.0, "B": 1.0}, self.soft,
                             pure_gauge_polarization(self.soft, self.gauge_vector))
        self.assertAlmostEqual(factor, 0)

    def test_nonuniversal_residual_equals_charge_difference_times_transfer(self):
        transfer = self.event.momenta[2] - self.event.momenta[0]
        residual = ward_residual(self.event, {"A": 1.0, "B": 1.25})
        np.testing.assert_allclose(residual, -0.25 * transfer, atol=1e-14)
        self.assertGreater(float(np.linalg.norm(residual)), 0.1)

    def test_gauge_variation_matches_ward_identity(self):
        for row in gauge_diagnostics()["rows"]:
            self.assertAlmostEqual(row["gauge_variation"], row["predicted_gauge_variation"], places=12)
        self.assertGreater(abs(gauge_diagnostics()["rows"][1]["gauge_variation"]), 0.01)

    def test_leading_soft_factor_scales_inversely_with_energy(self):
        couplings = {"A": 1.0, "B": 1.25}
        polarization = plus_polarization(self.soft)
        original = soft_factor(self.event, couplings, self.soft, polarization)
        rescaled = soft_factor(self.event, couplings, self.soft / 10, polarization)
        self.assertGreater(abs(original), 1.0)
        self.assertAlmostEqual(rescaled, 10 * original, places=10)


    def test_spin_one_allows_unequal_conserved_charges(self):
        diagnostics = spin_one_diagnostics()
        self.assertNotEqual(diagnostics["charges"]["A"], diagnostics["charges"]["B"])
        self.assertAlmostEqual(diagnostics["signed_charge_sum"], 0)
        self.assertAlmostEqual(diagnostics["photon_ward_contraction"], 0)

    def test_connected_channels_fix_ratios_but_not_common_normalization(self):
        events = interaction_events()
        diagnostics = design_summary(events, ("A", "B", "C", "D"))
        self.assertEqual((diagnostics["rank"], diagnostics["nullity"]), (3, 1))
        design = coupling_design(events, ("A", "B", "C", "D"))
        _, _, right_vectors = np.linalg.svd(design)
        np.testing.assert_allclose(abs(right_vectors[-1]), 0.5, atol=1e-13)

    def test_disconnected_channel_list_leaves_relative_sector_normalizations(self):
        diagnostics = interaction_diagnostics()
        self.assertEqual((diagnostics["disconnected"]["rank"], diagnostics["disconnected"]["nullity"]), (2, 2))
        self.assertLess(diagnostics["disconnected_residual_norm"], 1e-13)
        self.assertGreater(diagnostics["connected_residual_norm"], 0.1)

    def test_exact_forward_scattering_gives_no_coupling_constraint(self):
        diagnostics = design_summary((elastic_event(angle=0),), ("A", "B", "C", "D"))
        self.assertEqual((diagnostics["rank"], diagnostics["nullity"]), (0, 4))
        np.testing.assert_allclose(ward_residual(elastic_event(angle=0), {"A": 1.0, "B": 7.0}), 0)

    def test_any_common_coupling_obeys_this_identity(self):
        for coupling in (-0.3, 0.0, 0.7, 4.1):
            residual = ward_residual(self.event, {"A": coupling, "B": coupling})
            np.testing.assert_allclose(residual, 0, atol=1e-14)

    def test_one_gauge_choice_can_miss_a_nonuniversal_residual(self):
        couplings = {"A": 1.0, "B": 1.25}
        residual = ward_residual(self.event, couplings)
        spatial_gauge = np.cross(self.soft[1:], residual[1:])
        spatial_gauge /= np.linalg.norm(spatial_gauge)
        blind_gauge = np.concatenate(([0.0], spatial_gauge))
        self.assertAlmostEqual(float(minkowski_dot(self.soft, blind_gauge)), 0)
        self.assertGreater(float(np.linalg.norm(residual)), 0.1)
        self.assertAlmostEqual(soft_factor(self.event, couplings, self.soft,
                                          pure_gauge_polarization(self.soft, blind_gauge)), 0)

    def test_lorentz_boost_preserves_mass_shell_and_ward_covariance(self):
        transformation = lorentz_boost((0.17, -0.13, 0.21))
        metric = np.diag(SIGNATURE)
        np.testing.assert_allclose(transformation.T @ metric @ transformation, metric, atol=1e-14)
        boosted_momenta = self.event.momenta @ transformation.T
        np.testing.assert_allclose(minkowski_dot(boosted_momenta, boosted_momenta),
                                   minkowski_dot(self.event.momenta, self.event.momenta), atol=1e-14)
        self.assertLess(lorentz_diagnostics()["max_ward_vector_covariance_error"], 1e-13)

    def test_gauge_variation_is_lorentz_invariant(self):
        diagnostics = lorentz_diagnostics()
        self.assertGreater(abs(diagnostics["original_gauge_variation"]), 0.01)
        self.assertLess(abs(diagnostics["gauge_scalar_difference"]), 1e-13)

    def test_bounded_regular_term_gauge_variation_scales_with_soft_energy(self):
        diagnostics = regular_term_diagnostics()
        first = diagnostics["rows"][0]
        for row in diagnostics["rows"]:
            self.assertAlmostEqual(row["regular_gauge_variation"] / first["regular_gauge_variation"],
                                   row["soft_energy"] / first["soft_energy"])
            self.assertAlmostEqual(row["leading_gauge_variation"], first["leading_gauge_variation"])

    def test_one_energy_contact_cancellation_does_not_survive_soft_limit(self):
        rows = regular_term_diagnostics()["rows"]
        self.assertAlmostEqual(rows[0]["combined_gauge_variation"], 0)
        self.assertAlmostEqual(rows[-1]["combined_gauge_variation"],
                               0.99 * rows[-1]["leading_gauge_variation"])
        self.assertGreater(abs(rows[-1]["combined_gauge_variation"]), 0.3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(SoftGravitonTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures), "errors": len(outcome.errors)},
        "gauge_diagnostics": gauge_diagnostics(),
        "interaction_diagnostics": interaction_diagnostics(),
        "spin_one_diagnostics": spin_one_diagnostics(),
        "lorentz_diagnostics": lorentz_diagnostics(),
        "regular_term_diagnostics": regular_term_diagnostics(),
        "scope": "Conditional leading-soft spin-two consistency benchmark, not a derivation of spin two, spacetime, Einstein dynamics, or quantum gravity from cognition.",
    }
    if arguments.write_results:
        output = Path(__file__).with_name("soft_graviton_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()