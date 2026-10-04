"""Round 408: a calibrated binary effect anchors the residual H-commuting dictionary.

Only helper functions and frozen model data from 407 are reused.  No old test
suite is executed.  Exact claims use finite-field minors; estimated spectral
gaps below are numerical diagnostics, not certified interval lower bounds.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from spectral_tps_audit import tensor_basis, rank_certificate, PRIME

HERE = Path(__file__).resolve().parent
TARGET = HERE/"calibrated_effect_anchor_audit_results.json"
FROZEN = {
    "spectral_tps_audit.py": "c785ddaf92f2716bcd5f77d0bb806b3449b1728f340a51939b241d57e0e34fcf",
    "spectral_tps_audit_results.json": "17fff79f1b72b15daa6bc302e72749b064e63d3d1875ca65b4558a5885a7af15",
}


def weighted_laplacian(effect):
    weights = np.abs(effect)**2
    np.fill_diagonal(weights, 0)
    return np.diag(weights.sum(axis=1))-weights


def pure_distance(a, b):
    return float(np.sqrt(max(0., 1-abs(np.vdot(a, b))**2)))


@lru_cache(None)
def fixed_model():
    for name, sha in FROZEN.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == sha
    row = json.loads((HERE/"spectral_tps_audit_results.json").read_text(encoding="utf-8"))["exact_certificates"][1]
    assert row["qubits"] == 6
    assert row["hankel_certificate"]["rank"] == 64
    assert row["hankel_certificate"]["minor_determinant_mod_prime"] == 57
    strings = [tuple("IXYZ".index(x) for x in s) for s in row["pauli_strings"]]
    coeff = np.array(row["coefficients"], dtype=np.int64)
    h = np.einsum("a,aij->ij", coeff, tensor_basis(strings))
    hm = np.einsum("a,aij->ij", coeff, tensor_basis(strings, True)) % PRIME
    z = np.r_[np.ones(32, dtype=np.int64), -np.ones(32, dtype=np.int64)]
    effect = (np.eye(64)+np.diag(z))/2
    power = np.eye(64, dtype=np.int64)
    columns = []
    for k in range(1, 64):
        power = power@hm % PRIME
        columns.append(((power*z[None, :]-z[:, None]*power) % PRIME).reshape(-1))
    constraint = np.array(columns).T
    values, vectors = np.linalg.eigh(h)
    energy_effect = vectors.conj().T@effect@vectors
    lap = weighted_laplacian(energy_effect)
    return dict(h=h, effect=effect, values=values, vectors=vectors,
                constraint=constraint, energy_effect=energy_effect, lap=lap,
                gap=float(np.linalg.eigvalsh(lap)[1]))


def phase_checks():
    m = fixed_model()
    rng = np.random.default_rng(20260923408)
    rows = []
    for scale in (.0001, .01, .1, 1., 3.):
        phases = scale*rng.normal(size=64)
        w = np.exp(1j*phases)
        diff = (w.conj()[:, None]*w[None, :]-1)*m["energy_effect"]
        error_hs = float(np.linalg.norm(diff, "fro"))
        error_op = float(np.linalg.norm(diff, 2))
        mean = np.mean(w)
        phase = mean/abs(mean) if abs(mean) > 1e-14 else 1.
        distance_hs = float(np.linalg.norm(w-phase))
        distance_op = float(np.max(np.abs(w-phase)))
        variance = float(np.linalg.norm(w-mean)**2)
        bound = error_hs/np.sqrt(m["gap"])
        reference_distances = []
        for ref in (1, 3, 7):
            psi = rng.normal(size=(64, ref))+1j*rng.normal(size=(64, ref))
            psi /= np.linalg.norm(psi)
            moved = w[:, None]*psi
            reference_distances.append(pure_distance(psi, moved))
        rows.append(dict(
            phase_scale=scale, effect_error_hs=error_hs, effect_error_op=error_op,
            laplacian_identity_residual=abs(error_hs**2-2*np.vdot(w, m["lap"]@w).real),
            phase_distance_hs=distance_hs, phase_distance_op=distance_op,
            phase_distance_formula_residual=abs(distance_hs**2-128*(1-abs(mean))),
            laplacian_gap_margin=error_hs**2-2*m["gap"]*variance,
            theoretical_hs_bound_using_estimated_gap=float(bound),
            reference_dimensions=[1, 3, 7], reference_trace_distances=reference_distances))
    return rows


def ic_probabilities(delta):
    n = len(delta)
    diag = np.diag(delta).real
    out = list(diag)
    for j in range(n):
        for k in range(j+1, n):
            out += [(diag[j]+diag[k])/2+delta[j, k].real,
                    (diag[j]+diag[k])/2-delta[j, k].imag]
    return np.array(out)


def ic_reconstruct(values, n):
    a = np.diag(values[:n].astype(complex))
    index = n
    for j in range(n):
        for k in range(j+1, n):
            avg = (values[j]+values[k])/2
            a[j, k] = values[index]-avg+1j*(avg-values[index+1])
            a[k, j] = a[j, k].conjugate()
            index += 2
    assert index == n*n
    return a


def calibration_check():
    rng = np.random.default_rng(40877)
    n = 4
    a = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    a = (a+a.conj().T)/20
    exact = ic_probabilities(a)
    eta = 1e-4
    errors = rng.uniform(-eta, eta, n*n)
    reconstruction = ic_reconstruct(exact, n)
    noisy = ic_reconstruct(exact+errors, n)
    # Independently evaluate the indicated pure-state probabilities.
    vectors = list(np.eye(n, dtype=complex))
    for j in range(n):
        for k in range(j+1, n):
            for phase in (1, 1j):
                v = np.zeros(n, complex)
                v[j], v[k] = 1/np.sqrt(2), phase/np.sqrt(2)
                vectors.append(v)
    direct = np.array([np.vdot(v, a@v).real for v in vectors])
    return dict(dimension=n, preparations=n*n, probability_error_budget=eta,
                probability_formula_error=float(np.max(np.abs(exact-direct))),
                exact_reconstruction_error=float(np.linalg.norm(reconstruction-a)),
                noisy_reconstruction_error=float(np.linalg.norm(noisy-a, "fro")),
                analytic_error_bound=float(np.sqrt(8*n*n-7*n)*eta),
                six_qubit_preparation_count=64**2)


def counterexamples():
    # A connected graph can have an arbitrarily weak bottleneck.
    weak = []
    w = np.diag([-1., 1., 1.])
    psi = np.array([1., 1., 0.])/np.sqrt(2)
    h = np.diag([0., 1., 3.])
    for delta in (.5, .1, .01, .001):
        k = np.array([[0, delta, 0], [delta, 0, .5], [0, .5, 0]])
        effect = (np.eye(3)+k)/2
        gap = float(np.linalg.eigvalsh(weighted_laplacian(effect))[1])
        weak.append(dict(
            delta=delta, effect_min_eigenvalue=float(np.linalg.eigvalsh(effect)[0]),
            effect_max_eigenvalue=float(np.linalg.eigvalsh(effect)[-1]),
            positive_graph_gap=gap,
            anchor_error_op=float(np.linalg.norm(w@effect@w-effect, 2)),
            state_trace_distance=pure_distance(psi, w@psi),
            h_commutator_norm=float(np.linalg.norm(w@h-h@w))))

    disconnected_e = np.array([[.5, .5, 0], [.5, .5, 0], [0, 0, .5]])
    disconnected_w = np.diag([1., 1., -1.])
    probe = np.array([1., 0., 1.])/np.sqrt(2)

    # Degenerate H: a connected graph in one eigenbasis does not suffice.
    deg_h = np.diag([0., 0., 1.])
    u = np.ones(3)/np.sqrt(3)
    v = np.array([1., -1., 0.])/np.sqrt(2)
    deg_e = np.outer(u, u)
    deg_w = np.eye(3)-2*np.outer(v, v)

    # Simple H and connected effect graph, yet passive data miss populations.
    passive = []
    for t in (0., .2, 1.7, 5.):
        wt = np.diag(np.exp(-1j*t*np.diag(h)))
        et = wt.conj().T@deg_e@wt
        passive.append(dict(time=t, energy_eigenstate_probabilities=np.diag(et).real.tolist()))
    return dict(
        weak_links=weak,
        disconnected=dict(
            anchor_commutator=float(np.linalg.norm(disconnected_w@disconnected_e-disconnected_e@disconnected_w)),
            h_commutator=float(np.linalg.norm(disconnected_w@h-h@disconnected_w)),
            state_trace_distance=pure_distance(probe, disconnected_w@probe)),
        degenerate=dict(
            anchor_commutator=float(np.linalg.norm(deg_w@deg_e-deg_e@deg_w)),
            h_commutator=float(np.linalg.norm(deg_w@deg_h-deg_h@deg_w)),
            graph_edges=int(np.count_nonzero(np.triu(abs(deg_e)>1e-12, 1))),
            distance_from_scalar_span=float(np.linalg.norm(deg_w-np.trace(deg_w)*np.eye(3)/3))),
        passive_binary_measurement=passive)


@lru_cache(None)
def report():
    m = fixed_model()
    certificate = rank_certificate(m["constraint"])
    # A population-only calibration is blind to every energy-diagonal W.
    w = np.exp(-.13j*m["values"])
    changed = (w.conj()[:, None]*w[None, :]-1)*m["energy_effect"]
    return dict(
        date="2026-09-23", round=408, scientific_base_through_round=407,
        scope="Residual unitaries commuting with a fixed simple-spectrum H; an independently calibrated binary effect is fixed as an operator. Not a derivation of the effect, physical TPS, locality, or spatial dimension.",
        frozen_source_hashes=FROZEN, old_test_suites_executed=False,
        hilbert_dimension=64, fixed_binary_effect="(I + Z on the first factor)/2 of the frozen six-qubit witness",
        effect_projector_error=float(np.linalg.norm(m["effect"]@m["effect"]-m["effect"])),
        effect_trace=float(np.trace(m["effect"])),
        exact_commutant_certificate=dict(
            prime=PRIME, polynomial_powers=list(range(1, 64)),
            constraint_shape=list(m["constraint"].shape), **certificate),
        estimated_laplacian_gap=m["gap"], estimated_gap_is_certified_lower_bound=False,
        minimum_estimated_energy_basis_offdiagonal_magnitude=float(np.min(abs(m["energy_effect"])[~np.eye(64, dtype=bool)])),
        approximate_anchor_checks=phase_checks(), calibration=calibration_check(),
        boundary_cases=counterexamples(),
        eigenstate_only_calibration=dict(
            maximum_probability_difference=float(np.max(abs(np.diag(changed)))),
            full_effect_difference_op=float(np.linalg.norm(changed, 2))),
        scalar_joint_commutant_proved=True,
        original_h_commuting_unitary_freedom_removed_with_fixed_effect=True,
        robust_unknown_reference_bound_proved=True,
        one_binary_effect_is_one_measurement_shot=False,
        anchor_preparation_derived=False, actual_position_generated=False,
        global_unique_tps_proved=False, full_unitary_controllability_proved=False,
        spatial_dimension_generated=False, full_cognition_to_gr_refuted=False,
        infinite_resource_choice_required=False)


class Audit(unittest.TestCase):
    def test_exact_commutant_certificate_on_frozen_model(self):
        r = report()
        c = r["exact_commutant_certificate"]
        self.assertEqual(c["constraint_shape"], [4096, 63])
        self.assertEqual(c["rank"], 63)
        self.assertEqual(c["minor_determinant_mod_prime"], 26)
        self.assertEqual(r["effect_projector_error"], 0)
        self.assertEqual(r["effect_trace"], 32)

    def test_laplacian_identity_and_approximate_phase_rigidity(self):
        self.assertGreater(report()["estimated_laplacian_gap"], 0)
        for row in report()["approximate_anchor_checks"]:
            self.assertLess(row["laplacian_identity_residual"], 1e-11)
            self.assertLess(row["phase_distance_formula_residual"], 1e-11)
            self.assertGreaterEqual(row["laplacian_gap_margin"], -1e-11)
            self.assertLessEqual(row["phase_distance_hs"], row["theoretical_hs_bound_using_estimated_gap"]+1e-10)

    def test_bound_covers_unknown_reference_inputs(self):
        for row in report()["approximate_anchor_checks"]:
            bound = min(1., row["theoretical_hs_bound_using_estimated_gap"])
            self.assertLessEqual(max(row["reference_trace_distances"]), bound+1e-9)

    def test_binary_effect_calibration_requires_coherent_preparations(self):
        c = report()["calibration"]
        self.assertLess(c["probability_formula_error"], 1e-14)
        self.assertLess(c["exact_reconstruction_error"], 1e-14)
        self.assertLessEqual(c["noisy_reconstruction_error"], c["analytic_error_bound"])
        r = report()["eigenstate_only_calibration"]
        self.assertLess(r["maximum_probability_difference"], 1e-14)
        self.assertGreater(r["full_effect_difference_op"], .1)

    def test_disconnected_anchor_leaves_block_phase_freedom(self):
        r = report()["boundary_cases"]["disconnected"]
        self.assertLess(r["anchor_commutator"], 1e-14)
        self.assertLess(r["h_commutator"], 1e-14)
        self.assertAlmostEqual(r["state_trace_distance"], 1)

    def test_degenerate_spectrum_needs_a_different_criterion(self):
        r = report()["boundary_cases"]["degenerate"]
        self.assertEqual(r["graph_edges"], 3)
        self.assertLess(r["anchor_commutator"], 1e-14)
        self.assertLess(r["h_commutator"], 1e-14)
        self.assertGreater(r["distance_from_scalar_span"], 1)

    def test_connected_weak_anchor_has_no_uniform_noise_guarantee(self):
        rows = report()["boundary_cases"]["weak_links"]
        for row in rows:
            self.assertGreater(row["effect_min_eigenvalue"], 0)
            self.assertLess(row["effect_max_eigenvalue"], 1)
            self.assertGreater(row["positive_graph_gap"], 0)
            self.assertAlmostEqual(row["anchor_error_op"], row["delta"])
            self.assertAlmostEqual(row["state_trace_distance"], 1)
            self.assertEqual(row["h_commutator_norm"], 0)
        self.assertLess(rows[-1]["positive_graph_gap"], 1e-6)

    def test_no_symmetry_does_not_mean_passive_state_tomography(self):
        for row in report()["boundary_cases"]["passive_binary_measurement"]:
            self.assertTrue(np.allclose(row["energy_eigenstate_probabilities"], 1/3, atol=1e-14))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result["checks"] = dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors))
    result["runtime"] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))

