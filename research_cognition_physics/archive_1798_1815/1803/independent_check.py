"""Linked/empty-path fringe from the same low energy charged-current neutrino data."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def electron_weights(a13, solar, large_root=False):
    small = a13 / (2 * (1 + math.sqrt(1 - a13)))
    reactor = 1 - small if large_root else small
    return np.array([(1 - reactor) * (1 - solar), (1 - reactor) * solar, reactor])


def visibility_bound(weights):
    return max(0.0, 2 * float(max(weights)) - 1)


def pmns(weights, atmos, delta):
    c13, s13 = math.sqrt(sum(weights[:2])), math.sqrt(weights[2])
    c12, s12 = np.sqrt(weights[:2]) / c13
    c23, s23 = math.sqrt(1 - atmos), math.sqrt(atmos)
    phase = np.exp(1j * delta)
    return np.array([
        [c12*c13, s12*c13, s13/phase],
        [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
        [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13],
    ])


def point(saved):
    _, a13, solar, atmos = saved["low_ratio_A13_solar_atmospheric"]
    weights = electron_weights(a13, solar)
    U = pmns(weights, atmos, math.radians(saved["low"]["delta_deg"]))
    np.testing.assert_allclose(U.conj().T @ U, np.eye(3), atol=1e-12)
    # The common optical and rest-mass phases are taken out of the comparison.
    masses = np.array(saved["low"]["masses_eV"])
    delta_m2 = masses**2 - masses[0]**2
    k = 2 * 1.2669328
    baseline_a, baseline_b, energy = 100.0, 250.0, 5.0  # km, km, GeV
    def propagated(baseline):
        return U @ (U[0].conj() * np.exp(-1j * k * delta_m2 * baseline / energy))
    a, b = propagated(baseline_a), propagated(baseline_b)
    overlap = np.vdot(a, b)
    explicit = weights @ np.exp(-1j * k * delta_m2 * (baseline_b - baseline_a) / energy)
    np.testing.assert_allclose(explicit, overlap, atol=1e-12)
    bright = {}
    for holonomy in [1, -1]:
        ports = [float(np.linalg.norm((a + holonomy * np.exp(1j * phase) * b) / 2)**2)
                 for phase in [0, np.pi, np.pi/2, 3*np.pi/2]]
        complex_fringe = (ports[0] - ports[1]) - 1j*(ports[2] - ports[3])
        np.testing.assert_allclose(complex_fringe, holonomy * overlap, atol=1e-12)
        bright[str(holonomy)] = ports
    phases = np.where(weights == max(weights), 1.0, -1.0)
    np.testing.assert_allclose(abs(weights @ phases), visibility_bound(weights), atol=1e-12)
    return dict(weights=weights.tolist(), unconditional_kinematic_phase_visibility=visibility_bound(weights),
                energy_GeV=energy, baselines_km=[baseline_a, baseline_b],
                fringe_real_imag=[float(overlap.real), float(overlap.imag)],
                conditional_port_probabilities_by_holonomy=bright)


def singlet_check():
    eps = np.array([[0., 1.], [-1., 0.]])
    theta, phi, hyper = .41, .78, .67
    a, b = np.cos(theta), np.sin(theta) * np.exp(1j * phi)
    weak = np.array([[a, b], [-b.conjugate(), a]])
    H, lepton = np.array([.3 + .6j, -.17j]), np.array([.43, .1-.13j])
    source = H @ eps @ lepton
    ratios = []
    for holonomy in [1, -1]:
        new_h = np.exp(1j*hyper/2) * weak @ H
        new_l = holonomy * np.exp(-1j*hyper/2) * weak @ lepton
        ratio = (new_h @ eps @ new_l) / source
        np.testing.assert_allclose(ratio, holonomy, atol=1e-12)
        ratios.append(float(ratio.real))
    return ratios


def compute():
    witnesses = read("../1799/results.json")["same_two_tensor_mixed_type_I_and_II"]
    exp_data = read("../1782/inputs.json")
    j, d = exp_data["juno"], exp_data["daya_bay"]
    solar, ds = j["sin_squared_theta12"]
    n = j["marginal_widths"] + j["extra_solar_width_for_ordering_and_external_constraint"]
    solar_interval = [solar - n * ds, solar + n * ds]
    reactor, dr = d["sin_squared_two_theta13"]
    a13_interval = [reactor - d["widths"] * dr, reactor + d["widths"] * dr]
    # On the small reactor root the heaviest weight is (1-s13²)(1-s12²).
    min_first = max(electron_weights(a13_interval[1], solar_interval[1]))
    # On the large reactor root s13² itself is dominant and decreases with A13.
    min_large = max(electron_weights(a13_interval[1], solar_interval[1], True))
    vis = 2*min_first - 1
    # An illustrative *conditional* error budget, not an achieved experiment.
    worst_dominant_packet_loss, drift, error_per_contrast = .1, .03, .02
    min_partial = vis - min_first * worst_dominant_packet_loss
    assert min_partial > 2*error_per_contrast + drift
    return dict(source_parity_after_connected_action=singlet_check(),
                fitted_examples={name: point(val) for name, val in witnesses.items()},
                reactor_solar_two_marginal_box=dict(A13=a13_interval, solar=solar_interval,
                    lower_visibility_small_reactor_root=vis,
                    lower_visibility_large_reactor_root=2*min_large-1,
                    solar_interchange_leaves_max_weight=True),
                conditional_example=dict(max_dominant_packet_overlap_loss=worst_dominant_packet_loss,
                    guaranteed_contrast=min_partial, calibrated_link_insertion_drift_bound=drift,
                    complex_fringe_error_per_setting=error_per_contrast,
                    sum_two_errors_plus_drift=2*error_per_contrast+drift))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = compute()
    raw = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    dest = HERE / "independent_results.json"
    if args.check:
        assert data == read(dest.name)
    else:
        dest.write_bytes(raw.encode("utf-8"))
    print("1803 independent", hashlib.sha256(raw.encode()).hexdigest()[:12])

