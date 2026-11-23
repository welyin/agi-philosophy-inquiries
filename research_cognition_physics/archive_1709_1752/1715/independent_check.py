"""A direct-mass certificate for all up-row phases, independent of inverse orbits."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def published_frame():
    solver = runpy.run_path(str(HERE / "../1712/check.py"))
    flow = solver["FL"]
    original = flow["solve"]
    final = {}

    def retain(*args, **kwargs):
        rep = original(*args, **kwargs)
        final.update(rep)
        return rep

    flow["solve"] = retain
    try:
        point = solver["run_point"](**solver["CFG"]["points"]["interior_all_one_width"])
    finally:
        flow["solve"] = original
    original_yd = json.loads((HERE / "../1712/results.json").read_bytes())["cases"]["interior_all_one_width"]["y_d_MZ"]
    assert math.isclose(point["y_d_MZ"], original_yd, rel_tol=1e-9)

    def complex_frame(name):
        z = final[name]
        return np.array(z["real"]) + 1j * np.array(z["imag"])

    # On each open interval the inverse squared gauge couplings are affine in
    # log(mu). Both sides of the finite soft matching appear in these segments.
    gauge_peaks = np.sqrt(np.vstack([g[[0, -1], :2] for _, _, _, g in flow["segments"](.2)]).max(axis=0))
    return (np.array(final["up_texture_diag_10PeV"]), complex_frame("down_texture_10PeV"),
            complex_frame("charged_left_10PeV"), np.array(final["Y_charged_at_10PeV"]),
            flow["MZ"], flow["MS"], flow["END"], flow["TAN"], gauge_peaks)


def certificate(rows, electron_floor, transfer_bound):
    weights = np.sum(abs(rows)**2, axis=1)
    scalar_cap = weights.sum()
    gram = rows @ rows.conj().T
    electron = rows[:, 0]
    norm_lower = weights[0] - np.sum(abs(gram[1:, 0])**2) / weights[0]
    electron_cap = math.sqrt(abs(np.outer(electron, electron.conj()) * gram).sum())
    assert 0 < electron_cap / norm_lower < math.sqrt(electron_floor)
    tolerance = (math.sqrt(electron_floor) * norm_lower - electron_cap) / (1 + math.sqrt(electron_floor))
    threshold = math.log1p(tolerance / scalar_cap)
    assert threshold > transfer_bound
    correction = math.expm1(transfer_bound) * scalar_cap
    bound_after_running = ((electron_cap + correction) / (norm_lower - correction))**2
    # Samples only check conventions and the inequalities. The proof covers the whole torus.
    max_read = 0.
    for angle1 in np.arange(12) * (2 * np.pi / 12):
        for angle2 in np.arange(11) * (2 * np.pi / 11):
            phase = np.exp(1j * np.array([0., angle1, angle2]))
            mass = rows.T @ (phase[:, None] * rows)
            _, singular, vh = np.linalg.svd(mass)
            max_read = max(max_read, abs(vh[0, 0])**2)
            assert singular[0] >= norm_lower * (1 - 1e-12)
            assert np.linalg.norm(mass[:, 0]) <= electron_cap * (1 + 1e-12)
    assert max_read < (electron_cap / norm_lower)**2
    return {
        "phase_independent_row_norm_squared": weights.tolist(),
        "phase_independent_largest_singular_value_lower": norm_lower,
        "phase_independent_electron_column_norm_upper": electron_cap,
        "heaviest_electron_probability_upper": (electron_cap / norm_lower)**2,
        "general_additive_shape_error_needed_over_sum_weights": tolerance / scalar_cap,
        "diagonal_log_condition_needed": threshold,
        "heaviest_electron_probability_upper_after_declared_transport": bound_after_running,
        "electron_direction_min_rotation_deg_after_declared_transport": math.degrees(
            math.asin(math.sqrt(electron_floor)) - math.asin(math.sqrt(bound_after_running))),
        "grid_sanity_maximum_only": max_read,
    }


def compute():
    up, down, charge, charged_y, low, soft, high, tan_beta, gauge_peaks = published_frame()
    np.testing.assert_allclose(charge.conj().T @ charge, np.eye(3), atol=1e-14)
    windows = json.loads((HERE / "../1659/results.json").read_bytes())["joint_windows"]
    (rlo, rhi), (slo, shi) = windows[:2]
    lo13, hi13 = [(1 - math.sqrt(1 - r)) / 2 for r in (rlo, rhi)]
    floor = min(lo13, (1 - hi13) * slo, (1 - hi13) * (1 - shi))
    # Separate model class: below the matching no active singlet or other
    # non-universal source, diagonal charged basis, flavour-universal soft matching.
    # Ignoring positive Yukawa terms bounds the growth on running down.
    assert all(gauge_peaks < [.60, .70])
    g1, g2 = .60, .70
    high_tau = charged_y.max() * math.exp((9*g1**2/5 + 3*g2**2) * math.log(high/soft)/(16*math.pi**2))
    low_tau = high_tau / math.hypot(1, tan_beta) * math.exp(2.25*(g1**2+g2**2) * math.log(soft/low)/(16*math.pi**2))
    max_log = (1.5 * low_tau**2 * math.log(soft / low)
               + high_tau**2 * math.log(high / soft)) / (16 * math.pi**2)
    result = {
        "scale_GeV": [low, soft, high],
        "charged_yukawas_high": charged_y.tolist(),
        "endpoint_gauge_peaks_g1_GUT_and_g2": gauge_peaks.tolist(),
        "normal_small_branch_electron_floor_any_mass_label": floor,
        "one_loop_SM_MSSM_tau_envelopes": [low_tau, high_tau],
        "all_interval_conditional_log_spread": max_log,
    }
    for label, frame in (("registered_frame", charge), ("alternative_conjugate_frame", charge.conj())):
        rows = (down @ frame) / np.sqrt(up[:, None])
        result[label] = certificate(rows, floor, max_log)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    value = json.loads(json.dumps(compute()), parse_float=lambda x: float(f"{float(x):.10g}"))
    target = HERE / "independent_results.json"
    if args.check:
        assert json.loads(target.read_bytes()) == value
        print("1715 verified a direct all-phase bound on the heaviest electron content, including finite shape uncertainty")
    else:
        with target.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
