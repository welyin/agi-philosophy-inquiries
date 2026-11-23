"""Independent phase-lattice enclosure and neutral-interference witnesses for round 1661."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
OLD = runpy.run_path(str(ARCHIVE / "1659" / "independent_check.py"))
A = np.array([0, 1, 1])
B = np.array([1, 3, 1])
NULL = np.array([2, -1, 1]) / math.sqrt(6)
PHI = -math.pi / 6


def read(relative):
    return json.loads((ARCHIVE / relative).read_bytes())


def roots(ratio):
    scale = 6 * (ratio + 1 / ratio) + 16
    disc = np.sqrt(scale**2 - 1936)
    return 8 / (scale + disc), (scale + disc) / 242


def window(row):
    return np.array([row[0] - 2 * row[1], row[0] + 2 * row[2]])


def heavy(legs, sign, ratio, r):
    a = A @ legs
    b = B @ legs
    zero = abs(NULL @ legs)**2
    cross = np.real(np.exp(1j * sign * 2 * np.pi / 3) * a.conj() * b)
    return ((abs(a[:, None])**2 / (3 * r) + (4 / 3) * cross[:, None]
             + (11 / 6) * r * abs(b[:, None])**2 - (1 - zero[:, None]) / ratio)
            / (ratio - 1 / ratio))


def phase_ledger(e, ratio, angle_windows):
    lo, hi = angle_windows.T
    reactor = (1 - np.sqrt(1 - np.array([lo[0], hi[0]]))) / 2
    # Exclude the other root of A13 = 4 x (1-x), uniformly in P.
    null_weight_floor = ((2 * abs(e[0, 0]) - sum(abs(e[1:, 0])))**2) / 6
    assert null_weight_floor > .5
    lower = np.array([(1 - reactor[1]) * (1 - hi[1]), reactor[0], (1 - reactor[1]) * lo[2]])
    upper = np.array([(1 - reactor[0]) * (1 - lo[1]), reactor[1], (1 - reactor[0]) * hi[2]])
    scales = np.linalg.norm(e[1:, [0, 0, 1]], axis=0) + 1e-8
    radii = {"minimal": [], "charge_lifts": [], "signed_lifts": [], "all": []}
    rad_pad = .5 * (ratio[1] - ratio[0]) * np.array([[0], [2], [2]]) + 2e-8
    rlo, rhi = roots(ratio)
    rmax = max(rhi)
    min_disc = math.sqrt((6 * (ratio[0] + 1 / ratio[0]) + 16)**2 - 1936)
    deriv_bound = 44 * math.sqrt(6 * rmax * ratio[-1]) / ((ratio[0] - 1 / ratio[0]) * min_disc)
    assert deriv_bound < 2
    for k in range(12):
        for ell in range(12):
            phases = np.exp(1j * PHI * np.array([0, k, ell]))
            legs = phases[:, None] * e[:, :2]
            p0 = np.full(len(ratio), abs(NULL @ legs[:, 0])**2)
            entries = []
            for sign in (-1, 1):
                for label, r in (("small", rlo), ("large", rhi)):
                    p3 = heavy(legs, sign, ratio, r)
                    p = np.vstack([p0, p3])
                    floor = np.sqrt(np.clip(p - rad_pad, 0, 1))
                    ceil = np.sqrt(np.clip(p + rad_pad, 0, 1))
                    wrong = np.maximum(np.sqrt(lower[:, None]) - ceil, floor - np.sqrt(upper[:, None]))
                    chord_needed = np.max(wrong / scales[:, None], axis=0)
                    index = np.argmin(chord_needed)
                    # Chord > 2 already excludes this branch for any phase displacement.
                    degrees = math.degrees(2 * math.asin(min(1., chord_needed[index] / 2)))
                    entries.append(dict(k=k, ell=ell, eta_deg=sign * 120, root=label,
                                        R_center=float(ratio[index]), radius_bound_deg=degrees))
            winning = min(entries, key=lambda entry: entry["radius_bound_deg"])
            radii["all"].append(winning)
            if k % 2 == 0 and ell % 2 == 1:
                radii["signed_lifts"].append(winning)
            if k % 4 == 2 and ell % 4 == 3:
                radii["charge_lifts"].append(winning)
            if (k, ell) == (10, 11):
                radii["minimal"].append(winning)
    counts = [len(vals) for vals in radii.values()]
    assert counts == [1, 9, 36, 144]
    bounds = {name: min(vals, key=lambda entry: entry["radius_bound_deg"]) for name, vals in radii.items()}
    assert bounds["all"]["radius_bound_deg"] > .69 and bounds["charge_lifts"]["radius_bound_deg"] > 2.75
    return dict(R_full_cell_count=len(ratio)-1, R_cover=[ratio[0], ratio[-1]],
                spectral_derivative_ceiling=deriv_bound, high_e3_branch_excluded_by_p0_lower=null_weight_floor,
                probabilities_e1_e3_mu3_box=[lower.tolist(), upper.tolist()],
                phase_class_sizes=dict(zip(radii, counts)), uniform_max_phase_lower_bounds_deg=bounds)


def point(e, y, sign, ratio, ds, phases, windows, freqs):
    r = roots(ratio)[0]
    mass_modes, mass, masses = OLD["mass_source"](r, sign, ds)
    vals = OLD["rates"](e, (mass_modes, mass, masses), *phases)
    lo, hi = windows.T
    assert np.all(lo < vals["angles"]) and np.all(vals["angles"] < hi)
    gaps = np.array([ds, ds*(ratio**2-1)])
    assert np.all(freqs[:, 0] < gaps) and np.all(gaps < freqs[:, 1])
    phases_matrix = np.diag(np.exp(1j * np.r_[0, phases]))
    phys = phases_matrix @ mass @ phases_matrix
    hu = (e * y**2) @ e.conj().T
    def probabilities(charged, neutral):
        ec = np.linalg.eigh(charged)[1]
        w = np.linalg.eigh(neutral.conj().T @ neutral)[1]
        return abs(ec.conj().T @ w)**2
    origin = probabilities(hu, phys)
    d = np.diag(np.exp(.5j * np.pi * np.array([2, 0, 1])))
    other = probabilities(d @ hu @ d.conj().T, d.conj() @ phys @ d.conj().T)
    err = float(np.max(abs(origin - other)))
    assert err < 1e-12
    p_heavy = heavy(phases_matrix @ e[:, :2], sign, np.array([ratio]), np.array([r]))[:, 0]
    np.testing.assert_allclose(origin[:2, 2], p_heavy, atol=3e-14, rtol=0)
    return dict(R=ratio, r=r, gaps_eV2=gaps.tolist(), orbit_probability_error=err, **vals)


def solve():
    point_data = read("1657/independent_results.json")["extended_box_GJ_test"]
    e, y = OLD["charged"](point_data)
    angles = np.array(read("1659/results.json")["joint_windows"])
    osc = read("1484/inputs.json")
    freq = np.array([window(osc["KamLAND_dm21_eV2"]), window(osc["T2K_dm32_NO_eV2"])])
    observed_ratio = [math.sqrt(1+freq[1, 0]/freq[0, 1]), math.sqrt(1+freq[1, 1]/freq[0, 0])]
    ratio = np.linspace(5.49241, 6.22446, 4097)
    assert ratio[0] < observed_ratio[0] < observed_ratio[1] < ratio[-1]
    enclosure = phase_ledger(e, ratio, angles)
    phases = np.array([2, 3]) * PHI
    phases[1] += np.angle(3 + 2 * np.exp(4j * PHI))
    simple = point(e, y, 1, 5.63, 7.54e-5, phases, angles, freq)
    np.testing.assert_allclose(simple["angles"], [.08526891205, .30414859094, .55831001609], atol=1e-9)
    shift = math.radians(11.5)
    c1 = math.sin(shift) / math.sin(4*PHI)
    c0 = math.cos(shift) - c1 * math.cos(4*PHI)
    factor = c0 + c1 * np.exp(4j * PHI)
    np.testing.assert_allclose(factor, np.exp(1j*shift), atol=3e-16)
    other_phases = np.array([10, 7]) * PHI + np.array([np.angle(factor), 0])
    alternate = point(e, y, -1, 5.8804, 7.54e-5, other_phases, angles, freq)
    near = point(e, y, -1, 5.591, 7.8e-5, np.radians([-121.2, -58.9]), angles, freq)
    return dict(R_from_observed_gaps=observed_ratio, phase_lattice_certificate=enclosure,
                third_leg_3_to_2= simple, second_leg_counterpart=dict(real_terms=[c0, c1], **alternate),
                within_1p2_degrees_of_unrestricted_4_2=near)


def stable(x):
    if isinstance(x, dict):
        return {k: stable(v) for k, v in x.items()}
    if isinstance(x, list):
        return [stable(v) for v in x]
    return float(format(x, ".10g")) if isinstance(x, (float, np.floating)) else x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = stable(solve())
    output = HERE / "independent_results.json"
    if args.check:
        assert result == json.loads(output.read_bytes())
    else:
        with output.open("x", encoding="utf8") as f:
            f.write(json.dumps(result, indent=2) + "\n")
    print("1661 independent PASS: all lattice phases and frequency cells, common-basis covariance, two interference witnesses")
