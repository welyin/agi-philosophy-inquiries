"""Strong-angle torque with one or two real cubic operators on the weak phase circle."""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def stable_minimum(ratio, half_weak):
    # Test potential in units of chi, with positive real explicit amplitude.
    def slope(alpha):
        return 2 * math.sin(2 * alpha) - 3 * ratio * math.sin(3 * (alpha + half_weak))

    def energy(alpha):
        return 2 * math.sin(alpha)**2 + ratio * math.cos(3 * (alpha + half_weak))

    minima = []
    for center in (0, math.pi):
        lo, hi = center - 0.1, center + 0.1
        assert slope(lo) < 0 < slope(hi)
        for _ in range(65):
            mid = (lo + hi) / 2
            if slope(mid) > 0:
                hi = mid
            else:
                lo = mid
        minima.append((energy(mid), abs(2 * (mid - center))))
    # Away from these two wells, 2 sin^2(alpha) exceeds both explicit energies.
    assert 2 * math.sin(0.1)**2 > 2 * ratio
    assert 4 * math.cos(0.2) > 9 * ratio
    return min(minima)[1]


def calculate():
    prior = json.loads((HERE.parent / "1500/independent_results.json").read_text(encoding="utf-8"))
    chi = (prior["physical_lattice_chi_fourth_root_MeV"] / 1000)**4
    chi_high = prior["positive_two_width_rectangle_chi_MeV4"][1] / 1000**4
    theta_target, response_error = 1e-10, 0.01
    bounds = {}
    for degrees in (0, 5):
        radius = 1.5 * (math.radians(degrees) + theta_target)
        min_sine = math.sin(math.pi / 4 - radius)
        for branch in (-1, 1):
            for center in (0, math.pi / 2, math.pi, 3 * math.pi / 2):
                phase = 3 * (branch * math.pi / 4 + center)
                assert abs(abs(math.sin(phase)) - math.sqrt(0.5)) < 3e-15
                for sign in (-1, 1):
                    assert abs(math.sin(phase + sign * radius)) >= min_sine - 3e-15
        ratio = 2 * (1 + response_error) * theta_target / (3 * min_sine)
        bounds[str(degrees)] = {
            "min_abs_sine": min_sine, "B_over_chi_max": ratio,
            "B_GeV4_central_chi_max": chi * ratio,
            "B_fourth_root_MeV_upper_chi_max": 1000 * (chi_high * ratio)**0.25,
            "abs_mu3_GeV_max_at_x_GeV_central_chi": {
                str(v): chi * ratio / (2 * v**3) for v in (1000, 10**11)
            }
        }

    predicted_slope = 1.5 / math.sqrt(2)
    roots = {}
    for r in (1e-5, 1e-7):
        values = [stable_minimum(r, s * math.pi / 4) / r for s in (-1, 1)]
        assert all(abs(value - predicted_slope) < 4 * r for value in values)
        roots[str(r)] = values

    # At fixed weak phase either 0 or pi minimizes both QCD and the sum of equal cubics.
    pair = []
    for branch in (-1, 1):
        psi = branch * math.pi / 4
        for signed_ratio in (-4, 4):
            amplitude = 2 * signed_ratio * math.cos(3 * psi)
            alpha = math.pi if amplitude > 0 else 0
            term = signed_ratio * (math.cos(3 * (psi + alpha)) + math.cos(3 * (psi - alpha)))
            assert abs(term + abs(amplitude)) < 2e-14 and abs(math.sin(alpha)) < 1e-15
            pair.append({"weak_sign": branch, "each_B_over_chi": signed_ratio,
                         "alpha_over_pi": alpha / math.pi, "min_energy_over_chi": term})
    return {"strong_angle_design_target": theta_target, "local_response_margin": response_error,
            "chi_central_GeV4": chi, "bounds_by_weak_error_degrees": bounds,
            "toy_QCD_theta_over_ratio_by_B_over_chi": roots,
            "leading_abs_slope": predicted_slope, "equal_cubics_fixed_weak_circle": pair}


if __name__ == "__main__":
    report = calculate()
    if sys.argv[1:] == ["--write-new"]:
        with OUT.open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")
    elif sys.argv[1:] in ([], ["--check"]):
        assert report == json.loads(OUT.read_text(encoding="utf-8"))
        print("1503 independent: finite torque bounds, two conjugate wells and paired cubic exception verified.")
    else:
        raise SystemExit("Use --write-new or --check.")
