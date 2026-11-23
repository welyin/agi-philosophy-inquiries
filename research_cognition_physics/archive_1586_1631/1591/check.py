"""Two rank-one mass and dipole channels compared with tagged recoil branches."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def recoil_interval(settings, branch, widths):
    trap, trap_half = settings["offsets"]["trap_with_historical_SM"]
    atom, atom_half = settings["offsets"][branch]
    gap = atom - trap
    total_half = widths * (trap_half + atom_half)
    products = [delta * settings["offset_units"] * slope
                for delta in [gap-total_half, gap+total_half]
                for slope in settings["historical_abs_da_per_inverse_alpha"]]
    return [min(products), max(products)]


def matrix_checks(me, mmu):
    masses = np.diag([me, mmu])
    generator = np.random.default_rng(1591)
    for _ in range(160):
        columns = generator.normal(size=(2, 2)) + 1j*generator.normal(size=(2, 2))
        rows = np.linalg.solve(columns, masses)
        k = generator.uniform(0.1, 5, size=2) * 1e-7
        dipoles = columns @ np.diag(k) @ rows
        x = dipoles @ np.diag([1/me, 1/mmu])
        assert np.allclose(np.sort(np.linalg.eigvals(x).real), np.sort(k), rtol=2e-9, atol=1e-16)
        assert abs(np.prod(x[0, 0] - k) + x[0, 1]*x[1, 0]) < 1e-22
    # Three distinct rank-one channels. The first two contribute to the same mass.
    sources = [np.diag([-2*me, 0]), np.diag([3*me, 0]), np.diag([0, mmu])]
    k = np.array([3e-6, 1e-6, 5e-8])
    dipoles = sum(a*b for a, b in zip(k, sources))
    assert np.allclose(sum(sources), masses)
    assert all(np.linalg.matrix_rank(source) == 1 for source in sources)
    assert not dipoles[0, 1] and not dipoles[1, 0]
    return {"ae": 2*me*dipoles[0, 0], "amu": 2*mmu*dipoles[1, 1]}


def calculate(settings):
    me, mmu = settings["me_GeV"], settings["mmu_GeV"]
    mu_width = settings["hbar_GeV_s"]/settings["mu_lifetime_s"]
    phase = (1-(me/mmu)**2)**3
    rate_factor = settings["alpha"]*mmu**3*phase/(4*mu_width)
    windows = {atom: [recoil_interval(settings, atom, n) for n in (1, 2)]
               for atom in ["Cs", "Rb"]}
    scenarios = {}
    for label, eff, cancel, additive in [
        ("nominal", 1.0, 0.0, 0.0),
        ("finite_transfers", settings["stress_rate_floor"],
         settings["stress_unmodelled_amplitude_fraction"], settings["stress_ae_additive"])
    ]:
        rate = rate_factor*eff*(1-cancel)**2
        eta = math.sqrt(settings["meg_upper"]/(2*me*mmu*rate))
        worst_cs = max(0, -windows["Cs"][1][1]-additive)/(2*me**2)
        rb = (windows["Rb"][0][0]-additive)/(2*me**2) - eta
        scenarios[label] = {
            "nearest_kernel_cap_GeV_minus2": eta,
            "ae_pair_lower": -2*me**2*eta,
            "ae_including_other_lower": -2*me**2*eta-additive,
            "Cs_two_width_min_effective_Br_if_fit": rate*2*me*mmu*worst_cs**2,
            "Rb_one_width_one_mass_upper_GeV": 1/math.sqrt(6*rb) if rb > 0 else None
        }
    return {"recoil_ae_rectangles_one_two": windows, "scenarios": scenarios,
            "three_positive_channels_example": matrix_checks(me, mmu)}


def check(saved, actual):
    if isinstance(actual, dict):
        assert saved.keys() == actual.keys()
        for key in saved:
            check(saved[key], actual[key])
    elif isinstance(actual, list):
        assert len(saved) == len(actual)
        for left, right in zip(saved, actual):
            check(left, right)
    elif actual is None:
        assert saved is None
    else:
        assert math.isclose(saved, actual, rel_tol=1e-10, abs_tol=1e-27)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    opt = parser.parse_args()
    result = calculate(read("inputs.json"))
    if opt.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif opt.check:
        check(read("results.json"), result)
        print("1591: two-channel flavor identity, MEG and electron recoil branches verified")
    else:
        print(json.dumps(result, indent=2))
