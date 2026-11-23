"""Check a finite scalar-loop alternative against the absolute neutral Yukawa arm."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def endpoint():
    previous = runpy.run_path(str(HERE / "../1712/check.py"))
    model = previous["MODEL"]
    simulator = model["sim"]
    graph = simulator["graph"]
    measured = {}

    def record_graph(*args):
        value = graph(*args)
        # graph() returns the column in the mass basis D @ U.conj().
        family = value[0] @ model["state"]["U"].T
        measured["neutral_arm"] = float(abs(family[0, 0]) / simulator["vev"])
        return value

    simulator["graph"] = record_graph
    try:
        charged = previous["run_point"](**previous["CFG"]["points"]["interior_all_one_width"])
    finally:
        simulator["graph"] = graph
    saved = json.loads((HERE / "../1712/results.json").read_bytes())["cases"]["interior_all_one_width"]
    assert abs(charged["y_d_MZ"] / saved["y_d_MZ"] - 1) < 1e-9
    return (measured["neutral_arm"], np.array(charged["upper_probabilities"]),
            np.array(charged["upper_yukawas"]), model["spectra"]["neutral_poles_GeV"][2],
            simulator["vev"])


def loop_function(z, mass_squared):
    return math.log(z / mass_squared) / (1 - mass_squared / z)


def quadrature(z_low, z_high, mass_squared):
    # Substitute u=exp(t) in the finite Euclidean radial integral.
    t = np.linspace(math.log(mass_squared) - 32, math.log(z_high) + 32, 16001)
    u = np.exp(t)
    integrand = (z_high - z_low) * (u / (u + mass_squared)) * (u / (u + z_low)) / (u + z_high)
    return float(np.sum((integrand[1:] + integrand[:-1]) * (np.diff(t) / 2)))


def compute():
    target, weights, masses, lower, vev = endpoint()
    upper = 1e7
    norm_budget = math.sqrt(4 * math.pi)
    z_low, z_high = lower**2, upper**2
    spread = np.array([loop_function(z_high, (y * vev)**2)
                       - loop_function(z_low, (y * vev)**2) for y in masses])
    integrated = np.array([quadrature(z_low, z_high, (y * vev)**2) for y in masses])
    np.testing.assert_allclose(spread, integrated, rtol=3e-13)
    assert np.all((spread > 0) & (spread <= math.log(z_high / z_low)))
    amplitude_caps = (3 / (16 * math.pi**2)) * norm_budget**2 * masses * np.sqrt(weights) * spread / 2
    strange = float(amplitude_caps[1])
    all_flavours = float(np.linalg.norm(amplitude_caps))
    # Orthogonal endpoints in a two-level system at maximal mixing saturate the
    # range bound for all three internal flavours in the same mediator matrix.
    rotate = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    off_diagonals = []
    for y in masses:
        f = [loop_function(z, (y * vev)**2) for z in (z_low, z_high)]
        off_diagonals.append(abs((rotate @ np.diag(f) @ rotate.T)[0, 1]))
    np.testing.assert_allclose(2 * np.array(off_diagonals), spread, rtol=3e-15)
    optimal_right = norm_budget * amplitude_caps / all_flavours
    saturated = float(amplitude_caps @ optimal_right / norm_budget)
    assert math.isclose(saturated, all_flavours, rel_tol=2e-15)
    assert strange < all_flavours < target
    # Test the doublet-chain invariant in a complex basis unrelated to its source.
    rng = np.random.default_rng(17091714)
    frame, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    up = frame[:, 0]
    overlap = np.sqrt(weights[1]) * np.exp(.63j)
    down = up * overlap + frame[:, 2] * np.sqrt(1 - abs(overlap)**2)
    shell = np.zeros((4, 5), complex)
    shell[:3, 0], shell[:3, 1:4] = -down, np.eye(3)
    shell[3, 1:4], shell[3, 4] = -up.conj(), 1
    null = np.r_[1, down, overlap]
    np.testing.assert_allclose(shell @ null, 0, atol=5e-16)
    spectrum = np.linalg.svd(shell, compute_uv=False)**2
    inferred = 6**4 * spectrum.prod() / spectrum.sum()**4 - 2
    assert math.isclose(inferred, weights[1], abs_tol=9e-15)
    unequal = shell.copy()
    unequal[3] *= 1.2
    shifted = np.linalg.svd(unequal, compute_uv=False)**2
    unequal_readout = 6**4 * shifted.prod() / shifted.sum()**4 - 2
    assert abs(unequal_readout - weights[1]) > .01
    return {
        "order_and_scale": "finite one-loop Dirac self-energy at zero external momentum, couplings at the R1712 upper matching point",
        "mediator_mass_interval_GeV": [lower, upper],
        "maximum_endpoint_norm_each": norm_budget,
        "common_scale_CKM_u_row_moduli": np.sqrt(weights).tolist(),
        "common_scale_down_Yukawas": masses.tolist(),
        "matching_y_first_family_N0_required": target,
        "strange_only_maximum_y": strange,
        "all_down_flavours_maximum_y": all_flavours,
        "required_over_strange_cap": target / strange,
        "required_over_all_flavours_cap": target / all_flavours,
        "quadrature_max_absolute_error": float(np.max(abs(spread - integrated))),
        "saturating_norm_of_right_endpoint": float(np.linalg.norm(optimal_right)),
        "doublet_chain_equal_mass_spectrum_inferred_overlap_squared": float(inferred),
        "same_light_direction_with_second_bare_mass_20_percent_higher_spectral_readout": float(unequal_readout),
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    values = compute()
    result = json.loads(json.dumps(values), parse_float=lambda x: float(f"{float(x):.10g}"))
    out = HERE / "independent_results.json"
    if args.check:
        assert json.loads(out.read_bytes()) == result
        print("1714 verified absolute Dirac threshold, finite self-energy integral and sharp all-flavour bound")
    else:
        with out.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
