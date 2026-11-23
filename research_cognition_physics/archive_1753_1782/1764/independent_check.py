"""Test the matched double charged-current return with an independent neutral-current rate."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def read(path):
    return json.loads((STAGE / path).read_bytes())


def least_returns(gram_floor, plane_cap, target):
    if plane_cap >= target:
        return 0
    if gram_floor == 1:
        return None
    return math.ceil(math.log(target / plane_cap) / (-2 * math.log(gram_floor)))


def seven_state_test(target):
    masses, theta, mu_eV = runpy.run_path(str(STAGE / "1760/check.py"))["load_model"]()
    # Basis (active, two singlets S, two singlets N), as in round 1760.
    coupling = np.vstack([theta @ np.diag(masses), np.diag(masses)])
    majorana = np.zeros((5, 5), dtype=complex)
    majorana[3:, 3:] = mu_eV * 1e-9
    mass = np.block([[majorana, coupling], [coupling.T, np.zeros((2, 2))]])
    left, singular, _ = np.linalg.svd(mass)
    assert singular[3] > 5000 and singular[4] < 1e-8
    active_light = left[:3, 4:]
    gram = active_light @ active_light.conj().T
    exact_zero_lnv = np.linalg.inv(np.eye(3) + theta @ theta.conj().T)
    np.testing.assert_allclose(gram, exact_zero_lnv, atol=3e-13, rtol=0)
    eig = np.linalg.eigvalsh(gram)
    cap = read("1758/results.json")["profiled_comparison"]["normalized_two_source"][
        "all_symmetric_unitary_sources_mass_ratio_upper"]
    return dict(gram_eigenvalues=eig.tolist(), N_eff=float(eig @ eig),
                matched_massless_kinematics_e_mu_rate_ratio=float((gram[0, 0]/gram[1, 1]).real),
                cap_four=cap/eig[0]**8,
                depth_necessary=least_returns(eig[0], cap, target))


def solve():
    observations = read("1757/independent_results.json")
    solar, atmospheric = observations["independent_KamLAND"], observations["independent_SK"]
    target = math.sqrt(1 + (atmospheric[0]-atmospheric[1])/(solar[0]+solar[2]))
    caps = {name: model["normalized_two_source"]["all_symmetric_unitary_sources_mass_ratio_upper"]
            for name, model in read("1758/results.json").items()
            if name in ("primary_inclusive", "profiled_comparison")}
    plane_cap = max(caps.values())  # uniformly conservative for both parents
    observed_floor = 2.9963 - 2 * .0074  # Janot & Jadach, arXiv:1912.02067
    tests = []
    for additive_budget in (0., .05):
        gram_floor = math.sqrt(observed_floor - additive_budget - 2)
        n = least_returns(gram_floor, plane_cap, target)
        tests.append(dict(additive_budget=additive_budget, gram_floor=gram_floor,
                          mass_ratio_upper_four=plane_cap/gram_floor**8,
                          mass_ratio_upper_five=plane_cap/gram_floor**10,
                          least_returns_necessary=n))
        assert plane_cap / gram_floor**(2*(n-1)) < target <= plane_cap / gram_floor**(2*n)
    # Equal inclusive charged widths can hide arbitrary suppression of a coherent direction.
    v = np.ones(3)/math.sqrt(3)
    demo = []
    for q in (.1, .5, 1.):
        N = np.eye(3) - (1-math.sqrt(q))*np.outer(v, v)
        G = N @ N.T
        np.testing.assert_allclose(G.diagonal()/G[0, 0], np.ones(3), atol=1e-14)
        np.testing.assert_allclose(np.linalg.eigvalsh(G), [q, 1, 1], atol=1e-14)
        demo.append(dict(q=q, equal_charged_ratio=True, neutral_width=float(np.trace(G@G))))
    return dict(archived_plane_caps=caps, target_lower=target, pole_two_width_lower=observed_floor,
                pole_constraints=tests,
                additive_budget_needed_for_four=observed_floor-2-(plane_cap/target)**.25,
                four_step_counterexample_family=demo, canonical_tree=seven_state_test(target))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = json.loads(json.dumps(solve()), parse_float=lambda s: float(f"{float(s):.10g}"))
    path = HERE / "independent_results.json"
    if args.check:
        assert data == json.loads(path.read_bytes())
        print("1764 independent: Z-width singular-value bound and exact light-sector projector verified")
    else:
        with path.open("x", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, allow_nan=False)
            handle.write("\n")
