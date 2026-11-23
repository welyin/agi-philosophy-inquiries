"""Independent two-trial correlation, finite sampling and weak CP checks."""
import argparse
from itertools import combinations
import json
from math import ceil, exp, log, log1p, sqrt
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def run():
    embed = np.array([[1, 0, 0], [0, 1/sqrt(2), 0],
                      [0, 1/sqrt(2), 0], [0, 0, 1]])
    pair = (np.eye(9)[2] - np.eye(9)[4] + np.eye(9)[6]) / sqrt(3)
    gates = [np.eye(3), np.diag([1, -1, 1])]

    def probe(z):
        q = np.array([1, z]) / sqrt(1 + abs(z)**2)
        return embed.T @ np.kron(q, q)

    probabilities = {}
    for correlation in (0., 0.5, 1.):
        cases = []
        for j in range(5):
            state = probe(np.exp(2j*np.pi*j/5))
            probability = 0.
            for l, left in enumerate(gates):
                for k, right in enumerate(gates):
                    weight = (1 + (1 if l == k else -1)*correlation)/4
                    out = np.kron(left @ state, right @ state)
                    probability += weight * abs(np.vdot(pair, out))**2
            cases.append(probability)
        assert np.allclose(cases, (1-correlation)/6)
        probabilities[str(correlation)] = round(float(np.mean(cases)), 12)

    outputs = [g @ probe(1) for g in gates]
    density = sum(np.outer(v,v.conj()) for v in outputs)/2
    eigs = np.linalg.eigvalsh(density)
    assert np.allclose(eigs,[0,.5,.5])
    for pole in [np.eye(3)[0],np.eye(3)[2]]:
        assert all(abs(np.vdot(pole, g @ pole)) == 1 for g in gates)

    uniform_target = .01
    false_certification = .05
    n = ceil(log(false_certification) / log1p(-uniform_target/16))
    risk = lambda size: exp(size * log1p(-uniform_target/16))
    assert risk(n) <= false_certification < risk(n-1)
    rotation = np.array([[.8,.6],[-.6,.8]])
    matrix = np.diag(np.exp(1j*np.array([.73,-.28]))) @ rotation
    matrix = matrix @ np.diag(np.exp(1j*np.array([.45,1.13])))
    spin = embed.T @ np.kron(matrix,matrix) @ embed
    quartets = []
    for i,j in combinations(range(3), 2):
        for k,l in combinations(range(3), 2):
            quartets.append((spin[i,k]*spin[j,l]*spin[i,l].conjugate()*spin[j,k].conjugate()).imag)
    assert max(abs(x) for x in quartets) < 1e-15
    charm_max = -18.2 + 4 * (3.2 + .9)
    assert abs(charm_max+1.8) < 1e-12
    return {
        "probability_by_memory_correlation": probabilities,
        "one_use_output_eigenvalues": [round(float(x),12) for x in eigs],
        "mixed_input_output_distance_to_any_normalized_pure_state_at_least": .5,
        "single_global_test_zero_clicks_per_setting": n,
        "total_pair_attempts": 5*n,
        "single_global_test_false_certification_upper": round(risk(n),12),
        "common_left_frame_CP_quartets_all_zero": True,
        "prompt_charm_CP_difference_four_width_upper_in_units_1e-4": round(charm_max, 10),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    opts = parser.add_mutually_exclusive_group()
    opts.add_argument("--write-new", action="store_true")
    opts.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "independent_results.json"
    if args.write_new:
        with path.open("x", encoding="utf8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    elif args.check:
        assert result == json.loads(path.read_text(encoding="utf-8-sig"))
        print("1441 independent memory, finite confidence and weak CP checks passed")
    else:
        print(json.dumps(result,indent=2))
