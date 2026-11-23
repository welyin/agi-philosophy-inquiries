"""One routed visit comparing two independent response bits. Run with --check."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from joint_action import joint_window


def wrap(lanes):
    """Exclusive-or the public target into both private target registers."""
    public = lanes[0][1]
    return (lanes[0], *((setting, y ^ public) for setting, y in lanes[1:]))


def assembled_response(lanes, bit):
    encoded = wrap(lanes)
    returned = joint_window(encoded, (1, bit, 1, bit))
    return wrap(returned)


def check_window():
    rows = 0
    for table in product(range(2), repeat=2):
        for route in range(2):
            for bits in product(range(2), repeat=6):
                lanes = tuple(zip(bits[::2], bits[1::2]))
                output = assembled_response(lanes, table[route])
                expected = ((lanes[0][0], lanes[0][1] ^ table[route]), *lanes[1:])
                assert output == expected
                rows += 1
    return rows


def classically_addressed_scores():
    scores = []
    for route in range(2):
        for answer in product(range(2), repeat=2):
            wins = 0
            for table in product(range(2), repeat=2):
                observed = assembled_response(((0, 0),) * 3, table[route])[0][1]
                wins += (answer[observed] == (table[0] ^ table[1]))
            scores.append(wins / 4)
    assert set(scores) == {0.5}
    return len(scores)


def query_matrix(table):
    matrix = np.zeros((4, 4))
    for route, y in product(range(2), repeat=2):
        lanes = ((0, y), (1, 0), (0, 1))
        new_y = assembled_response(lanes, table[route])[0][1]
        matrix[2 * route + new_y, 2 * route + y] = 1
    return matrix


def parity_outputs(vector, visibility, memory, matrices):
    """Average the actual joint outputs, using independent uniform rule bits."""
    dimension = len(vector)
    split = dimension // 2
    states = [np.zeros((dimension, dimension), dtype=complex) for _ in range(2)]
    for (a, b), op in matrices.items():
        output = np.kron(op, np.eye(memory)) @ vector
        rho = np.outer(output, output.conj())
        rho[:split, split:] *= visibility
        rho[split:, :split] *= visibility
        states[a ^ b] += rho / 2
    return states


def helstrom(rho, sigma):
    return 0.5 + 0.25 * np.abs(np.linalg.eigvalsh(rho - sigma)).sum()


def run():
    window_rows = check_window()
    classical_count = classically_addressed_scores()
    matrices = {bits: query_matrix(bits) for bits in product(range(2), repeat=2)}
    plus = np.array([1., 1.]) / np.sqrt(2)
    minus = np.array([1., -1.]) / np.sqrt(2)
    probe = np.kron(plus, minus)
    x_read_agree = np.kron(np.outer(plus, plus), np.eye(2))
    rows = []
    rng = np.random.default_rng(220722)
    for visibility in (0., 0.25, 0.5, 0.8, 1.):
        even, odd = parity_outputs(probe, visibility, 1, matrices)
        success = (np.trace(x_read_agree @ even) + np.trace((np.eye(4) - x_read_agree) @ odd)).real / 2
        bound = (1 + visibility) / 2
        assert np.isclose(success, bound)
        assert np.isclose(helstrom(even, odd), bound)
        marker0 = np.array([1., 0.])
        marker1 = np.array([visibility, np.sqrt(1 - visibility**2)])
        path_success = helstrom(np.outer(marker0, marker0), np.outer(marker1, marker1))
        assert np.isclose((2*success-1)**2 + (2*path_success-1)**2, 1)
        for _ in range(48):
            v = rng.normal(size=12) + 1j*rng.normal(size=12)
            v /= np.linalg.norm(v)
            states = parity_outputs(v, visibility, 3, matrices)
            best = helstrom(*states)
            # Only the address-resolved |-> components carry the individual signs.
            tagged = v.reshape(2, 2, 3)
            weights = np.linalg.norm((tagged[:, 0] - tagged[:, 1]) / np.sqrt(2), axis=1)
            assert np.isclose(best, 0.5 + visibility * np.prod(weights))
            assert best < bound + 1e-12
        rows.append({"route_overlap": visibility, "best_parity": round(success, 10),
                     "best_path_guess": round(float(path_success), 10)})
    return {"window_basis_inputs": window_rows, "classical_strategies": classical_count,
            "joint_probe_checks": 48 * len(rows), "scores": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("single_window_agreement_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
