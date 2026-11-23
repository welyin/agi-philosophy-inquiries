"""A posteriori local or end-to-end trit audits on two fixed phase windows."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


OMEGA = np.exp(2j * np.pi / 3)


def density(ket):
    return np.outer(ket, ket.conj())


def scores_from_outputs(amplitudes):
    local = np.zeros((3, 2, 2), dtype=complex)
    end = np.zeros((3, 4, 4), dtype=complex)
    right = np.zeros_like(local)
    for a, b, c in itertools.product(range(3), repeat=3):
        d, e = (b-a) % 3, (c-b) % 3
        phase = np.kron([1, OMEGA**d], [1, OMEGA**e])
        state = density(amplitudes * phase)
        local[d] += np.trace(state.reshape(2, 2, 2, 2), axis1=1, axis2=3)/9
        right[e] += np.trace(state.reshape(2, 2, 2, 2), axis1=0, axis2=2)/9
        end[(c-a) % 3] += state/9
    total = np.zeros(3)
    for k in range(3):
        effect = density(np.array([1, OMEGA**k]))/3
        endpoint_effect = np.eye(4, dtype=complex)/3
        endpoint_effect[np.ix_([0, 3], [0, 3])] = effect
        for j, (matrix, measure) in enumerate(zip([local[k], right[k], end[k]],
                                                    [effect, effect, endpoint_effect])):
            total[j] += np.trace(matrix @ measure).real / 3
    return total


def trial_window(t):
    amplitudes = np.sqrt(np.array([t, 1-t, 1-t, t])/2)
    observed = scores_from_outputs(amplitudes)
    local = (1+2*np.sqrt(t*(1-t)))/3
    endpoint = (1+t)/3
    np.testing.assert_allclose(observed, [local, local, endpoint], atol=1e-12)
    return observed


def random_assisted(seed=290, trials=80):
    """Independent block-norm bounds; a reference may be supplied to the reader."""
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        # Include both an accessible reference of dimension three and discarded noise.
        columns = rng.normal(size=(12, 3)) + 1j*rng.normal(size=(12, 3))
        rho = columns @ columns.conj().T
        rho /= np.trace(rho).real
        block = rho.reshape(4, 3, 4, 3)
        norm = lambda x: np.linalg.svd(x, compute_uv=False).sum()
        populations = np.array([np.trace(block[i, :, i, :]).real for i in range(4)])
        t = populations[0]+populations[3]
        local1 = norm(block[0, :, 2, :] + block[1, :, 3, :])
        local2 = norm(block[0, :, 1, :] + block[2, :, 3, :])
        endpoint = norm(block[0, :, 3, :])
        assert (local1+local2)/2 <= np.sqrt(t*(1-t))+1e-13
        assert endpoint <= t/2+1e-13
    return trials


def longer_chain(n):
    size = 2**n
    bits = np.array(list(itertools.product(range(2), repeat=n)))
    states = np.zeros((3, size, size), dtype=complex)
    for increments in itertools.product(range(3), repeat=n):
        state = OMEGA**(bits @ increments)/np.sqrt(size)
        states[sum(increments) % 3] += density(state)/3**(n-1)
    score = 0.0
    for j in range(3):
        expected = np.eye(size, dtype=complex)/size
        expected[-1, 0] = OMEGA**j/size
        expected[0, -1] = OMEGA**(-j)/size
        np.testing.assert_allclose(states[j], expected, atol=1e-12)
        effect = np.eye(size, dtype=complex)/3
        effect[-1, 0], effect[0, -1] = OMEGA**j/3, OMEGA**(-j)/3
        score += np.trace(states[j] @ effect).real / 3
    np.testing.assert_allclose(score,(1+2**(1-n))/3,atol=1e-12)
    return round(score, 10)


def offline_preparation():
    hadamard = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    amplitudes = np.ones((4, 2))/np.sqrt(8)
    # Project onto "equal route bits or flag zero", then onto the plus flag.
    amplitudes[[1, 2], 1] = 0
    heralded = (amplitudes @ hadamard)[:, 0]
    probability = float(heralded @ heralded)
    np.testing.assert_allclose(probability, 5/8)
    np.testing.assert_allclose(heralded/np.sqrt(probability), np.array([2,1,1,2])/np.sqrt(10))
    return probability


def run():
    for t in np.linspace(0.5, 1, 101):
        trial_window(t)
    cases = {str(t): np.round(trial_window(t), 10).tolist() for t in (0.5, 0.8, 1.0)}
    return {"two_link_examples_local_local_end": cases,
            "balanced_probe_preparation_probability": round(offline_preparation(), 10),
            "pareto_points_recalculated": 101, "assisted_density_matrices": random_assisted(),
            "endpoint_if_every_link_optimal": {str(n):longer_chain(n) for n in range(1,6)}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = run()
    output = Path(__file__).with_name("local_end_audits_results.json")
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        with output.open("x",encoding="utf-8") as handle:
            handle.write(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))
