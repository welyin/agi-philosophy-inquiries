"""Finite matrix diagnostic for an unadopted local timing role; default read-only."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'local_timing_support_results_after1245.json'


def run():
    eye = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    h = 3 * np.eye(4) + np.kron(x, x)
    d = np.kron(z, eye) - np.kron(eye, z)
    swap = np.eye(4)[[0, 2, 1, 3]]
    rows = []
    for eta in (0., .25):
        ja, jb = h / 2 + eta * d, h / 2 - eta * d
        total = 1.5 * ja + .5 * jb
        energies, basis = np.linalg.eigh(total)
        unitary = (basis * np.exp(-1j * energies * np.pi / 2)) @ basis.T
        prob = float(abs(unitary[2, 1]) ** 2)
        freq = np.sqrt(1 + 4 * eta * eta)
        block_prob = float(np.sin(freq * np.pi / 2) ** 2 / (freq * freq))
        covariance = float(np.linalg.norm(swap @ ja @ swap.T - jb))
        uniform = float(np.linalg.norm(ja + jb - h))
        unitary_error = float(np.linalg.norm(unitary.conj().T @ unitary - np.eye(4)))
        energy_error = float(np.linalg.norm(unitary.conj().T @ total @ unitary - total))
        assert min(np.linalg.eigvalsh(ja)) >= .5 - 1e-12
        assert min(np.linalg.eigvalsh(jb)) >= .5 - 1e-12
        assert max(covariance, uniform, unitary_error, energy_error, abs(prob-block_prob)) < 1e-12
        rows.append(dict(eta=eta, JA_eigenvalues=np.linalg.eigvalsh(ja).tolist(),
            JB_eigenvalues=np.linalg.eigvalsh(jb).tolist(), total_eigenvalues=energies.tolist(),
            initial_complete_energy=float(total[1, 1]), probability=prob,
            independent_odd_block_probability=block_prob, covariance_residual=covariance,
            uniform_generator_residual=uniform, unitarity_residual=unitary_error,
            complete_energy_residual=energy_error))
    gap = rows[0]['probability'] - rows[1]['probability']
    assert gap >= .2 - 1e-12
    return dict(status='PASS_SCOPE_DIAGNOSTIC', scope='Mathematical timing-role test only; no actual clock or gravitational source realization.',
        N=[1.5, .5], time='pi/2', cases=rows, probability_gap=gap,
        analytic_gap_lower_bound='1/5', analytic_J_lower_bound='I/2',
        new_numbered_rounds=0, new_formal_scientific_groups=0,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def compare(actual, stored):
    if isinstance(actual, dict):
        assert actual.keys() == stored.keys()
        for key in actual:
            compare(actual[key], stored[key])
    elif isinstance(actual, list):
        assert len(actual) == len(stored)
        for a, b in zip(actual, stored):
            compare(a, b)
    elif isinstance(actual, float):
        assert abs(actual - stored) <= 1e-12 * max(1., abs(actual))
    else:
        assert actual == stored


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUTPUT.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        compare(result, json.loads(OUTPUT.read_text('utf-8')))
    print(json.dumps(dict(status=result['status'], mode='INITIAL_SAVE' if args.write else 'READ_ONLY_MATCH',
        gap=result['probability_gap'], analytic_lower_bound='1/5', new_numbered_rounds=0)))
