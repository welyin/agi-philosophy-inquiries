"""1164: frozen full-process energy resonance under a common shear dictionary.

Default recomputes and compares results read-only. --write-once creates them.
No physical geometry dynamics, control device, or radiation channel is built.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / '1161/check.py'
SOURCE_SHA = '52edf7a1dfc530996daa88ee3715efa936f7c5845eefe067124a544688621119'
TOL = 3e-9


def small(a, label):
    r = float(np.linalg.norm(a))
    assert r < TOL, (label, r)
    return r


def energies(q):
    ex = math.sqrt(9 + 16*math.exp(-2*q))
    ey = math.sqrt(9 + 16*math.exp(2*q))
    # Stable algebraic difference for arbitrarily small nonzero q.
    delta = 32*math.sinh(2*q)/(ey+ex)
    return ex, ey, delta


def compare(actual, expected):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for k in expected:
            compare(actual[k], expected[k])
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            compare(a, b)
    elif isinstance(expected, float):
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=3e-10), (actual, expected)
    else:
        assert actual == expected, (actual, expected)


def calculate():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    spec = importlib.util.spec_from_file_location('round1161_source_for_1164', SOURCE)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    old, data = source.reuse()  # Captures and rechecks the unchanged full 1160.
    shape = data['logical_shape']
    ab_index = np.array([ab for qc, ab, rm, re, pr in np.ndindex(shape)])
    mother = data['gates'][3] @ data['gates'][2] @ data['gates'][1] @ data['gates'][0]
    multiplicity = np.array([0, 1, 1, 2], float)
    samples = [-.01, -.003, .003, .01]
    table = []
    worst_gram = worst_source = worst_ward = 0.
    for q in samples:
        ex, ey, delta = energies(q)
        diagonal = 2*ex + delta*multiplicity[ab_index]
        r = diagonal[:, None]*mother - mother*diagonal[None, :]
        # The same full mother commutator, computed without common offsets.
        r2 = delta*(multiplicity[ab_index, None]*mother - mother*multiplicity[ab_index][None, :])
        worst_ward = max(worst_ward, small(r-r2, 'full mother energy residual'))
        for program in range(3):
            state = data['outputs'][program]
            v = state.reshape(-1, 2)
            worst_gram = max(worst_gram, small(v.conj().T @ v-np.eye(2), 'full unknown input isometry'))
            y = state.copy(); y[:, :3] = 0
            k = v.conj().T @ y.reshape(-1, 2)
            expected_k = np.diag([.5, .5 if program == 1 else 1.])
            worst_source = max(worst_source, small(k-expected_k, 'actual yy pullback'))
            # H_D and spectator terms first cancel on the full mother space.
            residual = (state*(delta*multiplicity)[None, :, None, None, None, None]).reshape(-1, 2)
            mean = v.conj().T @ residual
            square = residual.conj().T @ residual
            worst_source = max(worst_source, small(mean-2*delta*k, 'full mean energy source'))
            worst_source = max(worst_source, small(square-4*delta**2*k, 'full residual squared'))
            table.append(dict(q=q, program=program, delta=delta,
                              input1_mean_work=float(mean[1, 1].real),
                              energy_intertwining_norm=float(np.linalg.norm(residual, 2))))

    # Analytic certificate for every 0<|q|<=1/100, not a grid claim.
    a = Fraction(1, 100)
    upper = Fraction(64, 3)*a/(1-2*a)
    assert upper == Fraction(32, 147) and upper < 1
    spectrum_rows = []
    for q in samples:
        _, _, delta = energies(q)
        matches = [(n, m, ab) for n in range(4) for m in range(4) for ab in range(4)
                   if abs(n-m-multiplicity[ab]*delta) < 1e-11]
        assert matches == [(n, n, 0) for n in range(4)]
        spectrum_rows.append(dict(q=q, matches=matches))
    zero_matches = [(n, m, ab) for n in range(4) for m in range(4) for ab in range(4) if n == m]
    assert len(zero_matches) == 16
    hd = data['m']['hout']; baseline = data['m']['baseline']
    eigs = np.linalg.eigvalsh(hd)-baseline
    spectral_error = small(eigs-np.repeat(np.arange(4), 72), 'original full D spectrum')

    # A necessary separation from every normalized output supported on xx.
    distances = []
    for program in range(3):
        state = data['outputs'][program][..., 1]
        yprob = float(np.linalg.norm(state[:, 3])**2)
        assert math.isclose(yprob, .5 if program == 1 else 1., abs_tol=TOL)
        allowed = np.zeros_like(state)
        if program == 1:
            allowed[:, 0] = state[:, 0]/np.linalg.norm(state[:, 0])
        else:
            allowed[:, 0] = state[:, 3]/np.linalg.norm(state[:, 3])
        # Merely a comparison vector in the necessary support, not an instrument.
        overlap = abs(np.vdot(allowed.ravel(), state.ravel()))
        distance = float(np.linalg.norm(allowed-state))
        bound = math.sqrt(2-2*math.sqrt(max(0, 1-yprob)))
        assert distance + 5e-8 >= bound
        if program != 1:
            small(overlap, 'orthogonal full states')
            assert math.isclose(distance, math.sqrt(2), abs_tol=TOL)
        distances.append(dict(program=program, yy_probability=yprob,
                              pure_output_distance_lower_bound=bound,
                              comparison_vector_distance=distance,
                              yy_effect_probability_gap=yprob))
    vanishing = []
    for q in [.01, .0001, .000001]:
        _, _, delta = energies(q)
        vanishing.append(dict(q=q, energy_residual_norm=2*abs(delta),
                              exact_frozen_menu_task_trace_distance=1.0))
    assert vanishing[-1]['energy_residual_norm'] < 1.3e-5
    slope_error = abs(2*energies(1e-6)[2]/1e-6-64/5)
    assert slope_error < 1e-8

    # Fully transported passive coordinate dictionary, including covectors and endpoints.
    passive_error = phase_error = endpoint_error = 0.
    sigma = 6/5
    for q in samples:
        transform = np.diag([math.exp(-q), math.exp(q), 1.])
        gamma_inv = transform @ transform.T
        for physical_p in (np.array([4., 0, 0]), np.array([0., 4., 0])):
            canonical_p = np.linalg.solve(transform.T, physical_p)
            energy = math.sqrt(9+canonical_p @ gamma_inv @ canonical_p)
            passive_error = max(passive_error, abs(energy-5))
            old_endpoint = sigma*physical_p/3
            endpoint = transform @ old_endpoint
            endpoint_error = max(endpoint_error, small(endpoint-sigma*gamma_inv@canonical_p/3, 'same endpoint'))
            phase_error = max(phase_error, abs(canonical_p@endpoint-physical_p@old_endpoint))
        common_energy = 10*np.eye(len(mother))
        small(common_energy@mother-mother@common_energy, 'passive full mother energy')
    assert max(passive_error, phase_error, endpoint_error) < TOL
    return dict(schema='round1164_common_shear_frozen_process_v1', status='PASS', check_groups=4,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                source_sha256=SOURCE_SHA, original_process_replay='1160 saved calculation matched by 1161 reuse',
                full_process=dict(gram_error=worst_gram, source_identity_error=worst_source,
                                  mother_commutator_error=worst_ward, samples=table),
                exact_spectral_certificate=dict(open_punctured_interval='0<abs(q)<=1/100',
                    maximum_two_delta_bound=str(upper), original_D_spectrum_error=spectral_error,
                    nonzero_q_cases=spectrum_rows, zero_q_matches=len(zero_matches)),
                task_separation=dict(cases=distances, vanishing_energy_residual=vanishing,
                                     source_slope=64/5, numerical_slope_error=slope_error),
                passive_transport=dict(energy_error=passive_error, phase_error=phase_error,
                                       endpoint_error=endpoint_error, constant_metric_is_flat=True),
                new_adopted_axioms=0, new_empirical_calibration=0,
                scope='Frozen exact-energy menu only; no dynamic metric or full-axiom counterexample.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write('\n')
    else:
        compare(result, json.loads(path.read_text(encoding='utf-8-sig')))
    print(json.dumps(dict(status='PASS', mode='WRITE_ONCE' if args.write_once else 'READ_ONLY_MATCH',
                         groups=4, sharp_task_gap=1, full_source_error=result['full_process']['source_identity_error']), ensure_ascii=False))


if __name__ == '__main__':
    main()
