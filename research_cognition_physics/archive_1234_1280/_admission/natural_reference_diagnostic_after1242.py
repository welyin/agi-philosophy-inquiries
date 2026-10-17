"""Finite CAR check of an extra endpoint-only reference-return hypothesis.

This is a non-numbered diagnostic, not a sourced gravity solution. Normal replay
compares saved results; --write creates a result exactly once.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve()
RESULT = HERE.with_suffix('.json')


def run():
    ident = np.eye(4, dtype=complex)
    lower = np.array([[0, 1], [0, 0]], dtype=complex)
    parity = np.diag([1, -1]).astype(complex)
    a = np.kron(lower, np.eye(2))
    b = np.kron(parity, lower)
    number = a.conj().T @ a + b.conj().T @ b
    mass = number - ident
    pair = a.conj().T @ b.conj().T + b @ a

    def H(p):
        return 3 * ident + mass + p * pair

    def U(p, time):
        e, v = np.linalg.eigh(H(p))
        return (v * np.exp(-1j * time * e)) @ v.conj().T

    h0 = H(1)
    energies, modes = np.linalg.eigh(h0)
    beta = math.log(4) / math.sqrt(2)
    weights = np.exp(-beta * (energies - energies.min()))
    weights /= weights.sum()
    sigma = (modes * weights) @ modes.conj().T
    # Endpoint quasiparticle pair excitation, NOT the bare CAR |11> projector.
    effect = np.outer(modes[:, -1], modes[:, -1].conj())
    plateau = math.pi / (2 * math.sqrt(5))

    def metrics(transport):
        out = transport @ sigma @ transport.conj().T
        delta = out - sigma
        p0 = float(np.trace(effect @ sigma).real)
        p1 = float(np.trace(effect @ out).real)
        energy = float(np.trace(h0 @ delta).real)
        trace_distance = float(np.sum(np.abs(np.linalg.eigvalsh(delta))) / 2)
        state_entropy = lambda r: -sum(x * math.log(x) for x in np.linalg.eigvalsh(r) if x > 0)
        entropy_change = state_entropy(out) - state_entropy(sigma)
        e, v = np.linalg.eigh(out)
        log_out = (v * np.log(e)) @ v.conj().T
        log_sigma = (modes * np.log(weights)) @ modes.conj().T
        relative_entropy = float(np.trace(out @ (log_out - log_sigma)).real)
        assert abs(entropy_change) < 2e-11
        assert abs(relative_entropy - beta * energy) < 2e-11
        assert np.linalg.norm(transport.conj().T @ transport - ident) < 2e-11
        return dict(p_initial=p0, p_final=p1, probability_gap=p1-p0,
                    endpoint_energy_gain=energy, trace_distance=trace_distance,
                    relative_entropy=relative_entropy,
                    entropy_change=entropy_change)

    ideal = metrics(U(2, plateau))
    expected = dict(p_initial=1/25, p_final=1/10, probability_gap=3/50,
                    endpoint_energy_gain=3*math.sqrt(2)/25,
                    trace_distance=3/(5*math.sqrt(10)),
                    relative_entropy=6*math.log(2)/25)
    for k, value in expected.items():
        assert abs(ideal[k] - value) < 3e-14, (k, ideal[k], value)
    assert np.max(np.abs(weights - np.array([16,4,4,1])/25)) < 1e-15
    assert np.linalg.norm(a@a.conj().T+a.conj().T@a-ident) < 1e-15
    assert np.linalg.norm(a@b+b@a) < 1e-15

    # Smoothstep is C-infinity after attaching constants at 0 and 1.
    def smoothstep(s):
        if s <= 0: return 0.0
        if s >= 1: return 1.0
        x, y = math.exp(-1/s), math.exp(-1/(1-s))
        return x / (x+y)

    width = 1/1000
    smooth = []
    for steps in [256, 512, 1024]:
        dt = width/steps
        up = ident.copy()
        down = ident.copy()
        for i in range(steps):
            s = (i+.5)/steps
            up = U(1+smoothstep(s), dt) @ up
            down = U(2-smoothstep(s), dt) @ down
        transport = down @ U(2, plateau) @ up
        m = metrics(transport)
        reference_transport = U(1, width) @ U(2, plateau) @ U(1, width)
        numerical_norm_gap = float(np.linalg.norm(transport-reference_transport, ord=2))
        assert numerical_norm_gap < 2*width
        assert m['probability_gap'] > 29/500
        smooth.append(dict(steps_per_ramp=steps, **m,
                           numerical_unitary_gap=numerical_norm_gap))

    # Bound is analytical for ALL ramps in [1,2] of this width, not inferred from
    # midpoint convergence or floating-point integration.
    return dict(
        schema='reference_transport_non_numbered_diagnostic_v1', status='PASS',
        source_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),
        numbered_science_round=False,
        assumptions=['Fixed CAR pair model H(p)=3I+N_a+N_b-I+p(pair creation+annihilation)',
                     'Prescribed p=k/a history, not an autonomous sourced metric',
                     'Initial Gibbs reference and fixed endpoint energy effect',
                     'No extra restoring channel at the endpoint'],
        ideal=ideal, exact_expected=expected, smooth_diagnostics=smooth,
        analytic_bounds=dict(ramp_width='1/1000', integrated_operator_difference_le='1/500',
                             probability_gap_ge='29/500',
                             minimum_instantaneous_energy='3-sqrt(5)>0',
                             proof='Duhamel norm bound plus trace-distance contraction'),
        not_claimed=['Complete field state or all-mode stress', 'Actual preparation or detector',
                     'Internal pulse source, conservation and metric feedback',
                     'Einstein reconstruction', 'A counterexample to adopted TR or conjectured CT',
                     'Positive minimal reset work or dissipated heat'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with RESULT.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        status = 'WRITTEN'
    else:
        assert json.loads(RESULT.read_text('utf-8')) == result
        status = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status=status, ideal_probability_gap=result['ideal']['probability_gap'],
                          smooth_probability_lower_bound='29/500',
                          numbered_science_round=False), ensure_ascii=False))
