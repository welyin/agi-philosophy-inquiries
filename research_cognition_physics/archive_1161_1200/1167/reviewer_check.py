"""Independent read-only review of round 1167, not an additional research group.

Use spectral polynomials of the original 288D clock and a closed 2x2 Pauli
propagator. Do not use the author's eigensolver for the independent oracle.
The original full 1160 output and all its material factors are retained.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    author = load('review_author1167', HERE/'check.py')
    captured = {}
    previous = sys.getprofile()

    def retain(frame, event, arg):
        if event == 'return' and frame.f_code is author.calculate.__code__:
            for key in ('data', 'exact', 'perturbed'):
                captured[key] = frame.f_locals[key]

    sys.setprofile(retain)
    try:
        current = author.calculate()
    finally:
        sys.setprofile(previous)
    saved = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    author.compare(current, saved)
    leaves = floats = 0
    max_saved_difference = 0.

    def compare(a, b):
        nonlocal leaves, floats, max_saved_difference
        if isinstance(b, dict):
            assert a.keys() == b.keys()
            for key in b:
                compare(a[key], b[key])
        elif isinstance(b, list):
            assert len(a) == len(b)
            for x, y in zip(a, b):
                compare(x, y)
        else:
            leaves += 1
            if isinstance(b, float):
                floats += 1
                max_saved_difference = max(max_saved_difference, abs(a-b))
                assert math.isclose(a, b, rel_tol=2e-9, abs_tol=3e-10)
            else:
                assert a == b
    compare(current, saved)

    data = captured['data']
    hd = data['m']['hout']
    baseline = data['m']['baseline']
    identity = np.eye(288)
    x = hd-baseline*identity
    projectors = []
    for n in range(4):
        proj = identity.copy()
        for k in range(4):
            if k != n:
                proj = proj@(x-k*identity)/(n-k)
        projectors.append(proj)
    projector_residual = max(float(np.linalg.norm(sum(projectors)-identity)),
        *(float(np.linalg.norm(p@p-p)) for p in projectors),
        *(float(np.linalg.norm(x@p-n*p)) for n, p in enumerate(projectors)))
    assert projector_residual < 2e-10

    # Every 2x2 matrix is an effect on the complete unknown input; its transpose
    # is the unnormalized reference block of the normalized Choi state.
    joint_effects = []
    source_block_error = 0.
    weights = np.array([1., 3., 3., 1.])/8
    for program in range(3):
        blocks = np.zeros((4, 4, 2, 2), complex)
        for n, proj in enumerate(projectors):
            for ab in range(4):
                state = data['outputs'][program][:, ab].reshape(288, -1, 2)
                ps = (proj@state.reshape(288, -1)).reshape(state.shape)
                blocks[n, ab] = state.reshape(-1, 2).conj().T@ps.reshape(-1, 2)
                desired = np.zeros((2, 2))
                if ab == 0:
                    desired = weights[n]*np.diag([.5, .5 if program == 1 else 0.])
                if ab == 3:
                    desired = weights[n]*np.diag([.5, .5 if program == 1 else 1.])
                source_block_error = max(source_block_error, float(np.linalg.norm(blocks[n, ab]-desired)))
        joint_effects.append(blocks)
    assert source_block_error < 2e-11

    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.])
    eye = np.eye(2)
    paulis = [sx, sy, sz]
    g0 = np.ones(2)/math.sqrt(2)
    a = np.diag([.8, 1.2])
    hg = eye+.7*sx
    offset = np.diag([2*math.sqrt(16+9*.8**2), 2*math.sqrt(16+9*1.2**2)])
    direction = [-1., 0., 0., 1.]
    waits = [.2, .3]
    menus = [[(eye+sgn*p)/2 for sgn in (1, -1)] for p in (sy, sx)]
    unitarity_error = 0.

    def closed_u(h, t):
        nonlocal unitarity_error
        scalar = np.trace(h).real/2
        v = np.array([np.trace(h@p).real/2 for p in paulis])
        radius = float(np.linalg.norm(v))
        nonscalar = sum(vj*p for vj, p in zip(v, paulis))
        u = np.exp(-1j*scalar*t)*(math.cos(radius*t)*eye
             -1j*t*np.sinc(radius*t/math.pi)*nonscalar)
        unitarity_error = max(unitarity_error, float(np.linalg.norm(u.conj().T@u-eye)))
        return u

    def histories(program, perturbation):
        result = [np.zeros((4, 4), complex) for _ in range(4)]
        for n in range(4):
            for ab in range(4):
                h = hg+(baseline+n+3)*a+offset+perturbation*direction[ab]*sz
                kets = [g0]
                for t, menu in zip(waits, menus):
                    u = closed_u(h, t)
                    kets = [k@u@ket for ket in kets for k in menu]
                for index, ket in enumerate(kets):
                    result[index] += np.kron(np.outer(ket, ket.conj()),
                                             joint_effects[program][n, ab].T/2)
        return result

    independent_exact = [histories(program, 0.) for program in range(3)]
    independent_perturbed = [histories(program, .02) for program in (0, 1)]
    exact_matrix_error = max(float(np.linalg.norm(a-b))
        for oracle, actual in zip(independent_exact, captured['exact']) for a, b in zip(oracle, actual))
    perturbed_matrix_error = max(float(np.linalg.norm(a-b))
        for oracle, actual in zip(independent_perturbed, captured['perturbed']) for a, b in zip(oracle, actual))
    all_outcomes_error = max(abs(sum(np.trace(h).real for h in output)-1)
        for output in independent_exact+independent_perturbed)
    assert max(exact_matrix_error, perturbed_matrix_error, all_outcomes_error) < 2e-11
    distance = sum(np.sum(abs(np.linalg.eigvalsh(a-b)))/2
        for a, b in zip(*independent_perturbed))
    assert abs(distance-saved['relation_future']['perturbed_Choi_trace_distance']) < 2e-11
    assert distance <= saved['relation_future']['universal_half_diamond_upper_bound']

    # Independent central differentiation of all four AB branch energies.
    derivative_error = 0.
    for scale in (.8, 1., 1.2):
        def ab_energy(q):
            ex = math.sqrt(9*scale**2+16*math.exp(-2*q))
            ey = math.sqrt(9*scale**2+16*math.exp(2*q))
            return np.array([2*ex, ex+ey, ex+ey, 2*ey])
        epsilon = 1e-5
        actual = (ab_energy(epsilon)-ab_energy(-epsilon))/(2*epsilon)
        expected = 32/math.sqrt(16+9*scale**2)*np.array(direction)
        derivative_error = max(derivative_error, float(np.linalg.norm(actual-expected)))
    assert derivative_error < 2e-8
    print(json.dumps(dict(status='PASS', kind='independent_review_not_new_scientific_group',
        author_code_sha256=hashlib.sha256((HERE/'check.py').read_bytes()).hexdigest(),
        replay_comparison=dict(leaves=leaves, floats=floats, maximum_absolute_float_difference=max_saved_difference),
        clock_spectral_projector_residual=projector_residual,
        complete_HD_AB_input_effect_error=source_block_error,
        closed_pauli_unitarity_error=unitarity_error,
        exact_full_G_R_history_matrix_error=exact_matrix_error,
        perturbed_full_G_R_history_matrix_error=perturbed_matrix_error,
        complete_history_trace_error=all_outcomes_error,
        independent_perturbed_Choi_trace_distance=float(distance),
        finite_difference_shear_derivative_error=derivative_error), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
