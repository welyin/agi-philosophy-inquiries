"""1161: reuse the full 1160 process for a conditional leading soft-source map.

Default read-only. --write-once saves new results without replacing history.
No emitted radiation channel or physical graviton rate is implemented here.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OLD = BASE / 'archive_1159_1160/1160/check.py'
OLD_SHA = '8af83d78311b77c804693e055f2ece2c9f195e1d217b5b67e41fa6b3b6753170'
TOL = 1e-9


def small(a, label):
    r = float(np.linalg.norm(a))
    assert r < TOL, (label, r)
    return r


def reuse():
    assert hashlib.sha256(OLD.read_bytes()).hexdigest() == OLD_SHA
    spec = importlib.util.spec_from_file_location('r1160_frozen_source', OLD)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    captured = {}
    previous = sys.getprofile()

    def retain(frame, event, arg):
        if event == 'return' and frame.f_code is old.calculate.__code__:
            for name in ('outputs', 'm', 'gates', 'logical_charges', 'logical_shape'):
                captured[name] = frame.f_locals[name]

    # Observe the unchanged old calculation's final local data, rather than
    # copying its construction or modifying its frozen public return value.
    sys.setprofile(retain)
    try:
        result = old.calculate()
    finally:
        sys.setprofile(previous)
    saved = json.loads(OLD.with_name('results.json').read_text(encoding='utf-8-sig'))

    def compare(a, b):
        if isinstance(a, dict):
            assert a.keys() == b.keys()
            for k in a:
                compare(a[k], b[k])
        elif isinstance(a, list):
            assert len(a) == len(b)
            for x, y in zip(a, b):
                compare(x, y)
        elif isinstance(a, (int, float)):
            assert abs(a-b) <= TOL * max(1, abs(b)), (a, b)
        else:
            assert a == b, (a, b)
    compare(result, saved)
    return old, captured


def soft(particles, n):
    eta = np.diag([1., -1., -1., -1.])
    return sum(np.outer(p, p) / (p @ eta @ n) for p in particles)


def calculate():
    old, data = reuse()
    pa = np.array([[5, 4, 0, 0], [5, 0, 4, 0]], float)
    pb = np.array([[5, -4, 0, 0], [5, 0, -4, 0]], float)
    pair = [[pa[a], pb[b]] for a, b in np.ndindex(2, 2)]
    shape = data['logical_shape']
    ab_index = np.array([ab for qc, ab, rm, re, pr in np.ndindex(shape)])
    gates = data['gates']
    V = gates[3] @ gates[2] @ gates[1] @ gates[0]
    eta = np.diag([1., -1., -1., -1.])
    directions = [np.array([0., 0., 1.]), np.array([1., 0., 0.]),
                  np.array([1., 2., 3.])/math.sqrt(14)]
    worst_ward = worst_tensor = 0.
    old_recoil = old.load_old()
    boost = old_recoil.boost(np.array([.2, -.1, .3]))
    for direction in directions:
        n = np.r_[1., direction]
        tensors = np.array([soft(ps, n) for ps in pair])
        E = np.empty((4, 4, len(V), len(V)), complex)
        for mu, nu in np.ndindex(4, 4):
            d = tensors[ab_index, mu, nu]
            E[mu, nu] = d[:, None]*V - V*d[None, :]
        worst_ward = max(worst_ward, small(np.einsum('m,mnij->nij', eta@n, E), 'full-mother Ward'))
        boosted = np.array([soft([boost@p for p in ps], boost@n) for ps in pair])
        expected = np.einsum('am,kmn,bn->kab', boost, tensors, boost)
        worst_tensor = max(worst_tensor, small(boosted-expected, 'same null vector transported'))
    # Actual endpoints are the frozen 1160 free outgoing events, at dt=2.
    dt = 2.
    p0 = np.array([10., 0., 0., 0.])
    moments = []
    for ps in (pair[0], pair[3]):
        moment = sum(np.outer(dt*p[1:]/p[0], p[1:]) for p in ps)
        moments.append(moment)
    small(moments[0]-np.diag([64/5, 0, 0]), 'xx position-momentum moment')
    small(moments[1]-np.diag([0, 64/5, 0]), 'yy position-momentum moment')
    shifted = []
    displacement = np.array([.3, -.8, 1.2])
    for ps in (pair[0], pair[3]):
        shifted.append(sum(np.outer(dt*p[1:]/p[0]-displacement, p[1:]) for p in ps))
    origin = small((shifted[1]-shifted[0])-(moments[1]-moments[0]), 'origin-independent branch difference')
    small(sum(pair[0])-p0, 'incoming complete pair momentum')
    small(sum(pair[3])-p0, 'outgoing complete pair momentum')

    n = np.array([1., 0, 0, 1.])
    plus = np.diag([0., 1., -1., 0.])/math.sqrt(2)
    cross = np.zeros((4, 4)); cross[1, 2] = cross[2, 1] = 1/math.sqrt(2)
    source = np.array([soft(ps, n)-soft(pair[0], n) for ps in pair])
    s = np.einsum('mn,kmn->k', plus, source)
    sx = np.einsum('mn,kmn->k', cross, source)
    coefficient = -32*math.sqrt(2)/5
    small(source[3]-np.diag([0., -32/5, 32/5, 0.]), 'exact tensor branch difference')
    small(s[[0, 3]]-np.array([0., coefficient]), 'normalized TT-plus')
    rows = []
    max_iso = max_effect = max_source = max_jensen = 0.
    for program in range(3):
        state = data['outputs'][program]
        Vp = state.reshape(-1, 2)
        max_iso = max(max_iso, small(Vp.conj().T@Vp-np.eye(2), 'unknown input/reference isometry'))
        ys = state.copy(); ys[:, :3] = 0
        effect = Vp.conj().T @ ys.reshape(-1, 2)
        expected_effect = np.diag([.5, .5 if program == 1 else 1.])
        max_effect = max(max_effect, small(effect-expected_effect, 'whole output yy effect'))
        E = (state*s[None, :, None, None, None, None]).reshape(-1, 2)
        Ex = (state*sx[None, :, None, None, None, None]).reshape(-1, 2)
        small(Ex, 'TT cross zero in this orientation')
        intensity = E.conj().T@E
        amplitude_mean = Vp.conj().T@E
        max_source = max(max_source, small(intensity-coefficient**2*effect, 'full E dagger E'))
        gap = intensity-amplitude_mean.conj().T@amplitude_mean
        expected_gap = coefficient**2*(expected_effect-expected_effect@expected_effect)
        max_jensen = max(max_jensen, small(gap-expected_gap, 'omitted coherent-source variance'))
        # Entangled Q/R example: contraction uses the same operator certificate.
        bell = np.eye(2)/math.sqrt(2)
        fullR = np.einsum('oq,qr->or', E, bell)
        reference_value = float(np.vdot(fullR, fullR).real)
        small(reference_value-np.trace(intensity).real/2, 'arbitrary reference tensor identity check')
        rows.append(dict(program=program, yy_effect=np.diag(effect).real.tolist(),
                         normalized_plus_coefficient=coefficient,
                         intensity_on_Q0_Q1=np.diag(intensity).real.tolist(),
                         square_of_compressed_mean=np.diag(amplitude_mean.conj().T@amplitude_mean).real.tolist(),
                         variance_gap=np.diag(gap).real.tolist(),
                         bell_reference_intensity=reference_value))
    return dict(schema='round1161_complete_soft_source_v1', status='PASS', check_groups=4,
                date='2026-10-13', old_code_sha256=OLD_SHA,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                inherited_replay='1160 full calculation matches its frozen result',
                new_adopted_axioms=0, new_scientific_calibration=0,
                finite_moment=dict(dt=dt, xx=moments[0].tolist(), yy=moments[1].tolist(),
                                   difference=(moments[1]-moments[0]).tolist(), origin_residual=origin),
                full_mother=dict(logical_dimension=len(V), null_directions=len(directions),
                                 ward_residual=worst_ward, boost_residual=worst_tensor),
                actual_isometries=dict(rows=rows, gram_residual=max_iso, effect_residual=max_effect,
                                      intensity_residual=max_source, variance_residual=max_jensen),
                coefficient_error_budget=dict(isometry_norm_error_per_program=.01,
                     error_per_intensity=2*coefficient**2*.01,
                     program0_minus1_Q1_gap=coefficient**2/2,
                     remaining_coefficient_gap=coefficient**2*(.5-.04)),
                limitations=['spin2/soft factorization and a matching hard S-matrix remain conditions',
                             'soft coefficient map is not a radiation isometry or probability',
                             'no finite-frequency remainder or gravitational feedback is supplied',
                             'moment uses free separated rails and a common rest-frame time slice',
                             'no internal local stress tensor has been reconstructed'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        mode = 'CREATED_ONCE'
    else:
        saved = json.loads(path.read_text(encoding='utf-8'))
        def compare(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for k in a: compare(a[k], b[k])
            elif isinstance(a, list):
                assert len(a) == len(b)
                for x, y in zip(a, b): compare(x, y)
            elif isinstance(a, (int, float)):
                assert abs(a-b) <= TOL*max(1, abs(b))
            else: assert a == b
        compare(result, saved)
        mode = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status='PASS', mode=mode, groups=result['check_groups'],
         ward_residual=result['full_mother']['ward_residual'],
         coefficient_gap=result['coefficient_error_budget']['program0_minus1_Q1_gap']), ensure_ascii=False))


if __name__ == '__main__':
    main()
