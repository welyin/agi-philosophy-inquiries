"""775: auxiliary Feynman contact and the same-state cubic BV defect.

The finite symbols and exact polynomial diagnostic check the analytic map;
they do not compute all original loop coefficients or an interacting state.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_brst_relative_source as blocks
import joint_material_local_brst as poly

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_wick_source_anomaly_results.json'


def mx(x):
    return float(np.max(np.abs(x)))


def auxiliary_contact():
    _, _, k, ks, p, d0, ell, _, _, metric, _ = blocks.setup([1.3, .2, -.3, .4])
    d1 = p+k@ks
    inverse = np.linalg.inv(ell)
    ti = np.eye(111, dtype=complex)
    ti[63:79, :63] = -ks
    adj = lambda x: np.linalg.solve(metric, x.conj().T@metric)
    diagonal = ti@inverse@adj(ti)
    predicted = blocks.block_matrix([
        [np.linalg.inv(d1), None, None, None],
        [None, -np.eye(16), None, None],
        [None, None, None, np.linalg.inv(d0)],
        [None, None, np.linalg.inv(d0), None]], [63, 16, 16, 16])
    t = np.eye(111, dtype=complex)
    t[63:79, :63] = ks
    omitted = diagonal.copy()
    omitted[63:79, 63:79] = 0
    wrong_inverse = t@omitted@adj(t)
    inverse_defect = mx(ell@wrong_inverse-np.eye(111))
    errors = {'Schur_inverse': mx(diagonal-predicted),
              'auxiliary_contact_minus_identity': mx(diagonal[63:79, 63:79]+np.eye(16))}
    assert max(errors.values()) < 2e-13 and inverse_defect > .9

    # State algebra is tested at a null covector, not at the invertible symbol.
    spatial = np.array([1., 2., -3.])/np.sqrt(14.)
    _, _, kn, ksn, _, _, ln, _, _, mn, _ = blocks.setup(np.r_[1., spatial])
    wn = blocks.block_matrix([[np.eye(63), kn, None, None],
                              [ksn, None, None, None],
                              [None, None, None, np.eye(16)],
                              [None, None, np.eye(16), None]], [63, 16, 16, 16])
    tin = np.eye(111, dtype=complex)
    tin[63:79, :63] = -ksn
    adjn = lambda x: np.linalg.solve(mn, x.conj().T@mn)
    wd = tin@wn@adjn(tin)
    state_error = max(mx(ln@wn), mx(wd[63:79, :]), mx(wd[:, 63:79]))
    assert state_error < 2e-13
    return dict(errors=errors, omitted_Feynman_contact_inverse_defect=inverse_defect,
                null_symbol_state_auxiliary_zero_error=state_error,
                stripped_Feynman_auxiliary_contact=-1,
                physical_auxiliary_covariance_not_added=True,
                scope='Original 111-field symbols; noncharacteristic inverse and separate null state calibration. No actual Feynman distribution or interacting state is simulated.')


add, mul, scale, v = poly.add, poly.mul, poly.scale, poly.var


def derivative(p, name):
    return poly.derive(p, {name: poly.ONE}, int(name in poly.ODD))


def s0(p):
    # q1 is an inert physical variable w here, NOT a jet of q0 in this test.
    return poly.derive(p, {'q0': v('c0'), 'a': v('r'), 'z': v('p')}, 1)


def gamma0(p):
    return poly.derive(p, {'q0': v('c0')}, 1)


def wick(p, covariance, sign=-1):
    out, term, n = p, p, 0
    while term:
        n += 1
        term = add(*(scale(derivative(derivative(term, key), key), value)
                     for key, value in covariance.items()))
        term = scale(term, Q(sign, 2*n))
        out = add(out, term)
    return out


def constant(p):
    return p.get((poly.ZERO, 0), Q(0))


def first_vertex_defect():
    # S_2=r^2/2+p*c; S_3 is its compatible cubic BV completion.
    # w=q1 is a flat physical spectator with a nonzero Gaussian covariance.
    kappa, coupling = Q(2, 5), Q(3, 8)
    S3 = add(scale(poly.product(['q0', 'r', 'r']), -1),
             scale(poly.product(['r', 'q1', 'q1']), -kappa),
             scale(poly.product(['q1', 'q1', 'q1']), coupling/6),
             poly.product(['p', 'q0', 'c0']), poly.product(['a', 'r', 'c0']))
    assert not s0(S3)
    H = {'r': Q(3, 7), 'q1': Q(2, 9)}
    normalized = wick(S3, H)
    expected_defect = scale(v('c0'), H['r'])
    defect = add(s0(normalized), scale(wick(s0(S3), H), -1))
    assert defect == expected_defect
    records = []
    for sigma in (Q(5, 4), Q(7, 6)):
        W = {'q1': sigma}
        Nw = wick(S3, W)
        linear = add(normalized, scale(Nw, -1))
        j = {x: constant(wick(derivative(normalized, x), W, 1))
             for x in ('q0', 'r', 'q1')}
        expected_j = {'q0': H['r'], 'r': -kappa*(sigma-H['q1']),
                      'q1': coupling*(sigma-H['q1'])/2}
        assert j == expected_j
        assert linear == add(*(scale(v(x), value) for x, value in j.items()))
        assert not s0(Nw)
        assert s0(linear) == defect
        # Keep a prescribed physical finite part as well as the Ward repair.
        C1 = add(scale(v('q0'), -H['r']), scale(v('r'), Q(4, 11)),
                 scale(v('q1'), Q(-2, 13)))
        assert not add(defect, s0(C1))
        records.append(dict(state_covariance=str(sigma),
                            same_state_raw_source={x: str(y) for x, y in j.items()},
                            gauge_anomaly_coefficient=str(j['q0'])))
    assert records[0]['same_state_raw_source'] != records[1]['same_state_raw_source']
    assert records[0]['gauge_anomaly_coefficient'] == records[1]['gauge_anomaly_coefficient']
    # W intertwines the COMPLETE free BV differential, H only its field part.
    checked = 0
    for n in range(4):
        for names in itertools.combinations_with_replacement(('q0', 'r', 'q1', 'c0', 'p', 'a', 'z'), n):
            p = poly.product(names)
            if not p:
                continue
            assert s0(wick(p, {'q1': Q(5, 4)})) == wick(s0(p), {'q1': Q(5, 4)})
            assert gamma0(wick(p, H)) == wick(gamma0(p), H)
            checked += 1
    # Omitting antifield completion violates the cubic classical master identity.
    physical = add(scale(poly.product(['q0', 'r', 'r']), -1),
                   scale(poly.product(['r', 'q1', 'q1']), -kappa))
    assert s0(physical)
    # Making the EOM direction fluctuate breaks the bisolution prerequisite.
    badW = {'r': Q(1, 5), 'q1': Q(5, 4)}
    bad = add(s0(wick(S3, badW)), scale(wick(s0(S3), badW), -1))
    assert bad == scale(v('c0'), Q(1, 5))
    return dict(exact_monomials_checked=checked, two_same_background_states=records,
                full_free_BV_Wick_intertwining_exact=True,
                parametrix_field_BRST_but_not_full_BV_intertwining=True,
                tadpole_equals_single_vertex_anomaly_projection=True,
                prescribed_source_counterterm_cancels_linear_defect=True,
                wrong_non_bisolution_state_defect='1/5',
                antifield_completion_necessary=True,
                original_continuum_loop_coefficients_computed=False,
                scope='Exact rational cubic BV/Wick diagnostic with a genuine physical covariance parameter. It tests the map and failure controls, not original spacetime loop numbers.')


def run():
    dependencies = ['research_note_735.md', 'research_note_768.md', 'research_note_770.md',
                    'research_note_771.md', 'research_note_772.md', 'research_note_773.md',
                    'research_note_774.md', 'joint_brst_relative_source.py',
                    'joint_material_local_brst.py', 'round775_drafts/STATUS.md']
    return dict(round=775, tests_run=2, failures=0, errors=0,
                auxiliary=auxiliary_contact(), first_vertex=first_vertex_defect(),
                dependency_hashes={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                                   for p in dependencies},
                scope='Direct same-free-state Wick/BV identification of the local one-loop linear anomaly with the raw tadpole Ward defect, plus the original auxiliary Feynman contact. Removes the full in-in-functional prerequisite for this first-jet compatibility only; higher common normalization, interacting physical state, instrument and continuum map remain open.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result, ensure_ascii=False, indent=2))
