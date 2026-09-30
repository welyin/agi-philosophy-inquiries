"""Exact, unnumbered literature-interface checks; no spatial-dimension claim.

Source target: arXiv:2504.19272v1, Lemma 3.18(2), not an unseen journal version.
Python/Fraction certificates decide the checks; NumPy eigenvalues are diagnostics.
Default is read-only. --save uses exclusive creation and never overwrites evidence.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np


def matrix(rows):
    return [[Q(x) for x in row] for row in rows]


def ident(n):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def transpose(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def add(a, b, factor=Q(1)):
    return [[x+factor*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, factor):
    return [[factor*x for x in row] for row in a]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), Q(0))


def det2(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def tensor(a, b):
    return [[a[i][j]*b[k][ell] for j in range(len(a[0]))
             for ell in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def real_expectation_pair(rho_real, rho_imag, real_operator):
    """tr((rho_real + i rho_imag) real_operator), as two exact fractions."""
    return (trace(mul(rho_real, real_operator)),
            trace(mul(rho_imag, real_operator)))


def encoded(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encoded(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encoded(v) for v in value]
    return value


def verify():
    cases = {}
    x = matrix([[2, 0], [0, -1]])
    y = matrix([[2, 1], [1, -1]])
    b = mul(y, x)  # Full-rank pi_x = pi_y = I; A_xy = yx.
    i2 = ident(2)
    rr = scale(i2, Q(1, 2))
    ri = matrix([[0, Q(-1, 2)], [Q(1, 2), 0]])
    assert add(mul(rr, rr), mul(ri, ri), -1) == rr
    assert add(mul(rr, ri), mul(ri, rr)) == ri  # rho^2 = rho.
    assert trace(rr) == 1

    # (1) Exact signature, support, distinct eigenvalues and diagonalization.
    assert trace(x) == trace(y) == 1
    assert det2(x) == -2 and det2(y) == -3
    assert b == matrix([[4, -1], [2, 1]])
    assert trace(b) == 5 and det2(b) == 6
    v = matrix([[1, 1], [2, 1]])
    assert det2(v) == -1
    assert mul(b, v) == mul(v, matrix([[2, 0], [0, 3]]))
    cases['lawful_n1_timelike_diagonalizable_pair'] = {
        'x': x, 'y': y, 'Axy': b, 'spectrum': [2, 3],
        'x_signature': [1, 1], 'y_signature': [1, 1],
        'volume_constraint': 1, 'mean_trace_for_half_half_measure': 1,
        'spectral_weight_squared_double_integral': Q(31),
        'action_minimization_claimed': False,
    }
    # Double-integrated spectral-weight squared: x^2, y^2, xy, yx.
    weights = [trace(mul(x, x)), trace(mul(y, y)), trace(b)]
    assert weights == [5, 7, 5]
    assert (weights[0]**2+weights[1]**2+2*weights[2]**2)/4 == 31

    # (2) Non-real normalized single-state strength, basis-invariant.
    mean = real_expectation_pair(rr, ri, scale(b, Q(1, 2)))
    assert mean == (Q(5, 4), Q(-3, 4))
    o = matrix([[Q(3, 5), Q(-4, 5)], [Q(4, 5), Q(3, 5)]])
    assert mul(o, transpose(o)) == i2
    conjugate = lambda a: mul(mul(o, a), transpose(o))
    assert real_expectation_pair(conjugate(rr), conjugate(ri),
                                 scale(conjugate(b), Q(1, 2))) == mean
    assert b != transpose(b)
    cases['complex_strength_and_change_of_basis'] = {
        'b_u_real_imag': mean, 'basis_invariant': True,
        'Hilbert_self_adjoint': False,
    }

    # (3) Same obstruction at conventional spin dimension two, trace normalized.
    x4 = scale(tensor(x, i2), Q(1, 2))
    y4 = scale(tensor(y, i2), Q(1, 2))
    b4 = mul(y4, x4)
    assert b4 == scale(tensor(b, i2), Q(1, 4))
    assert trace(x4) == trace(y4) == 1
    e00 = matrix([[1, 0], [0, 0]])
    mean4 = real_expectation_pair(tensor(rr, e00), tensor(ri, e00),
                                  scale(b4, Q(1, 4)))
    assert mean4 == (Q(5, 32), Q(-3, 32))
    cases['n2_extension'] = {
        'x_signature': [2, 2], 'y_signature': [2, 2],
        'spectrum_Axy': [Q(1, 2), Q(1, 2), Q(3, 4), Q(3, 4)],
        'b_u_real_imag': mean4, 'each_point_trace': 1,
    }

    # (4) Exact parameter identity. General interval is proved in the note.
    family = []
    for s in [Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)]:
        ys = matrix([[2, s], [s, -1]])
        bs = mul(ys, x)
        assert det2(ys) == -2-s*s
        assert trace(bs) == 5 and det2(bs) == 4+2*s*s
        discriminant = trace(bs)**2-4*det2(bs)
        assert discriminant == 9-8*s*s and discriminant >= 1
        strength = real_expectation_pair(rr, ri, scale(bs, Q(1, 2)))
        assert strength == (Q(5, 4), -3*s/4)
        assert (bs == transpose(bs)) == (s == 0)
        family.append({'s': s, 'discriminant': discriminant,
                       'b_u_real_imag': strength})
    cases['parameter_identity_and_commuting_positive_control'] = family

    # (5) Trace spectral variance remains real; it is not Hermitian variance.
    a = scale(add(b, transpose(b)), Q(1, 2))
    # C = (B-B*)/(2i) = i D; D is real antisymmetric.
    d = scale(add(b, transpose(b), -1), Q(-1, 2))
    scalar_mean = trace(b)/2
    spectral_variance = trace(mul(b, b))/2-scalar_mean**2
    var_a = trace(mul(a, a))/2-(trace(a)/2)**2
    var_c = -trace(mul(d, d))/2
    assert (spectral_variance, var_a, var_c) == (Q(1, 4), Q(5, 2), Q(9, 4))
    assert trace(mul(a, d)) == 0
    assert spectral_variance == var_a-var_c
    assert 2*spectral_variance == Q(1, 2)  # n=1 causal Lagrangian for {2,3}.
    cases['trace_identity_survives_but_not_measurement_variance'] = {
        'mean_trace': scalar_mean, 'spectral_variance': spectral_variance,
        'Hermitian_real_part_variance': var_a,
        'Hermitian_imaginary_part_variance': var_c,
        'causal_lagrangian': Q(1, 2),
    }

    # (6) Two actual binary effects can read real and imaginary parts separately.
    # E_real=(17I+X+3Z)/24, E_imag=(4I-Y)/8.
    ea = matrix([[Q(20, 24), Q(1, 24)], [Q(1, 24), Q(14, 24)]])
    assert ea == add(scale(i2, Q(1, 2)), scale(a, Q(1, 12)))
    assert det2(ea) > 0 and det2(add(i2, ea, -1)) > 0
    assert ea[0][0] > 0 and 1-ea[0][0] > 0
    # E_imag has eigenvalues 3/8,5/8. Its complex expectation uses rho_imag.
    ec_real = scale(i2, Q(1, 2))
    ec_imag = scale(d, Q(1, 12))
    assert ec_imag == matrix([[0, Q(1, 8)], [Q(-1, 8), 0]])
    pa = trace(mul(rr, ea))
    pc = trace(mul(rr, ec_real))-trace(mul(ri, ec_imag))
    assert (pa, pc) == (Q(17, 24), Q(3, 8))
    assert (6*pa-3, 6*pc-3) == mean
    cases['separate_binary_measurement_interface'] = {
        'probabilities': [pa, pc], 'reconstructed_real_imag': mean,
        'measurements_are_separate_settings': True,
        'preparation_and_measurement_implementation_supplied': False,
        'single_joint_sharp_measurement_claimed': False,
    }

    # (7) The paper's complete position-effect family can miss the causal type.
    # Identify the two abstract event labels between the models; this is not
    # identity of the operator-valued point labels in the ambient F_n.
    ys = matrix([[Q(1, 2), 1], [1, Q(1, 2)]])
    bs = mul(ys, x)
    assert trace(ys) == 1 and det2(ys) == Q(-3, 4)
    assert trace(bs) == Q(1, 2) and det2(bs) == Q(3, 2)
    assert trace(bs)**2-4*det2(bs) == Q(-23, 4)
    # Roots (1 +/- i sqrt(23))/4 are conjugate and have equal moduli.
    # All three x,y,ys are invertible, so all support projections are I.
    event_effects = [scale(i2, Q(k, 2)) for k in [0, 1, 1, 2]]
    assert event_effects == [matrix([[0, 0], [0, 0]]),
                             scale(i2, Q(1, 2)), scale(i2, Q(1, 2)), i2]
    assert add(event_effects[1], event_effects[2]) == i2
    cases['same_position_effects_different_causal_type'] = {
        'y_spacelike': ys, 'Axy_spacelike': bs,
        'characteristic_polynomial': [Q(1), Q(-1, 2), Q(3, 2)],
        'discriminant': Q(-23, 4), 'same_effects_for_empty_a_b_ab': event_effects,
        'arbitrary_input_and_passive_reference_statistics_equal': True,
        'identified_abstract_event_labels_not_identical_ambient_subsets': True,
        'complete_dynamics_or_all_instruments_equivalent': False,
    }

    assert len(cases) == 7
    return encoded({
        'date': '2026-09-27', 'evidence_kind': 'unnumbered_literature_interface_audit',
        'source_version': 'arXiv:2504.19272v1', 'checks_run': len(cases),
        'all_exact_checks_passed': True, 'cases': cases,
        'numpy_diagnostic_eigenvalues': np.linalg.eigvals(np.array(b, dtype=float)).tolist(),
        'scope': {
            'frozen_science_through_round': 499, 'numbered_tests_unchanged': 2426,
            'new_numbered_round': None, 'preprint_lemma_3_18_2_counterexample': True,
            'journal_version_verified': False, 'action_minimizer_constructed': False,
            'causal_action_or_continuum_limit_refuted': False,
            'operational_position_or_dimension_generated': False,
            'full_cognition_to_GR_refuted': False, 'full_GR_goal_completed': False,
        },
    })


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save', action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.save:
        result['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        target = Path(__file__).with_name('correlation_geometry_bridge_results.json')
        with target.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    print(json.dumps({'checks_run': result['checks_run'],
                      'all_exact_checks_passed': result['all_exact_checks_passed'],
                      'scope': result['scope']}, ensure_ascii=False))
