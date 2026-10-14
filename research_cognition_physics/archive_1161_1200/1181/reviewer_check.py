"""Independent exact round-1181 audit; no author implementation imported.

The exact magnetic-index construction was written before reading author code.
--probe-only checks mathematics without recording/binding an unfinished draft.
--write records a frozen audit once. Default recomputes and matches it read-only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
INDEPENDENT = HERE / 'reviewer_algebra.py'


def load_independent():
    spec = importlib.util.spec_from_file_location('independent_spin_magnetic', INDEPENDENT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mathematics():
    m = load_independent()
    f, recouple = m.source_maps()
    ft = m.transpose(f)
    g = m.multiply(ft, f)
    vals = (Q(1, 3), Q(35, 81), Q(13, 27))
    projs = m.spectral_projectors(g, vals)
    multiplicities = {str(x): m.trace(projs[x]).text() for x in vals}
    assert multiplicities == {'1/3': '1', '35/81': '1', '13/27': '2'}
    gi = m.sum_matrix([m.scale(projs[x], 1/x) for x in vals], 4)
    isqrt = m.sum_matrix([m.scale(projs[x], m.sqrtq(1/x)) for x in vals], 4)
    v = m.multiply(f, isqrt)
    vt = m.transpose(v)
    projector = m.multiply(m.multiply(f, gi), ft)
    assert m.multiply(projector, f) == f
    assert m.multiply(projector, projector) == projector
    # Full-column left inverse certifies uniqueness algebraically: LF=I;
    # L(P_W F-F)=W G-I. No scalar/diagonal restriction on W is made.
    left = m.multiply(gi, ft)
    assert m.multiply(left, f) == m.identity(4)
    # Geometry effects transported using the same finite isometry.
    a12 = [[m.R(i*(i+1) if i == j else 0) for j in range(4)] for i in range(4)]
    a13 = m.multiply(m.multiply(recouple, a12), m.transpose(recouple))
    assert len({i*(i+1) for i in range(4)}) == 4
    assert all(a13[i][i+1] for i in range(3))
    comm12 = m.minus(m.multiply(g, a12), m.multiply(a12, g))
    comm13 = m.minus(m.multiply(g, a13), m.multiply(a13, g))
    assert not m.zero(comm12) and not m.zero(comm13)
    y = [[m.cg(2, mp, 1, mm, 3, mk) for mk in range(-3, 4, 2)]
         for mp in range(-2, 3, 2) for mm in (-1, 1)]
    assert m.multiply(m.transpose(y), y) == m.identity(4)
    effect_checks = 0
    for basis in (m.identity(4), recouple):
        for j in range(4):
            e = [[basis[i][j]*basis[k][j] for k in range(4)] for i in range(4)]
            ev = m.multiply(m.multiply(v, e), vt)
            assert m.multiply(m.multiply(vt, ev), v) == e
            assert m.multiply(ev, v) == m.multiply(v, e)
            effect_checks += 1
    # Scalar weight error on the original image, independently via endpoints.
    low, high = min(vals), max(vals)
    scalar = Q(2)/(low+high)
    err = (high-low)/(high+low)
    assert scalar == Q(27, 11) and err == Q(2, 11)
    assert scalar*low-1 == -err and scalar*high-1 == err
    # Fixed filter on the existing target boundary, distinct from the optimal
    # (noncontractive) scalar identity approximation.
    p0 = m.multiply(f, ft)
    gain = Q(27, 13)
    k = m.scale(p0, gain)
    compressed = m.multiply(m.multiply(vt, k), v)
    assert compressed == m.scale(g, gain)
    kap = [gain*x for x in vals]
    assert min(kap) == Q(9,13) and max(kap) == 1
    success = ((gain*low)**2+(gain*high)**2)/2
    success_and_return = ((gain*low+gain*high)/2)**2
    conditional_return = success_and_return/success
    assert success == Q(125,169)
    assert success_and_return == Q(121,169)
    assert conditional_return == Q(121,125)
    # Complete mathematical flagged map on all six target dimensions.
    fail = m.minus(m.identity(6), projector)
    for val in vals:
        pi_image = m.multiply(m.multiply(v, projs[val]), vt)
        fail = m.plus(fail, m.scale(pi_image, m.sqrtq(1-(gain*val)**2)))
    assert m.plus(m.multiply(m.transpose(k), k),
                  m.multiply(m.transpose(fail), fail)) == m.identity(6)
    assert m.multiply(vt, v) == m.identity(4)
    return {
        'group_1_exact_geometry': {
            'method': 'Factorial magnetic CG contraction in exact rational radicals, independently written before author code was read.',
            'F': m.texts(f), 'G': m.texts(g), 'A12': m.texts(a12), 'A13': m.texts(a13),
            'spectrum': multiplicities,
            'joint_commutant_dimension': 1,
            'commutant_proof': 'A12 simple spectrum forces a diagonal commuting matrix; nonzero adjacent A13 entries connect all four coordinates.',
            'adjacent_A13_entries': [a13[i][i+1].text() for i in range(3)],
            'G_A12_commutator': m.texts(comm12), 'G_A13_commutator': m.texts(comm13),
            'single_leg_Y_isometry_exact': True,
            'all_exact_checks_pass': True,
        },
        'group_2_identity_pairing': {
            'G_inverse': m.texts(gi),
            'left_inverse_identity_exact': True,
            'image_projector_identity_exact': True,
            'no_new_event_identity_weight': 'G^-1',
            'optimal_scalar': str(scalar), 'operator_error': str(err),
            'norm_scope': 'Image operator identity defect; not by itself an operational probability or channel norm.',
        },
        'group_3_full_quantum_interface': {
            'polar_isometry_exact': True,
            'sharp_geometry_effects_intertwined_exact': effect_checks,
            'filter_gain': str(gain), 'filter_image_eigenvalues': list(map(str,kap)),
            'success_probability': str(success),
            'joint_success_and_return_probability': str(success_and_return),
            'conditional_return_probability': str(conditional_return),
            'unconditional_failure_probability': str(1-success),
            'full_six_dimensional_coherent_failure_isometry_exact': True,
            'reference_scope': 'All stated isometry/intertwining identities tensor with arbitrary reference; the conditional successful branch alone is not deterministic identity.',
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--probe-only', action='store_true')
    args = parser.parse_args()
    data = mathematics()
    if args.probe_only:
        print('INDEPENDENT_EXACT_PROBE_PASS; groups=3; blockers=0')
        return
    author_results = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    sources = author_results['source_sha256']
    for relative, expected in sources.items():
        assert sha(STAGE/relative) == expected
    # This reads only the frozen result encoding after the independent
    # construction is finished. It does not import the author's implementation.
    m = load_independent()
    def decode_text(a):
        return m.texts([[m.R({int(d):Q(c) for d,c in x.items()}) for x in row] for row in a])
    pairs = [('fusion_exact','F'),('gram_exact','G'),('q12_exact','A12'),('q13_exact','A13')]
    for source, own in pairs:
        assert decode_text(author_results['g1_direct_CG'][source]) == data['group_1_exact_geometry'][own]
    assert decode_text(author_results['g2_complete_pairing_identity']['inverse_gram_exact']) == data['group_2_identity_pairing']['G_inverse']
    assert author_results['g3_complete_output_and_instruments']['success_probability_exact'] == '125/169'
    assert author_results['g3_complete_output_and_instruments']['conditional_original_psi_probability_exact'] == '121/125'
    core = [STAGE/'research_note_1181.md'] + [HERE/n for n in
            ('admission.md','proof.md','check.py','results.json','NEXT.md')]
    result = {
        'status': 'PASS_WITH_SCOPE', 'groups': 3, 'blockers': [],
        'independence': 'Prior scope advice disclosed. Exact CG contraction and all three groups implemented without reading/importing author source.',
        'reviewer_code_sha256': sha(Path(__file__)),
        'author_core_sha256': {p.relative_to(STAGE).as_posix():sha(p) for p in core},
        'author_five_sources_sha256_verified': sources,
        'frozen_author_exact_matrices_independently_equal': True,
        'independent_support_sha256': {
            INDEPENDENT.relative_to(STAGE).as_posix():sha(INDEPENDENT)},
        'mathematics': data,
        'adopted_new_axioms': 0, 'GR_hard_deliveries': 0,
    }
    stored = HERE/'review_checks.json'
    if args.write:
        with stored.open('x',encoding='utf-8',newline='\n') as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write('\n')
        print('WRITE_ONCE_PASS; groups=3; blockers=0')
    else:
        assert result == json.loads(stored.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH; groups=3; blockers=0')


if __name__ == '__main__':
    main()
