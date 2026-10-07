"""Working 793: complete BV doublets, physical remainder and quantum contact.

Finite rational calibration only. The finite Delta is NOT the continuum
renormalized anomaly, and the algebraic chart is not a new propagator gauge.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv793', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
bv.EVEN = bv.EVEN+('u', 'v')  # external commuting sources, with no BV partners
bv.NAMES = bv.EVEN+bv.ODD
bv.ZERO = (0,)*len(bv.NAMES)
bv.ONE = {bv.ZERO: Q(1)}
V = {x: bv.var(x) for x in bv.NAMES}
AUX = ('q', 'c', 'z', 'p', 'h', 'b', 'rho', 'a')


def project(f):
    return bv.substitute(f, {x: {} for x in AUX})


def h(f):
    out = {}
    for monomial, coefficient in f.items():
        degree = sum(monomial[bv.NAMES.index(x)] for x in AUX)
        if degree:
            term = {monomial: coefficient/Q(degree)}
            out = bv.add(out, bv.mul(V['q'], bv.diff(term, 'c')),
                         bv.mul(V['z'], bv.diff(term, 'p')),
                         bv.mul(V['h'], bv.diff(term, 'b')),
                         bv.scale(bv.mul(V['rho'], bv.diff(term, 'a')), -1))
    return out


def lap(f):
    return bv.add(*(bv.scale(bv.diff(bv.diff(f, anti), field), -1 if field in bv.ODD else 1)
                    for field, anti in bv.PAIRS))


def run():
    q, r, c, p, t, b, a, u, v = [V[x] for x in ('q', 'r', 'c', 'p', 't', 'b', 'a', 'u', 'v')]
    classical = bv.add(bv.scale(bv.mul(r, r), Q(3, 4)),
                       bv.scale(bv.prod(r, r, r), Q(5, 42)),
                       bv.mul(p, c), bv.scale(bv.mul(a, b), -1),
                       bv.mul(u, r), bv.scale(bv.prod(v, r, r), Q(1, 2)))
    s = lambda f: bv.bracket(classical, f)
    assert not bv.bracket(classical, classical) and not lap(classical)
    count = 0
    for degree in range(5):
        for names in combinations_with_replacement(AUX+('r', 't'), degree):
            f = bv.prod(*(V[x] for x in names))
            if not f:
                continue
            assert bv.add(s(h(f)), h(s(f))) == bv.add(f, bv.scale(project(f), -1))
            assert not s(s(f)) and not h(h(f))
            count += 1
    generator = bv.add(bv.prod(p, q, r), bv.scale(bv.prod(t, q, r, r), Q(2, 3)),
                       bv.scale(bv.prod(V['h'], r, r), Q(3, 7)),
                       bv.scale(bv.prod(V['z'], c, r), Q(1, 5)),
                       bv.scale(bv.prod(V['rho'], b, r), Q(2, 11)),
                       bv.scale(bv.prod(u, t, q, r), Q(2, 9)))
    physical = bv.add(bv.scale(r, Q(7, 13)), bv.scale(bv.mul(r, r), Q(8, 17)),
                      bv.scale(bv.prod(u, r, r), Q(3, 5)),
                      bv.scale(bv.prod(u, v, r, r, r), Q(2, 7)))
    first = bv.add(physical, s(generator))
    assert not s(first) and project(first) == physical
    primitive = h(bv.add(first, bv.scale(physical, -1)))
    assert s(primitive) == bv.add(first, bv.scale(physical, -1))
    remove_antifields = lambda f: bv.substitute(f, {anti: {} for _, anti in bv.PAIRS})
    wrong = remove_antifields(first)
    wrong_defect = remove_antifields(s(wrong))
    assert wrong_defect
    physical_source = bv.substitute(bv.diff(physical, 'r'), {x: {} for x in bv.NAMES})
    assert physical_source == bv.scale(bv.ONE, Q(7, 13))
    mixed_source = bv.diff(bv.diff(project(first), 'u'), 'v')
    assert mixed_source == bv.scale(bv.prod(r, r, r), Q(2, 7))
    # Positive-loop quantum canonical orbit of S + epsilon I. The -Delta F
    # term first occurs at epsilon^2 because the generator is epsilon F.
    second = bv.add(bv.bracket(physical, generator),
                    bv.scale(bv.bracket(s(generator), generator), Q(1, 2)),
                    bv.scale(lap(generator), -1))
    defect2 = bv.add(s(second), bv.scale(bv.bracket(first, first), Q(1, 2)), bv.scale(lap(first), -1))
    assert not defect2
    without_contact = bv.add(second, lap(generator))
    omitted = bv.add(s(without_contact), bv.scale(bv.bracket(first, first), Q(1, 2)), bv.scale(lap(first), -1))
    assert omitted and omitted == s(lap(generator))
    return dict(working_round=793, full_BV_homotopy_monomials=count,
                nonminimal_and_antifield_pairs_retained=True,
                physical_remainder=bv.display(physical),
                physical_one_loop_linear_coefficient='7/13',
                retained_mixed_source_contact=bv.display(mixed_source),
                dropping_antifields_breaks_closure=bv.display(wrong_defect),
                quantum_generator_laplacian=bv.display(lap(generator)),
                omitted_quantum_contact_two_loop_defect=bv.display(omitted),
                exact_decomposition_and_two_loop_QME_checks=True,
                all_checks_passed=True,
                original_continuum_quantum_flow_matched=False,
                original_actual_quantum_state_constructed=False,
                formal_round_completed=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('physical_representative_probe_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
