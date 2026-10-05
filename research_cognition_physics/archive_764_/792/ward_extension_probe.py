"""Working 792: Hom-complex Ward repair need not preserve field independence.

This is a finite unary-operator test of a proposed proof step, not a change
of the already fixed original T1 or a counterexample to continuum existence.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv792', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
bv.EVEN = ('w', 'q', 'b')
bv.ODD = ('c', 'h')
bv.NAMES = bv.EVEN+bv.ODD
bv.ZERO = (0,)*len(bv.NAMES)
bv.ONE = {bv.ZERO: Q(1)}
V = {x: bv.var(x) for x in bv.NAMES}


def d(f):
    return bv.add(bv.mul(V['c'], bv.diff(f, 'q')), bv.mul(V['b'], bv.diff(f, 'h')))


def project(f):
    return bv.substitute(f, {x: {} for x in ('q', 'b', 'c', 'h')})


def homotopy(f):
    out = {}
    for monomial, coeff in f.items():
        n = sum(monomial[bv.NAMES.index(x)] for x in ('q', 'b', 'c', 'h'))
        if n:
            term = {monomial: coeff/Q(n)}
            out = bv.add(out, bv.mul(V['q'], bv.diff(term, 'c')),
                          bv.mul(V['h'], bv.diff(term, 'b')))
    return out


def extension(f):
    return bv.diff(bv.diff(f, 'b'), 'q')


def defect(f):
    return bv.add(d(extension(f)), bv.scale(extension(d(f)), -1))


def hom_defect(f):
    # H(A)=h A + (-1)^|A| P A h for the odd anomaly A=delta E.
    return bv.add(homotopy(defect(f)), bv.scale(project(defect(homotopy(f))), -1))


def repaired(f):
    return bv.add(extension(f), bv.scale(hom_defect(f), -1))


def run():
    checks = 0
    for degree in range(7):
        for names in combinations_with_replacement(bv.NAMES, degree):
            mono = bv.prod(*(V[x] for x in names))
            if not mono:
                continue
            assert not d(d(mono)) and not homotopy(homotopy(mono))
            assert bv.add(d(homotopy(mono)), homotopy(d(mono))) == bv.add(mono, bv.scale(project(mono), -1))
            assert d(repaired(mono)) == repaired(d(mono))
            assert project(repaired(project(mono))) == project(extension(project(mono)))
            for x in bv.NAMES:
                assert bv.diff(extension(mono), x) == extension(bv.diff(mono, x))
            checks += 1
    witness = bv.prod(V['h'], V['q'], V['b'])
    lhs = bv.diff(repaired(witness), 'h')
    rhs = repaired(bv.diff(witness, 'h'))
    assert lhs == bv.scale(bv.ONE, 2) and rhs == bv.scale(bv.ONE, Q(1, 2))
    return dict(working_round=792, monomials_checked=checks,
                cochain_repair_is_Ward_compatible=True,
                cochain_repair_preserves_physical_block=True,
                unrepaired_operator_commutes_with_all_field_derivatives=True,
                witness='h*q*b', repaired_witness=bv.display(repaired(witness)),
                derivative_after_repair=bv.display(lhs), repair_after_derivative=bv.display(rhs),
                field_independence_defect=bv.display(bv.add(lhs, bv.scale(rhs, -1))),
                all_checks_passed=True,
                original_T1_modified=False, original_continuum_normalization_disproved=False,
                simultaneous_original_N1_N2_and_field_Ward_proven=False,
                formal_round_completed=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('ward_extension_probe_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
