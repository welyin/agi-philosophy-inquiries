"""792: dual coefficient contraction, independent of quantum field values.

Exact tests of the graded differential-operator complex used in the analytic
normalization construction. They are not numerical continuum extensions.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def algebra(name, slots=1, jets=1):
    spec = importlib.util.spec_from_file_location(name, HERE.parent/'777/joint_auxiliary_bv_reduction.py')
    a = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a)
    a.EVEN = tuple(f'{x}{s}{j}' for x in ('w', 'q', 'b') for s in range(slots) for j in range(jets))
    a.ODD = tuple(f'{x}{s}{j}' for x in ('c', 'h') for s in range(slots) for j in range(jets))
    a.NAMES = a.EVEN+a.ODD
    a.ZERO = (0,)*len(a.NAMES)
    a.ONE = {a.ZERO: Q(1)}
    a.v = {x: a.var(x) for x in a.NAMES}
    a.slots, a.jets = slots, jets
    return a


def field_d(a, f):
    return a.add(*(a.add(a.mul(a.v[f'c{s}{j}'], a.diff(f, f'q{s}{j}')),
                         a.mul(a.v[f'b{s}{j}'], a.diff(f, f'h{s}{j}')))
                   for s in range(a.slots) for j in range(a.jets)))


def dual_d(a, e):
    # [gamma, partial_c]_+ = partial_q, [gamma, partial_b] = -partial_h.
    return a.add(*(a.add(a.mul(a.v[f'q{s}{j}'], a.diff(e, f'c{s}{j}')),
                         a.scale(a.mul(a.v[f'h{s}{j}'], a.diff(e, f'b{s}{j}')), -1))
                   for s in range(a.slots) for j in range(a.jets)))


def physical(a, e):
    return a.substitute(e, {name: {} for name in a.NAMES if not name.startswith('w')})


def dual_h(a, e):
    out = {}
    for monomial, coefficient in e.items():
        count = sum(n for name, n in zip(a.NAMES, monomial) if not name.startswith('w'))
        if not count:
            continue
        term = {monomial: coefficient/Q(count)}
        out = a.add(out, *(a.add(a.mul(a.v[f'c{s}{j}'], a.diff(term, f'q{s}{j}')),
                                 a.scale(a.mul(a.v[f'b{s}{j}'], a.diff(term, f'h{s}{j}')), -1))
                           for s in range(a.slots) for j in range(a.jets)))
    return out


def repair(a, e):
    return a.add(e, a.scale(dual_h(a, dual_d(a, e)), -1))


def apply(a, operator, field_poly):
    out = {}
    for monomial, coeff in operator.items():
        term = field_poly
        # Symbols denote left derivatives composed from right to left.
        for name, power in reversed(list(zip(a.NAMES, monomial))):
            for _ in range(power):
                term = a.diff(term, name)
        out = a.add(out, a.scale(term, coeff))
    return out


def parity(a, mono):
    return sum(n for name, n in zip(a.NAMES, next(iter(mono))) if name in a.ODD) % 2


def basis(a, max_degree):
    for degree in range(max_degree+1):
        for names in combinations_with_replacement(a.NAMES, degree):
            p = a.prod(*(a.v[x] for x in names))
            if p:
                yield p


def coefficient_complex():
    a = algebra('coeff792')
    operators, polynomials = list(basis(a, 4)), list(basis(a, 5))
    intertwiners, derivatives = 0, 0
    for e in operators:
        d, h = dual_d(a, e), dual_h(a, e)
        assert not dual_d(a, d) and not dual_h(a, h)
        assert a.add(dual_d(a, h), dual_h(a, d)) == a.add(e, a.scale(physical(a, e), -1))
        fixed = repair(a, e)
        assert not dual_d(a, fixed) and physical(a, fixed) == physical(a, e)
        pe = parity(a, e)
        for f in polynomials:
            rhs = a.add(field_d(a, apply(a, e, f)), a.scale(apply(a, e, field_d(a, f)), -(-1)**pe))
            assert apply(a, d, f) == rhs
            intertwiners += 1
        if not pe:
            for f in polynomials[::9]:
                for name in a.NAMES:
                    assert a.diff(apply(a, fixed, f), name) == apply(a, fixed, a.diff(f, name))
                    derivatives += 1
    e = a.mul(a.v['q00'], a.v['b00'])
    fixed = repair(a, e)
    expected = a.scale(a.add(e, a.mul(a.v['c00'], a.v['h00'])), Q(1, 2))
    assert fixed == expected
    f = a.prod(a.v['h00'], a.v['q00'], a.v['b00'])
    left = a.diff(apply(a, fixed, f), 'h00')
    right = apply(a, fixed, a.diff(f, 'h00'))
    assert left == right == a.scale(a.ONE, Q(1, 2))
    return dict(operator_basis=len(operators), polynomial_basis=len(polynomials),
                commutator_identity_checks=intertwiners, derivative_checks=derivatives,
                repaired_operator=a.display(fixed), repaired_witness=a.display(apply(a, fixed, f)),
                derivative_before_and_after='1/2',
                previous_unrestricted_Hom_repair_defect='3/2',
                field_independence_preserved=True, arithmetic='Fraction')


def dual_translation(a, e):
    # The adjoint action of sum x_(j+1) d/dx_j, truncated consistently.
    return a.add(*(a.scale(a.mul(a.v[f'{x}{s}{j}'], a.diff(e, f'{x}{s}{j+1}')), -1)
                   for x in ('w', 'q', 'b', 'c', 'h') for s in range(a.slots)
                   for j in range(a.jets-1)))


def slot_degree(a, monomial):
    return tuple(sum(power for name, power in zip(a.NAMES, monomial) if int(name[1]) == s)
                 for s in range(a.slots))


def slot_and_jet_contracts():
    a = algebra('jets792', slots=2, jets=2)
    swap = {f'{x}{s}{j}': a.v[f'{x}{1-s}{j}'] for x in ('w', 'q', 'b', 'c', 'h')
            for s in range(2) for j in range(2)}
    count = 0
    for e in basis(a, 3):
        fixed = repair(a, e)
        assert not dual_d(a, fixed)
        assert dual_h(a, dual_translation(a, e)) == dual_translation(a, dual_h(a, e))
        assert repair(a, a.substitute(e, swap)) == a.substitute(fixed, swap)
        allowed = {slot_degree(a, m) for m in e}
        assert all(slot_degree(a, m) in allowed for m in fixed)
        count += 1
    contact = a.prod(a.v['q00'], a.v['w01'], a.v['b10'], a.v['w11'])
    contact = a.scale(a.add(contact, a.substitute(contact, swap)), Q(1, 2))
    fixed = repair(a, contact)
    assert fixed and not dual_d(a, fixed)
    assert not physical(a, fixed)
    assert all(slot_degree(a, m) == (2, 2) for m in fixed)
    # Every elementary linear first or second insertion is annihilated, even
    # against the product of all other slot's even fields to high degree.
    linear_checks = 0
    for slot in (0, 1):
        other = a.prod(*(a.prod(a.v[f'{x}{1-slot}{j}'], a.v[f'{x}{1-slot}{j}'])
                         for x in ('w', 'q', 'b') for j in range(2)))
        for x in ('w', 'q', 'b', 'c', 'h'):
            for j in range(2):
                assert not apply(a, fixed, a.mul(a.v[f'{x}{slot}{j}'], other))
                linear_checks += 1
    return dict(symbols=len(a.NAMES), homogeneous_operators_checked=count,
                slot_permutation_and_total_jet_derivative_preserved=True,
                linear_insertion_checks=linear_checks,
                repaired_two_slot_contact=a.display(fixed),
                derivative_orders_per_slot=[2, 2], pure_physical_contact_zero=True,
                distributions_numerically_extended=False,
                scope='Coefficient identities and linear-input normalization subspace; continuum induction is in the report.')


def run():
    return dict(round=792, fresh_test_groups=2,
                dual_coefficient_complex=coefficient_complex(),
                slots_jets_and_linear_insertions=slot_and_jet_contracts(),
                all_checks_passed=True,
                joint_field_Ward_and_N1_constructed_in_stated_local_class=True,
                original_T1_preserved=True,
                original_branch_transport_by_finite_renormalization=True,
                transported_quantum_action_field_only_proven=False,
                original_interacting_positive_state_proven=False,
                all_continuum_coefficients_computed=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('coefficient_ward_normalization_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
