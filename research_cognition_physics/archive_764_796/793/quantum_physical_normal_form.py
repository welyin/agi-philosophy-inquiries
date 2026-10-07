"""793: full BV contraction and positive-loop quantum normal forms.

Exact finite rational diagnostics. The continuum proof uses the renormalized
Ward/L-infinity brackets, NOT this finite dimensional BV Laplacian.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('probe793_nf', HERE/'physical_representative_probe.py')
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)
bv, V, h, project, lap = probe.bv, probe.V, probe.h, probe.project, probe.lap
ORDER = 4


def base():
    r, p, c, a, b, u, v = [V[x] for x in ('r', 'p', 'c', 'a', 'b', 'u', 'v')]
    return bv.add(bv.scale(bv.mul(r, r), Q(3, 4)),
                  bv.scale(bv.prod(r, r, r), Q(5, 42)),
                  bv.mul(p, c), bv.scale(bv.mul(a, b), -1),
                  bv.mul(u, r), bv.scale(bv.prod(v, r, r), Q(1, 2)))


def contraction_in_original_gauge():
    original_probe = probe.run()
    s0 = base()
    q, r, b, hc = [V[x] for x in ('q', 'r', 'b', 'h')]
    psi = bv.mul(hc, bv.add(q, bv.scale(r, Q(2, 5)),
                 bv.scale(bv.mul(q, r), Q(3, 7)),
                 bv.scale(bv.mul(r, r), Q(1, 3)), bv.scale(b, Q(-2, 9))))
    maps = [{x: bv.add(V[x], bv.scale(bv.bracket(V[x], psi), sign))
             for x in bv.NAMES} for sign in (1, -1)]
    transform = lambda f: bv.substitute(f, maps[0])
    inverse = lambda f: bv.substitute(f, maps[1])
    fixed = transform(s0)
    assert fixed == bv.add(s0, bv.bracket(s0, psi))
    assert not bv.bracket(fixed, fixed) and not lap(fixed)
    pairs = 0
    for x in bv.NAMES:
        assert inverse(transform(V[x])) == V[x]
        for y in bv.NAMES:
            assert bv.bracket(transform(V[x]), transform(V[y])) == transform(bv.bracket(V[x], V[y]))
            pairs += 1
    sf = lambda f: bv.bracket(fixed, f)
    hf = lambda f: transform(h(inverse(f)))
    pf = lambda f: transform(project(inverse(f)))
    count = 0
    for degree in range(4):
        for names in combinations_with_replacement(probe.AUX+('r', 't'), degree):
            f = bv.prod(*(V[x] for x in names))
            if f:
                assert bv.add(sf(hf(f)), hf(sf(f))) == bv.add(f, bv.scale(pf(f), -1))
                assert not hf(hf(f))
                for source in ('u', 'v'):
                    assert bv.diff(hf(f), source) == hf(bv.diff(f, source))
                count += 1
    return dict(unfixed_full_bv_probe=original_probe,
                original_gauge_canonical_generator_pairs=pairs,
                conjugated_contraction_monomials=count,
                gauge_fermion=bv.display(psi),
                ghost_and_nonminimal_antifields_all_kept=True,
                physical_nonlinear_equations_and_sources_kept=True,
                action_unchanged_by_proof_chart_roundtrip=(inverse(fixed) == s0),
                scope='Algebraic antifield shift only; no change of propagator, state or physical gauge.')


def series_add(*series):
    return [bv.add(*(s[n] for s in series)) for n in range(ORDER+1)]


def qme(series):
    return [bv.add(*(bv.scale(bv.bracket(series[j], series[n-j]), Q(1, 2))
                     for j in range(n+1)), bv.scale(lap(series[n-1]), -1) if n else {})
            for n in range(ORDER+1)]


def flow(series, generator, loop_order, quantum=True):
    """Exact truncated solution dS/dt=(S,eps^n F)-eps^(n+1) Delta F."""
    out = [{} for _ in range(ORDER+1)]
    for start, poly in enumerate(series):
        term = poly
        for k in range((ORDER-start)//loop_order+1):
            end = start+k*loop_order
            out[end] = bv.add(out[end], bv.scale(term, Q(1, factorial(k))))
            term = bv.bracket(term, generator)
    if quantum:
        term = lap(generator)
        for k in range((ORDER-1)//loop_order):
            end = 1+(k+1)*loop_order
            out[end] = bv.add(out[end], bv.scale(term, Q(-1, factorial(k+1))))
            term = bv.bracket(term, generator)
    return out


def normal_form_recursion():
    q, r, p, t, c, z, b, hc, rho, u, v = [V[x] for x in ('q','r','p','t','c','z','b','h','rho','u','v')]
    physical1 = bv.add(bv.scale(r, Q(7, 13)), bv.scale(bv.mul(r, r), Q(8, 17)),
                      bv.scale(bv.prod(u, r, r), Q(3, 5)),
                      bv.scale(bv.prod(u, v, r, r, r), Q(2, 7)))
    physical2 = bv.add(bv.scale(bv.mul(r, r), Q(3, 11)),
                      bv.scale(bv.prod(u, v, r), Q(4, 23)),
                      bv.scale(bv.prod(v, v, r, r), Q(2, 29)))
    target = [base(), physical1, physical2,
              bv.scale(bv.prod(r, r, r), Q(5, 31)),
              bv.scale(bv.prod(u, v, r, r), Q(7, 37))]
    generator = bv.add(bv.prod(p, q, r), bv.scale(bv.prod(t, q, r, r), Q(2, 3)),
                       bv.scale(bv.prod(hc, r, r), Q(3, 7)),
                       bv.scale(bv.prod(z, c, r), Q(1, 5)),
                       bv.scale(bv.prod(rho, b, r), Q(2, 11)),
                       bv.scale(bv.prod(u, t, q, r), Q(2, 9)), bv.scale(hc, Q(5, 19)))
    original = flow(target, generator, 1)
    assert not any(qme(original))
    assert flow(original, bv.scale(generator, -1), 1) == target
    naive_defect = qme(flow(target, generator, 1, quantum=False))
    assert naive_defect[2]
    current = original
    steps, inverse_steps = [], []
    for n in range(1, ORDER+1):
        discrepancy = current[n]
        # All already fixed lower coefficients are physical QME solutions.
        assert not bv.bracket(base(), discrepancy)
        physical = project(discrepancy)
        primitive = h(discrepancy)
        assert discrepancy == bv.add(physical, bv.bracket(base(), primitive))
        before = current
        current = flow(current, bv.scale(primitive, -1), n)
        assert current[:n] == before[:n] and current[n] == physical
        assert not any(qme(current))
        assert all(project(current[j]) == current[j] for j in range(1, n+1))
        assert flow(current, primitive, n) == before
        inverse_steps.append((n, primitive))
        steps.append(dict(loop=n, input_terms=len(discrepancy),
                          exact_primitive_terms=len(primitive), physical_terms=len(physical),
                          physical_coefficient=bv.display(physical),
                          qme_residual_terms=0, lower_loop_coefficients_unchanged=True))
    reconstructed = current
    for n, primitive in reversed(inverse_steps):
        reconstructed = flow(reconstructed, primitive, n)
    assert reconstructed == original
    # Joint source partners come from the whole series, including UV contacts.
    source_data = []
    for n in range(1, ORDER+1):
        mixed = bv.diff(bv.diff(current[n], 'u'), 'v')
        assert mixed == bv.diff(bv.diff(current[n], 'v'), 'u')
        source_data.append(bv.display(mixed))
    assert source_data[0] == bv.display(bv.scale(bv.prod(r, r, r), Q(2, 7)))
    at_zero = lambda f: bv.substitute(f, {x: {} for x in bv.NAMES})
    assert at_zero(bv.diff(current[1], 'r')) == bv.scale(bv.ONE, Q(7, 13))
    # Check the reduced linear auxiliary generator causes no physical field shift.
    linear_primitive = {m: a for m, a in inverse_steps[0][1].items() if sum(m) == 1}
    assert linear_primitive == bv.scale(hc, Q(5, 19))
    assert not bv.bracket(r, linear_primitive)
    return dict(hbar_order=ORDER, source_variables=['u', 'v'],
                original_action_terms_by_loop=[len(f) for f in original],
                steps=steps, inverse_joint_flow_recovers_original_action=True,
                first_physical_linear_coefficient='7/13',
                linear_auxiliary_generator=bv.display(linear_primitive),
                mixed_source_coefficients=source_data,
                omitted_quantum_contact_qme_defect_order2=bv.display(naive_defect[2]),
                physical_corrections_kept_not_set_to_zero=True,
                scope='Finite BV quantum-flow and normal-form calibration, not continuum anomaly values or a star algebra/state map.')


def run():
    return dict(round=793, fresh_test_groups=2,
                original_gauge_full_bv_contraction=contraction_in_original_gauge(),
                positive_loop_normal_form=normal_form_recursion(), all_checks_passed=True,
                local_original_branch_quantum_normal_form_proven=True,
                physical_finite_terms_and_joint_sources_retained=True,
                original_interacting_star_algebra_dictionary_proven=False,
                original_interacting_positive_state_proven=False,
                original_continuum_anomaly_coefficients_computed=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('quantum_physical_normal_form_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
