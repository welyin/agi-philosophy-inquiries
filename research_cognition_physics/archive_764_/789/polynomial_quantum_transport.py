"""789: fixed-linear-part coordinate transport, with exact finite diagnostics.

The finite BV Laplacian and Gaussian integrals are diagnostic only. The report
uses the ORIGINAL renormalized QAP, not a continuum functional Laplacian.
"""
from fractions import Fraction as Q
from math import factorial, comb
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bv = load('bv777_transport', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
probe = load('coordinate_probe789', HERE/'coordinate_endpoint_probe.py')
V = {name: bv.var(name) for name in bv.NAMES}
WEIGHT = 8


def lap(poly):
    return bv.add(*(bv.scale(bv.diff(bv.diff(poly, anti), field), -1 if field in bv.ODD else 1)
                    for field, anti in bv.PAIRS))


def clean(series):
    out = {}
    for key, poly in series.items():
        loop, parameter, source_a, source_b = key
        p = {m: a for m, a in poly.items()
             if a and sum(m)+2*loop+3*(source_a+source_b) <= WEIGHT}
        if p:
            out[key] = p
    return out


def add(*series):
    out = {}
    for s in series:
        for key, p in s.items():
            out[key] = bv.add(out.get(key, {}), p)
    return clean(out)


def bracket(a, b):
    out = {}
    for i, p in a.items():
        for j, q in b.items():
            key = tuple(x+y for x, y in zip(i, j))
            out[key] = bv.add(out.get(key, {}), bv.bracket(p, q))
    return clean(out)


def defect(s):
    terms = {k: bv.scale(p, Q(1, 2)) for k, p in bracket(s, s).items()}
    delta = {(ell+1, t, a, b): bv.scale(lap(p), -1)
             for (ell, t, a, b), p in s.items()}
    return add(terms, delta)


def at_endpoint(series):
    out = {}
    for (ell, t, a, b), p in series.items():
        key = (ell, 0, a, b)
        out[key] = bv.add(out.get(key, {}), p)
    return clean(out)


def source_derivative(series, source_axis):
    out = {}
    for key, p in series.items():
        if key[source_axis]:
            new = list(key)
            new[source_axis] -= 1
            out[tuple(new)] = bv.scale(p, key[source_axis])
    # No extra truncation: keep the same originally computed source budget.
    return out


def finite_quantum_bv_flow():
    q, r, p, anti_r, c, aa, b = [V[x] for x in ('q', 'r', 'p', 't', 'c', 'a', 'b')]
    mass, cubic, kappa = Q(3, 2), Q(5, 7), Q(2, 3)
    original = bv.add(bv.scale(bv.prod(r, r), mass/2), bv.scale(bv.prod(r, r, r), cubic/6),
                      bv.mul(p, c), bv.scale(bv.mul(aa, b), -1))
    generator = bv.add(bv.scale(bv.prod(p, q, q), Q(1, 2)),
                       bv.scale(bv.prod(anti_r, q, r), kappa),
                       bv.scale(bv.prod(anti_r, r, r), Q(1, 2)))
    delta_g = lap(generator)
    assert delta_g == bv.add(bv.scale(q, 1+kappa), r)
    assert not bv.bracket(original, original) and not lap(original)
    classical = {}
    inputs = ((0, 0, original), (1, 0, r), (0, 1, bv.scale(bv.mul(r, r), Q(1, 2))))
    for sa, sb, term in inputs:
        current = term
        for power in range(WEIGHT):
            classical[(0, power, sa, sb)] = bv.scale(current, Q(1, factorial(power)))
            current = bv.bracket(current, generator)
    classical = clean(classical)
    # Integrate exp((t-u) ad_G) (-hbar Delta G), rather than omitting density.
    quantum = {}
    current = delta_g
    for power in range(1, WEIGHT):
        quantum[(1, power, 0, 0)] = bv.scale(current, Q(-1, factorial(power)))
        current = bv.bracket(current, generator)
    quantum = clean(quantum)
    full = add(classical, quantum)
    assert not defect(full)
    assert not defect(at_endpoint(full))
    wrong = defect(classical)
    assert wrong and (1, 1, 0, 0) in wrong
    # The path equation itself is checked coefficientwise, including sources.
    dt = {(ell, tau-1, a, b): bv.scale(poly, tau)
          for (ell, tau, a, b), poly in full.items() if tau}
    rhs = {key: bv.bracket(poly, generator) for key, poly in full.items()}
    rhs = add(rhs, {(1, 0, 0, 0): bv.scale(delta_g, -1)})
    assert not add(dt, {key: bv.scale(poly, -1) for key, poly in rhs.items()})
    # The repaired action and both source partners share one QME. Differentiating
    # its zero polynomial is stronger than separately picking closed operators.
    assert source_derivative(source_derivative(full, 2), 3) == source_derivative(source_derivative(full, 3), 2)
    original_at_anchor = {k: v for k, v in full.items() if k[1] == 0}
    assert original_at_anchor == clean({(0, 0, a, b): p for a, b, p in inputs})
    return dict(weight_cutoff=WEIGHT, field_weight=1, loop_weight=2, source_weight=3,
                source_menu=['r', 'r^2/2'], action_includes_physical_cubic=True,
                generator=bv.display(generator), quantum_contact_driver=bv.display(delta_g),
                quantum_counterterm_leading=bv.display(quantum[(1, 1, 0, 0)]),
                omitted_contact_leading_defect=bv.display(wrong[(1, 1, 0, 0)]),
                qme_residual_terms=0, endpoint_qme_residual_terms=0,
                path_equation_residual_terms=0, anchored_source_action_preserved=True,
                parameter_source_coefficients=len(full),
                maximum_parameter_degree=max(key[1] for key in full),
                scope='Exact finite BV flow. No original continuous anomaly coefficients computed.')


def double_factorial_odd(k):
    value = 1
    for j in range(1, k+1, 2):
        value *= j
    return value


def s_add(*series):
    out = {}
    for s in series:
        for k, v in s.items():
            out[k] = out.get(k, Q(0))+v
    return {k: v for k, v in out.items() if v}


def s_mul(p, q, order=3):
    out = {}
    for (a, t), v in p.items():
        for (b, u), w in q.items():
            if a+b <= order:
                key = (a+b, t+u)
                out[key] = out.get(key, Q(0))+v*w
    return {k: v for k, v in out.items() if v}


def s_inverse(p, order=3):
    assert p.get((0, 0)) == 1 and not any(a == 0 and t for a, t in p)
    minus_tail = {k: -v for k, v in p.items() if k != (0, 0)}
    out, power = {(0, 0): Q(1)}, {(0, 0): Q(1)}
    for _ in range(order):
        power = s_mul(power, minus_tail, order)
        out = s_add(out, power)
    return out


def finite_source_integral():
    mass, alpha, order = Q(3, 2), Q(2, 3), 3

    def integral(moment, jacobian=True, transform_observable=True):
        # y=x+t*alpha*x^2; m*y^2/2 has cubic and quartic interactions.
        # Exact centered Gaussian moments sum all relevant finite diagrams.
        result = {}
        for extra in range(moment+1 if transform_observable else 1):
            observed_degree = moment+extra
            observed_coefficient = Q(comb(moment, extra))*alpha**extra
            for jac_power in range(2 if jacobian else 1):
                for n3 in range(2*order+1):
                    for n4 in range(order+1):
                        degree = observed_degree+jac_power+3*n3+4*n4
                        if degree % 2:
                            continue
                        pairs = degree//2
                        loop_power = pairs-n3-n4
                        if not 0 <= loop_power <= order:
                            continue
                        coefficient = observed_coefficient*(2*alpha)**jac_power
                        coefficient *= (-mass*alpha)**n3/Q(factorial(n3))
                        coefficient *= (-mass*alpha*alpha/2)**n4/Q(factorial(n4))
                        coefficient *= Q(double_factorial_odd(degree-1), 1)/mass**pairs
                        key = (loop_power, extra+jac_power+n3+2*n4)
                        result[key] = result.get(key, Q(0))+coefficient
        return {k: v for k, v in result.items() if v}

    norm = integral(0)
    assert norm == {(0, 0): Q(1)}
    moments = [integral(i) for i in range(1, 5)]
    assert moments == [{}, {(1, 0): 1/mass}, {}, {(2, 0): 3/mass**2}]
    mixed = s_add(moments[2], {k: -v for k, v in s_mul(moments[0], moments[1]).items()})
    variance_square = s_add(moments[3], {k: -v for k, v in s_mul(moments[1], moments[1]).items()})
    assert not mixed and variance_square == {(2, 0): 2/mass**2}
    no_jacobian = s_mul(integral(1, jacobian=False), s_inverse(integral(0, jacobian=False)))
    linear_observable = integral(1, transform_observable=False)
    assert no_jacobian[(1, 1)] == -2*alpha/mass
    assert linear_observable[(1, 1)] == -alpha/mass
    return dict(expansion_through_hbar=order, mass=str(mass), alpha=str(alpha),
                normalized_moments=['0', 'hbar/m', '0', '3*hbar^2/m^2'],
                connected_mixed_source_covariance='0',
                connected_square_source_variance_coefficient=str(2/mass**2),
                wrong_mean_without_jacobian_hbar_t=str(no_jacobian[(1, 1)]),
                wrong_mean_without_observable_transport_hbar_t=str(linear_observable[(1, 1)]),
                all_positive_parameter_powers_cancel_for_transported_menu=True,
                scope='Formal finite Gaussian change of coordinates; no global invertibility or Lorentzian state asserted.')


def run():
    return dict(round=789, fresh_test_groups=3,
                fixed_linear_coordinate_filtration=probe.run(),
                quantum_bv_action_and_menu=finite_quantum_bv_flow(),
                shared_source_moments=finite_source_integral(),
                all_checks_passed=True,
                formal_polynomial_coordinate_endpoint_repair_proven=True,
                original_continuum_contact_coefficients_computed=False,
                prescribed_788_source_contact_identified_with_endpoint_action=False,
                auxiliary_gauge_endpoint_matching_proven=False,
                original_interacting_positive_state_proven=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('polynomial_quantum_transport_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
