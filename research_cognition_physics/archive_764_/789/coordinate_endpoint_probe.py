"""Unpublished 789: polynomial parameter degree for a nonlinear coordinate path.

This probes a possible route beyond 788, not a quantum matching theorem.
All series remain formal in field degree. No singular gauge endpoint is solved.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

ORDER = 9
X = {(1, 0): Q(1)}
ONE = {(0, 0): Q(1)}


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, Q(0))+v
    return {k: v for k, v in out.items() if v}


def scale(p, a):
    return {k: a*v for k, v in p.items() if a*v}


def mul(p, q):
    out = {}
    for (a, b), v in p.items():
        for (c, d), w in q.items():
            if a+c <= ORDER:
                key = (a+c, b+d)
                out[key] = out.get(key, Q(0))+v*w
    return {k: v for k, v in out.items() if v}


def power(p, n):
    out = ONE
    for _ in range(n):
        out = mul(out, p)
    return out


def times_t(p):
    return {(a, b+1): v for (a, b), v in p.items()}


def degree(p, n):
    return {k: v for k, v in p.items() if k[0] == n}


def endpoint(p):
    out = {}
    for (a, b), v in p.items():
        out[a] = out.get(a, Q(0))+v
    return {str(a): str(v) for a, v in out.items() if v}


def run():
    a, b, mass, cubic = Q(2, 3), Q(-1, 5), Q(3, 2), Q(5, 7)

    def q(h):
        return add(scale(power(h, 2), a), scale(power(h, 3), b))

    h = X
    for n in range(2, ORDER+1):
        error = degree(add(h, times_t(q(h)), scale(X, -1)), n)
        h = add(h, scale(error, -1))
    assert not add(h, times_t(q(h)), scale(X, -1))
    action = add(scale(power(h, 2), mass/2), scale(power(h, 3), cubic/6))
    generator = q(h)  # Nonlinear field part; its antifield adds one field slot.
    log_derivative = {}
    jac_defect = times_t(add(scale(h, 2*a), scale(power(h, 2), 3*b)))
    for k in range(1, ORDER+1):
        log_derivative = add(log_derivative, scale(power(jac_defect, k), Q((-1)**(k+1), k)))
    rows = []
    for n in range(2, ORDER+1):
        inverse_deg = max(k[1] for k in degree(h, n))
        action_deg = max(k[1] for k in degree(action, n))
        generator_deg = max(k[1] for k in degree(generator, n))
        assert inverse_deg <= n-1 and action_deg <= n-2 and generator_deg <= n-2
        rows.append(dict(field_degree=n, inverse_parameter_degree=inverse_deg,
                         action_parameter_degree=action_deg,
                         generator_field_part_parameter_degree=generator_deg))
    assert all(b <= a for a, b in log_derivative)
    # Arbitrarily long chains of quadratic vertices have E=2, L=0.
    chains = [dict(vertices=k, internal_edges=k-1, external_legs=2,
                   loops=(k-1)-k+1, parameter_degree=k) for k in (1, 2, 4, 8, 16)]
    return dict(round=789, status='working_probe_not_signed_off', field_order=ORDER,
                nonlinear_inverse_residual='0', fixed_linear_part=True,
                polynomial_degree_checks=rows,
                inverse_at_parameter_one=endpoint(h),
                finite_jacobian_log_at_parameter_one=endpoint(log_derivative),
                varying_quadratic_part_negative_control=chains,
                endpoint_evaluation_only_for_finite_parameter_polynomials=True,
                continuous_quantum_coordinate_matching_proven=False,
                auxiliary_gauge_endpoint_matching_proven=False,
                original_interacting_positive_state_proven=False)


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('coordinate_endpoint_probe_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
