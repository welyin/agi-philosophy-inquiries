"""Exact radial identities and limited source/signal matching diagnostics."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent


def clean(p):
    return {k: v for k, v in p.items() if v}


def plus(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, F(0)) + v
    return clean(out)


def scale(p, a):
    return clean({k: a*v for k, v in p.items()})


def shift(p, n):
    return {k+n: v for k, v in p.items()}


def diff(p):
    return clean({k-1: k*v for k, v in p.items()})


def lap(p):
    return plus(diff(diff(p)), scale(shift(diff(p), -1), 2))


def at(p, r):
    return sum((v*r**k for k, v in p.items()), F(0))


def kernel(w, u):
    # Hess W=(W'/r) delta +(W''/r^2-W'/r^3) x_i x_j.
    a = plus(u, scale(shift(diff(w), -1), -1))
    b = plus(scale(shift(diff(diff(w)), -2), -1), shift(diff(w), -3))
    return a, b


def divergence_coefficient(a, b):
    return plus(shift(diff(a), -1), shift(diff(b), 1), scale(b, 4))


def calculate():
    # Units a=1; homogeneity restores all radius powers in the note.
    ui, uo = {0:F(3,2), 2:F(-1,2)}, {-1:F(1)}
    wi = {0:F(3,4), 2:F(1,2), 4:F(-1,20)}
    wo = {1:F(1), -1:F(1,5)}
    assert lap(wi) == scale(ui, 2)
    assert lap(wo) == scale(uo, 2)
    boundary = []
    di, do = wi.copy(), wo.copy()
    for order in range(4):
        assert at(di,F(1)) == at(do,F(1))
        boundary.append(str(at(di,F(1))))
        di, do = diff(di), diff(do)
    ai, bi = kernel(wi, ui)
    ao, bo = kernel(wo, uo)
    assert ai == {0:F(1,2), 2:F(-3,10)} and bi == {0:F(2,5)}
    assert ao == {-3:F(1,5)} and bo == {-3:F(1), -5:F(-3,5)}
    for a,b,u in [(ai,bi,ui),(ao,bo,uo)]:
        assert plus(scale(a,3),shift(b,2)) == u
        assert divergence_coefficient(a,b) == scale(shift(diff(u),-1),-1)
    original_trace_defect = plus({2:F(1)},scale(ui,-1))
    original_div_defect = plus({0:F(4)}, shift(diff(ui),-1))
    assert original_trace_defect == {0:F(-3,2),2:F(3,2)}
    assert original_div_defect == {0:F(3)}
    # External correction eigenvalues relative to point-kernel norm: +a²/5r² twice, -2a²/5r² once.
    ratio = 2*F(1,100)**2/5
    assert ratio == F(1,25000)
    # Exact slow-mass bound samples where E/m is rational; analytic all-z proof is in the note.
    examples = []
    for q in [F(1),F(10001,10000),F(3,2),F(2),F(10)]:
        z = q*q-1
        f = (1+2*z)/q
        error = 1+3*z/2-f
        upper = 5*z*z/8
        assert 0 <= error <= upper
        examples.append(dict(E_over_m=str(q),z=str(z),error_over_m=str(error),upper_over_m=str(upper)))
    # Equilateral opposite-edge primitive: with L=1, h²=3/4, endpoint distance=1.
    assert F(1,4)+F(3,4) == 1
    # Integral factor = -m/h, curvature factor = 2m/h = 4m/(sqrt(3)L).
    return dict(status='PASS_SOURCE_MATCHING_DIAGNOSTICS_ONLY',
        original_uniform_source_identification='REJECTED_BY_TRACE_AND_DIVERGENCE',
        uniform_sphere=dict(W_boundary_derivatives_a1=boundary,
            exact_trace_and_divergence=True,
            outer_relative_kernel_correction_at_a_over_R_1percent=str(ratio)),
        optical_triangle=dict(first_order_equilateral_coefficient='4 G sum(m)/(sqrt(3) L)',
            actual_geodesic_error_bound=False),
        test_source_signal=dict(massless_coefficient='2 |p|', slow_mass_bound='5 p^4/(8 m^3)',
            exact_rational_samples=examples, finite_recoil_or_quantum_certificate=False),
        scope=dict(new_scientific_rounds=0, original1260_modified=False,
            actual_terminal_instrument=False, complete_source_process=False,
            CT_adopted=False, goal_completed=False))


def main():
    result = calculate()
    result['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['note_sha256'] = hashlib.sha256((HERE/'README.md').read_bytes()).hexdigest()
    path=HERE/'results.json'
    if sys.argv[1:] == ['--write-once']:
        with path.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('WROTE_SOURCE_MATCHING_DIRECTION_RESULT')
    else:
        assert not sys.argv[1:]
        assert json.loads(path.read_text('utf-8')) == result
        print('READ_ONLY_MATCH_SOURCE_MATCHING_DIRECTION_RESULT')


if __name__ == '__main__':
    main()
