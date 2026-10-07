"""796: exact mixed Wick algebra and full-BV mean displacement diagnostics.

The oscillator/CAR algebra is a calibration of star, state and translation
identities. Continuum wavefront/domain claims are proved in the report.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import factorial, comb
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'795'))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


gc = load('gaussian_rational796', HERE.parent/'795/real_cutoff_comparison.py')
probe = load('formal_state796', HERE/'formal_state_probe.py')
C, matrix, scale, dagger = gc.C, gc.matrix, gc.scale, gc.dagger
ZERO, IDENTITY = gc.ZERO, gc.IDENTITY
K = ((C(Q(3, 4)), C(0, Q(1, 2))),
     (C(0, -Q(1, 2)), C(Q(5, 4))))


def cpow(x, n):
    out = C(1)
    for _ in range(n):
        out *= x
    return out


def fall(n, k):
    return factorial(n)//factorial(n-k) if k <= n else 0


def clean(p):
    return {k: v for k, v in p.items() if any(v.flat)}


def add(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out.get(k, ZERO)+v
    return clean(out)


def smul(p, c):
    return clean({k: scale(v, c) for k, v in p.items()})


def star(p, q):
    out = {}
    for (lq, lp, le), x in p.items():
        for (rq, rp, re), y in q.items():
            for a in range(min(lq, rq)+1):
                for b in range(min(lq-a, rp)+1):
                    for c in range(min(lp, rq-a)+1):
                        for d in range(min(lp-c, rp-b)+1):
                            f = Q(fall(lq, a+b)*fall(lp, c+d)*
                                  fall(rq, a+c)*fall(rp, b+d),
                                  factorial(a)*factorial(b)*factorial(c)*factorial(d))
                            weight = C(f)*cpow(K[0][0], a)*cpow(K[0][1], b)
                            weight *= cpow(K[1][0], c)*cpow(K[1][1], d)
                            key = (lq+rq-2*a-b-c, lp+rp-b-c-2*d, le+re)
                            out[key] = out.get(key, ZERO)+scale(x@y, weight)
    return clean(out)


def adj(p):
    return {k: dagger(v) for k, v in p.items()}


def eq(p, q):
    return all(gc.equal(p.get(k, ZERO), q.get(k, ZERO)) for k in p.keys() | q.keys())


def translate(p, hq, hp):
    # q -> q + epsilon*hq, p -> p + epsilon*hp.
    out = {}
    for (a, b, n), v in p.items():
        for i in range(a+1):
            for j in range(b+1):
                key = (i, j, n+a+b-i-j)
                coefficient = Q(comb(a, i)*comb(b, j))*hq**(a-i)*hp**(b-j)
                out[key] = out.get(key, ZERO)+scale(v, coefficient)
    return clean(out)


def expect(p):
    # Positive CAR density, independent of the oscillator.
    out = {}
    for (a, b, n), v in p.items():
        if a == b == 0:
            out[n] = out.get(n, C())+Q(3, 5)*v[0, 0]+Q(2, 5)*v[1, 1]
    return {n: v for n, v in out.items() if v}


def wick_probe():
    one = {(0, 0, 0): IDENTITY}
    q, p = {(1, 0, 0): IDENTITY}, {(0, 1, 0): IDENTITY}
    annihilation = matrix([[0, 1], [0, 0]])
    creation = dagger(annihilation)
    number = creation@annihilation
    assert gc.equal(annihilation@creation+creation@annihilation, IDENTITY)
    assert eq(add(star(q, p), smul(star(p, q), -1)), smul(one, C(0, 1)))
    assert K[0][0].re > 0
    determinant = K[0][0]*K[1][1]-K[0][1]*K[1][0]
    assert determinant == C(Q(11, 16))
    a = {(1, 0, 0): annihilation, (0, 1, 0): scale(annihilation, C(0, 1)),
         (2, 0, 1): creation, (0, 2, 1): scale(creation, -1),
         (1, 1, 0): number, (1, 1, 1): number,
         (0, 0, 2): matrix([[1, Q(2, 3)], [-Q(1, 3), 2]])}
    b = {(1, 0, 0): number, (0, 1, 1): creation,
         (0, 0, 0): annihilation}
    c = {(1, 1, 0): IDENTITY, (0, 0, 1): creation}
    assert eq(star(star(a, b), c), star(a, star(b, c)))
    assert eq(adj(star(a, b)), star(adj(b), adj(a)))
    displacement = (Q(2, 7), -Q(3, 11))
    assert eq(translate(star(a, b), *displacement),
              star(translate(a, *displacement), translate(b, *displacement)))
    assert eq(adj(translate(a, *displacement)), translate(adj(a), *displacement))
    before, after = expect(star(adj(a), a)), expect(star(adj(translate(a, *displacement)),
                                                                translate(a, *displacement)))
    for result in (before, after):
        assert all(not v.im for v in result.values())
        assert result[min(result)].re > 0
    for field, mean, variance in ((q, displacement[0], K[0][0]),
                                  (p, displacement[1], K[1][1])):
        shifted = translate(field, *displacement)
        assert expect(shifted) == {1: C(mean)}
        second = expect(star(shifted, shifted))
        assert second == {0: variance, 2: C(mean*mean)}
    return dict(CCR_exact=True, CAR_exact=True,
                covariance_determinant=str(determinant),
                mixed_star_associativity_and_adjoint=True,
                coherent_translation_star_automorphism=True,
                centered_covariance_unchanged=True,
                norm_before={str(k): str(v) for k, v in sorted(before.items())},
                norm_after={str(k): str(v) for k, v in sorted(after.items())},
                fermion_occupation="2/5",
                continuum_distribution_limit_tested=False)


def displacement_bv_probe():
    bv = load('bv796_displacement', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
    bv.EVEN += ('w', 'epsilon')
    bv.ODD += ('wstar',)
    bv.PAIRS += (('w', 'wstar'),)
    bv.NAMES = bv.EVEN+bv.ODD
    bv.ZERO = (0,)*len(bv.NAMES)
    bv.ONE = {bv.ZERO: Q(1)}
    v = {name: bv.var(name) for name in bv.NAMES}
    q, r, w, b, hc, c = [v[x] for x in ('q', 'r', 'w', 'b', 'h', 'c')]
    physical = bv.scale(bv.mul(bv.add(r, bv.scale(w, -1)),
                               bv.add(r, bv.scale(w, -1))), Q(1, 2))
    ungauged = bv.add(physical, bv.mul(v['p'], c), bv.scale(bv.mul(v['a'], b), -1))
    psi = bv.mul(hc, bv.add(bv.scale(q, 2), r, bv.scale(w, 2), bv.scale(b, Q(3, 2))))
    # Canonical antifield shift; the engine fixes the sign convention.
    shift = {anti: bv.add(v[anti], bv.diff(psi, field))
             for field, anti in bv.PAIRS}
    s2 = bv.substitute(ungauged, shift)
    assert not bv.bracket(s2, s2)
    differential = lambda f: bv.bracket(s2, f)
    translation = {'r': bv.add(r, v['epsilon']),
                   'w': bv.add(w, v['epsilon']),
                   'q': bv.add(q, bv.scale(v['epsilon'], -Q(3, 2)))}
    tau = lambda f: bv.substitute(f, translation)
    assert tau(s2) == s2
    count = 0
    for degree in range(4):
        for names in combinations_with_replacement(tuple(n for n in bv.NAMES if n != 'epsilon'), degree):
            f = bv.prod(*(v[name] for name in names))
            if f:
                assert tau(differential(f)) == differential(tau(f))
                count += 1
    bad = lambda f: bv.substitute(f, {'r': bv.add(r, v['epsilon']),
                                    'w': bv.add(w, v['epsilon'])})
    defects = [bv.add(bad(differential(v[name])), bv.scale(differential(bad(v[name])), -1))
               for name in bv.NAMES if name != 'epsilon']
    assert any(defects)
    return dict(complete_BV_translation_checks=count,
                physical_mean_shift=["1", "1"],
                gauge_partner_mean_shift="-3/2", auxiliary_b_mean_shift="0",
                full_free_action_unchanged=True,
                full_BV_translation_intertwines=True,
                physical_only_shift_rejected=True,
                physical_only_defects=[bv.display(f) for f in defects if f],
                continuum_Cauchy_solvability_tested=False)


def run():
    formal = probe.run()
    assert formal == json.loads((HERE/'formal_state_probe_results.json').read_text(encoding='utf-8'))
    return dict(round=796, fresh_test_groups=3, all_checks_passed=True,
                closed_formal_state=formal, mixed_Wick_state=wick_probe(),
                compatible_mean_displacement=displacement_bv_probe(),
                arithmetic="Exact Fraction / Gaussian rationals",
                scope="Finite algebra calibration only; continuous domain and original local-state existence are addressed analytically in the report.",
                physical_experiment_or_finite_coupling_convergence=False)


if __name__ == '__main__':
    result = run()
    (HERE/'wick_state_initial_data_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
