"""797: exact diagnostics of the source-algebra resolvent instrument.
Finite matrices calibrate Ward, cohomology, source transport and record order.
They do not implement an internal detector in the original continuum theory.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'795'))
spec = importlib.util.spec_from_file_location('rational797', HERE.parent/'795/real_cutoff_comparison.py')
old = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = old
spec.loader.exec_module(old)
C, matrix, scale = old.C, old.matrix, old.scale
I, Z0 = old.IDENTITY, old.ZERO
X = matrix([[0, 1], [1, 0]])
Y = matrix([[0, C(0, -1)], [C(0, 1), 0]])
Z = matrix([[1, 0], [0, -1]])
MAX_H, MAX_T = 2, 4

def dag(a):
    return np.array([[a[j, i].conjugate() for j in range(a.shape[0])]
                     for i in range(a.shape[1])], dtype=object)

def same(a, b):
    return not any((a-b).flat)

def clean(p):
    return {k: v for k, v in p.items() if any(v.flat)}

def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, Z0)+v
    return clean(out)

def smul(p, z):
    return clean({k: scale(v, z) for k, v in p.items()})

def mul(p, q):
    out = {}
    for (h, t), a in p.items():
        for (k, u), b in q.items():
            key = h+k, t+u
            if key[0] <= MAX_H and key[1] <= MAX_T:
                out[key] = out.get(key, Z0)+a@b
    return clean(out)

def adj(p):
    return {k: dag(v) for k, v in p.items()}

def eq(p, q):
    return all(same(p.get(k, Z0), q.get(k, Z0)) for k in p.keys()|q.keys())

ONE = {(0, 0): I}

def instrument(a):
    x = {(h, t+1): scale(v, C(0, 1)) for (h, t), v in a.items()
         if t < MAX_T}
    out, power = ONE, ONE
    for _ in range(MAX_T):
        power = smul(mul(power, x), -1)
        out = add(out, power)
    return out, add(ONE, smul(out, -1))

def operation(ks, b):
    return add(*(mul(mul(adj(k), b), k) for k in ks))

def trace_series(p, density):
    return {f'{h},{t}': str(np.trace(density@v)) for (h, t), v in sorted(p.items())
            if np.trace(density@v)}

def coeff(p, h, t):
    return p.get((h, t), Z0)

def inv_complex(z):
    d = z.re*z.re+z.im*z.im
    assert d
    return C(z.re/d, -z.im/d)

def inv2(a):
    d = inv_complex(a[0, 0]*a[1, 1]-a[0, 1]*a[1, 0])
    return scale(matrix([[a[1, 1], -a[0, 1]], [-a[1, 0], a[0, 0]]]), d)

def exact_instrument(a, t):
    k0 = inv2(I+scale(a, C(0, t)))
    return k0, I-k0

def value(rho, a):
    v = np.trace(rho@a)
    assert v.im == 0
    return v.re

def source_history_probe():
    # Noncommuting higher-order source correction is intentionally retained.
    a = {(0, 0): Z, (1, 0): X, (2, 0): scale(I+Y, Q(1, 3))}
    b = {(0, 0): X+scale(Z, Q(1, 2)), (1, 0): scale(Y, Q(2, 3))}
    ka, kb = instrument(a), instrument(b)
    assert eq(operation(ka, ONE), ONE)
    assert eq(operation(kb, ONE), ONE)
    # One common source-independent unitary comparison, through h^2.
    u = {(0, 0): I, (1, 0): scale(Y, C(0, 1)), (2, 0): scale(I, -Q(1, 2))}
    assert eq(mul(adj(u), u), ONE)
    alpha = lambda p: mul(mul(adj(u), p), u)
    for p, ks in ((a, ka), (b, kb)):
        transported = instrument(alpha(p))
        assert all(eq(alpha(k), kt) for k, kt in zip(ks, transported))
    histories = [mul(second, first) for first in ka for second in kb]
    assert eq(operation(histories, ONE), ONE)
    # Matrix-amplified Gram identity, with noncommuting test entries.
    entries = (a, b)
    test = (ONE, {(0, 0): Y, (1, 0): X})
    left, right = {}, {}
    for k in ka:
        v = add(*(mul(mul(entries[j], k), test[j]) for j in range(2)))
        right = add(right, mul(adj(v), v))
    for i in range(2):
        for j in range(2):
            gij = operation(ka, mul(adj(entries[i]), entries[j]))
            left = add(left, mul(mul(adj(test[i]), gij), test[j]))
    assert eq(left, right)
    # Transport both state and observable; leaving state fixed is not equivalent.
    rho = scale(I+scale(Y, Q(1, 3))+scale(Z, Q(1, 5)), Q(1, 2))
    rho_new = alpha({(0, 0): rho})
    for d in histories:
        e = mul(adj(d), d)
        assert eq(alpha(e), mul(adj(alpha(d)), alpha(d)))
        assert trace_series(mul(rho_new, alpha(e)), I) == trace_series(e, rho)
    untransported = add(alpha(b), smul(b, -1))
    assert trace_series(untransported, rho)
    # Removing source corrections makes an actual, nonzero h-dependent change.
    bare = instrument({(0, 0): Z})
    defect = add(mul(adj(ka[1]), ka[1]), smul(mul(adj(bare[1]), bare[1]), -1))
    assert any(h for h, t in defect)
    # First nonselective terms: i[A,B] and 2ABA - {A^2,B}.
    a0, b0 = coeff(a, 0, 0), coeff(b, 0, 0)
    channel = operation(ka, b)
    assert same(coeff(channel, 0, 1), scale(a0@b0-b0@a0, C(0, 1)))
    assert same(coeff(channel, 0, 2), scale(a0@b0@a0, 2)-a0@a0@b0-b0@a0@a0)
    # A normalized finite-dimensional law tests history order independently.
    aa = matrix([[2, 1], [1, 0]])
    bb = matrix([[1, C(0, -1)], [C(0, 1), -2]])
    e1, e2 = exact_instrument(aa, Q(2, 3)), exact_instrument(bb, Q(3, 5))
    forward, backward = [], []
    for first in e1:
        for second in e2:
            d, dr = second@first, first@second
            forward.append(value(rho, dag(d)@d))
            backward.append(value(rho, dag(dr)@dr))
    assert sum(forward) == sum(backward) == 1
    assert all(p > 0 for p in forward+backward)
    assert forward != backward
    return dict(formal_orders=dict(h=MAX_H, read_strength=MAX_T),
                two_instruments_and_four_histories_normalized=True,
                matrix_amplified_Gram_identity=True,
                common_cutoff_instrument_transport=True,
                common_cutoff_record_state_transport=True,
                untransported_state_source_difference=trace_series(untransported, rho),
                nonselective_first_two_coefficients_verified=True,
                omitted_quantum_source_effect=trace_series(defect, rho),
                exact_two_record_probabilities=[str(p) for p in forward],
                reverse_order_probabilities=[str(p) for p in backward],
                finite_probabilities_are_calibration_only=True)

def ward_probe():
    spec = importlib.util.spec_from_file_location('quartet797', HERE.parent/'784/local_charge_boundary.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    m5, _, d5, h5, p5 = mod.quartet()
    to_c = lambda a: np.array([[old.cast(x) for x in row] for row in a], dtype=object)
    m = to_c(np.kron(m5, np.eye(2, dtype=object)))
    d = to_c(np.kron(d5, np.eye(2, dtype=object)))
    id10 = to_c(np.eye(10, dtype=object))
    z10 = scale(id10, 0)
    delta = lambda a, odd=False: d@a-(-1 if odd else 1)*a@d
    sharp = lambda a: m@dag(a)@m
    physical = lambda a: a[8:10, 8:10]
    g = to_c(np.kron(mod.elementary(5, 0, 4)+mod.elementary(5, 4, 3), X))
    exact = delta(g, True)
    exact = exact+sharp(exact)
    aa = to_c(np.kron(p5, Z+scale(X, Q(1, 3))))+exact
    assert any(exact.flat) and not any(physical(exact).flat)
    assert same(sharp(aa), aa) and not any(delta(aa).flat)
    order = 4
    power, k0 = id10, [id10]
    for n in range(1, order+1):
        power = scale(power@aa, C(0, -1))
        k0.append(power)
    k1 = [z10]+[scale(c, -1) for c in k0[1:]]
    identities = 0
    for n in range(order+1):
        norm = z10.copy()
        post_exact = z10.copy()
        for ks in (k0, k1):
            for i in range(n+1):
                norm += sharp(ks[i])@ks[n-i]
                post_exact += sharp(ks[i])@exact@ks[n-i]
            assert not any(delta(ks[n]).flat)
        assert same(norm, id10 if n == 0 else z10)
        assert not any(physical(post_exact).flat)
        identities += 1
    # Nonclosed raw readout fails Ward already at first reading order.
    raw = aa+to_c(np.kron(mod.elementary(5, 1, 4)-mod.elementary(5, 2, 4), I))
    violation = delta(raw)
    assert any(violation.flat)
    return dict(nonzero_exact_source_representative=True,
                resolvent_coefficients_closed=True,
                complete_BV_normalization_orders_checked=identities,
                exact_insertions_still_annihilated=True,
                raw_nonclosed_readout_Ward_defect_entries=sum(bool(x) for x in violation.flat))

def projection_bridge_probe():
    n = matrix([[0, 0], [0, 1]])
    rho = scale(I+scale(X, Q(1, 3))+scale(Z, Q(1, 5)), Q(1, 2))
    rows = []
    for t in (Q(1, 2), Q(1), Q(2), Q(7)):
        k0, k1 = exact_instrument(n, t)
        assert same(dag(k0)@k0+dag(k1)@k1, I)
        assert same(dag(k1)@k1, scale(n, t*t/(1+t*t)))
        rows.append(dict(t=str(t), outcome1=str(value(rho, dag(k1)@k1)),
                         retained_coherence_factor=str(inv_complex(C(1, -t)))))
    # A four-block unitary is an algebraic dilation, not an apparatus model.
    k0, k1 = exact_instrument(n, Q(2))
    v = np.block([[k0, -dag(k1)], [k1, dag(k0)]])
    id4 = matrix([[1 if i == j else 0 for j in range(4)] for i in range(4)])
    assert same(dag(v)@v, id4) and same(v@dag(v), id4)
    # The labels 0,1 reveal |n| but a general A's effects do not distinguish +/-A.
    a = matrix([[2, 1], [1, -1]])
    kp, km = exact_instrument(a, Q(2, 3)), exact_instrument(-a, Q(2, 3))
    assert all(same(dag(p)@p, dag(q)@q) for p, q in zip(kp, km))
    assert not same(sum((dag(k)@X@k for k in kp), Z0),
                    sum((dag(k)@X@k for k in km), Z0))
    return dict(original_free_CAR_projection_formula_recovered=True,
                projection_readout_calibrations=rows,
                algebraic_dilation_unitary=True,
                dilation_is_original_internal_apparatus=False,
                effect_loses_sign_but_backaction_retains_it=True,
                strong_readout_limit_used_in_formal_theory=False)

def run():
    return dict(round=797, all_checks_passed=True, fresh_test_groups=2,
                source_records_and_Ward=dict(source=source_history_probe(), ward=ward_probe()),
                original_CAR_interface=projection_bridge_probe(),
                original_continuum_proof_is_in_report=True,
                finite_probes_prove_continuum_or_finite_coupling=False,
                internal_preparation_or_autonomous_detector_implemented=False)

if __name__ == '__main__':
    result = run()
    (HERE/'source_record_instrument_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
