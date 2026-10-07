"""Working 796: formal state after cohomological observable transport.

Reuses the old quartet. New diagnostics concern a zero leading physical vector,
higher-order positivity and the need to retain both cross terms in a square.
This is not a continuum Wick-extension or interacting-state existence proof.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('quartet796', HERE.parent/'784/local_charge_boundary.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def run():
    m5, _, d5, h5, p5 = old.quartet()
    identity2 = np.eye(2, dtype=object)
    metric = np.kron(m5, identity2).astype(object)
    d = np.kron(d5, identity2).astype(object)
    p = np.kron(p5, identity2).astype(object)
    h = np.kron(h5, identity2).astype(object)
    zero = np.zeros((10, 10), dtype=object)
    assert np.array_equal(d@h+h@d, np.eye(10, dtype=object)-p)
    delta = lambda a, odd=False: d@a-(-1 if odd else 1)*a@d
    star = lambda a: metric@a.T@metric
    physical = lambda a: a[8:10, 8:10]
    density = np.diag([Q(2, 5), Q(3, 5)]).astype(object)
    value = lambda a: np.trace(density@physical(a))
    x = np.array([[0, 1], [1, 0]], dtype=object)
    z = np.array([[1, 0], [0, -1]], dtype=object)
    yreal = x@z
    lift = lambda a: np.kron(p5, a).astype(object)
    k0 = np.kron(old.elementary(5, 0, 4)+old.elementary(5, 4, 3), x).astype(object)
    k1 = np.kron(old.elementary(5, 0, 2), z).astype(object)
    k2 = np.kron(old.elementary(5, 2, 3), x).astype(object)
    exact = [delta(k, True) for k in (k0, k1, k2)]
    a = [exact[0], lift(x)+exact[1], lift(z)+exact[2], lift(yreal)+exact[0]]
    b = [physical(coefficient) for coefficient in a]
    assert np.any(a[0]) and not np.any(b[0])
    for coefficient in a:
        assert not np.any(delta(coefficient))
        assert np.array_equal(physical(star(coefficient)), physical(coefficient).T)
    for left in a:
        for right in a:
            assert np.array_equal(physical(left@right), physical(left)@physical(right))
    coefficients = []
    for n in range(2*len(a)-1):
        square = zero.copy()
        expected = np.zeros((2, 2), dtype=object)
        for i in range(len(a)):
            j = n-i
            if 0 <= j < len(a):
                square += star(a[i])@a[j]
                expected += b[i].T@b[j]
        assert np.array_equal(physical(square), expected)
        assert value(square) == np.trace(density@expected)
        coefficients.append(value(square))
    first = next(n for n, c in enumerate(coefficients) if c)
    assert first == 2 and coefficients[first] > 0
    # Outside ker D, the same extended functional is genuinely indefinite.
    bad = np.kron(old.elementary(5, 1, 4)-old.elementary(5, 2, 4), identity2)
    assert np.any(delta(bad)) and value(star(bad)@bad) == -2
    # A finite truncation of a formal norm is not a nonnegative numerical law.
    # ||(1-epsilon) psi||^2 = 1-2 epsilon+epsilon^2; dropping the last term fails.
    assert 1-2*Q(1) == -1 and (1-Q(1))**2 == 0
    return dict(working_round=796, all_checks_passed=True,
                exact_leading_observable_nonzero=True,
                physical_leading_observable_zero=True,
                closed_coefficient_products_checked=16,
                norm_coefficients=[str(c) for c in coefficients],
                first_nonzero_norm_order=first,
                first_nonzero_norm_coefficient=str(coefficients[first]),
                negative_square_outside_closed_algebra="-2",
                positivity_requires_closed_observables=True,
                formal_positivity_is_not_truncated_numeric_positivity=True,
                original_continuum_Wick_state_extension_proven=False,
                original_interacting_state_constructed=False,
                formal_round_completed=False)


if __name__ == '__main__':
    result = run()
    (HERE/'formal_state_probe_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
