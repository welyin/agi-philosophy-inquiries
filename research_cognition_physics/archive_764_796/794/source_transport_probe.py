"""Working 794: source-dependent quantum flows and insertion products.

The finite Delta calibrates contacts; it does not define a continuum product
or prove the original interacting star algebra/state correspondence.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('nf794_probe', HERE.parent/'793/quantum_physical_normal_form.py')
nf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nf)
bv, V, lap, N = nf.bv, nf.V, nf.lap, nf.ORDER


def q_action(action, insertion):
    return [bv.add(*(bv.bracket(action[j], insertion[n-j]) for j in range(n+1)),
                   bv.scale(lap(insertion[n-1]), -1) if n else {})
            for n in range(N+1)]


def sd(series, source):
    return [bv.diff(f, source) for f in series]


def connection(generator, source):
    term = bv.diff(generator, source)
    result = [{} for _ in range(N+1)]
    for k in range(N):
        result[k+1] = bv.scale(term, Q(1, factorial(k+1)))
        term = bv.bracket(term, generator)
    return result


def run():
    q, r, p, t, hc, u, v = [V[x] for x in ('q', 'r', 'p', 't', 'h', 'u', 'v')]
    physical = bv.add(bv.scale(r, Q(7, 13)), bv.scale(bv.prod(u, v, r, r, r), Q(2, 7)))
    start = [nf.base(), physical]+[{} for _ in range(N-1)]
    generator = bv.add(bv.scale(bv.prod(p, q, q), Q(1, 2)),
                       bv.prod(t, q, r), bv.prod(u, t, q, r),
                       bv.scale(bv.prod(v, hc, r, r), Q(1, 3)),
                       bv.scale(bv.prod(u, v, t, q), Q(2, 5)))
    endpoint = nf.flow(start, generator, 1)
    assert not any(nf.qme(endpoint))
    partners = {}
    for source in ('u', 'v'):
        actual = sd(endpoint, source)
        homogeneous = nf.flow(sd(start, source), generator, 1, quantum=False)
        correction = q_action(endpoint, connection(generator, source))
        assert actual == nf.series_add(homogeneous, correction)
        assert any(correction)
        assert not any(q_action(endpoint, actual))
        partners[source] = dict(derivative_equals_homogeneous_plus_q_exact_connection=True,
                               quantum_closed_through_order=N,
                               leading_omitted_connection=bv.display(next(f for f in correction if f)),
                               nonzero_connection_terms_by_loop=[len(f) for f in correction])
    assert sd(sd(endpoint, 'u'), 'v') == sd(sd(endpoint, 'v'), 'u')
    # This is the BV insertion product, NOT the physical interacting star product.
    base = [nf.base()]+[{} for _ in range(N)]
    f = [r]+[{} for _ in range(N)]
    primitive = [bv.mul(t, q)]+[{} for _ in range(N)]
    g = q_action(base, primitive)
    assert not any(q_action(base, f)) and not any(q_action(base, g))
    bare = [bv.mul(r, term) for term in g]
    defect = q_action(base, bare)
    assert defect[1] == V['c'] and not any(defect[j] for j in (0, 2, 3, 4))
    corrected = q_action(base, [bv.prod(r, t, q)]+[{} for _ in range(N)])
    contact = nf.series_add(corrected, [bv.scale(term, -1) for term in bare])
    assert contact[1] == bv.scale(q, -1)
    assert not any(q_action(base, corrected))
    return dict(working_round=794, hbar_order=N,
                joint_source_transport=partners,
                mixed_source_derivatives_commute=True,
                closed_times_exact_bare_product_defect='hbar*c',
                restored_exact_insertion_contact='-hbar*q',
                all_checks_passed=True,
                scope='Source-connection and bare insertion-product diagnostics only; original star algebra, boundary completion, conjugation and state remain unproved.',
                original_interacting_star_dictionary_proven=False,
                original_interacting_positive_state_proven=False,
                formal_round_completed=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('source_transport_probe_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
