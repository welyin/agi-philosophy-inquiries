"""794: compact-completion algebra and noncommutative Ward homotopies.

The switch parameter test is finite BV, not a spacetime distribution test.
Matrix tests calibrate identities in a differential star algebra, not the
original continuum state or its cutoff-identification map.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


nf = load('nf794_compact', HERE.parent/'793/quantum_physical_normal_form.py')
probe = load('source794_compact', HERE/'source_transport_probe.py')
bv, V, N = nf.bv, nf.V, nf.ORDER


def compact_completion():
    q, r, p, t, c, a, b, hc, g, z = [V[x] for x in ('q','r','p','t','c','a','b','h','u','v')]
    free = bv.add(bv.scale(bv.mul(r, r), Q(3, 4)), bv.mul(p, c), bv.scale(bv.mul(a, b), -1))
    interaction = bv.add(bv.scale(bv.prod(r, r, r), Q(5, 42)), bv.mul(z, r))
    sg = bv.add(free, bv.mul(g, interaction))
    phys1 = bv.add(bv.scale(r, Q(7, 13)), bv.scale(bv.prod(z, z, r, r), Q(2, 7)))
    target = [sg, bv.mul(g, phys1), bv.scale(bv.prod(g, r, r), Q(3, 11)), {}, {}]
    f = bv.add(bv.prod(p, q, r), bv.scale(bv.prod(t, q, r, r), Q(2, 3)),
               bv.scale(bv.prod(hc, r, r), Q(3, 7)), bv.scale(hc, Q(5, 19)))
    completed = nf.flow(target, bv.mul(g, f), 1)
    assert not any(nf.qme(target)) and not any(nf.qme(completed))
    at = lambda s, value: [bv.substitute(x, {'u': bv.scale(bv.ONE, value)}) for x in s]
    full = nf.flow(at(target, Q(1)), f, 1)
    assert at(completed, Q(1)) == full
    assert at(completed, Q(0)) == [free, {}, {}, {}, {}]
    assert nf.flow(completed, bv.scale(bv.mul(g, f), -1), 1) == target
    naive = [sg]+[bv.mul(g, x) for x in full[1:]]
    defect = nf.qme(naive)
    assert any(defect) and not any(at(defect, Q(0))) and not any(at(defect, Q(1)))
    first = next(i for i, x in enumerate(defect) if x)
    middle = at(defect, Q(1, 2))
    assert middle[first]
    return dict(loop_order=N, switched_classical_action_qme_exact=True,
                inverse_local_flow_qme_exact_for_all_switch_values=True,
                platform_value_recovers_full_original_quantum_orbit=True,
                zero_switch_recovers_original_free_action=True,
                inverse_dictionary_recovers_switched_physical_target=True,
                naive_multiply_counterterms_first_defect_loop=first,
                naive_half_switch_defect=bv.display(middle[first]),
                completed_terms_by_loop=[len(x) for x in completed],
                scope='Finite coupling-value diagnostic; spacetime cutoff-jet/support statements are proved analytically in the report.')


def matrix_homotopy():
    old = load('quartet794', HERE.parent/'784/local_charge_boundary.py')
    metric5, parity5, charge5, h5, physical5 = old.quartet()
    eye2 = np.eye(2, dtype=object)
    ident = np.eye(10, dtype=object)
    zero = np.zeros((10, 10), dtype=object)
    metric = np.kron(metric5, eye2)
    charge = np.kron(charge5, eye2)
    h = np.kron(h5, eye2)
    physical = np.kron(physical5, eye2)
    parity = np.kron(parity5, eye2)
    x = np.array([[0, 1], [1, 0]], dtype=object)
    z = np.array([[1, 0], [0, -1]], dtype=object)
    k0 = np.kron(old.elementary(5, 0, 4)+old.elementary(5, 4, 3), x)
    k1 = np.kron(old.elementary(5, 0, 4), z)
    k2 = np.kron(old.elementary(5, 4, 3), x@z)
    d = lambda m, odd=False: charge@m-(-1 if odd else 1)*m@charge
    star = lambda m: metric@m.T@metric
    contract = lambda m: h@m+physical@m@h  # even closed input
    assert not np.any(charge@charge)
    assert np.array_equal(star(charge), charge)
    assert np.array_equal(charge@h+h@charge, ident-physical)
    for k in (k0, k1, k2):
        assert np.array_equal(parity@k@parity, -k)
    # Keys: formal interaction order, one source power, homotopy parameter.
    limit = 3

    def add(*terms):
        out = {}
        for term in terms:
            for key, m in term.items():
                out[key] = out.get(key, zero)+m
        return {key: m for key, m in out.items() if np.any(m)}

    def scale(series, scalar):
        return {key: scalar*m for key, m in series.items() if np.any(scalar*m)}

    def mul(left, right):
        out = {}
        for k, a in left.items():
            for l, b in right.items():
                key = tuple(i+j for i, j in zip(k, l))
                if key[0] <= limit:
                    out[key] = out.get(key, zero)+a@b
        return {key: m for key, m in out.items() if np.any(m)}

    def deriv(series, axis):
        out = {}
        for key, m in series.items():
            if key[axis]:
                reduced = list(key)
                reduced[axis] -= 1
                out[tuple(reduced)] = key[axis]*m
        return out

    def at(series, axis, value):
        out = {}
        for key, m in series.items():
            reduced = list(key)
            reduced[axis] = 0
            key2 = tuple(reduced)
            out[key2] = out.get(key2, zero)+(value**key[axis])*m
        return {key: m for key, m in out.items() if np.any(m)}

    def delta(series, odd=False):
        return {key: d(m, odd) for key, m in series.items() if np.any(d(m, odd))}

    def inverse(series):
        tail = scale(add(series, {(0,0,0): -ident}), -1)
        out = power = {(0,0,0): ident}
        for _ in range(limit):
            power = mul(power, tail)
            out = add(out, power)
        assert not add(mul(out, series), {(0,0,0): -ident})
        return out

    # Nontrivial physical matrix source plus an exact, source-dependent orbit.
    base = {(0,0,0): ident, (1,0,0): np.kron(physical5, z),
            (1,1,0): np.kron(physical5, x), (1,2,0): np.kron(physical5, z)}
    primitive = {(1,0,1): k0, (1,1,1): k1, (2,2,2): k2}
    st = add(base, delta(primitive, odd=True))
    assert not delta(st)
    inv0 = inverse(at(st, 1, Q(0)))
    relative = mul(inv0, st)
    assert not delta(relative)
    dotk = deriv(primitive, 2)
    # D is odd; inverse and both S factors are even.
    primitive_relative = add(
        scale(mul(mul(mul(inv0, at(dotk, 1, Q(0))), inv0), st), -1),
        mul(inv0, dotk))
    assert not add(deriv(relative, 2), scale(delta(primitive_relative, odd=True), -1))
    before, after = at(relative, 2, Q(0)), at(relative, 2, Q(1))
    difference = add(after, scale(before, -1))
    assert difference
    for m in difference.values():
        assert not np.any(d(m))
        assert not np.any(physical@m@physical)
        assert np.array_equal(d(contract(m), odd=True), m)
        assert np.array_equal(d(star(contract(m)), odd=True), star(m))
    first = at(deriv(relative, 1), 1, Q(0))
    second = scale(at(deriv(deriv(relative, 1), 1), 1, Q(0)), Q(1, 2))
    products = [mul(first, second), mul(second, first)]
    noncomm = add(products[0], scale(products[1], -1))
    assert any(np.any(physical@m@physical) for m in noncomm.values())
    product_coefficients = 0
    for product in products:
        change = add(at(product, 2, Q(1)), scale(at(product, 2, Q(0)), -1))
        for m in change.values():
            assert np.array_equal(d(contract(m), odd=True), m)
            product_coefficients += 1
    # An arbitrary cutoff-comparison similarity is not automatically a star map.
    u_bad = np.diag([Q(2), Q(1)]).astype(object)
    ui_bad = np.diag([Q(1, 2), Q(1)]).astype(object)
    image_x = ui_bad@x@u_bad
    assert np.any(image_x-image_x.T)
    imaginary_expectation = (image_x[0,1]-image_x[1,0])/2
    assert imaginary_expectation == Q(-3, 4)
    return dict(quartet_times_physical_dimension=10, physical_matrix_dimension=2,
                formal_order=limit, inverse_and_source_homotopy_exact=True,
                relative_S_homotopy_primitive_coefficients=len(primitive_relative),
                endpoint_exact_difference_coefficients=len(difference),
                endpoint_difference_nonzero_but_exact=True,
                conjugates_of_differences_exact=True,
                ordered_product_difference_coefficients=product_coefficients,
                physical_source_commutator_nonzero=True,
                nonunitary_similarity_selfadjointness_defect=True,
                nonunitary_similarity_X_expectation_in_vector_1_i_over_sqrt2=str(imaginary_expectation)+'*i',
                scope='Noncommutative differential star-algebra identity calibration; no original continuum state/charge-domain assertion.')


def run():
    return dict(round=794, fresh_test_groups=3,
                source_connection_and_insertion_contact=probe.run(),
                compatible_compact_completion=compact_completion(),
                relative_S_and_noncommutative_products=matrix_homotopy(),
                all_checks_passed=True,
                compact_QME_completion_of_original_local_germ_proven=True,
                constructed_completion_relative_S_cohomology_identity_proven=True,
                original_arbitrary_cutoff_star_dictionary_proven=False,
                original_interacting_positive_state_proven=False,
                global_or_UV_result=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('compact_ward_homotopy_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
