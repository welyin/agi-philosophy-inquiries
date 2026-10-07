"""784: exact boundary-charge and relative-completion diagnostics.

Integer polynomial matrices test the algebraic inference, not the original
continuum Ward kernels or a Hilbert-space domain. No new physical model is used
in place of the original coupled theory.
"""
from pathlib import Path
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def elementary(n, row, col):
    out = np.zeros((n, n), dtype=np.int64)
    out[row, col] = 1
    return out


def quartet():
    # Ordered basis: a(-1), u(0), v(0), b(+1), p(0 physical).
    a, u, v, b, p = range(5)
    metric = elementary(5, a, b)+elementary(5, b, a)
    metric += elementary(5, u, v)+elementary(5, v, u)+elementary(5, p, p)
    parity = np.diag([-1, 1, 1, -1, 1])
    charge = elementary(5, u, a)+elementary(5, b, v)
    homotopy = elementary(5, a, u)+elementary(5, v, b)
    physical = elementary(5, p, p)
    return metric, parity, charge, homotopy, physical


def multiply_series(left, right):
    out = [np.zeros((left[0].shape[0], right[0].shape[1]), dtype=np.int64)
           for _ in range(len(left)+len(right)-1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j] += a@b
    return out


def assert_series_equal(left, right):
    zero = np.zeros_like(left[0])
    for k in range(max(len(left), len(right))):
        assert np.array_equal(left[k] if k < len(left) else zero,
                              right[k] if k < len(right) else zero), k


def causal_boundary_identity():
    metric, parity, q0, _, _ = quartet()
    r = elementary(5, 4, 0)+elementary(5, 3, 4)
    classical_past = q0
    quantum_past = r
    past = classical_past+quantum_past
    future = 2*q0-r
    z = np.diag([1, 2, 3, 4, 5])
    assert np.array_equal(parity@z, z@parity)
    master = future-past
    # Even observable: future insertions stand left, past insertions right.
    ordered = future@z-z@past
    contact_difference = ordered-master@z
    assert np.array_equal(contact_difference, past@z-z@past)
    after_complete_subtraction = contact_difference-(past@z-z@past)
    after_classical_only = contact_difference-(classical_past@z-z@classical_past)
    assert not np.any(after_complete_subtraction)
    assert np.array_equal(after_classical_only, quantum_past@z-z@quantum_past)
    assert np.any(after_classical_only)
    return dict(graded_causal_split_identity_exact=True,
                complete_boundary_defect_cancels=True,
                omit_quantum_boundary_defect_nonzero_entries=int(np.count_nonzero(after_classical_only)),
                omit_quantum_boundary_defect_squared_frobenius=int(np.sum(after_classical_only**2)),
                continuum_anomaly_or_time_ordering_computed=False)


def relative_completion():
    g, parity, q, h, physical = quartet()
    identity = np.eye(5, dtype=np.int64)
    assert not np.any(q@q)
    assert np.array_equal(q.T@g, g@q)
    assert np.array_equal(q@h+h@q, identity-physical)
    # A second quartet is the unobserved/boundary factor. It is already in Q0.
    q0 = np.kron(q, identity)+np.kron(parity, q)
    metric = np.kron(g, g)
    pi = np.kron(physical, physical)
    contraction = np.kron(h, identity)+np.kron(parity@physical, h)
    assert not np.any(q0@q0)
    assert np.array_equal(q0.T@metric, metric@q0)
    assert np.array_equal(q0@contraction+contraction@q0, np.eye(25, dtype=np.int64)-pi)
    assert np.array_equal(metric@pi, pi)  # Positive physical line, exact null complement.
    r = elementary(5, 4, 0)+elementary(5, 3, 4)
    b = elementary(5, 2, 0)+elementary(5, 3, 1)
    x = elementary(5, 4, 1)-elementary(5, 2, 4)
    assert np.array_equal(r.T@g, g@r)
    assert np.array_equal(b.T@g, g@b)
    assert not np.any(q@r+r@q)
    assert np.array_equal(r@r, elementary(5, 3, 0))
    assert np.array_equal(q@b+b@q, 2*(r@r))
    assert not np.any(x@x@x) and not np.any(x.T@g+g@x)
    # Relative endomorphism contraction is available in this tensor example.
    endomorphism_checks = 0
    for row in range(5):
        for col in range(5):
            a = elementary(5, row, col)
            sign = int(parity[row, row]*parity[col, col])
            da = q@a-sign*a@q
            ha = h@a+sign*physical@a@h
            dha = q@ha+sign*ha@q
            hda = h@da-sign*physical@da@h
            assert np.array_equal(dha+hda, a-physical@a@physical)
            endomorphism_checks += 1
    raw_primitive = h@(4*r@r)+physical@(4*r@r)@h
    # The real part preserves the equation because the curvature is symmetric.
    assert np.array_equal(-(raw_primitive+g@raw_primitive.T@g), -4*b)
    # kappa=2*t removes fractions, while retaining exact polynomial identities.
    charge = [q0, np.kron(parity, 2*r), np.kron(parity, -2*b)]
    raw = charge[:2]
    raw_square = multiply_series(raw, raw)
    assert not np.any(raw_square[0]) and not np.any(raw_square[1])
    assert np.array_equal(raw_square[2], np.kron(identity, 4*r@r))
    assert np.any(raw_square[2])
    for coefficient in multiply_series(charge, charge):
        assert not np.any(coefficient)
    for coefficient in charge:
        assert np.array_equal(coefficient.T@metric, metric@coefficient)
    # All matrix units of the local factor, including odd ones.
    local_checks = 0
    nontrivial_free_actions = 0
    for row in range(5):
        for col in range(5):
            a = elementary(5, row, col)
            sign = int(parity[row, row]*parity[col, col])
            aa = np.kron(a, identity)
            dq = q@a-sign*a@q
            assert np.array_equal(q0@aa-sign*aa@q0, np.kron(dq, identity))
            assert not np.any(q@dq+sign*dq@q)
            for boundary in charge[1:]:
                assert not np.any(boundary@aa-sign*aa@boundary)
            assert not np.any(raw_square[2]@aa-aa@raw_square[2])
            nontrivial_free_actions += int(np.any(dq))
            local_checks += 1
    outside = np.kron(identity, elementary(5, 0, 0))
    assert np.any(raw_square[2]@outside-outside@raw_square[2])
    # Exact Krein-unitary completion entirely inside the relative commutant.
    ub = [identity, 2*x, 2*x@x]
    ubi = [identity, -2*x, 2*x@x]
    assert_series_equal(multiply_series(ub, ubi), [identity])
    assert_series_equal(multiply_series([v.T@g for v in ub], ub), [g])
    lifted_q = multiply_series(multiply_series(ub, [q]), ubi)
    assert_series_equal(lifted_q, [q, 2*r, -2*b])
    # The original physical vector must be lifted, even in the completed model.
    p = np.zeros((5, 1), dtype=np.int64); p[4, 0] = 1
    lifted_vector = [v@p for v in ub]
    assert_series_equal(multiply_series(lifted_q, lifted_vector), [np.zeros((5, 1), dtype=np.int64)])
    assert_series_equal(multiply_series([v.T@g for v in lifted_vector], lifted_vector),
                        [np.ones((1, 1), dtype=np.int64)])
    fixed_exact = elementary(5, 0, 4)  # odd, ghost degree -1
    at_one = sum(lifted_q)
    exact_observable = at_one@fixed_exact+fixed_exact@at_one
    old_expectation = int((p.T@g@exact_observable@p)[0, 0])
    lifted_at_one = sum(lifted_vector)
    new_expectation = int((lifted_at_one.T@g@exact_observable@lifted_at_one)[0, 0])
    assert old_expectation == 2 and new_expectation == 0
    assert np.any(at_one@p)
    # Replacing the global free charge by its spatial/local part is also unsafe.
    negative_boundary = np.zeros((5, 1), dtype=np.int64)
    negative_boundary[1, 0] = 1; negative_boundary[2, 0] = -1
    negative = np.kron(p, negative_boundary)
    local_free_charge = np.kron(q, identity)
    assert not np.any(local_free_charge@negative) and np.any(q0@negative)
    negative_norm = int((negative.T@metric@negative)[0, 0])
    assert negative_norm == -2
    return dict(matrix_dimension=25, exact_integer_polynomial_arithmetic=True,
                free_nilpotency_metric_and_positive_quotient_verified=True,
                local_matrix_units_checked=local_checks,
                nontrivial_local_free_BRST_actions=nontrivial_free_actions,
                raw_charge_square='4*t^2 I_local tensor |b><a|',
                local_derivation_nilpotent_but_raw_charge_not_nilpotent=True,
                curvature_in_relative_commutant_but_not_global_center=True,
                relative_endomorphism_homotopy_units_checked=endomorphism_checks,
                exact_relative_correction='-2*t^2 parity_local tensor (|v><a|+|b><u|)',
                completed_charge_nilpotent_and_Krein_symmetric=True,
                completion_Krein_unitarily_conjugate_to_same_Q0=True,
                old_vector_expectation_of_exact_observable=old_expectation,
                lifted_vector_expectation_of_exact_observable=new_expectation,
                spatially_truncated_free_charge_closed_negative_norm=negative_norm,
                actual_coupled_theory_relative_cohomology_proven=False)


def run():
    return dict(round=784, test_groups=2,
                causal_boundary=causal_boundary_identity(),
                relative_charge=relative_completion(),
                continuum_charge_operator_domain_proven=False,
                original_interacting_positive_state_proven=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'local_charge_boundary_results.json'
    if args.write:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    if args.check:
        assert result == json.loads(target.read_text(encoding='utf-8'))
    print(json.dumps(result, ensure_ascii=False, indent=2))
