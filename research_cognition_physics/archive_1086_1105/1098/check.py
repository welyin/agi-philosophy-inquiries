"""Finite algebra checks for the conditional observer/influence bridge.

No entries certify actual apparatus, observer preparation or an axiom source.
The default run compares saved results without modifying files.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import comb
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def coefficients(n):
    # chi0 = log(2). Exact values, independent of floating hyperbolic functions.
    r = F(2**n)
    return (r + 1/r)/2, (r - 1/r)/2


def calculate():
    groups = []

    def add(name, condition, **data):
        assert bool(condition), name
        groups.append(dict(name=name, status='PASS', **data))

    b = np.eye(4)
    b[:2, :2] = [[1.25, -.75], [-.75, 1.25]]
    rows = []
    maximum_matrix_error = 0.
    for n in range(1, 7):
        ch, sh = coefficients(n)
        assert ch*ch - sh*sh == 1
        bn = np.linalg.matrix_power(b, n)
        target = np.eye(4)
        target[:2, :2] = [[float(ch), -float(sh)], [-float(sh), float(ch)]]
        maximum_matrix_error = max(maximum_matrix_error, float(np.max(abs(bn-target))))
        rows.append(dict(n=n, beta=str(sh/ch), prefix_speed_bound=str(ch/sh),
                         null_future_time=str(ch-sh)))
    thresholds = []
    for w in [F(6, 5), F(101, 100), F(1000001, 1000000)]:
        n = 1
        while coefficients(n)[0] - coefficients(n)[1]*w >= 0:
            n += 1
        prior = F(1) if n == 1 else coefficients(n-1)[0]-coefficients(n-1)[1]*w
        current = coefficients(n)[0]-coefficients(n)[1]*w
        assert prior >= 0 and current < 0
        thresholds.append(dict(speed=str(w), first_negative_n=n,
                               previous_time=str(prior), witness_time=str(current)))
    read_rows = [dict(repetitions=k, completion_time=str(k+F(3)),
                     distance=str(k*F(4, 5)),
                     observed_speed=str(k*F(4, 5)/(k+F(3))))
                 for k in [1, 2, 10, 100]]
    assert all(F(row['observed_speed']) < F(4, 5) for row in read_rows)
    add('finite_powers_and_null_boundary', maximum_matrix_error < 1e-12,
        prefix_table=rows, first_negative_witnesses=thresholds,
        maximum_matrix_error=maximum_matrix_error,
        future_null_kept_for_every_finite_n_by_exact_formula=True,
        finite_readout_overhead=3, illustrative_serial_readout=read_rows,
        fixed_readout_overhead_permission_certified=False,
        actual_finite_iteration_permission_certified=False)

    w = F(6, 5)
    finite_time = 1-F(3, 5)*w
    diagonal = np.array([1., 1., 1.])
    six_axes = [s*e for e in np.eye(3) for s in [-1., 1.]]
    axis_min = min(1-.6*float(e@diagonal) for e in six_axes)
    aligned = 1-.6*np.linalg.norm(diagonal)
    add('finite_reference_and_direction_coverage',
        finite_time == F(7, 25) and axis_min > 0 and aligned < 0,
        speed_outside_unit_cone=str(w), all_direction_beta_cap='3/5',
        all_direction_time_bracket=str(finite_time), remaining_speed_bound='5/3',
        six_axis_minimum=axis_min, same_event_aligned_axis_bracket=float(aligned),
        finite_axes_do_not_certify_all_directions=True)

    n = 2
    ch, sh = coefficients(n)
    scale = F(3, 2)
    t, x = F(1), F(6, 5)
    ideal = scale*(ch*t-sh*x)
    input_contrast = F(1, 5)
    eta_total = 2*F(1, 100)
    remaining = input_contrast-2*eta_total
    at = ax = ao = F(1, 1000)
    tau = scale*(ch*at+sh*ax)+ao
    upper = ideal+tau
    corners = [scale*(ch*(t+st*at)-sh*(x+sx*ax))+so*ao
               for st, sx, so in product([-1, 1], repeat=3)]
    p0, p1 = F(2, 5), F(3, 5)
    q0, q1 = p0+eta_total, p1-eta_total
    # Negative controls test whether the certificate (not the phenomenon) fails.
    large_probability_error = input_contrast/2
    large_clock_error = F(1, 4)
    add('finite_record_and_timestamp_certificate',
        ideal == -F(3, 16) and remaining == F(4, 25)
        and tau == F(7, 1000) and upper == -F(361, 2000)
        and max(corners) == upper and q1-q0 == remaining
        and input_contrast-2*large_probability_error == 0
        and ideal+large_clock_error > 0,
        chain_length=n, scale=str(scale), ideal_time=str(ideal),
        per_input_transport_error=str(eta_total), input_contrast=str(input_contrast),
        output_contrast_lower_bound=str(remaining), timestamp_error_bound=str(tau),
        exact_box_corners=len(corners), maximum_actual_time=str(upper),
        declared_output_probabilities=[str(q0), str(q1)],
        large_probability_error_certifies_signal=False,
        large_clock_error_certifies_negative_time=False,
        empirical_calibration=False)

    events = [(F(0), F(0)), (F(1), F(6, 5)), (F(2), F(0))]
    times = [ch*te-sh*xe for te, xe in events]
    edges = [(0, 1), (1, 2)]
    add('dependence_identity_is_not_every_clock_order',
        times == [0, -F(1, 8), F(17, 4)]
        and all(i < j for i, j in edges) and ch > 0,
        directed_chain=edges, transformed_times=[str(q) for q in times],
        physical_edge_identity_reversed=False, directed_cycle=False,
        stationary_local_clock_derivative=str(ch),
        all_remote_timestamp_orders_forward=False,
        full_physical_countermodel=False)

    rotation = np.diag([1., -1., -1., 1.])
    mixed = rotation@b
    residual = float(np.max(abs(mixed@mixed-np.eye(4))))
    eta = np.diag([1., -1., -1., -1.])
    metric_error = float(np.max(abs(mixed.T@eta@mixed-eta)))
    signs = [F(1, 10)*(ch-sh*F(6, 5)),
             F(7)*(ch-sh*F(6, 5))]
    add('mixed_involution_and_positive_scale',
        residual == 0 and metric_error == 0 and mixed[0, 1] != 0
        and all(s < 0 for s in signs),
        mixed_matrix=mixed.tolist(), involution_residual=residual,
        Lorentz_metric_residual=metric_error,
        time_space_mixing_without_unbounded_powers=True,
        positive_scale_time_examples=[str(s) for s in signs],
        physical_scale_removed=False)

    p, samples = F(1, 10), 3
    e = sum(F(comb(samples, k))*p**k*(1-p)**(samples-k)
            for k in range((samples+1)//2, samples+1))
    total_error = 1-(1-e)**2
    weights = [(1-e)**2, e*(1-e), e*(1-e), e*e]
    identity = np.eye(2, dtype=complex)
    xx = np.array([[0, 1], [1, 0]], complex)
    zz = np.diag([1., -1.]).astype(complex)
    yy = np.array([[0, -1j], [1j, 0]])
    bell = np.array([1., 0., 0., 1.])/math.sqrt(2)
    target = np.outer(bell, bell.conj())
    noisy = np.zeros((4, 4), complex)
    for weight, pauli in zip(weights, [identity, xx, zz, yy]):
        v = np.kron(pauli, identity)@bell
        noisy += float(weight)*np.outer(v, v.conj())
    choi_distance = float(np.sum(abs(np.linalg.eigvalsh(noisy-target)))/2)
    add('fixed_noisy_strategy_is_not_exact_quantum_delivery',
        e == F(7, 250) and total_error == F(3451, 62500)
        and abs(choi_distance-float(total_error)) < 1e-12,
        binary_noise=str(p), repetitions=samples, majority_error=str(e),
        half_diamond_error_from_Pauli_identity=str(total_error),
        normalized_choi_half_trace_distance=choi_distance,
        theorem_about_all_possible_error_corrected_strategies=False,
        physical_completion_latency_assumed_zero=False)

    return dict(schema='round1098_actual_observer_influence_bridge_v1',
        round=1098, status='PASS', numpy_version=np.__version__, groups=groups,
        scope=dict(conditional_full_influence_bridge=True,
            new_OC_or_MC_adopted=False, actual_observer_resources_certified=False,
            full_influence_cone_unconditionally_derived=False,
            physical_clock_scale_fixed=False, all_matter_dynamics_certified=False,
            entire_conjecture_decided=False, new_adopted_axioms=0,
            scientific_count_increment=0, scientific_count_total=3860))


def compare(a, b):
    if isinstance(a, dict):
        assert set(a) == set(b)
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-11), (a, b)
    else:
        assert a == b, (a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    out = calculate()
    path = HERE/'results.json'
    if args.write:
        with path.open('x', encoding='utf8') as f:
            json.dump(out, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        compare(out, json.loads(path.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1098, status='PASS', groups=len(out['groups']),
                         saved_result_matches=not args.write), ensure_ascii=False))


if __name__ == '__main__':
    main()
