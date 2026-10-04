"""Round 530: a finite prepared reference patch in the unchanged ordinary torus H.

The lattice, three reference copies, thermal source and addressed Weyl preparation
are supplied. FFTs evolve exact quantum common means, not the full relative state.
The non-Gaussian relative readout bound is inherited from round 529. No geometry
generation, autonomous preparation, all-time readout or dimension selection claim.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from thermal_reference_dimension_model import thermal_q

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'finite_patch_reference_source_results.json'


def lap(values):
    out = 2 * values.ndim * values
    for axis in range(values.ndim):
        out = out - np.roll(values, 1, axis) - np.roll(values, -1, axis)
    return out


def source(n, d, r, s, axis):
    radius = r + s
    assert n % 2 == 1 and n >= 2 * radius + 3 and s >= 1
    coordinates = np.meshgrid(*([np.arange(n) - n // 2] * d), indexing='ij')
    norm = np.maximum.reduce([abs(v) for v in coordinates])
    return coordinates[axis] * (norm <= radius), norm <= r, coordinates


def eigenvalues(n, d):
    one = 4 * np.sin(np.pi * np.arange(n) / n)**2
    return np.sum(np.meshgrid(*([one] * d), indexing='ij'), axis=0)


def propagate(profile, t):
    omega = np.sqrt(eigenvalues(profile.shape[0], profile.ndim))
    transformed = np.fft.fftn(profile)
    mean = np.fft.ifftn(np.cos(t * omega) * transformed).real
    momentum = np.fft.ifftn(-omega * np.sin(t * omega) * transformed).real
    return mean, momentum


def energy(profile, momentum=None):
    kinetic = 0 if momentum is None else float(np.sum(momentum**2))
    return float((kinetic + np.sum(profile * lap(profile))) / 2)


def buffer_bound(d, r, s, time=F(1)):
    x = 4 * d * time**2
    q = x / ((2 * s + 3) * (2 * s + 4))
    assert q < 1
    return (r + s) * x**(s + 1) / math.factorial(2 * s + 2) / (1 - q)


def rational(value):
    value = F(value)
    return dict(numerator=value.numerator, denominator=value.denominator,
                decimal=float(value))


def source_and_energy_checks():
    rows = []
    for d, r, s in ((1, 1, 2), (2, 1, 3), (3, 1, 5)):
        radius = r + s
        for n in (2 * radius + 3, 2 * radius + 9):
            for axis in range(d):
                profile, core, _ = source(n, d, r, s, axis)
                assert int(np.sum(profile)) == 0
                work = profile.copy()
                zero_checks = []
                for k in range(1, s + 1):
                    work = lap(work)
                    exact_zero = int(np.max(abs(work[core])))
                    assert exact_zero == 0
                    zero_checks.append(exact_zero)
                next_value = int(np.max(abs(lap(work)[core])))
                assert next_value > 0
                expected = F((2 * radius + 1)**(d - 1) * radius *
                             (radius + 1) * (d + 2), 3)
                actual = int(np.sum(profile * lap(profile)))
                assert F(actual, 2) == expected
                # Both boundary faces and all internal edges are retained.
                edge_sum = sum(int(np.sum((profile - np.roll(profile, 1, a))**2))
                               for a in range(d))
                assert edge_sum == actual
                rows.append(dict(d=d, n=n, r=r, s=s, axis=axis,
                    total_mean=0, core_vanishing_powers=zero_checks,
                    first_unprotected_power_maximum=next_value,
                    preparation_energy_for_unit_gradient=rational(expected)))
    return rows


def evolution_checks():
    rows = []
    bound = float(buffer_bound(3, 1, 5))
    for n in (15, 21, 31):
        profile, core, _ = source(n, 3, 1, 5, 0)
        initial_energy = energy(profile)
        max_bias, max_energy_residual, max_zero_residual = 0., 0., 0.
        for t in np.linspace(0, 1, 33):
            mean, momentum = propagate(profile, float(t))
            max_bias = max(max_bias, float(np.max(abs(mean[core] - profile[core]))))
            max_energy_residual = max(max_energy_residual,
                abs(energy(mean, momentum) / initial_energy - 1))
            max_zero_residual = max(max_zero_residual, abs(float(np.mean(mean))))
        assert max_bias < bound
        assert max_energy_residual < 5e-14 and max_zero_residual < 1e-14
        rows.append(dict(n=n, time_points=33, maximum_sampled_core_bias=max_bias,
            analytic_all_times_bias_bound=bound,
            relative_energy_residual=max_energy_residual,
            zero_mode_mean_residual=max_zero_residual))
    return rows


def covariance_checks():
    rows = []
    for n in (15, 21, 31):
        ell = eigenvalues(n, 3)
        spectral = np.zeros_like(ell)
        active = ell > 0
        spectral[active] = thermal_q(ell[active], 1.)
        covariance = np.fft.ifftn(spectral).real
        site = float(covariance[(0, 0, 0)])
        assert site < float(F(143, 96))
        pairs = [2 * (site - float(covariance[tuple(k % n for k in v)]))
                 for v in itertools.product(range(-2, 3), repeat=3) if any(v)]
        assert min(pairs) >= 0 and max(pairs) < 12 * (1 + math.sqrt(3))
        # A Weyl shift changes means only. Nonzero common modes stay thermal;
        # the separately prepared common zero mode cancels from any difference.
        rows.append(dict(n=n, nonzero_common_site_variance=site,
            shell_site_upper=rational(F(143, 96)),
            maximum_core_common_difference_variance=float(max(pairs)),
            inherited_local_difference_bound=12 * (1 + math.sqrt(3)),
            zero_mode_site_variance_at_time_one=1 / n**3))
    return rows


def certificate():
    d, r, s, g, nu = 3, 1, 5, 512, F(1, 4)
    radius = r + s
    bound = buffer_bound(d, r, s)
    assert bound == F(936, 23375) and bound < F(1, 16)
    sites, reads = (2 * r + 1)**d, d * (2 * r + 1)**d
    pairs = sites * (sites - 1) // 2
    spacing, tolerance = F(1, 64), F(1, 4)
    noise_allowance = tolerance - F(1, 8) - spacing
    assert noise_allowance == F(7, 64)
    pair_proxy, site_proxy = F(347, 4), F(223, 8)
    exponent = g**2 * noise_allowance**2 / (2 * pair_proxy)
    assert exponent == F(6272, 347) and exponent > 18
    pair_failure = F(2 * d * pairs, 2**18)
    overflow_exponent = g**2 * F(15, 16)**2 / (2 * site_proxy)
    assert overflow_exponent == F(921600, 223)
    # e > 2, 2*reads < 2**8, exponent > 28 imply overflow < 2**(-20).
    assert 2 * reads < 2**8 and overflow_exponent > 28
    total_failure = pair_failure + F(1, 2**20)
    assert total_failure == F(8425, 1048576) and total_failure < F(1, 100)
    symbols = int(4 / spacing) + 1 + 1
    bits = (symbols - 1).bit_length()
    assert (symbols, bits, reads * bits) == (258, 9, 729)
    prep_energy = F(d * g**2 * (2 * radius + 1)**(d - 1) *
                    radius * (radius + 1) * (d + 2), 3)
    meter_energy = reads / (8 * nu**2)
    assert prep_energy == 9303490560 and meter_energy == 162
    return dict(d=d, r=r, s=s, R=radius, maximum_time=1, minimum_odd_n=15,
        gradient=g, beta=1, nu=rational(nu),
        core_sites=sites, reference_copies=d, actual_reads=reads,
        pairs=pairs, pair_component_constraints=d*pairs,
        all_times_common_mean_bias= rational(bound),
        pair_subgaussian_proxy_upper=rational(pair_proxy),
        site_subgaussian_proxy_upper=rational(site_proxy),
        random_difference_allowance=rational(noise_allowance),
        pair_failure_exponent=rational(exponent),
        overflow_failure_exponent=rational(overflow_exponent),
        conservative_total_failure_upper=rational(total_failure),
        success_lower=rational(1-total_failure),
        numeric_symbols=symbols-1, symbols_including_overflow=symbols,
        bits_per_record=bits, raw_payload_bits=reads*bits,
        total_mean_preparation_energy=rational(prep_energy),
        batch_meter_energy_increment=rational(meter_energy))


def quantize(values):
    values = np.asarray(values)
    # Code 257 is explicit overflow. Valid codes 0..256 decode to [-2,2].
    bad = (values < -2) | (values > 2)
    codes = np.rint(np.clip(values, -2, 2) * 64).astype(int) + 128
    codes[bad] = 257
    return codes


def decode(codes):
    assert np.all((codes >= 0) & (codes <= 256)), 'overflow: no decoded coordinate'
    return (codes - 128) / 64.


def record_checks():
    rows = []
    for n, t in ((15, 0.), (15, 1.), (21, .5), (31, 1.)):
        means, targets = [], []
        for axis in range(3):
            profile, core, coordinates = source(n, 3, 1, 5, axis)
            mean, _ = propagate(profile, t)
            means.append(mean[core])
            targets.append(coordinates[axis][core])
        means, targets = np.array(means), np.array(targets)
        # Deliberately correlated bounded records, NOT Monte Carlo draws from
        # the quantum state. They test decoding on the theorem's good event.
        signs = np.where(np.arange(means.size).reshape(means.shape) % 2, 1., -1.)
        noise = F(7, 128) * signs.astype(float)
        analog = means + np.asarray(noise, dtype=float) + .5
        codes = quantize(analog)
        decoded = decode(codes)
        assert np.max(abs(decoded - analog)) <= 1/128 + 1e-15
        error, minimum_separation = 0., float('inf')
        for i in range(27):
            for j in range(i):
                actual_difference = decoded[:, i] - decoded[:, j]
                error = max(error, float(np.max(abs(actual_difference -
                                                    targets[:, i] + targets[:, j]))))
                minimum_separation = min(minimum_separation,
                                         float(np.max(abs(actual_difference))))
        assert error < .25 and minimum_separation >= .75
        rows.append(dict(n=n, t=t, decoded_pair_maximum_error=error,
            minimum_record_linf_separation=minimum_separation,
            raw_payload_bits=int(codes.size * 9),
            bounded_test_records_not_quantum_samples=True))
    probe = np.array([-2.001, -2., 0., 2., 2.001])
    assert quantize(probe).tolist() == [257, 0, 128, 256, 257]
    rejected = False
    try:
        decode(quantize(probe))
    except AssertionError:
        rejected = True
    assert rejected
    return dict(cases=rows, overflow_explicitly_rejected=True,
                quantum_confidence_from_analytic_bound_not_sampling=True)


def scope_and_resource_checks():
    # The source-energy identity is for extra coherent means, not background.
    rows = []
    for n in (15, 21, 31):
        ell = eigenvalues(n, 3)
        nonzero = ell > 0
        background = float(np.sum(ell[nonzero] * thermal_q(ell[nonzero], 1.)) + .25)
        rows.append(dict(n=n, one_copy_common_background_energy=background,
                         three_copy_extra_preparation_energy=9303490560))
    assert all(a['one_copy_common_background_energy'] < b['one_copy_common_background_energy']
               for a, b in zip(rows, rows[1:]))
    return rows


def run():
    source_rows = source_and_energy_checks()
    evolution = evolution_checks()
    covariance = covariance_checks()
    exact_certificate = certificate()
    records = record_checks()
    resources = scope_and_resource_checks()
    dependencies = ('research_note_529.md', 'nonlinear_thermal_reference.py',
        'nonlinear_thermal_reference_results.json', 'research_note_522.md',
        'thermal_reference_dimension_model.py', 'research_note_287.md', 'research_note_290.md')
    return dict(date='2026-09-30', round=530, scientific_base_through_round=529,
        tests_run=6, failures=0, errors=0, finite_source_and_exact_energy=source_rows,
        exact_common_mean_diagnostics=evolution, common_covariance_diagnostics=covariance,
        finite_record_certificate=exact_certificate, finite_decoder_diagnostics=records,
        resource_scope_diagnostics=resources,
        dependency_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                           for name in dependencies},
        scope=dict(same_relative_gate_and_ordinary_periodic_H=True,
            permanent_twist_or_boundary_force_required=False,
            finite_addressed_Weyl_preparation_assumed=True,
            graph_geometry_and_coordinate_axes_supplied=True,
            thermal_source_and_readout_controls_supplied=True,
            full_relative_quantum_state_not_replaced_by_Gaussian=True,
            old_infrared_dimension_theorem_not_rerun=True,
            finite_patch_single_selected_time_guarantee=True,
            strict_lightcone_claimed=False, source_autonomously_generated=False,
            original_degree_three_graph_to_lattice_map_proved=False,
            unique_dimension_derived=False, global_atlas_proved=False,
            continuous_geometry_or_GR_derived=False, stage_complete=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text('utf8')) == json.loads(json.dumps(result))
    else:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    cert = result['finite_record_certificate']
    print(json.dumps(dict(round=530, tests=6, failures=0,
        maximum_sampled_bias=max(row['maximum_sampled_core_bias'] for row in
                                 result['exact_common_mean_diagnostics']),
        success_lower=cert['success_lower']['decimal'],
        raw_payload_bits=cert['raw_payload_bits']), ensure_ascii=False))
