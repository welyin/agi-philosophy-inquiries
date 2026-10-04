"""Working entry after 761: exact scale accounting and initial packet checks.

These are entry diagnostics, not a new completed scientific round.
The quadrature is one gauge-invariant local marginal, not full graph dynamics.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_loop_scale_transport as prev
import joint_background_source_closure as alg
TARGET = HERE / 'resource_scale_entry_results.json'


def exact_resource_identity():
    I = alg.I
    mass = np.array([[Q(1), Q(2)], [Q(2), Q(-1)]], dtype=object)
    hopping = np.array([[Q(0), Q(1)], [Q(1), Q(0)]], dtype=object)
    B = {(0, 1): mass, (0, 0): hopping}
    laplace = {(2, 0): I, (2, 2): I / Q(6), (1, 1): I * Q(5, 6)}
    potential = {(0, 4): I * Q(3, 8), (0, 2): I * Q(1, 9)}
    rows = []
    for M in (Q(1), Q(5), Q(11)):
        epsilon = 1 / M
        physical = alg.add(alg.scaled(laplace, -1 / (2 * M)),
                           alg.scaled(potential, M))
        physical = alg.add(physical, B)
        scaled = alg.add(alg.scaled(laplace, -epsilon**2 / 2), potential)
        scaled = alg.add(scaled, alg.scaled(B, epsilon))
        assert not alg.add(alg.scaled(physical, epsilon), scaled, Q(-1))
        rows.append(dict(M_exact=str(M), epsilon_exact=str(epsilon),
                         same_physical_time_generator=True))
    return dict(rows=rows, exact_rational_operator_identity=True,
                physical_hbar_fixed_to_one=True,
                bosonic_action_prefactor_is_added_model_input=True,
                not_derived_by_unit_change_or_by_adding_nodes=True)


def smooth_cutoff(scaled):
    """C-infinity cutoff: 1 inside 0.6, 0 outside 1."""
    z = np.abs(scaled)
    ans = np.ones_like(z)
    ans[z >= 1] = 0
    middle = (z > .6) & (z < 1)
    u = (z[middle] - .6) / .4
    a = np.exp(-1 / (1-u))
    b = np.exp(-1 / u)
    ans[middle] = a / (a+b)
    return ans


def original_marginal_check():
    phi, x, F, B, M, vp, vm, psi = prev.original_data()
    s0 = float(phi[4])
    A = 2 - float(phi[:4] @ phi[:4]) / 6
    radius = min(.5, .4 * (np.sqrt(6*A) - abs(s0)))
    assert radius > .2
    nodes, weights = np.polynomial.legendre.leggauss(512)
    u = radius * nodes
    weights = radius * weights
    cutoff = smooth_cutoff(u/radius)
    s = s0 + u
    Fs = A - s*s/6
    assert np.min(Fs) > 0
    ys = abs(prev.old.matter.Y['s'])
    record = .5 + np.sin(s) / 4
    energy = ys * s / np.sqrt(Fs)
    c1, c2 = .5, 1.5
    predicted_record = -(c2-c1)*np.sin(s0)/8
    predicted_energy = (c2-c1) * ys * A*s0/(4*F**2.5)
    rows = []
    for eps in (1/1024, 1/2048, 1/4096, 1/8192):
        measurements = []
        for c in (c1, c2):
            density = weights * cutoff**2 * np.exp(-u*u/(2*eps*c))
            density /= np.sum(density)
            measurements.append(dict(
                variance_coefficient=c,
                measured_variance_over_epsilon=float(density @ (u*u))/eps,
                record_probability=float(density @ record),
                normalized_mass_energy=float(density @ energy),
                canonical_momentum_variance_coefficient=1/(4*c),
                leading_uncertainty_product=c/(4*c)))
        first, second = measurements
        dr = (second['record_probability']-first['record_probability'])/eps
        de = (second['normalized_mass_energy']-first['normalized_mass_energy'])/eps
        gaussian_finite = np.sin(s0)/4 * np.exp(-eps*c1/2) * \
            np.expm1(-eps*(c2-c1)/2)/eps
        remainder_bound = abs(np.sin(s0))/4 * eps*(c1*c1+c2*c2)/8
        assert abs(dr-gaussian_finite) < 1e-10
        assert abs(dr-predicted_record) <= remainder_bound + 1e-10
        assert abs(de-predicted_energy) < 5e-5
        for row in measurements:
            assert abs(row['measured_variance_over_epsilon'] -
                       row['variance_coefficient']) < 1e-7
            assert row['leading_uncertainty_product'] == .25
        rows.append(dict(epsilon=eps, preparations=measurements,
                         record_difference_over_epsilon=dr,
                         physical_mass_difference_over_epsilon_squared=de))
    assert abs(predicted_record) > .01
    assert abs(rows[-1]['record_difference_over_epsilon']-predicted_record) < \
           abs(rows[0]['record_difference_over_epsilon']-predicted_record)
    return dict(original_background_s=s0, original_positive_F=F,
                original_Ys_modulus=float(ys), compact_s_radius=float(radius),
                limiting_record_difference_over_epsilon=float(predicted_record),
                limiting_physical_mass_difference_over_epsilon_squared=float(predicted_energy),
                rows=rows,
                scope='Only an initial gauge-invariant s marginal in the original positive-F patch. '
                      'A smooth compact cutoff and half-density preparation permit this marginal '
                      'inside a full Gauss packet. No full-graph evolution, bosonic loop, or gravity '
                      'solution is numerically calculated.')


def run():
    dependencies = ['research_note_353.md', 'research_note_357.md',
                    'research_note_525.md', 'research_note_574.md',
                    'research_note_741.md', 'research_note_758.md',
                    'research_note_761.md', 'joint_loop_scale_transport.py']
    return dict(kind='round762_working_entry', latest_completed_round=761,
                new_numbered_scientific_tests=0, diagnostic_groups=2,
                resource_identity=exact_resource_identity(),
                original_initial_marginal=original_marginal_check(),
                dependency_hashes={n: hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest()
                                   for n in dependencies},
                all_diagnostics_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('kind', 'diagnostic_groups',
                                             'all_diagnostics_passed')}))
