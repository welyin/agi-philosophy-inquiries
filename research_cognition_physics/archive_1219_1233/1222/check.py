"""Round 1222: finite-budget/finite-switching certificates, not a detector simulation."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]


def f(energy, beta):
    return 1.0 / beta if energy == 0 else energy / math.expm1(beta * energy)


def curvature(energy, beta):
    t = abs(beta * energy / 2)
    if t < 1e-3:
        z = beta * energy
        return beta * (1 / 6 - z*z / 60 + z**4 / 1008)
    return beta / 2 * (t * math.cosh(t) - math.sinh(t)) / math.sinh(t)**3


def certificate(energy, duration, acceleration=1.0):
    beta = 2 * math.pi / acceleration
    plus, minus = f(energy, beta), f(-energy, beta)
    moment = 4 * math.pi**2 / (3 * duration**2)
    k = beta * moment / 12
    bias = math.log1p(k / plus) - math.log1p(k / minus)
    return dict(gap=energy, duration=duration, second_moment=moment,
                delta_upper=k, f_plus=plus, f_minus=minus,
                log_ratio_bias_upper=bias, passes_switching=bias <= 0.1)


def threshold(energy, tolerance, acceleration=1.0):
    beta = 2 * math.pi / acceleration
    assert 0 < tolerance < beta * energy
    plus, minus = f(energy, beta), f(-energy, beta)
    k = plus * minus * math.expm1(tolerance) / (minus - math.exp(tolerance) * plus)
    best = math.pi * math.sqrt(beta / (12 * k))
    return dict(gap=energy, k_star=k, optimal_H1_duration=best,
                cos_squared_duration=2 * best / math.sqrt(3))


def budget_at(duration, mass, remainder_mass, acceleration=1.0):
    rapidity = acceleration * duration / 2
    energy = mass * math.cosh(rapidity)
    momentum = mass * math.sinh(rapidity)
    return energy + math.hypot(remainder_mass, momentum)


def maximum_duration(budget, mass, remainder_mass, acceleration=1.0):
    assert budget >= mass + remainder_mass
    return 2 / acceleration * math.acosh(
        (budget**2 + mass**2 - remainder_mass**2) / (2 * budget * mass))


def run_checks():
    beta = 2 * math.pi
    grid = np.linspace(-10, 10, 2001)
    identities = [abs(f(-float(x), beta) - f(float(x), beta) - float(x)) for x in grid]
    curves = np.array([curvature(float(x), beta) for x in grid])
    assert max(identities) < 4e-15
    assert curves.min() >= 0 and curves.max() <= beta / 6 + 1e-12
    # Directly average several even spectral measures, independently of the envelope formula.
    max_violation = 0.0
    for gap in (.1, .2, .7):
        for offset in (.01, .1, 1.0):
            avg = .2*f(gap, beta) + .4*(f(gap+offset, beta)+f(gap-offset, beta))
            delta = avg-f(gap, beta)
            k = beta*(.8*offset**2)/12
            actual_bias = beta*gap-math.log((avg+gap)/avg)
            envelope = math.log1p(k/f(gap,beta))-math.log1p(k/f(-gap,beta))
            assert delta >= -1e-14 and delta <= k + 1e-14
            assert actual_bias >= -1e-14 and actual_bias <= envelope + 1e-14
            max_violation = max(max_violation, actual_bias-envelope)

    # Actual quadrature of switching and derivative norms (not a response quadrature).
    nodes, weights = np.polynomial.legendre.leggauss(128)
    norm_checks = []
    for duration in (20., 30.):
        times = nodes * duration / 2
        chi = np.cos(math.pi*times/duration)**2
        derivative = -math.pi/duration*np.sin(2*math.pi*times/duration)
        norm = float(duration/2*np.dot(weights, chi**2))
        derivative_norm = float(duration/2*np.dot(weights, derivative**2))
        assert abs(norm-3*duration/8) < 1e-12
        assert abs(derivative_norm/norm-4*math.pi**2/(3*duration**2)) < 1e-14
        norm_checks.append(dict(duration=duration, norm=norm, derivative_norm=derivative_norm))
    thresholds = [threshold(gap, .1) for gap in (.1,.2)]
    certified_min = max(x['cos_squared_duration'] for x in thresholds)
    optimal_min = max(x['optimal_H1_duration'] for x in thresholds)
    for row in thresholds:
        boundary = certificate(row['gap'], row['cos_squared_duration'])
        assert abs(boundary['log_ratio_bias_upper']-.1) < 1e-14
    rows = [certificate(gap,duration) for duration in (20.,30.) for gap in (.1,.2)]

    mass, remainder = 100., 500.
    budgets = []
    for value in (1e5,1e8):
        tmax = maximum_duration(value,mass,remainder)
        assert abs(budget_at(tmax,mass,remainder)/value-1) < 1e-13
        budgets.append(dict(budget=value, maximum_symmetric_duration=tmax,
            maximum_rest_start_duration=tmax/2,
            cos_squared_certificate_interval=[certified_min,tmax] if tmax>=certified_min else None,
            any_H1_second_moment_certificate_possible=tmax>=optimal_min,
            fixed_windows=[dict(duration=t, energy_required=budget_at(t,mass,remainder),
                passes_budget=budget_at(t,mass,remainder)<=value,
                passes_both_gaps_and_budget=all(certificate(e,t)['passes_switching'] for e in (.1,.2))
                    and budget_at(t,mass,remainder)<=value) for t in (20.,30.)]))
    massless_checks = []
    for rapidity in (.5,1.,3.):
        total = budget_at(2*rapidity,mass,0.)
        alone = mass*math.cosh(rapidity)
        assert abs(total/(mass*math.exp(rapidity))-1)<1e-14 and total>alone
        massless_checks.append(dict(rapidity=rapidity, complete_budget=total,
                                   material_only=alone, missing_recoil_energy=total-alone))
    # Endpoint minimax check: any shifted interval needs at least the symmetric energy.
    for length in (1.,5.,20.):
        lower = budget_at(length,mass,remainder)
        for shift in np.linspace(-5,5,21):
            endpoint = max(abs(shift-length/2),abs(shift+length/2))
            needed = budget_at(2*endpoint,mass,remainder)
            assert needed >= lower-1e-9

    coupling = .02
    lower_probabilities = [dict(gap=row['gap'],duration=row['duration'],
        leading_probability_lower=coupling**2*(3*row['duration']/8)/(2*math.pi)*row['f_plus']) for row in rows]
    lower_probability = min(row['leading_probability_lower'] for row in lower_probabilities)
    relative_error = math.tanh(.1/2)
    absolute_error_contract = relative_error*lower_probability
    assert abs(math.log((1+relative_error)/(1-relative_error))-.1)<1e-14
    assert not any(w['passes_both_gaps_and_budget'] for b in budgets for w in b['fixed_windows'])
    assert budgets[0]['cos_squared_certificate_interval'] is None
    assert budgets[1]['cos_squared_certificate_interval'] is not None
    sources = [HERE/'proof.md',HERE/'check.py',
        BASE/'archive_1106_1119/1114/proof.md',
        BASE/'archive_990_1008/research_note_1004.md',
        HERE.parent/'_admission/finite_source_unruh_screen_after1221.md']
    return dict(schema='round1222_window_certificate_v1',status='PASS',check_groups=4,
        scope='Analytic certificate and conditional resource feasibility, not a full Unruh implementation',
        inputs=dict(acceleration=1.,mass=mass,remainder_mass_lower=remainder,gaps=[.1,.2],
            fixed_durations=[20.,30.],budgets=[1e5,1e8],coupling=coupling,
            switching_log_tolerance=.1,actual_probability_log_tolerance=.1),
        spectral_checks=dict(identity_max_error=max(identities),curvature_min=float(curves.min()),
            curvature_max=float(curves.max()),envelope_positive_violation=max_violation,
            grid_is_not_continuous_proof=True),switching_norms=norm_checks,
        thresholds=thresholds,finite_switching_certificates=rows,budgets=budgets,
        recoil_deletion_witness=massless_checks,probability_contract=dict(
            leading_probability_lowers=lower_probabilities,common_probability_lower=lower_probability,
            relative_error_max=relative_error,required_full_absolute_error=absolute_error_contract,
            full_error_achieved=False,leading_coefficient_is_not_full_probability=True),
        fixed_windows_jointly_pass=False,
        higher_budget_continuous_certificate_intersection_nonempty=True,
        actual_common_source_clock_field_instrument_constructed=False,stage_complete=False,
        source_sha256={p.relative_to(BASE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run_checks()
    output=HERE/'results.json'
    if args.verify:
        saved=json.loads(output.read_text(encoding='utf-8'))
        assert saved==result,'Saved output differs from recomputation'
    else:
        if output.exists():
            assert json.loads(output.read_text(encoding='utf-8'))==result,'Do not overwrite different results'
        else:
            output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PASS',check_groups=4,
        cos_squared_minimum=max(x['cos_squared_duration'] for x in result['thresholds']),
        high_budget_interval=result['budgets'][1]['cos_squared_certificate_interval'],
        fixed_windows_jointly_pass=False,
        required_absolute_error=result['probability_contract']['required_full_absolute_error']),ensure_ascii=False))


if __name__=='__main__':
    main()
