"""Ordinary charged-current pair source and a three-null, untagged spectroscopy test."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TAGS = ("d", "s", "b")
P = 100_000.0                  # GeV, common magnitude of momentum
T = 0.04                       # GeV^-1, uniformly randomized single delay
TICKS = 256                    # equal weights at the clock-bin midpoints
LIGHT_MAX = 5.0                # GeV, explicit low-spectrum envelope at the colour threshold
PHASE_BUDGET = 1.0e-3          # absolute relative phase error at each dark port, radians
DELTA = 0.05


def energy(mass):
    return mass**2 / (math.hypot(P, mass) + P)


def transmission(frequency, roots, duration=T):
    """Uniform time average of the product of three dark-port probabilities."""
    terms = []
    for powers in itertools.product((-1, 0, 1), repeat=3):
        coefficient = math.prod(0.5 if k == 0 else -0.25 for k in powers)
        angular = math.fsum([frequency * sum(powers), *(-k * r for k, r in zip(powers, roots))])
        x = angular * duration
        terms.append(coefficient * (math.sin(x) / x if x else 1.0))
    return math.fsum(terms)


def polarization_shape(top_ratio, w_ratio):
    return math.sqrt((1-top_ratio-w_ratio)**2 - 4*top_ratio*w_ratio) * (
        (1-top_ratio)**2 + w_ratio*(1+top_ratio) - 2*w_ratio**2
    )


def exposure(click):
    return math.ceil(math.log(DELTA) / math.log1p(-click))


def assess(record, heavy, vev, top=172.5, w_mass=80.379):
    weights = {}
    admixture = dict.fromkeys(TAGS, 0.0)
    max_shift = 0.0
    for tag, data in record["heavy_modes"].items():
        pole = heavy + data["heavy_excess_GeV"]
        max_shift = max(max_shift, data["heavy_excess_GeV"])
        eps = data["left_doublet_probability"]
        kinematic = polarization_shape((top/pole)**2, (w_mass/pole)**2)
        weights[tag] = data["partial_W_u_c_t_GeV"][2] / (pole**3/(16*math.pi*vev**2)*eps*kinematic)
        admixture[tag] = eps
    omitted, = set(TAGS) - weights.keys()
    weights[omitted] = 1 - math.fsum(weights.values())
    assert min(weights.values()) > 0

    p_one, p_two, missing_overlap = [], [], []
    for i, j in itertools.combinations(TAGS, 2):
        k, = set(TAGS) - {i, j}
        ordinary_pair_weight = weights[k]
        a, b = admixture[i], admixture[j]
        leakage = a + b - a*b
        p_one.append(ordinary_pair_weight*(a*(1-b) + (1-a)*b))
        p_two.append(ordinary_pair_weight*a*b)
        missing_overlap.append(ordinary_pair_weight*leakage/(1+math.sqrt(1-leakage)))
    one, two, lost = map(math.fsum, (p_one, p_two, missing_overlap))

    # All ordinary pair roots lie in [0, 2*energy(LIGHT_MAX)] after subtracting 2p.
    # For one/two heavy poles, bound the variation of each sine-squared factor.
    gap = energy(heavy)
    energy_shift = max_shift  # |dE/dm|<=1, including rounding of the saved common M below.
    rounding = 1e-6
    avg, errs = [], []
    for n_heavy in (1, 2):
        shift = n_heavy*(energy_shift + rounding) + (2-n_heavy)*energy(LIGHT_MAX)
        delays = T*(np.arange(TICKS)+0.5)/TICKS
        finite = float(np.mean(np.sin(n_heavy*gap*delays/2)**6))
        assert abs(finite - transmission(n_heavy*gap, (0.0,)*3)) < 3e-5
        avg.append(finite)
        errs.append(1.5*T*(shift+2*energy(LIGHT_MAX)))
    central = one*avg[0] + two*avg[1]
    bound = one*errs[0] + two*errs[1]
    lo, hi = central-bound, central+bound
    assert 0 < lo <= hi < one+two
    guaranteed = lo - 1.5*PHASE_BUDGET*(one+two)
    equal_prior_min = math.ceil(math.log(4*DELTA*(1-DELTA))/(2*math.log1p(-lost)))
    return dict(top_row_moduli_squared=weights, left_heavy_probabilities=admixture,
                prob_one_heavy=one, prob_two_heavy=two, dark_spectrum_probability=one+two,
                one_minus_pure_overlap=lost,
                copies_lower_bound_any_measurement_equal_prior_error_5pct=equal_prior_min,
                copies_lower_bound_any_zero_background_test_95pct_power=math.ceil(
                    math.log(DELTA)/(2*math.log1p(-lost))),
                transmission_one_two_at_zero_light_mass=avg, transmission_abs_envelopes=errs,
                click_probability_envelope=[lo, hi],
                trials_for_95pct_click_with_phase_budget=exposure(guaranteed),
                independent_trials_for_95pct_click=[exposure(hi),exposure(lo)])


def run():
    data = json.loads((HERE/"../1768/independent_results.json").read_bytes())
    # Quadrature, alias at a fixed delay, and stable analytic integration across close light roots.
    roots = (0.0001, 0.0002, 0.00022)
    grid, quad_weights = np.polynomial.legendre.leggauss(64)
    delays = (grid+1)*T/2
    for omega in (*roots, 2.1, 278.2, 600.0):
        by_quad = np.dot(quad_weights/2, np.prod(np.sin(
            (omega-np.array(roots)[:,None])*delays[None,:]/2)**2,axis=0))
        assert abs(transmission(omega,roots)-by_quad) < 4e-15
    alias = roots[0]+2*math.pi/T
    assert math.sin((alias-roots[0])*T/2)**2<1e-28 and transmission(alias,roots)>0.1
    by_point = {
        source: {key:assess(row,data["common_mass_GeV"],data["vev_GeV"])
                 for key,row in rows.items()} for source,rows in data["by_calibration"].items()
    }
    mixing = [0.1,0.01,0.001,0.0001]
    return dict(common_source="fixed weak up/charm antisymmetric pair; ordinary three-line controls",
                common_momentum_GeV=P, max_single_delay_GeV_inverse=T, clock_bins=TICKS,
                max_light_working_mass_GeV=LIGHT_MAX, ordinary_phase_error_example=PHASE_BUDGET,
                false_click_per_trial_bound=(math.sin(PHASE_BUDGET/2)**2 *
                   math.sin((PHASE_BUDGET+2*energy(LIGHT_MAX)*T)/2)**4),
                adaptive_weak_calls_needed_for_equal_prior_error_5pct=[
                    [theta,math.ceil(math.asin(1-2*DELTA)/theta)] for theta in mixing],
                by_calibration=by_point)


def rounded(value):
    if isinstance(value, dict):
        return {k:rounded(v) for k,v in value.items()}
    if isinstance(value, list):
        return [rounded(v) for v in value]
    return float(f"{value:.11g}") if isinstance(value,float) else value


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=rounded(run())
    payload=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == result
    else:
        path.write_text(payload,encoding="utf-8")
    print("1774 fixed-source weak reach verified",hashlib.sha256(payload.encode()).hexdigest()[:12])
