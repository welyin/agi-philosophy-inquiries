"""Read-only count0 checks of the already adopted 1055 force and moment map."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys


def poly_integral(coefficients, left=F(0), right=F(1)):
    return sum(c*(right**(i+1)-left**(i+1))/F(i+1)
               for i, c in enumerate(coefficients))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, q):
    return tuple(q*x for x in a)


def norm_squared(a):
    return sum(x*x for x in a)


def calculate():
    mu = F(10**20, 2)
    sigma = F(1, 10**6)
    force = (6*mu, -12*mu)  # coefficient times each component of d
    momentum_increment = poly_integral(force)
    weighted_force = (force[0], force[1]-force[0], -force[1])
    position_increment_coefficient = poly_integral(weighted_force)/mu
    impulse_coefficient = (poly_integral(force, F(0), F(1, 2))
                           - poly_integral(force, F(1, 2), F(1)))
    assert momentum_increment == 0
    assert position_increment_coefficient == 1
    assert impulse_coefficient == 3*mu

    zero = (F(0),)*3
    d = (F(3, 10**6), F(4, 10**6), F(0))
    e = (F(-1, 10**6), F(0), F(1, 10**6))
    radius = F(1, 10**5)
    assert norm_squared(d) < radius**2
    assert norm_squared(scale(d, F(1, 2))) < radius**2
    assert add(d, scale(d, -1)) == zero
    assert add(scale(d, F(1, 2)), scale(d, F(1, 2))) == d
    assert scale(add(d, d), F(1, 2)) == d
    assert norm_squared(scale(d, F(1, 2))) == norm_squared(d)/4
    assert add(add(d, e), scale(d, -1)) == e

    norm_d = F(5, 10**6)
    assert norm_d**2 == norm_squared(d)
    canonical_cost = impulse_coefficient*norm_d
    half_cost = impulse_coefficient*norm_d/2
    assert half_cost == canonical_cost/2
    assert 2*half_cost == canonical_cost
    loop_cost = 2*canonical_cost
    assert loop_cost > 0

    # Zero endpoint does not identify whole histories or their resource totals.
    phase = 3*mu*norm_squared(d)/5
    two_halves_phase = 2*3*mu*norm_squared(scale(d, F(1, 2)))/5
    assert two_halves_phase == phase/2
    p_variance = 1/(4*sigma**2)
    r_variance_1 = sigma**2 + p_variance/mu**2
    r_variance_2 = sigma**2 + 4*p_variance/mu**2
    assert r_variance_2-r_variance_1 == 3*p_variance/mu**2 > 0

    # Convex mixtures preserve the unconditional moment contract, not each
    # reference-conditioned component if momenta cancel only in the mixture.
    r_a, r_b = F(1, 10**6), F(-2, 10**6)
    weight = F(2, 3)
    mixed = weight*r_a+(1-weight)*r_b
    assert mixed == 0
    displacement = F(1, 10**6)
    assert weight*(r_a+displacement)+(1-weight)*(r_b+displacement) == mixed+displacement
    velocity = F(1, 1000)
    reference_conditioned_means = (velocity+displacement, -velocity+displacement)
    assert sum(reference_conditioned_means)/2 == displacement
    assert reference_conditioned_means[0] != reference_conditioned_means[1]

    return {
        "scope": "count0 reuse of 1055; unconditional mean endpoint quotient, not full-history state quotient",
        "new_scientific_groups": 0,
        "new_empirical_groups": 0,
        "new_cognitive_principles": 0,
        "mu": str(mu),
        "sigma": str(sigma),
        "pulse_duration": "1",
        "net_momentum_increment_per_d": str(momentum_increment),
        "mean_position_increment_per_d": str(position_increment_coefficient),
        "absolute_impulse_per_norm_d": str(impulse_coefficient),
        "sample_d": [str(x) for x in d],
        "sample_norm_d": str(norm_d),
        "canonical_cost": str(canonical_cost),
        "half_cost": str(half_cost),
        "zero_endpoint_loop_actual_cost": str(loop_cost),
        "canonical_phase": str(phase),
        "two_half_pulses_phase": str(two_halves_phase),
        "one_wait_r_variance": str(r_variance_1),
        "two_waits_r_variance": str(r_variance_2),
        "variance_difference": str(r_variance_2-r_variance_1),
        "conditional_reference_means": [str(x) for x in reference_conditioned_means],
        "unconditional_reference_mean": str(displacement),
        "all_assertions_passed": True,
        "direction_effect_descends_from_all_histories": False,
        "new_optical_or_dynamics_simulation": False,
        "complete_M3A_or_dimension_origin_claimed": False
    }


if __name__ == "__main__":
    result = calculate()
    output = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("Use no arguments, or --save-exclusive for first creation.")
    else:
        assert result == json.loads(output.read_text(encoding="utf-8"))
    print("PASS: force moments, endpoint compositions, impulse budget and retained-history boundaries; count0.")
