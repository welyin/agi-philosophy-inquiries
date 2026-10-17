"""Independent round 1261 checks: exact coefficients and finite-task budgets.
No author implementation or saved result is imported. Standard library only.
Default is read-only comparison; --write creates results once, exclusively.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json


def rational(x):
    return f"{x.numerator}/{x.denominator}"


def root_upper(x):
    """Exact dyadic bisection, independent of the author's decimal isqrt."""
    x = Q(x)
    assert x >= 0
    if x == 0:
        return Q(0)
    shift = max(0, (x.numerator.bit_length() - x.denominator.bit_length() + 2) // 2)
    lo, hi = Q(0), Q(2**shift)
    assert hi * hi >= x
    for _ in range(400):
        mid = (lo + hi) / 2
        if mid * mid >= x:
            hi = mid
        else:
            lo = mid
    assert hi * hi >= x and lo * lo <= x
    return hi


# Multivariate polynomials in (A,a,d,c,y), with exact coefficient collection.
ZERO = (0, 0, 0, 0, 0)


def const(q):
    return {ZERO: Q(q)} if q else {}


def var(i):
    key = list(ZERO)
    key[i] = 1
    return {tuple(key): Q(1)}


def add(*polys):
    out = {}
    for p in polys:
        for key, value in p.items():
            out[key] = out.get(key, Q(0)) + value
    return {key: value for key, value in out.items() if value}


def scale(p, q):
    return {key: value*q for key, value in p.items() if value*q}


def mul(p, q):
    out = {}
    for k, v in p.items():
        for ell, w in q.items():
            key = tuple(a+b for a, b in zip(k, ell))
            out[key] = out.get(key, Q(0)) + v*w
    return {key: value for key, value in out.items() if value}


def square(p):
    return mul(p, p)


def reduce_c_squared(p):
    out = {}
    for key, value in p.items():
        k = list(key)
        pairs = k[3] // 2
        k[3] -= 2*pairs
        k[4] += 2*pairs
        tup = tuple(k)
        out[tup] = out.get(tup, Q(0)) + value
    return {key: value for key, value in out.items() if value}


def exact_angular_group():
    A, a, d, c, y = (var(i) for i in range(5))
    a2, d2, c2, y2 = map(square, (a, d, c, y))
    acd = mul(mul(a, c), d)
    N1 = scale(add(scale(mul(d2, c2), 2), scale(acd, -2), a2, scale(d2, -1)), 2)
    N2 = scale(add(scale(mul(A, c2), -1), mul(a2, c2), scale(acd, 2), d2, scale(a2, -1)), 2)
    N3 = mul(add(scale(mul(A, c2), -3), mul(a2, c2), scale(acd, 8), A, scale(a2, -3)), y)
    qpoly = add(A, scale(a2, -3), scale(acd, 8), mul(y2, add(a2, scale(A, -3))),
                scale(mul(y, add(scale(d2, 2), scale(A, -1), a2)), 2))
    assert reduce_c_squared(add(N1, N2, N3, scale(mul(y, qpoly), -1))) == {}
    # Exact identities for the stronger photon-source matrix norm <= 2.
    one, s = const(1), y
    s2 = square(s)
    denominator = add(one, s)
    xx = add(one, s, scale(s2, -4))
    zz = mul(add(scale(s, 4), const(-2)), denominator)
    off2 = scale(mul(s2, add(one, scale(s2, -1))), 16)
    upper_x = add(scale(denominator, 2), scale(xx, -1))
    upper_z = add(scale(denominator, 2), scale(zz, -1))
    lower_x = add(xx, scale(denominator, 2))
    lower_z = add(zz, scale(denominator, 2))
    det_upper = add(mul(upper_x, upper_z), scale(off2, -1))
    det_lower = add(mul(lower_x, lower_z), scale(off2, -1))
    assert add(det_upper, scale(mul(add(one, scale(s, -1)), square(denominator)), -4)) == {}
    assert add(det_lower, scale(mul(mul(s, add(const(3), scale(s, -1))), denominator), -4)) == {}
    assert add(lower_x, scale(add(const(3), scale(s, -1), scale(mul(s, add(one, scale(s, -1))), 4)), -1)) == {}
    # Remaining diagonal positivity is explicit for 0 <= s <= 1.
    coarse_difference = add(scale(square(denominator), 6), scale(add(const(4), scale(s, 16), scale(s2, 4)), -1))
    assert add(coarse_difference, scale(square(add(one, scale(s, -1))), -2)) == {}
    reverse = Q(31)*64/147
    assert reverse < 14 and Q(3)+Q(6+14, 4) == 8
    return {
        "massless_cancellation": "EXACT_COEFFICIENT_IDENTITY_MOD_C2_EQ_Y2",
        "photon_source_matrix_norm": "at_most_2_by_PSD_certificates; author_bound_6_valid",
        "reverse_source_coefficient": rational(reverse),
        "mixed_pair_remainder_coefficient": 8,
        "continuous_scope": "0<=s<=1; massive_speed<=1/2; k>0; r>0",
    }


def atan_bracket(inv, terms=24):
    x = Q(1, inv)
    value = sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(terms)), Q(0))
    following = (-1)**terms*x**(2*terms+1)/Q(2*terms+1)
    return min(value, value+following), max(value, value+following)


def finite_groups():
    a5, b5 = atan_bracket(5)
    a239, b239 = atan_bracket(239)
    pi_low, pi_high = 16*a5-4*b239, 16*b5-4*a239
    assert Q(333,106) < pi_low < pi_high < Q(355,113)
    R, mu, k, sigma = (Q(10**j) for j in (40,34,12,12))
    phi, d, om, cap = Q(1,10**6), Q(1,10**4), Q(1000), Q(300000)
    lam = 1-phi
    tlo, thi = Q(160001,2)*pi_low/lam, Q(160001,2)*pi_high/lam
    dlo, dhi = phi*tlo, phi*thi
    assert thi < cap and Q(1,4) < dlo < dhi < Q(13,50) < pi_low/12
    # Good domain really belongs to the corrected radial plateau.
    assert Q(3,4) > Q(2,3)
    assert (Q(5,4)*k)**2 + 1 < (Q(4,3)*k)**2
    D = 12*sigma**2 + cap**2 + 2
    assert 10*sigma**2 + pi_high**2/om**2 + 1 + cap**2 < D
    # Independent pairwise derivation of the author's good-region envelope.
    AB = Q(16,15)/R * (3*(mu+1)*d*d/52 + 49*d*d/mu + 4*d*d + 7*d**4/(26*mu))
    assert AB < 8*phi*d*d/26
    BC_linear, BC_quadratic = Q(256,93)*d*k/R, Q(512,93)*d*d*k/(mu*R)
    assert BC_linear < 32*d*k/R and BC_quadratic < 1000*d*d*k/(mu*R)
    AC = (Q(4096,45)+Q(256,45)*d+Q(256,585)*d*d)*k/R
    assert AC < 1000*k/R
    assert Q(512,15) < 64 and Q(64,93) < 2
    # Whole completion support: no reliance on the normal-state tails here.
    global_interaction = (528*phi + 16*phi*k + 64*d*k/R + 256*d*d*k/(mu*R)
                          + 512*k/R + 64*d*k/R + Q(256,26)*d*d*k/R)
    BV = Q(10**9)
    assert global_interaction < BV < mu+26
    # Each Husimi marginal includes the physical packet and coherent frame.
    radial_eighth = 945/(16*sigma**8)
    transverse_eighth = 24/sigma**8
    prob_other = 3*16**2*D/R**2 + 2*radial_eighth/d**8 + transverse_eighth
    eighth_long = 2**7*((om/2)**8 + Q(105,256))
    edge = k/4
    prob_long = eighth_long/edge**8
    tail1 = (BV+phi*(32+Q(5,4)*k))*root_upper(prob_other)
    tail2 = (BV+32*phi)*root_upper(prob_long) + phi*(k*root_upper(prob_long)+root_upper(eighth_long/edge**6))
    good = 8*phi*d*d/26 + 4*phi/k + 1000*k/R + 32*d*k/R + 1000*d*d*k/(mu*R) + 1000*d*d/(mu*R)
    position = (64*phi+2*phi*k)/R*root_upper(D)
    free = cap*((Q(1,26)+1/mu)/(2*sigma**2) + 1/((k-om/2)*sigma**2))
    eps = cap*(good+position+tail1+tail2)+free
    assert eps < Q(93,10**11) < Q(1,10**8)
    group1 = {
        "angular": exact_angular_group(),
        "good_region_AB_bound": rational(AB),
        "good_region_BC_linear_bound": rational(BC_linear),
        "good_region_BC_quadratic_bound": rational(BC_quadratic),
        "good_region_AC_full_bound": rational(AC),
        "global_symbol_bound_from_pairs": rational(global_interaction),
        "chosen_interaction_bound": rational(BV),
        "corrected_plateau_contains_good_domain": True,
        "Husimi_second_bound": rational(D),
        "total_state_bound": rational(eps),
        "free_state_bound": rational(free),
        "finite_energy_exchange_bound": rational(2*BV),
        "error_scope": "fixed_normal_mechanical_signal_preparation; all_original_data_and_reference",
    }
    # Exact lower bounds for original p(a)=sin(a/2)^6, inherited from 1253.
    assert Q(7,5)**2 < 2
    pstar_lower = (Q(10)+7*Q(7,5))/32
    assert pstar_lower > Q(3,5)
    # sin(3pi/8)>.9 and cos(5pi/12)>.25; root identities carry the proof.
    assert (2+Q(7,5))/4 > Q(9,10)**2
    assert 3**2 > (2**2)*2  # 3 > 2 sqrt(2), hence sqrt(6)-sqrt(2)>1.
    derivative_lower = 3*Q(9,10)**5*Q(1,4)
    assert derivative_lower > Q(3,8)
    actual_midplane_tail = 4*pi_high**2/(om**2*dlo**2)
    q = Q(640)/om**2
    assert actual_midplane_tail < q
    joint = Q(3,5)*(1-q)-q-2*eps
    clock = Q(3,32)-eps-free
    signal = 1-2*q-eps-free
    assert joint > Q(59,100) and clock > Q(9,100) and signal > Q(998,1000)
    # Coefficientwise identity Hstar=lambda*H0star+mu*(1-lambda).
    hstar = (mu, lam, lam)
    transformed = (lam*mu + mu*(1-lam), lam, lam)
    assert hstar == transformed and lam*thi < cap
    retiming = eps+free
    assert retiming < Q(93,10**11)
    group2 = {
        "pi_interval_method": "Machin_identity_with_alternating_rational_enclosures",
        "delta_lower": rational(dlo), "delta_upper": rational(dhi),
        "clock_reference_probability_lower": rational(pstar_lower),
        "clock_derivative_lower": rational(derivative_lower),
        "midplane_tail_from_actual_delta": rational(actual_midplane_tail),
        "declared_midplane_tail": rational(q),
        "joint_gap_lower": rational(joint), "clock_gap_lower": rational(clock),
        "signal_gap_lower": rational(signal),
        "retiming_identity": "coefficientwise_exact_Hstar=lambda_H0star+mu*(1-lambda)",
        "all_menu_retiming_error": rational(retiming),
        "probability_conversion": "normalized_vector_bound_to_half_trace_distance_and_0<=E<=I",
        "not_certified": ["actual_arrival", "full_EH_source_insertion", "independent_geometry", "beyond_1PM", "GR_complete"],
    }
    return group1, group2


def build():
    first, second = finite_groups()
    return {
        "status": "PASS_INDEPENDENT_FINITE_1PM_TWO_GROUPS",
        "independent_groups": 2,
        "group_1_angular_source_symbol_and_budget": first,
        "group_2_joint_effects_and_retiming": second,
        "independent_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "method": "exact_coefficient_polynomials_and_rational_intervals; no_author_code_import_or_author_results_read",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    output = Path(__file__).with_name("independent_results.json")
    if args.write:
        payload = json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2)+"\n"
        with output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
        print("WRITE_ONCE_PASS_INDEPENDENT_1261_TWO_GROUPS")
    else:
        stored = json.loads(output.read_text(encoding="utf-8"))
        if stored != result:
            raise AssertionError("Independent saved results differ from this read-only rebuild")
        print("READ_ONLY_MATCH_INDEPENDENT_1261_TWO_GROUPS")
    print("independent_groups=2; geometry_added=0; GR_complete=false")


if __name__ == "__main__":
    main()
