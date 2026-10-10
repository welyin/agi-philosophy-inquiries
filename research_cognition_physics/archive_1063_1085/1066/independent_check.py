"""1066 independent exact checks; default mode never writes.

Only finite rational circle and quaternion identities are computed.  Compact
neighborhoods, uniform convergence, inverse-limit paths and the Lie conclusion
must be proved analytically.  No author implementation is imported.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def demand(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def circle(x: F) -> F:
    shifted = x + F(1, 2)
    return x - shifted.numerator // shifted.denominator


def circle_root(x: F) -> F:
    return circle(x) / 2


def circle_distance(x: F, y: F) -> F:
    return abs(circle(x - y))


def circle_domain_boundary() -> dict:
    radius = F(9, 20)
    g = F(2, 5)
    square = circle(2*g)
    back = circle_root(square)
    demand(abs(g) < radius and abs(square) < radius and abs(back) < radius, "all points in the same declared U")
    demand(circle(2*back) == square, "right-root identity remains true")
    demand(back != g and circle_distance(back, g) == F(1, 2), "same-U left identity fails maximally")
    demand(square == -F(1, 5) and back == -F(1, 10), "exact branch witnesses")

    # On V=(-1/5,1/5), t, s=2t mod1 in V imply
    # |2t-s|<3/5<1, so the intervening integer is zero.
    small_radius = F(1, 5)
    demand(3*small_radius < 1, "analytic no-wrap bound on the small domain")
    tested = 0
    for numerator in range(-89, 90):
        h = F(numerator, 200)
        demand(abs(h) < radius, "grid belongs to U")
        root = circle_root(h)
        demand(abs(root) < radius and circle(2*root) == h, "finite right-root controls")
        if abs(h) < small_radius and abs(circle(2*h)) < small_radius:
            demand(circle_root(circle(2*h)) == h, "left identity in the smaller no-wrap domain")
            tested += 1
    return {
        "passed": True,
        "large_U_radius": str(radius),
        "g": str(g),
        "g_squared": str(square),
        "R_g_squared": str(back),
        "return_error": "1/2",
        "small_V_radius": str(small_radius),
        "right_root_grid_points": 179,
        "small_domain_left_identity_points": tested,
        "scope": "Continuity of the coordinate half-map and the slow-log strict cost decrease are analytic facts. A smaller-domain root identity is not asserted on the original U.",
    }


def circle_dyadic_chains() -> dict:
    seeds = (F(0), F(2, 5), -F(2, 5), F(1, 3), -F(1, 7), F(3, 7))
    integer_powers = (-9, -3, -1, 0, 1, 2, 7, 9)
    total = 0
    for seed in seeds:
        chain = [seed]
        for _ in range(12):
            chain.append(circle_root(chain[-1]))
        for n, value in enumerate(chain):
            demand(circle(2**n * value) == seed, "repeated roots reconstruct the original point")
        for m in range(len(chain)):
            for n in range(m, len(chain)):
                for k in integer_powers:
                    left = circle(k * 2**(n-m) * chain[n])
                    right = circle(k * chain[m])
                    demand(left == right, "dyadic descriptions agree without a left-root axiom")
                    total += 1
    return {
        "passed": True,
        "seeds": [str(seed) for seed in seeds],
        "maximum_root_depth": 12,
        "exact_refinement_identities": total,
        "scope": "Only the coherent chain h_(n+1)^2=h_n is used. Integer powers may leave U; no root selection is invoked outside U.",
    }


Quaternion = tuple[F, F, F, F]
ONE: Quaternion = (F(1), F(0), F(0), F(0))


def multiply(q: Quaternion, r: Quaternion) -> Quaternion:
    a, b, c, d = q
    e, f, g, h = r
    return (
        a*e - b*f - c*g - d*h,
        a*f + b*e + c*h - d*g,
        a*g - b*h + c*e + d*f,
        a*h + b*g - c*f + d*e,
    )


def power(q: Quaternion, exponent: int) -> Quaternion:
    if exponent < 0:
        q = (q[0], -q[1], -q[2], -q[3])
        exponent = -exponent
    answer = ONE
    while exponent:
        if exponent & 1:
            answer = multiply(answer, q)
        q = multiply(q, q)
        exponent //= 2
    return answer


def quaternion_chains() -> dict:
    # Rational parametrization of a unit quaternion.  The deepest root has
    # angle 2 atan(1/64); its 32nd power stays below angle 1<pi/2.
    s = F(1, 64)
    scalar = (1-s*s)/(1+s*s)
    vector = 2*s/(1+s*s)
    bases = ((scalar, vector, F(0), F(0)), (scalar, F(0), vector, F(0)))
    depth = 5
    chains = [[power(base, 2**(depth-n)) for n in range(depth+1)] for base in bases]
    checks = 0
    for chain in chains:
        for n, q in enumerate(chain):
            demand(sum((component*component for component in q), F(0)) == 1, "quaternion is exactly unit")
            demand(q[0] > 0, "positive-scalar local root branch")
            demand(power(q, 2**n) == chain[0], "full-root reconstruction")
            if n:
                demand(multiply(q, q) == chain[n-1], "adjacent coherent roots")
        for m in range(depth+1):
            for n in range(m, depth+1):
                for k in (-3, -1, 0, 1, 2, 5):
                    demand(power(chain[n], k*2**(n-m)) == power(chain[m], k), "noncommutative-group dyadic refinement")
                    checks += 1
    demand(multiply(chains[0][0], chains[1][0]) != multiply(chains[1][0], chains[0][0]), "ambient group is genuinely noncommutative")
    return {
        "passed": True,
        "rational_parameter": str(s),
        "root_depth": depth,
        "chains": 2,
        "exact_refinement_identities": checks,
        "two_seed_elements_commute": False,
        "scope": "Quaternion SU(2) controls use exact rational multiplication, no matrices, logarithms or floating roots. Each one-point root chain commutes internally; different seed elements need not commute.",
    }


def discontinuous_product_root_controls() -> dict:
    rows = []
    for denominator in (16, 64, 256, 1024):
        epsilon = F(1, denominator)
        left, right = F(1, 2)-epsilon, -F(1, 2)+epsilon
        demand(circle_distance(left, right) == 2*epsilon, "inputs coalesce at the circle cut")
        root_gap = circle_distance(circle_root(left), circle_root(right))
        demand(root_gap == F(1, 2)-epsilon, "root jump stays macroscopically nonzero")
        for coordinate in (1, 4, 11):
            weight = F(1, 2**coordinate)
            demand(weight*abs(circle_root(left)) == weight*abs(left)/2, "weighted cost halves on the left")
            demand(weight*abs(circle_root(right)) == weight*abs(right)/2, "weighted cost halves on the right")
        rows.append({
            "epsilon": str(epsilon),
            "input_circle_distance": str(2*epsilon),
            "root_circle_distance": str(root_gap),
        })

    # An arbitrary finite prefix and the uniform infinite-tail bound check
    # the ingredients, without treating a truncation as the infinite group.
    coordinates = (F(0), F(1, 2), -F(1, 3), F(2, 7), -F(4, 9), F(1, 8))
    cost = sum((F(1, 2**j)*abs(circle(x)) for j, x in enumerate(coordinates, 1)), F(0))
    root_cost = sum((F(1, 2**j)*abs(circle_root(x)) for j, x in enumerate(coordinates, 1)), F(0))
    demand(root_cost == cost/2 and cost > 0, "finite prefix exact weighted contraction")
    tail_bound = F(1, 2**(len(coordinates)+1))
    return {
        "passed": True,
        "cut_controls": rows,
        "finite_prefix_cost": str(cost),
        "finite_prefix_root_cost": str(root_cost),
        "uniform_tail_cost_bound": str(tail_bound),
        "analytic_scope": "On T^N, sum_(j>=1) 2^-j distance(x_j,0) is continuous and positive by uniform tail control; the coordinate principal root halves it exactly but is discontinuous at cuts in every identity neighborhood. Compactness, connectedness and arbitrarily small tail subgroups are analytic, not finite-prefix conclusions.",
    }


def run() -> dict:
    return {
        "round": 1066,
        "passed": True,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest().upper(),
        "method": "standard-library exact rational circle coordinates and unit quaternions",
        "circle_domain_boundary": circle_domain_boundary(),
        "circle_dyadic_chains": circle_dyadic_chains(),
        "quaternion_chains": quaternion_chains(),
        "discontinuous_product_root": discontinuous_product_root_controls(),
        "not_computed": [
            "Existence of an invariant compact neighborhood and uniform convergence of roots.",
            "Lie quotient compatibility and actual-group joint continuity of the constructed contraction.",
            "Local contractibility to Lie, or any spatial dimension/cognitive necessity.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusively create independent_results.json")
    args = parser.parse_args()
    payload = run()
    destination = Path(__file__).with_name("independent_results.json")
    if args.write:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        mode = "created"
    elif destination.exists():
        demand(json.loads(destination.read_text(encoding="utf-8")) == payload, "saved independent results match exactly")
        mode = "read_only_recomputed_and_matched"
    else:
        mode = "read_only_recomputed_no_saved_result"
    print(json.dumps({
        "passed": True,
        "mode": mode,
        "circle_same_U_return_error": payload["circle_domain_boundary"]["return_error"],
        "circle_refinement_identities": payload["circle_dyadic_chains"]["exact_refinement_identities"],
        "quaternion_refinement_identities": payload["quaternion_chains"]["exact_refinement_identities"],
        "cut_pairs": len(payload["discontinuous_product_root"]["cut_controls"]),
        "floating_point_used": False,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
