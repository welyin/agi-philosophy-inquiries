"""Independent exact finite witnesses for 1067; no author-code imports.

Default: recompute and, if present, compare the saved independent result.
--write: exclusively create independent_results.json; never overwrite it.
Uses rational quaternions, rational planar fractional maps and exact ranks.
The global topological exclusion and all-domain inequalities require proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path


def vec(*xs):
    return tuple(F(x) for x in xs)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def qmul(q, p):
    return (q[0]*p[0]-dot(q[1:], p[1:]),) + add(add(scale(q[0], p[1:]), scale(p[0], q[1:])), cross(q[1:], p[1:]))


def rotate(q, v):
    # q need not be normalized: division by its exact squared norm suffices.
    qbar = (q[0],) + scale(-1, q[1:])
    result = qmul(qmul(q, (F(0),) + v), qbar)
    assert result[0] == 0
    return scale(1 / dot(q, q), result[1:])


def rank(rows):
    rows = [list(row) for row in rows]
    pivot = 0
    for col in range(len(rows[0])):
        found = next((i for i in range(pivot, len(rows)) if rows[i][col]), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        divisor = rows[pivot][col]
        rows[pivot] = [x / divisor for x in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot:
                factor = rows[i][col]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def solve(columns, value):
    rows = [[columns[j][i] for j in range(3)] + [value[i]] for i in range(3)]
    for col in range(3):
        k = next(i for i in range(col, 3) if rows[i][col])
        rows[col], rows[k] = rows[k], rows[col]
        divisor = rows[col][col]
        rows[col] = [x/divisor for x in rows[col]]
        for i in range(3):
            if i != col:
                c = rows[i][col]
                rows[i] = [x-c*y for x, y in zip(rows[i], rows[col])]
    return tuple(row[3] for row in rows)


def lincomb(columns, coefficients):
    return tuple(sum((columns[j][i]*coefficients[j] for j in range(3)), F(0)) for i in range(3))


def cmul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def cdiv(z, w):
    return scale(1/dot(w, w), cmul(z, (w[0], -w[1])))


def mobius(z, r):
    return cdiv((z[0]-r, z[1]), (1-r*z[0], -r*z[1]))


def effect(v, t=F(0), variable_trace=False):
    # E = c I + p dot sigma; all four coefficients are exact rationals.
    c = F(1, 2) + (F(1, 10)*t if variable_trace else 0)
    return (c,) + scale(F(3, 10), v)


def admissible(e):
    c, p = e[0], e[1:]
    return 0 <= c <= 1 and dot(p, p) <= c*c and dot(p, p) <= (1-c)*(1-c)


def diffrank(effects, traceless=False):
    differences = [add(e, scale(-1, effects[0])) for e in effects[1:]]
    return rank([row[1:] if traceless else row for row in differences])


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, (tuple, list)):
        return [encode(x) for x in value]
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    return value


def run():
    checks = []
    def verify(name, condition, data):
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "passed": True, "data": encode(data)})

    ex, ey, ez = vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)
    axes = [ex, ey, ez]
    qa, qb = vec(1, 0, 0, 1), vec(1, 1, 0, 0)
    ab, ba = rotate(qmul(qa, qb), ez), rotate(qmul(qb, qa), ez)
    visibility = F(3, 5)
    verify("S2_same_preparation_order_witness", ab == ex and ba == scale(-1, ey),
           {"ab_ez": ab, "ba_ez": ba, "preparation": "+(ex+ey)/sqrt(2)",
            "gap_exact": "(3/5)/sqrt(2)", "gap_squared": visibility**2/F(2),
            "gap_numeric": float(visibility)/math.sqrt(2)})

    # Continuous paths q_z(u), q_x(u), u in [0,1], yield rational rotations.
    parameters = [F(i, 8) for i in range(9)]
    path_quats = [vec(1, 0, 0, u) for u in parameters] + [vec(1, u, 0, 0) for u in parameters]
    identities = 0
    for q in path_quats + [qmul(qa, qb), qmul(qb, qa)]:
        images = [rotate(q, v) for v in axes]
        assert all(dot(images[i], images[j]) == F(i == j) for i in range(3) for j in range(3))
        assert cross(images[0], images[1]) == images[2]
        assert all(rotate(qmul(q, qb), v) == rotate(q, rotate(qb, v)) for v in axes)
        identities += 1
    verify("quaternion_paths_and_composition", identities == 20,
           {"exact_rotations": identities, "path_parameterization": "q=(1,u*axis), normalized implicitly; 0<=u<=1"})

    v1, v2 = ex, vec(F(3, 5), F(4, 5), 0)
    frame = [v1, v2, cross(v1, v2)]
    for q in path_quats:
        f1, f2 = rotate(q, v1), rotate(q, v2)
        image_frame = [f1, f2, cross(f1, f2)]
        for v in axes:
            assert lincomb(image_frame, solve(frame, v)) == rotate(q, v)
    verify("two_vector_plus_cross_representation_reconstruction", True,
           {"path_settings": len(path_quats), "probe_vectors": [v1, v2]})

    sphere_effects = [effect(scale(sign, axis)) for axis in axes for sign in [-1, 1]]
    verify("nonzero_fixed_trace_orbit_difference_rank", all(admissible(e) for e in sphere_effects) and diffrank(sphere_effects) == 3,
           {"full_coefficient_difference_rank": diffrank(sphere_effects), "trace": F(1)})

    radius, t = F(3, 5), F(4, 5)
    assert radius*radius+t*t == 1
    verify("S3_nonminimal_positive_contract", rotate(qmul(qa, qb), scale(radius, ez)) == scale(radius, ex),
           {"invariant_t": t, "radius": radius, "gap_exact": "(9/25)/sqrt(2)",
            "gap_squared": (visibility*radius)**2/F(2), "pole_orbits": "two fixed singleton orbits; action is not minimal"})
    north, south = effect(vec(0, 0, 0), F(1)), effect(vec(0, 0, 0), F(-1))
    verify("S3_constant_trace_antipodal_failure", north == south,
           {"north_effect": north, "south_effect": south, "antipodal_gap_at_poles": F(0)})

    variable_effects = sphere_effects + [effect(vec(0, 0, 0), t, True) for t in [-1, 1]]
    full_rank, zero_rank = diffrank(variable_effects), diffrank(variable_effects, True)
    verify("S3_variable_trace_changes_full_difference_space", full_rank == 4 and zero_rank == 3,
           {"full_difference_rank": full_rank, "traceless_projection_rank": zero_rank,
            "north_trace": F(6, 5), "south_trace": F(4, 5), "global_fixed_trace": False})

    antipodal_gaps = []
    for u in [F(i, 16) for i in range(-32, 33)]:
        r, t = abs(2*u/(1+u*u)), (1-u*u)/(1+u*u)
        assert r*r+t*t == 1
        for axis in axes:
            v = scale(r, axis)
            e = effect(v, t, True)
            anti = effect(scale(-1, v), -t, True)
            assert admissible(e) and admissible(anti)
            # Exact operator norm of aI+p.sigma is |a|+||p||.
            gap = F(1, 5)*abs(t) + F(3, 5)*r
            assert gap >= F(1, 5)
            antipodal_gaps.append(gap)
    verify("S3_variable_trace_antipodal_certificate", min(antipodal_gaps) == F(1, 5),
           {"samples": len(antipodal_gaps), "sample_minimum": min(antipodal_gaps),
            "global_exact_minimum": F(1, 5), "global_proof": "gap=(|t|+3r)/5 >= sqrt(t^2+r^2)/5=1/5; equality at r=0"})

    qreflection = vec(0, 1, 0, 0)
    oab, oba = rotate(qmul(qa, qreflection), ex), rotate(qmul(qreflection, qa), ex)
    verify("O2_circle_without_identity_reachability", oab == ey and oba == scale(-1, ey),
           {"ab_ex": oab, "ba_ex": oba, "same_preparation": "+Y", "probability_gap": visibility,
            "reflection_normal_sign": -1, "identity_component_normal_sign": 1,
            "missing_condition": "reflection is not reachable from identity inside O(2)"})
    assert rotate(qb, ey) == ez  # Ambient SO3 path leaves the equatorial interface.

    one, ii = vec(1, 0), vec(0, 1)
    r = F(1, 2)
    mab, mba = cmul(ii, mobius(one, r)), mobius(cmul(ii, one), r)
    gap = visibility*(mab[0]-mba[0])/2
    circle_points = [vec((1-u*u)/(1+u*u), 2*u/(1+u*u)) for u in parameters] + [vec(-1, 0)]
    for rr in [F(i, 16) for i in range(9)]:
        for z in circle_points:
            assert dot(mobius(z, rr), mobius(z, rr)) == 1
            assert mobius(mobius(z, rr), -rr) == z
    verify("Mobius_circle_without_global_unitary_covariance", mab == ii and mba == vec(F(-4, 5), F(3, 5)) and gap == F(6, 25),
           {"ab_1": mab, "ba_1": mba, "same_preparation": "+X", "probability_gap": gap,
            "inverse_identity_checks": 9*len(circle_points), "initial_direction_inner_product": F(0),
            "after_b_inner_product": dot(mobius(one, r), mobius(ii, r)),
            "missing_condition": "one common unitary must preserve every Bloch inner product"})

    return {"schema": "1067-independent-v1", "arithmetic": "exact fractions; irrational witness stored symbolically and squared exactly",
            "algorithm": "unnormalized rational quaternions and rational Mobius transformations; no author-code imports",
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "checks": checks, "passed": len(checks), "failed": 0,
            "scope": "Finite model witnesses and exact ranks only; not proof of maximum equicontinuous factor, circle quotient, or dimension theorem."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusively create independent_results.json")
    args = parser.parse_args()
    result = run()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        state = "exclusively saved"
    elif target.exists():
        assert json.loads(target.read_text(encoding="utf-8")) == result, "Saved independent result differs from recomputation"
        state = "read-only recomputation matches saved result"
    else:
        state = "read-only computation; no result file created"
    print(json.dumps({"passed": result["passed"], "failed": result["failed"], "status": state}, ensure_ascii=False))


if __name__ == "__main__":
    main()
