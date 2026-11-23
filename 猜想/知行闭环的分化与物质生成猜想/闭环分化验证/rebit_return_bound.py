"""Rational polygon certificate for the delayed-return task at p = 4/5.

For four planar Bloch vectors a,b,c,d the optimized score is 7/10 + G/80.
G is a sum of six Euclidean norms given in polygon_upper_bound. It is
jointly convex and positively homogeneous. Thus if r D is in the polygon
P and D is the unit disk, max_(D^4) G <= max_(vertices(P)^4) G / r.
"""

import argparse
import json
from fractions import Fraction
from math import isqrt, lcm
from pathlib import Path


def rational_polygon():
    # The rational parametrization at t=j/8 enumerates the first quadrant
    # in angle order. Append its exact rotations to cover the full circle.
    quadrant = []
    for j in range(8):
        denominator = 64 + j * j
        quadrant.append((Fraction(64 - j * j, denominator),
                         Fraction(16 * j, denominator)))
    points = []
    for _ in range(4):
        points.extend(quadrant)
        quadrant = [(-y, x) for x, y in quadrant]
    scale = lcm(*(v.denominator for p in points for v in p))
    vertices = [(int(scale * x), int(scale * y)) for x, y in points]
    assert all(x * x + y * y == scale * scale for x, y in vertices)
    # Distance from origin to a chord of unit vectors u,v is
    # sqrt((1 + u.v)/2). Check an explicit rational inscribed radius.
    radius = Fraction(127, 128)
    for (x, y), (u, v) in zip(vertices, vertices[1:] + vertices[:1]):
        assert x * v - y * u > 0
        squared_distance = Fraction(scale * scale + x * u + y * v, 2 * scale * scale)
        assert squared_distance >= radius * radius
    return vertices, scale, radius


def ceil_norm(x, y):
    squared = x * x + y * y
    return 1 + isqrt(squared - 1) if squared else 0


def polygon_upper_bound(vertices):
    # All coordinates are integers in units 1/scale. Rounding each norm
    # upward gives a bound on each tuple, independently of floating point.
    distances = [[ceil_norm(x - u, y - v) for u, v in vertices] for x, y in vertices]
    maximum = 0
    checked = 0
    # Simultaneously rotating four vectors by 90 degrees preserves all terms.
    for i, (ax, ay) in enumerate(vertices[:len(vertices) // 4]):
        for j, (bx, by) in enumerate(vertices):
            for k, (cx, cy) in enumerate(vertices):
                common_edges = distances[i][j] + distances[i][k]
                for m, (dx, dy) in enumerate(vertices):
                    first = ceil_norm(ax + bx - cx - dx, ay + by - cy - dy)
                    second = ceil_norm(ax + cx - bx - dx, ay + cy - by - dy)
                    edges = common_edges + distances[j][m] + distances[k][m]
                    maximum = max(maximum, first + second + 2 * edges)
                    checked += 1
    return maximum, checked


def compute():
    vertices, scale, radius = rational_polygon()
    maximum, checked = polygon_upper_bound(vertices)
    score = Fraction(7, 10) + Fraction(maximum, 80 * scale) / radius
    threshold = Fraction(9235, 10000)
    assert score < threshold
    # The three-dimensional witness gives 7/10 + sqrt(5)/10.
    assert (10 * (threshold - Fraction(7, 10))) ** 2 < 5
    return {
        "return_probability": "4/5",
        "rational_polygon_vertices": len(vertices),
        "coordinate_denominator": scale,
        "certified_radius_lower": str(radius),
        "vertex_tuples_checked": checked,
        "polygon_G_upper_scaled": maximum,
        "rebit_score_upper": str(score),
        "rebit_score_upper_ceiling_millionths": -(-score.numerator * 10**6 // score.denominator),
        "separating_threshold": str(threshold),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = compute()
    path = Path(__file__).with_name("rebit_return_bound_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == results
    else:
        with path.open("x", encoding="utf-8") as out:
            out.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
