"""Adaptive circle partitions, local tangent triangles, and conditional certificates."""

import argparse
import copy
import hashlib
import json
import math
import sys
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, ceil_div, pi_interval, sin_interval, cos_interval
from compensated_angle_kernel import preconvolution_density
from outcome_updates import WINDOW_ATTENUATION as ALPHA


CASES = {
    "0144739": ("0.144739", 248, "181d68b7596fd6946d2b1c6c2a75b506fae8548fe09d3e643b86439b9ab5ae9a"),
    "014473923": ("0.14473923", 296, "c0afa5bd0a210a7fd9efbfc8553b659bde11d0167708b94329dec85aa0472e6e"),
    "0144739232": ("0.144739232", 324, "8a6c89646e7e800ce49639928a777cecd352a43a8754466dc2d5c727ccb1c1d4"),
}
BEST_CASE = "0144739232"


def candidate_path(case):
    return Path(__file__).with_name("adaptive_arc_candidate_" + case + ".json")


def read_candidate(case=BEST_CASE):
    return json.loads(candidate_path(case).read_text(encoding="utf-8"))


def adaptive_quarter(minimum=Fraction(1, 8192), growth=Fraction(6, 5), maximum_step=Fraction(1, 25)):
    minimum, growth, maximum_step = map(Fraction, (minimum, growth, maximum_step))
    if not 0 < minimum <= maximum_step < 1 or not 1 < growth <= 2:
        raise ValueError("Invalid adaptive spacing parameters.")
    distance = [Fraction(0), minimum]
    while distance[-1] < 1:
        value = min(distance[-1] * growth, distance[-1] + maximum_step, Fraction(1))
        distance.append(value)
        if len(distance) > 1024:
            raise ValueError("The requested partition is too large.")
    return tuple(sorted((1 - value) / 2 for value in distance))


def reflection(angle, element):
    return (angle, -angle, 1 - angle, 1 + angle)[element] % 2


@lru_cache(maxsize=8)
def partition(quarter):
    quarter = tuple(map(Fraction, quarter))
    if len(quarter) < 3 or quarter[0] != 0 or quarter[-1] != Fraction(1, 2):
        raise ValueError("Quarter boundaries must include 0 and 1/2 in units of pi.")
    if any(a >= b for a, b in zip(quarter, quarter[1:])):
        raise ValueError("Quarter boundaries must strictly increase.")
    boundary = tuple(sorted({reflection(a, g) for a in quarter for g in range(4)}))
    right = boundary[1:] + (Fraction(2),)
    midpoint = tuple((a + b) / 2 for a, b in zip(boundary, right))
    half_width = tuple((b - a) / 2 for a, b in zip(boundary, right))
    size = len(boundary)
    source_angles = boundary + midpoint
    source_lookup = {(i // size, angle): i for i, angle in enumerate(source_angles)}
    target_lookup = {angle: j for j, angle in enumerate(midpoint)}
    representatives = tuple(i for i, angle in enumerate(source_angles) if 0 <= angle <= Fraction(1, 2))
    target_representatives = tuple(j for j, angle in enumerate(midpoint) if 0 < angle < Fraction(1, 2))
    source_maps = [[source_lookup[(i // size, reflection(angle, g))] for i, angle in enumerate(source_angles)] for g in range(4)]
    target_maps = [[target_lookup[reflection(angle, g)] for angle in midpoint] for g in range(4)]
    canonical_target = {}
    for column, j in enumerate(target_representatives):
        for g in range(4):
            canonical_target[target_maps[g][j]] = column
    return {"quarter": quarter, "size": size, "boundary": boundary, "right": right,
            "midpoint": midpoint, "half_width": half_width, "mass": half_width,
            "source_angles": source_angles, "representatives": representatives,
            "target_representatives": target_representatives, "source_maps": source_maps,
            "target_maps": target_maps, "canonical_target": canonical_target}


def reduced_angle(angle):
    angle = Fraction(angle) % 2
    sin_sign, cos_sign = 1, 1
    if angle > 1:
        angle = 2 - angle
        sin_sign = -1
    if angle > Fraction(1, 2):
        angle = 1 - angle
        cos_sign = -1
    return angle, cos_sign, sin_sign


@lru_cache(maxsize=16384)
def trig_pi(angle):
    small, cos_sign, sin_sign = reduced_angle(angle)
    if small == 0:
        cosine, sine = I.exact(1), I.exact(0)
    elif small == Fraction(1, 2):
        cosine, sine = I.exact(0), I.exact(1)
    else:
        cosine, sine = cos_interval(pi_interval() * small), sin_interval(pi_interval() * small)
    return cos_sign * cosine, sin_sign * sine


def trig_pi_float(angle):
    small, cos_sign, sin_sign = reduced_angle(angle)
    if small == 0:
        return float(cos_sign), 0.0
    if small == Fraction(1, 2):
        return 0.0, float(sin_sign)
    value = math.pi * float(small)
    return cos_sign * math.cos(value), sin_sign * math.sin(value)


def sinc_interval(value):
    value = I.exact(value)
    if value.abs_upper().hi > 2 * SCALE:
        raise ValueError("Use |argument|<=2 for this sinc certificate.")
    total = I.exact(0)
    for k in range(24):
        total += Fraction((-1)**k, math.factorial(2 * k + 1)) * value**(2 * k)
    return total.widen(value.abs_upper()**48 / math.factorial(49))


@lru_cache(maxsize=4096)
def triangle_fourier_moment(width_pi, frequency):
    """Integral_0^1 q_A(u) exp(i*n*Delta*u) du, with a stable local series."""
    if type(frequency) is not int or frequency < 0 or frequency > 24 or frequency % 2:
        raise ValueError("Use an even integer frequency from 0 through 24.")
    delta = pi_interval() * width_pi
    if delta.lo <= 0 or delta.hi > I.rational(1, 10).lo:
        raise ValueError("Source cells must have angular width in (0, 0.1].")
    argument = max(1, frequency) * delta
    real, imaginary = I.exact(0), I.exact(0)
    for degree in range(36):
        power_sum = sum(frequency**(degree - 2 * j) for j in range(degree // 2 + 1))
        coefficient = Fraction(power_sum, max(1, frequency)**degree * math.factorial(degree + 3))
        term = coefficient * argument**degree
        if degree % 4 == 0:
            real += term
        elif degree % 4 == 1:
            imaginary += term
        elif degree % 4 == 2:
            real -= term
        else:
            imaginary -= term
    tail = Fraction(4, 3) * argument**36 / math.factorial(39) / (1 - argument / 40)
    prefactor = 2 / sinc_interval(delta / 2)**2
    return (prefactor * real).widen(prefactor * tail), (prefactor * imaginary).widen(prefactor * tail)


def geometry_float(quarter, contrast):
    p = partition(quarter)
    size = p["size"]
    points = np.array([trig_pi_float(a) for a in p["source_angles"]]) * ALPHA
    for j, width in enumerate(p["half_width"]):
        points[size + j] /= math.cos(math.pi * float(width))
    targets = np.array([trig_pi_float(a) for a in p["midpoint"]])
    widths = np.array(list(map(float, p["half_width"])))
    targets *= np.sinc(widths)[:, None]
    weights = np.zeros(2 * size)
    nodes, gauss_weights = np.polynomial.legendre.leggauss(32)
    u = (nodes + 1) / 2
    for j, (a, width) in enumerate(zip(p["boundary"], widths)):
        d = math.pi * width
        left = np.sin(d * (1 - u)) / math.sin(d)
        right = np.sin(d * u) / math.sin(d)
        density = preconvolution_density(math.pi * float(a) + 2 * d * u, contrast)
        weights[j] += width * np.sum(gauss_weights * left**2 * density) / 2
        weights[(j + 1) % size] += width * np.sum(gauss_weights * right**2 * density) / 2
        weights[size + j] = width * np.sum(gauss_weights * 2 * math.cos(d) * left * right * density) / 2
    for i in p["representatives"]:
        orbit = sorted({p["source_maps"][g][i] for g in range(4)})
        weights[orbit] = weights[orbit].mean()
    weights /= weights.sum()
    return dict(p, weights=weights, sources=points, targets=targets, masses=widths)


@lru_cache(maxsize=4)
def geometry_interval(quarter, contrast):
    p = partition(quarter)
    size = p["size"]
    pi, alpha = pi_interval(), 4 * sin_interval(I.rational(1, 4))
    eta = I.exact(contrast)
    if eta.lo < 0 or eta.hi >= I.rational(1, 5).lo:
        raise ValueError("Use 0<=eta<0.2 in this certificate.")
    q = (1 - eta**2).sqrt()
    ratio = (eta / (1 + q))**2
    coefficients = [ratio**k * (1 + 2 * k * q) for k in range(1, 13)]
    tail = ratio**13 * ((1 + 26 * q) / (1 - ratio) + 2 * q * ratio / (1 - ratio)**2)
    weights = [I.exact(0) for _ in range(2 * size)]
    sources = [(alpha * c, alpha * s) for c, s in map(trig_pi, p["source_angles"])]
    targets, radius_margins = [], []
    for j, (a, b, middle, half_width) in enumerate(zip(p["boundary"], p["right"], p["midpoint"], p["half_width"])):
        d = pi * half_width
        beta = sinc_interval(d)
        cosine, sine = trig_pi(middle)
        targets.append((beta * cosine, beta * sine))
        tangent_radius = alpha / cos_interval(d)
        sources[size + j] = (tangent_radius * cosine, tangent_radius * sine)
        radius_margins.append(1 - tangent_radius)
        zeroth, _ = triangle_fourier_moment(b - a, 0)
        left, right, mass = zeroth, zeroth, I.exact(1)
        for k, coefficient in enumerate(coefficients, 1):
            frequency = 2 * k
            real, imaginary = triangle_fourier_moment(b - a, frequency)
            ca, sa = trig_pi(frequency * a)
            cb, sb = trig_pi(frequency * b)
            cm, _ = trig_pi(frequency * middle)
            left += 2 * coefficient * (ca * real - sa * imaginary)
            right += 2 * coefficient * (cb * real + sb * imaginary)
            mass += 2 * coefficient * sinc_interval(frequency * d) * cm
        remainder = 2 * tail * half_width
        wa, wb = (left * half_width).widen(remainder), (right * half_width).widen(remainder)
        # q_A+q_B+q_C=1; only one Fourier tail is needed for the third weight.
        wc = ((mass - left - right) * half_width).widen(remainder)
        weights[j] += wa
        weights[(j + 1) % size] += wb
        weights[size + j] = wc
    cx = sum(mass * target[0]**2 for mass, target in zip(p["mass"], targets))
    cy = sum(mass * target[1]**2 for mass, target in zip(p["mass"], targets))
    return dict(p, weights=weights, sources=sources, targets=targets, variance_x=cx, variance_y=cy,
                minimum_radius_margin=I(min(value.lo for value in radius_margins)), fourier_tail=tail)


def expand_representatives(p, rows):
    size = p["size"]
    output = np.zeros((2 * size, size))
    for row, i in zip(rows, p["representatives"]):
        stabilizer = 4 / len({p["source_maps"][g][i] for g in range(4)})
        for g in range(4):
            output[p["source_maps"][g][i], p["target_maps"][g]] += row / stabilizer
    return output


def solve_adaptive(quarter, contrast, time_limit=180):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / ".research_runtime"))
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    geometry = geometry_float(quarter, contrast)
    size = geometry["size"]
    representatives = geometry["representatives"]
    count, columns_count = len(representatives), len(geometry["target_representatives"])
    index = np.arange(count * size).reshape(count, size)
    row_indices, column_indices, values = [], [], []
    rhs = []
    constraint = 0
    for k, i in enumerate(representatives):
        row_indices.extend([constraint] * (size + 1))
        column_indices.extend(index[k].tolist() + [count * size])
        values.extend([1.] * (size + 1))
        rhs.append(1.)
        constraint += 1
        for coordinate in (0, 1):
            if geometry["sources"][i, coordinate] == 0:
                continue  # The stabilizer average forces this coordinate to zero.
            row_indices.extend([constraint] * size)
            column_indices.extend(index[k].tolist())
            values.extend(geometry["targets"][:, coordinate].tolist())
            rhs.append(geometry["sources"][i, coordinate])
            constraint += 1
    column_start = constraint
    for k, i in enumerate(representatives):
        orbit_size = len({geometry["source_maps"][g][i] for g in range(4)})
        for j in range(size):
            row_indices.append(column_start + geometry["canonical_target"][j])
            column_indices.append(int(index[k, j]))
            values.append(geometry["weights"][i] * orbit_size / (4 * geometry["masses"][j]))
    row_indices.extend(range(column_start, column_start + columns_count))
    column_indices.extend([count * size] * columns_count)
    values.extend([1.] * columns_count)
    rhs.extend([1.] * columns_count)
    matrix = coo_matrix((values, (row_indices, column_indices)), shape=(len(rhs), count * size + 1)).tocsr()
    objective = np.zeros(count * size + 1)
    objective[-1] = -1
    result = linprog(objective, A_eq=matrix, b_eq=rhs,
                     bounds=[(0, None)] * (count * size) + [(0, 1)], method="highs-ipm",
                     options={"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10,
                              "ipm_optimality_tolerance": 1e-12, "time_limit": time_limit})
    report = {"contrast": contrast, "target_labels": size, "source_points": 2 * size,
              "variables": count * size + 1, "constraints": len(rhs),
              "status": int(result.status), "message": result.message}
    if not result.success:
        return report, None
    raw = result.x[:-1].reshape(count, size)
    report.update(reserve=float(result.x[-1]), maximum_scaled_constraint_residual=float(np.max(np.abs(matrix @ result.x - rhs))),
                  negative_Q_entries_clipped=int(np.count_nonzero(raw < 0)))
    full = expand_representatives(geometry, np.maximum(raw, 0))
    bits = 100
    candidate = {"schema": "adaptive_conditional_v1", "report": report,
                 "quarter_boundary_pi": list(map(str, quarter)), "denominator_bits": bits,
                 "reserve_numerator": int(result.x[-1] * (1 << bits)),
                 "conditional_sparse": [[int(i), int(j), int(full[i, j] * (1 << bits))] for i, j in zip(*np.nonzero(full))]}
    return report, candidate


def unpack_candidate(data):
    if data["schema"] != "adaptive_conditional_v1" or data["denominator_bits"] != 100:
        raise ValueError("Unexpected conditional candidate format.")
    quarter = tuple(map(Fraction, data["quarter_boundary_pi"]))
    p = partition(quarter)
    size = p["size"]
    if (data["report"]["target_labels"], data["report"]["source_points"]) != (size, 2 * size):
        raise ValueError("Candidate dimensions do not match its partition.")
    reserve = data["reserve_numerator"]
    if type(reserve) is not int or not 0 < reserve < SCALE:
        raise ValueError("Expected a strictly positive dyadic reserve below one.")
    rows, columns, seen = [[] for _ in range(2 * size)], [[] for _ in range(size)], set()
    for entry in data["conditional_sparse"]:
        if len(entry) != 3 or any(type(value) is not int for value in entry):
            raise ValueError("Invalid sparse entry.")
        i, j, value = entry
        if not 0 <= i < 2 * size or not 0 <= j < size or value < 0 or (i, j) in seen:
            raise ValueError("Out of range, negative or repeated sparse entry.")
        rows[i].append((j, value))
        columns[j].append((i, value))
        seen.add((i, j))
    return quarter, rows, columns


def interval_integer_dot(entries, points, coordinate):
    lower = sum(value * points[j][coordinate].lo for j, value in entries)
    upper = sum(value * points[j][coordinate].hi for j, value in entries)
    return I(lower // SCALE, ceil_div(upper, SCALE))


def verify_adaptive(data, target_contrast=None):
    quarter, rows, columns = unpack_candidate(data)
    contrast = Fraction(str(data["report"]["contrast"])) if target_contrast is None else Fraction(target_contrast)
    geometry = geometry_interval(quarter, contrast)
    reserve = I(data["reserve_numerator"])
    weights, targets = geometry["weights"], geometry["targets"]
    if min(weight.lo for weight in weights) <= 0 or geometry["minimum_radius_margin"].lo <= 0:
        return {"verified": False, "reason": "Source geometry not certified strictly positive."}
    cx, cy = geometry["variance_x"], geometry["variance_y"]
    if cx.lo <= 0 or cy.lo <= 0:
        return {"verified": False, "reason": "Target covariance not certified positive."}
    mx = I(max(target[0].abs_upper().hi for target in targets))
    my = I(max(target[1].abs_upper().hi for target in targets))
    row_bounds = []
    for row, source in zip(rows, geometry["sources"]):
        mass_error = 1 - reserve - I(sum(value for _, value in row))
        bx = source[0] - interval_integer_dot(row, targets, 0)
        by = source[1] - interval_integer_dot(row, targets, 1)
        bound = mass_error.abs_upper() + bx.abs_upper() * mx / I(cx.lo) + by.abs_upper() * my / I(cy.lo)
        row_bounds.append(I(bound.hi))
    column_relative_errors = []
    for column, mass in zip(columns, geometry["mass"]):
        lower = sum(value * weights[i].lo for i, value in column)
        upper = sum(value * weights[i].hi for i, value in column)
        total = I(lower // SCALE, ceil_div(upper, SCALE))
        error = ((1 - reserve) * mass - total).abs_upper() / mass
        column_relative_errors.append(I(error.hi))
    average_row_bound = sum(weight * bound for weight, bound in zip(weights, row_bounds))
    row_max = I(max(value.hi for value in row_bounds))
    column_max = I(max(value.hi for value in column_relative_errors))
    margin = reserve - row_max - column_max - average_row_bound
    if margin.lo > 0:
        minimum_joint = Fraction(margin.lo, SCALE) * Fraction(min(w.lo for w in weights), SCALE) * min(geometry["mass"])
    else:
        minimum_joint = Fraction(margin.lo, SCALE) * Fraction(max(w.hi for w in weights), SCALE) * max(geometry["mass"])
    return {"verified": margin.lo > 0, "contrast": str(contrast), "target_labels": geometry["size"],
            "source_points": len(weights), "minimum_target_mass": str(min(geometry["mass"])),
            "reserve_interval": reserve.floats(), "maximum_normalized_row_repair": row_max.hi / SCALE,
            "maximum_relative_column_residual": column_max.hi / SCALE,
            "weighted_row_repair_bound": average_row_bound.hi / SCALE,
            "corrected_product_reserve_lower_bound_exact": str(Fraction(margin.lo, SCALE)),
            "corrected_product_reserve_lower_bound_float": margin.lo / SCALE,
            "minimum_joint_entry_lower_bound_exact": str(minimum_joint),
            "minimum_joint_entry_lower_bound_float": float(minimum_joint),
            "minimum_tangent_radius_margin": geometry["minimum_radius_margin"].lo / SCALE,
            "scope": ("Local tangent-triangle quantization, nonuniform uniform-arc lifting, and exact conditional/covariance repair certify the continuous canonical instrument."
                      if margin.lo > 0 else "This sufficient repair bound failed. This is not an infeasibility proof.")}


def repaired_conditional_float(data, target_contrast=None):
    quarter, rows, _ = unpack_candidate(data)
    contrast = float(data["report"]["contrast"] if target_contrast is None else target_contrast)
    geometry = geometry_float(quarter, contrast)
    reserve = data["reserve_numerator"] / SCALE
    matrix = np.tile(reserve * geometry["masses"], (2 * geometry["size"], 1))
    for i, row in enumerate(rows):
        for j, value in row:
            matrix[i, j] += value / SCALE
    mass_error = 1 - matrix.sum(axis=1)
    moment_error = geometry["sources"] - matrix @ geometry["targets"]
    covariance = (geometry["targets"].T * geometry["masses"]) @ geometry["targets"]
    matrix += geometry["masses"] * (mass_error[:, None] + moment_error @ np.linalg.solve(covariance, geometry["targets"].T))
    column_error = geometry["masses"] - geometry["weights"] @ matrix
    matrix += column_error
    return geometry, matrix


def quantize_angles(angles, geometry):
    theta = np.mod(np.asarray(angles), 2 * math.pi)
    boundary = math.pi * np.array(list(map(float, geometry["boundary"])))
    cell = np.searchsorted(boundary, theta, side="right") - 1
    width = math.pi * np.array(list(map(float, geometry["half_width"])))[cell]
    offset = theta - boundary[cell]
    u = offset / (2 * width)
    left = np.sin(width * (1 - u)) / np.sin(width)
    right = np.sin(width * u) / np.sin(width)
    return cell, left**2, right**2, 2 * np.cos(width) * left * right


@lru_cache(maxsize=1)
def endpoint_bracket():
    """Bracket alpha*G(eta)=1; every bisection sign uses outward intervals."""
    from polygon_limit_analysis import factor_interval
    alpha = 4 * sin_interval(I.rational(1, 4))
    lower, upper = Fraction(0), Fraction(1, 5)
    equation = lambda value: alpha * factor_interval(value) - 1
    if equation(lower).hi >= 0 or equation(upper).lo <= 0:
        raise ArithmeticError("Missing strict endpoint bracket.")
    for _ in range(64):
        middle = (lower + upper) / 2
        value = equation(middle)
        if value.hi < 0:
            lower = middle
        elif value.lo > 0:
            upper = middle
        else:
            raise ArithmeticError("Increase interval precision before continuing bisection.")
    return lower, upper


@lru_cache(maxsize=3)
def saved_float_model(case=BEST_CASE):
    return repaired_conditional_float(read_candidate(case))


class AdaptiveArcBridgeTests(unittest.TestCase):
    def test_partition_is_a_reflection_symmetric_probability_partition(self):
        for case in CASES:
            quarter, _, _ = unpack_candidate(read_candidate(case))
            p = partition(quarter)
            self.assertEqual(sum(p["mass"]), 1)
            self.assertEqual(p["size"], CASES[case][1])
            self.assertTrue(all(axis in p["boundary"] for axis in (0, Fraction(1, 2), 1, Fraction(3, 2))))
            for g in range(4):
                self.assertEqual(sorted(p["source_maps"][g]), list(range(2 * p["size"])))
                self.assertEqual(sorted(p["target_maps"][g]), list(range(p["size"])))
                self.assertTrue(all(p["mass"][j] == p["mass"][p["target_maps"][g][j]] for j in range(p["size"])))
        with self.assertRaises(ValueError):
            adaptive_quarter(growth=1)
        with self.assertRaises(ValueError):
            partition((Fraction(0), Fraction(1, 2), Fraction(1, 2)))

    def test_local_triangle_preserves_mass_vector_and_both_absolute_coordinates(self):
        geometry, _ = saved_float_model()
        # Include every boundary and points inside the smallest cells.
        boundaries = np.array(list(map(float, geometry["boundary"])))
        widths = np.array(list(map(float, geometry["half_width"])))
        angles = math.pi * (boundaries[:, None] + 2 * widths[:, None] * np.array([0, .001, .27, .5, .999])).ravel()
        angles = np.r_[angles, -2 * math.pi, 2 * math.pi, 4 * math.pi]
        cell, qa, qb, qc = quantize_angles(angles, geometry)
        self.assertGreaterEqual(min(qa.min(), qb.min(), qc.min()), -2e-15)
        np.testing.assert_allclose(qa + qb + qc, 1, atol=2e-15, rtol=0)
        points = geometry["sources"]
        size = geometry["size"]
        for absolute in (False, True):
            vertices = np.abs(points) if absolute else points
            actual = qa[:, None] * vertices[cell] + qb[:, None] * vertices[(cell + 1) % size] + qc[:, None] * vertices[size + cell]
            expected = ALPHA * np.column_stack((np.cos(angles), np.sin(angles)))
            np.testing.assert_allclose(actual, np.abs(expected) if absolute else expected, atol=3e-15, rtol=0)

    def test_stable_fourier_intervals_contain_independent_decimal_closed_form(self):
        # This reference integrates the exponential explicitly; it does not use
        # the local moment series or the Fourier density used by the certificate.
        with localcontext() as context:
            context.prec = 90
            def atan_inverse(n):
                x = Decimal(1) / n
                return sum((-1)**k * x**(2 * k + 1) / (2 * k + 1) for k in range(160))
            pi = 16 * atan_inverse(5) - 4 * atan_inverse(239)
            def trig(x):
                return (sum((-1)**k * x**(2 * k) / Decimal(math.factorial(2 * k)) for k in range(80)),
                        sum((-1)**k * x**(2 * k + 1) / Decimal(math.factorial(2 * k + 1)) for k in range(80)))
            def multiply(a, b):
                return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]
            for width in (Fraction(1, 50), Fraction(1, 81920), Fraction(1, 2**24)):
                delta = pi * Decimal(width.numerator) / width.denominator
                c, s = trig(delta)
                def exponential_integral(n):
                    if n == 0:
                        return Decimal(1), Decimal(0)
                    cn, sn = trig(n * delta)
                    return sn / (n * delta), (1 - cn) / (n * delta)
                for n in (0, 2, 24):
                    first = exponential_integral(n)
                    left = multiply((c, s), exponential_integral(n - 1))
                    right = multiply((c, -s), exponential_integral(n + 1))
                    answer = [(first[k] - (left[k] + right[k]) / 2) / (1 - c) for k in (0, 1)]
                    for interval, expected in zip(triangle_fourier_moment(width, n), answer):
                        self.assertLessEqual(Decimal(interval.lo) / SCALE, expected)
                        self.assertGreaterEqual(Decimal(interval.hi) / SCALE, expected)

    def test_interval_source_weights_agree_with_independent_positive_quadrature(self):
        from convex_order_interface import absolute_factor
        for case, (eta, _, _) in CASES.items():
            quarter, _, _ = unpack_candidate(read_candidate(case))
            interval = geometry_interval(quarter, Fraction(eta))
            geometry, _ = saved_float_model(case)
            weights = np.array([(w.lo + w.hi) / (2 * SCALE) for w in interval["weights"]])
            np.testing.assert_allclose(weights, geometry["weights"], atol=1e-16, rtol=2e-13)
            self.assertAlmostEqual(weights.sum(), 1, places=14)
            np.testing.assert_allclose(weights @ geometry["sources"], 0, atol=3e-16)
            self.assertAlmostEqual(float(weights @ np.abs(geometry["sources"][:, 0])),
                                   2 * ALPHA * absolute_factor(float(eta)) / math.pi, places=14)

    def test_axis_split_target_has_no_absolute_projection_loss(self):
        for case in CASES:
            geometry, _ = saved_float_model(case)
            np.testing.assert_allclose(geometry["masses"] @ np.abs(geometry["targets"]), [2 / math.pi] * 2, atol=3e-16)
            np.testing.assert_allclose(geometry["masses"] @ geometry["targets"], 0, atol=2e-16)
            covariance = (geometry["targets"].T * geometry["masses"]) @ geometry["targets"]
            self.assertLess(abs(covariance[0, 1]), 1e-16)
            self.assertGreater(np.linalg.eigvalsh(covariance).min(), .49)

    def test_all_archived_candidates_pass_strict_independent_certification(self):
        for case, (eta, labels, digest) in CASES.items():
            self.assertEqual(hashlib.sha256(candidate_path(case).read_bytes()).hexdigest(), digest)
            report = verify_adaptive(read_candidate(case))
            self.assertTrue(report["verified"])
            self.assertEqual(Fraction(report["contrast"]), Fraction(eta))
            self.assertEqual(report["target_labels"], labels)
            self.assertGreater(Fraction(report["minimum_joint_entry_lower_bound_exact"]), 0)

    def test_damaged_and_malformed_candidates_are_rejected(self):
        data = read_candidate()
        damaged = copy.deepcopy(data)
        damaged["conditional_sparse"][0][2] += 1 << 85
        self.assertFalse(verify_adaptive(damaged)["verified"])
        for kind in ("negative", "duplicate", "dimensions", "reserve"):
            damaged = copy.deepcopy(data)
            if kind == "negative":
                damaged["conditional_sparse"][0][2] = -1
            elif kind == "duplicate":
                damaged["conditional_sparse"].append(damaged["conditional_sparse"][0])
            elif kind == "dimensions":
                damaged["report"]["source_points"] += 1
            else:
                damaged["reserve_numerator"] = 0
            with self.assertRaises(ValueError):
                unpack_candidate(damaged)

    def test_repaired_float_models_have_positive_entries_and_complete_marginals(self):
        for case in CASES:
            geometry, conditional = saved_float_model(case)
            self.assertGreater(conditional.min(), 0)
            np.testing.assert_allclose(conditional.sum(axis=1), 1, atol=2e-15, rtol=0)
            np.testing.assert_allclose(geometry["weights"] @ conditional, geometry["masses"], atol=1e-16, rtol=2e-13)
            np.testing.assert_allclose(conditional @ geometry["targets"], geometry["sources"], atol=2e-15, rtol=0)

    def test_full_selective_density_for_nonuniform_arcs(self):
        from compensated_angle_kernel import warped_angle
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        outgoing = np.r_[np.linspace(-3, 9, 81), 0, math.pi / 2, math.pi, 3 * math.pi / 2]
        for case, (eta_string, _, _) in CASES.items():
            eta = float(eta_string)
            geometry, conditional = saved_float_model(case)
            size = geometry["size"]
            for setting in (0, .29, -1.7):
                rotation = np.array([[math.cos(setting), -math.sin(setting)], [math.sin(setting), math.cos(setting)]])
                targets = geometry["targets"] @ rotation.T
                for mark in (0, 1):
                    relative = outgoing - setting
                    cell, qa, qb, qc = quantize_angles(warped_angle(relative, mark, eta), geometry)
                    probabilities = qa[:, None] * conditional[cell] + qb[:, None] * conditional[(cell + 1) % size] + qc[:, None] * conditional[size + cell]
                    response = (1 + (2 * mark - 1) * eta * np.cos(relative)) / 2
                    for summary in ([1, ALPHA, 0], [1, ALPHA / math.sqrt(2), ALPHA / math.sqrt(2)], [.4, .1, -.2]):
                        arc_mass = geometry["masses"] * (summary[0] + (targets @ summary[1:]) / ALPHA)
                        observed = response / (2 * math.pi) * (probabilities @ (arc_mass / geometry["masses"]))
                        expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, outgoing)
                        np.testing.assert_allclose(observed, expected, atol=3e-14, rtol=0)

    def test_strict_endpoint_gap_and_old_fixed_grid_exclusion(self):
        from polygon_limit_analysis import grid_certificate
        lower, upper = endpoint_bracket()
        eta = Fraction(CASES[BEST_CASE][0])
        self.assertGreater(lower, eta)
        self.assertLess(upper - eta, Fraction(1, 10**10))
        self.assertLess(upper - lower, Fraction(1, 10**19))
        self.assertTrue(grid_certificate(2048, eta)["grid_excluded"])
        # This comparison is only for the old, pole-centered equal-arc family.
        self.assertTrue(grid_certificate(600688, eta)["grid_excluded"])
        self.assertGreater(int(grid_certificate(600692, eta)["budget_integer_interval"][0]), 0)


def results_report(checks):
    from finite_bridge_cost import necessary_label_count
    from polygon_limit_analysis import grid_certificate
    lower, upper = endpoint_bracket()
    reports = []
    for case, (eta, _, digest) in CASES.items():
        report = verify_adaptive(read_candidate(case))
        report.update(candidate=candidate_path(case).name, candidate_sha256=digest,
                      endpoint_gap_bracket_exact=[str(lower - Fraction(eta)), str(upper - Fraction(eta))],
                      endpoint_gap_diagnostic=float((lower + upper) / 2 - Fraction(eta)),
                      necessary_label_count_from_round_42=necessary_label_count(upper - Fraction(eta)))
        reports.append(report)
    eta = Fraction(CASES[BEST_CASE][0])
    return {"certificates": reports, "endpoint_bracket_exact": list(map(str, (lower, upper))),
            "canonical_family_lower_bound_exact": str(Fraction(CASES[BEST_CASE][0])),
            "old_pole_centered_grid_projection_comparison": {
                "excluded": grid_certificate(600688, eta), "first_not_excluded_multiple_of_four": grid_certificate(600692, eta),
                "scope": "Monotonicity of cos(pi/N) excludes every smaller admissible old grid; a positive projection budget is not sufficient for existence."},
            "endpoint_feasibility_proven": False, "logarithmic_sufficient_label_bound_proven": False,
            "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--contrast", type=str)
    parser.add_argument("--minimum", type=str, default="1/8192")
    parser.add_argument("--growth", type=str, default="6/5")
    parser.add_argument("--maximum-step", type=str, default="1/25")
    parser.add_argument("--candidate-output", type=Path)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    if sum((args.solve, args.verify is not None, args.write_results)) > 1:
        parser.error("Choose a single operation.")
    if args.solve:
        if args.candidate_output is None or args.candidate_output.exists():
            parser.error("Choose a new candidate output file.")
        quarter = adaptive_quarter(Fraction(args.minimum), Fraction(args.growth), Fraction(args.maximum_step))
        report, candidate = solve_adaptive(quarter, float(args.contrast or CASES["0144739"][0]))
        print(json.dumps(report), flush=True)
        if candidate is None:
            raise SystemExit(1)
        certificate = verify_adaptive(candidate)
        print(json.dumps(certificate), flush=True)
        if not certificate["verified"]:
            raise SystemExit(1)
        args.candidate_output.write_text(json.dumps(candidate, separators=(",", ":")) + "\n", encoding="utf-8")
    elif args.verify:
        candidate = json.loads(args.verify.read_text(encoding="utf-8"))
        report = verify_adaptive(candidate, Fraction(args.contrast) if args.contrast is not None else None)
        print(json.dumps(report, indent=2))
        if not report["verified"]:
            raise SystemExit(1)
    else:
        checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AdaptiveArcBridgeTests))
        if not checks.wasSuccessful():
            raise SystemExit(1)
        report = results_report(checks)
        if args.write_results:
            Path(__file__).with_name("adaptive_arc_bridge_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
