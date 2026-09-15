"""Exact covariance repair inside a half-circle block, preserving forbidden zeros."""

import argparse
import hashlib
import json
import math
import unittest
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, ceil_div, pi_interval, cos_interval
from polygon_martingale_certificate import certified_geometry
from polygon_limit_analysis import factor_interval
from half_polygon_reduction import half_geometry, expand_half_coupling, solve_half_polygon


CASES = {"512": (512, .1445, .989635), "2048": (2048, .14473, .98961702)}
RECERTIFIED_CONTRAST = Fraction("0.14473104667")


def verification_cases():
    return [(case, None) for case in CASES] + [("2048", RECERTIFIED_CONTRAST)]


def candidate_path(case):
    return Path(__file__).with_name("half_face_candidate_" + case + ".json")


def pack_candidate(report, half, bits=80):
    baselines, exceptions = [], []
    for i, row in enumerate(half):
        values = [int(value * (1 << bits)) for value in row]
        # Pole columns have a smaller product background than interior columns.
        # The most frequent value keeps both backgrounds losslessly compact.
        baseline = Counter(values).most_common(1)[0][0]
        baselines.append(baseline)
        exceptions.extend([i, j, value] for j, value in enumerate(values) if value != baseline)
    return {"report": report, "entry_denominator_bits": bits,
            "storage": "half_row_baselines_v1", "baselines": baselines, "exceptions": exceptions}


def unpack_candidate(data):
    size = data["report"]["size"]
    bits = data["entry_denominator_bits"]
    if type(size) is not int or size < 32 or size % 4 or type(bits) is not int or not 1 <= bits <= 100:
        raise ValueError("Invalid polygon size or denominator.")
    if data["storage"] != "half_row_baselines_v1":
        raise ValueError("Unknown format.")
    baselines = data["baselines"]
    if len(baselines) != size // 2 - 1 or any(type(v) is not int or v < 0 for v in baselines):
        raise ValueError("Invalid half-block rows.")
    rows = [[v] * (size // 2 + 1) for v in baselines]
    seen = set()
    for entry in data["exceptions"]:
        if len(entry) != 3 or any(type(v) is not int for v in entry):
            raise ValueError("Invalid exception.")
        i, j, value = entry
        if not 0 <= i < len(rows) or not 0 <= j < len(rows[0]) or value < 0 or (i, j) in seen:
            raise ValueError("Out-of-range, negative or repeated exception.")
        rows[i][j] = value
        seen.add((i, j))
    return rows


def read_candidate(case):
    data = json.loads(candidate_path(case).read_text(encoding="utf-8"))
    actual = tuple(data["report"][key] for key in ("size", "contrast", "source_radius"))
    if actual != CASES[case] or data["entry_denominator_bits"] != 80:
        raise ValueError("Unexpected archived parameters.")
    return data


@lru_cache(maxsize=4)
def half_interval_geometry(size, eta, radius):
    geo = certified_geometry(size, eta, radius)
    beta, weights = geo["beta"], geo["weights"]
    cosine = cos_interval(pi_interval() / size)
    attenuation = geo["alpha"] * factor_interval(eta) / cosine
    epsilon = 1 - attenuation
    pole = size // 4
    source_indices = [i % size for i in range(-pole + 1, pole)]
    target_indices = [j % size for j in range(-pole, pole + 1)]
    wi = [weights[i] for i in source_indices]
    total = (1 - weights[-1] - 2 * weights[pole]) / 2
    pole_mass = I.rational(1, 2 * size) - weights[pole] / 2 - weights[-1] / 4
    targets = [geo["targets"][j] for j in target_indices]
    # The pole x coordinates and the mean y coordinate vanish identically.
    targets[0], targets[-1] = (I.exact(0), -beta), (I.exact(0), beta)
    tj = [pole_mass] + [I.rational(1, size)] * (size // 2 - 1) + [pole_mass]
    probabilities = [mass / total for mass in tj]
    sources = [(geo["sources"][i][0] / attenuation, geo["sources"][i][1] / attenuation)
               for i in source_indices]
    mean_x = cosine / (pi_interval() * total)
    variance_x = beta**2 / (4 * total) - mean_x**2
    variance_y = beta**2 * (1 - 1 / (4 * total))
    centered = [(point[0] - mean_x, point[1]) for point in targets]
    return {"size": size, "weights": wi, "masses": tj, "total": total,
            "probabilities": probabilities, "sources": sources, "targets": targets,
            "mean_x": mean_x, "variance_x": variance_x, "variance_y": variance_y,
            "centered": centered, "epsilon": epsilon,
            "axis_radius_ratio": I.exact(radius) / (attenuation * beta),
            "pole_mass": pole_mass, "containment": geo["polygon_containment_margin"],
            "center_mass": weights[-1]}


def integer_dot(row, targets, coordinate, denominator):
    lo = sum(p * point[coordinate].lo for p, point in zip(row, targets))
    hi = sum(p * point[coordinate].hi for p, point in zip(row, targets))
    return I(lo // denominator, ceil_div(hi, denominator))


def near_axis_capacity_certificate(size, contrast, radius):
    geo = half_interval_geometry(size, Fraction(contrast), Fraction(radius))
    gap = geo["pole_mass"] - geo["weights"][-1] * (1 - geo["axis_radius_ratio"])
    elementary = (geo["pole_mass"].lo > 0 and 0 <= geo["epsilon"].lo <= geo["epsilon"].hi < SCALE
                  and 0 <= geo["axis_radius_ratio"].lo <= geo["axis_radius_ratio"].hi <= SCALE
                  and geo["containment"].lo > 0 and geo["center_mass"].lo >= 0
                  and min(w.lo for w in geo["weights"]) > 0)
    return {"size": size, "contrast": str(contrast), "source_radius": str(radius),
            "elementary_conditions_certified": elementary,
            "pole_capacity_minus_nearest_source_requirement_interval": gap.floats(),
            "gap_upper_bound_exact": str(gap.hi) + "/" + str(SCALE),
            "equality_face_excluded": elementary and gap.hi < 0,
            "scope": "Necessary condition r_N >= w_(N/4-1)*(1-t). A violation excludes saturation for this auxiliary grid, not all continuous kernels."}


def verify_half_candidate(data, target_contrast=None):
    rows = unpack_candidate(data)
    report = data["report"]
    size = report["size"]
    eta, radius = Fraction(str(report["contrast"])), Fraction(str(report["source_radius"]))
    if target_contrast is not None:
        eta = Fraction(target_contrast)
    if not 0 <= eta < 1 or not 0 < radius < 1:
        raise ValueError("Invalid strength or radius.")
    geo = half_interval_geometry(size, eta, radius)
    return _verify_half_geometry(data, geo, str(eta), str(radius), rows=rows)


def _verify_half_geometry(data, geo, contrast_label, radius_label, rows=None):
    """Certify positivity; callers must justify the exact geometry identities of note 40."""
    rows = unpack_candidate(data) if rows is None else rows
    size = data["report"]["size"]
    structural = (geo["total"].lo > 0 and geo["pole_mass"].lo > 0
                  and geo["containment"].lo > 0 and geo["center_mass"].lo >= 0
                  and geo["epsilon"].lo >= 0 and geo["epsilon"].hi < SCALE
                  and geo["axis_radius_ratio"].hi <= SCALE
                  and geo["variance_x"].lo > 0 and geo["variance_y"].lo > 0
                  and min(w.lo for w in geo["weights"]) > 0)
    if not structural:
        return {"verified": False, "reason": "Required exact half geometry is not certified positive."}
    denominator = 1 << data["entry_denominator_bits"]
    max_vx = I(max(v[0].abs_upper().hi for v in geo["centered"]))
    max_vy = I(max(v[1].abs_upper().hi for v in geo["centered"]))
    row_bounds, mass_errors, moment_errors = [], [], []
    for row, wi, point in zip(rows, geo["weights"], geo["sources"]):
        a = wi - I.rational(sum(row), denominator)
        bx = wi * point[0] - integer_dot(row, geo["targets"], 0, denominator)
        by = wi * point[1] - integer_dot(row, geo["targets"], 1, denominator)
        cx = bx - a * geo["mean_x"]
        bound = (a.abs_upper() + cx.abs_upper() * max_vx / I(geo["variance_x"].lo)
                 + by.abs_upper() * max_vy / I(geo["variance_y"].lo))
        row_bounds.append(I(bound.hi))
        mass_errors.append(a.abs_upper().hi)
        moment_errors.append((bx.abs_upper() + by.abs_upper()).hi)
    col_errors = [(mass - I.rational(sum(row[j] for row in rows), denominator)).abs_upper()
                  for j, mass in enumerate(geo["masses"])]
    max_pi = I(max(p.hi for p in geo["probabilities"]))
    sum_row_bounds = sum(row_bounds)
    max_col_repair = I(max(e.hi for e in col_errors)) + max_pi * sum_row_bounds
    max_row_fraction = I(max(w.hi for w in geo["weights"])) / I(geo["total"].lo)
    correction = max_pi * I(max(e.hi for e in row_bounds)) + max_row_fraction * max_col_repair
    minimum = I.rational(min(min(row) for row in rows), denominator)
    direct_margin = minimum - correction
    # Divide by the known target probability to avoid treating tiny pole columns
    # as though they had the same correction scale as the interior columns.
    scaled_col_error = I(max((e / I(p.lo)).hi for e, p in zip(col_errors, geo["probabilities"])))
    scaled_col_repair = scaled_col_error + sum_row_bounds
    normalized_margins = []
    for row, wi, bound in zip(rows, geo["weights"], row_bounds):
        # All interior columns have the same target probability.
        normalized_minimum = min(
            value * SCALE**2 // (denominator * probability.hi)
            for value, probability in ((row[0], geo["probabilities"][0]),
                                       (min(row[1:-1]), geo["probabilities"][1]),
                                       (row[-1], geo["probabilities"][-1])))
        normalized_margins.append(I(normalized_minimum) - bound - wi / geo["total"] * scaled_col_repair)
    normalized_lower = min(v.lo for v in normalized_margins)
    probability_endpoint = (min(p.lo for p in geo["probabilities"]) if normalized_lower >= 0
                            else max(p.hi for p in geo["probabilities"]))
    weighted_margin = I(normalized_lower) * I(probability_endpoint)
    margin = I(max(direct_margin.lo, weighted_margin.lo))
    return {"verified": margin.lo > 0, "size": size, "contrast": contrast_label, "source_radius": radius_label,
            "uniform_reserve_cap_interval": geo["epsilon"].floats(),
            "pole_mass_each_interval": geo["pole_mass"].floats(),
            "axis_radius_ratio_interval": geo["axis_radius_ratio"].floats(),
            "half_target_variance_x_interval": geo["variance_x"].floats(),
            "half_target_variance_y_interval": geo["variance_y"].floats(),
            "maximum_row_mass_error_bound": max(mass_errors) / SCALE,
            "maximum_row_moment_l1_error_bound": max(moment_errors) / SCALE,
            "maximum_original_column_error_bound": max(e.hi for e in col_errors) / SCALE,
            "maximum_half_entry_correction_bound": correction.hi / SCALE,
            "minimum_original_half_entry": minimum.lo / SCALE,
            "repaired_half_entry_lower_bound": margin.lo / SCALE,
            "repaired_half_entry_lower_bound_exact": str(margin.lo) + "/" + str(SCALE),
            "direct_repaired_half_entry_lower_bound": direct_margin.lo / SCALE,
            "probability_scaled_repaired_half_entry_lower_bound": weighted_margin.lo / SCALE,
            "scope": "Correction acts only inside a positive half block. Reflection and the exact pole rows preserve every forbidden zero in the core; the analytic full uniform reserve equals its proven projection upper bound."}


@lru_cache(maxsize=3)
def repaired_half_float(case, target_contrast=None):
    data = read_candidate(case)
    size, contrast, radius = CASES[case]
    if target_contrast is not None:
        contrast = float(target_contrast)
    geo = half_geometry(size, contrast, radius)
    matrix = np.array(unpack_candidate(data), dtype=float) / (1 << data["entry_denominator_bits"])
    weights, masses, targets = geo["source_mass"], geo["target_mass"], geo["targets"]
    total = geo["total_mass"]
    probabilities = masses / total
    mean = probabilities @ targets
    centered = targets - mean
    covariance = (centered.T * probabilities) @ centered
    a = weights - matrix.sum(axis=1)
    b = weights[:, None] * geo["sources"] - matrix @ targets
    c = b - a[:, None] * mean
    matrix += probabilities * (a[:, None] + c @ np.linalg.solve(covariance, centered.T))
    d = masses - matrix.sum(axis=0)
    matrix += weights[:, None] / total * d
    return geo, matrix


class HalfFaceCertificateTests(unittest.TestCase):
    def test_archived_half_faces_have_strict_certificates(self):
        for case in CASES:
            self.assertTrue(verify_half_candidate(read_candidate(case))["verified"])

    def test_corrupted_matrix_is_rejected(self):
        data = read_candidate("512")
        data["exceptions"][0][2] += 1 << 65
        self.assertFalse(verify_half_candidate(data)["verified"])

    def test_retargeted_strength_is_certified_by_probability_scaled_bound(self):
        certificate = verify_half_candidate(read_candidate("2048"), RECERTIFIED_CONTRAST)
        self.assertTrue(certificate["verified"])
        self.assertEqual(certificate["contrast"], str(RECERTIFIED_CONTRAST))
        self.assertLess(certificate["direct_repaired_half_entry_lower_bound"], 0)

    def test_positive_pole_capacity_can_still_fail_near_axis_requirement(self):
        certificate = near_axis_capacity_certificate(64, "0.135", "0.991")
        self.assertTrue(certificate["elementary_conditions_certified"])
        self.assertTrue(certificate["equality_face_excluded"])
        for n, eta, radius in CASES.values():
            certificate = near_axis_capacity_certificate(n, str(eta), str(radius))
            self.assertGreater(certificate["pole_capacity_minus_nearest_source_requirement_interval"][0], 0)

    def test_analytic_covariance_matches_direct_weighted_sum(self):
        for case, (size, eta, radius) in CASES.items():
            g = half_geometry(size, eta, radius)
            p = g["target_mass"] / g["total_mass"]
            mean = p @ g["targets"]
            v = g["targets"] - mean
            covariance = (v.T * p) @ v
            certified = half_interval_geometry(size, Fraction(str(eta)), Fraction(str(radius)))
            expected = np.diag([sum(certified["variance_x"].floats()) / 2,
                                sum(certified["variance_y"].floats()) / 2])
            np.testing.assert_allclose(covariance, expected, atol=5e-15)

    def test_repaired_half_rows_columns_and_vector_means(self):
        for case, target in verification_cases():
            g, h = repaired_half_float(case, target)
            self.assertGreater(float(h.min()), 0)
            np.testing.assert_allclose(h.sum(axis=1), g["source_mass"], atol=3e-16)
            np.testing.assert_allclose(h.sum(axis=0), g["target_mass"], atol=3e-16)
            np.testing.assert_allclose(h @ g["targets"], g["source_mass"][:, None] * g["sources"], atol=3e-16)

    def test_full_reconstruction_keeps_exact_forbidden_support(self):
        for case, target in verification_cases():
            g, half = repaired_half_float(case, target)
            size = g["size"]
            # Setting epsilon to zero exposes the core without subtractive cancellation.
            core_geo = dict(g, epsilon=0)
            core = expand_half_coupling(core_geo, half)
            cosine = np.cos(np.arange(size) * 2 * math.pi / size)
            crossing = cosine[:, None] * cosine < -1e-12
            self.assertEqual(np.count_nonzero(core[:-1][crossing]), 0)
            for row in (size // 4, 3 * size // 4, size):
                self.assertEqual(np.flatnonzero(core[row]).tolist(), [size // 4, 3 * size // 4])

    def test_full_matrix_and_projection_reserve_agree(self):
        from convex_order_interface import absolute_factor
        for case, target in verification_cases():
            g, half = repaired_half_float(case, target)
            size = g["size"]
            j = expand_half_coupling(g, half)
            angles = np.arange(size) * 2 * math.pi / size
            unit = np.stack((np.cos(angles), np.sin(angles)), axis=1)
            source_points = np.vstack((g["radius"] * unit, np.zeros(2)))
            np.testing.assert_allclose(j.sum(axis=1), g["weights"], atol=4e-16)
            np.testing.assert_allclose(j.sum(axis=0), 1 / size, atol=4e-16)
            np.testing.assert_allclose(j @ (g["beta"] * unit), g["weights"][:, None] * source_points, atol=4e-16)
            from outcome_updates import WINDOW_ATTENUATION as alpha
            expected = 1 - alpha * absolute_factor(g["contrast"]) / math.cos(math.pi / size)
            self.assertAlmostEqual(g["epsilon"], expected, places=14)

    def test_full_selective_density_after_half_face_repair(self):
        from polygon_martingale_certificate import quantize_source
        from compensated_angle_kernel import warped_angle
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        from outcome_updates import WINDOW_ATTENUATION as alpha
        outgoing = np.linspace(-1, 7, 41)
        setting = .37
        for case, target in verification_cases():
            g, half = repaired_half_float(case, target)
            size, eta = g["size"], g["contrast"]
            j = expand_half_coupling(g, half)
            conditional = j / g["weights"][:, None]
            centers = np.arange(size) * 2 * math.pi / size + setting
            for mark in (0, 1):
                relative = outgoing - setting
                left, lower, upper, center = quantize_source(warped_angle(relative, mark, eta), size, g["radius"])
                probabilities = (lower[:, None] * conditional[left]
                                 + upper[:, None] * conditional[(left + 1) % size]
                                 + center[:, None] * conditional[-1])
                f = (1 + (2 * mark - 1) * eta * np.cos(relative)) / 2
                for z in ([1, alpha / math.sqrt(2), alpha / math.sqrt(2)], [.4, .1, -.2]):
                    arc_mass = (z[0] + g["beta"] / alpha * (z[1] * np.cos(centers) + z[2] * np.sin(centers))) / size
                    observed = f * size / (2 * math.pi) * (probabilities @ arc_mass)
                    expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ z, outgoing)
                    np.testing.assert_allclose(observed, expected, atol=2e-13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve-case", choices=tuple(CASES))
    parser.add_argument("--candidate-output", type=Path)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    if args.solve_case:
        if args.candidate_output is None or args.candidate_output.exists():
            parser.error("Choose a new --candidate-output path.")
        report, half = solve_half_polygon(*CASES[args.solve_case], return_half=True)
        print(json.dumps(report), flush=True)
        if half is None:
            raise SystemExit(1)
        data = pack_candidate(report, half)
        args.candidate_output.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")
        certificate = verify_half_candidate(data)
        print(json.dumps(certificate), flush=True)
        if not certificate["verified"]:
            raise SystemExit(1)
        return
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(HalfFaceCertificateTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"certificates": {case: dict(verify_half_candidate(read_candidate(case)),
                                         candidate_sha256=hashlib.sha256(candidate_path(case).read_bytes()).hexdigest(),
                                         candidate_bytes=candidate_path(case).stat().st_size)
                                for case in CASES},
              "recertified_strength": dict(verify_half_candidate(read_candidate("2048"), RECERTIFIED_CONTRAST),
                                            candidate_sha256=hashlib.sha256(candidate_path("2048").read_bytes()).hexdigest(),
                                            candidate_generated_at_contrast="0.14473"),
              "near_axis_capacity_counterexample": near_axis_capacity_certificate(64, "0.135", "0.991"),
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("half_face_certificate_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
