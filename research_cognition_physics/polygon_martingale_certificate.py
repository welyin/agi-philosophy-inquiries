"""Finite coupling with analytic source quantization and continuous target lifting."""

import argparse
import hashlib
import json
import math
import platform
import sys
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA
from certified_intervals import Interval as I, SCALE, ceil_div, pi_interval, sin_interval, cos_interval


CANDIDATE_PATH = Path(__file__).with_name("polygon_coupling_candidate.json")


def source_weights(size, contrast, radius, terms=12):
    step = 2 * math.pi / size
    angles = np.arange(size) * step
    q = math.sqrt(1 - contrast**2)
    r = (contrast / (1 + q))**2
    # Product identities avoid subtracting nearly equal cosines on fine grids.
    integral = np.full(size, 4 * math.sin(step / 2)**2)
    for k in range(1, terms + 1):
        ck = r**k * (1 + 2 * k * q)
        ik = (4 * math.sin((2 * k + 1) * step / 2)
              * math.sin((2 * k - 1) * step / 2) / (4 * k**2 - 1))
        integral += 2 * ck * ik * np.cos(2 * k * angles)
    weights = ALPHA * integral / (2 * math.pi * radius * math.sin(step))
    return np.append(weights, 1 - weights.sum())


def solve_polygon(size=256, contrast=0.144, radius=0.9897):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / ".research_runtime"))
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    angles = np.arange(size) * 2 * math.pi / size
    unit = np.stack((np.cos(angles), np.sin(angles)), axis=-1)
    beta = np.sinc(1 / size)
    source = source_weights(size, contrast, radius)
    if min(source) < 0 or radius * math.cos(math.pi / size) < ALPHA:
        raise ValueError("The source polygon must contain the entire inner circle.")
    points = np.vstack((radius * unit, np.zeros(2)))
    targets = beta * unit
    count = size + 1
    ri, ci = np.indices((count, size))
    variables = np.arange(count * size)
    rows = np.concatenate((ri.ravel(), count + ci.ravel(), count + size + ri.ravel(),
                           2 * count + size + ri.ravel(), np.arange(count), count + np.arange(size)))
    columns = np.concatenate((variables, variables, variables, variables,
                              np.full(count + size, count * size)))
    values = np.concatenate((np.ones(count * size), np.ones(count * size),
                             np.tile(targets[:, 0], count), np.tile(targets[:, 1], count),
                             source, np.full(size, 1 / size)))
    matrix = coo_matrix((values, (rows, columns)), shape=(3 * count + size, count * size + 1)).tocsr()
    right = np.concatenate((source, np.full(size, 1 / size), source * points[:, 0], source * points[:, 1]))
    objective = np.zeros(count * size + 1)
    objective[-1] = -1
    result = linprog(objective, A_eq=matrix, b_eq=right,
                     bounds=[(0, None)] * (count * size) + [(0, 1)], method="highs",
                     options={"dual_feasibility_tolerance": 1e-9, "primal_feasibility_tolerance": 1e-9})
    report = {"size": size, "contrast": contrast, "source_radius": radius,
              "target_arc_mean_radius": float(beta), "source_center_mass": float(source[-1]),
              "status": int(result.status), "message": result.message}
    if not result.success:
        return report, None
    reserve = float(result.x[-1])
    coupling = result.x[:-1].reshape(count, size) + reserve * source[:, None] / size
    report.update({"uniform_reserve": reserve, "minimum_entry": float(coupling.min()),
                   "maximum_equality_residual": float(np.max(np.abs(matrix @ result.x - right)))})
    return report, coupling


@lru_cache(maxsize=4)
def certified_geometry(size, contrast, radius, terms=12):
    """All arithmetic here is outward-rounded integer interval arithmetic."""
    pi = pi_interval()
    alpha = 4 * sin_interval(I.rational(1, 4))
    eta, rho = I.exact(contrast), I.exact(radius)
    step, half = 2 * pi / size, pi / size
    beta = sin_interval(half) / half
    cosine, sine = [], []
    for j in range(size):
        index = j if j <= size // 2 else j - size
        angle = 2 * pi * index / size
        cosine.append(cos_interval(angle))
        sine.append(sin_interval(angle))
    q = (1 - eta**2).sqrt()
    r = (eta / (1 + q))**2
    integral0 = 2 * (1 - cosine[1])
    prefactor = alpha / (2 * pi * rho * sin_interval(step))
    tail = r**(terms + 1) * ((1 + 2 * (terms + 1) * q) / (1 - r) + 2 * q * r / (1 - r)**2)
    tail_weight = (2 * prefactor * integral0 * tail).abs_upper()
    coefficients = []
    for k in range(1, terms + 1):
        ck = r**k * (1 + 2 * k * q)
        ik = 2 * (cosine[(2 * k) % size] - cosine[1]) / (1 - 4 * k**2)
        coefficients.append(2 * ck * ik)
    weights = []
    for j in range(size):
        integral = integral0
        for k, coefficient in enumerate(coefficients, 1):
            integral += coefficient * cosine[(2 * k * j) % size]
        weights.append((prefactor * integral).widen(tail_weight))
    weights.append(1 - sum(weights))
    return {"alpha": alpha, "beta": beta, "weights": weights,
            "targets": [(beta * x, beta * y) for x, y in zip(cosine, sine)],
            "sources": [(rho * x, rho * y) for x, y in zip(cosine, sine)] + [(I.exact(0), I.exact(0))],
            "polygon_containment_margin": rho * cos_interval(half) - alpha,
            "fourier_tail_weight_bound": tail_weight}


def read_candidate():
    return json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))


def verify_candidate(data=None):
    if data is None:
        data = read_candidate()
    report = data["report"]
    size = int(report["size"])
    eta, radius = Fraction(str(report["contrast"])), Fraction(str(report["source_radius"]))
    bits = data["entry_denominator_bits"]
    if size != 256 or eta != Fraction(18, 125) or radius != Fraction(9897, 10000) or bits != 60:
        raise ValueError("This saved certificate has fixed N=256, eta=0.144, rho=0.9897 and 60-bit entries.")
    return verify_polygon_data(data)


def verify_polygon_data(data):
    """Verify a dyadic candidate for an even polygon, independently of its solver."""
    report = data["report"]
    size, bits = report["size"], data["entry_denominator_bits"]
    eta = Fraction(str(report["contrast"]))
    radius = Fraction(str(report["source_radius"]))
    if (type(size) is not int or size < 32 or size % 4 != 0
            or type(bits) is not int or not 1 <= bits <= 100
            or not 0 <= eta < 1 or not 0 < radius < 1):
        raise ValueError("Invalid certified polygon parameters.")
    numerators = data["coupling_numerators"]
    if len(numerators) != size + 1 or any(len(row) != size for row in numerators):
        raise ValueError("Invalid matrix shape.")
    if any(not isinstance(p, int) or p < 0 for row in numerators for p in row):
        raise ValueError("Entries must be nonnegative integers.")
    geometry = certified_geometry(size, eta, radius)
    weights, sources, targets = geometry["weights"], geometry["sources"], geometry["targets"]
    denominator = 1 << bits
    correction_bounds = []
    row_mass_errors, row_mean_errors = [], []
    for row, weight, source_point in zip(numerators, weights, sources):
        mass_error = weight - I.rational(sum(row), denominator)
        mean_errors = []
        for coordinate in (0, 1):
            lower = sum(p * target[coordinate].lo for p, target in zip(row, targets))
            upper = sum(p * target[coordinate].hi for p, target in zip(row, targets))
            moment = I(lower // denominator, ceil_div(upper, denominator))
            mean_errors.append(weight * source_point[coordinate] - moment)
        moment_bound = mean_errors[0].abs_upper() + mean_errors[1].abs_upper()
        bound = (mass_error.abs_upper() + 2 * moment_bound / I(geometry["beta"].lo)) / size
        correction_bounds.append(I(bound.hi))
        row_mass_errors.append(mass_error.abs_upper().hi)
        row_mean_errors.append(moment_bound.hi)
    column_errors = [
        (I.rational(1, size) - I.rational(sum(row[j] for row in numerators), denominator)).abs_upper()
        for j in range(size)]
    largest_column_error = I(max(e.hi for e in column_errors))
    column_repair_bound = largest_column_error + sum(correction_bounds)
    largest_weight = I(max(w.hi for w in weights))
    largest_correction = I(max(b.hi for b in correction_bounds)) + largest_weight * column_repair_bound
    minimum_entry = I.rational(min(min(row) for row in numerators), denominator)
    positivity_margin = minimum_entry - largest_correction
    success = (geometry["polygon_containment_margin"].lo > 0
               and min(w.lo for w in weights) > 0 and positivity_margin.lo > 0)
    return {"verified": success, "size": size, "contrast": str(eta), "source_radius": str(radius),
            "interval_denominator_bits": 100,
            "polygon_containment_margin_interval": geometry["polygon_containment_margin"].floats(),
            "source_center_mass_interval": weights[-1].floats(),
            "fourier_tail_weight_bound": geometry["fourier_tail_weight_bound"].hi / SCALE,
            "maximum_row_mass_error_bound": max(row_mass_errors) / SCALE,
            "maximum_row_mean_l1_error_bound": max(row_mean_errors) / SCALE,
            "maximum_original_column_error_bound": largest_column_error.hi / SCALE,
            "maximum_entry_correction_bound": largest_correction.hi / SCALE,
            "minimum_original_entry": minimum_entry.lo / SCALE,
            "repaired_entry_lower_bound": positivity_margin.lo / SCALE,
            "repaired_entry_lower_bound_exact": f"{positivity_margin.lo}/{SCALE}"}


def repaired_coupling_float(data=None):
    """Numerical view of the exact algebraic correction, for independent QA."""
    data = read_candidate() if data is None else data
    size = data["report"]["size"]
    eta, radius = data["report"]["contrast"], data["report"]["source_radius"]
    matrix = np.asarray(data["coupling_numerators"], dtype=float) / (1 << data["entry_denominator_bits"])
    angles = np.arange(size) * 2 * math.pi / size
    unit = np.stack((np.cos(angles), np.sin(angles)), axis=1)
    beta = np.sinc(1 / size)
    targets = beta * unit
    points = np.vstack((radius * unit, np.zeros(2)))
    source = source_weights(size, eta, radius)
    mass_error = source - matrix.sum(axis=1)
    mean_error = source[:, None] * points - matrix @ targets
    row_repair = mass_error[:, None] / size + 2 * mean_error @ targets.T / (size * beta**2)
    updated = matrix + row_repair
    columns = 1 / size - updated.sum(axis=0)
    updated += source[:, None] * columns
    return updated, source, points, targets


def quantize_source(angles, size=256, radius=0.9897):
    angles = np.mod(np.asarray(angles), 2 * math.pi)
    step = 2 * math.pi / size
    left = np.floor(angles / step).astype(int) % size
    offset = angles - left * step
    lower = ALPHA / radius * np.sin(step - offset) / math.sin(step)
    upper = ALPHA / radius * np.sin(offset) / math.sin(step)
    return left, lower, upper, 1 - lower - upper


def reverse_arc_probabilities(angles, selected_columns=None):
    matrix, source, _, _ = repaired_coupling_float()
    if selected_columns is not None:
        matrix = matrix[:, selected_columns]
    conditional = matrix / source[:, None]
    left, lower, upper, center = quantize_source(angles)
    return (lower[..., None] * conditional[left] + upper[..., None] * conditional[(left + 1) % 256]
            + center[..., None] * conditional[-1])


def pushed_selective_density(summary, outgoing, outcome, setting=0.0):
    """Integrate the input density exactly over each target arc."""
    from compensated_angle_kernel import warped_angle
    eta, size = 0.144, 256
    relative = np.asarray(outgoing) - setting
    psi = warped_angle(relative, outcome, eta)
    probabilities = reverse_arc_probabilities(psi)
    centers = np.arange(size) * 2 * math.pi / size + setting
    beta = np.sinc(1 / size)
    input_arc_mass = (summary[0] + beta / ALPHA * (summary[1] * np.cos(centers) + summary[2] * np.sin(centers))) / size
    f = (1 + (2 * outcome - 1) * eta * np.cos(relative)) / 2
    return f * size / (2 * math.pi) * (probabilities @ input_arc_mass)


class PolygonMartingaleCertificateTests(unittest.TestCase):
    def test_saved_dyadic_candidate_has_a_strict_integer_certificate(self):
        self.assertTrue(verify_candidate()["verified"])

    def test_certificate_rejects_a_substantially_corrupted_matrix(self):
        data = read_candidate()
        data["coupling_numerators"][0][0] += 1 << 45
        self.assertFalse(verify_candidate(data)["verified"])

    def test_quantization_has_positive_weights_and_exact_barycenters(self):
        size, radius = 256, 0.9897
        angles = np.linspace(-3, 9, 1007)
        left, lower, upper, center = quantize_source(angles)
        self.assertGreaterEqual(float(min(lower.min(), upper.min(), center.min())), -1e-12)
        actual = radius * (lower * np.exp(2j * math.pi * left / size)
                           + upper * np.exp(2j * math.pi * (left + 1) / size))
        np.testing.assert_allclose(actual, ALPHA * np.exp(1j * angles), atol=1e-13)

    def test_analytic_source_weights_match_independent_sector_quadrature(self):
        from compensated_angle_kernel import preconvolution_density
        size, radius, eta = 256, 0.9897, 0.144
        nodes, weights = np.polynomial.legendre.leggauss(48)
        step = 2 * math.pi / size
        expected = source_weights(size, eta, radius)
        for index in (0, 1, 17, 63, 127):
            integral = 0.0
            for sign in (-1, 1):
                offset = sign * (nodes + 1) * step / 2
                integrand = preconvolution_density(index * step + offset, eta) * np.sin(step - np.abs(offset))
                integral += float(weights @ integrand) * step / 2
            actual = ALPHA * integral / (2 * math.pi * radius * math.sin(step))
            self.assertAlmostEqual(actual, expected[index], places=14)

    def test_corrected_matrix_has_the_required_rows_columns_and_vector_means(self):
        matrix, source, points, targets = repaired_coupling_float()
        self.assertGreater(float(matrix.min()), 0)
        np.testing.assert_allclose(matrix.sum(axis=1), source, atol=1e-14)
        np.testing.assert_allclose(matrix.sum(axis=0), 1 / 256, atol=1e-14)
        np.testing.assert_allclose(matrix @ targets, source[:, None] * points, atol=1e-14)

    def test_arc_lift_and_composed_reverse_kernel_have_correct_conditional_means(self):
        matrix, source, _, targets = repaired_coupling_float()
        eta_angles = np.linspace(-2, 8, 211)
        left, lower, upper, center = quantize_source(eta_angles)
        conditional = matrix / source[:, None]
        probabilities = (lower[:, None] * conditional[left]
                         + upper[:, None] * conditional[(left + 1) % 256]
                         + center[:, None] * conditional[-1])
        np.testing.assert_allclose(probabilities.sum(axis=1), 1, atol=1e-12)
        expected = ALPHA * np.stack((np.cos(eta_angles), np.sin(eta_angles)), axis=1)
        np.testing.assert_allclose(probabilities @ targets, expected, atol=1e-12)

    def test_whole_selective_output_density_matches_original_instrument(self):
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        grid = np.linspace(-2, 8, 43)
        for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
            for setting in (0, 0.7):
                for mark in (0, 1):
                    expected = summary_density(fresh_branch_matrix(setting, mark, 0.144) @ summary, grid)
                    np.testing.assert_allclose(pushed_selective_density(summary, grid, mark, setting), expected, atol=1e-12)

    def test_forward_normalization_by_independent_piecewise_output_quadrature(self):
        from atomic_reset_kernel import atom_target
        nodes, weights = np.polynomial.legendre.leggauss(16)
        step = 2 * math.pi / 256
        psi = (np.arange(256)[:, None] + (nodes + 1) / 2) * step
        psi = psi.ravel()
        integration_weights = np.tile(weights * step / 2, 256)
        probabilities = reverse_arc_probabilities(psi, [0, 31, 99])
        total = np.zeros(3)
        eta, q = 0.144, math.sqrt(1 - 0.144**2)
        for mark in (0, 1):
            sign = 2 * mark - 1
            y = atom_target(psi, mark, eta)
            f = (1 + sign * eta * np.cos(y)) / 2
            dy_dpsi = q / (1 - sign * eta * np.cos(psi))
            total += ((integration_weights * f * dy_dpsi)[:, None] * probabilities).sum(axis=0) * 256 / (2 * math.pi)
        np.testing.assert_allclose(total, 1, atol=1e-11)

    def test_actual_mixture_keeps_identical_full_outputs(self):
        grid = np.linspace(0, 6, 31)
        first, second = np.array([1, ALPHA, 0]), np.array([1, 0, ALPHA])
        mixed = (pushed_selective_density(first, grid, 1) + pushed_selective_density(second, grid, 1)) / 2
        np.testing.assert_allclose(mixed, pushed_selective_density((first + second) / 2, grid, 1), atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=256)
    parser.add_argument("--contrast", type=float, default=0.144)
    parser.add_argument("--radius", type=float, default=0.9897)
    parser.add_argument("--write-results", action="store_true")
    parser.add_argument("--solve", action="store_true", help="Regenerate the numerical candidate using optional SciPy.")
    args = parser.parse_args()
    if not args.solve:
        checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PolygonMartingaleCertificateTests))
        if not checks.wasSuccessful():
            raise SystemExit(1)
        certificate = verify_candidate()
        certificate.update({"environment": {"python": platform.python_version(), "numpy": np.__version__},
                            "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
                            "candidate_sha256": hashlib.sha256(CANDIDATE_PATH.read_bytes()).hexdigest(),
                            "scope": "An independently checked integer interval certificate proves that the exact analytic correction of the saved dyadic matrix is positive. Source barycentric quantization and uniform-arc target lifting give an exact continuous martingale coupling at eta=0.144. Strength degradation extends it to every lower contrast; this is not merely finite-grid feasibility."})
        if args.write_results:
            Path(__file__).with_name("polygon_martingale_certificate_results.json").write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(certificate, indent=2))
        return
    report, coupling = solve_polygon(args.size, args.contrast, args.radius)
    print(json.dumps(report, indent=2), flush=True)
    if args.write_results and coupling is not None:
        # Dyadic entries are exact rationals; a separate certificate will repair residuals.
        bits = 60
        integers = np.floor(coupling * (1 << bits)).astype(np.int64)
        CANDIDATE_PATH.write_text(json.dumps({"report": report, "entry_denominator_bits": bits,
                                    "coupling_numerators": integers.tolist()}, separators=(",", ":")) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
