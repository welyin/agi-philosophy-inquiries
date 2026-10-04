"""New polygon candidates, independent certification, and full-instrument checks."""

import argparse
import hashlib
import json
import math
import platform
import unittest
from functools import lru_cache
from pathlib import Path

import numpy as np

import polygon_martingale_certificate as polygon
from outcome_updates import WINDOW_ATTENUATION as ALPHA


CASES = {
    "512": (512, 0.1445, 0.989635, "polygon_coupling_512_01445.json"),
    "1024": (1024, 0.1447, 0.9896206, "polygon_coupling_1024_01447.json"),
}


def candidate_path(case):
    return Path(__file__).with_name(CASES[case][3])


def read_candidate(case):
    data = json.loads(candidate_path(case).read_text(encoding="utf-8"))
    expected = CASES[case][:3]
    actual = tuple(data["report"][key] for key in ("size", "contrast", "source_radius"))
    if actual != expected or data["entry_denominator_bits"] != 60:
        raise ValueError("Saved candidate parameters do not match this research case.")
    return data


@lru_cache(maxsize=2)
def numerical_coupling(case):
    return polygon.repaired_coupling_float(read_candidate(case))


def reverse_probabilities(case, angles, selected_columns=None):
    size, _, radius, _ = CASES[case]
    matrix, weights, _, _ = numerical_coupling(case)
    if selected_columns is not None:
        matrix = matrix[:, selected_columns]
    conditional = matrix / weights[:, None]
    left, lower, upper, center = polygon.quantize_source(angles, size, radius)
    return (lower[..., None] * conditional[left]
            + upper[..., None] * conditional[(left + 1) % size]
            + center[..., None] * conditional[-1])


def selective_density(case, summary, outgoing, mark, setting=0):
    from compensated_angle_kernel import warped_angle
    size, eta, _, _ = CASES[case]
    relative = np.asarray(outgoing) - setting
    probabilities = reverse_probabilities(case, warped_angle(relative, mark, eta))
    centers = np.arange(size) * 2 * math.pi / size + setting
    input_masses = (summary[0] + np.sinc(1 / size) / ALPHA
                    * (summary[1] * np.cos(centers) + summary[2] * np.sin(centers))) / size
    f = (1 + (2 * mark - 1) * eta * np.cos(relative)) / 2
    return f * size / (2 * math.pi) * (probabilities @ input_masses)


def certificate_report(case):
    certificate = polygon.verify_polygon_data(read_candidate(case))
    certificate["candidate_sha256"] = hashlib.sha256(candidate_path(case).read_bytes()).hexdigest()
    return certificate


class FinePolygonTests(unittest.TestCase):
    def test_both_new_candidates_have_strict_integer_certificates(self):
        for case in CASES:
            with self.subTest(case=case):
                self.assertTrue(certificate_report(case)["verified"])

    def test_invalid_parameters_are_rejected_before_geometry(self):
        data = read_candidate("512")
        for key, value in (("size", 513), ("contrast", 1), ("source_radius", 0)):
            original = data["report"][key]
            data["report"][key] = value
            with self.assertRaises(ValueError):
                polygon.verify_polygon_data(data)
            data["report"][key] = original

    def test_altered_matrix_is_not_a_certificate(self):
        data = read_candidate("1024")
        data["coupling_numerators"][0][0] += 1 << 45
        self.assertFalse(polygon.verify_polygon_data(data)["verified"])

    def test_stable_weights_match_independent_piecewise_quadrature(self):
        from compensated_angle_kernel import preconvolution_density
        nodes, gauss = np.polynomial.legendre.leggauss(32)
        for size, eta, radius, _ in CASES.values():
            step = 2 * math.pi / size
            angle = (np.arange(size)[:, None] + (nodes + 1) / 2) * step
            left, lower, upper, center = polygon.quantize_source(angle, size, radius)
            mass = preconvolution_density(angle, eta) * gauss * step / (4 * math.pi)
            observed = np.zeros(size + 1)
            np.add.at(observed, left.ravel(), (mass * lower).ravel())
            np.add.at(observed, ((left + 1) % size).ravel(), (mass * upper).ravel())
            observed[-1] = np.sum(mass * center)
            np.testing.assert_allclose(observed, polygon.source_weights(size, eta, radius), atol=8e-16)

    def test_repaired_rows_columns_and_barycenters(self):
        for case in CASES:
            matrix, weights, sources, targets = numerical_coupling(case)
            self.assertGreater(float(matrix.min()), 0)
            np.testing.assert_allclose(matrix.sum(axis=1), weights, atol=2e-17)
            np.testing.assert_allclose(matrix.sum(axis=0), 1 / matrix.shape[1], atol=2e-17)
            np.testing.assert_allclose(matrix @ targets, weights[:, None] * sources, atol=2e-17)

    def test_reverse_distribution_has_exact_required_first_moment(self):
        angles = np.linspace(-3, 9, 137)
        expected = ALPHA * np.stack((np.cos(angles), np.sin(angles)), axis=1)
        for case in CASES:
            probabilities = reverse_probabilities(case, angles)
            targets = numerical_coupling(case)[3]
            self.assertGreater(float(probabilities.min()), 0)
            np.testing.assert_allclose(probabilities.sum(axis=1), 1, atol=5e-13)
            np.testing.assert_allclose(probabilities @ targets, expected, atol=5e-13)

    def test_entire_selective_output_matches_original_fresh_instrument(self):
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        angles = np.linspace(-2, 8, 53)
        for case, (_, eta, _, _) in CASES.items():
            for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
                for setting in (0, 0.71):
                    for mark in (0, 1):
                        expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, angles)
                        np.testing.assert_allclose(selective_density(case, summary, angles, mark, setting),
                                                   expected, atol=2e-13)

    def test_forward_normalization_by_output_quadrature(self):
        from atomic_reset_kernel import atom_target
        nodes, gauss = np.polynomial.legendre.leggauss(12)
        for case, (size, eta, _, _) in CASES.items():
            step = 2 * math.pi / size
            psi = ((np.arange(size)[:, None] + (nodes + 1) / 2) * step).ravel()
            weights = np.tile(gauss * step / 2, size)
            probabilities = reverse_probabilities(case, psi, [0, size // 4, size // 3])
            total = np.zeros(3)
            q = math.sqrt(1 - eta**2)
            for mark in (0, 1):
                sign = 2 * mark - 1
                y = atom_target(psi, mark, eta)
                factor = (1 + sign * eta * np.cos(y)) / 2 * q / (1 - sign * eta * np.cos(psi))
                total += ((weights * factor)[:, None] * probabilities).sum(axis=0) * size / (2 * math.pi)
            np.testing.assert_allclose(total, 1, atol=2e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    parser.add_argument("--solve-case", choices=tuple(CASES))
    parser.add_argument("--candidate-output", type=Path)
    parser.add_argument("--symmetric", action="store_true", help="Use the optional reflection-reduced IPM generator.")
    args = parser.parse_args()
    if args.solve_case:
        if args.candidate_output is None:
            parser.error("--solve-case requires a new --candidate-output path.")
        if args.candidate_output.exists():
            parser.error("Choose a new output path; archived candidates are not overwritten.")
        solver = polygon.solve_polygon
        if args.symmetric:
            from symmetric_polygon_solver import solve_symmetric_polygon
            solver = solve_symmetric_polygon
        report, matrix = solver(*CASES[args.solve_case][:3])
        print(json.dumps(report, indent=2), flush=True)
        if matrix is None:
            raise SystemExit(1)
        data = {"report": report, "entry_denominator_bits": 60,
                "coupling_numerators": np.floor(matrix * (1 << 60)).astype(np.int64).tolist()}
        args.candidate_output.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")
        certificate = polygon.verify_polygon_data(data)
        print(json.dumps(certificate, indent=2))
        if not certificate["verified"]:
            raise SystemExit(1)
        return
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FinePolygonTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"certificates": {case: certificate_report(case) for case in CASES},
              "environment": {"python": platform.python_version(), "numpy": np.__version__},
              "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
              "scope": "Exact analytic repair is strictly positive for both archived candidates. Barycentric source quantization and uniform target arcs give full continuous instrument kernels. Round 34 extends eta=0.1447 to all smaller strengths; endpoint attainability remains unproved."}
    if args.write_results:
        Path(__file__).with_name("fine_polygon_certificate_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
