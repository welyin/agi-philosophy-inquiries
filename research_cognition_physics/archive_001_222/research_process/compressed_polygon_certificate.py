"""Lossless row-baseline storage for fine dyadic polygon certificates."""

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


CANDIDATE_PATH = Path(__file__).with_name("polygon_coupling_2048_014473_compressed.json")


def compress_candidate(report, coupling, bits=80):
    baselines, exceptions = [], []
    for i, row in enumerate(coupling):
        integers = [int(value * (1 << bits)) for value in row]
        baseline = min(integers)
        baselines.append(baseline)
        exceptions.extend([i, j, value] for j, value in enumerate(integers) if value != baseline)
    return {"report": report, "entry_denominator_bits": bits,
            "storage": "row_baseline_and_absolute_exceptions_v1",
            "row_baselines": baselines, "exceptions": exceptions}


def expand_candidate(data):
    size = data["report"]["size"]
    if data["storage"] != "row_baseline_and_absolute_exceptions_v1":
        raise ValueError("Unknown storage format.")
    baselines = data["row_baselines"]
    if len(baselines) != size + 1 or any(type(v) is not int or v < 0 for v in baselines):
        raise ValueError("Invalid row baselines.")
    rows = [[value] * size for value in baselines]
    seen = set()
    for item in data["exceptions"]:
        if len(item) != 3 or any(type(v) is not int for v in item):
            raise ValueError("Invalid exception.")
        i, j, value = item
        if not 0 <= i <= size or not 0 <= j < size or value < 0 or (i, j) in seen:
            raise ValueError("Invalid or duplicate exception.")
        seen.add((i, j))
        rows[i][j] = value
    return {"report": data["report"], "entry_denominator_bits": data["entry_denominator_bits"],
            "coupling_numerators": rows}


def read_candidate():
    data = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    expected = (2048, .14473, .98961702)
    actual = tuple(data["report"][key] for key in ("size", "contrast", "source_radius"))
    if actual != expected or data["entry_denominator_bits"] != 80:
        raise ValueError("Unexpected saved research case.")
    return expand_candidate(data)


@lru_cache(maxsize=1)
def numerical_coupling():
    return polygon.repaired_coupling_float(read_candidate())


def reverse_probabilities(angles, selected_columns=None):
    matrix, weights, _, _ = numerical_coupling()
    if selected_columns is not None:
        matrix = matrix[:, selected_columns]
    conditional = matrix / weights[:, None]
    left, lower, upper, center = polygon.quantize_source(angles, 2048, .98961702)
    return (lower[..., None] * conditional[left] + upper[..., None] * conditional[(left + 1) % 2048]
            + center[..., None] * conditional[-1])


class CompressedCandidateTests(unittest.TestCase):
    def test_storage_preserves_every_quantized_entry(self):
        matrix = np.array([[.1, .1, .2, .1], [.01, .03, .01, .01], [.2] * 4,
                           [.04] * 4, [.11, .12, .11, .11]])
        data = compress_candidate({"size": 4}, matrix, 80)
        restored = expand_candidate(data)["coupling_numerators"]
        self.assertEqual(restored, [[int(value * (1 << 80)) for value in row] for row in matrix])
        self.assertEqual(len(data["exceptions"]), 3)

    def test_duplicate_exceptions_are_rejected(self):
        data = compress_candidate({"size": 4}, np.full((5, 4), .1))
        data["exceptions"] = [[0, 1, 0], [0, 1, 2]]
        with self.assertRaises(ValueError):
            expand_candidate(data)

    def test_saved_fine_candidate_has_positive_exact_repair(self):
        self.assertTrue(polygon.verify_polygon_data(read_candidate())["verified"])

    def test_corrupted_compressed_entry_is_rejected_by_exact_certificate(self):
        data = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
        data["exceptions"][0][2] += 1 << 65
        self.assertFalse(polygon.verify_polygon_data(expand_candidate(data))["verified"])

    def test_repaired_matrix_has_required_marginals_and_vector_means(self):
        matrix, weights, sources, targets = numerical_coupling()
        self.assertGreater(float(matrix.min()), 0)
        np.testing.assert_allclose(matrix.sum(axis=1), weights, atol=3e-17)
        np.testing.assert_allclose(matrix.sum(axis=0), 1 / 2048, atol=3e-17)
        np.testing.assert_allclose(matrix @ targets, weights[:, None] * sources, atol=3e-17)

    def test_complete_selective_output_at_new_strength(self):
        from compensated_angle_kernel import warped_angle
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        from outcome_updates import WINDOW_ATTENUATION as alpha
        grid = np.linspace(-1, 8, 47)
        eta, size = .14473, 2048
        for setting in (0, .61):
            centers = np.arange(size) * 2 * math.pi / size + setting
            for mark in (0, 1):
                relative = grid - setting
                probabilities = reverse_probabilities(warped_angle(relative, mark, eta))
                f = (1 + (2 * mark - 1) * eta * np.cos(relative)) / 2
                for summary in ([1, alpha, 0], [1, 0, alpha], [.6, .2, -.1], [1, 0, 0]):
                    arc_mass = (summary[0] + np.sinc(1 / size) / alpha
                                * (summary[1] * np.cos(centers) + summary[2] * np.sin(centers))) / size
                    observed = f * size / (2 * math.pi) * (probabilities @ arc_mass)
                    expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, grid)
                    np.testing.assert_allclose(observed, expected, atol=2e-13)

    def test_forward_kernel_normalization_by_independent_quadrature(self):
        from atomic_reset_kernel import atom_target
        size, eta = 2048, .14473
        nodes, gauss = np.polynomial.legendre.leggauss(12)
        step = 2 * math.pi / size
        psi = ((np.arange(size)[:, None] + (nodes + 1) / 2) * step).ravel()
        weights = np.tile(gauss * step / 2, size)
        probabilities = reverse_probabilities(psi, [0, 512, 731])
        total = np.zeros(3)
        q = math.sqrt(1 - eta**2)
        for mark in (0, 1):
            sign = 2 * mark - 1
            y = atom_target(psi, mark, eta)
            factor = (1 + sign * eta * np.cos(y)) / 2 * q / (1 - sign * eta * np.cos(psi))
            total += ((weights * factor)[:, None] * probabilities).sum(axis=0) * size / (2 * math.pi)
        np.testing.assert_allclose(total, 1, atol=3e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--candidate-output", type=Path)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    if args.solve:
        if args.candidate_output is None or args.candidate_output.exists():
            parser.error("Supply a new --candidate-output path; preserve archived files.")
        from symmetric_polygon_solver import solve_symmetric_polygon
        report, coupling = solve_symmetric_polygon(2048, .14473, .98961702)
        print(json.dumps(report), flush=True)
        if coupling is None:
            raise SystemExit(1)
        data = compress_candidate(report, coupling)
        args.candidate_output.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")
        certificate = polygon.verify_polygon_data(expand_candidate(data))
        print(json.dumps(certificate), flush=True)
        if not certificate["verified"]:
            raise SystemExit(1)
        return
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CompressedCandidateTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    certificate = polygon.verify_polygon_data(read_candidate())
    certificate.update({"candidate_sha256": hashlib.sha256(CANDIDATE_PATH.read_bytes()).hexdigest(),
                        "environment": {"python": platform.python_version(), "numpy": np.__version__},
                        "automated_checks": {"run": checks.testsRun, "errors": len(checks.errors), "failures": len(checks.failures)},
                        "scope": "Lossless storage of an 80-bit dyadic matrix; exact analytic repair verified with 100-bit intervals. Source quantization and target arc lift give a continuous martingale. Round 34 gives every eta<=0.14473; endpoint remains open."})
    if args.write_results:
        Path(__file__).with_name("compressed_polygon_certificate_results.json").write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()
