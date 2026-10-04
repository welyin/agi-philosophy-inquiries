"""Round 406: variable completion time in a fixed finite quantum processor.

Finite matrix checks support, but do not prove, the analytic-set dimension
theorem in research_note_406.md. No earlier scientific test is rerun.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name("variable_completion_time_audit_results.json")
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
P0, P1 = np.diag([1., 0.]), np.diag([0., 1.])
# Tensor ordering: data, program.
H = np.kron(Z, P0)+np.kron(X, P1)


def exp_h(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*time*values))@vectors.conj().T


def rotation(axis, time):
    return np.cos(time)*I-1j*np.sin(time)*axis


def channel_on_reference(w, program, rho):
    refdim = rho.shape[0]//2
    source = np.einsum("arbs,ij->airbjs", rho.reshape(2, refdim, 2, refdim),
                       program).reshape(4*refdim, 4*refdim)
    whole = np.kron(w, np.eye(refdim))
    output = (whole@source@whole.conj().T).reshape(2, 2, refdim, 2, 2, refdim)
    return np.einsum("airbis->arbs", output).reshape(2*refdim, 2*refdim)


def unitary_choi(u):
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    vector = np.kron(u, I)@bell
    return np.outer(vector, vector.conj())


def random_density(dim, rng):
    a = rng.normal(size=(dim, dim))+1j*rng.normal(size=(dim, dim))
    rho = a@a.conj().T
    return rho/np.trace(rho)


def euler(a, b, c):
    return rotation(Z, a)@rotation(X, b)@rotation(Z, c)


def euler_tangents(a, b, c):
    u = euler(a, b, c)
    derivatives = [-1j*Z@u,
                   rotation(Z, a)@(-1j*X@rotation(X, b))@rotation(Z, c),
                   u@(-1j*Z)]
    return np.array([[np.trace(axis@(1j*du@u.conj().T)).real/2
                      for du in derivatives] for axis in (X, Y, Z)])


@lru_cache(None)
def report():
    rng = np.random.default_rng(406)
    identities, purity, fidelity = [], [], []
    bell_rho = unitary_choi(I)
    target = rotation(Y, np.pi/4)
    target_choi = unitary_choi(target)
    times = (0., .17, .61, np.pi/2, np.pi, 2.31, 10.7)
    for time in times:
        w = exp_h(H, time)
        uz, ux = rotation(Z, time), rotation(X, time)
        expected = np.kron(uz, P0)+np.kron(ux, P1)
        for j in range(3):
            program = random_density(2, rng)
            state = random_density(6, rng)  # Arbitrary data/reference correlation.
            p = float(program[0, 0].real)
            actual = channel_on_reference(w, program, state)
            wz, wx = np.kron(uz, np.eye(3)), np.kron(ux, np.eye(3))
            mixture = p*wz@state@wz.conj().T+(1-p)*wx@state@wx.conj().T
            identities.append(dict(time=time, case=j,
                global_unitary_error=float(np.linalg.norm(w.conj().T@w-np.eye(4))),
                block_formula_error=float(np.linalg.norm(w-expected)),
                full_reference_channel_error=float(np.linalg.norm(actual-mixture)),
                unchanged_reference_error=float(np.linalg.norm(
                    np.einsum("aras->rs", actual.reshape(2, 3, 2, 3))-
                    np.einsum("aras->rs", state.reshape(2, 3, 2, 3))))))
        for p in (0., .2, .5, .8, 1.):
            program_vector = np.array([np.sqrt(p), np.exp(.7j)*np.sqrt(1-p)])
            program = np.outer(program_vector, program_vector.conj())
            choi = channel_on_reference(w, program, bell_rho)
            expected_purity = 1-2*p*(1-p)*(1-np.cos(time)**4)
            observed_purity = float(np.trace(choi@choi).real)
            overlap = float(np.trace(target_choi@choi).real)
            trace_distance = float(np.abs(np.linalg.eigvalsh(choi-target_choi)).sum()/2)
            purity.append(dict(time=time, program_weight=p, choi_purity=observed_purity,
                expected_purity=float(expected_purity),
                exact_unitary_expected=bool(p in (0., 1.) or abs(np.sin(time)) < 1e-12)))
            fidelity.append(dict(time=time, program_weight=p, target_choi_overlap=overlap,
                predicted_overlap=float(np.cos(time)**2/2),
                witness_lower_bound=1-overlap, choi_trace_distance=trace_distance))
    curve = []
    for time in (.1, .3, .5, .7):
        p = P0
        obtained = channel_on_reference(exp_h(H, time), p, bell_rho)
        curve.append(dict(time=time, exact_channel_error=float(np.linalg.norm(
            obtained-unitary_choi(rotation(Z, time))))))
    distinct = min(float(np.linalg.norm(unitary_choi(rotation(Z, a))-
                                       unitary_choi(rotation(Z, b))))
                   for a in (.1, .3, .5, .7) for b in (.1, .3, .5, .7) if a < b)
    eulers = []
    for angles in ((.2, .4, .7), (.13, 0., .21), (.3, np.pi/2, .5)):
        tangents = euler_tangents(*angles)
        singular = np.linalg.svd(tangents, compute_uv=False)
        eulers.append(dict(angles=list(map(float, angles)),
                           tangent_singular_values=singular.tolist(),
                           rank=int(np.sum(singular > 1e-10)),
                           determinant=float(np.linalg.det(tangents)),
                           predicted_abs_determinant=float(abs(np.sin(2*angles[1])))))
    # These are three externally chosen stage durations, not a single run.
    external = euler(np.pi/4, np.pi/4, 3*np.pi/4)
    external_error = float(np.linalg.norm(unitary_choi(external)-target_choi))
    dense = []
    phase1, phase2 = np.pi/3, np.pi/5
    for cutoff in (128, 2048, 32768):
        ns = np.arange(cutoff+1)
        ts = phase1+2*np.pi*ns
        residuals = np.angle(np.exp(1j*(np.sqrt(2)*ts-phase2)))
        index = int(np.argmin(np.abs(residuals)))
        residual = float(residuals[index])
        dense.append(dict(integer_search_limit=cutoff, selected_integer=index,
            time=float(ts[index]), second_phase_residual=residual,
            half_diamond_error=float(abs(np.sin(residual/2)))))
    return dict(
        round=406,
        scope="For a fixed finite-dimensional deterministic processor, real-analytic dependence on m external continuous settings permits an exactly implemented unitary-channel set of Hausdorff dimension at most m. A single variable completion time is insufficient for all qubit unitaries. This is an implementation qualification, not a spatial dimension theorem or a general approximate-control no-go.",
        processor_reference_checks=identities, complete_choi_purity_checks=purity,
        missing_rotation_witness_checks=fidelity, variable_time_exact_curve_checks=curve,
        minimum_sampled_curve_separation=distinct, euler_parameter_checks=eulers,
        externally_switched_target_channel_error=external_error,
        irrational_time_approximation_checks=dense,
        variable_completion_time_included=True,
        fixed_finite_processor_single_run_exact_universality_excluded=True,
        continuous_program_state_coordinates_counted_as_external_settings=False,
        countably_many_finite_hardware_types_extension_proved=True,
        arbitrary_accuracy_at_unbounded_times_excluded=False,
        dense_reorientation_contract_excluded=False,
        finite_group_position_certificate_excluded=False,
        all_internal_resource_models_excluded=False,
        control_parameter_dimension_is_spatial_dimension=False,
        spatial_dimension_generated=False, full_cognition_to_gr_refuted=False)


class Audit(unittest.TestCase):
    def test_full_processor_and_unknown_reference_channel(self):
        for row in report()["processor_reference_checks"]:
            for key in ("global_unitary_error", "block_formula_error",
                        "full_reference_channel_error", "unchanged_reference_error"):
                self.assertLess(row[key], 1e-12)

    def test_unitary_extremality_with_coherent_programs(self):
        for row in report()["complete_choi_purity_checks"]:
            self.assertAlmostEqual(row["choi_purity"], row["expected_purity"], places=12)
            self.assertEqual(abs(row["choi_purity"]-1) < 1e-12, row["exact_unitary_expected"])

    def test_variable_time_really_gives_distinct_exact_gates(self):
        for row in report()["variable_time_exact_curve_checks"]:
            self.assertLess(row["exact_channel_error"], 1e-12)
        self.assertGreater(report()["minimum_sampled_curve_separation"], .1)

    def test_missing_y_rotation_has_complete_reference_witness(self):
        for row in report()["missing_rotation_witness_checks"]:
            self.assertAlmostEqual(row["target_choi_overlap"], row["predicted_overlap"], places=12)
            self.assertGreaterEqual(row["witness_lower_bound"], .5-1e-12)
            self.assertGreaterEqual(row["choi_trace_distance"], row["witness_lower_bound"]-1e-12)

    def test_three_external_parameters_and_chart_singularities(self):
        rows = report()["euler_parameter_checks"]
        self.assertEqual([row["rank"] for row in rows], [3, 2, 2])
        for row in rows:
            self.assertAlmostEqual(abs(row["determinant"]), row["predicted_abs_determinant"], places=12)

    def test_external_stage_switching_changes_the_contract(self):
        self.assertLess(report()["externally_switched_target_channel_error"], 1e-12)

    def test_exact_dimension_bound_does_not_forbid_dense_approximation(self):
        rows = report()["irrational_time_approximation_checks"]
        errors = [row["half_diamond_error"] for row in rows]
        self.assertTrue(all(a >= b for a, b in zip(errors, errors[1:])))
        self.assertLess(errors[-1], 1e-4)
        self.assertTrue(all(error > 0 for error in errors))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result["checks"] = dict(run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors))
    result["runtime"] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))

