"""Round 369: modular inclusion prerequisites and finite affine approximations."""
import unittest
import numpy as np
from growing_stream_audit import main

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def exp_hermitian(h, coefficient):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(coefficient*values)) @ vectors.conjugate().T


def thermal(h):
    rho = exp_hermitian(h, -1.)
    return rho/np.trace(rho)


def matrix_log(rho):
    values, vectors = np.linalg.eigh(rho)
    if values.min() <= 0:
        raise ValueError("Faithful reference state required")
    return (vectors*np.log(values)) @ vectors.conjugate().T


def project_first_algebra(a):
    # Trace-preserving Hilbert--Schmidt projection onto M_2 tensor I_2.
    reduced = np.trace(a.reshape(2, 2, 2, 2), axis1=1, axis2=3)/2
    return np.kron(reduced, I2)


def modular_leakage(rho, time, observable):
    u = exp_hermitian(matrix_log(rho), 1j*time)
    evolved = u @ observable @ u.conjugate().T
    return float(np.linalg.norm(evolved-project_first_algebra(evolved), "fro")/
                 np.linalg.norm(observable, "fro"))


def cyclic_ranks(rho):
    omega = exp_hermitian(matrix_log(rho), .5)
    full, small = [], []
    for i in range(4):
        for j in range(4):
            e = np.zeros((4, 4), complex); e[i, j] = 1
            full.append((e @ omega).ravel())
    for i in range(2):
        for j in range(2):
            e = np.zeros((2, 2), complex); e[i, j] = 1
            small.append((np.kron(e, I2) @ omega).ravel())
    return {"Hilbert_space_dimension": 16,
            "full_algebra_cyclic_rank": int(np.linalg.matrix_rank(np.column_stack(full))),
            "proper_subalgebra_cyclic_rank": int(np.linalg.matrix_rank(np.column_stack(small)))}


def affine_matrix_witness(seed=369, dimension=5, time=.7):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
    k = (raw+raw.conjugate().T)/2
    raw = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
    p = raw @ raw.conjugate().T+np.eye(dimension)
    p /= np.linalg.norm(p, 2)
    comm = k @ p-p @ k
    residual = comm-1j*p
    u = exp_hermitian(k, 1j*time)
    diff = u @ p @ u.conjugate().T-np.exp(-time)*p
    return {"dimension": dimension, "time": time,
            "P_norm": float(np.linalg.norm(p, 2)),
            "P_HS_norm": float(np.linalg.norm(p, "fro")),
            "commutator_HS_norm": float(np.linalg.norm(comm, "fro")),
            "affine_residual_HS_norm": float(np.linalg.norm(residual, "fro")),
            "squared_orthogonality_residual": float(abs(np.linalg.norm(residual, "fro")**2-
                np.linalg.norm(comm, "fro")**2-np.linalg.norm(p, "fro")**2)),
            "finite_flow_operator_error": float(np.linalg.norm(diff, 2)),
            "finite_flow_operator_lower_bound": float(abs(1-np.exp(-time))*np.linalg.norm(p, 2))}


def window_code(halfwidth, spacing=1/16):
    n = int(round(2*halfwidth/spacing))
    q = -halfwidth+spacing*np.arange(n)
    gaussian = np.exp(-(q+.5)**2/(4*.65**2))*np.exp(.4j*q)
    raw = np.column_stack([gaussian, (q+.5)*gaussian])
    code, _ = np.linalg.qr(raw)
    return q, code


def shift_fft(values, time, spacing):
    k = 2*np.pi*np.fft.fftfreq(len(values), d=spacing)
    multiplier = np.exp(-1j*time*k)
    if values.ndim == 2:
        multiplier = multiplier[:, None]
    return np.fft.ifft(np.fft.fft(values, axis=0)*multiplier, axis=0)


def window_row(halfwidth, time=.75, translation=.4):
    spacing = 1/16
    q, code = window_code(halfwidth, spacing)
    steps = int(round(time/spacing))
    if not np.isclose(time, steps*spacing):
        raise ValueError("Exact grid translation required for the tail certificate")
    p = np.exp(q)
    actual = np.roll(np.exp(1j*translation*p)[:, None]*
                     np.roll(code, -steps, axis=0), steps, axis=0)
    expected = np.exp(1j*translation*np.exp(-time)*p)[:, None]*code
    wrap = q < -halfwidth+time
    tail = float(np.linalg.norm(code[wrap, :], 2))
    residual = actual-expected
    a, b = actual.ravel()/np.sqrt(2), expected.ravel()/np.sqrt(2)
    # Stable pure-state trace distance; avoids 1-|overlap|^2 cancellation.
    overlap = np.vdot(a, b)
    reference_distance = float(np.linalg.norm(b-overlap*a))
    return {"L": halfwidth, "N": len(q), "grid_spacing": spacing,
            "dilation_parameter_t": time, "translation_parameter_a": translation,
            "window_code_dimension": 2,
            "uniform_reference_half_diamond_bound": min(1., 2*tail),
            "actual_input_code_operator_error": float(np.linalg.norm(residual, 2)),
            "entangled_reference_witness_distance": reference_distance,
            "P_norm": float(p.max()),
            "generator_covariance_on_code_error": float(np.linalg.norm(
                (np.roll(p, steps)-np.exp(-time)*p)[:, None]*code, 2)),
            "full_generator_operator_error_lower": float((1-np.exp(-time))*p.max()),
            "K_norm": float(np.pi/spacing),
            "log_P_dynamic_range": float(q[-1]-q[0])}


def report():
    eta = .7
    correlated = thermal(eta*np.kron(Z, Z))
    product = np.kron(thermal(.3*Z), thermal(-.2*Z))
    times = (.2, .7, 1.1)
    return {"round": 369,
            "scope": "Finite matrix prerequisite audit plus a specified positive affine-group representation and windowed unitary approximation. Not a construction of a half-sided modular inclusion, a derived vacuum, physical dimension, or Einstein dynamics.",
            "hypothesis_tested": "Do nested finite quantum algebras and a faithful shared state already provide the modular inclusion prerequisites that generate nontrivial translations?",
            "finite_common_vector_dimension_obstruction": cyclic_ranks(correlated),
            "modular_inclusion_leakage": [
                {"t": t, "product_reference_leakage": modular_leakage(product, t, np.kron(X, I2)),
                 "correlated_reference_leakage": modular_leakage(correlated, t, np.kron(X, I2)),
                 "correlated_exact": float(abs(np.sin(2*eta*t)))} for t in times],
            "finite_affine_obstructions": [affine_matrix_witness(369+n, n) for n in (2, 5, 8)],
            "windowed_affine_group": [window_row(L) for L in (3, 4, 5, 6, 8)],
            "analytically_proved": [
                "A finite algebra inclusion with a common cyclic and separating vector must be equality.",
                "A finite-dimensional subalgebra invariant by inclusion under a reversible flow for one half of time is invariant by equality for both halves.",
                "For finite Hermitian K,P, ||[K,P]-iP||_HS^2=||[K,P]||_HS^2+||P||_HS^2.",
                "Bounded nonzero P cannot obey exact dilation covariance; its full operator error is at least |1-exp(-t)| ||P||.",
                "Grid affine-unitary error on an input code is bounded by twice its wrap-region amplitude, uniformly for every reference."],
            "additional_inputs": ["nested observable algebras and their representation",
                                  "a faithful reference state; no vacuum selection",
                                  "for the positive comparison: L2(R,dq), P=exp(q), K=i*d/dq",
                                  "finite-grid boundary condition and a selected two-mode input window"],
            "not_proved": ["the supplied affine representation comes from common modular data",
                           "the grid generators are modular logarithms of the required algebras",
                           "the full generator norm converges",
                           "a fixed-spacing table is a continuum convergence theorem",
                           "one translation and dilation generate 3+1 geometry or gravity"],
            "sources": ["https://arxiv.org/html/math/0412061",
                        "https://www.numdam.org/article/AIHPA_1995__63_4_331_0.pdf"]}


class Checks(unittest.TestCase):
    def test_01_gns_common_cyclic_vector_dimension_gap(self):
        row = cyclic_ranks(thermal(.7*np.kron(Z, Z)))
        self.assertEqual(row, {"Hilbert_space_dimension": 16,
                              "full_algebra_cyclic_rank": 16,
                              "proper_subalgebra_cyclic_rank": 4})

    def test_02_algebra_projection_is_idempotent_and_trace_preserving(self):
        rng = np.random.default_rng(1369)
        a = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        projected = project_first_algebra(a)
        np.testing.assert_allclose(project_first_algebra(projected), projected, atol=1e-14)
        self.assertAlmostEqual(abs(np.trace(projected)-np.trace(a)), 0., places=13)
        for local in (I2, X, Y, Z):
            self.assertLess(abs(np.trace(np.kron(local, I2).conjugate().T @
                                         (a-projected))), 1e-13)

    def test_03_product_state_preserves_subalgebra_both_time_directions(self):
        rho = np.kron(thermal(.3*Z+.1*X), thermal(-.2*Y))
        for t in (-1.1, -.2, .2, 1.1):
            for a in (I2, X, Y, Z):
                self.assertLess(modular_leakage(rho, t, np.kron(a, I2)), 1e-13)

    def test_04_correlated_reference_fails_half_sided_inclusion(self):
        eta = .7; rho = thermal(eta*np.kron(Z, Z))
        for t in (-1.1, -.2, .2, 1.1):
            self.assertAlmostEqual(modular_leakage(rho, t, np.kron(X, I2)),
                                   abs(np.sin(2*eta*t)), places=13)

    def test_05_faithfulness_does_not_imply_modular_inclusion(self):
        rho = thermal(.7*np.kron(Z, Z))
        self.assertGreater(np.linalg.eigvalsh(rho).min(), 0)
        self.assertGreater(modular_leakage(rho, .7, np.kron(X, I2)), .8)
        marginal = np.trace(rho.reshape(2, 2, 2, 2), axis1=1, axis2=3)
        np.testing.assert_allclose(marginal, I2/2, atol=1e-14)

    def test_06_finite_affine_residual_exact_pythagoras(self):
        for dimension in (2, 5, 8):
            row = affine_matrix_witness(369+dimension, dimension)
            self.assertLess(row["squared_orthogonality_residual"], 1e-12)
            self.assertGreaterEqual(row["affine_residual_HS_norm"], row["P_HS_norm"])

    def test_07_bounded_generator_norm_obstruction(self):
        for t in (-.4, .2, .7):
            row = affine_matrix_witness(2369, 7, t)
            self.assertGreaterEqual(row["finite_flow_operator_error"]+1e-13,
                                    row["finite_flow_operator_lower_bound"])

    def test_08_commuting_generator_attains_HS_lower_bound(self):
        p = np.diag([.1, .3, 1.]); k = np.diag([.7, -.5, .2])
        comm = k @ p-p @ k
        self.assertEqual(np.linalg.norm(comm-1j*p, "fro"), np.linalg.norm(p, "fro"))

    def test_09_exact_grid_shift_matches_independent_fourier_unitary(self):
        q, code = window_code(4)
        for steps in (-12, 0, 7, 12):
            np.testing.assert_allclose(shift_fft(code, steps/16, 1/16),
                                       np.roll(code, steps, axis=0), atol=2e-15)

    def test_10_group_covariance_error_is_supported_only_on_wrap(self):
        q, code = window_code(4); time=.75; steps=12; a=.4
        actual = np.roll(np.exp(1j*a*np.exp(q))[:, None]*np.roll(code,-steps,axis=0),
                         steps,axis=0)
        desired = np.exp(1j*a*np.exp(q-time))[:, None]*code
        np.testing.assert_allclose((actual-desired)[q >= -4+time], 0, atol=2e-15)

    def test_11_window_bound_and_reference_witness(self):
        for L in (3, 4, 5, 6, 8):
            row = window_row(L)
            self.assertLessEqual(row["actual_input_code_operator_error"],
                                 2 if row["uniform_reference_half_diamond_bound"]==1 else
                                 row["uniform_reference_half_diamond_bound"]+3e-15)
            self.assertLessEqual(row["entangled_reference_witness_distance"],
                                 row["uniform_reference_half_diamond_bound"]+3e-15)
        self.assertLess(window_row(8)["actual_input_code_operator_error"], 1e-10)
        self.assertGreater(window_row(8)["full_generator_operator_error_lower"], 1000.)

    def test_12_compact_window_is_exact_at_allowed_grid_times(self):
        q, _ = window_code(4)
        rng = np.random.default_rng(3369)
        raw = np.zeros((len(q), 3), complex)
        mask = abs(q) < 1
        raw[mask] = rng.normal(size=(mask.sum(), 3))+1j*rng.normal(size=(mask.sum(), 3))
        code, _ = np.linalg.qr(raw)
        for steps in (2, 12, 24):
            t=steps/16
            for a in (-1., .2, 3.):
                actual=np.roll(np.exp(1j*a*np.exp(q))[:,None]*
                               np.roll(code,-steps,axis=0),steps,axis=0)
                desired=np.exp(1j*a*np.exp(q-t))[:,None]*code
                np.testing.assert_allclose(actual,desired,atol=3e-15)

    def test_13_full_positive_P_and_selfadjoint_K_are_supplied_not_modular(self):
        q, _ = window_code(3)
        k = 2*np.pi*np.fft.fftfreq(len(q), d=1/16)
        eye = np.eye(len(q), dtype=complex)
        generator = np.fft.ifft(-k[:,None]*np.fft.fft(eye,axis=0),axis=0)
        np.testing.assert_allclose(generator, generator.conjugate().T, atol=2e-14)
        self.assertGreater(np.exp(q).min(), 0)
        self.assertAlmostEqual(np.linalg.norm(generator, 2), 16*np.pi, places=12)

    def test_14_no_uniform_improvement_from_window_success(self):
        rows = [window_row(L) for L in (3, 4, 5, 6, 8)]
        self.assertTrue(all(b["uniform_reference_half_diamond_bound"] <
                            a["uniform_reference_half_diamond_bound"]
                            for a,b in zip(rows,rows[1:])))
        self.assertTrue(all(b["full_generator_operator_error_lower"] >
                            a["full_generator_operator_error_lower"]
                            for a,b in zip(rows,rows[1:])))


if __name__ == "__main__":
    main(__name__, "modular_inclusion_limit_audit", report)
