"""Round 409: regularity bounds for exact autonomous programs on separable spaces."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "separable_program_regularity_audit_results.json"
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def op(a):
    return float(np.linalg.norm(a, 2))


def evolve(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def embed(p, d=2):
    return np.kron(np.eye(d), np.asarray(p, complex).reshape(-1, 1))


def ket(rng, d):
    a = rng.normal(size=d) + 1j * rng.normal(size=d)
    return a / np.linalg.norm(a)


def trace_distance_pure(a, b):
    return float(np.sqrt(max(0., 1. - abs(np.vdot(a, b))**2)))


def choi(u):
    v = u.reshape(-1) / np.sqrt(u.shape[0])
    return np.outer(v, v.conj())


def factorization_examples():
    p = np.array([np.sqrt(.6), np.sqrt(.4)], complex)
    q = np.array([np.cos(.72), np.sin(.72)], complex)
    a = np.vdot(p, q)
    examples = []
    h = np.kron(Y, I2) + np.kron(I2, Z)
    for label, hh, t, s, v, w, ep, eq in [
        ("short_separable_evolution", h, .2, .205,
         evolve(Y, .2), evolve(Y, .205), evolve(Z, .2) @ p, evolve(Z, .205) @ q),
        ("interacting_recurrence", np.kron(Z, Z), np.pi/2, np.pi,
         Z, I2, -1j * Z @ p, -q),
    ]:
        ut, us = evolve(hh, t), evolve(hh, s)
        tp, sq = ut @ embed(p), us @ embed(q)
        delta = op((us - ut) @ embed(q))
        b = np.vdot(ep, eq)
        residual = op(a * I2 - b * v.conj().T @ w)
        # Phase picked by the proof, not a numerical phase minimization.
        phase = (a / abs(a)) * np.conj(b / abs(b)) if abs(b) > 1e-14 else 1.
        examples.append(dict(
            name=label, time_difference=abs(t-s), program_overlap=float(abs(a)),
            factorization_error=max(op(tp - embed(ep) @ v), op(sq - embed(eq) @ w)),
            overlap_equation_residual=residual, elapsed_isometry_difference=delta,
            phase_operator_distance=op(v.conj().T @ w-phase*I2),
            phase_bound=2*delta/abs(a),
            choi_trace_distance=trace_distance_pure(v.reshape(-1)/np.sqrt(2), w.reshape(-1)/np.sqrt(2)),
            bounded_generator_upper=op(hh)*abs(s-t)))
    hc = np.kron(Z, np.diag([1., 0.])) + np.kron(X, np.diag([0., 1.]))
    u = evolve(hc, .23)
    e0, e1 = I2[:, 0], I2[:, 1]
    return dict(nonorthogonal=examples, same_time_orthogonal=dict(
        overlap=float(abs(np.vdot(e0, e1))),
        factorization_error=max(op(u@embed(e0)-embed(e0)@evolve(Z,.23)),
                                op(u@embed(e1)-embed(e1)@evolve(X,.23))),
        distinct_channel_choi_distance=trace_distance_pure(
            evolve(Z,.23).reshape(-1)/np.sqrt(2), evolve(X,.23).reshape(-1)/np.sqrt(2))))


def scalar_lemma_checks():
    rng = np.random.default_rng(40901)
    margins = []
    for d in (2, 3, 5):
        for k in range(16):
            z = np.diag(np.exp(1j*rng.uniform(-np.pi, np.pi, d)))
            a = rng.uniform(.2, 1)*np.exp(1j*rng.uniform(-np.pi, np.pi))
            b = 0j if k == 0 else rng.uniform(0, 1)*np.exp(1j*rng.uniform(-np.pi, np.pi))
            delta = op(a*np.eye(d)-b*z)
            phase = a/abs(a)*np.conj(b/abs(b)) if b else 1.
            margins.append(2*delta/abs(a)-op(z-phase*np.eye(d)))
    return dict(cases=len(margins), minimum_bound_margin=float(min(margins)),
                zero_final_program_overlap_included=True)


def holder_checks():
    rng = np.random.default_rng(40902)
    z = rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
    v, _ = np.linalg.qr(z)
    energies = np.array([0., .1, 1., 10., 100., 1000.])
    h = (v*energies)@v.conj().T
    q = ket(rng, 3)
    b = embed(q, 2)
    rows = []
    for alpha in (.5, 1.):
        ha = (v*(energies**alpha))@v.conj().T
        moment = float(np.linalg.norm(ha@b, "fro"))
        for tau in (1e-5, .003, .1):
            observed = op((evolve(h, tau)-np.eye(6))@b)
            bound = 2**(1-alpha)*moment*tau**alpha
            rows.append(dict(alpha=alpha, time_difference=tau,
                             basis_moment=moment, isometry_difference=observed,
                             spectral_theorem_bound=bound, ratio=observed/bound))
    return rows


def reference_checks():
    rng = np.random.default_rng(40903)
    g = np.array([[.2, .3j, .1],[-.3j,-.4,.2j],[.1,-.2j,.7]],complex)
    v, w = evolve(g, .12), evolve(g, .125)
    # Same program is allowed; overlap a=1. The displayed bound is conservative.
    delta = op(w-v)
    rows = []
    for rd in (1, 3, 7):
        psi = ket(rng, 3*rd)
        lhs = trace_distance_pure(np.kron(v,np.eye(rd))@psi, np.kron(w,np.eye(rd))@psi)
        rows.append(dict(reference_dimension=rd, pure_joint_trace_distance=lhs,
                         overlap_channel_bound=min(1.,2*delta)))
    return rows


def euler(a, b, c):
    def rotation(p, x):
        return np.cos(x)*I2-1j*np.sin(x)*p
    return rotation(Z,a)@rotation(Y,b)@rotation(Z,c)


def positive_channel_log(u):
    """A positive generator of this SU(2) channel, allowing an overall phase."""
    co = float(np.clip(np.trace(u).real/2, -1, 1))
    theta = float(np.arccos(co))
    anti = (u-u.conj().T)/(-2j)
    si = float(np.sin(theta))
    if abs(si) < 1e-10:
        return np.zeros((2,2),complex) if co > 0 else np.pi*I2
    return theta*I2+(theta/si)*anti


def countable_catalogue_checks():
    rng = np.random.default_rng(40904)
    targets = [rng.uniform(0,2*np.pi,3) for _ in range(4)]
    rows = []
    for m in (16, 64):
        for angles in targets:
            rounded = (np.floor(angles*m/(2*np.pi)+.5).astype(int)%m)
            grid = 2*np.pi*rounded/m
            target, approx = euler(*angles), euler(*grid)
            h = positive_channel_log(approx)
            actual = evolve(h,1)
            gap = op(target-approx)
            rows.append(dict(
                denominator=m, integer_angles=rounded.tolist(),
                positive_generator_min=float(np.linalg.eigvalsh(h)[0]),
                generator_norm=op(h),
                exact_program_choi_error=op(choi(actual)-choi(approx)),
                unitary_grid_error=gap, approximation_bound=3*np.pi/m))
    # A finite restriction checks the direct-sum wiring, not the infinite theorem.
    hs = [positive_channel_log(euler(2*np.pi*a/7,2*np.pi*b/7,2*np.pi*c/7))
          for a,b,c in ((1,2,3),(2,1,4),(4,3,1))]
    full = np.zeros((6,6),complex)
    for j,h in enumerate(hs):
        full += np.kron(h, np.diag([float(k==j) for k in range(3)]))
    u = evolve(full,1)
    wiring = max(op(u@embed(np.eye(3)[:,j])-embed(np.eye(3)[:,j])@evolve(hs[j],1))
                 for j in range(3))
    return dict(samples=rows, finite_direct_sum_wiring_error=wiring,
                infinite_construction="countable direct sum over all rational Euler triples",
                global_generator_norm_upper=2*np.pi,
                arbitrary_accuracy_from_dense_catalogue_proved_analytically=True,
                exact_universal_unitary_catalogue_claimed=False)


def infinite_energy_example():
    # p_n = n^(-3/2)/sqrt(zeta(3)); H=|1><1| tensor diag(n), n>=1.
    # Integral remainders enclose the omitted positive series. Floating endpoints
    # are diagnostics, not interval-arithmetic certificates.
    nmax = 4096
    n = np.arange(1,nmax+1,dtype=float)
    z3prefix = float(np.sum(n**-3))
    z2prefix = float(np.sum(n**-2))
    z3lo = z3prefix+1/(2*(nmax+1)**2)
    z3hi = z3prefix+1/(2*nmax**2)
    energylo = (z2prefix+1/(nmax+1))/z3hi
    energyhi = (z2prefix+1/nmax)/z3lo
    rows = []
    for tau in (.001,.01,.2):
        observed_prefix = float(np.sum(4*np.sin(n*tau/2)**2/n**3))
        upper = (observed_prefix+4/(2*nmax**2))/z3lo
        rows.append(dict(time_difference=tau,
                         squared_isometry_difference_lower=observed_prefix/z3hi,
                         squared_isometry_difference_upper=upper,
                         holder_squared_bound_using_energy_lower=2*energylo*tau))
    return dict(prefix_terms=nmax, normalization_interval=[z3lo,z3hi],
                mean_energy_interval=[energylo,energyhi],
                second_moment_partial_sums=[float(np.sum(np.arange(1,k+1,dtype=float)**-1)/z3hi)
                                           for k in (64,nmax)],
                second_moment_diverges_by_harmonic_series=True,
                inequalities=rows,
                floating_endpoints_are_interval_certificates=False)


@lru_cache(maxsize=1)
def report():
    return dict(
        date="2026-09-23", round=409, scientific_base_through_round=408,
        scope="One fixed autonomous self-adjoint generator, finite data, separable normal independent programs, exact deterministic unitary output; no intermediate outside control.",
        old_test_suites_executed=False,
        factorization=factorization_examples(),
        scalar_phase_lemma=scalar_lemma_checks(),
        spectral_holder=holder_checks(),
        arbitrary_reference=reference_checks(),
        countable_dense_catalogue=countable_catalogue_checks(),
        unbounded_finite_mean_energy=infinite_energy_example(),
        bounded_generator_exact_channel_hausdorff_dimension_upper=1,
        finite_mean_energy_exact_channel_hausdorff_dimension_upper=2,
        uniform_program_energy_cap_required=False,
        finite_program_dimension_required=False,
        program_selection_continuity_required=False,
        approximation_excluded=False, spatial_dimension_generated=False,
        full_cognition_to_gr_refuted=False,
        strengthened_exact_autonomous_control_adopted_as_axiom=False,
        initial_resource_amount_choice_required=False,
        infinite_dimension_theorem_established_by_finite_numerics=False,
        analytic_proof_in_note=True)


class Audit(unittest.TestCase):
    def test_distinct_times_and_orthogonal_exception(self):
        r=report()["factorization"]
        for row in r["nonorthogonal"]:
            self.assertGreater(row["program_overlap"],.5)
            self.assertLess(row["factorization_error"],1e-12)
            self.assertLessEqual(row["overlap_equation_residual"],row["elapsed_isometry_difference"]+1e-12)
            self.assertLessEqual(row["elapsed_isometry_difference"],row["bounded_generator_upper"]+1e-12)
            self.assertLessEqual(row["phase_operator_distance"],row["phase_bound"]+1e-12)
        self.assertEqual(r["same_time_orthogonal"]["overlap"],0)
        self.assertLess(r["same_time_orthogonal"]["factorization_error"],1e-12)
        self.assertGreater(r["same_time_orthogonal"]["distinct_channel_choi_distance"],.1)

    def test_complex_scalar_to_phase_lemma(self):
        self.assertEqual(report()["scalar_phase_lemma"]["cases"],48)
        self.assertGreaterEqual(report()["scalar_phase_lemma"]["minimum_bound_margin"],-1e-12)

    def test_fractional_energy_moment_bounds(self):
        for row in report()["spectral_holder"]:
            self.assertLessEqual(row["isometry_difference"],row["spectral_theorem_bound"]+1e-10)

    def test_unknown_reference_bounds(self):
        for row in report()["arbitrary_reference"]:
            self.assertLessEqual(row["pure_joint_trace_distance"],row["overlap_channel_bound"]+1e-11)

    def test_countable_exact_catalogue_still_allows_approximation(self):
        r=report()["countable_dense_catalogue"]
        self.assertLess(r["finite_direct_sum_wiring_error"],1e-12)
        for row in r["samples"]:
            self.assertGreaterEqual(row["positive_generator_min"],-1e-12)
            self.assertLessEqual(row["generator_norm"],2*np.pi+1e-12)
            self.assertLess(row["exact_program_choi_error"],1e-12)
            self.assertLessEqual(row["unitary_grid_error"],row["approximation_bound"]+1e-12)

    def test_infinite_normal_program_with_finite_mean_not_variance(self):
        r=report()["unbounded_finite_mean_energy"]
        self.assertLess(r["mean_energy_interval"][1]-r["mean_energy_interval"][0],1e-6)
        self.assertGreater(r["second_moment_partial_sums"][1],r["second_moment_partial_sums"][0]+3)
        for row in r["inequalities"]:
            self.assertLessEqual(row["squared_isometry_difference_upper"],
                                 row["holder_squared_bound_using_energy_lower"])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    tests=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    result=dict(report())
    result["checks"]=dict(run=tests.testsRun,failures=len(tests.failures),errors=len(tests.errors))
    result["runtime"]=dict(python=platform.python_version(),numpy=np.__version__)
    if args.write_results:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
