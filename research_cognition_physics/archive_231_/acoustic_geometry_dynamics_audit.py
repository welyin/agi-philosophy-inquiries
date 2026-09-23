"""Round 355: autonomous GP flow, acoustic geometry and Einstein-source gap.

NumPy only. The order parameter is a mean-field variable, not a fundamental
nonlinear quantum state. All geometric and laboratory units are explicit.
"""
from dataclasses import dataclass
import unittest
import numpy as np
from growing_stream_audit import main


@dataclass(frozen=True)
class Parameters:
    mass: float = 1.
    hbar: float = 1.
    coupling: float = .01
    density: float = 100.
    expansion: float = .002

    @property
    def c0(self):
        return np.sqrt(self.coupling*self.density/self.mass)


DEFAULT = Parameters()


def background(t, p=DEFAULT):
    a = 1+p.expansion*np.asarray(t)
    if np.any(a <= 0):
        raise ValueError("The expanding solution requires a>0.")
    return a, p.expansion/a, p.density/a**3, p.c0/a**1.5


def order_parameter(t, x, p=DEFAULT):
    a, hubble, density, _ = background(t, p)
    r2 = np.dot(x, x)
    # Integral_0^t n(s) ds; this form also has the correct H0=0 limit.
    integral_density = p.density*t*(2+p.expansion*t)/(2*a*a)
    phase = p.mass*hubble*r2/(2*p.hbar)-p.coupling*integral_density/p.hbar
    return np.sqrt(density)*np.exp(1j*phase)


def analytic_gp_residual(t, x, p=DEFAULT):
    _, hubble, density, _ = background(t, p)
    psi = order_parameter(t, x, p)
    phase_t = -p.mass*hubble*hubble*np.dot(x, x)/(2*p.hbar)
    phase_t -= p.coupling*density/p.hbar
    psi_t = psi*(-1.5*hubble+1j*phase_t)
    laplacian = psi*(3j*p.mass*hubble/p.hbar -
                     (p.mass*hubble/p.hbar)**2*np.dot(x, x))
    return (1j*p.hbar*psi_t+p.hbar**2*laplacian/(2*p.mass) -
            p.coupling*abs(psi)**2*psi)/psi


def derivative_5(function, at, step):
    return (function(at-2*step)-8*function(at-step)+
            8*function(at+step)-function(at+2*step))/(12*step)


def finite_gp_residual(t, x, step, p=DEFAULT):
    psi = order_parameter(t, x, p)
    psi_t = derivative_5(lambda s: order_parameter(s, x, p), t, step)
    laplacian = 0j
    for axis in np.eye(3):
        laplacian += (-order_parameter(t, x+2*step*axis, p) +
                      16*order_parameter(t, x+step*axis, p)-30*psi +
                      16*order_parameter(t, x-step*axis, p) -
                      order_parameter(t, x-2*step*axis, p))/(12*step**2)
    return (1j*p.hbar*psi_t+p.hbar**2*laplacian/(2*p.mass) -
            p.coupling*abs(psi)**2*psi)/psi


def laboratory_energy(t, x, p=DEFAULT):
    _, hubble, density, _ = background(t, p)
    velocity = hubble*np.asarray(x)
    kinetic = .5*p.mass*density*np.dot(velocity, velocity)
    pressure = .5*p.coupling*density*density
    energy = kinetic+pressure
    return energy, (energy+pressure)*velocity


def local_balance_residuals(t, x, step, p=DEFAULT):
    """Finite differences on the actual density, velocity and energy fields."""
    def density_at(s):
        return background(s, p)[2]

    _, hubble, density, _ = background(t, p)
    continuity = derivative_5(density_at, t, step)+3*hubble*density
    velocity_t = derivative_5(
        lambda s: background(s, p)[1]*x, t, step)
    euler = velocity_t+hubble*hubble*x
    energy_t = derivative_5(
        lambda s: laboratory_energy(s, x, p)[0], t, step)
    flux_divergence = 0.
    for i, axis in enumerate(np.eye(3)):
        flux_divergence += derivative_5(
            lambda s: laboratory_energy(t, x+s*axis, p)[1][i], 0., step)
    return continuity, np.linalg.norm(euler), energy_t+flux_divergence


def acoustic_metric(t, x, p=DEFAULT, normalized=True):
    _, hubble, density, sound = background(t, p)
    velocity = hubble*np.asarray(x)
    metric = np.eye(4)
    metric[0, 0] = -sound**2+np.dot(velocity, velocity)
    metric[0, 1:] = metric[1:, 0] = -velocity
    factor = density/sound
    if normalized:
        factor /= p.density/p.c0
    return factor*metric


def comoving_metric_jets(t, p=DEFAULT):
    a, hubble, _, _ = background(t, p)
    powers = np.array([-4.5, .5, .5, .5])
    values = np.array([-p.c0**2, 1., 1., 1.])*a**powers
    return (np.diag(values), np.diag(values*powers*hubble),
            np.diag(values*powers*(powers-1)*hubble*hubble))


def curvature_from_metric_jets(metric, first, second):
    """Contract Christoffel symbols and their derivatives, no FLRW formula.

    The metric only depends on coordinate 0. R_{mn} uses the convention
    partial_a Gamma^a_mn - partial_n Gamma^a_ma
      + Gamma^a_ab Gamma^b_mn - Gamma^a_nb Gamma^b_ma.
    """
    inverse = np.linalg.inv(metric)
    inverse_first = -inverse @ first @ inverse
    connection = np.zeros((4, 4, 4))
    connection_first = np.zeros_like(connection)
    for upper in range(4):
        for i in range(4):
            for j in range(4):
                for ell in range(4):
                    term = ((i == 0)*first[ell, j]+(j == 0)*first[ell, i]
                            -(ell == 0)*first[i, j])
                    term_first = ((i == 0)*second[ell, j]+
                                  (j == 0)*second[ell, i] -
                                  (ell == 0)*second[i, j])
                    connection[upper, i, j] += .5*inverse[upper, ell]*term
                    connection_first[upper, i, j] += .5*(
                        inverse_first[upper, ell]*term +
                        inverse[upper, ell]*term_first)
    ricci = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            ricci[i, j] = connection_first[0, i, j]
            if j == 0:
                ricci[i, j] -= sum(connection_first[k, i, k] for k in range(4))
            for k in range(4):
                for ell in range(4):
                    ricci[i, j] += (connection[k, k, ell]*connection[ell, i, j]
                                    -connection[k, j, ell]*connection[ell, i, k])
    scalar = float(np.einsum("ij,ij", inverse, ricci))
    einstein = ricci-.5*metric*scalar
    orthonormal = np.diag(einstein)/abs(np.diag(metric))
    return {"ricci": ricci, "scalar": scalar,
            "einstein": einstein, "orthonormal_einstein": orthonormal}


def acoustic_proper_time(t, p=DEFAULT):
    a = background(t, p)[0]
    if p.expansion == 0:
        return p.c0*t
    return 4*p.c0/(5*p.expansion)*(1-a**(-1.25))


def acoustic_hubble(t, p=DEFAULT):
    a = background(t, p)[0]
    return p.expansion*a**1.25/(4*p.c0)


def mode_window(a, comoving_k, p=DEFAULT):
    sound = p.c0/a**1.5
    physical_k = comoving_k/a
    eta = p.hbar*physical_k/(2*p.mass*sound)
    wkb = (p.expansion/a)/(sound*physical_k)
    scattering_length = p.coupling*p.mass/(4*np.pi*p.hbar**2)
    density = p.density/a**3
    return {"a": float(a), "comoving_k": comoving_k,
            "eta_quantum_pressure": float(eta),
            "instantaneous_frequency_relative_shift": float(np.sqrt(1+eta**2)-1),
            "H_over_hydro_frequency": float(wkb),
            "gas_parameter": float(density*scattering_length**3)}


def mode_rhs(s, vector, k, include_quantum_pressure, p=DEFAULT):
    """s=H0*t, vector=(density contrast, k*velocity potential/c0)."""
    a = 1+s
    eta0 = p.hbar*k/(2*p.mass*p.c0)
    quantum = eta0**2/a**2 if include_quantum_pressure else 0.
    prefactor = p.c0*k/p.expansion
    return prefactor*np.array([vector[1]/a**2,
                               -(a**(-3)+quantum)*vector[0]])


def evolve_mode(k=.1, steps=2048, include_quantum_pressure=True, p=DEFAULT):
    """Classical linear perturbations; RK4 on 0<=H0*t<=1."""
    vector = np.array([1.+0j, -1j])
    step = 1/steps
    max_wronskian_error = 0.
    for index in range(steps):
        s = index*step
        k1 = mode_rhs(s, vector, k, include_quantum_pressure, p)
        k2 = mode_rhs(s+step/2, vector+step*k1/2, k, include_quantum_pressure, p)
        k3 = mode_rhs(s+step/2, vector+step*k2/2, k, include_quantum_pressure, p)
        k4 = mode_rhs(s+step, vector+step*k3, k, include_quantum_pressure, p)
        vector += step*(k1+2*k2+2*k3+k4)/6
        max_wronskian_error = max(max_wronskian_error,
                                  abs(np.imag(vector[0].conjugate()*vector[1])+1))
    return vector, max_wronskian_error


def report():
    gp_parameters = Parameters(.7, .8, .4, 3., .2)
    gp_point = np.array([1.2, .6, -.4])
    gp_rows = [{"step": h, "absolute_relative_GP_residual":
                float(abs(finite_gp_residual(.3, gp_point, h, gp_parameters)))}
               for h in (.04, .02, .01)]
    geometry_rows = []
    for a in (1., 1.25, 1.5, 2.):
        t = (a-1)/DEFAULT.expansion
        h = acoustic_hubble(t)
        actual = curvature_from_metric_jets(*comoving_metric_jets(t))
        geometry_rows.append({
            "a": a, "laboratory_time": t, "acoustic_proper_time": acoustic_proper_time(t),
            "acoustic_scale_factor": a**.25, "acoustic_H": h,
            "Ricci_scalar_from_connection": actual["scalar"],
            "Ricci_scalar_predicted": 42*h*h,
            "G_orthonormal_from_connection": actual["orthonormal_einstein"].tolist(),
            "G_null_contraction": float(actual["orthonormal_einstein"][:2].sum()),
            "atom_density": DEFAULT.density/a**3,
            "number_in_unit_comoving_volume": DEFAULT.density,
        })
    comparison = []
    for k in (.1, .5):
        full, w_full = evolve_mode(k, 4096, True)
        hydro, w_hydro = evolve_mode(k, 4096, False)
        finer, _ = evolve_mode(k, 8192, True)
        comparison.append({
            "comoving_k": k, "final_a": 2.,
            "full_vs_hydro_final_vector_difference": float(np.linalg.norm(full-hydro)),
            "full_RK4_refinement_difference": float(np.linalg.norm(full-finer)),
            "max_full_wronskian_error": w_full, "max_hydro_wronskian_error": w_hydro,
            "final_window": mode_window(2., k),
        })
    convergence = []
    reference, _ = evolve_mode(steps=4096)
    for steps in (256, 512, 1024, 2048):
        value, wronskian = evolve_mode(steps=steps)
        convergence.append({"steps": steps, "difference_to_4096":
                            float(np.linalg.norm(value-reference)),
                            "wronskian_error": wronskian})
    return {
        "round": 355,
        "scope": ("Exact non-normalizable homogeneous expanding solution of a specified "
                  "3D repulsive Gross-Pitaevskii mean-field model; low-k acoustic metric, "
                  "independent curvature contraction, finite perturbation checks. "
                  "Not a full FUCP countermodel, closed finite universe, fundamental "
                  "nonlinear quantum theory, or derivation/exclusion of all sourced GR."),
        "parameters": DEFAULT.__dict__,
        "coupling_convention": "g=4*pi*hbar^2*a_s/m; n=|Psi|^2 is number density",
        "GP_difference_parameters": gp_parameters.__dict__,
        "GP_residual_convergence": gp_rows,
        "geometry": geometry_rows,
        "mode_windows": [mode_window(a, .1) for a in (1., 2., 4., 16., 400.)],
        "mode_comparison": comparison,
        "mode_RK4_convergence": convergence,
        "metric_normalization": "(n/c)/(n0/c0); proper time has initial d_tau=c0*d_t",
        "source_obstruction": {
            "G00_over_Hac_squared": 3., "Gii_over_Hac_squared": -13.,
            "G_null_over_Hac_squared": -10.,
            "posthoc_effective_w": -13/3,
            "claim": ("No vacuum Einstein+constant Lambda solution for H0!=0; "
                      "no positive-kappa Einstein source obeying acoustic null NEC. "
                      "Atom stress is NOT assumed to be an acoustic-metric source."),
        },
        "added_inputs": ["3D laboratory space and time", "GP approximation and fixed g,m,hbar",
                         "homogeneous density and expanding initial phase",
                         "linear low-k probe sector", "fixed acoustic conformal normalization"],
        "sources": ["https://doi.org/10.1103/PhysRevLett.46.1351",
                    "https://arxiv.org/html/gr-qc/9712010",
                    "https://arxiv.org/html/gr-qc/0011026",
                    "https://homepages.ecs.vuw.ac.nz/~visser/Articles/Journals/ijmpd-frw.pdf"],
    }


class Checks(unittest.TestCase):
    def test_01_exact_GP_substitution_and_static_limit(self):
        for p in (DEFAULT, Parameters(.7, .8, .4, 3., .2), Parameters(expansion=0.)):
            for t in (0., .3, 2.):
                self.assertLess(abs(analytic_gp_residual(t, np.array([1.2, -.4, .7]), p)),
                                2e-14)
        self.assertEqual(acoustic_proper_time(2., Parameters(expansion=0.)), 2.)

    def test_02_actual_order_parameter_finite_differences_converge(self):
        p = Parameters(.7, .8, .4, 3., .2)
        x = np.array([1.2, .6, -.4])
        errors = [abs(finite_gp_residual(.3, x, h, p)) for h in (.04, .02, .01)]
        self.assertGreater(errors[0]/errors[1], 12.)
        self.assertGreater(errors[1]/errors[2], 10.)
        self.assertLess(errors[-1], 2e-8)

    def test_03_continuity_Euler_and_local_energy_balance(self):
        p = Parameters(.7, .8, .4, 3., .2)
        residuals = local_balance_residuals(.3, np.array([1.2, .6, -.4]), 1e-3, p)
        self.assertLess(max(abs(x) for x in residuals), 5e-11)

    def test_04_material_volume_number_is_conserved_not_density(self):
        for a in (1., 1.5, 2., 10.):
            t = (a-1)/DEFAULT.expansion
            self.assertAlmostEqual(background(t)[2]*a**3, DEFAULT.density, places=12)
        self.assertNotEqual(background(1/DEFAULT.expansion)[2], DEFAULT.density)

    def test_05_metric_determinant_inverse_give_the_quadratic_action(self):
        for t in (0., 7., 200.):
            x = np.array([1., -.3, .7])
            _, hubble, density, sound = background(t)
            velocity = hubble*x
            actual = acoustic_metric(t, x, normalized=False)
            density_inverse = np.sqrt(-np.linalg.det(actual))*np.linalg.inv(actual)
            expected = np.zeros((4, 4))
            expected[0, 0] = -density/sound**2
            expected[0, 1:] = expected[1:, 0] = -density*velocity/sound**2
            expected[1:, 1:] = density*(np.eye(3)-np.outer(velocity, velocity)/sound**2)
            np.testing.assert_allclose(density_inverse, expected, atol=3e-13)
            self.assertEqual(np.count_nonzero(np.linalg.eigvalsh(actual) < 0), 1)

    def test_06_full_coordinate_pullback_removes_the_flow(self):
        for t in (0., 31., 500.):
            a = background(t)[0]
            comoving = np.array([.3, -.4, .7])
            jacobian = np.eye(4)
            jacobian[1:, 0] = DEFAULT.expansion*comoving
            jacobian[1:, 1:] = a*np.eye(3)
            pulled = jacobian.T @ acoustic_metric(t, a*comoving) @ jacobian
            np.testing.assert_allclose(pulled, comoving_metric_jets(t)[0], atol=3e-15)

    def test_07_proper_time_and_scale_derivatives(self):
        for t in (0., 20., 500.):
            a, hubble, _, _ = background(t)
            tau_dot = derivative_5(acoustic_proper_time, t, .01)
            self.assertAlmostEqual(tau_dot, DEFAULT.c0*a**(-2.25), places=10)
            actual_hubble = .25*hubble/tau_dot
            self.assertAlmostEqual(actual_hubble, acoustic_hubble(t), places=12)
            h_prime = derivative_5(acoustic_hubble, t, .01)/tau_dot
            self.assertAlmostEqual(h_prime, 5*actual_hubble**2, places=12)

    def test_08_direct_connection_contraction_reproduces_all_tensor_components(self):
        for t in (0., 100., 500.):
            metric, first, second = comoving_metric_jets(t)
            actual = curvature_from_metric_jets(metric, first, second)
            h2 = acoustic_hubble(t)**2
            np.testing.assert_allclose(actual["orthonormal_einstein"],
                                       h2*np.array([3., -13., -13., -13.]), rtol=2e-14)
            self.assertAlmostEqual(actual["scalar"]/h2, 42., places=11)
            self.assertLess(np.max(abs(actual["ricci"]-actual["ricci"].T)), 1e-15)
        zero = curvature_from_metric_jets(*comoving_metric_jets(1., Parameters(expansion=0.)))
        np.testing.assert_array_equal(zero["einstein"], np.zeros((4, 4)))

    def test_09_no_Lambda_can_cancel_the_null_projection(self):
        at = 100.
        metric, first, second = comoving_metric_jets(at)
        tensor = curvature_from_metric_jets(metric, first, second)["einstein"]
        null = np.array([1/np.sqrt(-metric[0, 0]), 1/np.sqrt(metric[1, 1]), 0., 0.])
        self.assertAlmostEqual(null @ metric @ null, 0., places=14)
        for lam in (-.3, 0., .7):
            actual = null @ (tensor+lam*metric) @ null
            self.assertAlmostEqual(actual, -10*acoustic_hubble(at)**2, places=14)

    def test_10_Bianchi_identity_does_not_select_the_source(self):
        for t in (0., 100., 500.):
            a = background(t)[0]
            h = acoustic_hubble(t)
            rho_prime = derivative_5(lambda s: 3*acoustic_hubble(s)**2, t, .01)
            rho_prime /= DEFAULT.c0*a**(-2.25)
            self.assertAlmostEqual(rho_prime+3*h*(3*h*h-13*h*h), 0., places=15)

    def test_11_acoustic_wave_equation_matches_hydrodynamic_mode(self):
        k = .1
        for t in (0., 100., 500.):
            a, hubble, _, sound = background(t)
            tau_dot = DEFAULT.c0*a**(-2.25)
            scale = a**.25
            time_friction = 3*(.25*hubble)-(-2.25*hubble)
            self.assertAlmostEqual(time_friction, 3*hubble)
            self.assertAlmostEqual(tau_dot**2*k*k/scale**2, sound**2*k*k/a**2)

    def test_12_linearized_GP_matrix_has_Bogoliubov_frequency(self):
        for a in (1., 2., 4.):
            for k in (.1, .5, 2.):
                eta0 = DEFAULT.hbar*k/(2*DEFAULT.mass*DEFAULT.c0)
                matrix = DEFAULT.c0*k*np.array([[0., a**(-2)],
                                                [-a**(-3)-eta0**2*a**(-2), 0.]])
                actual = np.sort(abs(np.linalg.eigvals(matrix)))
                expected = np.sqrt(DEFAULT.c0**2*k*k*a**(-5) +
                                   DEFAULT.hbar**2*k**4/(4*DEFAULT.mass**2)*a**(-4))
                np.testing.assert_allclose(actual, expected, rtol=2e-15)

    def test_13_diluteness_improves_while_both_probe_windows_deteriorate(self):
        start, end = mode_window(1., .1), mode_window(4., .1)
        self.assertAlmostEqual(end["eta_quantum_pressure"]/start["eta_quantum_pressure"], 2.)
        self.assertAlmostEqual(end["H_over_hydro_frequency"]/start["H_over_hydro_frequency"], 8.)
        self.assertAlmostEqual(end["gas_parameter"]/start["gas_parameter"], 1/64)
        self.assertLess(start["gas_parameter"], 1e-6)

    def test_14_RK4_fourth_order_refinement(self):
        coarse, _ = evolve_mode(steps=256)
        middle, _ = evolve_mode(steps=512)
        fine, _ = evolve_mode(steps=1024)
        ratio = np.linalg.norm(coarse-middle)/np.linalg.norm(middle-fine)
        self.assertGreater(ratio, 14.)
        self.assertLess(ratio, 18.)

    def test_15_linear_symplectic_Wronskian_is_conserved(self):
        for quantum in (False, True):
            _, error = evolve_mode(steps=2048, include_quantum_pressure=quantum)
            self.assertLess(error, 2e-9)

    def test_16_quantum_pressure_changes_modes_above_numerical_error(self):
        hydro, _ = evolve_mode(steps=2048, include_quantum_pressure=False)
        full, _ = evolve_mode(steps=2048, include_quantum_pressure=True)
        finer, _ = evolve_mode(steps=4096, include_quantum_pressure=True)
        difference = np.linalg.norm(full-hydro)
        numerical = np.linalg.norm(full-finer)
        self.assertGreater(difference, 1e-3)
        self.assertLess(difference, .2)
        self.assertGreater(difference, 1000*numerical)


if __name__ == "__main__":
    main(__name__, "acoustic_geometry_dynamics_audit", report)
