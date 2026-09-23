"""Round 337: mixed scalar characteristics, independent action Hessians, and limits."""
import unittest
import numpy as np
from growing_stream_audit import main
from disformal_occupied_attractor_audit import solve_phase


def background(fraction, q, coupling=.6):
    if not (0 <= fraction < 1 and 0 < q < 1 and coupling > 0):
        raise ValueError('Use the positive timelike homogeneous branch.')
    v = np.sqrt(q/coupling)
    rho_b = .5*v*v*fraction/(1-fraction)
    w = np.sqrt(2*rho_b*(1-q)**1.5)
    return v, w


def action_density(derivatives, coupling=.6):
    """Full scalar action at one inertial point; order is v,w,phi_x,chi_x."""
    v, w, px, cx = derivatives
    xx = -v*v+px*px
    yy = -w*w+cx*cx
    zz = -v*w+px*cx
    determinant_ratio = 1+coupling*xx
    return (-xx/2-np.sqrt(determinant_ratio)*yy/2
            +coupling*zz*zz/(2*np.sqrt(determinant_ratio)))


def matrices(fraction, q, coupling=.6):
    v, w = background(fraction, q, coupling)
    c = 1-q
    r = coupling*w*w/(2*c**1.5)
    mixing = coupling*v*w/c
    kinetic = np.array([[1+r+3*r*q/c, mixing/np.sqrt(c)],
                        [mixing/np.sqrt(c), 1/np.sqrt(c)]])
    gradient = np.array([[1+r*(2*q-1), mixing*np.sqrt(c)],
                         [mixing*np.sqrt(c), np.sqrt(c)]])
    return kinetic, gradient


def finite_hessian(point, h=2e-4):
    """Differentiate the full covariant action, independent of matrix formulas."""
    point = np.asarray(point, dtype=float)
    eye = np.eye(4)*h
    result = np.zeros((4, 4))
    center = action_density(point)
    for i in range(4):
        result[i, i] = (action_density(point+eye[i])-2*center
                        +action_density(point-eye[i]))/h**2
        for j in range(i):
            result[i, j] = result[j, i] = (
                action_density(point+eye[i]+eye[j])
                -action_density(point+eye[i]-eye[j])
                -action_density(point-eye[i]+eye[j])
                +action_density(point-eye[i]-eye[j]))/(4*h*h)
    return result


def spectrum(fraction, q):
    kinetic, gradient = matrices(fraction, q)
    chol = np.linalg.cholesky(kinetic)
    inverse = np.linalg.inv(chol)
    return np.linalg.eigvalsh(inverse@gradient@inverse.T)


def analytic_squares(fraction, q):
    r = fraction*q/(2*(1-fraction))
    c = 1-q
    return np.array([(1-r)/(1+r/c), c])


def simultaneous_basis(fraction, q):
    v, w = background(fraction, q)
    mixing = .6*v*w/(1-q)
    return np.array([[1., 0.], [-mixing, 1.]])


def row(fraction, q):
    kinetic, gradient = matrices(fraction, q)
    speeds = analytic_squares(fraction, q)
    return {'fraction':float(fraction), 'q':float(q),
            'r_D_rho_chi':float(fraction*q/(2*(1-fraction))),
            'phi_mode_speed_squared':float(speeds[0]),
            'theta_mode_speed_squared':float(speeds[1]),
            'smallest_kinetic_eigenvalue':float(np.linalg.eigvalsh(kinetic)[0]),
            'smallest_gradient_eigenvalue':float(np.linalg.eigvalsh(gradient)[0]),
            'generalized_spectrum_error':float(np.max(abs(spectrum(fraction,q)-np.sort(speeds))))}


def report():
    rows = []
    for fraction in (.01, .1, .5):
        path = solve_phase(fraction)
        for index in (0, 4096):
            data = row(*path[index, :2])
            data['initial_fraction'] = fraction
            data['log_scale_factor'] = index*100/4096
            rows.append(data)
    v,w = background(.1,.4)
    exact_k,exact_g = matrices(.1,.4)
    approximate = finite_hessian([v,w,0.,0.])
    return {'round':337,
            'scope':'Principal (two-derivative) scalar perturbations of the stipulated classical Einstein plus canonical phi plus disformal massless chi theory, on homogeneous timelike backgrounds. The gravitational scalar constraints do not change these principal coefficients; lower-derivative mixing, long-wavelength growth, nonlinear stability and quantum completion are not proved.',
            'literature_interface':{'source':'https://arxiv.org/html/1510.01650',
                                    'location':'Section 2.2.2 and 3.1, equations 39-45; set their conformal C to one and potentials to zero.'},
            'action_hessian_error':{'kinetic_max':float(np.max(abs(approximate[:2,:2]-exact_k))),
                                    'gradient_max':float(np.max(abs(-approximate[2:,2:]-exact_g)))},
            'attractor_trajectories':rows,
            'lorentzian_but_gradient_unstable':row(.9,.4),
            'metric_limit_without_fraction_control':[row(1/(1+q),q) for q in (.01,.001,.0001)],
            'analytic_region':{'fraction_upper':.5,'q_upper':.4,
                               'r_upper':.2,'all_scalar_speed_squared_lower':.6,
                               'maximum_speed_squared_deficit_bound':'q',
                               'mixing_upper_bound':'q/(1-q)^(1/4)'},
            'next_interface':'Determine whether attraction and hyperbolicity persist when the specified exponential coupling slope changes; common principal cones do not imply complete operation equivalence.'}


class Checks(unittest.TestCase):
    def test_01_full_covariant_action_reduces_to_homogeneous_lagrangian(self):
        for f,q in ((0.,.4),(.1,.4),(.5,.1)):
            v,w = background(f,q)
            self.assertAlmostEqual(action_density([v,w,0.,0.]),.5*v*v+.5*w*w/np.sqrt(1-q),places=14)

    def test_02_independent_action_hessians_and_refinement(self):
        v,w = background(.1,.4)
        k,g = matrices(.1,.4)
        exact = np.zeros((4,4)); exact[:2,:2]=k; exact[2:,2:]=-g
        coarse = np.max(abs(finite_hessian([v,w,0.,0.],4e-4)-exact))
        fine = np.max(abs(finite_hessian([v,w,0.,0.],2e-4)-exact))
        self.assertLess(fine,2e-7)
        self.assertLess(fine,coarse/2)

    def test_03_mixing_basis_simultaneously_diagonalizes_principal_action(self):
        for f,q in ((.1,.4),(.5,.4),(.9,.4),(.01,.001)):
            k,g = matrices(f,q); basis=simultaneous_basis(f,q)
            r=f*q/(2*(1-f)); c=1-q
            np.testing.assert_allclose(basis.T@k@basis,np.diag([1+r/c,1/np.sqrt(c)]),atol=2e-14)
            np.testing.assert_allclose(basis.T@g@basis,np.diag([1-r,np.sqrt(c)]),atol=2e-14)

    def test_04_generalized_eigenproblem_matches_known_speeds(self):
        for f in (0.,.01,.1,.5,.9):
            for q in (.001,.1,.4,.7):
                np.testing.assert_allclose(spectrum(f,q),np.sort(analytic_squares(f,q)),atol=3e-14,rtol=2e-14)

    def test_05_fraction_region_has_ordered_positive_subluminal_cones(self):
        for f in np.linspace(0,.5,31):
            for q in np.linspace(.0001,.4,31):
                phi,theta = analytic_squares(f,q)
                self.assertGreaterEqual(theta,.6)
                self.assertGreaterEqual(phi,theta-1e-14)
                self.assertLessEqual(phi,1+1e-14)
                self.assertGreater(np.linalg.eigvalsh(matrices(f,q)[0])[0],0.)

    def test_06_zero_matter_preserves_two_distinct_cones(self):
        np.testing.assert_allclose(analytic_squares(0.,.4),[1.,.6])

    def test_07_energy_balance_point_has_degenerate_cones_without_singular_basis(self):
        k,g = matrices(.5,.4)
        np.testing.assert_allclose(g,.6*k,atol=2e-15)
        self.assertAlmostEqual(np.linalg.det(simultaneous_basis(.5,.4)),1.)

    def test_08_positive_time_kinetic_and_lorentzian_metrics_allow_gradient_instability(self):
        k,g=matrices(.9,.4)
        self.assertGreater(np.linalg.eigvalsh(k)[0],0.)
        self.assertLess(np.linalg.eigvalsh(g)[0],0.)
        np.testing.assert_allclose(analytic_squares(.9,.4),[-.2,.6],atol=1e-14)

    def test_09_stability_boundary_is_not_metric_degeneracy(self):
        phi,theta=analytic_squares(5/6,.4)
        self.assertAlmostEqual(phi,0.,places=14)
        self.assertAlmostEqual(theta,.6)

    def test_10_occupied_evolution_keeps_stable_characteristics(self):
        for f in (.01,.1,.5):
            path=solve_phase(f)
            for fraction,q,_ in path[::32]:
                phi,theta=analytic_squares(fraction,q)
                self.assertGreaterEqual(phi,theta-2e-14)
                self.assertGreaterEqual(theta,.6-2e-14)
                self.assertLessEqual(1-phi,q+2e-14)

    def test_11_mixing_and_all_principal_coefficients_vanish_in_controlled_limit(self):
        for f in (.01,.1,.5):
            errors=[]
            for q in (.01,.001,.0001):
                k,g=matrices(f,q)
                errors.append(max(np.max(abs(k-np.eye(2))),np.max(abs(g-np.eye(2)))))
                mix=abs(simultaneous_basis(f,q)[1,0])
                self.assertLessEqual(mix,q/(1-q)**.25+1e-15)
            self.assertLess(errors[-1],errors[0]/90)

    def test_12_q_alone_does_not_control_other_characteristic(self):
        q=.0001; f=1/(1+q)
        phi,theta=analytic_squares(f,q)
        self.assertLess(phi,.34)
        self.assertGreater(theta,.999)


if __name__=='__main__':
    main(__name__,'disformal_characteristics_audit',report)
