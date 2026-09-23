"""Round 341: perturbative scale windows on the frozen round-336 branch.

M_D = D**(-1/4) is an operator power-counting scale, not a proved UV cutoff.
The contact amplitude is about a constant-field Minkowski reference; its use
for a slowly varying, weak-gradient background is an additional approximation.
"""
from functools import lru_cache
import itertools
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate
from disformal_occupied_attractor_audit import BETA, D0, Q0, phase_rhs, solve_phase

ETA = np.diag([-1., 1., 1., 1.])
FIXED_GROWTH_LOWER = 3/(2*np.sqrt(2))-1


def dot(p, q):
    return float(np.asarray(p)@ETA@np.asarray(q))


def lagrangian_exact(phi_gradient, chi_gradient, coupling):
    x = dot(phi_gradient, phi_gradient)
    y = dot(chi_gradient, chi_gradient)
    z = dot(phi_gradient, chi_gradient)
    root = np.sqrt(1+coupling*x)
    return -.5*x-.5*root*y+.5*coupling*z*z/root


def quartic(phi_gradient, chi_gradient, coupling):
    x = dot(phi_gradient, phi_gradient)
    y = dot(chi_gradient, chi_gradient)
    z = dot(phi_gradient, chi_gradient)
    return coupling*(.5*z*z-.25*x*y)


def scattering_momenta(energy, cosine):
    sine = np.sqrt(1-cosine*cosine)
    # All incoming: phi(p1), chi(p2), phi(p3), chi(p4).
    return np.array([[energy, 0., 0., energy],
                     [energy, 0., 0., -energy],
                     [-energy, -energy*sine, 0., -energy*cosine],
                     [-energy, energy*sine, 0., energy*cosine]])


def mandelstam(p):
    return tuple(-dot(p[0]+p[j], p[0]+p[j]) for j in (1, 2, 3))


def contact_contraction(p, coupling):
    return coupling*(dot(p[0], p[1])*dot(p[2], p[3])
                     +dot(p[0], p[3])*dot(p[2], p[1])
                     -dot(p[0], p[2])*dot(p[1], p[3]))


def contact_polynomial_coefficient(p, coupling):
    # Independent coefficient of z1*z2*z3*z4 in the local quartic density.
    # Inclusion-exclusion is exact here since every surviving monomial has
    # total degree four and must contain each of the four independent labels.
    total = 0.
    for bits in itertools.product((0, 1), repeat=4):
        phi = bits[0]*p[0]+bits[2]*p[2]
        chi = bits[1]*p[1]+bits[3]*p[3]
        total += (-1)**(4-sum(bits))*quartic(phi, chi, coupling)
    return total


def augmented_rhs(y):
    f, q, logh, phi = y
    return np.r_[phase_rhs(y[:3]), np.sqrt(6*(1-f))]


@lru_cache(maxsize=12)
def background(fraction=.1, end=100., steps=8192):
    values = integrate(augmented_rhs, np.array([fraction, Q0, 0., 0.]), end, steps)
    n = np.linspace(0., end, steps+1)
    logh = -np.log(3*np.sqrt(1-fraction))+values[:, 2]
    logd = np.log(D0)+BETA*values[:, 3]
    logm = -.25*logd
    # Independent recovery of D from q = D*6*H**2*(1-f), all in log space.
    recovered_logd = np.log(values[:, 1])-np.log(6.)-2*logh-np.log1p(-values[:, 0])
    return n, values, logh, logm, recovered_logd


def log_rate_bounds(f, q):
    logm_prime = -1.5*np.sqrt(1-f)
    logh_prime = -3+1.5*q*f
    return -1-logm_prime, logh_prime-logm_prime


def first_crossing(n, values, threshold=0.):
    candidates = np.flatnonzero(values >= threshold)
    if not len(candidates):
        return None
    i = int(candidates[0])
    if i == 0:
        return float(n[0])
    return float(n[i-1]+(n[i]-n[i-1])*(threshold-values[i-1])/(values[i]-values[i-1]))


def window_protocol(fraction=.1, start=50., end=150., steps=12288):
    n, y, logh, logm, _ = background(fraction, end, steps)
    start_index = int(round(start/end*steps))
    if abs(n[start_index]-start) > 1e-10:
        raise ValueError('Choose a grid containing the protocol starting time.')
    # A fixed comoving mode chosen to have E = 0.01 M_D at N=start.
    logk = start+logm[start_index]+np.log(.01)
    log_fixed_over_m = logk-n-logm
    log_fixed_over_h = logk-n-logh
    sl = slice(start_index, None)
    exit_small = first_crossing(n[sl], log_fixed_over_m[sl], np.log(.1))
    reaches_m = first_crossing(n[sl], log_fixed_over_m[sl], 0.)
    # Rechosen local geometric-mean probe: H/E = E/M_D = sqrt(H/M_D).
    return {'start_N':start, 'end_N':end, 'log_fixed_comoving_k':float(logk),
            'q_at_start':float(y[start_index, 1]),
            'initial_E_over_M_D':float(np.exp(log_fixed_over_m[start_index])),
            'initial_H_over_E':float(np.exp(-log_fixed_over_h[start_index])),
            'fixed_mode_reaches_E_over_M_D_0_1_at_N':exit_small,
            'fixed_mode_reaches_E_over_M_D_1_at_N':reaches_m,
            'end_log_E_fixed_over_M_D':float(log_fixed_over_m[-1]),
            'end_log_E_fixed_over_H':float(log_fixed_over_h[-1]),
            'end_reselected_H_over_E_equals_E_over_M_D':float(np.exp(.5*(logh[-1]-logm[-1]))),
            'end_q':float(y[-1, 1])}


def report():
    n, y, logh, logm, recovered_logd = background()
    rows = []
    for target in (0., 10., 50., 100.):
        i = int(round(target/100*8192))
        rows.append({'N':float(n[i]), 'fraction':float(y[i, 0]), 'q':float(y[i, 1]),
                     'log_H':float(logh[i]), 'log_M_D':float(logm[i]),
                     'log_H_over_M_D':float(logh[i]-logm[i]),
                     'd_log_fixed_k_phys_over_M_D_dN':float(log_rate_bounds(y[i, 0], y[i, 1])[0]),
                     'reselected_H_over_E_equals_E_over_M_D':float(np.exp(.5*(logh[i]-logm[i])))})
    amplitudes = []
    for cosine in (-.8, 0., .7):
        p = scattering_momenta(.2, cosine)
        s, t, u = mandelstam(p)
        amplitudes.append({'cos_theta':cosine, 'energy':.2,
                           'direct_quartic_coefficient':contact_polynomial_coefficient(p, D0),
                           'mandelstam_amplitude':float(-D0*s*u/2),
                           's_plus_t_plus_u':float(s+t+u)})
    return {'round':341,
            'scope':'Frozen round-336 classical homogeneous positive-occupation branch; local dimension-eight scalar contact amplitude is evaluated on a constant-field Minkowski reference. M_D=D^(-1/4) is a Wilson-coefficient power-counting scale, not a proved cutoff. Applying this count to late-time local probes assumes weak background gradients, slow variation, small probe energy/backreaction, and no lower omitted-operator scale.',
            'analytic_bounds':{'fixed_mode_log_ratio_growth_lower':float(FIXED_GROWTH_LOWER),
                               'H_over_M_D_log_rate_upper':-1.2,
                               'physical_coefficient_drift_rate_over_H_upper':6.,
                               'proper_time_curvature_scale_order':'H; |dot H|/H^2 <= 3',
                               'fixed_mode_result':'Any fixed nonzero comoving k eventually leaves k/a << M_D.',
                               'reselected_window_result':'For every sufficiently late N, E=sqrt(H*M_D) obeys H/E=E/M_D -> 0, with a different comoving k(N).'},
            'contact_amplitude':{'convention':'S=1+i(2pi)^4 delta^4(sum p) M; signature (-+++), all incoming phi1 chi2 phi3 chi4. Disformal contact only; gravitational exchange is excluded.',
                                 'formula':'M=D*(s*s+u*u-t*t)/4=-D*s*u/2; generic-angle M scales as D*E^4.',
                                 'rows':amplitudes},
            'scale_history':rows,
            'log_D_two_reconstructions_max_difference':float(np.max(abs(-4*logm-recovered_logd))),
            'fixed_vs_reselected_protocol':window_protocol(),
            'limitations':['No exact strong-coupling threshold or UV completion is established.',
                           'No Minkowski S-matrix is asserted for the full time-dependent occupied background.',
                           'The contact calculation alone neither computes loop corrections nor ensures absence of all lower scales.',
                           'The existence of a scale interval does not construct its preparation or detector resource budget.',
                           'No quantum instability, universal equivalence principle, or cognitive origin of the action is inferred.'],
            'next_interface':'Audit background-canonically-normalized interaction vertices and omitted operator scales before a quantum radiative-stability claim; preserve the distinction between a fixed comoving mode and separately prepared late-time probes.'}


class Checks(unittest.TestCase):
    def test_01_complete_disformal_density_has_the_independent_quartic_expansion(self):
        p = np.array([.3, .1, -.04, .02]); c = np.array([.2, -.03, .06, .02])
        base = -.5*(dot(p, p)+dot(c, c))
        errors = []
        for coupling in (.04, .02, .01):
            errors.append(abs(lagrangian_exact(p, c, coupling)-base-quartic(p, c, coupling)))
        self.assertGreater(errors[0]/errors[1], 3.98)
        self.assertGreater(errors[1]/errors[2], 3.98)

    def test_02_on_shell_conservation_and_mandelstam_identity(self):
        for energy in (.05, .2, 2.):
            for cosine in (-.8, 0., .7):
                p = scattering_momenta(energy, cosine)
                np.testing.assert_allclose(np.sum(p, axis=0), 0., atol=1e-14)
                self.assertLess(max(abs(dot(v, v)) for v in p), 1e-14)
                self.assertAlmostEqual(sum(mandelstam(p)), 0., places=13)

    def test_03_lagrangian_coefficient_and_tensor_vertex_agree(self):
        for cosine in (-.8, 0., .7):
            p = scattering_momenta(.2, cosine)
            self.assertAlmostEqual(contact_polynomial_coefficient(p, D0),
                                   contact_contraction(p, D0), places=15)

    def test_04_mandelstam_amplitude_and_fourth_power_energy_scaling(self):
        for cosine in (-.8, 0., .7):
            p = scattering_momenta(.2, cosine); s, t, u = mandelstam(p)
            amp = contact_contraction(p, D0)
            self.assertAlmostEqual(amp, D0*(s*s+u*u-t*t)/4, places=15)
            self.assertAlmostEqual(amp, -D0*s*u/2, places=15)
            self.assertAlmostEqual(contact_contraction(2*p, D0)/amp, 16., places=13)

    def test_05_augmented_flow_reuses_the_frozen_background(self):
        _, y, _, _, _ = background(end=10., steps=1024)
        old = solve_phase(.1, end=10., steps=1024)
        np.testing.assert_array_equal(y[:, :3], old)

    def test_06_phi_integral_and_constraint_recover_the_same_coupling(self):
        _, _, _, logm, logd_constraint = background()
        self.assertLess(np.max(abs(-4*logm-logd_constraint)), 3e-8)

    def test_07_uniform_rate_bounds_cover_the_invariant_rectangle(self):
        for f in np.linspace(0., .5, 21):
            for q in np.linspace(0., .4, 21):
                fixed, curvature = log_rate_bounds(f, q)
                self.assertGreaterEqual(fixed, FIXED_GROWTH_LOWER-1e-14)
                self.assertLessEqual(curvature, -1.2+1e-14)
                self.assertLessEqual(6*np.sqrt(1-f), 6.)

    def test_08_integrated_bounds_hold_on_independent_histories(self):
        for fraction in (.01, .1, .5):
            n, _, logh, logm, _ = background(fraction, end=100., steps=4096)
            fixed_change = -n-(logm-logm[0])
            curvature_change = (logh-logm)-(logh[0]-logm[0])
            self.assertTrue(np.all(fixed_change >= FIXED_GROWTH_LOWER*n-1e-10))
            self.assertTrue(np.all(curvature_change <= -1.2*n+1e-10))

    def test_09_fixed_mode_leaves_the_power_counting_window_while_staying_subhorizon(self):
        row = window_protocol()
        self.assertLess(row['q_at_start'], .001)
        self.assertLess(row['initial_H_over_E'], 1e-20)
        self.assertGreater(row['fixed_mode_reaches_E_over_M_D_0_1_at_N'], 50.)
        self.assertLess(row['fixed_mode_reaches_E_over_M_D_1_at_N'], 65.)
        self.assertGreater(row['end_log_E_fixed_over_M_D'], 0.)
        self.assertGreater(row['end_log_E_fixed_over_H'], 100.)

    def test_10_reselected_window_is_wide_and_slowly_varying_at_late_times(self):
        n, y, logh, logm, _ = background()
        late = n >= 50.
        ratio = np.exp(.5*(logh[late]-logm[late]))
        self.assertLess(np.max(ratio), 1e-15)
        self.assertTrue(np.all(np.diff(ratio) < 0.))
        # |dot log D|/E <= 6H/E, an adiabatic background-coefficient check.
        self.assertLess(np.max(6*ratio), 1e-14)
        self.assertLess(np.max(y[late, 1]), .001)

    def test_11_window_crossing_is_resolved_under_step_refinement(self):
        fine = window_protocol(); coarse = window_protocol(steps=6144)
        key = 'fixed_mode_reaches_E_over_M_D_1_at_N'
        self.assertLess(abs(fine[key]-coarse[key]), 2e-6)

    def test_12_initial_planck_normalization_does_not_fake_a_large_window(self):
        _, _, logh, logm, _ = background()
        initial = np.exp(logh[0]-logm[0])
        self.assertGreater(initial, .1)
        self.assertLess(initial, .5)
        # Early-time weak-background and hierarchy assumptions must not be
        # silently inferred from the late-time theorem.
        self.assertEqual(Q0, .4)


if __name__ == '__main__':
    main(__name__, 'disformal_validity_window_audit', report)
