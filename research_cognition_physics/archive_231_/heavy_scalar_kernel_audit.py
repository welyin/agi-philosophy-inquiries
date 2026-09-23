"""Round 328: static heavy-scalar elimination and a controlled derivative expansion."""
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import gauss


def kernel(k, mass=8., order=None, omega=0.):
    if mass <= 0:
        raise ValueError('The expansion requires a positive heavy mass.')
    z = (np.asarray(k)**2-omega**2)/mass**2
    if order is None:
        return 1/(mass**2*(1+z))
    return sum((-z)**j for j in range(order+1))/mass**2


def source_modes(size=64):
    modes = np.zeros(size, dtype=complex)
    modes[0] = .6*size
    for n, value in ((1, .2), (2, .15j), (3, .1)):
        modes[n] = size*value
        modes[-n] = size*np.conj(value)
    return modes


def periodic_audit(mass=8., size=64):
    k = np.fft.fftfreq(size)*size
    modes = source_modes(size)
    source = np.fft.ifft(modes).real
    exact = np.fft.ifft(-kernel(k,mass)*modes).real
    identity = np.eye(size)
    operator = np.fft.ifft((mass**2+k**2)[:,None]*np.fft.fft(identity,axis=0),axis=0).real
    direct = np.linalg.solve(operator,-source)
    reduced_energy = .5*np.mean(source*exact)
    full_energy = (.5*exact@operator@exact+source@exact)/size
    rows=[]
    for order in range(4):
        approx = np.fft.ifft(-kernel(k,mass,order)*modes).real
        rows.append({'order':order,
                     'relative_field_error':float(np.linalg.norm(approx-exact)/np.linalg.norm(exact)),
                     'spectral_bound':float((3/mass)**(2*(order+1)))})
    return {'direct_solve_error':float(np.max(abs(direct-exact))),
            'stationary_equation_error':float(np.max(abs(operator@exact+source))),
            'full_minus_reduced_energy':float(full_energy-reduced_energy),
            'reduced_energy':float(reduced_energy),'expansion':rows}


def separated_sources(mass=8., distance=1., width=.2):
    # A 1-D infinite-line diagnostic, distinct from the band-limited periodic example.
    nodes, weights = gauss(80)
    bump = np.exp(-1/(1-nodes**2))
    measures = width*weights*bump
    x, y = width*nodes, distance+width*nodes
    green = np.exp(-mass*abs(x[:,None]-y[None,:]))/(2*mass)
    energy = -measures@green@measures
    bound = np.sum(measures)**2*np.exp(-mass*(distance-2*width))/(2*mass)
    # Finite-difference local operators have no overlapping source support here.
    grid=np.linspace(-1.,3.,1024,endpoint=False); dx=4/len(grid)
    def compact(center):
        coordinate=(grid-center)/width
        result=np.zeros_like(grid); inside=abs(coordinate)<1
        result[inside]=np.exp(-1/(1-coordinate[inside]**2))
        return result
    left,right=compact(0.),compact(distance)
    term=right/mass**2; local=term.copy()
    for _ in range(3):
        term=(np.roll(term,1)-2*term+np.roll(term,-1))/(dx*dx*mass**2)
        local+=term
    return {'exact_cross_energy':float(energy), 'absolute_tail_bound':float(bound),
            'finite_local_expansion_cross_energy':float(-dx*left@local)}


def report():
    return {'round':328,
            'scope':'Classical/tree-level elimination of the quadratic heavy scalar with prescribed small source J. Static positive operator and band-limited error bounds are exact in this linear sector. No loop determinant, generic real-time Green-function choice or all-amplitude exponential-mass elimination is included.',
            'periodic_static_solution':periodic_audit(),
            'disjoint_compact_sources':[{'mass':m,**separated_sources(m)} for m in (4.,8.,12.)],
            'outside_expansion_window': [{'order':n,'relative_error':float(abs(kernel(16.,8.,n)/kernel(16.,8.)-1))} for n in range(4)],
            'near_timelike_pole_order_2_error':float(abs(kernel(0.,8.,2,omega=7.92)/kernel(0.,8.,omega=7.92)-1)),
            'next_interface':'Restore the exponential mass functions and test amplitude control, positivity and the apparent negative-quartic instability.'}


class Checks(unittest.TestCase):
    def test_01_spectral_solution_matches_direct_matrix_solve(self):
        r=periodic_audit()
        self.assertLess(r['direct_solve_error'],2e-15)
        self.assertLess(r['stationary_equation_error'],2e-14)

    def test_02_stationary_energy_matches_eliminated_action_sign(self):
        r=periodic_audit()
        self.assertLess(abs(r['full_minus_reduced_energy']),2e-16)
        self.assertLess(r['reduced_energy'],0.)

    def test_03_spectral_relative_bound(self):
        for mass in (4.,8.,12.):
            for row in periodic_audit(mass)['expansion']:
                self.assertLessEqual(row['relative_field_error'],row['spectral_bound']+2e-15)

    def test_04_exact_remainder_for_a_single_mode(self):
        z=(3/8)**2
        for n in range(4):
            error=kernel(3.,8.,n)/kernel(3.,8.)-1
            self.assertAlmostEqual(float(error),-(-z)**(n+1),places=14)

    def test_05_high_spatial_frequency_series_fails(self):
        errors=[abs(kernel(16.,8.,n)/kernel(16.,8.)-1) for n in range(4)]
        np.testing.assert_allclose(errors,[4.,16.,64.,256.],atol=1e-12)

    def test_06_static_low_momentum_does_not_control_time_dependence(self):
        self.assertEqual(float(kernel(0.,8.,2)/kernel(0.,8.)-1),0.)
        error=abs(kernel(0.,8.,2,omega=7.92)/kernel(0.,8.,omega=7.92)-1)
        self.assertGreater(error,.9)

    def test_07_nonzero_separated_tail_obeys_absolute_bound(self):
        for mass in (4.,8.,12.):
            r=separated_sources(mass)
            self.assertGreater(abs(r['exact_cross_energy']),0.)
            self.assertLessEqual(abs(r['exact_cross_energy']),r['absolute_tail_bound'])
            self.assertEqual(r['finite_local_expansion_cross_energy'],0.)

    def test_08_separation_has_exact_exponential_dependence_in_one_dimension(self):
        first=separated_sources(8.,1.)['exact_cross_energy']
        second=separated_sources(8.,1.5)['exact_cross_energy']
        self.assertAlmostEqual(second/first,float(np.exp(-4.)),places=14)

    def test_09_induced_mixed_species_term_remains(self):
        j1,j2,mass=.2,.3,8.
        mixed=((j1+j2)**2-j1**2-j2**2)/(2*mass**2)
        self.assertAlmostEqual(mixed,j1*j2/mass**2,places=15)
        self.assertGreater(mixed,0.)

    def test_10_mass_gap_required_for_this_expansion(self):
        with self.assertRaises(ValueError):
            kernel(1.,0.,2)


if __name__=='__main__':
    main(__name__,'heavy_scalar_kernel_audit',report)
