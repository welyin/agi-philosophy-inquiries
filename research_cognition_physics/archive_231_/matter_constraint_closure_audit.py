"""Round 352: continuum scalar constraint brackets and common propagation.

Fourier quadrature probes continuum functional derivatives on smooth periodic
configurations. This does not construct a first-class finite lattice theory.
"""
import unittest
import numpy as np
from growing_stream_audit import main


def derivative(f):
    n = len(f)
    k = np.fft.fftfreq(n, d=1.0/n)
    shape = (n,) + (1,)*(np.ndim(f)-1)
    return np.fft.ifft(1j*k.reshape(shape)*np.fft.fft(f, axis=0), axis=0).real


def integral(f):
    return float(2*np.pi*np.mean(f))


def potential(phi, masses, quartic):
    norm = np.sum(phi**2, axis=1)
    return (0.5*np.sum(masses*phi**2, axis=1)+quartic*norm**2/4,
            masses*phi+quartic*norm[:, None]*phi)


def data(n=96):
    x = 2*np.pi*np.arange(n)/n
    phi = np.column_stack((np.sin(x)+0.2*np.cos(2*x),
                           0.3*np.cos(x)-0.1*np.sin(3*x)))
    pi = np.column_stack((0.4+0.3*np.sin(2*x), 0.2*np.cos(x)+0.3*np.sin(x)))
    lapse = 1.2+0.2*np.cos(x)
    other = 0.4+np.sin(x)+0.1*np.cos(2*x)
    return x, phi, pi, lapse, other


def hamiltonian(phi, pi, lapse, a, b, masses, quartic, volume=None, inverse=None):
    if volume is None:
        volume = np.ones(len(phi))
    if inverse is None:
        inverse = np.ones(len(phi))
    grad = derivative(phi)
    v, _ = potential(phi, masses, quartic)
    density = (np.einsum('ni,ij,nj->n', pi, a, pi)/(2*volume)
               +volume*inverse*np.einsum('ni,ij,nj->n', grad, b, grad)/2
               +volume*v)
    return integral(lapse*density)


def functional_derivatives(phi, pi, lapse, a, b, masses, quartic,
                           volume=None, inverse=None):
    if volume is None:
        volume = np.ones(len(phi))
    if inverse is None:
        inverse = np.ones(len(phi))
    _, dv = potential(phi, masses, quartic)
    df = (-derivative((lapse*volume*inverse)[:, None]*(derivative(phi)@b))
          +(lapse*volume)[:, None]*dv)
    dp = (lapse/volume)[:, None]*(pi@a)
    return df, dp


def bracket(phi, pi, lapse, other, a, b, masses, quartic,
            volume=None, inverse=None):
    fn, pn = functional_derivatives(phi, pi, lapse, a, b, masses, quartic, volume, inverse)
    fm, pm = functional_derivatives(phi, pi, other, a, b, masses, quartic, volume, inverse)
    return integral(np.sum(fn*pm-pn*fm, axis=1))


def expected(phi, pi, lapse, other, matrix, inverse=None):
    if inverse is None:
        inverse = np.ones(len(phi))
    v = inverse*(lapse*derivative(other)-other*derivative(lapse))
    return integral(v*np.einsum('ni,ij,nj->n', pi, matrix, derivative(phi)))


def witness(matrix, row, col, n=64):
    """Flat torus, N=1, M=sin x, pi_a=1, phi_b=sin x."""
    x = 2*np.pi*np.arange(n)/n
    fields = len(matrix)
    phi = np.zeros((n, fields))
    pi = np.zeros_like(phi)
    phi[:, col] = np.sin(x)
    pi[:, row] = 1
    return expected(phi, pi, np.ones(n), np.sin(x), matrix)


def speeds_squared(a, b):
    w, v = np.linalg.eigh(a)
    root = (v*np.sqrt(w))@v.T
    return np.linalg.eigvalsh(root@b@root)


def nonlinear_legendre(p, s, alpha):
    """Regular branch P(X)=X+alpha X^2, alpha*s<1."""
    if np.any(alpha*s >= 1) or alpha < 0:
        raise ValueError('This audit uses the monotone regular branch.')
    v=np.asarray(p,dtype=float).copy()
    for _ in range(30):
        v-=((1+alpha*(v*v-s))*v-p)/(1+alpha*(3*v*v-s))
    x=(v*v-s)/2
    px=1+2*alpha*x
    return p*v-x-alpha*x*x, v, px/2


def report():
    x, phi, pi, lapse, other = data()
    a = np.array([[2., .3], [.3, 1.]])
    b_good = np.linalg.inv(a)
    b_bad = np.array([[1., .2], [.2, 1.3]])
    masses = np.array([.4, 1.7])
    cases = []
    for name, b in (('common', b_good), ('different', b_bad)):
        value = bracket(phi, pi, lapse, other, a, b, masses, .8)
        rhs = expected(phi, pi, lapse, other, a@b)
        common = expected(phi, pi, lapse, other, np.eye(2))
        cases.append({'case':name, 'bracket':value, 'matrix_rhs':rhs,
                      'common_rhs':common, 'identity_error':abs(value-rhs),
                      'common_closure_residual':value-common,
                      'squared_speeds':speeds_squared(a,b).tolist()})
    matrix = a@b_bad-np.eye(2)
    witness_values = [[witness(matrix, i, j) for j in range(2)] for i in range(2)]
    masses_k = [np.array([.1, .9]), np.array([4., 7.])]
    arbitrary_potentials = [
        {'masses_squared':m.tolist(), 'quartic':q,
         'bracket':bracket(phi,pi,lapse,other,a,b_good,m,q)}
        for m in masses_k for q in (-3., 0., 2.)]
    return {'round':352,
            'scope':'Conditional continuum scalar ansatz, arbitrary lapse and off-shell common deformation algebra; not a quantum reconstruction, lattice constraint proof, or derivation of GR from FUCP.',
            'inputs':{'constant_positive_momentum_matrix':a.tolist(),
                      'common_gradient_matrix':b_good.tolist(),
                      'unequal_gradient_matrix':b_bad.tolist(),
                      'gravity_structure_coefficient':1.0,
                      'spatial_dimension_selected':False,
                      'potential_selected':False},
            'cases':cases,
            'entrywise_necessity_witness':{'pi_times_defect':(np.pi*matrix).tolist(),
                                          'quadrature':witness_values},
            'arbitrary_nonderivative_potentials':arbitrary_potentials,
            'homogeneous_probe':{'gradient':0, 'bracket':0,
                                 'distinguishes_speeds':False},
            'outside_quadratic_ansatz':{
                'lagrangian':'P(X)=X+alpha X^2',
                'alpha':.5, 'timelike_X':.4,
                'time_kinetic_coefficient':2.2,
                'spatial_kinetic_coefficient':1.4,
                'sound_speed_squared':1.4/2.2,
                'deformation_bracket_coefficient':1.,
                'closure_alone_implies_common_characteristics':False},
            'sources':['https://arxiv.org/html/gr-qc/0012089',
                       'https://arxiv.org/html/hep-th/9904176']}


class Checks(unittest.TestCase):
    def setUp(self):
        self.x, self.phi, self.pi, self.n, self.m = data()
        self.a = np.array([[2., .3], [.3, 1.]])
        self.b = np.array([[1., .2], [.2, 1.3]])
        self.mass = np.array([.4, 1.7])

    def test_01_derivative_and_integration_by_parts(self):
        self.assertLess(np.max(abs(derivative(np.sin(3*self.x))-3*np.cos(3*self.x))), 1e-12)
        self.assertAlmostEqual(integral(self.n*derivative(self.m)),
                               -integral(derivative(self.n)*self.m), places=12)

    def test_02_functional_directional_variation(self):
        u = np.column_stack((np.cos(2*self.x), np.sin(self.x)))
        v = np.column_stack((np.sin(3*self.x), np.cos(self.x)))
        df, dp = functional_derivatives(self.phi,self.pi,self.n,self.a,self.b,self.mass,.8)
        eps = 1e-5
        def h(sign):
            return hamiltonian(self.phi+sign*eps*u,self.pi+sign*eps*v,
                               self.n,self.a,self.b,self.mass,.8)
        self.assertAlmostEqual((h(1)-h(-1))/(2*eps),
                               integral(np.sum(df*u+dp*v,axis=1)), places=8)

    def test_03_continuum_bracket_identity(self):
        for n in (32,64,96):
            _,f,p,l,m=data(n)
            self.assertAlmostEqual(bracket(f,p,l,m,self.a,self.b,self.mass,.8),
                                   expected(f,p,l,m,self.a@self.b), places=11)

    def test_04_curved_volume_cancels_in_bracket(self):
        # A smooth diagonal 3-metric depending only on x:
        # qxx=e^(2f), qyy=qzz=e^(-f), hence sqrt(q)=1.
        inverse=np.exp(-.3*np.cos(self.x))
        volume=np.ones(len(self.x))
        self.assertAlmostEqual(bracket(self.phi,self.pi,self.n,self.m,
                                       self.a,self.b,self.mass,.8,volume,inverse),
                               expected(self.phi,self.pi,self.n,self.m,
                                        self.a@self.b,inverse), places=11)
        # Another diagonal metric, qxx=e^(2f), qyy=qzz=1.
        volume=np.exp(.15*np.cos(self.x))
        self.assertAlmostEqual(bracket(self.phi,self.pi,self.n,self.m,
                                       self.a,self.b,self.mass,.8,volume,inverse),
                               expected(self.phi,self.pi,self.n,self.m,
                                        self.a@self.b,inverse), places=11)

    def test_05_potential_and_interaction_cancel(self):
        values=report()['arbitrary_nonderivative_potentials']
        self.assertLess(max(v['bracket'] for v in values)-min(v['bracket'] for v in values),1e-12)

    def test_06_common_closure_matrix(self):
        for sigma in (.4,1.,3.):
            b=sigma*np.linalg.inv(self.a)
            value=bracket(self.phi,self.pi,self.n,self.m,self.a,b,self.mass,.8)
            rhs=expected(self.phi,self.pi,self.n,self.m,sigma*np.eye(2))
            self.assertAlmostEqual(value,rhs,places=11)

    def test_07_offdiagonal_necessity_witness(self):
        defect=self.a@self.b-np.eye(2)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(witness(defect,i,j),np.pi*defect[i,j],places=12)
        self.assertGreater(abs(witness(defect,0,1)),.1)

    def test_08_mixing_does_not_remove_unequal_characteristics(self):
        change=np.array([[1.,.3],[.1,1.4]])
        inv=np.linalg.inv(change)
        transformed_a=inv@self.a@inv.T
        transformed_b=change.T@self.b@change
        np.testing.assert_allclose(speeds_squared(transformed_a,transformed_b),
                                   speeds_squared(self.a,self.b),atol=1e-13)

    def test_09_canonical_normalization_for_common_cone(self):
        np.testing.assert_allclose(speeds_squared(self.a,np.linalg.inv(self.a)),[1.,1.],atol=1e-13)
        self.assertGreater(np.ptp(speeds_squared(self.a,self.b)),.2)

    def test_10_homogeneous_fields_hide_failure(self):
        phi=np.tile([.3,.7],(len(self.x),1))
        self.assertAlmostEqual(bracket(phi,self.pi,self.n,self.m,self.a,self.b,self.mass,.8),0,places=12)

    def test_11_projectable_lapses_hide_failure(self):
        self.assertAlmostEqual(bracket(self.phi,self.pi,np.ones(len(self.x)),
                                       2*np.ones(len(self.x)),self.a,self.b,self.mass,.8),0,places=12)

    def test_12_ultralocal_branch_is_different_algebra(self):
        b=np.zeros((2,2))
        value=bracket(self.phi,self.pi,self.n,self.m,self.a,b,self.mass,.8)
        self.assertAlmostEqual(value,0,places=12)
        self.assertGreater(abs(expected(self.phi,self.pi,self.n,self.m,np.eye(2))),.01)

    def test_13_random_positive_matrices_and_field_numbers(self):
        rng=np.random.default_rng(352)
        for count in (1,2,4):
            r=rng.normal(size=(count,count))
            a=r@r.T+np.eye(count)
            b=np.linalg.inv(a)
            np.testing.assert_allclose(speeds_squared(a,b),np.ones(count),atol=1e-12)
            for i in range(count):
                for j in range(count):
                    self.assertAlmostEqual(witness(a@b-np.eye(count),i,j),0,places=11)

    def test_14_closure_does_not_guarantee_potential_stability(self):
        f=np.ones((len(self.x),2))*10
        p=np.zeros_like(f)
        h=hamiltonian(f,p,np.ones(len(self.x)),self.a,np.linalg.inv(self.a),
                      np.zeros(2),-1)
        self.assertLess(h,-1e4)
        # Same matrix closure condition despite the unbounded-below potential.
        np.testing.assert_allclose(self.a@np.linalg.inv(self.a),np.eye(2),atol=1e-13)

    def test_15_nonlinear_legendre_and_bracket_coefficient(self):
        p=np.array([.3,1.,1.7])
        s=np.array([.1,.2,.3])
        alpha=.5
        h,v,hs=nonlinear_legendre(p,s,alpha)
        np.testing.assert_allclose(2*v*hs,p,atol=1e-13)
        e=1e-5
        hp=(nonlinear_legendre(p+e,s,alpha)[0]-nonlinear_legendre(p-e,s,alpha)[0])/(2*e)
        hs_fd=(nonlinear_legendre(p,s+e,alpha)[0]-nonlinear_legendre(p,s-e,alpha)[0])/(2*e)
        np.testing.assert_allclose(hp,v,atol=1e-9,rtol=1e-9)
        np.testing.assert_allclose(hs_fd,hs,atol=1e-9,rtol=1e-9)

    def test_16_closed_nonlinear_bracket_with_distinct_sound_cone(self):
        x=2*np.pi*np.arange(128)/128
        phi=.2*np.sin(x)
        p=1+.1*np.cos(x)
        n=1+.2*np.sin(2*x)
        m=.3+np.cos(x)
        grad=derivative(phi)
        _,v,hs=nonlinear_legendre(p,grad**2,.5)
        fn=-derivative(n*2*hs*grad)
        fm=-derivative(m*2*hs*grad)
        actual=integral(fn*m*v-n*v*fm)
        target=integral((n*derivative(m)-m*derivative(n))*p*grad)
        self.assertAlmostEqual(actual,target,places=11)
        c=report()['outside_quadratic_ansatz']
        self.assertGreater(c['time_kinetic_coefficient'],0)
        self.assertGreater(c['spatial_kinetic_coefficient'],0)
        self.assertLess(c['sound_speed_squared'],1.)


if __name__=='__main__':
    main(__name__,'matter_constraint_closure_audit',report)
