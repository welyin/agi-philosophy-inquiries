"""Round 339: independent finite species, aggregate stability and common cones.

This audits the principal action at homogeneous timelike states of a stipulated
classical theory. No multi-species cosmological attraction is assumed or solved.
"""
import unittest
import numpy as np
from growing_stream_audit import main
from disformal_characteristics_audit import matrices as single_matrices


def state_from_q_r(q, r, velocity=1.):
    """Construct one local homogeneous state; r_i=D_i rho_i, q_i=D_i v^2."""
    q, r = np.asarray(q, dtype=float), np.asarray(r, dtype=float)
    if (q.ndim != 1 or len(q) == 0 or q.shape != r.shape or velocity == 0
            or np.any(q <= 0) or np.any(q >= 1) or np.any(r < 0)):
        raise ValueError('Use finite nonempty arrays, 0<q_i<1, r_i>=0 and v!=0.')
    coupling = q/velocity**2
    matter_velocity = np.sqrt(2*r*(1-q)**1.5/coupling)
    return float(velocity), matter_velocity, coupling


def state_from_fractions(q, fractions, velocity=1.):
    """Fractions share one total density, including the canonical phi energy."""
    q, fractions = np.asarray(q, dtype=float), np.asarray(fractions, dtype=float)
    total = np.sum(fractions)
    if q.shape != fractions.shape or np.any(fractions < 0) or total >= 1:
        raise ValueError('Nonnegative matter fractions must sum to less than one.')
    return state_from_q_r(q, q*fractions/(2*(1-total)), velocity)


def parameters(velocity, matter_velocity, coupling):
    matter_velocity, coupling = np.asarray(matter_velocity), np.asarray(coupling)
    q = coupling*velocity**2
    c = 1-q
    if np.any(c <= 0) or np.any(coupling < 0):
        raise ValueError('The chosen branch has D_i>=0 and C_i>0.')
    r = coupling*matter_velocity**2/(2*c**1.5)
    mixing = coupling*velocity*matter_velocity/c
    return q, c, r, mixing


def action_density(derivatives, coupling):
    """Full covariant local action, not an expansion.

    Derivative rows are (time,x,y,z), columns (phi,chi_1,...,chi_m).
    """
    derivatives = np.asarray(derivatives)
    coupling = np.asarray(coupling)
    phi = derivatives[:, 0]
    matter = derivatives[:, 1:]
    xx = -phi[0]**2+np.dot(phi[1:], phi[1:])
    yy = -matter[0]**2+np.sum(matter[1:]**2, axis=0)
    zz = -phi[0]*matter[0]+phi[1:]@matter[1:]
    root = np.sqrt(1+coupling*xx)
    return float(-xx/2+np.sum(-root*yy/2+coupling*zz**2/(2*root)))


def matrices(velocity, matter_velocity, coupling):
    q, c, r, mixing = parameters(velocity, matter_velocity, coupling)
    kinetic = np.diag(np.r_[1.+np.sum(r+3*r*q/c), 1/np.sqrt(c)])
    gradient = np.diag(np.r_[1.+np.sum(r*(2*q-1)), np.sqrt(c)])
    kinetic[0, 1:] = kinetic[1:, 0] = mixing/np.sqrt(c)
    gradient[0, 1:] = gradient[1:, 0] = mixing*np.sqrt(c)
    return kinetic, gradient


def simultaneous_basis(velocity, matter_velocity, coupling):
    _, _, _, mixing = parameters(velocity, matter_velocity, coupling)
    result = np.eye(1+len(matter_velocity))
    result[1:, 0] = -mixing
    return result


def analytic_squares(velocity, matter_velocity, coupling):
    _, c, r, _ = parameters(velocity, matter_velocity, coupling)
    return np.r_[(1-np.sum(r))/(1+np.sum(r/c)), c]


def spectrum(velocity, matter_velocity, coupling):
    kinetic, gradient = matrices(velocity, matter_velocity, coupling)
    inverse = np.linalg.inv(np.linalg.cholesky(kinetic))
    return np.linalg.eigvalsh(inverse@gradient@inverse.T)


def finite_hessian(velocity, matter_velocity, coupling, step=2e-4):
    """Differentiate all four spacetime components of the complete action."""
    size = 1+len(matter_velocity)
    point = np.zeros((4, size))
    point[0] = np.r_[velocity, matter_velocity]
    flat = point.ravel()
    shifts = np.eye(len(flat))*step
    def evaluate(x):
        return action_density(x.reshape(4, size), coupling)
    center = evaluate(flat)
    result = np.zeros((len(flat), len(flat)))
    for i in range(len(flat)):
        result[i, i] = (evaluate(flat+shifts[i])-2*center
                        +evaluate(flat-shifts[i]))/step**2
        for j in range(i):
            result[i, j] = result[j, i] = (
                evaluate(flat+shifts[i]+shifts[j])
                -evaluate(flat+shifts[i]-shifts[j])
                -evaluate(flat-shifts[i]+shifts[j])
                +evaluate(flat-shifts[i]-shifts[j]))/(4*step**2)
    return result


def hessian_audit(step):
    state = state_from_fractions([.12, .27, .4], [.08, .11, .16], velocity=.8)
    size = 4
    kinetic, gradient = matrices(*state)
    exact = np.zeros((4*size, 4*size))
    exact[:size, :size] = kinetic
    for axis in range(1, 4):
        exact[axis*size:(axis+1)*size, axis*size:(axis+1)*size] = -gradient
    actual = finite_hessian(*state, step=step)
    return {'step':step, 'all_components_max_error':float(np.max(abs(actual-exact))),
            'kinetic_max_error':float(np.max(abs(actual[:size,:size]-kinetic))),
            'gradient_max_error':float(np.max(abs(actual[size:2*size,size:2*size]+gradient))),
            'time_space_cross_max':float(np.max(abs(actual[:size,size:])))}


def row(state):
    velocity, matter_velocity, coupling = state
    q, c, r, mixing = parameters(*state)
    rho = matter_velocity**2/(2*c**1.5)
    fractions = rho/(.5*velocity**2+np.sum(rho))
    kinetic, gradient = matrices(*state)
    squares = analytic_squares(*state)
    return {'species_count':int(len(q)), 'D_i':coupling.tolist(),
            'q_i':q.tolist(), 'r_i':r.tolist(), 'total_R':float(np.sum(r)),
            'total_S':float(np.sum(r/c)), 'matter_fractions':fractions.tolist(),
            'total_matter_fraction':float(np.sum(fractions)),
            'phi_mode_speed_squared':float(squares[0]),
            'matter_mode_speeds_squared':squares[1:].tolist(),
            'maximum_squared_speed_deficit':float(np.max(1-squares)),
            'mixing_norm':float(np.linalg.norm(mixing)),
            'smallest_kinetic_eigenvalue':float(np.linalg.eigvalsh(kinetic)[0]),
            'smallest_gradient_eigenvalue':float(np.linalg.eigvalsh(gradient)[0]),
            'generalized_spectrum_error':float(np.max(abs(spectrum(*state)-np.sort(squares))))}


def report():
    q, r = np.array([.2, .4]), np.array([.6, .6])
    controlled = []
    for scale in (1., .1, .01, .001):
        state = state_from_fractions(scale*np.array([.12,.27,.4]), [.12,.18,.2])
        data = row(state)
        kinetic, gradient = matrices(*state)
        data['coefficient_operator_error'] = float(max(
            np.linalg.norm(kinetic-np.eye(4), 2),
            np.linalg.norm(gradient-np.eye(4), 2)))
        controlled.append(data)
    return {'round':339,
            'scope':'Finite independent massless matter scalars, each disformally coupled to the same canonical phi by a possibly different nonnegative D_i(phi), on homogeneous timelike states of a stipulated classical Einstein theory. Only the local principal action is derived; neither multi-species background attraction nor nonlinear, low-frequency or quantum stability is assumed or proved.',
            'literature_interface':{'source':'https://arxiv.org/html/1510.01650',
                'location':'Sections 2.1 and 3.1, equations 1-4 and 39-45. Their single-sector action is summed over independent matter fields; their conformal C is set to one. The multi-sector proof and counterexample are independently checked here.'},
            'analytic_summary':{
                'strict_gradient_stability':'R=sum_i D_i rho_i < 1, with C_i=1-D_i phi_dot^2>0 and D_i>=0',
                'phi_speed_squared':'(1-R)/(1+S), S=sum_i r_i/C_i',
                'matter_speeds_squared':'C_i',
                'common_gravity_cone_iff':'max_i q_i -> 0 and R -> 0, for fixed finite species on the positive branch',
                'sufficient_uniform_certificate':'If total matter energy fraction F<=1/2 and Q=max_i q_i<1, all scalar speed squares lie in [1-Q,1], R<=Q/2 and norm(A)<=Q/(1-Q)^(1/4).',
                'attraction_not_assumed':True},
            'full_action_hessian':[hessian_audit(4e-4), hessian_audit(2e-4)],
            'separately_stable_but_jointly_unstable':{
                'each_alone':[row(state_from_q_r([qi],[ri])) for qi,ri in zip(q,r)],
                'combined':row(state_from_q_r(q,r)),
                'comparison':'Same local phi velocity, couplings and corresponding matter velocities. Each separate theory adjusts its own Friedmann H; no common full background is claimed.'},
            'aggregate_stability_boundary':row(state_from_q_r([.2,.4],[.4,.6])),
            'controlled_distinct_coupling_limit':controlled,
            'metric_limit_without_aggregate_occupation_control':[
                row(state_from_q_r([epsilon,2*epsilon],[.2,.3]))
                for epsilon in (.01,.001,.0001)],
            'zero_occupation_with_distinct_matter_cones':row(state_from_q_r([.2,.4],[0.,0.])),
            'shared_scalar_cone_is_not_gravity_cone':row(
                state_from_fractions([.3]*4,[.125]*4)),
            'next_interface':'Derive closed multi-species background equations and prove or refute an invariant region controlling both Q and the total occupation F, without requiring all D_i to be identical. This report itself supplies only local principal conditions.'}


class Checks(unittest.TestCase):
    def test_01_full_action_has_correct_homogeneous_reduction(self):
        state = state_from_q_r([.2,.4,.6],[.1,.2,.3])
        v,w,d = state
        point=np.zeros((4,4)); point[0]=np.r_[v,w]
        expected=.5*v*v+np.sum(.5*w*w/np.sqrt(1-d*v*v))
        self.assertAlmostEqual(action_density(point,d),expected,places=14)

    def test_02_independent_full_spacetime_hessian_and_refinement(self):
        coarse,fine=hessian_audit(4e-4),hessian_audit(2e-4)
        self.assertLess(fine['all_components_max_error'],4e-7)
        self.assertLess(fine['all_components_max_error'],
                        coarse['all_components_max_error']/2)
        self.assertLess(fine['time_space_cross_max'],1e-9)

    def test_03_one_species_reproduces_frozen_round_337(self):
        for f,q in ((0.,.4),(.1,.4),(.5,.2),(.9,.4)):
            state=state_from_fractions([q],[f],velocity=np.sqrt(q/.6))
            k,g=matrices(*state); old_k,old_g=single_matrices(f,q)
            np.testing.assert_allclose(k,old_k,atol=2e-14)
            np.testing.assert_allclose(g,old_g,atol=2e-14)

    def test_04_one_triangular_basis_diagonalizes_both_forms(self):
        for q,r in (([.2,.4],[.6,.6]),([.1,.3,.5],[.02,.1,.15])):
            state=state_from_q_r(q,r)
            _,c,r,_=parameters(*state)
            k,g=matrices(*state); basis=simultaneous_basis(*state)
            np.testing.assert_allclose(basis.T@k@basis,
                np.diag(np.r_[1+sum(r/c),1/np.sqrt(c)]),atol=3e-14)
            np.testing.assert_allclose(basis.T@g@basis,
                np.diag(np.r_[1-sum(r),np.sqrt(c)]),atol=3e-14)
            self.assertAlmostEqual(np.linalg.det(basis),1.)

    def test_05_independent_generalized_eigensolver_for_distinct_species(self):
        rng=np.random.default_rng(339)
        for count in (1,2,4,7):
            for _ in range(20):
                state=state_from_q_r(rng.uniform(.01,.8,count),rng.uniform(0,.6,count))
                np.testing.assert_allclose(spectrum(*state),
                    np.sort(analytic_squares(*state)),atol=8e-14,rtol=8e-14)

    def test_06_characteristic_determinant_factors_for_all_modes(self):
        state=state_from_q_r([.1,.27,.45],[.08,.11,.09])
        k,g=matrices(*state); _,c,r,_=parameters(*state)
        speed=analytic_squares(*state)
        for lam in (-1.,0.,.5,1.2,2.):
            expected=(1+sum(r/c))*np.prod(1/np.sqrt(c))*np.prod(speed-lam)
            self.assertAlmostEqual(np.linalg.det(g-lam*k),expected,places=12)

    def test_07_separate_stability_does_not_survive_unbudgeted_combination(self):
        for q in (.2,.4):
            self.assertGreater(min(analytic_squares(*state_from_q_r([q],[.6]))),0.)
        state=state_from_q_r([.2,.4],[.6,.6])
        k,g=matrices(*state)
        self.assertGreater(np.linalg.eigvalsh(k)[0],0.)
        self.assertLess(np.linalg.eigvalsh(g)[0],0.)
        self.assertAlmostEqual(analytic_squares(*state)[0],-4/55,places=14)

    def test_08_total_stability_boundary_occurs_before_any_metric_degeneracy(self):
        state=state_from_q_r([.2,.4],[.4,.6])
        self.assertAlmostEqual(analytic_squares(*state)[0],0.,places=14)
        self.assertGreater(min(parameters(*state)[1]),.59)
        self.assertAlmostEqual(np.linalg.eigvalsh(matrices(*state)[1])[0],0.,places=14)

    def test_09_total_energy_certificate_is_uniform_in_species_count(self):
        rng=np.random.default_rng(9339)
        for count in (1,2,5,17):
            for total in (.01,.2,.5):
                weights=rng.random(count); fractions=total*weights/sum(weights)
                q=rng.uniform(.001,.7,count); largest=max(q)
                state=state_from_fractions(q,fractions)
                speeds=analytic_squares(*state)
                _,c,r,mix=parameters(*state)
                self.assertGreaterEqual(min(speeds),1-largest-2e-14)
                self.assertLessEqual(max(speeds),1+2e-14)
                self.assertLessEqual(sum(r),largest/2+2e-14)
                self.assertLessEqual(np.linalg.norm(mix),largest/(1-largest)**.25+2e-14)
                self.assertGreater(np.linalg.eigvalsh(matrices(*state)[1])[0],0.)

    def test_10_metric_convergence_without_occupation_control_is_insufficient(self):
        state=state_from_q_r([1e-5,2e-5],[.2,.3])
        speeds=analytic_squares(*state)
        self.assertLess(abs(speeds[0]-1/3),1e-5)
        self.assertGreater(min(speeds[1:]),.9999)
        self.assertGreater(row(state)['total_matter_fraction'],.999)

    def test_11_zero_occupation_does_not_erase_matter_characteristics(self):
        state=state_from_q_r([.2,.4],[0.,0.])
        np.testing.assert_allclose(analytic_squares(*state),[1.,.8,.6])
        np.testing.assert_allclose(matrices(*state)[0],
                                  np.diag([1.,1/np.sqrt(.8),1/np.sqrt(.6)]))

    def test_12_all_principal_coefficients_converge_in_controlled_unequal_limit(self):
        errors=[]
        for scale in (1.,.1,.01,.001):
            state=state_from_fractions(scale*np.array([.12,.27,.4]),[.12,.18,.2])
            k,g=matrices(*state)
            errors.append(max(np.linalg.norm(k-np.eye(4),2),
                              np.linalg.norm(g-np.eye(4),2)))
            self.assertLessEqual(max(1-analytic_squares(*state)),.4*scale+2e-14)
        self.assertLess(errors[-1],errors[0]/900)

    def test_13_equal_scalar_cones_do_not_equal_the_tensor_cone(self):
        state=state_from_fractions([.3]*4,[.125]*4)
        k,g=matrices(*state)
        np.testing.assert_allclose(g,.7*k,atol=2e-14)
        np.testing.assert_allclose(analytic_squares(*state),[.7]*5,atol=2e-14)

    def test_14_velocity_sign_changes_preserve_characteristics_and_positivity(self):
        state=state_from_q_r([.1,.3,.4],[.05,.1,.2])
        v,w,d=state
        signed=(-v,w*np.array([-1.,1.,-1.]),d)
        np.testing.assert_allclose(spectrum(*signed),spectrum(*state),atol=2e-14)
        self.assertGreater(np.linalg.eigvalsh(matrices(*signed)[0])[0],0.)


if __name__ == '__main__':
    main(__name__, 'multispecies_characteristics_audit', report)
