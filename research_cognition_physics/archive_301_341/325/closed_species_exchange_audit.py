"""Round 325: closed scalar/matter exchange, including homogeneous Einstein backreaction."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import ETA, integrate

ALPHA = np.array([.3, .8])
MASSES = np.array([1., 1.3])
MU = .4


def potential(q, alpha=ALPHA):
    q = np.asarray(q)
    mass2 = MASSES**2*np.exp(2*alpha*q[..., :1])
    return .5*MU**2*q[..., 0]**2+.5*np.sum(mass2*q[..., 1:]**2, axis=-1)


def gradient(q, alpha=ALPHA):
    mass2 = MASSES**2*np.exp(2*alpha*q[0])
    return np.r_[MU**2*q[0]+np.sum(alpha*mass2*q[1:]**2), mass2*q[1:]]


def densities(q, v):
    mass2 = MASSES**2*np.exp(2*ALPHA*q[..., :1])
    sector = np.concatenate([(.5*v[..., :1]**2+.5*MU**2*q[..., :1]**2),
                             .5*v[..., 1:]**2+.5*mass2*q[..., 1:]**2], axis=-1)
    return sector, .5*np.sum(v*v, axis=-1)-potential(q)


def stress(q, dq, sector=None):
    if sector is None:
        chosen, local_potential = dq, potential(q)
    elif sector == 0:
        chosen, local_potential = dq[:1], .5*MU**2*q[0]**2
    else:
        chosen = dq[sector:sector+1]
        local_potential = .5*MASSES[sector-1]**2*np.exp(2*ALPHA[sector-1]*q[0])*q[sector]**2
    return chosen.T@chosen-ETA*(.5*np.einsum('ia,ab,ib', chosen, ETA, chosen)+local_potential)


def jet_check(on_shell=False, sector=None):
    rng = np.random.default_rng(325)
    q0 = np.array([.2, .4, .35])
    dq = rng.normal(size=(3, 4))*.2
    raw = rng.normal(size=(3, 4, 4))*.1
    hessian = (raw+raw.swapaxes(1, 2))/2
    if on_shell:
        hessian[:, 0, 0] = np.trace(hessian[:, 1:, 1:], axis1=1, axis2=2)-gradient(q0)
    def value(x):
        q = q0+dq@x+.5*np.einsum('a,iab,b->i', x, hessian, x)
        return stress(q, dq+np.einsum('iab,b->ia', hessian, x), sector)
    h = 1e-4
    divergence = np.zeros(4)
    for a in range(4):
        offset = np.eye(4)[a]*h
        derivative = (-value(2*offset)+8*value(offset)-8*value(-offset)+value(-2*offset))/(12*h)
        divergence += ETA[a, a]*derivative[a]
    if sector is None:
        field_equation = np.einsum('ab,iab->i', ETA, hessian)-gradient(q0)
        expected = field_equation@dq
    else:
        sources = ALPHA*MASSES**2*np.exp(2*ALPHA*q0[0])*q0[1:]**2
        expected = (np.sum(sources) if sector == 0 else -sources[sector-1])*dq[0]
    return float(np.max(abs(divergence-expected))), float(np.linalg.norm(divergence))


@lru_cache(maxsize=16)
def solve(steps=1024, gravity=True, reciprocal=True):
    q0 = np.array([.2, .4, .35]); v0 = np.array([.1, .2, -.15])
    rho0 = np.sum(densities(q0, v0)[0])
    h0 = np.sqrt(rho0/3) if gravity else 0.
    initial = np.r_[q0, v0, 0., h0, 0.]
    def rhs(y):
        q, v = y[:3], y[3:6]
        force = gradient(q)
        if not reciprocal:
            force[0] = MU**2*q[0]
        hubble = y[7]
        pressure = densities(q, v)[1]
        return np.r_[v, -3*hubble*v-force, hubble,
                     -.5*np.dot(v, v) if gravity else 0.,
                     3*hubble*np.exp(3*y[6])*pressure]
    return integrate(rhs, initial, 6., steps)


def audit(steps=1024, gravity=True, reciprocal=True):
    data = solve(steps, gravity, reciprocal)
    q, v = data[:, :3], data[:, 3:6]
    sectors, pressure = densities(q, v)
    rho = np.sum(sectors, axis=1)
    ratio = MASSES[0]/MASSES[1]*np.exp((ALPHA[0]-ALPHA[1])*q[:, 0])
    constraint = 3*data[:, 7]**2-rho
    budget = np.exp(3*data[:, 6])*rho+data[:, 8]-rho[0]
    return {'energy_density_initial': float(rho[0]), 'energy_density_final': float(rho[-1]),
            'flat_energy_drift': float(np.max(abs(rho-rho[0]))) if not gravity else None,
            'Einstein_constraint_drift': float(np.max(abs(constraint-constraint[0]))) if gravity else None,
            'pressure_work_budget_error': float(np.max(abs(budget))),
            'mass_ratio_min': float(np.min(ratio)), 'mass_ratio_max': float(np.max(ratio)),
            'sector_energy_changes': (sectors[-1]-sectors[0]).tolist(),
            'scale_factor_final': float(np.exp(data[-1, 6])), 'H_final': float(data[-1, 7])}


def report():
    return {'round': 325,
            'scope': 'An added classical covariant Einstein plus three canonical scalar action with masses m_i(phi)=m_i0 exp(alpha_i phi), alpha=(.3,.8). The conserved total stress includes reciprocal scalar sources. A homogeneous finite-time solution includes metric backreaction; no global conserved cosmological energy, full nonlinear stability, microscopic origin or uniqueness of this action is claimed.',
            'off_shell_Ward_error': jet_check()[0],
            'on_shell_total_divergence': jet_check(True)[1],
            'sector_exchange_errors': [jet_check(True, i)[0] for i in range(3)],
            'flat_closed_system': audit(gravity=False),
            'joint_Einstein_matter_solution': audit(),
            'omitted_reciprocal_source_control': audit(reciprocal=False),
            'constraint_convergence': [{'steps': n, 'error': audit(n)['Einstein_constraint_drift']} for n in (128, 256, 512, 1024)],
            'next_interface': 'Total conservation does not set equal scalar charges. Audit the distinct spin-2 soft consistency condition and its assumptions.'}


class Checks(unittest.TestCase):
    def test_01_action_potential_derivative(self):
        q = np.array([.2, .4, -.3]); h = 1e-5
        numerical = [(potential(q+h*d)-potential(q-h*d))/(2*h) for d in np.eye(3)]
        np.testing.assert_allclose(numerical, gradient(q), atol=2e-11)

    def test_02_off_shell_covariant_Ward_identity(self):
        error, norm = jet_check()
        self.assertLess(error, 2e-10)
        self.assertGreater(norm, .01)

    def test_03_on_shell_total_and_sector_exchange(self):
        self.assertLess(jet_check(True)[1], 2e-10)
        for i in range(3):
            error, norm = jet_check(True, i)
            self.assertLess(error, 2e-10)
            self.assertGreater(norm, .001)

    def test_04_flat_closed_energy(self):
        self.assertLess(audit(gravity=False)['flat_energy_drift'], 2e-9)
        self.assertGreater(abs(audit(gravity=False)['sector_energy_changes'][0]), .005)

    def test_05_full_Einstein_constraint(self):
        self.assertLess(audit()['Einstein_constraint_drift'], 2e-9)
        self.assertGreater(audit()['scale_factor_final'], 2.)

    def test_06_pressure_work_not_constant_comoving_energy(self):
        self.assertLess(audit()['pressure_work_budget_error'], 2e-9)
        data = solve()
        self.assertGreater(abs(data[-1, 8]), .005)

    def test_07_mass_ratios_really_change(self):
        for gravity in (False, True):
            a = audit(gravity=gravity)
            self.assertGreater(a['mass_ratio_max']-a['mass_ratio_min'], .05)

    def test_08_missing_backreaction_is_detected(self):
        a = audit(reciprocal=False)
        self.assertGreater(a['Einstein_constraint_drift'], .005)
        self.assertGreater(a['pressure_work_budget_error'], .005)

    def test_09_constraint_converges_at_fourth_order(self):
        errors = [audit(n)['Einstein_constraint_drift'] for n in (128, 256, 512)]
        self.assertGreater(errors[0]/errors[1], 12)
        self.assertGreater(errors[1]/errors[2], 12)

    def test_10_positive_potential_and_finite_branch(self):
        for gravity in (False, True):
            data = solve(gravity=gravity)
            self.assertTrue(np.all(np.isfinite(data)))
            self.assertGreaterEqual(float(np.min(potential(data[:, :3]))), 0.)
        self.assertGreater(float(np.min(solve()[:, 7])), 0.)


if __name__ == '__main__':
    main(__name__, 'closed_species_exchange_audit', report)
