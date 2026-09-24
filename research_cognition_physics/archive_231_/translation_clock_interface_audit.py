"""Round 422: exact outgoing interfaces of a normal translation clock.

Only quadrature of characteristic solutions is discretized; P is never replaced
by a finite matrix.  General statements are proved in the companion note.
Python/NumPy only. --check writes nothing; --output creates a new result file.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
Y0 = -2 * math.pi
OBS: dict[str, object] = {}


def exp_h(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * t * values)) @ vectors.conj().T


def norm(a):
    return float(np.linalg.norm(a, 2))


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh((a-b + (a-b).conj().T)/2))) / 2)


def sinc(x):
    return float(np.sinc(x / math.pi))


def characteristic(w):
    """Integral of exp(-4iz) against the compact cosine-squared density."""
    return sinc(2*w) + (sinc(2*w + math.pi) + sinc(2*w - math.pi))/2


@lru_cache(maxsize=None)
def gauss(n):
    return np.polynomial.legendre.leggauss(n)


def packet_quadrature(w, center=Y0, n=96):
    x, weights = gauss(n)
    z = w*x/2
    density = 2/w * np.cos(math.pi*z/w)**2
    return center+z, weights*w/2, density


def characteristic_propagator(h, segments, y, t):
    """Exact ordered exponentials along one ray x=y+s for step potentials."""
    assert t >= 0
    end = y+t
    cuts = sorted({y, end} | {float(c) for a,b,_ in segments for c in (a,b) if y<c<end})
    out = np.eye(h.shape[0], dtype=complex)
    for left, right in zip(cuts[:-1], cuts[1:]):
        middle = (left+right)/2
        local = h.copy()
        for a,b,v in segments:
            if a < middle < b:
                local += v
        out = exp_h(local, right-left) @ out
    return out


def scattering_matrix(h, segments):
    a = min(item[0] for item in segments)
    b = max(item[1] for item in segments)
    s = characteristic_propagator(h, segments, a, b-a)
    return exp_h(h, -b) @ s @ exp_h(h, a)


def outgoing_propagator(h, m, y, t):
    return exp_h(h, t) @ exp_h(h, y) @ m @ exp_h(h, -y)


def channel_from_rays(rho, h, segments, t, w, center=Y0, reference_dim=1, n=96):
    y, quadrature, density = packet_quadrature(w, center, n)
    out = np.zeros_like(rho, dtype=complex)
    for yy, weight in zip(y, quadrature*density):
        u = characteristic_propagator(h, segments, float(yy), t)
        ur = np.kron(u, np.eye(reference_dim))
        out += weight * (ur @ rho @ ur.conj().T)
    return out


def pure_density(vector):
    vector = vector / np.linalg.norm(vector)
    return np.outer(vector, vector.conj())


def choi_state(d=2):
    return pure_density(np.eye(d).reshape(-1).astype(complex))


def target_dephased_x(rho, chi, t, reference_dim=1):
    x = np.kron(X, np.eye(reference_dim))
    z = np.kron(Z, np.eye(reference_dim))
    u = np.kron(exp_h(Z, t), np.eye(reference_dim))
    r = x @ rho @ x
    r = (1+chi)/2*r + (1-chi)/2*(z @ r @ z)
    return u @ r @ u.conj().T


EXAMPLE_SEGMENTS = [(0., math.pi, -Z + (I2+X)/2)]


class TranslationClockTests(unittest.TestCase):
    def test_01_characteristics_and_scattering_factorization(self):
        h = np.array([[.4,.1j,.2],[-.1j,-.7,.05],[.2,.05,.9]], complex)
        v1 = np.array([[.2,.3,0],[.3,0,.2j],[0,-.2j,-.1]], complex)
        v2 = np.array([[0,.1,.2j],[.1,-.3,.2],[-.2j,.2,.4]], complex)
        segments = [(.1,.7,v1),(.7,1.6,v2)]
        m = scattering_matrix(h, segments)
        errors = []
        for y in (-2.3,-1.7,-.4):
            for t in (4.2,6.1):
                u = characteristic_propagator(h, segments, y, t)
                errors.append(norm(u-outgoing_propagator(h,m,y,t)))
        y, t, s = -1.2, 1.6, 2.1
        direct = characteristic_propagator(h, segments, y, t+s)
        split = characteristic_propagator(h, segments, y+t, s) @ characteristic_propagator(h, segments, y, t)
        errors.append(norm(direct-split))
        # Local Schrödinger derivative at an endpoint strictly inside a segment.
        y, t, dt = -.2, .6, 1e-6
        plus = characteristic_propagator(h, segments, y, t+dt)
        minus = characteristic_propagator(h, segments, y, t-dt)
        u = characteristic_propagator(h, segments, y, t)
        derivative_error = norm((plus-minus)/(2*dt) + 1j*(h+v1) @ u)
        self.assertLess(max(errors), 2e-13)
        self.assertLess(derivative_error, 2e-9)
        OBS['characteristics'] = {'maximum_factorization_and_cocycle_error': max(errors), 'local_derivative_error': derivative_error}

    def test_02_commuting_gate_complete_release(self):
        h = np.diag([-.6,.2,.8]).astype(complex)
        k = np.diag([.3,-.5,.7]).astype(complex)
        segments = [(0.,1.2,k/1.2)]
        target = exp_h(k,1.)
        max_error = 0.
        for y in (-4.1,-2.2,-.1):
            for t in (6.,9.7):
                u = characteristic_propagator(h,segments,y,t)
                max_error = max(max_error,norm(u-exp_h(h,t)@target))
        # An arbitrary future operation need not commute with the old H.
        g = exp_h(np.array([[0,1,.2j],[1,.2,0],[-.2j,0,-.1]],complex),.7)
        t, delta = 6., .53
        for y in (-4.1,-2.2,-.1):
            first = characteristic_propagator(h,segments,y,t)
            later = characteristic_propagator(h,segments,y+t,delta) @ g @ first
            expected = exp_h(h,delta) @ g @ exp_h(h,t) @ target
            max_error = max(max_error,norm(later-expected))
        self.assertLess(max_error, 2e-13)
        OBS['commuting_release'] = {'dimension':3, 'all_ray_operator_error': max_error, 'future_noncommuting_operation_checked':True}

    def test_03_zero_drift_all_gate_positive_example(self):
        k = np.array([[.2,.3j,.4],[-.3j,-.1,.25],[.4,.25,.5]], complex)
        h = np.zeros((3,3),complex)
        segments = [(0.,.4,k*.75),(.4,1.4,k*.7)]
        target = exp_h(k,1.)
        error = max(norm(characteristic_propagator(h,segments,y,t)-target)
                    for y in (-2.4,-.8) for t in (4.1,8.))
        self.assertLess(error,2e-13)
        OBS['zero_drift'] = {'operator_error':error, 'non_diagonal_logarithm':True}

    def test_04_normal_clock_moments_and_characteristic(self):
        rows = []
        for w in (.1,.4,1.,2.):
            y, q, density = packet_quadrature(w,n=128)
            z = y-Y0
            phi = np.sqrt(2/w)*np.cos(math.pi*z/w)
            derivative = -np.sqrt(2/w)*(math.pi/w)*np.sin(math.pi*z/w)
            p_mean = np.sum(q * phi * (-1j*derivative))
            p_second = float(np.sum(q*derivative**2))
            chi = np.sum(q*density*np.exp(-4j*z))
            zz = np.linspace(-w/2,w/2,8193)
            dd = 2/w*np.cos(math.pi*zz/w)**2
            trap = np.trapezoid(dd*np.exp(-4j*zz),zz)
            self.assertAlmostEqual(float(np.sum(q*density)),1.,places=12)
            self.assertLess(abs(p_mean),2e-11)
            self.assertLess(abs(p_second/(math.pi**2/w**2)-1),3e-13)
            self.assertLess(abs(chi-characteristic(w)),3e-13)
            self.assertLess(abs(trap-chi),3e-12)
            rows.append({'width':w,'normalization':float(np.sum(q*density)), 'mean_P_abs':float(abs(p_mean)), 'P_second_moment':p_second,'chi4':float(chi.real),'independent_trapezoid_difference':float(abs(trap-chi))})
        OBS['clock_moments'] = rows

    def test_05_noncommuting_pulse_entire_channel_and_free_future(self):
        m = scattering_matrix(Z,EXAMPLE_SEGMENTS)
        self.assertLess(norm(m-X),3e-14)
        rho = choi_state()
        rows = []
        for w in (.2,.5,1.,2.):
            t = 3*math.pi+w/2+.7
            out = channel_from_rays(rho,Z,EXAMPLE_SEGMENTS,t,w,reference_dim=2)
            predicted = target_dephased_x(rho,characteristic(w),t,reference_dim=2)
            error = norm(out-predicted)
            delta = .91
            later = channel_from_rays(rho,Z,EXAMPLE_SEGMENTS,t+delta,w,reference_dim=2)
            free = np.kron(exp_h(Z,delta),I2)
            future_error = norm(later-free@out@free.conj().T)
            refined = channel_from_rays(rho,Z,EXAMPLE_SEGMENTS,t,w,reference_dim=2,n=192)
            eigen = np.linalg.eigvalsh(out)
            self.assertLess(error,3e-13)
            self.assertLess(future_error,3e-13)
            self.assertLess(norm(refined-out),3e-13)
            self.assertGreater(eigen[-2],1e-4)
            rows.append({'width':w,'complete_channel_choi_error':error,'future_free_channel_error':future_error,'quadrature_refinement_error':norm(refined-out),'choi_nonzero_eigenvalues':[float(eigen[-2]),float(eigen[-1])]})
        OBS['noncommuting_pulse'] = {'M_minus_X':norm(m-X),'rows':rows}

    def test_06_exact_diamond_error_and_unknown_references(self):
        rng = np.random.default_rng(422)
        rows=[]
        for w in (.2,.5,1.,2.):
            eps=(1-characteristic(w))/2
            t=3*math.pi+w/2+.3
            target=exp_h(Z,t)@X
            plus=pure_density(np.array([1.,1.],complex))
            actual=channel_from_rays(plus,Z,EXAMPLE_SEGMENTS,t,w)
            witness=trace_distance(actual,target@plus@target.conj().T)
            bell=choi_state()
            u2=np.kron(target,I2)
            bell_distance=trace_distance(channel_from_rays(bell,Z,EXAMPLE_SEGMENTS,t,w,reference_dim=2),u2@bell@u2.conj().T)
            max_random=0.
            for _ in range(8):
                a=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
                r=a@a.conj().T
                r/=np.trace(r)
                u3=np.kron(target,np.eye(3))
                actual=channel_from_rays(r,Z,EXAMPLE_SEGMENTS,t,w,reference_dim=3)
                d=trace_distance(actual,u3@r@u3.conj().T)
                self.assertLessEqual(d,eps+3e-13)
                max_random=max(max_random,d)
            self.assertLess(abs(witness-eps),3e-13)
            self.assertLess(abs(bell_distance-eps),3e-13)
            rows.append({'width':w,'half_diamond_analytic':eps,'plus_witness':witness,'bell_witness':bell_distance,'maximum_random_reference_distance':max_random})
        OBS['reference_error'] = rows

    def test_07_normal_mixed_clock_and_clock_phase(self):
        # A normal mixture of two wave packets has the mixture of their outgoing channels.
        r=choi_state()
        t=4*math.pi
        mixture=.37*channel_from_rays(r,Z,EXAMPLE_SEGMENTS,t,.6,center=Y0,reference_dim=2)
        mixture+=.63*channel_from_rays(r,Z,EXAMPLE_SEGMENTS,t,.9,center=Y0+.23,reference_dim=2)
        eigen=np.linalg.eigvalsh(mixture)
        self.assertGreater(eigen[-2],.02)
        # Changing the wave function's position-dependent phase does not change a partial trace over position.
        y,q,density=packet_quadrature(.8,n=64)
        weights=q*density
        vectors=[]
        bell=np.eye(2).reshape(-1).astype(complex)/np.sqrt(2)
        for yy,weight in zip(y,weights):
            u=np.kron(characteristic_propagator(Z,EXAMPLE_SEGMENTS,float(yy),t),I2)
            vectors.append(np.sqrt(weight)*(u@bell))
        amplitudes=np.asarray(vectors)
        phased=amplitudes*np.exp(1j*.7*y*y)[:,None]
        reduced=np.einsum('ki,kj->ij',amplitudes,amplitudes.conj())
        reduced_phased=np.einsum('ki,kj->ij',phased,phased.conj())
        error=norm(reduced-reduced_phased)
        self.assertLess(error,3e-14)
        OBS['normal_mixed_clock']={'choi_second_eigenvalue':float(eigen[-2]),'coherent_clock_phase_independence_error':error}

    def test_08_phase_ambiguity_requires_finite_trace_argument(self):
        # A finite cyclic shift looks like an energy ladder shift except at the wrap-around.
        d=7
        h=np.diag(np.arange(-3,4)).astype(complex)
        shift=np.roll(np.eye(d),1,axis=0).astype(complex)
        delta=shift.conj().T@h@shift-h
        self.assertLess(abs(np.trace(delta)),1e-13)
        self.assertLess(norm(delta-np.diag([1.]*6+[-6.])),1e-13)
        y=.17
        orbit=exp_h(h,y)@shift@exp_h(h,-y)
        phase_target=np.exp(-1j*y)*shift
        interior_error=norm((orbit-phase_target)[:,:-1])
        full_error=norm(orbit-phase_target)
        self.assertLess(interior_error,2e-14)
        self.assertGreater(full_error,1.)
        OBS['finite_trace_phase_check']={'dimension':d,'trace_of_energy_difference':float(np.trace(delta).real),'interior_phase_error':interior_error,'full_cyclic_shift_phase_error':full_error,'boundary_energy_change':-6.}

    def test_09_width_energy_tradeoff_is_only_this_packet_family(self):
        coefficient=1/3-2/math.pi**2
        limit=math.pi**2/3-2
        rows=[]
        errors=[]
        for w in (.5,.25,.125,.0625,.03125):
            eps=(1-characteristic(w))/2
            second=math.pi**2/w**2
            errors.append(abs(eps/w**2-coefficient))
            rows.append({'width':w,'half_diamond_error':eps,'P_second_moment':second,'epsilon_over_width_squared':eps/w**2,'epsilon_times_P_second':eps*second})
        self.assertTrue(all(b<a for a,b in zip(errors[:-1],errors[1:])))
        self.assertLess(errors[-1],2e-5)
        self.assertLess(abs(rows[-1]['epsilon_times_P_second']-limit),2e-4)
        OBS['packet_family_tradeoff']={'coefficient_limit':coefficient,'product_limit':limit,'rows':rows,'not_a_universal_optimality_bound':True}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.check and args.output:
        parser.error('--check and --output are exclusive')
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(TranslationClockTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    report={'round':422, 'baseline_round':420, 'status':'verified_conditional_result',
            'tests_run':result.testsRun, 'failures':len(result.failures), 'errors':len(result.errors),
            'python':platform.python_version(), 'numpy':np.__version__,
            'scope':dict(normal_clock_required=True,
                         infinite_clock_and_translation_are_model_input=True,
                         free_future_after_exit_proved=True,
                         noncommuting_exact_unitary_in_this_model=False,
                         semibounded_energy_required=False,
                         permanent_outputs_required=False,
                         finite_object_hardware_closure_proved=False,
                         three_dimensional_space_derived=False,
                         full_cognitive_countermodel_completed=False,
                         full_GR_goal_completed=False, phase_closure_triggered=False),
            'numerics':'Gaussian quadrature and independent trapezoid integration of exact characteristic solutions; no finite-matrix discretization of P and no plots.',
            'observations':OBS}
    if args.output:
        with args.output.open('x',encoding='utf-8') as handle:
            json.dump(report,handle,ensure_ascii=False,indent=2,allow_nan=False)
            handle.write('\n')
        print(f'Created {args.output}')
    else:
        print(f'{result.testsRun} checks passed; no files written.')


if __name__=='__main__':
    main()
