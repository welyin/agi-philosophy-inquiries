"""Round 403: passive dependence versus arbitrary intervening controls.

An exactly soluble constant outside control removes detuning. Drift and control
norms are uniformly bounded; the required duration and integrated action grow.
Zero-error equivalence and resource-dependent bounds are proved in the note.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('passive_control_locality_audit_results.json')
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
PS = (I, X, Y, Z)


def kron(*factors):
    out = np.ones((1, 1), complex)
    for a in factors:
        out = np.kron(out, a)
    return out


def norm(a):
    return float(np.linalg.norm(a, 2))


def unitary(h, t):
    eig, vec = np.linalg.eigh(h)
    return (vec*np.exp(-1j*t*eig)) @ vec.conj().T


def td(a, b):
    d = a-b
    return float(np.abs(np.linalg.eigvalsh((d+d.conj().T)/2)).sum()/2)


def images(protocol, qubits):
    answer = []
    for p in PS:
        a = kron(p, np.eye(2**(qubits-1)))
        for operators in reversed(protocol):
            a = sum(k.conj().T@a@k for k in operators)
        answer.append(a)
    return answer


def choi(protocol, qubits):
    din = 2**qubits
    return sum(np.kron(p, a.T) for p,a in zip(PS, images(protocol, qubits)))/(2*din)


def receive(rho, environment_dim, reference_dim=1):
    tensor = rho.reshape(2, environment_dim, reference_dim, 2, environment_dim, reference_dim)
    return np.einsum('aerbes->arbs', tensor).reshape(2*reference_dim, 2*reference_dim)


def signal(u):
    a = np.array([1., 0.])
    outputs = []
    for sign in (1., -1.):
        b = np.array([1., sign])/np.sqrt(2)
        state = u@np.kron(a, b)
        outputs.append(receive(np.outer(state, state.conj()), 2))
    return td(*outputs)


def control_case(g):
    omega = math.sqrt(1-g*g)
    h = omega*kron(I, Z)+g*kron(X, X)
    control = -omega*kron(I, Z)
    t = math.pi/(4*g)
    ideal = unitary(g*kron(X, X), t)
    controlled = unitary(h+control, t)
    eta = math.sqrt(2*(1-omega))
    residue = [norm(g*kron(X, X)@kron(a, I)-kron(a, I)@ (g*kron(X, X))) for a in PS[1:]]
    bare = unitary(h, t)
    signal_gap_lower_bound = (signal(controlled)-signal(bare))/2
    integral_sine = 2*math.floor(t/math.pi)+1-math.cos(t % math.pi)
    row = dict(g=g, omega=omega, drift_operator_norm=norm(h),
               h_squared_error=norm(h@h-np.eye(4)),
               control_operator_norm=norm(control),
               controlled_generator_error=norm(h+control-g*kron(X, X)),
               control_commutators_on_receiver=[norm(control@kron(a,I)-kron(a,I)@control) for a in PS],
               passive_all_time_half_diamond_upper_bound=eta,
               passive_canonical_dependency_upper_bound=2*eta,
               duration=t, integrated_control_norm=omega*t,
               controlled_unitary_error=norm(controlled-ideal),
               controlled_signal=signal(controlled),
               uncontrolled_signal_same_time=signal(bare),
               all_time_uncontrolled_signal_upper_bound=g,
               control_vs_bare_channel_lower_bound=signal_gap_lower_bound,
               uniform_channel_gap_lower_bound=(1-g)/2,
               control_vs_bare_choi_distance=td(choi([[controlled]],2),choi([[bare]],2)),
               continuous_control_resource_bound=min(1.,2*eta*omega*integral_sine),
               generator_residual_paulis=residue,
               round399_uncapped_receiver_bound=float(.5*t*sum(residue)))
    return row, h, eta


def hybrid_protocol(h, times, final_time):
    p = kron(I, X)
    rot = kron(I, unitary(.3*Y+.2*Z, .8))
    damp = [kron(I, np.diag([1., np.sqrt(.37)])),
            kron(I, np.array([[0., np.sqrt(.63)], [0., 0.]]))]
    reset = [kron(I, np.array([[1., 0.], [0., 0.]])),
             kron(I, np.array([[0., 1.], [0., 0.]]))]
    choices = ([p], damp, [rot], reset)
    protocol = []
    previous = 0.
    for j, time in enumerate(times):
        protocol.extend(([unitary(h,time-previous)], choices[j % len(choices)]))
        previous = time
    protocol.append([unitary(h,final_time-previous)])
    return protocol


def exact_protected_model():
    h = .4*kron(X,X,I)+.2*kron(Z,I,I)+.7*kron(I,X,Z)+.1*kron(I,I,X)
    algebra = [kron(a,b,I) for a in PS for b in (I,X)]
    def project(a):
        return sum(np.trace(b.conj().T@a)/8*b for b in algebra)
    invariant = max(norm((h@a-a@h)-project(h@a-a@h)) for a in algebra)
    cross = .7*kron(I,X,Z)
    cross_error = max(norm(cross@a-a@cross) for a in algebra)
    rot = [kron(I,I,unitary(.3*X+.6*Y,1.1))]
    damp = [kron(I,I,np.diag([1.,np.sqrt(.2)])),
            kron(I,I,np.array([[0.,np.sqrt(.8)],[0.,0.]]))]
    reset = [kron(I,I,np.array([[1.,0.],[0.,0.]])),
             kron(I,I,np.array([[0.,1.],[0.,0.]]))]
    cases = []
    for t in (.2, 1., 3.7):
        protocol = [[unitary(h,.17*t)],rot,[unitary(h,.23*t)],damp,
                    [unitary(h,.37*t)],reset,[unitary(h,.23*t)]]
        original = [[unitary(h,t)]]
        a, b = choi(protocol,3), choi(original,3)
        cases.append(dict(time=t, full_channel_difference=norm(a-b),
                           largest_adjoint_difference=max(norm(c-d) for c,d in zip(images(protocol,3),images(original,3))),
                           controlled_choi_min_eigenvalue=float(np.linalg.eigvalsh(a)[0]),
                           controlled_choi_trace_error=float(abs(np.trace(a)-1))))
    evolved = images([[unitary(h,1.)]],3)[3]
    a_basis = [kron(p,I,I) for p in PS]
    a_projection = sum(np.trace(p.conj().T@evolved)/8*p for p in a_basis)
    return dict(invariant_algebra_residual=invariant, cross_interaction_norm=norm(cross),
                cross_commutator_on_invariant_algebra=cross_error,
                receiver_dependence_beyond_A=norm(evolved-a_projection), protocols=cases)


@lru_cache(None)
def report():
    rng = np.random.default_rng(403)
    a = rng.normal(size=(12,5))+1j*rng.normal(size=(12,5))
    rho = a@a.conj().T
    rho /= np.trace(rho)
    base_reference = receive(rho, 2, 3)
    cases, passive, hybrids = [], [], []
    for g in (.08, .035, .015):
        row, h, eta = control_case(g)
        cases.append(row)
        for t in (0., math.pi/4, math.pi/2, 3., row['duration']):
            u = unitary(h,t)
            analytic = np.cos(t)*np.eye(4)-1j*np.sin(t)*h
            comparator = kron(I,unitary(Z,t))
            v = kron(u,np.eye(3))
            out = receive(v@rho@v.conj().T,2,3)
            passive.append(dict(g=g,time=t,analytic_propagator_error=norm(u-analytic),
                                unitary_local_comparator_error=norm(u-comparator),
                                unitary_error_formula=eta*abs(math.sin(t)),
                                mixed_reference_output_distance=td(out,base_reference),
                                complete_channel_choi_distance=td(choi([[u]],2),choi([[np.eye(4)]],2)),
                                uniform_channel_upper_bound=eta,
                                signal=signal(u), signal_formula=abs(g*math.sin(2*t))))
        times = (.11,.29,.56,.79)
        final_time = 1.1
        protocol = hybrid_protocol(h,times,final_time)
        cp = choi(protocol,2)
        bare = choi([[unitary(h,final_time)]],2)
        bound = 2*eta*sum(abs(math.sin(final_time-time)) for time in times)
        hybrids.append(dict(g=g,steps=len(times),
                            half_diamond_upper_bound=min(1.,bound),
                            full_choi_distance=td(cp,bare),
                            choi_min_eigenvalue=float(np.linalg.eigvalsh(cp)[0]),
                            choi_trace_error=float(abs(np.trace(cp)-1))))
    return dict(round=403,
                scope='Fixed finite closed drift, full input/reference domain and declared outside controls. Exact passive locality equals intervention immunity; uniformly small passive influence alone gives no control-independent error-only guarantee.',
                controlled_cases=cases, passive_all_input_checks=passive,
                finite_protocol_hybrid_checks=hybrids, exact_protection=exact_protected_model(),
                exact_zero_equivalence_proved=True,
                uniform_passive_smallness_sufficient_for_unrestricted_control=False,
                drift_norm_uniformly_bounded=True,
                propagation_time_uniform_in_small_coupling=False,
                control_amplitude_uniform_in_small_coupling=True,
                integrated_control_action_uniform_in_small_coupling=False,
                constant_control_exactly_solved=True,
                outside_control_cost_ignored=False,
                full_autonomous_controller_derived=False,
                spatial_dimension_generated=False, full_cognition_to_gr_refuted=False)


class Audit(unittest.TestCase):
    def test_bounded_drift_analytic_propagator_and_full_reference_bound(self):
        for row in report()['controlled_cases']:
            self.assertAlmostEqual(row['drift_operator_norm'],1.)
            self.assertLess(row['h_squared_error'],1e-12)
        for row in report()['passive_all_input_checks']:
            self.assertLess(row['analytic_propagator_error'],1e-12)
            self.assertAlmostEqual(row['unitary_local_comparator_error'],row['unitary_error_formula'])
            self.assertLessEqual(row['mixed_reference_output_distance'],row['unitary_error_formula']+1e-12)
            self.assertLessEqual(row['complete_channel_choi_distance'],row['unitary_error_formula']+1e-12)

    def test_uncontrolled_signal_formula_is_valid_at_all_sampled_times(self):
        for row in report()['passive_all_input_checks']:
            self.assertAlmostEqual(row['signal'],row['signal_formula'])

    def test_constant_control_only_acts_outside_receiver(self):
        for row in report()['controlled_cases']:
            self.assertLess(row['controlled_generator_error'],1e-12)
            self.assertLessEqual(row['control_operator_norm'],1.)
            self.assertLess(max(row['control_commutators_on_receiver']),1e-12)
            self.assertLess(row['controlled_unitary_error'],1e-12)

    def test_continuous_control_bound_and_channel_gap_witness(self):
        for row in report()['controlled_cases']:
            bound = row['continuous_control_resource_bound']
            self.assertLessEqual(row['control_vs_bare_choi_distance'],bound+1e-12)
            self.assertLessEqual(row['control_vs_bare_channel_lower_bound'],bound+1e-12)
            self.assertGreaterEqual(row['control_vs_bare_channel_lower_bound']+1e-12,row['uniform_channel_gap_lower_bound'])

    def test_controlled_signal_gap_and_resource_accounting(self):
        rows = report()['controlled_cases']
        for row in rows:
            self.assertAlmostEqual(row['controlled_signal'],1.)
            self.assertLessEqual(row['uncontrolled_signal_same_time'],row['g']+1e-12)
            self.assertAlmostEqual(row['duration']*row['g'],math.pi/4)
            self.assertAlmostEqual(row['integrated_control_norm'],row['omega']*row['duration'])
        self.assertLess(rows[-1]['passive_all_time_half_diamond_upper_bound'],rows[0]['passive_all_time_half_diamond_upper_bound'])
        self.assertGreater(rows[-1]['integrated_control_norm'],rows[0]['integrated_control_norm'])
        self.assertGreater(rows[-1]['duration'],rows[0]['duration'])

    def test_finite_nonunitary_protocol_hybrid_bound(self):
        for row in report()['finite_protocol_hybrid_checks']:
            self.assertLessEqual(row['full_choi_distance'],row['half_diamond_upper_bound']+1e-12)
            self.assertGreaterEqual(row['choi_min_eigenvalue'],-1e-12)
            self.assertLess(row['choi_trace_error'],1e-12)

    def test_exact_immunity_with_nonzero_cross_interaction(self):
        row = report()['exact_protection']
        self.assertGreater(row['cross_interaction_norm'],.5)
        self.assertGreater(row['receiver_dependence_beyond_A'],.1)
        self.assertLess(row['invariant_algebra_residual'],1e-12)
        self.assertLess(row['cross_commutator_on_invariant_algebra'],1e-12)
        for case in row['protocols']:
            self.assertLess(case['full_channel_difference'],1e-12)
            self.assertLess(case['largest_adjoint_difference'],1e-12)
            self.assertGreaterEqual(case['controlled_choi_min_eigenvalue'],-1e-12)
            self.assertLess(case['controlled_choi_trace_error'],1e-12)

    def test_generator_residual_does_not_certify_false_isolation(self):
        for row in report()['controlled_cases']:
            self.assertAlmostEqual(row['generator_residual_paulis'][0],0.)
            for residual in row['generator_residual_paulis'][1:]:
                self.assertAlmostEqual(residual,2*row['g'])
            self.assertAlmostEqual(row['round399_uncapped_receiver_bound'],math.pi/2)
            self.assertFalse(report()['uniform_passive_smallness_sufficient_for_unrestricted_control'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checks.testsRun, failures=len(checks.failures), errors=len(checks.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
