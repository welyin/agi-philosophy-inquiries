"""Round 454: an enlarged common relation supports an exact static interaction.

Six raw qubits, seven positive unit SWAP couplings. A pre-existing singlet
between the two gauge factors is an explicit input resource, not free joining.
"""
import argparse
from functools import lru_cache
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np

import exchange_relation_audit as old
import relational_subject_composition_audit as previous

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'shared_relation_revival_audit_results.json'
OBS = {}
T = np.pi/(2*np.sqrt(2))
DELTA = 1/32


def short(x):
    return float(f'{float(x):.12g}')


@lru_cache(None)
def model():
    local, p3, iv, v3, w, p9, cross, comms, aa, cg = previous.model()
    singlet = np.array([0., 1., -1., 0.])/np.sqrt(2)
    v = w @ np.kron(singlet[:, None], np.eye(4))
    s = cross[8]
    edges = [(0,1), (0,2), (1,2), (3,4), (3,5), (4,5)]
    h0 = sum(previous.swap(6, a, b) for a, b in edges)
    b = (s @ v[:, 3]-v[:, 3]/3)*3/(2*np.sqrt(2))
    q = np.column_stack((v, b))
    h = h0+s
    m = np.diag([-1., 1., 1., 1/3, 17/3])
    m[3,4] = m[4,3] = 2*np.sqrt(2)/3
    x = old.operators()[4]
    effect = (np.eye(64)+np.kron(np.eye(8), x))/2
    return v, b, q, h0, s, h, m, effect, w, p9/9


def gate():
    return np.diag([np.exp(1j*T), np.exp(-1j*T),
                    np.exp(-1j*T), -np.exp(-3j*T)])


def analytic_columns(t):
    v, b, *_ = model()
    out = v.astype(complex).copy()
    out[:,0] *= np.exp(1j*t)
    out[:,1:3] *= np.exp(-1j*t)
    angle = 2*np.sqrt(2)*t
    survival = np.exp(-3j*t)*(np.cos(angle)+1j*(2*np.sqrt(2)/3)*np.sin(angle))
    transfer = -1j*np.exp(-3j*t)*np.sin(angle)/3
    out[:,3] = survival*v[:,3]+transfer*b
    return out


class Audit(unittest.TestCase):
    def test_01_common_singlet_and_invariant_five_state_interface(self):
        v, b, q, h0, s, h, m, effect, w, p = model()
        np.testing.assert_allclose(q.T@q, np.eye(5), atol=2e-15)
        np.testing.assert_allclose(h0@v, 0, atol=2e-15)
        np.testing.assert_allclose(h0@b, 6*b, atol=2e-15)
        np.testing.assert_allclose(h@q, q@m, atol=3e-15)
        for mu in range(3):
            total = sum(previous.site_pauli(6, a, mu) for a in range(6))
            np.testing.assert_allclose(total@q, 0, atol=3e-15)
        # Direct raw matrix multiplication checks the analytic restricted
        # matrix independently of the spectral exponential below.
        restricted = q.T@s@q
        expected = np.diag([-1.,1.,1.,1/3,-1/3])
        expected[3,4] = expected[4,3] = 2*np.sqrt(2)/3
        np.testing.assert_allclose(restricted, expected, atol=2e-15)
        # Independent exact integer certificate, with no eigensolver.
        iv = previous.model()[2]
        nn = np.column_stack([np.kron(iv[:,a],iv[:,2+c])-np.kron(iv[:,2+a],iv[:,c])
                              for a in range(2) for c in range(2)])
        np.testing.assert_array_equal(nn.T@nn,np.diag([8,24,24,72]))
        n = nn[:,3]
        bn = (3*s-np.eye(64,dtype=np.int64))@n
        np.testing.assert_array_equal(h@nn[:,:3],nn[:,:3]@np.diag([-1,1,1]))
        np.testing.assert_array_equal(3*h@n,n+bn)
        np.testing.assert_array_equal(3*h@bn,17*bn+8*n)
        np.testing.assert_array_equal(h0@bn,6*bn)
        self.assertEqual(int(n@bn),0)
        self.assertEqual(int(bn@bn),576)
        np.testing.assert_allclose(bn/24,b,atol=2e-15)
        self.assertGreater(np.linalg.norm((np.eye(64)-v@v.T)@h@v), .9)
        OBS['enlarged_interface'] = dict(raw_dimension=64, original_logic_dimension=4,
            common_invariant_dimension=5, additional_joint_mode=1,
            all_time_original_code_invariance=False,
            seven_raw_positive_equal_edges=[[0,1],[0,2],[1,2],[3,4],[3,5],[4,5],[2,5]],
            common_gauge_singlet_prepared=True,
            restricted_hamiltonian_diagonal=['-1','1','1','1/3','17/3'],
            restricted_hamiltonian_off_diagonal='H[3,4]=H[4,3]=2*sqrt(2)/3',
            total_SU2_invariance_checked=True, exact_integer_intertwiner_certificate=True,
            integer_logic_column_norms_squared=[8,24,24,72], integer_bridge_norm_squared=576)

    def test_02_exact_solution_and_unknown_input_gate(self):
        v, b, q, h0, s, h, m, effect, *_ = model()
        errors = []
        for t in [0., .17, T/2, T, T+DELTA, 2*T]:
            actual = old.evolve(h,t)@v
            error = np.linalg.norm(actual-analytic_columns(t), 2)
            self.assertLess(error, 9e-15)
            errors.append(short(error))
        u = old.evolve(h,T)
        np.testing.assert_allclose(u@v, v@gate(), atol=4e-15)
        # A gate-level identity proves all inputs and all references, including
        # complex amplitudes; this separate reference input catches indexing.
        rng = np.random.default_rng(454)
        state = rng.normal(size=12)+1j*rng.normal(size=12)
        state /= np.linalg.norm(state)
        actual = np.kron(u@v, np.eye(3))@state
        ideal = np.kron(v@gate(), np.eye(3))@state
        self.assertLess(np.linalg.norm(actual-ideal), 6e-15)
        OBS['exact_endpoint'] = dict(time='pi/(2*sqrt(2))', time_numeric=short(T),
            diagonal_gate=['exp(i*T)','exp(-i*T)','exp(-i*T)','-exp(-3i*T)'],
            arbitrary_logic_and_reference_operator_identity=True,
            analytic_evolution_max_numeric_error=short(max(errors)),
            reference_dimension_checked=3,
            input_dependent_control_or_postselection=False)

    def test_03_entanglement_and_joint_mode_occupation(self):
        v, b, q, h0, s, h, m, effect, *_ = model()
        diagonal = np.diag(gate())
        ratio = diagonal[0]*diagonal[3]/(diagonal[1]*diagonal[2])
        self.assertLess(abs(ratio+1), 3e-15)
        output = gate()@np.ones(4)/2
        concurrence = 2*abs(output[0]*output[3]-output[1]*output[2])
        self.assertAlmostEqual(concurrence, 1, places=14)
        midway = old.evolve(h,T/2)@v[:,3]
        probability = abs(np.vdot(b,midway))**2
        self.assertAlmostEqual(probability, 1/9, places=14)
        OBS['genuine_joint_interaction'] = dict(conditional_phase_ratio='-1',
            plus_plus_output_concurrence=short(concurrence),
            maximum_joint_mode_probability='1/9',
            joint_mode_probability_at_half_period=short(probability),
            deleting_joint_mode_changes_the_dynamics=True,
            noninteraction_restriction_of_round453_not_violated=True)

    def test_04_actual_local_readout_and_uniform_time_window(self):
        v, b, q, h0, s, h, m, effect, *_ = model()
        eig = np.linalg.eigvalsh(effect)
        self.assertGreaterEqual(eig.min(), -1e-15)
        self.assertLessEqual(eig.max(), 1+1e-15)
        np.testing.assert_allclose(h0@effect, effect@h0, atol=2e-15)
        derivative_bound = np.linalg.norm(h@effect-effect@h, 2)
        self.assertLessEqual(derivative_bound, 1+2e-15)
        states = [v@np.kron(np.eye(2)[:,a], np.ones(2)/np.sqrt(2)) for a in [0,1]]
        center = []
        for st in states:
            y = old.evolve(h,T)@st
            center.append(float(np.vdot(y,effect@y).real))
        np.testing.assert_allclose(center, [(1+np.cos(2*T))/2,(1-np.cos(2*T))/2], atol=3e-15)
        # Exact rational certificate: theta=pi*(1-1/sqrt(2))<93/100,
        # since pi<22/7 and 1/sqrt(2)>70/99.
        # -cos(pi/sqrt(2))=cos(theta)>=1-theta^2/2>11351/20000.
        from fractions import Fraction as F
        central_lower = 1-F(93,100)**2/2
        whole_lower = central_lower-2*F(1,32)
        self.assertEqual(whole_lower, F(10101,20000))
        self.assertGreater(whole_lower, F(1,2))
        self.assertLess(F(22,7)*(1-F(70,99)), F(93,100))
        self.assertLess(2*F(70,99)**2, 1)
        self.assertGreater(whole_lower-2*F(37,32)*F(1,1024), F(1,2))
        whole = []
        for t in [T-DELTA,T,T+DELTA]:
            u = old.evolve(h,t)
            pp = [float(np.vdot(u@st,effect@(u@st)).real) for st in states]
            whole.append(short(pp[1]-pp[0]))
            self.assertGreater(pp[1]-pp[0], float(whole_lower))
        OBS['physical_readout'] = dict(effect='(I+X_L_on_B)/2 on all 64 states',
            center_probabilities=[short(x) for x in center],
            center_gap=short(center[1]-center[0]),
            single_probability_derivative_bound='1', gap_derivative_bound='2',
            time_window='abs(t-pi/(2*sqrt(2)))<=1/32',
            rigorous_gap_lower='10101/20000 > 1/2',
            equal_prior_success_lower='30101/40000 > 3/4',
            raw_three_time_spot_checks=whole, postselection=False,
            physical_effect_hardware_autonomously_built=False)

    def test_05_unknown_reference_window_error_and_coefficient_error(self):
        v, b, q, h0, s, h, m, effect, *_ = model()
        self.assertAlmostEqual(np.linalg.norm(m-3*np.eye(5),2), 4, places=13)
        worst = []
        for dt in [-DELTA, DELTA]:
            actual = old.evolve(h,T+dt)@v
            ideal = np.exp(-3j*dt)*v@gate()
            error = np.linalg.norm(actual-ideal,2)
            self.assertLessEqual(error, 4*abs(dt)+2e-14)
            leakage = np.linalg.norm((np.eye(64)-v@v.T)@actual,2)**2
            self.assertLessEqual(leakage, 8*dt*dt/9+2e-14)
            worst.append(dict(operator_error=short(error),leakage=short(leakage)))
        # Perturb one raw coupling by eta; Duhamel is an analytic bound for
        # any Hermitian perturbation with this norm, not a random-noise model.
        eta = 1/256
        actual = old.evolve(h+eta*previous.swap(6,0,4),T)@v
        error = np.linalg.norm(actual-v@gate(),2)
        self.assertLessEqual(error, eta*T+2e-14)
        OBS['finite_precision'] = dict(operator_error_up_to_global_phase='<=4*abs(t-T)',
            reference_uniform_trace_distance_bound_at_window_edge='1/8',
            maximum_old_code_leakage_at_window_edge='<=1/1152',
            perturbation_Hermitian_norm_budget='eta',
            added_operator_error_bound='eta*abs(t)',
            seven_coupling_errors_individually_bounded_by_e='eta<=7e',
            endpoint_spot_checks=worst, perturbation_spot_check_error=short(error),
            uniform_readout_gap_with_eta_1_over_1024='>1/2',
            exact_gate_claim_requires_exact_coefficients_and_time=True)

    def test_06_common_resource_and_unknown_gauge_boundary(self):
        v, b, q, h0, s, h, m, effect, w, p = model()
        u = old.evolve(h,T)
        # The singlet gauge is entangled across the two actors and cannot be
        # replaced by arbitrary independent gauge preparation.
        vg = w@np.kron(np.array([[1.],[0.],[0.],[0.]]),np.eye(4))
        leakage = np.linalg.norm((np.eye(64)-p)@u@vg,2)**2
        self.assertGreater(leakage, .39)
        singlet = np.array([0.,1.,-1.,0.])/np.sqrt(2)
        gs = np.outer(singlet,singlet)
        np.testing.assert_allclose(old.partial(gs,[2,2],[0]), np.eye(2)/2, atol=1e-15)
        # Closed unitary exchange preserves global spin sector and cannot
        # turn a nonzero spin-one input into this total-spin-zero resource.
        all_pairs = sum(previous.swap(6,a,b) for a in range(6) for b in range(a+1,6))
        j2 = all_pairs-3*np.eye(64)
        np.testing.assert_allclose(j2@v, 0, atol=4e-15)
        np.testing.assert_allclose(j2@vg, 2*vg, atol=4e-15)
        np.testing.assert_array_equal(h@j2,j2@h)
        OBS['resource_boundary'] = dict(shared_gauge_singlet_entanglement_ebits=1,
            unspecified_independent_gauge_endpoint_max_leakage=short(leakage),
            full_unknown_gauge_joining_not_implemented=True,
            closed_six_spin_exchange_cannot_prepare_singlet_from_spin_one=True,
            joint_mode_carries_temporary_state_and_energy=True,
            resources=['six qubits','prepared common gauge singlet','unknown logical input',
                       'seven fixed specified exchange edges','local logical effect access',
                       'specified finite reading window'],
            network_topology_or_resource_origin_derived=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=454, baseline_round=453, date='2026-09-24',
        runtime=dict(python=platform.python_version(),numpy=np.__version__),
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        observations=OBS, scope=dict(
            enlarged_common_interface_positive_construction=True,
            fixed_equal_positive_pair_exchanges_only=True,
            exact_unknown_logic_endpoint_entangling_map=True,
            finite_window_raw_physical_readout_certified=True,
            shared_gauge_preparation_input_explicit=True,
            arbitrary_unknown_gauge_joining_implemented=False,
            autonomous_resource_preparation_or_permanent_halting=False,
            complete_SoCA_or_quantum_necessity_derived=False,
            physical_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    answer=run()
    data=json.dumps(answer,ensure_ascii=False,indent=2)+'\n'
    if not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:
            f.write(data)
    print(data)
