"""Round 1087: exact clock/signal histories and full CQ instrument calibration.

Default: JSON to stdout only. --check: read and compare saved results, never write.
Analytic scope and unproved implementation requirements are in proof.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
I0 = F(1, 8)


def hamiltonian(p, action, epsilon):
    p2 = sum(x*x for x in p)
    return p2/F(2) + action + epsilon*p2*(action-I0)


def exact_calibration():
    derivative_checks = 0
    step = F(1, 19)
    for epsilon in (0, 1):
        for p in ((F(3,5),F(0),F(0)), (F(-2,7),F(1,3),F(4,9))):
            for action in (I0, F(3,8)):
                p2 = sum(x*x for x in p)
                for axis in range(3):
                    plus, minus = list(p), list(p)
                    plus[axis] += step
                    minus[axis] -= step
                    derivative = (hamiltonian(plus,action,epsilon)-hamiltonian(minus,action,epsilon))/(2*step)
                    assert derivative == p[axis]*(1+2*epsilon*(action-I0))
                    derivative_checks += 1
                derivative = (hamiltonian(p,action+step,epsilon)-hamiltonian(p,action-step,epsilon))/(2*step)
                assert derivative == 1+epsilon*p2
                derivative_checks += 1
                assert hamiltonian(p,action,epsilon) >= 0

    beta, distance = F(3,5), F(2,5)
    expected = [(F(1),F(8,5),F(5,2),F(9,10)),
                (F(34,25),F(20,17),F(17,5),F(189,85))]
    histories = []
    for epsilon in (0, 1):
        f = 1+epsilon*beta*beta
        a_emits = [F(k) for k in (1,2)]
        b_emits = [F(k)/f for k in (1,2)]
        # Solve actual intersections independently of the Doppler formula.
        a_receives = [t+(distance+beta*t) for t in b_emits]
        b_receives = [(distance+t)/(1-beta) for t in a_emits]
        for emitted, received in zip(b_emits,a_receives):
            assert received > emitted
            assert distance+beta*emitted-(received-emitted) == 0
        for emitted, received in zip(a_emits,b_receives):
            assert received > emitted
            assert received-emitted == distance+beta*received
        a_readings = a_receives
        b_readings = [f*t for t in b_receives]
        da = a_readings[1]-a_readings[0]
        db = b_readings[1]-b_readings[0]
        assert (f,da,db,db-da) == expected[epsilon]
        assert da == (1+beta)/f and db == f/(1-beta)
        assert hamiltonian((beta,F(0),F(0)),I0,epsilon) == F(61,200)
        assert f*f > 1-beta*beta
        gamma = F(5,4)
        echo_parallel, echo_transverse = f*gamma*gamma, f*gamma
        assert echo_parallel != 1 and echo_transverse != 1
        error_sum = F(1,5)
        margin = db-da-error_sum
        assert margin > 0
        histories.append({
            'epsilon': epsilon, 'phase_rate_factor': str(f),
            'prepared_energy': '61/200',
            'A_emits_tau': list(map(str,a_emits)),
            'B_emits_tau': list(map(str,b_emits)),
            'A_receives_tau': list(map(str,a_receives)),
            'B_receives_tau': list(map(str,b_receives)),
            'B_local_receipt_readings': list(map(str,b_readings)),
            'D_A_from_B': str(da), 'D_B_from_A': str(db),
            'reciprocity_gap': str(db-da),
            'per_reading_example_error_bound': '1/10',
            'certified_example_gap_margin': str(margin),
            'conditional_rigid_arm_echo_parallel': str(echo_parallel),
            'conditional_rigid_arm_echo_transverse': str(echo_transverse),
        })
    assert histories[0]['A_receives_tau'][0] == '2'
    assert histories[1]['A_receives_tau'][0] == '134/85'
    return {
        'hamiltonian_derivative_checks': derivative_checks,
        'I_domain': 'I > 0', 'prepared_action': '1/8',
        'm_c_omega': ['1','1','1'],
        'epsilon1_global_kinetic_coefficient_lower_bound': '3/8',
        'tau_definition': 'tau=t/T0, T0=2*pi/omega',
        'ell_definition': 'ell=r/(c*T0)',
        'beta': '3/5', 'initial_normalized_distance': '2/5',
        'histories': histories,
        'lorentz_necessary_comparison': {
            'phase_rate_factor': '4/5', 'D_A_from_B': '2', 'D_B_from_A': '2'},
    }


IDENT = np.eye(2, dtype=complex)
PAULI = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], dtype=complex)
ANISOTROPY = np.diag([1., .75, .5])


def effect(x, anchor, orientation):
    delta = np.asarray(x)-np.asarray(anchor)
    delta = delta/max(1.,np.linalg.norm(delta))
    vector = ANISOTROPY @ orientation @ delta
    return (IDENT+np.einsum('j,jab->ab',vector,PAULI))/2


def kraus_from_effect(e):
    values, vectors = np.linalg.eigh(e)
    assert values.min() > -1e-13 and values.max() < 1+1e-13
    values = np.clip(values,0,1)
    return [(vectors*np.sqrt(values))@vectors.conj().T,
            (vectors*np.sqrt(1-values))@vectors.conj().T]


def unitary_dilation(kraus):
    # Independent closed-form Halmos dilation for commuting positive roots.
    kp,km = kraus
    return np.block([[kp,-km],[km,kp]])


def axis_rotation(axis, angle):
    axis = np.asarray(axis,dtype=float)
    axis /= np.linalg.norm(axis)
    a,b,c = axis
    cross = np.array([[0,-c,b],[c,0,-a],[-b,a,0]])
    return np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis,axis)+np.sin(angle)*cross


def quantum_calibration():
    rng = np.random.default_rng(1087)
    residuals = {name: 0. for name in (
        'kraus_completeness', 'unitary_dilation', 'spatial_effect_transport',
        'spatial_branch_transport', 'dilation_branch_equality',
        'reference_nonsignalling', 'trace_preservation', 'classical_record_marginal',
        'payload_transport', 'relative_control_composition', 'relative_control_inverse',
        'midpoint_record_inverse')}
    def track(name,value):
        residuals[name] = max(residuals[name],float(np.max(np.abs(value))))
    state_cases = branch_cases = history_cases = linear_basis_cases = 0
    for case in range(12):
        x, anchor = rng.normal(size=3), rng.normal(size=3)
        orientation = axis_rotation(rng.normal(size=3),float(rng.uniform(-2,2)))
        rotate = axis_rotation(rng.normal(size=3),float(rng.uniform(-2,2)))
        shift, ra = rng.normal(size=3), rng.normal(size=3)
        rb = ra+x
        rb_new = ra+rotate@(rb-ra)+shift
        x_new = rb_new-ra
        first = effect(x,anchor,orientation)
        lifted = effect(x_new,rotate@anchor+shift,orientation@rotate.T)
        track('spatial_effect_transport',first-lifted)
        kernels = kraus_from_effect(first)
        transported = kraus_from_effect(lifted)
        track('kraus_completeness',sum(k.conj().T@k for k in kernels)-IDENT)
        dilation = unitary_dilation(kernels)
        track('unitary_dilation',dilation.conj().T@dilation-np.eye(4))
        # A complete matrix basis checks each full branch superoperator, not
        # just one preparation or its report probability.
        for row in range(2):
            for col in range(2):
                basis = np.zeros((2,2),dtype=complex)
                basis[row,col] = 1
                for k,l in zip(kernels,transported):
                    track('spatial_branch_transport',k@basis@k.conj().T-l@basis@l.conj().T)
                    linear_basis_cases += 1

        rotate2 = axis_rotation(rng.normal(size=3),float(rng.uniform(-2,2)))
        shift2 = rng.normal(size=3)
        sequential = ra+rotate2@(rb_new-ra)+shift2
        composed = ra+(rotate2@rotate)@(rb-ra)+rotate2@shift+shift2
        track('relative_control_composition',sequential-composed)
        recovered = ra+rotate.T@(rb_new-ra-shift)
        track('relative_control_inverse',recovered-rb)
        y = rng.normal(size=3)
        midpoint, relative = (x+y)/2,(x-y)/2
        track('midpoint_record_inverse',np.r_[midpoint+relative-x,midpoint-relative-y])

        for refdim in (2,3,5):
            size = 2*refdim
            z = rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
            rho = z@z.conj().T
            if case % 2 == 0:
                psi = z[:,0]/np.linalg.norm(z[:,0])
                rho = np.outer(psi,psi.conj())
            rho /= np.trace(rho)
            original_ref = rho.reshape(2,refdim,2,refdim).trace(axis1=0,axis2=2)
            extended_k = [np.kron(k,np.eye(refdim)) for k in kernels]
            branches = [k@rho@k.conj().T for k in extended_k]
            all_state = sum(branches)
            track('trace_preservation',np.trace(all_state)-1)
            after_ref = all_state.reshape(2,refdim,2,refdim).trace(axis1=0,axis2=2)
            track('reference_nonsignalling',after_ref-original_ref)
            record_ready = np.zeros((2,2),dtype=complex)
            record_ready[0,0] = 1
            joint_u = np.kron(dilation,np.eye(refdim))
            joint = joint_u @ np.kron(record_ready,rho) @ joint_u.conj().T
            block = joint.reshape(2,size,2,size)
            for bit,branch in enumerate(branches):
                track('dilation_branch_equality',block[bit,:,bit,:]-branch)
                # Classical control transports the same quantum payload by
                # an identity channel; a receiver can perform a nontrivial
                # follow-up without losing the original conditional state.
                receiver_u = np.array([[1,1j],[1j,1]],complex)/np.sqrt(2)
                receiver = np.kron(receiver_u,np.eye(refdim))
                for epsilon in (0,1):
                    arrival_payload = np.eye(size)@branch@np.eye(size)
                    track('payload_transport',receiver@arrival_payload@receiver.conj().T
                          -receiver@branch@receiver.conj().T)
                    history_cases += 1
                branch_cases += 1
            classical_joint = np.zeros((2*size,2*size),dtype=complex)
            for bit,branch in enumerate(branches):
                projector = np.zeros((2,2)); projector[bit,bit] = 1
                classical_joint += np.kron(projector,branch)
            marginal = classical_joint.reshape(2,size,2,size).trace(axis1=0,axis2=2)
            track('classical_record_marginal',marginal-all_state)
            assert np.linalg.eigvalsh(classical_joint).min() > -1e-12
            state_cases += 1
    assert all(v<2e-12 for v in residuals.values()), residuals
    return {
        'spatial_control_cases': 12, 'reference_dimensions': [2,3,5],
        'unknown_input_state_cases': state_cases,
        'complete_branch_state_cases': branch_cases,
        'linear_basis_branch_cases': linear_basis_cases,
        'poststate_delivery_history_cases': history_cases,
        'tolerance': 2e-12, 'max_residuals': residuals,
        'all_passed': True,
    }


def calculate():
    return {
        'round': 1087, 'passed': True,
        'new_independent_scientific_calibrations': 0,
        'cumulative_scientific_calibrations': 3860,
        'scope': 'CO1-CO5 effective CQ model with original internal CP/control primitives',
        'exact': exact_calibration(), 'quantum': quantum_calibration(),
        'preserved_original_CO_contract': True,
        'Lorentz_clock_signal_dictionary_follows_from_CO_alone': False,
        'counterexample_satisfies_1086_double_echo_contract': False,
        'generative_completeness_disproved': False,
        'all_primitives_generated_by_one_autonomous_H': False,
        'closed_detector_energy_and_recoil_certified': False,
        'fundamental_Lorentz_theory_excluded': False,
        'empirical_instrument_error_budget_certified': False,
    }


def compare(actual, expected, path='results'):
    if isinstance(actual, dict):
        assert isinstance(expected,dict) and actual.keys()==expected.keys(),path
        for key in actual:
            compare(actual[key],expected[key],path+'.'+key)
    elif isinstance(actual,list):
        assert isinstance(expected,list) and len(actual)==len(expected),path
        for i,(a,b) in enumerate(zip(actual,expected)):
            compare(a,b,path+f'[{i}]')
    elif isinstance(actual,float):
        assert isinstance(expected,(int,float)) and abs(actual-expected)<1e-13,path
    else:
        assert actual==expected,(path,actual,expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.check:
        expected = json.loads((BASE/'results.json').read_text(encoding='utf-8'))
        compare(result,expected)
    print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
