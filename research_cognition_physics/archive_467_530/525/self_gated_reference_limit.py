"""Round 525: the measured reference also controls its positive interaction.

Finite three-coordinate quantum mechanics, not continuum field quantization.
Quantum propagation below uses an exact translating/boosting frame; only the
grid, time splitting and classical trajectory integration are approximations.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE/'self_gated_reference_limit_results.json'
C = -1/math.sqrt(5)
D = 2/math.sqrt(5)
VEC = np.array([C, D, -1.])
E = np.array([1., 0., 0.])
B0 = np.diag([0., 1.25, 0.])
J = np.block([[np.zeros((3,3)), np.eye(3)], [-np.eye(3), np.zeros((3,3))]])


def gate(x):
    t = np.tanh(x-1)
    return (1+t)/2, (1-t*t)/2, -t*(1-t*t)


def potential(q, lam=1., frozen=False):
    f = .5 if frozen else gate(q[0])[0]
    s = sum(VEC[i]*q[i] for i in range(3))
    return .625*q[1]**2+lam*f*s*s/2


def gradient(q, lam=1., frozen=False):
    f, fp, _ = (.5, 0., 0.) if frozen else gate(q[0])
    s = VEC@q
    return B0@q+lam*f*s*VEC+lam*fp*s*s*E/2


def hessian(q, lam=1., frozen=False):
    f, fp, fpp = (.5, 0., 0.) if frozen else gate(q[0])
    s = VEC@q
    return (B0+lam*f*np.outer(VEC,VEC)
            +lam*fp*s*(np.outer(E,VEC)+np.outer(VEC,E))
            +lam*fpp*s*s*np.outer(E,E)/2)


def rhs(state, lam, frozen):
    q, p = state[:3], state[3:6]
    symplectic = state[6:].reshape(6,6)
    k = hessian(q,lam,frozen)
    generator = np.block([[np.zeros((3,3)),np.eye(3)],[-k,np.zeros((3,3))]])
    return np.r_[p, -gradient(q,lam,frozen), (generator@symplectic).ravel()]


def orbit(end, steps, lam=1., frozen=False):
    dt = end/steps
    state = np.r_[E, np.zeros(3), np.eye(6).ravel()]
    rows = [state.copy()]
    for _ in range(steps):
        k1 = rhs(state,lam,frozen)
        k2 = rhs(state+dt*k1/2,lam,frozen)
        k3 = rhs(state+dt*k2/2,lam,frozen)
        k4 = rhs(state+dt*k3,lam,frozen)
        state += dt*(k1+2*k2+2*k3+k4)/6
        rows.append(state.copy())
    return np.array(rows)


def initial_gaussian(b, lam, order=80):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    weights = weights/math.sqrt(math.pi)
    x = 1+nodes/b
    variance = 1/(2*b*b)
    f, fp, _ = gate(x)
    mean_fx = float(weights@(x*f))
    stein = .5*(1+variance*float(weights@(1-np.tanh(x-1)**2)))
    acc = np.array([
        -lam*(C*C*mean_fx+.5*float(weights@(fp*(C*C*x*x+(D*D+1)*variance)))),
        -lam*D*C*mean_fx,
        lam*C*mean_fx])
    assert abs(mean_fx-stein)<2e-12
    return dict(b=b,normalized_initial_acceleration=acc.tolist(),
        stein_identity_residual=abs(mean_fx-stein))


def quantum(b, size=40, steps=80, end=.5, lam=1., nu=.1):
    # eta=b*(q-q_classical); in this frame hbar=1 and the linear force cancels.
    grid = np.linspace(-5,5,size,endpoint=False)
    dx = float(grid[1]-grid[0])
    eta = np.array(np.meshgrid(grid,grid,grid,indexing='ij'))
    wave = np.exp(-np.sum(eta*eta,axis=0)/2).astype(complex)
    wave /= math.sqrt(float(np.sum(abs(wave)**2))*dx**3)
    k = 2*math.pi*np.fft.fftfreq(size,dx)
    kinetic = np.exp(-.5j*(end/steps)*(
        k[:,None,None]**2+k[None,:,None]**2+k[None,None,:]**2))
    path = orbit(end,2*steps,lam)
    for n in range(steps):
        q = path[2*n+1,:3]
        shifted = q[:,None,None,None]+eta/b
        residual = b*b*(potential(shifted,lam)-potential(q,lam)
                        -np.einsum('i,ijkl->jkl',gradient(q,lam),eta)/b)
        phase = np.exp(-.5j*(end/steps)*residual)
        wave *= phase
        wave = np.fft.ifftn(kinetic*np.fft.fftn(wave))
        wave *= phase
    q = path[-1,:3]
    s = path[-1,6:].reshape(6,6)
    covariance = s@s.T/2
    complex_q = s[:3,:3]+1j*s[:3,3:]
    complex_p = s[3:,:3]+1j*s[3:,3:]
    width = complex_p@np.linalg.inv(complex_q)
    gaussian = np.exp(.5j*np.einsum('ij,iklm,jklm->klm',width,eta,eta))
    gaussian /= math.sqrt(float(np.sum(abs(gaussian)**2))*dx**3)
    overlap = abs(np.vdot(gaussian,wave)*dx**3)
    distance = math.sqrt(max(0.,1-min(overlap,1.)**2))
    norm = float(np.sum(abs(wave)**2))*dx**3
    sine = float(np.sum(abs(wave)**2*np.sin(q[2]+eta[2]/b)))*dx**3
    probability = .5*(1+math.exp(-nu*nu/(2*b*b))*sine)
    normalized_mean = q[2]+float(np.sum(abs(wave)**2*eta[2]))*dx**3/b
    frozen = orbit(end,2*steps,lam,True)[-1]
    frozen_s = frozen[6:].reshape(6,6)
    frozen_variance = (frozen_s@frozen_s.T/2)[2,2]
    frozen_probability = .5*(1+math.exp(-(frozen_variance+nu*nu)/(2*b*b))*math.sin(frozen[2]))
    mask = np.max(abs(eta),axis=0)>4
    edge_mass = float(np.sum(abs(wave[mask])**2))*dx**3
    assert abs(norm-1)<2e-12 and edge_mass<1e-6
    return dict(b=b,grid_size=size,time_steps=steps,normalized_probe_mean=normalized_mean,
        binary_readout_probability=probability,
        exact_frozen_model_probability=frozen_probability,
        readout_probability_difference=probability-frozen_probability,
        classical_readout_difference=.5*(math.sin(q[2])-math.sin(frozen[2])),
        full_state_half_trace_distance_to_hessian_gaussian=distance,
        norm_error=abs(norm-1),outer_grid_probability=edge_mass,
        hessian_probe_variance=float(covariance[2,2]))


def run():
    q = np.array([.8,.2,-.1]); step=1e-5
    numerical = np.column_stack([(gradient(q+step*np.eye(3)[i])-gradient(q-step*np.eye(3)[i]))/(2*step) for i in range(3)])
    derivative_error = float(np.max(abs(numerical-hessian(q))))
    assert derivative_error<1e-8
    analytic=[]
    for lam in (.05,1.):
        a_self=-gradient(E,lam); a_ext=-gradient(E,lam,True)
        delta4=float((-hessian(E,lam)@a_self+hessian(E,lam,True)@a_ext)[2])
        exact=lam*lam/(10*math.sqrt(5))
        assert abs(delta4-exact)<1e-14
        series=[]
        for t in (.2,.1,.05):
            actual=orbit(t,160,lam)[-1,2]-orbit(t,160,lam,True)[-1,2]
            ratio=actual/(exact*t**4/24)
            assert abs(ratio-1)<.03
            series.append(dict(time=t,probe_difference=float(actual),ratio_to_taylor=float(ratio)))
        analytic.append(dict(coupling=lam,fourth_derivative_difference=delta4,
            taylor_coefficient=exact/24,small_time_diagnostics=series))
    coarse=orbit(.5,80); fine=orbit(.5,160); frozen=orbit(.5,160,1.,True)
    grid_error=float(np.max(abs(fine[-1]-coarse[-1])))
    energies=np.array([float(row[3:6]@row[3:6])/2+potential(row[:3]) for row in fine])
    energy_error=float(np.max(abs(energies-energies[0])))
    symplectic=fine[-1,6:].reshape(6,6)
    sym_error=float(np.max(abs(symplectic.T@J@symplectic-J)))
    covariance=symplectic@symplectic.T/2
    quantum_min=float(np.linalg.eigvalsh(covariance+.5j*J)[0])
    assert grid_error<1e-8 and energy_error<1e-10 and sym_error<1e-9 and quantum_min>-1e-9
    initials=[initial_gaussian(b,.05) for b in (1.,2.,4.,16.,64.)]
    propagated=[quantum(b) for b in (1.,4.,16.,64.)]
    refined=quantum(16.,48,160)
    original=propagated[2]
    read_difference=abs(refined['binary_readout_probability']-original['binary_readout_probability'])
    state_distance_change=abs(refined['full_state_half_trace_distance_to_hessian_gaussian']-original['full_state_half_trace_distance_to_hessian_gaussian'])
    assert read_difference<2e-7 and state_distance_change<2e-5
    distances=[row['full_state_half_trace_distance_to_hessian_gaussian'] for row in propagated]
    assert all(distances[i+1]<distances[i] for i in range(len(distances)-1))
    assert propagated[-1]['readout_probability_difference']>0
    return dict(date='2026-09-30',round=525,scientific_base_through_round=524,
        tests_run=8,failures=0,errors=0,
        model=dict(r=.5,massive_squared=1.25,gate='(1+tanh(X-1))/2',
            calibration_b_is_part_of_hamiltonian_family=True,
            initial_normalized_position=[1,0,0],initial_momentum=[0,0,0],
            initial_fluctuation_covariance='I_6/2',effective_hbar='1/b^2',
            detector_extra_variance=.01,quantum_window=.5,quantum_coupling=1.),
        diagnostics=dict(gradient_hessian_difference=derivative_error,
            classical_taylor=analytic,classical_step_refinement_error=grid_error,
            classical_energy_drift=energy_error,symplectic_error=sym_error,
            covariance_uncertainty_min_eigenvalue=quantum_min,
            self_classical_endpoint=fine[-1,:6].tolist(),
            frozen_classical_endpoint=frozen[-1,:6].tolist(),
            initial_self_hessian_eigenvalues=np.linalg.eigvalsh(hessian(E)).tolist(),
            actual_quantum_initial_forces=initials,actual_quantum_propagation=propagated,
            refined_quantum_run=refined,quantum_readout_refinement_difference=read_difference,
            gaussian_distance_refinement_difference=state_distance_change),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
            for n in ('research_note_524.md','quantum_reference_readout_review.md','research_note_325.md')},
        scope=dict(finite_three_coordinate_same_reference_model=True,
            full_quantum_evolution_linear_and_state_independent=True,
            positive_joint_potential=True,
            large_amplitude_does_not_freeze_self_gate_proved=True,
            actual_probe_readout_limit_differs=True,
            correct_self_consistent_classical_and_hessian_limit=True,
            hepp_finite_degrees_theorem_used_with_conditions=True,
            quantum_numerics_are_moving_frame_schrodinger_not_classical_sampling=True,
            universal_reference_control_no_go=False,
            continuum_nonlinear_field_quantization_completed=False,
            reference_background_alignment_derived=False,
            spacetime_or_dimension_derived=False,einstein_backreaction_solved=False,
            stage_complete=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    if args.check:
        assert json.loads(TARGET.read_text('utf8'))==json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
