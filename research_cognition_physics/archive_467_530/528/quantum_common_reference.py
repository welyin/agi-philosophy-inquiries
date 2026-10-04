"""Round 528: quantum feedback on a classically flat common reference.

Full finite-coordinate Schrodinger propagation, not a continuum quantum field
calculation. A second, explicitly different relative-only gate tests a repair.
The analytic identities are established in research_note_528.md; FFT grids are
diagnostics, never interval certificates or proofs of continuum convergence.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from self_gated_reference_limit import C, D, gate

HERE = Path(__file__).resolve().parent
TARGET = HERE/'quantum_common_reference_results.json'
A = math.sqrt(1+C*C)
MU2 = 1.25
KAPPA = (9-math.sqrt(33))/8


def relative_gate(s):
    return .5+.25*np.tanh(s)**2


def initial_force(b, x=1., order=100):
    eta, weight = np.polynomial.hermite.hermgauss(order)
    fp = gate(x+eta/b)[1]
    force = -float(weight@(fp*(.9+.2*eta*eta)))/math.sqrt(math.pi)/(2*A*b)
    return dict(b=b, x=x, initial_common_force=force,
                b_times_force=b*force, asymptotic_coefficient=-float(gate(x)[1])/(2*A))


def propagate(b, relative=False, size=40, steps=128, halfwidth=6., end=.5, nu=.1):
    """Centered coordinates (u=U-A*b, M, W); center is classically stationary."""
    grid = np.linspace(-halfwidth, halfwidth, size, endpoint=False)
    dx = float(grid[1]-grid[0]); volume = dx**3
    u, m, w = np.meshgrid(grid, grid, grid, indexing='ij')
    s = A*w+D*m
    f, fp, _ = gate(1+(u+C*w)/(A*b))
    potential = MU2*m*m/2+(relative_gate(s/b) if relative else f)*s*s/2
    force = np.zeros_like(u) if relative else -fp*s*s/(2*A*b)
    wave = np.exp(-(u*u+m*m+w*w)/2).astype(complex)
    wave /= math.sqrt(float(np.sum(abs(wave)**2))*volume)
    k = 2*math.pi*np.fft.fftfreq(size, dx)
    k2 = k[:,None,None]**2+k[None,:,None]**2+k[None,None,:]**2
    dt = end/steps
    kinetic = np.exp(-.5j*dt*k2)
    phase = np.exp(-.5j*dt*potential)

    def diagnostics(psi):
        prob = abs(psi)**2
        ft = np.fft.fftn(psi, norm='ortho')
        p_u = float(np.sum(abs(ft)**2*k[:,None,None]))*volume
        energy = float(np.sum(abs(ft)**2*k2/2)+np.sum(prob*potential))*volume
        acc = float(np.sum(prob*force))*volume
        return p_u, energy, acc

    p0, energy0, acc0 = diagnostics(wave)
    force_integral = 0.; old_acc = acc0
    max_energy_error = 0.
    for _ in range(steps):
        wave *= phase
        wave = np.fft.ifftn(kinetic*np.fft.fftn(wave))
        wave *= phase
        p_u, energy, acc = diagnostics(wave)
        force_integral += dt*(old_acc+acc)/2
        old_acc = acc
        max_energy_error = max(max_energy_error, abs(energy-energy0))
    prob = abs(wave)**2
    def mean(value):
        return float(np.sum(prob*value))*volume
    def variance(value):
        return mean(value*value)-mean(value)**2
    u_mean, w_mean = mean(u), mean(w)
    estimate = u-w/C
    estimate_mean = mean(estimate)
    estimate_variance = variance(estimate)+A*A*nu*nu/(C*C)
    covariance_uw = mean(u*w)-u_mean*w_mean
    chi_centered = (C*u-w)/A
    # A genuine bounded effect on the noisy chi readout, evaluated by convolution.
    read_probability = .5*(1+math.exp(-nu*nu/2)*mean(np.sin(chi_centered)))
    w_var = variance(w)
    w_fourth_cumulant = mean((w-w_mean)**4)-3*w_var*w_var
    flatten = wave.reshape(size,-1)*dx**1.5
    rho_u = flatten@flatten.conj().T
    purity = float(np.trace(rho_u@rho_u).real)
    product_error = None
    relative_energy = None
    if relative:
        free = np.exp(-grid*grid/2).astype(complex)
        free /= math.sqrt(float(np.sum(abs(free)**2))*dx)
        free = np.fft.ifft(np.exp(-.5j*end*k*k)*np.fft.fft(free))
        rel = np.exp(-(m[0]**2+w[0]**2)/2).astype(complex)
        rel /= math.sqrt(float(np.sum(abs(rel)**2))*dx*dx)
        rel_phase = np.exp(-.5j*dt*potential[0])
        rel_kinetic = np.exp(-.5j*dt*(k[:,None]**2+k[None,:]**2))
        for _ in range(steps):
            rel *= rel_phase
            rel = np.fft.ifftn(rel_kinetic*np.fft.fftn(rel))
            rel *= rel_phase
        product = free[:,None,None]*rel[None,:,:]
        product_error = math.sqrt(float(np.sum(abs(wave-product)**2))*volume)
        relative_energy = energy0-.25
    norm_error = abs(mean(np.ones_like(u))-1)
    outer_mass = mean((np.maximum.reduce([abs(u),abs(m),abs(w)])>halfwidth-1).astype(float))
    return dict(b=b, relative_gate=relative, grid_size=size, time_steps=steps,
        halfwidth=halfwidth, end=end, norm_error=norm_error, outer_grid_probability=outer_mass,
        energy_initial=energy0, max_energy_drift=max_energy_error,
        initial_force=acc0, common_momentum=p_u, integrated_force=force_integral,
        ehrenfest_residual=abs(p_u-p0-force_integral), common_displacement=u_mean,
        relative_mean=w_mean, common_variance=variance(u), relative_variance=w_var,
        common_relative_covariance=covariance_uw, relative_fourth_cumulant=w_fourth_cumulant,
        calibrated_readout_bias=estimate_mean-u_mean,
        calibrated_readout_variance=estimate_variance,
        separated_variance_formula=variance(u)+w_var/(C*C)+A*A*nu*nu/(C*C),
        centered_binary_readout_probability=read_probability,
        common_reduced_purity=purity, full_vs_factorized_wave_norm=product_error,
        relative_energy_initial=relative_energy,
        relative_second_moment_bound=None if relative_energy is None else 2*relative_energy/KAPPA)


def graph_identity():
    rng = np.random.default_rng(528)
    n=4
    lap = np.array([[1,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,1.]])
    q = rng.normal(size=(n,3)); p = rng.normal(size=(n,3))
    # Original order R,M,chi; new order U,M,W.
    rotation = np.array([[1/A,0,C/A],[0,1,0],[C/A,0,-1/A]])
    z = q@rotation.T; pi = p@rotation.T
    old_quad = float(np.sum(p*p)/2+np.sum(q*(lap@q))/2)
    new_quad = float(np.sum(pi*pi)/2+np.sum(z*(lap@z))/2)
    s = C*q[:,0]+D*q[:,1]-q[:,2]
    b=2.
    old_potential = MU2*q[:,1]**2/2+gate(q[:,0]/b)[0]*s*s/2
    rotated_potential = MU2*z[:,1]**2/2+gate((z[:,0]+C*z[:,2])/(A*b))[0]*s*s/2
    local_a = gate(q[:,0]/b)[1]*s*s/(2*A*b)
    total_p_force = -lap@z[:,0]-local_a
    energies = old_quad+float(np.sum(old_potential))
    # Repair: derivatives along U vanish pointwise; no reference to a Gaussian limit.
    shifted = q+1.3*rotation[0]
    shifted_s = C*shifted[:,0]+D*shifted[:,1]-shifted[:,2]
    repair_v = MU2*q[:,1]**2/2+relative_gate(s/b)*s*s/2
    repair_shift_v = MU2*shifted[:,1]**2/2+relative_gate(shifted_s/b)*shifted_s**2/2
    stiff = np.array([[0.,0.],[0.,MU2]])+.5*np.outer([A,D],[A,D])
    return dict(orthogonal_error=float(np.max(abs(rotation@rotation.T-np.eye(3)))),
        graph_kinetic_split_error=abs(old_quad-new_quad),
        rotated_potential_error=float(np.max(abs(old_potential-rotated_potential))),
        total_force_sum_error=abs(float(np.sum(total_p_force)+np.sum(local_a))),
        positive_local_force_coefficients=local_a.tolist(),
        force_energy_bound_slack=2*energies/(A*b)-float(np.sum(local_a)),
        repair_common_shift_error=float(np.max(abs(repair_v-repair_shift_v))),
        coercivity_eigenvalues=np.linalg.eigvalsh(stiff).tolist(), analytic_kappa=KAPPA)


def run():
    identities = graph_identity()
    for name in ('orthogonal_error','graph_kinetic_split_error','rotated_potential_error',
                 'total_force_sum_error','repair_common_shift_error'):
        assert identities[name]<2e-14, (name,identities[name])
    assert min(identities['positive_local_force_coefficients'])>0
    assert identities['force_energy_bound_slack']>0
    assert abs(identities['coercivity_eigenvalues'][0]-KAPPA)<1e-14
    initial = [initial_force(b) for b in (1.,4.,16.,64.,256.)]
    refined_initial = initial_force(1.,order=160)
    assert abs(initial[0]['initial_common_force']-refined_initial['initial_common_force'])<1e-11
    assert abs(initial[-1]['b_times_force']+1/(4*A))<3e-6
    other_x = initial_force(4.,x=.75)
    assert abs(other_x['initial_common_force']-initial[1]['initial_common_force'])>1e-3
    old = [propagate(b) for b in (1.,4.,16.)]
    for row in old:
        assert row['common_momentum']<0 and row['common_displacement']<0
        assert row['ehrenfest_residual']<2e-6
        assert abs(row['initial_force']-initial_force(row['b'])['initial_common_force'])<2e-10
    repaired = [propagate(b,True) for b in (1.,4.)]
    for row in repaired:
        assert row['full_vs_factorized_wave_norm']<1e-11
        assert abs(row['common_momentum'])<1e-9 and abs(row['relative_mean'])<1e-9
        assert abs(row['common_reduced_purity']-1)<1e-11
        assert abs(row['calibrated_readout_bias'])<1e-9
        assert abs(row['calibrated_readout_variance']-row['separated_variance_formula'])<1e-9
        assert row['relative_variance']<row['relative_second_moment_bound']
    assert abs(repaired[0]['relative_fourth_cumulant'])>1e-5
    fine_old = propagate(1.,size=56,steps=256,halfwidth=7.)
    fine_repair = propagate(1.,True,size=56,steps=256,halfwidth=7.)
    errors={}
    for name,coarse,fine in (('old',old[0],fine_old),('repair',repaired[0],fine_repair)):
        errors[name] = {key:abs(coarse[key]-fine[key]) for key in (
            'common_momentum','common_displacement','calibrated_readout_variance',
            'centered_binary_readout_probability','relative_fourth_cumulant')}
        assert max(errors[name].values())<2e-5,errors[name]
    for row in old+repaired+[fine_old,fine_repair]:
        assert row['norm_error']<2e-12
        assert row['outer_grid_probability']<1e-8
        assert row['max_energy_drift']<2e-5
    dependencies = ('self_gated_reference_limit.py','research_note_525.md',
        'research_note_527.md','reference_alignment_completion_review.md',
        'material_reference_geometry_review.md')
    return dict(date='2026-09-30',round=528,scientific_base_through_round=527,
        tests_run=8,failures=0,errors=0,
        model=dict(c=C,d=D,a=A,mu_squared=MU2,repair_gate='.5+.25*tanh(S/b)^2',
                   finite_closed_graph=True,shared_reference_source_not_blank_probe=True),
        graph_and_force_identities=identities,initial_force_quadrature=initial,
        refined_initial_force=refined_initial,nonuniform_source_force=other_x,
        original_rule_quantum=old,relative_rule_quantum=repaired,
        refined_original=fine_old,refined_repair=fine_repair,refinement_errors=errors,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in dependencies},
        scope=dict(exact_finite_energy_common_force_proved=True,
            old_classical_flat_valley_not_exact_free_quantum_factor=True,
            macroscopic_free_approximation_not_excluded=True,
            repair_is_explicit_rule_change=True,repair_exact_common_factor=True,
            actual_relative_noise_retained=True,full_schrodinger_not_hessian_only=True,
            spatially_separating_source_generated=False,continuum_field_limit_proved=False,
            quantum_einstein_feedback_closed=False,three_dimensions_derived=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert json.loads(TARGET.read_text('utf8'))==json.loads(json.dumps(result))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(round=528,tests=result['tests_run'],failures=0,
        old_common_momentum=result['original_rule_quantum'][0]['common_momentum'],
        repair_factorization_error=result['relative_rule_quantum'][0]['full_vs_factorized_wave_norm']),ensure_ascii=False))
