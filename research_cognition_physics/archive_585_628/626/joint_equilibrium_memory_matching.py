"""626: conserved thermal sources constrain source-only local effective actions.

The full fixed-graph Gauss result is analytic. Numerical spectra are the
unchanged 598/602 frozen 32-mode mass sector, with an independent neutral
16-dimensional Fock calculation. They are not full Gauss spectra.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quantum_response_matching as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_equilibrium_memory_matching_results.json'
BETA=2.
HBAR=1.


def data():
    bdg=old.matrices(old.S0)[0][0]
    eigen=np.linalg.eigvalsh(bdg)
    assert np.max(abs(eigen+eigen[::-1]))<1e-12
    eps=eigen[32:]
    assert eps.min()>0
    f=1/(1+np.exp(BETA*eps))
    e0=-.5*eps.sum()  # tr(h)=0 in the original mass matrices.
    mean=float(e0+eps@f)
    var=float(np.sum(eps**2*f*(1-f)))
    k3=float(np.sum(eps**3*f*(1-f)*(1-2*f)))
    return dict(eps=eps,f=f,e0=float(e0),mean=mean,variance=var,kappa3=k3)


def characteristic(d,tau):
    return np.exp(-1j*d['e0']*tau/HBAR)*np.prod(
        1-d['f']+d['f']*np.exp(-1j*d['eps']*tau/HBAR))


def encoded(z):
    return [float(z.real),float(z.imag)]


def bump(order=128):
    x,w=np.polynomial.legendre.leggauss(order)
    t=(x+1)/2;w=w/2
    raw=np.exp(-1/(t*(1-t)))
    normal=float(w@raw)
    value=raw/normal
    # C-infinity bump, extended by zero outside (0,1).
    derivative=value*(1-2*t)/(t*t*(1-t)**2)
    return dict(t=t,w=w,value=value,derivative=derivative,normal=normal,
                area=float(w@value),square=float(w@(value*value)))


def spectral_check(d):
    previous=old.response(beta=BETA)
    err=max(abs(d['mean']-previous['mean'][0]),
            abs(d['variance']-previous['D'][0,0]))
    assert err<2e-13
    _,pairs=old.matrices(old.S0,[24,25,30,31])
    H=old.fock(pairs[0])
    e,v,p,rho=old.thermal(H,BETA)
    q=np.linalg.eigvalsh(old.matrices(old.S0,[24,25,30,31])[0][0])[4:]
    neutral=dict(eps=q,f=1/(1+np.exp(BETA*q)),e0=float(-.5*q.sum()))
    rows=[]
    for tau in (0.,.08,.4,1.1):
        U=(v*np.exp(-1j*e*tau))@v.conj().T
        direct=np.trace(U@rho)
        factored=characteristic(neutral,tau)
        error=float(abs(direct-factored))
        assert error<3e-13
        rows.append(dict(tau=tau,Fock_trace=encoded(direct),
                         quasiparticle_product=encoded(factored),error=error))
    # Symmetric differences are diagnostics of the exact characteristic function.
    eps=.001
    lplus=np.log(characteristic(d,2*eps));lminus=np.log(characteristic(d,-2*eps))
    coefficient=(lplus+lminus)/(2*eps*eps)
    expected=-2*d['variance']/HBAR**2
    assert abs(coefficient-expected)<2e-6
    return dict(CAR_modes=32,BdG_dimension=64,beta=BETA,
                mean_energy=d['mean'],energy_variance=d['variance'],
                third_energy_cumulant=d['kappa3'],previous_source_identity_error=float(err),
                independent_neutral_Fock=rows,
                two_unit_area_pulses_log_quadratic=encoded(coefficient),
                exact_log_quadratic=expected,
                finite_difference_error=float(abs(coefficient-expected)))


def locality_check(d):
    a=bump(128);b=bump(256)
    quad_error=abs(a['normal']-b['normal'])
    assert quad_error<1e-15
    assert abs(a['w']@a['derivative'])<2e-13
    # f supported in (0,1), g=f(t-2) in (2,3). Supports stay disjoint
    # under all finite derivatives. The theorem covers arbitrary local
    # diagonal-supported Hessians, not just this illustrative white kernel.
    V=d['variance']
    q_plus=4*V*a['area']**2
    q_minus=0.
    white_strength=V/a['square']
    local_both=2*white_strength*a['square']
    error=max(abs(local_both-q_plus),abs(local_both-q_minus))
    lower_bound=2*V*a['area']**2
    assert abs(error-lower_bound)<2e-13
    # Keeping this fixed local coefficient, broaden each pulse and reduce
    # its amplitude to keep its area. Slowness does not remove the zero mode.
    slow=[]
    for scale in (1.,4.,16.):
        local=2*white_strength*a['square']/scale
        slow.append(dict(width_scale=scale,exact_plus_variance=q_plus,
                         exact_compensated_variance=q_minus,
                         fixed_white_plus_and_minus_variance=local,
                         plus_error=abs(q_plus-local)))
    return dict(unit_area=a['area'],bump_square_integral=a['square'],
                quadrature_normalization_difference=quad_error,
                endpoint_preserving_derivative_area=float(a['w']@a['derivative']),
                true_plus_variance=q_plus,true_compensated_variance=q_minus,
                one_pulse_calibrated_local_variance_for_both=local_both,
                any_local_quadratic_kernel_minimax_variance_error=lower_bound,
                any_local_log_quadratic_minimax_error=lower_bound/(2*HBAR**2),
                calibrated_white_example_saturates_two_probe_bound=True,
                slow_rows=slow,local_model_coefficients_not_refitted_across_scales=True)


def shared_initial_memory_check(d):
    rows=[]
    for epsilon in (.05,.1,.2):
        shared_plus=characteristic(d,2*epsilon)
        shared_minus=characteristic(d,0.)
        fresh_plus=characteristic(d,epsilon)**2
        fresh_minus=characteristic(d,epsilon)*characteristic(d,-epsilon)
        assert abs(shared_minus-1)<1e-13 and abs(fresh_minus)<1
        rows.append(dict(epsilon=epsilon,shared_plus=encoded(shared_plus),
                         shared_compensated=encoded(shared_minus),
                         independently_redrawn_plus=encoded(fresh_plus),
                         independently_redrawn_compensated=encoded(fresh_minus),
                         compensated_coherence_error=float(abs(shared_minus-fresh_minus))))
    gaussian=[]
    for tau in (.2,.1,.05):
        exact=characteristic(d,tau)
        approx=np.exp(-1j*d['mean']*tau-.5*d['variance']*tau*tau)
        error=float(abs(exact-approx))
        gaussian.append(dict(tau=tau,characteristic_error=error,
                             error_over_tau_cubed=error/tau**3))
    assert d['kappa3']>1e-5
    assert abs(gaussian[-1]['error_over_tau_cubed']-abs(d['kappa3'])/6)<.002
    # Existing mixed conserved source matrix is retained, not a new theorem.
    conserved=old.response(beta=BETA)['D']
    area=np.array([1.,.7])
    shared=float(area@conserved@area)
    wrong=float((area*area)@np.diag(conserved))
    assert abs(shared-wrong)>1e-3
    return dict(shared_energy_sample_is_a_spectral_representation_not_an_added_bath=True,
                rows=rows,Gaussian_constant_source_matches_only_up_to_order_two=gaussian,
                third_cumulant_over_six=d['kappa3']/6,
                inherited_joint_conserved_matrix=conserved.tolist(),
                mixed_source_area=area.tolist(),joint_conserved_variance=shared,
                deleting_mixed_covariance_variance=wrong)


def run():
    d=data()
    dependencies=('research_note_600.md','research_note_601.md','research_note_602.md',
                  'research_note_603.md','research_note_620.md','research_note_622.md',
                  'research_note_623.md','research_note_625.md',
                  'joint_quantum_response_matching.py','joint_fermion_gauss_completion.py')
    return dict(round=626,tests_run=3,failures=0,errors=0,
                spectrum_and_actual_process=spectral_check(d),
                source_only_locality=locality_check(d),
                retained_common_initial_data=shared_initial_memory_check(d),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                                   for n in dependencies},
                scope='Analytic full fixed-graph Gauss Gibbs uniform-lapse sector, with no intermediate instruments: nonzero conserved variance obstructs source-only temporally local quadratic influence matching; shared energy spectral data restores this sector exactly. Numerics use unchanged frozen 32-mode CAR masses and neutral Fock factor, not the full Gauss spectrum. Does not rule out local microscopic or Wilsonian EFT with retained fields, nonlocal influence terms, different declared states or windows, and does not generate gravity.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
