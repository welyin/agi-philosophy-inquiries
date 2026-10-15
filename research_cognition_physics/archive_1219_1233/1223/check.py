"""1223: exact rational continuum bound and same-524-material diagnostics.

No saved historical output is modified. PDE/spectral numbers are diagnostics;
the exact Fraction inequalities, not grids or floats, certify the 3.32 constant.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]


def rational_certificate():
    pi_low, pi_high = F(157,50), F(22,7)
    tau, lam, T, nu = F(1,4), F(1,20), F(11,20), F(1,100)
    kap = lam*T*T
    root2_upper, H_upper = F(283,200), F(29,5)
    assert root2_upper**2 > 2
    assert H_upper**2 > F(32,3)*pi_high
    cq_upper = tau*lam*root2_upper*T*H_upper/(1-kap)**2
    cq_loose = F(59,1000)
    assert cq_upper < cq_loose
    c0_upper = F(113,200)
    # C0 = 3/4 (4/(3 pi))^(1/3).
    assert (F(4,3)*c0_upper)**3 > F(4,3)/pi_low
    mean2_upper = F(839,1000)
    assert F(4,15)*pi_high < mean2_upper
    amp_upper = F(351,250)
    assert amp_upper**2 > mean2_upper+2*c0_upper
    denominator_lower = F(1,40)
    assert denominator_lower**2 < 2*pi_low*nu**2
    coeff = cq_loose*amp_upper/denominator_lower
    asserted = F(83,25)
    assert coeff < asserted
    assert F(925,3)*pi_high < 1000
    assert F(1882,12)*pi_high < 500
    return dict(kappa=str(kap),Q_coefficient_rational_upper=str(cq_upper),
        Q_coefficient_loose_upper=str(cq_loose),variance_C0_upper=str(c0_upper),
        mean_squared_coefficient_upper=str(mean2_upper),
        combined_amplitude_upper=str(amp_upper),
        TV_coefficient_intermediate_upper=str(coeff),TV_coefficient_upper=str(asserted),
        fine_delta_rho='1/10000000',fine_TV_upper=str(asserted/F(10000000)),
        reference_excitation_energy_upper=1000,reference_mean_number_upper=500,
        exact_rational_checks=True,pi_bounds_are_classical_analytic_inputs=True)


def constants():
    V=4*math.pi/3; H=math.sqrt(32*math.pi/3)
    kap=.05*.55**2
    cq=.25*.05*math.sqrt(2)*.55*H/(1-kap)**2
    c0opt=.75*(V/(2*math.pi**2))**(1/3)
    c0=.75*(V/math.pi**2)**(1/3)
    c1=3*(V/(4*math.pi**2))**(2/3)
    tv=cq*math.sqrt(V/5+2*c0)/(math.sqrt(2*math.pi)*.01)
    # Independent algebraic stationary points of the low/high split bounds.
    L=(2*math.pi**2/V)**(1/3)
    assert abs(V*L*L/(8*math.pi**2)+1/(2*L)-c0opt)<1e-14
    Lth=(4*math.pi**2/V)**(1/3)
    assert abs(V*Lth/(2*math.pi**2)+1/Lth**2-c1)<1e-14
    return dict(spatial_volume=V,h_L1_L2_upper=H,kappa=kap,Q_coefficient=cq,
        vacuum_split_optimal=c0opt,vacuum_split_used=c0,thermal_extra_coefficient=c1,
        TV_coefficient=tv,figures_are_float_diagnostics=True)


def gaussian_checks():
    # Canonical vacuum Q,P are noncommuting; their symmetric covariance is I/2.
    # Rotated linear combinations test the auxiliary classical covariance only.
    rows=[]
    for angle in (.1,.7,1.3):
        x=np.array([1.,0.]); y=np.array([math.cos(angle),math.sin(angle)])
        S=.5*np.array([[x@x,x@y],[y@x,y@y]])
        eig=np.linalg.eigvalsh(S)
        difference=.5*float((x-y)@(x-y))
        assert eig.min()>0 and abs(np.array([1.,-1.])@S@np.array([1.,-1.])-difference)<1e-14
        rows.append(dict(angle=angle,symmetric_covariance_eigenvalues=eig.tolist(),
                         quantum_difference_variance=difference))
    maximum_excess=0.
    for shift in np.linspace(0,5,101):
        exact=math.erf(shift/(2*math.sqrt(2)))
        upper=shift/math.sqrt(2*math.pi)
        maximum_excess=max(maximum_excess,exact-upper)
        assert exact<=upper+1e-15
    return dict(noncommuting_QP_examples=rows,Gaussian_shift_envelope_max_excess=maximum_excess,
        auxiliary_classical_coupling_not_joint_quantum_measurement=True,
        general_nonGaussian_states_covered=False)


def compact_preparation():
    z,w=np.polynomial.legendre.leggauss(160)
    r=1.75+.25*z; weights=.25*w
    x=2*(2-r)
    # Stable logistic form of the smooth compact plateau.
    exponent=1/x-1/(1-x)
    s=np.empty_like(x)
    positive=exponent>=0
    temp=np.exp(-exponent[positive]); s[positive]=temp/(1+temp)
    temp=np.exp(exponent[~positive]); s[~positive]=1/(1+temp)
    derivative=s*(1-s)*(1/x**2+1/(1-x)**2)
    radial=-2*derivative
    energy=2*math.pi*float(weights@(r*r*radial*radial))
    norm2=4*math.pi*1.5**3/3+4*math.pi*float(weights@(r*r*s*s))
    assert max(abs(radial))<10 and energy<1000
    assert (norm2+2*energy)/4<500
    assert 1+.25<1.5
    return dict(Cauchy_plateau_radius=1.5,Cauchy_support_radius=2.,
        maximum_past_domain_radius=1.25,excitation_energy_quadrature_diagnostic=energy,
        L2_squared_quadrature_diagnostic=norm2,
        number_upper_from_quadrature_diagnostic=(norm2+2*energy)/4,
        independent_analytic_energy_upper=1000,independent_analytic_number_upper=500,
        global_affine_reference_used=False,normal_Fock_coherent_preparation=True)


def load_legacy():
    # Exact old entry points in the migrated layout, imported without invoking main.
    direction_path=BASE/'archive_467_530/523/thermal_direction_bridge_model.py'
    spec=importlib.util.spec_from_file_location('thermal_direction_bridge_model',direction_path)
    direction=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=direction;spec.loader.exec_module(direction)
    path=BASE/'archive_467_530/524/compact_probe_reference_model.py'
    spec=importlib.util.spec_from_file_location('round524_readonly',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def same_model_diagnostics(c):
    old=load_legacy();original=old.LAMBDA
    datasets=[]
    try:
        for delta in (0.,.001,1e-7):
            old.LAMBDA=original*(1-delta)
            data=old.solve(.025)
            datasets.append((delta,data))
    finally:
        old.LAMBDA=original
    zero=datasets[0][1]
    k,kw=old.direction.quad(192,0.,20.)
    hf=old.transform(zero,k,0.,'h')
    rows=[]
    for delta,data in datasets:
        qr=old.transform(data,k,0.)
        qm=old.transform(data,k,old.MU2)
        # Exact vacuum weights, but finite-k integral and PDE remain diagnostic.
        var_R=float((kw*k/(4*math.pi**2))@abs(qr)**2)
        var_M=float((kw*k*k/(4*math.pi**2*np.sqrt(k*k+old.MU2)))@abs(qm)**2)
        var_chi=float((kw*k/(4*math.pi**2))@abs(hf-qr)**2)
        variance=var_R/5+4*var_M/5+var_chi+.01**2
        mean=-data['gain']/math.sqrt(5)
        prob=.5*(1+math.erf(mean/math.sqrt(2*variance)))
        dq=data['q']-zero['q']
        spatial=4*math.pi*data['r']**2*old.trap_weights(data['r'])
        Q=float(old.trap_weights(data['t'])@np.sqrt((dq*dq)@spatial))
        if delta:
            assert Q<c['Q_coefficient']*delta
        else:
            assert Q==0
        assert abs(data['summary']['forward_probe_mean']-mean)<1e-10
        rows.append(dict(delta_rho_upper=delta,effective_coupling_for_legacy_evaluation=original*(1-delta),
            gain=data['gain'],mean=mean,variance_truncated=variance,
            fixed_event_Y_ge_0_probability_truncated=prob,Q_diagnostic=Q,
            analytic_Q_upper=c['Q_coefficient']*delta,analytic_TV_rational_upper=3.32*delta,
            forward_advanced_residual=data['summary']['forward_advanced_reciprocity_residual']))
    for row in rows:
        row['fixed_event_difference_truncated']=abs(row['fixed_event_Y_ge_0_probability_truncated']-rows[0]['fixed_event_Y_ge_0_probability_truncated'])
        assert row['fixed_event_difference_truncated']<=row['analytic_TV_rational_upper']+1e-15
    return dict(spatial_step=.025,spectral_endpoint=20.,spectral_order=192,
        inherited_PDE_modified_on_disk=False,rows=rows,
        numerical_probability_is_not_interval_or_full_continuum_certification=True,
        finite_k_endpoint_is_not_physical_UV_cutoff=True)


def run():
    rat=rational_certificate();const=constants();gau=gaussian_checks();prep=compact_preparation()
    diag=same_model_diagnostics(const)
    sources=[HERE/'check.py',HERE/'proof.md',
        BASE/'archive_467_530/research_note_524.md',
        BASE/'archive_467_530/524/compact_probe_reference_model.py',
        BASE/'archive_467_530/524/compact_probe_reference_results.json',
        BASE/'archive_467_530/523/thermal_direction_bridge_model.py']
    return dict(schema='round1223_continuum_Green_probability_v1',status='PASS',check_groups=4,
        group_names=['continuous_resolvent_and_rational_envelope','Gaussian_probability_and_covariance',
                     'normal_compact_preparation_and_resource_bounds','same_524_PDE_and_fixed_event_diagnostics'],
        analytic_rational_certificate=rat,floating_constants=const,Gaussian_checks=gau,
        preparation=prep,same_material_diagnostics=diag,
        scope=dict(full_single_readout_distribution_bound=True,nonperturbative=True,
            normally_prepared_finite_excitation_energy_example=True,actual_old_Gaussian_readout_reused=True,
            continuous_history_bound_independent_of_mesh=True,quantum_source_history_internalized=False,
            co_moving_clock_constructed=False,full_instrument_poststate_distance_controlled=False,
            directly_fills_1222_all_material_errors=False,general_nonGaussian_state_bound=False,
            new_adopted_axioms=0,empirical_calibrations=0,TG2_complete=False,stage_complete=False),
        source_sha256={p.relative_to(BASE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();result=run();path=HERE/'results.json'
    if args.write_results:
        with path.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        saved=json.loads(path.read_text(encoding='utf8'))
        assert saved==result,'Read-only recomputation differs from saved results'
    print(json.dumps(dict(status='PASS',check_groups=4,TV_coefficient_upper='3.32',
        fixed_small_history_difference='1e-7',full_TV_upper='3.32e-7',
        source_clock_internalized=False),ensure_ascii=False))


if __name__=='__main__':
    main()
