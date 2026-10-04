"""713: actual quotient-group loop instruments and their full-H injection.

Spatially averaging disjoint loop characters is a new specified readout,
not an autonomous detector, a flow field, or a continuum quantum-state proof.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_region_energy_gluing as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_holonomy_readout_scale_results.json'
GAUGE=original.gauge
BC=GAUGE.parameters()['b']
BW=float(BC[1]);B0=float(BC[2]);CSTAR=max(BW/4,9*B0)
CASIMIR=3*BW/4+9*B0
IDENTITY=(np.eye(3,dtype=complex),np.eye(2,dtype=complex),1.+0j)


def multiply(items):
    result=IDENTITY
    for g in items:result=original.product(result,g)
    return result


def f(g):
    return float(np.trace(original.rep(g,'L')).real/2)


def gamma(g):
    a=float(np.trace(g[1]).real/2);phase=3*np.angle(g[2])
    return BW*(1-a*a)*np.cos(phase)**2/4+9*B0*a*a*np.sin(phase)**2


def shift(g,kind,a,t):
    if kind=='weak':
        coords=np.zeros(3);coords[a]=t
        return g[0],GAUGE.group_exp(coords,2)@g[1],g[2]
    return g[0],g[1],np.exp(1j*t)*g[2]


def topology_check():
    rng=np.random.default_rng(7131)
    loops=[[original.sample(rng) for _ in range(5)] for _ in range(3)]
    before=[f(multiply(loop)) for loop in loops]
    errors=[]
    center=(np.exp(2j*np.pi/3)*np.eye(3),-np.eye(2),np.exp(1j*np.pi/3))
    for loop,value in zip(loops,before):
        frames=[original.sample(rng) for _ in loop]
        transformed=[original.product(original.product(frames[i],g),original.inverse(frames[(i+1)%len(loop)])) for i,g in enumerate(loop)]
        errors.append(abs(f(multiply(transformed))-value))
        for i in range(len(loop)):
            changed=list(loop);changed[i]=original.product(center,changed[i])
            errors.append(abs(f(multiply(changed))-value))
    assert max(errors)<1e-12
    rows=[]
    for n in (4,8,16,32,64):
        W=np.diag(np.exp(1j*np.pi/n*np.array([1,-1])))
        link=(np.eye(3),W,1.)
        end=multiply([link]*n)
        raw=original.matter.representation(*link)
        endpoint=original.matter.representation(*end)
        rows.append(dict(edges=n,loop_readout=f(end),single_32_mode_link_distance=float(np.linalg.norm(raw-np.eye(32))),
                         full_32_mode_holonomy_distance=float(np.linalg.norm(endpoint-np.eye(32))),
                         plaquette_error=float(np.linalg.norm(W@np.eye(2)@W.conj().T-np.eye(2)))))
        assert abs(f(end)+1)<1e-13 and abs(rows[-1]['full_32_mode_holonomy_distance']-8)<1e-12
    assert all(rows[i+1]['single_32_mode_link_distance']<rows[i]['single_32_mode_link_distance'] for i in range(len(rows)-1))
    return dict(original_L_quotient_character=True,local_Gauss_and_Z6_errors=errors,flat_center_refinement=rows,
        noncontractible_loop_not_instanton_charge=True,original_699_kernel_not_transported=True)


def readout_check():
    rng=np.random.default_rng(7132);n=4;m=3;eps=.73;sigma=.12;eta=.61;step=1e-5
    loops=[[original.sample(rng) for _ in range(n)] for _ in range(m)]
    def amplitudes(data):
        val=np.mean([f(multiply(loop)) for loop in data])
        return np.sqrt((1+np.array([-1,1])*eta*val)/2)
    values=[f(multiply(loop)) for loop in loops]
    analytic=n/eps*np.exp(-2*sigma)*eta**2/(4*m*m)*(sum(gamma(multiply(loop)) for loop in loops))/(1-eta**2*np.mean(values)**2)
    numerical=0.
    for l in range(m):
        for e in range(n):
            for kind,dim,b in (('weak',3,BW),('abelian',1,B0)):
                for a in range(dim):
                    plus=[list(loop) for loop in loops];minus=[list(loop) for loop in loops]
                    plus[l][e]=shift(loops[l][e],kind,a,step);minus[l][e]=shift(loops[l][e],kind,a,-step)
                    derivative=(amplitudes(plus)-amplitudes(minus))/(2*step)
                    numerical+=b/eps*np.exp(-2*sigma)*float(derivative@derivative)
    error=abs(numerical-analytic);assert error<1e-10
    norm=n/(eps*m)*np.exp(-2*sigma)*eta**2*CSTAR/4
    maximizer=(np.eye(3),np.diag([1j,-1j]),1.) if BW/4>=9*B0 else (np.eye(3),np.eye(2),np.exp(1j*np.pi/6))
    exact=n/(eps*m)*np.exp(-2*sigma)*eta**2*gamma(maximizer)/4/(1-eta**2*f(maximizer)**2)
    assert abs(exact-norm)<1e-15
    # Bilinear bound that yields the exact essential supremum, not a mesh estimate.
    x=rng.random(100);y=rng.random(100)
    gam=BW*(1-x)*y/4+9*B0*x*(1-y)
    assert np.all(gam<=CSTAR*(1-x*y)+1e-15)
    derivative=(analytic*np.exp(-2*step)-analytic*np.exp(2*step))/(2*step)
    source_error=abs(derivative+2*analytic);assert source_error<1e-10
    return dict(loop_edges=n,parallel_loops=m,eta=eta,epsilon=eps,sigma=sigma,
        original_b_weak=BW,original_b_U1=B0,original_character_Cstar=CSTAR,
        analytic_full_H_injection_multiplier=float(analytic),matrix_directional_injection=float(numerical),
        direct_instrument_derivative_error=error,exact_operator_norm=norm,
        attained_norm_error=abs(exact-norm),same_geometry_source_error=source_error,
        scalar_multiplication_commutes_all_original_mass_and_hopping=True,
        fixed_gain_instrument_is_additional_input=True)


def scale_check():
    eta=.6;length=1.;area=.25;rows=[]
    for n in (4,8,16,32,64):
        a=length/n;m=(n//2)**2
        thin=eta**2*CSTAR*n/(4*a)
        thick=thin/m
        assert abs(m*a*a-area)<1e-15
        assert abs(thick-eta**2*CSTAR*length/(4*area))<1e-15
        rows.append(dict(edges=n,physical_spacing=a,parallel_loops=m,physical_cross_section=m*a*a,
                         thin_worst_state_budget=thin,averaged_worst_state_budget=thick,
                         flat_holonomy_effect_contrast=eta,pointwise_injection_at_both_exact_centers=0.))
    return dict(eta=eta,physical_length=length,physical_area=area,refinements=rows,
        normal_quantum_state_preparation_uniform_energy_not_proved=True,
        continuum_renormalized_Wilson_observable_not_proved=True,
        state_independent_injection_bound_not_total_resource_cost=True)


def state_energy_check():
    # Haar on SU2: a0=cos(theta), (2/pi)sin(theta)^2 dtheta.
    # The quotient-invariant integrand can be integrated on the covering group.
    def integrals(r,n):
        theta=np.pi*np.arange(1,n+1)/(n+1)
        a=np.cos(theta)[:,None];weights=(2/(n+1)*np.sin(theta)**2)[:,None]
        phase=2*np.pi*np.arange(2*n)/(2*n)
        val=a*np.cos(phase)[None,:]
        gam=BW*(1-a*a)*np.cos(phase)[None,:]**2/4+9*B0*a*a*np.sin(phase)[None,:]**2
        density=weights*np.exp(r*val)/(2*n)
        z=np.sum(density)
        return np.array([z,np.sum(density*val)/z,np.sum(density*gam)/z])
    rows=[]
    for r in (1.,3.,6.):
        low=integrals(r,48);high=integrals(r,96)
        error=float(max(abs(low-high)))
        z,mean,grad=high
        ward=abs(CASIMIR*mean-r*grad)
        energy_per_NM_over_a=r*r*grad/4
        bound_per_NM_over_a=CASIMIR**2*mean**2/(4*CSTAR)
        assert max(error,ward)<1e-11 and energy_per_NM_over_a>=bound_per_NM_over_a
        rows.append(dict(tilt=r,normalizer=float(z),signal_mean=float(mean),
            weighted_gradient_average=float(grad),exact_integration_by_parts_error=float(ward),
            two_quadrature_orders_difference=error,
            actual_electric_energy_per_NM_over_a=float(energy_per_NM_over_a),
            universal_electric_lower_bound_per_NM_over_a=float(bound_per_NM_over_a)))
    scale=[]
    m=rows[1]['signal_mean'];factor=rows[1]['actual_electric_energy_per_NM_over_a']
    for n in (4,8,16,32,64):
        a=1/n;loops=(n//2)**2
        scale.append(dict(edges=n,parallel_loops=loops,mean=m,
            exact_electric_expectation=factor*n*loops/a,
            necessary_electric_expectation=CASIMIR**2*m*m*n*loops/(4*CSTAR*a)))
    return dict(original_character_weighted_Casimir=CASIMIR,normal_Gauss_states=True,
        full_H_still_contains_all_original_terms=True,tilted_Haar_checks=rows,
        fixed_signal_refinement=scale,bound_applies_to_any_normal_finite_electric_energy_state=True,
        electric_expectation_is_not_preparation_work_above_a_common_reference=True,
        reference_subtracted_continuum_cost_not_determined=True)


def run():
    a=topology_check();b=readout_check();c=scale_check();d=state_energy_check()
    deps=('research_note_568.md','research_note_574.md','research_note_598.md','research_note_617.md',
          'research_note_637.md','research_note_699.md','research_note_700.md','research_note_701.md',
          'research_note_712.md','joint_region_energy_gluing.py','joint_quotient_gauge_completion.py',
          'round713_drafts/local_interaction_entry.md')
    return dict(round=713,tests_run=4,failures=0,errors=0,physical_holonomy=a,actual_instrument=b,common_scale=c,actual_state_energy=d,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original full finite Hamiltonian with specified loop instruments on a given 3D isotropic refinement family; exact Gauss, injection and source identities. No continuum state, autonomous device or transported699 negative form claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as out:json.dump(result,out,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
