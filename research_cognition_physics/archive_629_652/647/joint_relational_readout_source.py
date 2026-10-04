"""647: original573 material reference, original readout and moving location.

Classical relational observables; numerics test a spatial diffeomorphism of
the same573 initial data. No quantum instrument or BV quantization claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gravity_material_coordinates as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_relational_readout_source_results.json'
PROFILE=np.array([3.,-2.,1.5,.8])


def transformed(N,epsilon):
    old=original.fields(N);q=old['q'];x,y,z=np.moveaxis(q['grid'],-1,0)
    zeta=1+.2*np.cos(2*z);yp=y+epsilon*zeta
    # Exact shift of the collocation trigonometric polynomial along y for each z.
    modes=np.fft.fftfreq(N,1/N)
    phase=np.exp(1j*modes[None,:,None]*epsilon*zeta)
    psi=np.fft.ifft(np.fft.fft(old['psi'],axis=1)*phase,axis=1).real
    hs,ss=old['hstar'],old['sstar'];ah=.05*hs;bs=.06*ss
    h=hs*(1+.05*np.sin(x+yp));s=ss*(1+.06*np.cos(yp))
    ph=-.0056/ah*np.cos(x+yp);ps=.025*np.cos(x)-.0056/bs*np.sin(yp)
    F=2-(h*h+s*s)/6;dot=h*ph+s*ps
    vh=psi**-6*F*(ph-h*dot/12);vs=psi**-6*F*(ps-s*dot/12)
    hy=ah*np.cos(x+yp);sy=-bs*np.sin(yp)
    shear=-.4*epsilon*np.sin(2*z)
    dh=np.stack([hy,hy,shear*hy],axis=-1)
    ds=np.stack([np.zeros_like(s),sy,shear*sy],axis=-1)
    # Full pulled-back inverse gamma, including the new yz off-diagonal entries.
    gi=np.zeros((N,N,N,3,3))
    gi[...,0,0]=1;gi[...,1,1]=1+shear**2;gi[...,2,2]=1
    gi[...,1,2]=gi[...,2,1]=-shear
    gi*=psi[...,None,None]**-4
    A=np.einsum('...i,...ij,...j->...',dh,gi,dh)-vh*vh
    B=np.einsum('...i,...ij,...j->...',ds,gi,ds)-vs*vs
    expectedA=psi**-4*2*hy**2-vh*vh
    expectedB=psi**-4*sy**2-vs*vs
    pullback_error=max(float(np.max(abs(A-expectedA))),float(np.max(abs(B-expectedB))))
    menu=np.stack([h,s,A,B],axis=-1)
    profile=np.exp(menu@PROFILE)
    mu=psi**6
    effect=.5+.25*np.sin(s)
    return dict(menu=menu,profile=profile,mu=mu,effect=effect,psi=psi,zeta=zeta,
                h_y=hy,s_y=sy,pullback_error=pullback_error,old=old)


def region_value(fields,profile=None):
    w=fields['mu']*(fields['profile'] if profile is None else profile)
    return float(np.sum(w*fields['effect'])/np.sum(w))


def original_reference_check():
    f=original.fields(32);q=f['q'];base=transformed(32,0.)
    target=np.stack([q['f']['h'],q['f']['s'],f['rh'],f['rs']],axis=-1)
    err=float(np.max(abs(base['menu']-target)));assert err<2e-14
    # These are references to the original constrained source, not new initial data.
    idx=(0,8,4);s=float(target[idx][1]);sy=float(f['ds'][idx][1])
    fixed_coordinate=float(.25*np.cos(s)*sy)
    for eps in (1e-3,5e-4):
        plus=transformed(32,eps);minus=transformed(32,-eps)
        fd=(plus['effect'][idx]-minus['effect'][idx])/(2*eps)
        assert abs(fd-fixed_coordinate)<2e-9
    # Fixed relational X^1=s gives e(X^1), so its change includes the location term.
    location=-.25*np.cos(s)*sy
    assert abs(fixed_coordinate+location)<1e-16 and abs(fixed_coordinate)>.001
    errors=[]
    for eps in (-.19,.11,.23):
        shifted=transformed(32,eps)
        assert shifted['pullback_error']<1e-14
        errors.append(shifted['pullback_error'])
    return dict(original_menu_match_error=err,spatial_point=[0.,np.pi/2,np.pi/4],
                original_s=s,original_s_y=sy,
                fixed_coordinate_effect_response=fixed_coordinate,
                fixed_relational_label_location_response=location,
                fixed_relational_label_total_response=fixed_coordinate+location,
                finite_diffeomorphism_tensor_errors=errors,
                on_shell_reference_existence_from573_not_reproved=True,
                last_two_clock_time_derivatives_not_invented=True,
                no_statement_that_quantum_clock_fluctuations_vanish=True)


def moving_region_check():
    rows=[]
    for N in (24,32,48):
        base=transformed(N,0.);old=base['old'];zeta=base['zeta']
        w=base['mu']*base['profile'];w=w/np.sum(w)
        P=region_value(base);center=base['effect']-P
        dx=np.stack([base['h_y'],base['s_y'],old['grh'][...,1],old['grs'][...,1]],axis=-1)
        dx*=zeta[...,None]
        measure=6*zeta*original.geo.derivative(base['psi'],1)/base['psi']
        effect=.25*np.cos(base['menu'][...,1])*dx[...,1]
        location=dx@PROFILE
        components=np.array([np.sum(w*effect),np.sum(w*center*measure),np.sum(w*center*location)])
        total=float(np.sum(components));frozen=float(np.sum(components[:2]))
        steps=[]
        for step in (2e-3,1e-3):
            plus=transformed(N,step);minus=transformed(N,-step)
            full_fd=(region_value(plus)-region_value(minus))/(2*step)
            wrong_fd=(region_value(plus,base['profile'])-region_value(minus,base['profile']))/(2*step)
            steps.append(dict(step=step,full_response=full_fd,frozen_location_response=wrong_fd,
                              frozen_response_difference=abs(wrong_fd-frozen)))
        finite=[]
        for eps in (-.19,.11,.23):
            shifted=transformed(N,eps)
            finite.append(dict(epsilon=eps,complete_difference=region_value(shifted)-P,
                frozen_location_difference=region_value(shifted,base['profile'])-P))
        rows.append(dict(N=N,base_region_readout=P,
            effect_volume_location_responses=components.tolist(),total_pure_gauge_response=total,
            omitted_location_response=frozen,steps=steps,finite_changes=finite))
    final=rows[-1]
    assert abs(final['total_pure_gauge_response'])<2e-12
    assert abs(final['omitted_location_response'])>1e-6
    assert final['steps'][-1]['frozen_response_difference']<2e-9
    assert max(abs(x['complete_difference']) for x in final['finite_changes'])<2e-12
    frozen_gap=max(abs(x['frozen_location_difference']) for x in final['finite_changes'])
    numerical_scale=max(abs(x['complete_difference']) for x in final['finite_changes'])
    numerical_scale+=final['steps'][-1]['frozen_response_difference']+np.finfo(float).eps
    assert frozen_gap>1000*numerical_scale
    return dict(profile_coefficients=PROFILE.tolist(),rows=rows,
                region_is_smooth_weighted_initial_slice_not_a_quantum_instrument=True,
                periodic_spatial_diffeomorphism_generates_same_geometric_initial_data=True,
                spacetime_Ward_identity_is_analytic_not_an_evolution_simulation=True,
                no_global_unique_material_chart_claim=True)


def run():
    deps=('research_note_548.md','research_note_573.md','research_note_624.md','research_note_639.md',
          'research_note_646.md','joint_gravity_material_coordinates.py',
          'joint_gravity_material_coordinates_results.json','joint_gauss_einstein_initial_data.py')
    return dict(round=647,tests_run=2,failures=0,errors=0,
                original_reference_and_readout=original_reference_check(),
                same_field_region_sources=moving_region_check(),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(original573_classical_material_reference_and_constraint_branch=True,
                           original_readout_function_reused_but_no_quantum_instrument_map=True,
                           full_classical_location_variation_and_limited_spatial_diagnostics=True,
                           no_new_basic_clocks_or_cognition_axioms=True,
                           no_BV_quantization_or_full_Gauss_GR_derivation_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
