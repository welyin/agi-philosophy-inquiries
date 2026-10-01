"""639: original curved singlet read, fixed support and shared geometry.

Collective read is a candidate mathematical instrument, not an instantaneous
causal implementation. Configurations below are original field samples, NOT
Gibbs averages or simulations of an emerging continuum.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_record_relative_entropy as prior

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_regional_read_scale_results.json'
original=prior.prior.model.original
HBAR=prior.prior.HBAR


def kappa(q):
    return np.cos(q)**2/(16*(1-np.sin(q)**2/4))


def kappa_prime(q):
    return -3*np.sin(q)*np.cos(q)/(32*(1-np.sin(q)**2/4)**2)


def quantities(phi,w,f):
    alpha=w*f/np.sum(w*f)
    q=float(np.sum(alpha*phi[...,4]))
    Kss=original.inverse(phi)[...,4,4]
    B=float(np.sum(alpha*alpha*Kss/w))
    injection=float(HBAR**2/2*kappa(q)*B)
    inv_volume=float(np.sum(alpha*alpha/w))
    bound=float(HBAR**2*original.M/32*inv_volume)
    return dict(alpha=alpha,q=q,Kss=Kss,B=B,injection=injection,
                inv_volume=inv_volume,bound=bound)


def product_rule_check():
    rng=np.random.default_rng(639);errors=[]
    for _ in range(30):
        phi=rng.normal(size=(3,5))*.3
        w=rng.uniform(.3,1.2,3);f=rng.uniform(.2,1,3)
        data=quantities(phi,w,f)
        L,dL,_=prior.instrument(data['q'])
        psi=rng.normal()+1j*rng.normal()
        grad=rng.normal(size=(3,5))+1j*rng.normal(size=(3,5))
        K=original.inverse(phi)
        before=sum(np.vdot(grad[v],K[v]@grad[v]).real/w[v] for v in range(3))
        after=0.
        for r in range(2):
            new=L[r]*grad.copy()
            new[:,4]+=psi*dL[r]*data['alpha']
            after+=sum(np.vdot(new[v],K[v]@new[v]).real/w[v] for v in range(3))
        direct=HBAR**2/2*(after-before)
        expected=data['injection']*abs(psi)**2
        errors.append(float(abs(direct-expected)))
    assert max(errors)<2e-13
    # The exact multiplication-operator norm is approached at original phi=0.
    v=np.array([[.2,-.1,.3,.1,.4],[.3,.2,-.2,.1,-.3],[.1,.2,.1,.1,.2]])
    w=np.array([.4,.7,1.1]);f=np.array([.8,1.2,.5])
    rows=[]
    for scale in (1.,.3,.1,.03,0.):
        d=quantities(scale*v,w,f)
        rows.append(dict(field_scale=scale,injection=d['injection'],
                         exact_operator_norm=d['bound'],ratio=d['injection']/d['bound']))
    assert all(rows[k+1]['ratio']>rows[k]['ratio'] for k in range(len(rows)-1))
    assert abs(rows[-1]['ratio']-1)<1e-14
    return dict(original_three_node_Dirichlet_error=max(errors),
        original_target_origin_saturation=rows,
        saturation_is_operator_norm_not_Gibbs_expectation=True)


def sample(N):
    q=original.shared_source(N)
    x,y,z=np.moveaxis(q['grid'],-1,0)
    psi=1.1+.05*np.cos(x)+.02*np.sin(y)
    w=q['eps']**3*psi**6
    # Smooth periodic profile with proper fixed support cos(x)>0.
    a=np.cos(x)
    f=np.zeros_like(a);support=a>0
    f[support]=np.exp(-1/a[support]**2)*(1+.2*np.cos(y[support]))
    return q,w,f


def scale_check():
    rows=[]
    # All N are multiples of 8, permitting a simple uniform support count.
    for N in (8,16,24,32):
        q,w,f=sample(N);d=quantities(q['phi'],w,f)
        volume=float(np.sum(w));support_volume=float(np.sum(w[f>0]))
        rows.append(dict(N=N,nodes=N**3,total_physical_volume=volume,
            support_volume=support_volume,effective_volume=1/d['inv_volume'],
            regional_read_operator_norm=d['bound'],
            maximum_point_read_operator_norm=float(HBAR**2*original.M/(32*w.min())),
            original_configuration_injection=d['injection'],
            original_configuration_read_value=d['q']))
        assert d['injection']<=d['bound']*(1+1e-13)
        assert 1/d['inv_volume']<=support_volume*(1+1e-13)
    # On |x|<=pi/4, at least N/8 layers: f>=.8e^-2 and psi>=1.03.
    box=(2*np.pi)**3
    I1_lower=box/8*1.03**6*.8*np.exp(-2)
    I2_upper=box*1.17**6*(1.2/np.e)**2
    uniform_volume_lower=I1_lower**2/I2_upper
    uniform_norm_upper=HBAR**2*original.M/(32*uniform_volume_lower)
    assert min(r['effective_volume'] for r in rows)>uniform_volume_lower
    assert max(r['regional_read_operator_norm'] for r in rows)<uniform_norm_upper
    ratio=rows[-1]['maximum_point_read_operator_norm']/rows[0]['maximum_point_read_operator_norm']
    assert abs(ratio-64)<1e-10
    assert abs(rows[-1]['effective_volume']/rows[-2]['effective_volume']-1)<.01
    return dict(rows=rows,analytic_uniform_effective_volume_lower=float(uniform_volume_lower),
        analytic_uniform_injection_norm_upper=float(uniform_norm_upper),
        fixed_beta_2_unconditional_relative_entropy_upper=float(2*uniform_norm_upper),
        chosen_fixed_support_profile_is_extra_input=True,
        no_thermal_averages_or_cross_graph_state_convergence_computed=True)


def geometry_check():
    q,w,f=sample(8);x,y,z=np.moveaxis(q['grid'],-1,0)
    phi=q['phi'];d=quantities(phi,w,f)
    rows=[]
    for label,c in (('nonuniform',.2+np.sin(x+y)+.4*np.cos(z)),
                    ('uniform',np.ones_like(x))):
        alpha=d['alpha'];cbar=float(np.sum(alpha*c))
        dalpha=6*alpha*(c-cbar)
        dq=float(np.sum(dalpha*phi[...,4]))
        weights=alpha**2*d['Kss']/w
        Bprime=float(np.sum(weights*(6*c-12*cbar)))
        total=HBAR**2/2*(kappa_prime(d['q'])*dq*d['B']+kappa(d['q'])*Bprime)
        source_fixed=HBAR**2/2*kappa(d['q'])*float(np.sum(-6*c*weights))
        operation_change=HBAR**2/2*(kappa_prime(d['q'])*dq*d['B']+
                                   2*kappa(d['q'])*float(np.sum(alpha*dalpha*d['Kss']/w)))
        step=2e-6
        plus=quantities(phi,w*np.exp(6*step*c),f)['injection']
        minus=quantities(phi,w*np.exp(-6*step*c),f)['injection']
        fd=(plus-minus)/(2*step)
        error=abs(fd-total)
        assert error<1e-11 and abs(total-source_fixed-operation_change)<2e-17
        if label=='uniform':
            assert abs(dq)<1e-13
            assert abs(total+6*d['injection'])<1e-16
            assert abs(operation_change)<1e-16
        else:assert abs(operation_change)>1e-8
        rows.append(dict(geometry_direction=label,read_value_derivative=dq,
            total_geometry_derivative=float(total),
            original_source_change_at_fixed_instrument=float(source_fixed),
            extra_instrument_variation=float(operation_change),
            finite_difference_error=float(error)))
    return dict(rows=rows,actual_original_field_metric_volume_and_profile=True,
        original_uniform_conformal_identity_recovered=True,
        varying_instrument_not_confused_with_fixed_source=True)


def run():
    deps=('research_note_593.md','research_note_598.md','research_note_603.md',
          'research_note_625.md','research_note_633.md','research_note_638.md',
          'joint_record_relative_entropy.py','joint_curved_quantum_source.py',
          'joint_full_spatial_metric.py','joint_geometry_work_noise.py')
    return dict(round=639,tests_run=3,failures=0,errors=0,
        original_curved_read=product_rule_check(),
        physical_support_and_scale=scale_check(),shared_geometry=geometry_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_full_model_record_injection_proved_analytically=True,
            original_H_not_replaced_by_free_scalar=True,
            collective_singlet_read_is_an_additional_operation_family=True,
            fixed_physical_support_controls_uniform_information_upper_bound=True,
            source_response_and_operation_variation_separated=True,
            no_instantaneous_causal_implementation_or_nontrivial_signal_lower_bound=True,
            no_continuum_state_map_area_entropy_or_GR_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
