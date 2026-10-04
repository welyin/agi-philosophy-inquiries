"""Preliminary 557 link-readout diagnostics; no completed round is registered.

t is an instrument resolution parameter, not time evolution under the full
matter-link Hamiltonian. D is only the link square-root instrument's cost.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'link_heat_kernel_probe_results.json'


def kernel(theta,t,terms):
    theta=np.asarray(theta)
    n=np.arange(1,terms+1)[:,None]
    angle=theta[None,:]
    coefficients=n*np.exp(-t*(n*n-1)/4)
    char=np.sin(n*angle)/np.sin(angle)
    derivative=(n*np.cos(n*angle)*np.sin(angle)-np.sin(n*angle)*np.cos(angle))/np.sin(angle)**2
    return np.sum(coefficients*char,axis=0),np.sum(coefficients*derivative,axis=0)


def evaluate(t,quadrature,terms):
    x,w=np.polynomial.legendre.leggauss(quadrature)
    theta=np.pi*(x+1)/2
    haar=w*np.sin(theta)**2
    k,kprime=kernel(theta,t,terms)
    assert np.min(k)>0
    return dict(t=t,normalization=float(haar@k),
                fundamental_attenuation=float(haar@(k*np.cos(theta))),
                predicted_attenuation=math.exp(-3*t/4),
                link_D=float(haar@(kprime*kprime/k)/16),
                minimum_sampled_kernel=float(min(k)),
                singlet_survival_for_constant_link_input=float(haar@np.sqrt(k))**2)


def run():
    rows=[]
    for t in (.5,.75,1.,2.,4.):
        row=evaluate(t,180,55);other=evaluate(t,260,80)
        assert abs(row['normalization']-1)<1e-12
        assert abs(row['fundamental_attenuation']-row['predicted_attenuation'])<1e-12
        residual=abs(row['link_D']-other['link_D'])
        assert residual<1e-9
        assert 0<row['singlet_survival_for_constant_link_input']<1
        row['independent_quadrature_difference']=residual
        rows.append(row)
    assert all(rows[i]['link_D']>rows[i+1]['link_D'] for i in range(len(rows)-1))
    return dict(status='preliminary_557_not_a_completed_round',rows=rows,
                model_input='single_dynamic_SU2_link_with_Casimir_j_jplus1',
                scope=dict(instrument_t_is_not_dynamical_time=True,
                           moments_for_probabilities_only=True,
                           square_root_group_reading_not_Gauss_preserving=True,
                           D_excludes_matter_instrument_energy=True,
                           no_small_t_asymptotic_claim_verified_numerically=True,
                           no_unified_model_completion=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists():
        assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
