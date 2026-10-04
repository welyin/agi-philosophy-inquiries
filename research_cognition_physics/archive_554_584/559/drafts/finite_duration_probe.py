"""Preliminary 559: finite rectangular probe in a fixed matter fiber.

The duration, switching and independent neutral pointer are explicit inputs.
No autonomous full-matter detector or universal low-energy instrument is proved.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'finite_duration_probe_results.json'


def data(alpha,b,N):
    n=np.arange(N+1,dtype=float)
    C=np.diag(np.ones(N)/2,1)+np.diag(np.ones(N)/2,-1)
    K=b*np.diag(n*(n+2)/4)-alpha*C
    E,V=np.linalg.eigh(K)
    return K,alpha*(np.eye(N+1)-C),V[:,0],float(E[0])


def propagate(p,tau,K,F,psi):
    E,V=np.linalg.eigh(K+p*F/tau)
    phase=np.exp(-1j*tau*E)
    source=V.T@psi
    output=V@(phase*source)
    Fe=V.T@F@V
    average=(E[:,None]+E[None,:])/2
    difference=(E[:,None]-E[None,:])/2
    derivative_e=-1j*Fe*np.exp(-1j*tau*average)*np.sinc(tau*difference/np.pi)
    derivative=V@(derivative_e@source)
    return output,derivative


def evaluate(alpha,b,tau,nu=.5,N=20,nq=40):
    K,F,psi,lam=data(alpha,b,N)
    nodes,w=np.polynomial.hermite.hermgauss(nq); w=w/math.sqrt(math.pi)
    sigma_p=1/(2*nu)
    mean=0.;second=nu*nu;energy=0.;outside=0.;normalization=0.
    for p,weight in zip(math.sqrt(2)*sigma_p*nodes,w):
        out,der=propagate(p,tau,K,F,psi)
        normalization+=weight*np.vdot(out,out).real
        mean+=weight*(1j*np.vdot(out,der)).real
        second+=weight*np.vdot(der,der).real
        energy+=weight*(np.vdot(out,K@out).real-lam)
        outside+=weight*(1-abs(np.vdot(psi,out))**2)
    point=.17;eps=1e-5
    _,der=propagate(point,tau,K,F,psi)
    fd=(propagate(point+eps,tau,K,F,psi)[0]-propagate(point-eps,tau,K,F,psi)[0])/(2*eps)
    residual=float(np.linalg.norm(der-fd))
    assert abs(normalization-1)<1e-12 and residual<1e-7
    assert energy>-1e-11 and second-mean*mean>=nu*nu-1e-10
    return dict(alpha=alpha,b=b,duration=tau,noise_variance=nu*nu,
        actual_mean=float(mean),actual_variance=float(second-mean*mean),
        matched_ground_mean=float(psi@F@psi),
        source_electric_plus_edge_energy_increase=float(energy),
        probability_outside_initial_fast_ground=float(outside),
        pointer_derivative_finite_difference_residual=residual)


def run():
    rows=[evaluate(.3,b,tau) for b in (1.,10.) for tau in (.01,.2,1.,5.,20.,100.)]
    for row in rows:
        other=evaluate(.3,row['b'],row['duration'],N=28,nq=56)
        row['independent_cutoff_quadrature_difference']=max(
            abs(row[key]-other[key]) for key in ('actual_mean','actual_variance','source_electric_plus_edge_energy_increase'))
        assert row['independent_cutoff_quadrature_difference']<1e-9
    return dict(status='preliminary_559_not_a_completed_round',rows=rows,
        scope=dict(frozen_matter_fiber_only=True,
            rectangular_finite_time_coupling_is_new_control_input=True,
            pointer_free_evolution_not_included=True,
            actual_pointer_moments_not_inferred_only_from_effective_H=True,
            no_analytic_uniform_error_bound_completed=True,
            no_full_source_or_autonomous_apparatus_claim=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists(): assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],cases=len(result['rows']),
          examples=[result['rows'][0],result['rows'][5],result['rows'][6],result['rows'][-1]]),ensure_ascii=False))

