"""Finite-window state-dependent response: 1+1 derivative-coupling calibration.

The four-dimensional HH/Unruh comparison is analytic in the report. This
script neither computes 4D radial modes nor solves semiclassical backreaction.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'finite_horizon_response_results.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def rule(n,left,right):
    x,w=leggauss(n)
    return (right+left)/2+(right-left)*x/2,(right-left)*w/2


def switch(x,L):
    """C-infinity compact support, normalized to peak one."""
    x=np.asarray(x);out=np.zeros_like(x,dtype=float)
    mask=np.abs(x)<L
    out[mask]=np.exp(1-1/(1-(x[mask]/L)**2))
    return out


def smooth_kernel(z):
    """D(s)/alpha**2, z=alpha*s; stable coincidence series."""
    z=np.asarray(z,dtype=float);out=np.empty_like(z)
    small=np.abs(z)<.05;u=z[small]**2
    out[small]=(1/12-u/240+u*u/6048-u*u*u/172800)/(2*np.pi)
    v=z[~small]
    out[~small]=(1/(v*v)-1/(4*np.sinh(v/2)**2))/(2*np.pi)
    return out


def time_difference(n,L,e):
    x,w=rule(n,-L,L);weighted=w*switch(x,L)
    z=x[:,None]-x[None,:]
    return float(weighted@(smooth_kernel(z)*np.cos(e*z))@weighted)


def spectral_difference(n,L,e,Q=10.0):
    x,w=rule(n,-L,L);weighted=w*switch(x,L)
    q,v=rule(n,0,Q)
    # The compact switch is real and even, so its transform is a cosine integral.
    plus=np.cos((e+q[:,None])*x[None,:])@weighted
    minus=np.cos((e-q[:,None])*x[None,:])@weighted
    weight=q/np.expm1(2*np.pi*q)
    return float(v@(weight*(plus*plus+minus*minus))/(2*np.pi))


def tail_bound(L,Q):
    """Analytic upper bound using ||chi||_1 <= 2 L in dimensionless variables."""
    b=2*math.pi
    return (4*L*L/math.pi)*math.exp(-b*Q)/(1-math.exp(-b*Q))*(Q/b+1/(b*b))


def calculate():
    M=1.0;R=8.0;N=math.sqrt(1-2*M/R)
    kappa=1/(4*M);alpha=kappa/N
    L=8.0;Q=10.0;rows=[]
    for e in (.2,.7,1.5):
        evaluations=[]
        for n in (96,192,384):
            td=time_difference(n,L,e);sp=spectral_difference(n,L,e,Q)
            evaluations.append(dict(nodes_per_integral=n,time_H_minus_B=td,
                spectrum_H_minus_B=sp,absolute_algorithm_difference=abs(td-sp)))
        last=evaluations[-1];before=evaluations[-2]
        assert last['time_H_minus_B']>0 and last['spectrum_H_minus_B']>0
        assert last['absolute_algorithm_difference']<1e-9
        stability=max(abs(last['time_H_minus_B']-before['time_H_minus_B']),
                      abs(last['spectrum_H_minus_B']-before['spectrum_H_minus_B']))
        assert stability<1e-8
        # Independent sign reversal uses the full integrand, not stored output.
        assert abs(time_difference(192,L,-e)-before['time_H_minus_B'])<1e-14
        rows.append(dict(gap_over_alpha=e,gap=e*alpha,evaluations=evaluations,
            H_minus_U=last['time_H_minus_B']/2,refinement_difference=stability))
    zero=float(smooth_kernel(np.array(0.0)))
    assert abs(zero-1/(24*np.pi))<1e-16
    assert np.all(smooth_kernel(np.array([0,.001,.01,1,10]))>0)
    evidence=[STAGE/'1003/research_round_1003_checks.json',
        HERE/'drafts/resource_adoption_checks.json',HERE/'drafts/NEXT.md',
        STAGE/'1003/overall_operation_hypothesis_v2.md',
        STAGE.parent/'archive_301_341/research_note_304.md',
        STAGE.parent/'archive_301_341/research_note_312.md',
        STAGE.parent/'archive_301_341/research_note_317.md',
        STAGE.parent/'archive_301_341/research_note_319.md',
        STAGE.parent/'archive_629_652/research_note_635.md',
        STAGE.parent/'archive_629_652/research_note_648.md',Path(__file__)]
    return dict(round=1004,date='2026-10-07',all_scientific_checks_passed=True,
        kind='one_plus_one_finite_window_response_calibration',
        parameters=dict(M=M,R=R,redshift_lapse=N,kappa=kappa,alpha=alpha,
            detector_proper_acceleration=M/(R*R*N),support_half_width=L/alpha,
            dimensionless_half_width=L,dimensionless_frequency_cut=Q,
            switch='exp(1-1/(1-(alpha*tau/L)^2)) on |alpha*tau|<L'),
        smooth_difference_at_zero=alpha*alpha*zero,
        spectral_H_minus_B_tail_upper_bound=tail_bound(L,Q),
        spectral_H_minus_U_tail_upper_bound=tail_bound(L,Q)/2,
        rows=rows,finite_time_support=True,
        rigorous_quadrature_error_certified=False,
        numerical_calibration_is_four_dimensional=False,
        finite_window_certifies_finite_HH_reservoir=False,
        exact_finite_window_KMS_ratio_claimed=False,
        same_background_is_two_self_consistent_backreaction_solutions=False,
        autonomous_physical_detector_constructed=False,
        hawking_spectrum_derived_from_cognitive_principles=False,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence})


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'/{i}')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-12),path
    else:assert a==b,path


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=calculate()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
