"""Independent checks of round 1223's old-material probability-stability bound.

No author code is imported. Decimal calculations check constants; rational
inequalities below certify the deliberately rounded analytic bound. Gaussian
tests check probability transport, not a substitute for the continuum proof.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
PI='3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647'


def constants(precision):
    with localcontext() as context:
        context.prec=precision
        pi=D(PI)
        volume=4*pi/3
        tau, lam, duration, nu=D('.25'),D('.05'),D('.55'),D('.01')
        kappa=lam*duration**2
        hnorm=(32*pi/3).sqrt()
        q=tau*lam*D(2).sqrt()*duration*hnorm/(1-kappa)**2
        c0=D('.75')*((volume/(pi*pi)).ln()/3).exp()
        c0_opt=D('.75')*((volume/(2*pi*pi)).ln()/3).exp()
        rms=(2*c0+volume/5).sqrt()
        coef=rms*q/((2*pi).sqrt()*nu)
        energy=925*pi/3
        particle=1882*pi/12
        assert energy<1000 and particle<500
        assert coef<D('3.32') and c0_opt<c0
        data=dict(kappa=kappa,h_L1L2_upper=hnorm,spatial_volume=volume,
            delta_q_L1L2_coefficient=q,vacuum_covariance_conservative=c0,
            vacuum_covariance_optimal_split=c0_opt,
            mean_and_covariance_rms_coefficient=rms,
            full_readout_TV_coefficient=coef,
            coherent_preparation_energy_upper=energy,
            coherent_preparation_number_upper=particle)
        return {k:format(v,'.24g') for k,v in data.items()}


def rational_certificate():
    # 3.14 < pi < 22/7. Every inequality uses exact rational arithmetic.
    pi_lo,pi_hi=F(314,100),F(22,7)
    assert F(2)<F(1415,1000)**2
    assert 32*pi_hi/3<F(58,10)**2
    assert F(1,20)*F(55,100)**2<F(16,1000)
    volume_hi=F(42,10)
    assert 4*pi_hi/3<volume_hi
    # C0=.75*(V/pi^2)^(1/3), using V/pi^2 = 4/(3*pi).
    c0_hi=F(57,100)
    assert (F(3,4)**3)*4/(3*pi_lo)<c0_hi**3
    q_hi_raw=F(1,4)*F(1,20)*F(1415,1000)*F(55,100)*F(58,10)/F(984,1000)**2
    q_hi=F(585,10000)
    assert q_hi_raw<q_hi
    rms_hi=F(1408,1000)
    assert 2*c0_hi+volume_hi/5<rms_hi**2
    assert F(25,10)**2<2*pi_lo
    tv_hi=rms_hi*q_hi/F(25,1000)
    assert tv_hi<F(332,100)
    return dict(status='PASS',TV_coefficient_strictly_below='3.32',
        intermediate_rational_TV_upper=str(tv_hi),
        no_floating_point_rounding_used=True)


def normal_cdf(x,mean,sd):
    if x==math.inf:return 1.
    if x==-math.inf:return 0.
    return .5*math.erfc(-(x-mean)/(math.sqrt(2)*sd))


def gaussian_tv(m1,s1,m2,s2):
    if s1==s2:
        return math.erf(abs(m1-m2)/(2*math.sqrt(2)*s1))
    a=.5*(1/s2**2-1/s1**2)
    b=m1/s1**2-m2/s2**2
    c=math.log(s2/s1)+.5*(m2*m2/s2**2-m1*m1/s1**2)
    roots=sorted(float(x.real) for x in np.roots([a,b,c]) if abs(x.imag)<1e-10)
    cuts=[-math.inf,*roots,math.inf]
    total=0.
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        p=normal_cdf(hi,m1,s1)-normal_cdf(lo,m1,s1)
        q=normal_cdf(hi,m2,s2)-normal_cdf(lo,m2,s2)
        total+=abs(p-q)
    return total/2


def probability_checks():
    tests=[]
    cases=[(np.array([1.,.2]),np.array([.98,.23]),.02,.07,.1),
           (np.array([.01,.02]),np.array([.012,.019]),-.03,.01,.2),
           (np.array([1.,0.]),np.array([0.,1.]),0.,0.,.1),
           (np.zeros(2),np.zeros(2),0.,.001,.01)]
    for a,b,m1,m2,nu in cases:
        variance_delta=float(np.dot(a-b,a-b))
        gram=np.array([[np.dot(a,a),np.dot(a,b)],[np.dot(a,b),np.dot(b,b)]])
        assert np.linalg.eigvalsh(gram).min()>-1e-14
        s1=math.sqrt(float(np.dot(a,a))+nu*nu)
        s2=math.sqrt(float(np.dot(b,b))+nu*nu)
        tv=gaussian_tv(m1,s1,m2,s2)
        bound=math.sqrt((m1-m2)**2+variance_delta)/(math.sqrt(2*math.pi)*nu)
        assert tv<=min(1,bound)+1e-12
        tests.append(dict(actual_Gaussian_TV=tv,rms_envelope=bound))
    # Equal marginals do not imply equal joint source or post-measurement states.
    assert tests[2]['actual_Gaussian_TV']==0 and tests[2]['rms_envelope']>0
    return tests


def calculate():
    low,high=constants(70),constants(95)
    assert low==high
    sources=[HERE/'independent_check.py',BASE/'archive_467_530/research_note_524.md']
    return dict(schema='round1223_independent_stability_v1',status='PASS',check_groups=3,
        scope='Old three-field material, Gaussian input and fixed noisy marginal readout; no autonomous source or joint-output bound',
        constants=low,rational_envelope=rational_certificate(),
        exact_Gaussian_probability_tests=probability_checks(),
        normal_finite_energy_preparation=True,
        full_probability_given_classical_history=True,
        quantum_source_and_poststate_error_certified=False,
        continuum_PDE_grid_not_used_as_certificate=True,
        source_sha256={p.relative_to(BASE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    payload=calculate()
    out=HERE/'independent_results.json'
    if args.write_results:
        with out.open('x',encoding='utf-8') as stream:
            json.dump(payload,stream,ensure_ascii=False,indent=2);stream.write('\n')
    else:
        assert json.loads(out.read_text(encoding='utf-8'))==payload
    print(json.dumps(dict(status='PASS',check_groups=3,
        full_probability_TV_coefficient=payload['constants']['full_readout_TV_coefficient'],
        certified_rounded_coefficient='3.32'),ensure_ascii=False))


if __name__=='__main__':main()
