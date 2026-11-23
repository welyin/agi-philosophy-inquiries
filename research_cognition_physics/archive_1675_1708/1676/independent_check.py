"""Positive triangle-integral certificate for the finite scalar H->Z gamma kernel."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def series(h,t,terms=78):
    # Triangle beta integration. Coefficients of A0 are positive;
    # their successive ratio after factoring out h is < 1.
    weight=1/3
    a0=b0=0.
    for n in range(terms):
        multiplier=sum((1-k/(n+1))*t**k for k in range(n+1))
        a0+=weight
        b0+=weight*multiplier
        weight*=h*4*(n+2)**2/((2*n+6)*(2*n+5))
    error=weight/((1-h)*(1-t))
    assert error<1e-16
    return a0,b0


def quadrature(h,t):
    # Map the simplex to a square; independent of both series and closed I1.
    x,w=np.polynomial.legendre.leggauss(36)
    x=(x+1)/2
    w=w/2
    y=(1-x[:,None])*x
    v=1-x[:,None]-y
    xy=x[:,None]*y
    measure=8*w[:,None]*w*(1-x[:,None])*xy
    return [np.sum(measure/(1-4*h*xy)),
            np.sum(measure/(1-4*h*xy-4*h*t*x[:,None]*v))]


def solve():
    cfg=read('../1667/inputs.json')
    low=cfg['mass_floor_GeV']
    mh=cfg['higgs_GeV']
    prior=runpy.run_path(str(HERE.parent/'1504/check.py'))
    old=read('../1373/inputs.json')
    mw,mz,fs=prior['prior'].common_inputs(old)
    cw=mw/mz
    sw=np.sqrt(1-cw*cw)
    sc=sw*cw
    sw2=sw*sw
    def pair(mass,zpole):
        h=(mh/(2*mass))**2
        t=(zpole/mh)**2
        p=series(h,t)
        np.testing.assert_allclose(p,quadrature(h,t),rtol=2e-14)
        i1=prior['integrals'](mh,zpole,mass)[0]
        a0=prior['prior'].loop_functions(mh,mass)[0]
        np.testing.assert_allclose(p,[a0,2*i1],rtol=3e-10)
        return np.array(p)
    def standard(zpole):
        w2=sw2/(cw*cw)
        tau=(2*mw/mh)**2
        i1,i2=prior['integrals'](mh,zpole,mw)
        z=cw/sw*(4*(3-w2)*i2+((1+2/tau)*w2-(5+2/tau))*i1)
        for mf,nc,charge in fs:
            a,b=prior['integrals'](mh,zpole,mf)
            t3=.5*np.sign(charge)
            z+=2*nc*charge*(t3-2*charge*sw2)/sc*(a-b)
        return z
    # Finite mass comparisons use a constant electroweak scheme. Separately
    # quote sensitivity of the dimensionless kernel to the external Z mass.
    ratios={}
    for name,zpole in [('tree_level_1504',mz),('Z_kinematics_91_1876',low)]:
        pts=[pair(m,zpole) for m in (low,120.,240.,480.)]
        vals=[b/a for a,b in pts]
        assert all(x>y for x,y in zip(vals,vals[1:])) and vals[-1]>1
        ratios[name]=vals
    rmax=ratios['tree_level_1504'][0]
    smz=standard(mz)
    oldres=read('../1667/results.json')
    xmax=max(x[1] for x in oldres['gluon_amplitude_shift_sets']['production'])
    max_delta=rmax*((.5-sw2)*4/3+1/3)*xmax/sc
    minrate=abs(smz-max_delta)**2/abs(smz)**2
    masses=[120,480,240,240]
    forms=np.array([pair(m,mz) for m in masses])
    a0,b0=forms.T
    fA=np.array([.1,0.,0.,0.])
    fB=np.array([0.,.1*a0[0]/a0[1],0.,0.])
    def get(f):
        g=2*f@a0
        d=(5/3)*f@a0
        k=-(3/sc)*(.5-5*sw2/9)*(f@b0)
        return np.array([g,d,k])
    A,B,N=[get(f) for f in (fA,fB,fA-fB)]
    np.testing.assert_allclose(A[:2],B[:2],atol=1e-14)
    assert abs(N[2])>1e-3 and N[0]**2+N[1]**2<1e-26
    assert A[0]<xmax and abs(np.prod(masses)**.25-240)<1e-11
    def rounded(x):
        return np.asarray(x).round(11).tolist()
    return {'reference_Z_GeV':round(float(mz),10),'sin2thetaW':round(float(sw2),10),
            'r_2I1_over_A0_at_masses':[low,120,240,480],
            'ratios':{s:rounded(v) for s,v in ratios.items()},
            'baseline_complex_Z_amplitude':rounded([smz.real,smz.imag]),
            'nonnegative_portal_Z_amplitude_abs_ceiling':round(float(max_delta),11),
            'positive_partial_width_outer_range':[round(float(minrate),11),1],
            'fixed_mass_positive_A':rounded(A),'fixed_mass_positive_B':rounded(B),
            'zero_calibrations_signed_shift':rounded(N),
            'masses_GeV':masses,'positive_portal_fractions_A':rounded(fA),
            'positive_portal_fractions_B':rounded(fB)}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1676 independent PASS: positive integral, sharp finite-mass factor, split-port and signed witnesses')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as output:
            output.write(json.dumps(result,indent=2)+'\n')
