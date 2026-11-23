"""Cross-ratios from the full heavy resolvents and tests with independent pole scales."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def static_kernel(z,x):
    return z*(1+x)/(x+abs(z)**2)


def point(x2,r,q,s):
    def weights(x):
        w=(1-static_kernel(r,x2))/(static_kernel(r,x)-static_kernel(r,x2))
        return np.array([w,1-w])
    def test(x):
        return weights(x)@static_kernel(q,np.array([x,x2]))-s
    lo,hi=0.,.1
    assert test(lo)>0>test(hi)
    for _ in range(60):
        mid=(lo+hi)/2
        if test(mid)>0:
            lo=mid
        else:
            hi=mid
    x=(lo+hi)/2
    xs=np.array([x,x2])
    ws=weights(x)
    assert np.all(ws>0) and np.all(xs>0)
    return xs,ws


def score(xs,ws,g,r,q):
    c=ws*(1+xs)
    port_products=c*g
    def full(z,Q2):
        value=0j
        for x,scale,ends in zip(xs,g,port_products):
            M=scale*np.sqrt(x)
            B=np.array([[M,scale*z],[-scale*z.conjugate(),M]])
            K=np.block([[np.zeros((2,2)),B],[B.conj().T,np.zeros((2,2))]])
            value-=ends*np.linalg.inv(K-1j*np.sqrt(Q2)*np.eye(4))[2,1]
        simple=z*np.sum(c/(xs+abs(z)**2+Q2/g**2))
        np.testing.assert_allclose(value,simple,rtol=1e-13)
        return value
    phi=-np.pi/6
    zn=(3+2*np.exp(4j*phi))/5
    ze=q*np.exp(1j*phi)
    def amplitude_ratio(Q2):
        return abs(full(ze,Q2)/full(zn,Q2))
    def cross(a,b):
        return (amplitude_ratio(b)/amplitude_ratio(a))**2
    a,b=q**2,r**2
    def at(z2):
        return np.sum(c/(xs+z2))
    coeff=c*(1/((xs+b)**2*at(b))-1/((xs+a)**2*at(a)))
    derivative=float(coeff@(1/g**2))
    # Tiny central differences use positive virtuality; regularity at zero
    # also permits an accurate three-forward-point slope.
    step=1e-7
    f=np.log([amplitude_ratio(k*step) for k in (0,1,2)])
    np.testing.assert_allclose((-3*f[0]+4*f[1]-f[2])/(2*step),derivative,atol=2e-7)
    ref=max(g**2*(xs+b))
    return {'x':xs,'reference_fractions':ws,'g_over_g1':g,
            'physical_end_products_in_common_units':port_products,
            'zero_charged_over_neutral':amplitude_ratio(0),
            'coefficients_of_inverse_g_squared_in_zero_log_slope':coeff,
            'zero_log_slope':derivative,'cross_rate_zero_to_0_01':cross(0,.01),
            'nonzero_pair_cross_0_001_to_0_01':cross(.001,.01),
            'neutral_max_unmixed_pole_mass_squared':ref,
            'cross_rate_to_fraction_half_of_neutral_max_squared':cross(0,.5*ref),
            'cross_rate_to_neutral_max_squared':cross(0,ref)}


def tidy(val):
    if isinstance(val,dict):
        return {k:tidy(v) for k,v in val.items()}
    return np.asarray(val).round(11).tolist()


def solve():
    saved=read('../1682/positive_spectrum_results.json')
    q=.3
    r=math.sqrt(7)/5
    s=saved['required_charged_gain']
    xs,ws=point(5.,r,q,s)
    np.testing.assert_allclose(xs,saved['masses_squared_over_g2'],atol=2e-12)
    results={name:score(xs,ws,g,r,q) for name,g in
             [('original_common_g',np.array([1.,1.])),('original_fractions_different_g',np.array([1.,.05]))]}
    x2,w2=point(20.,r,q,s)
    physical=w2*(1+x2)
    equal_ends=physical[0]/physical
    results['both_poles_same_physical_end_product']=score(x2,w2,equal_ends,r,q)
    coeff=results['original_common_g']['coefficients_of_inverse_g_squared_in_zero_log_slope']
    critical=math.sqrt(-coeff[0]/coeff[1])
    assert 11.9<critical<12.1
    cases=list(results.values())
    assert cases[0]['zero_log_slope']<0 and all(k['zero_log_slope']>0 for k in cases[1:])
    for i,case in enumerate(cases):
        assert (case['nonzero_pair_cross_0_001_to_0_01']>1)==(i>0)
        assert case['cross_rate_to_fraction_half_of_neutral_max_squared']<1
        assert case['cross_rate_to_neutral_max_squared']<1/s**2
    np.testing.assert_allclose(cases[2]['physical_end_products_in_common_units'][0],cases[2]['physical_end_products_in_common_units'][1])
    floor=read('../1682/results.json')['required_large_s_strict_lower']
    xcap=max(read('../1682/results.json')['necessary_x_upper_if_unchanged_relative_gains'])
    qa=1/math.sqrt(3)
    bound_A=(1+(qa**2-r**2)/(2*xcap+qa**2+r**2))**2
    return tidy({'cases':results,'g1_over_g2_required_for_original_zero_slope_flip':critical,
                 'sufficient_Q2_over_neutral_max_squared_for_B_decrease':(floor-math.sqrt(floor**2-1))**2,
                 'A_rate_cross_lower_bound_at_own_neutral_max':bound_A,
                 'B_rate_cross_upper_bound_at_own_neutral_max':1/floor**2})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1683 independent PASS: complete propagators, ratio reversal with equal ports, certified finite holdout')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
