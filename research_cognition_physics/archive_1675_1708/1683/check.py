"""Finite spacelike form-factor cross ratio for the two realized common flavor mechanisms."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(file):
    return json.loads((HERE/file).read_bytes())


def amplitudes(q, r, x, g, w, momenta):
    # Reference zero-momentum weights w_i. The residue in physical momentum units
    # is proportional to w_i*(1+x_i)*g_i^2 for both third-port propagators.
    x,g,w = map(np.asarray,(x,g,w))
    mass_parameter = x*g*g
    residue = w*(1+x)*g*g
    def form(y,t):
        return y*sum(residue/(mass_parameter+g*g*y*y+t))
    return np.array([[form(y,t) for y in (q,r)] for t in momenta])


def squared_cross(q,r,x,g,w,t):
    a,b=amplitudes(q,r,x,g,w,[0,t])
    return float((b[0]*a[1]/(b[1]*a[0]))**2)


def compute():
    cfg=read('inputs.json')
    previous=read(cfg['baseline'])
    pos=read(cfg['unit_neutral_witness'])
    one=read(cfg['free_neutral_witness'])['witnesses']
    r=math.sqrt(7)/5
    qe=1/math.sqrt(3)
    s0=previous['required_large_s_strict_lower']
    xupper=previous['necessary_x_upper_if_unchanged_relative_gains'][1]
    lowerA=(1+(qe**2-r**2)/(2*xupper+qe**2+r**2))**2
    upperB=s0**-2
    single=[]
    for ww in one:
        x=ww['x']
        t=x+r*r
        d=squared_cross(qe,r,[x],[1.],[1.],t)
        np.testing.assert_allclose(d, (1+(qe*qe-r*r)/(2*x+qe*qe+r*r))**2)
        assert d>lowerA>1
        single.append(d)
    q=.3
    x=np.array(pos['masses_squared_over_g2'])
    w=np.array(pos['pole_fractions_at_reference'])
    s=pos['required_charged_gain']
    two=[]
    for g in (np.array([1.,1.]),np.array([1.,.05]),np.array([.05,1.]),np.array([.004,3.])):
        t=max(g*g*(x+r*r))
        high=amplitudes(q,r,x,g,w,[t])[0]
        cross=squared_cross(q,r,x,g,w,t)
        assert high[0]/high[1]<1 and cross < 1/s**2 < upperB
        two.append({'g_over_reference':g.tolist(),'Q2_over_reference_squared':t,'double_ratio':cross})
    cA=one[0]['x']
    monotoneA=[squared_cross(qe,r,[cA],[1],[1],v) for v in (0,.01,.1,1.,10.)]
    monotoneB=[squared_cross(q,r,x,[1,1],w,v) for v in (0,.01,.1,1.,10.)]
    assert min(np.diff(monotoneA))>0 and max(np.diff(monotoneB))<0
    reversed_low=squared_cross(q,r,x,[1.,.05],w,.01)
    assert reversed_low>1
    return {'single_lower_when_Q_at_least_largest_neutral_mass':round(lowerA,12),
            'positive_multipole_upper_at_the_same_scale':round(upperB,12),
            'single_pole_two_witnesses':list(map(lambda z:round(z,12),single)),
            'two_pole_measured_scale_examples':two,'fixed_g_A':monotoneA,'fixed_g_B':monotoneB,
            'B_with_unequal_couplings_low_Q2_point_01':reversed_low,
            'untagged_neutral_weight_max_for_benchmark_sign':(s-1)/(r/q-1)}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--check',action='store_true')
    out=compute()
    file=HERE/'results.json'
    if p.parse_args().check:
        assert out==json.loads(file.read_bytes())
        print('PASS: spectral-scale double ratio, unequal-pole examples and finite-tail allowance')
    else:
        with file.open('x',encoding='utf-8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
