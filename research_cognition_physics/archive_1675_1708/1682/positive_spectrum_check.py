"""Normalized multipole transfer with common positive residues for neutral and charged flavors."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def kernel(q,x):
    return q*(1+x)/(q*q+x)


def compute():
    r = math.sqrt(7)/5
    qc = 1/math.sqrt(3)
    witness = json.loads((HERE.parent/'1678/independent_results.json').read_bytes())['joint_witness']
    s = witness['common_abs_t']
    all_single_max = (1+r)/(2*math.sqrt(r))
    qcut = 2*r / (s*(1+r) + math.sqrt((s*(1+r))**2-4*r))
    def remap(q,u):
        return q*(1-r*r)*u / (r*(1-q*q)+(q*q-r*r)*u)
    for q in [.1,.3,r,qc,.8,1]:
        for x in [.0001,.03,.1,r,1.,30.]:
            assert abs(remap(q,kernel(r,x))-kernel(q,x))<2e-14
    # All channels use the same two poles, portal signs, weights and reference.
    q,x2 = .3,5.
    def weight(x1):
        return (1-kernel(r,x2))/(kernel(r,x1)-kernel(r,x2))
    def gain(x1):
        w = weight(x1)
        return w*kernel(q,x1)+(1-w)*kernel(q,x2)
    a,b = .02,.05
    assert gain(a)>s>gain(b)
    for _ in range(70):
        mid = (a+b)/2
        if gain(mid)>s:
            a = mid
        else:
            b = mid
    x1 = (a+b)/2
    w = weight(x1)
    xs = np.array([x1,x2])
    ws = np.array([w,1-w])
    residues = ws*(1+xs)
    probs = np.array([1+math.sqrt(3)*q,1.,1-math.sqrt(3)*q])/3
    phi = -np.pi/6
    z = np.exp(4j*phi)
    ze = probs @ np.array([1,z,z*z])
    zn = (3+2*z)/5
    assert min(probs)>0 and abs(ze-q*np.exp(1j*phi))<1e-15
    def response(value,momentum2_over_g2):
        p = momentum2_over_g2
        denominator = sum(residues/(p+xs+1))
        return value*sum(residues/(p+xs+abs(value)**2))/denominator
    np.testing.assert_allclose(response(ze,0),s*np.exp(1j*phi),rtol=1e-11)
    np.testing.assert_allclose(response(zn,0),zn/abs(zn),rtol=1e-13)
    def record(x):
        return round(float(x),12)
    return {'r_neutral':record(r),'original_charged_modulus':record(qc),
            'sharp_ceiling_at_original_charged_modulus':record(kernel(qc,r)),
            'all_q_at_or_above_neutral_uniform_ceiling':record(all_single_max),
            'necessary_q_upper_for_alternative':record(qcut),'masses_squared_over_g2':list(map(record,xs)),
            'pole_fractions_at_reference':list(map(record,ws)),
            'three_charged_path_probabilities':list(map(record,probs)),
            'required_charged_gain':record(s),
            'finite_spacelike_amplitudes':[{'Q_squared_over_g_squared':p,'charged_modulus':record(abs(response(ze,p))),
                 'neutral_modulus':record(abs(response(zn,p)))} for p in [0.,.1,1.,10.]]}


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    result=compute()
    target=HERE/'positive_spectrum_results.json'
    if parser.parse_args().check:
        assert result==json.loads(target.read_bytes())
        print('PASS: sharp positive-spectrum envelopes and finite two-pole common-messenger witness')
    else:
        with target.open('x',encoding='utf-8') as f:
            f.write(json.dumps(result,indent=2)+'\n')
