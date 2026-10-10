"""Finite calibration of the explicit clock/rod/echo contract, not a spacetime proof."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import math
import numpy as np

BASE = Path(__file__).resolve().parent


def matmul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def boost(b, g, lam):
    return [[g/lam, -g*b/lam, F(0), F(0)],
            [-g*b/lam, g/lam, F(0), F(0)],
            [F(0), F(0), 1/lam, F(0)],
            [F(0), F(0), F(0), 1/lam]]


def compute():
    exact = []
    for b, g in [(F(0), F(1)), (F(3,5), F(5,4)), (F(5,13), F(13,12)), (F(8,17), F(17,15))]:
        assert g*g*(1-b*b) == 1
        for lam in ([F(1)] if b == 0 else [F(1), F(6,5), F(4,5)]):
            ap, at, rate = lam/g, lam, 1/(lam*g)
            echo_parallel = rate*2*ap/(1-b*b)
            echo_transverse = rate*2*at*g
            # Independent intersection calculation for successive emitted pulses.
            b_arrivals = [(1+b)*F(k)/rate for k in range(1,5)]
            a_arrivals = [rate*F(k)/(1-b) for k in range(1,5)]
            dab = b_arrivals[1]-b_arrivals[0]
            dba = a_arrivals[1]-a_arrivals[0]
            assert echo_parallel == echo_transverse == 2
            assert dab/dba == lam*lam
            assert dab == lam*g*(1+b) and dba == g*(1+b)/lam
            m = boost(b,g,lam)
            inv_rule = boost(-b,g,lam)
            comp = matmul(inv_rule,m)
            for i in range(4):
                for j in range(4):
                    assert comp[i][j] == (1/(lam*lam) if i == j else 0)
                    quadratic = sum((1 if k == 0 else -1)*m[k][i]*m[k][j] for k in range(4))
                    assert quadratic == ((1 if i == 0 else -1)/(lam*lam) if i == j else 0)
            exact.append({'beta':str(b),'lambda':str(lam),'a_parallel':str(ap),'a_perp':str(at),
                          'clock_rate':str(rate),'echo_parallel':str(echo_parallel),
                          'echo_transverse':str(echo_transverse),'doppler_A_from_B':str(dab),
                          'doppler_B_from_A':str(dba),'doppler_ratio':str(dab/dba)})
    # Three-dimensional ray/moving-mirror intersections: the full affine lattice
    # and finite window are inputs. No random sampling is a proof of that input.
    rng = np.random.default_rng(1086)
    residuals = []
    min_time = float('inf')
    for beta in [0., .2, .6, .8]:
        gamma = 1/math.sqrt(1-beta*beta)
        for lam in ([1.] if beta == 0 else [.8, 1., 1.2]):
            for _ in range(24):
                n = rng.normal(size=3)
                n /= np.linalg.norm(n)
                d = lam*np.array([n[0]/gamma,n[1],n[2]])
                disc = beta*beta*d[0]*d[0]+(1-beta*beta)*float(d@d)
                out = (beta*d[0]+math.sqrt(disc))/(1-beta*beta)
                back = (-beta*d[0]+math.sqrt(disc))/(1-beta*beta)
                residuals.extend([abs(np.linalg.norm(d+np.array([beta*out,0,0]))-out),
                                  abs(np.linalg.norm(np.array([beta*back,0,0])-d)-back),
                                  abs((out+back)/(lam*gamma)-2)])
                min_time = min(min_time,out,back)
    max_residual = max(residuals)
    assert min_time > 0 and max_residual < 1e-12
    # Analytic log bound and its equality cases; all are finite-window examples.
    for eta, xi in [(0.,0.),(.02,.01),(.1,.2)]:
        lp, lm = (eta+xi)/2, (eta-xi)/2
        assert abs(lp) <= (eta+xi)/2+1e-15
        assert abs(lp+lm) <= eta+1e-15 and abs(lp-lm) <= xi+1e-15
        assert math.isclose(.5*math.log(math.exp(lp)**2),lp,abs_tol=1e-15)
    witness = next(r for r in exact if r['beta']=='3/5' and r['lambda']=='6/5')
    assert F(witness['doppler_A_from_B'])-F(witness['doppler_B_from_A']) == F(11,15)
    return {'round':1086,'passed':True,'scope':'declared finite ray-clock-rod kinematics only',
            'exact_calibration_cases':len(exact),'exact_cases':exact,
            'three_dimensional_ray_cases':240,'max_ray_residual':max_residual,
            'minimum_leg_time':min_time,'finite_doppler_gap':'11/15',
            'nonunit_scale_selectable_by_double_echo':False,
            'lorentz_from_CO1_CO5_proved':False,'complete_quantum_material_model':False,
            'numpy_version':np.__version__}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    data = compute()
    if args.write:
        (BASE/'results.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old = json.loads((BASE/'results.json').read_text(encoding='utf-8'))
        assert old['exact_cases'] == data['exact_cases']
        assert old['finite_doppler_gap'] == data['finite_doppler_gap']
        assert old['three_dimensional_ray_cases'] == data['three_dimensional_ray_cases']
        assert abs(old['max_ray_residual']-data['max_ray_residual']) < 1e-12
    print(json.dumps({k:v for k,v in data.items() if k!='exact_cases'},ensure_ascii=False))
