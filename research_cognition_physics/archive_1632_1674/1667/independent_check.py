"""Independent one-loop scalar and Higgs coupling test with the CMS Run 2 fit."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def background(mh):
    v = read('../1373/inputs.json')['v_MZ_GeV']
    light = read('../1367/inputs.json')
    yd, ys, ye, ymu = np.multiply(light['points']['MZ']['central'], light['units'])
    up = read('../1368/inputs.json')['MZ']
    heavy = read('../1369/inputs.json')
    fermions = [(up['up'][0], 3, 2/3), (up['charm'][0], 3, 2/3),
                (heavy['yt'][0], 3, 2/3), (yd, 3, -1/3), (ys, 3, -1/3),
                (heavy['yb'][0], 3, -1/3), (ye, 1, -1), (ymu, 1, -1),
                (heavy['ytau'][0], 1, -1)]
    w = v*heavy['g2'][0]/2
    x = (mh/(2*w))**2
    photon = -(2*x*x+3*x+3*(2*x-1)*math.asin(math.sqrt(x))**2)/x**2
    gluon = 0j
    for yukawa, colours, charge in fermions:
        tau = (math.sqrt(2)*v*yukawa/mh)**2
        continuation = cmath.asin(1/math.sqrt(tau))**2
        fermion = 2*tau*(1+(1-tau)*continuation)
        photon += colours*charge**2*fermion
        if colours == 3:
            gluon += fermion
    return v, gluon, photon


def scalar(mh, mass):
    # Exact convergent series for 8*integral_triangle xy/(1-4zxy).
    z = (mh/(2*mass))**2
    assert 0 <= z < 0.5
    term = 1/3
    total = term
    for j in range(80):
        term *= 4*z*(j+2)**2/((2*j+5)*(2*j+6))
        total += term
    return total


def preimage(amp, coefficient, lower, upper):
    def horizontal(radius):
        return math.sqrt((radius * abs(amp))**2-amp.imag**2)
    near, far = horizontal(lower), horizontal(upper)
    return [[(-amp.real-far)/coefficient, (-amp.real-near)/coefficient],
            [(-amp.real+near)/coefficient, (-amp.real+far)/coefficient]]


def intersections(first, second):
    return [[max(a,c), min(b,d)] for a,b in first for c,d in second
            if max(a,c) <= min(b,d)]


def compute():
    # CMS Fig. 4 left, https://arxiv.org/pdf/2207.00043#page=11 .
    # Use the collaboration's nominal mass convention, independent of the ATLAS point.
    mh = read('../1373/inputs.json')['higgs_mass_GeV']
    v,g,p = background(mh)
    gauge = read('../1663/results.json')
    cap = gauge['four_copies_geometric_mean_upper_GeV']
    point = gauge['central_four_copy_geometric_mean_upper_GeV']
    ratio = (3*((2/3)**2+(-1/3)**2))/2
    assert abs(ratio-5/6) < 2e-15
    g_window, p_window = [0.92-2*.08,0.92+2*.08], [1.10-2*.08,1.10+2*.08]
    gs = preimage(g,1,*g_window)
    ps = preimage(p,ratio,*p_window)
    joint = intersections(gs,ps)
    assert len(joint)==1 and joint[0][0]<0<joint[0][1]
    upper = joint[0][1]
    signed_masses = [240]*4
    signed_f = np.array([.25,-.25,0.,0.])
    lam = 2*(np.array(signed_masses)/v)**2*signed_f
    higgs_quartic = mh**2/(2*v*v)
    sufficient_quartic = max(abs(lam))**2/(4*higgs_quartic)
    assert sufficient_quartic<1
    # Local soft two-Higgs response via an independent finite difference in h/v.
    def sum_log(t):
        return float(np.log1p(signed_f*(2*t+t*t)).sum()/2)
    step = 0.0001
    numeric = (sum_log(step)+sum_log(-step))/step**2
    analytic = float(sum(signed_f-2*signed_f**2))
    assert abs(analytic-numeric)<1e-7
    endpoints = {}
    for label, geo in [('all_inputs',cap),('charged_witness',point)]:
        min_shift = 8*scalar(mh,geo)
        gamma_upper = abs(p+ratio*min_shift)/abs(p)
        assert 8*scalar(mh,91.1876) < -p.real/ratio
        endpoints[label] = {'all_ew_kappa_g_min':abs(g+min_shift)/abs(g),
                            'all_ew_kappa_gamma_max':gamma_upper,
                            'common_positive_lambda_upper':upper*geo**2/(4*v*v*scalar(mh,geo)),
                            'common_positive_fraction_upper':upper/min_shift}
    return {'cms_inputs':{'mh_GeV':mh,'v_GeV':v,'kappa_g_window':g_window,
                          'kappa_gamma_window':p_window},
            'sm_one_loop':{'g':[g.real,g.imag],'gamma':[p.real,p.imag]},
            'gluon_shift_sets':{'from_g':gs,'from_photon':ps,'common':joint},
            'nonnegative_fractions':{'sum_upper':1.5*upper,
                                     'kappa_gamma_min_from_g_only':abs(p+ratio*upper)/abs(p)},
            'fixed_gain_comparisons':endpoints,
            'signed_competitor':{'portal':lam.tolist(),'min_sufficient_scalar_quartic':sufficient_quartic,
                                 'soft_second_derivative':analytic}}


def compact(obj):
    if isinstance(obj,dict):
        return {k:compact(v) for k,v in obj.items()}
    if isinstance(obj,list):
        return [compact(x) for x in obj]
    return float(f'{obj:.10g}')


if __name__ == '__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    result=compact(compute())
    if args.check:
        assert result == read('independent_results.json')
        print('1667 independent PASS: CMS finite-loop intervals, common portal bounds and signed potential')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
