"""Reconstruct the charged/neutral semileptonic isospin contrast without supplying V_us."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def data():
    # Seng et al. 2107.14708, Eqs. 3, 4, 11, 12, 16, 17, 29, 33.
    # Widths share the unit 1e-14 MeV. Order is charged then neutral within each block.
    vals = np.array([.0507,5.3167,.4055,1.2866,.15915,.15470,
                     .0021,.0116,493.677,497.611])
    sig = np.array([.0004,.0086,.0011,.0053,.00015,.00015,.0005,.0003,.016,.013])
    corr = np.eye(len(vals))
    for i,j,rho in [(0,1,.01425396),(2,3,-.2676291),(4,5,1.0)]:
        corr[i,j] = corr[j,i] = rho
    cov = corr*np.outer(sig,sig)
    e2 = 4*math.pi/137.035999
    cov[6,7] = cov[7,6] = -4*e2**2*(.0005)**2
    def measure(v):
        bp,gp,bn,gn,ip,ineut,emp,emn,mp,mn = v
        ratio = 2*bp*gp/(bn*gn)*(mn/mp)**5 * ineut/ip
        return ratio*(1+emn)-1-emp
    delta = measure(vals)
    grad = np.array([(measure(vals+s)-measure(vals-s))/(2*s[i])
                     for i,s in enumerate(np.diag(sig*1e-3))])
    width = math.sqrt(grad@cov@grad)
    return delta,width


def coeff(s):
    mesons = (494.2/134.8)**2
    return 6*(mesons+.219*(1+s)/2)/(s**2-1)


def predict(y,s):
    x = 2*y/s-1
    if not (0<x<1):
        raise ValueError("requires ordered positive u,d,s masses")
    q2 = y*(s**2-1)/(4*(s-y))
    return dict(y=y, mu_over_md=x, Q=math.sqrt(q2), rate_correction=coeff(s)*(s/y-1))


def mass_from_charge(delta,s):
    return s/(1+delta/coeff(s))


def solve():
    old = json.loads((HERE.parent/'1734/inputs.json').read_bytes())
    s = old['BMW_S'][0]
    es = math.hypot(*old['BMW_S'][1:])
    delta,err = data()
    refits = json.loads((HERE.parent/'1766/independent_results.json').read_bytes())
    rows = {}
    for tag in ['refits','exploratory_64_deg']:
        rows[tag] = {}
        for key,row in refits[tag].items():
            val = predict(1/row['ds'],s)
            val['charge_difference_in_experimental_widths'] = (val['rate_correction']-delta)/err
            val['eta_Q_difference_in_published_widths'] = (val['Q']-old['eta_Q'][0])/old['eta_Q'][1]
            assert abs(mass_from_charge(val['rate_correction'],s)-val['y'])<1e-12
            rows[tag][key] = val
    yrange = [mass_from_charge(delta+k*err,s+j*es) for k in (-1,1) for j in (-1,1)]
    return dict(electron_isospin_rate_contrast=delta, local_width=err, S=s,
                central_y_from_kaons=mass_from_charge(delta,s),
                inverted_one_width_corners=[min(yrange),max(yrange)], branches=rows)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(solve()),parse_float=lambda v:float(f'{float(v):.9g}'))
    path = HERE/'kaon_charge_results.json'
    if args.check:
        assert result == json.loads(path.read_bytes())
        print("1768 charge contrast: correlated raw Ke3 ratio, mass inverse, 12 B-channel points")
    else:
        with path.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
