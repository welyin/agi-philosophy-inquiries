"""Isolate the supersymmetric gauge part of the two-loop projection above the highest threshold."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def solve():
    main = runpy.run_path(str(HERE/'check.py'))
    adj = np.array([0.,2.,3.])
    beta_mssm = np.array([[199/25,27/5,88/5],[9/5,25,24],[11/5,9,14]])
    index_d = np.array([2/5,0,1])
    index_l = np.array([3/5,1,0])
    delta = (2*np.diag(adj)+4*np.outer(index_d,[1/15,0,4/3])
             +4*np.outer(index_l,[3/20,3/4,0]))
    np.testing.assert_allclose(delta,[[7/15,9/5,32/15],[3/5,7,0],[4/15,0,34/3]])
    b = main['SUSY']+4
    projected = main['W']@(beta_mssm+4*delta)
    table = {}
    levels = main['floor_spectrum'](main['CFG']['step_mass_floor_GeV'],main['MC'])
    floor = 2e5
    assert np.hypot(max(*main['COLOR'],*levels),main['SOFT']) < floor
    for label,n in [('central',0),('favourable_two_width_corner',2)]:
        gs = main['electroweak'](-n,n)
        mt, et = main['CFG']['top_MS_matching_GeV']
        one = main['gauge'](levels,gs=gs,mt=mt+n*et)
        upper, ainv = one['crossing_GeV'],one['inverse_alpha_U']
        def inverse(q):
            return ainv+b*math.log(upper/q)/(2*math.pi)
        def integrate(lo,hi):
            integrals = 2*math.pi/b*np.log(inverse(lo)/inverse(hi))
            nodes,weights = np.polynomial.legendre.leggauss(32)
            dt = math.log(hi/lo)
            values = [1/inverse(lo*math.exp(dt*(u+1)/2)) for u in nodes]
            np.testing.assert_allclose(dt/2*(weights@np.array(values)),integrals,rtol=2e-13)
            return float(projected@integrals/(8*math.pi**2))
        low = integrate(floor,1e7)
        high = integrate(1e7,upper)
        total = integrate(floor,upper)
        np.testing.assert_allclose(total,low+high)
        observed, sigma = main['CFG']['strong_lattice_5f']['ALPHA_2026']
        needed = 1/(observed-n*sigma)-one['a5_QCD']
        table[label] = dict(MU_GeV=upper, alpha_U=1/ainv, target_inverse_shift=needed,
                            through_flavour_cutoff=low, continuation_above_cutoff=high,
                            gauge_piece_sum=total)
    return dict(projected_two_loop_gauge_coefficients=projected.tolist(), rows=table)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(solve()),parse_float=lambda v:float(f'{float(v):.10g}'))
    path = HERE/'gauge_two_loop_results.json'
    if args.check:
        assert result == json.loads(path.read_bytes())
        print("1769 two-loop gauge projection: analytic running integral and independent quadrature agree")
    else:
        with path.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
