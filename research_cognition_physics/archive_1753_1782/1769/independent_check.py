"""Invert the common gauge threshold and the subtracted photon-polarization slope."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def bisect_increasing(fun, lo, hi):
    if not (fun(lo) <= 0 <= fun(hi)):
        raise ValueError("root is not bracketed")
    for _ in range(90):
        mid = (lo+hi)/2
        if fun(mid) <= 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def threshold(mass, soft):
    return 2*math.log(mass)/3 + math.log(mass*mass+soft*soft)/6


def em_moment(mass, soft):
    return 1/mass**2 + .25/(mass**2+soft**2)


def delta_inverse_alpha(mass, soft, low, high, order=48):
    x, w = np.polynomial.legendre.leggauss(order)
    x, w = (x+1)/2, w/2
    z = x*(1-x)
    def difference(m2):
        return np.log1p((high-low)*z/(m2+low*z))
    fermion = 2*np.dot(w, z*difference(mass**2))/math.pi
    two_scalars = np.dot(w, (1-2*x)**2*difference(mass**2+soft**2))/(2*math.pi)
    return float(fermion+two_scalars)


def solve():
    old = json.loads((HERE.parent/'1698/results.json').read_bytes())['shared_witness']
    mz = json.loads((HERE.parent/'1369/inputs.json').read_bytes())['mu_ref_GeV']
    soft = old['soft_GeV']
    level = (sum(threshold(m,soft) for m in old['colored_fermions_GeV'])
             - old['required_effective_split'])
    equal = bisect_increasing(lambda m: 2*threshold(m,soft)-level, mz, 300)
    def partner(m):
        return bisect_increasing(lambda b: threshold(m,soft)+threshold(b,soft)-level,
                                 mz*(1-1e-9), 2000)
    np.testing.assert_allclose(equal, old['weak_fermions_GeV'][0], rtol=1e-8)
    def moment_sum(m):
        return em_moment(m,soft)+em_moment(partner(m),soft)
    rows = {}
    for tag, m in [('lower_moment',equal), ('mass_reconstruction',107.), ('upper_moment',mz)]:
        masses = [m, partner(m)]
        sig = sum(delta_inverse_alpha(p, soft, 12.25, 3434.) for p in masses)
        refined = sum(delta_inverse_alpha(p, soft, 12.25, 3434., 80) for p in masses)
        moment = moment_sum(m)
        linear = moment*(3434.-12.25)/(15*math.pi)
        assert linear/(1+3434./(4*min(masses)**2)) < sig < linear
        np.testing.assert_allclose(refined,sig,rtol=1e-12)
        rows[tag] = dict(masses_GeV=masses, moment_per_GeV2=moment,
                         delta_inverse_alpha_L3=sig, fraction_of_L3_width=sig/1.29)
        if tag=='mass_reconstruction':
            recovered = bisect_increasing(lambda a:moment-moment_sum(a), mz, equal)
            np.testing.assert_allclose([recovered,partner(recovered)],masses,rtol=1e-12)
    target = 2*em_moment(225, soft)
    return dict(witnesses=rows, upper_moment_if_both_at_least_225=target,
                minimal_moment_over_225_bound=moment_sum(equal)/target)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(solve()),parse_float=lambda v:float(f'{float(v):.10g}'))
    path = HERE/'independent_results.json'
    if args.check:
        assert result == json.loads(path.read_bytes())
        print("1769 independent: two masses from joint gauge and QED response, L3 kernel verified")
    else:
        with path.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
