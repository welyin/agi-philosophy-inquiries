"""Read-only replay by default. Exact bounds and numerical diagnostics are separate."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

D = Path(__file__).resolve().parent
S = D.parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check_sources():
    b = json.loads((D/'source_baseline.json').read_text(encoding='utf-8-sig'))
    for name, expected in b['source_sha256'].items():
        assert sha(S/name) == expected, name
    for name, expected in b['external_source_sha256'].items():
        assert sha(Path(name)) == expected, name
    return len(b['source_sha256']) + len(b['external_source_sha256'])

def wightman(q, r, rp, dy, deta):
    alpha = math.acosh((r*r+rp*rp+dy*dy)/(2*r*rp))
    assert abs(deta) < alpha
    beta = 2*math.pi/q
    return math.sinh(q*alpha)/(4*math.pi*beta*r*rp*math.sinh(alpha)
            *(math.cosh(q*alpha)-math.cosh(q*deta)))

def quadrature(n, logarithmic):
    z, v = np.polynomial.legendre.leggauss(n)
    if logarithmic:
        s = (z+1)*math.log(2)/2
        r = np.exp(s)
        weights = v*math.log(2)*r/2
    else:
        r = (z+3)/2
        weights = v/2
    bump = np.exp(-1/((r-1)*(2-r)))
    return float(np.sum(weights*bump/r**4)/np.sum(weights*bump))

def compute():
    bindings = check_sources()
    lo, hi = F(99749,100000), F(1002491,1000000)
    assert lo**4 < F(99,100) and hi**4 > F(101,100)
    reject, retained = F(11,10), F(10001,10000)
    d_reject, d_retained = reject**4-1, retained**4-1
    reject_ratio_lower = 100*abs(d_reject)
    retained_ratio_upper = 1600*abs(d_retained)
    assert reject_ratio_lower == F(4641,100) > 1
    assert 0 < retained_ratio_upper < 1
    # Expand both the accepted source band and rejected prediction by epsilon/2.
    separation = reject_ratio_lower-F(1,2)-(1+F(1,2))
    assert separation == F(4441,100) > 0
    js = [quadrature(n, log) for n in (96,192) for log in (False, True)]
    assert all(1/16 < j < 1 for j in js)
    spread = max(js)-min(js)
    assert spread < 5e-12
    j = js[-1]
    eps = 1/(768000*math.pi**2)
    delta = 480*math.pi**2*eps/j
    assert 0 < delta < .01
    # Eq.(2) tensor evaluated in an inertial orthonormal frame; check invariant.
    eta = np.diag([-1.,1.,1.,1.])
    invariant_residuals = []
    for rapidity in (0., .2, 1.3):
        u = np.array([math.cosh(rapidity),math.sinh(rapidity),0.,0.])
        tensor = 4*np.outer(u,u)+eta  # coefficient c=1, energy density=3
        energy = float((eta@u)@tensor@(eta@u))
        invariant = float(np.trace((tensor@eta)@(tensor@eta)))
        invariant_residuals.append(max(abs(energy-3),abs(invariant-12)))
    assert max(invariant_residuals) < 1e-11
    wm_residuals, altered = [], []
    for r,rp,dy,deta in [(1.1,1.7,.4,.2),(1.4,1.8,.6,.1),(.8,1.3,.2,.15)]:
        wm = 1/(4*math.pi**2*(r*r+rp*rp+dy*dy-2*r*rp*math.cosh(deta)))
        wm_residuals.append(abs(wightman(1,r,rp,dy,deta)-wm))
        altered.append(wightman(1.1,r,rp,dy,deta)-wm)
    assert max(wm_residuals) < 1e-12
    assert all(x>0 for x in altered)
    return {
        'status':'PASS_WITH_SCOPE', 'source_bindings':bindings,
        'exact': {
            'q_outward_lower':str(lo), 'q_outward_upper':str(hi),
            'lower_fourth_minus_99_100':str(lo**4-F(99,100)),
            'upper_fourth_minus_101_100':str(hi**4-F(101,100)),
            'rejected_q':str(reject), 'rejection_ratio_lower':str(reject_ratio_lower),
            'retained_nonunit_q':str(retained), 'retention_ratio_upper':str(retained_ratio_upper),
            'conditional_readout_gap_in_epsilon_units':str(separation)},
        'diagnostics': {
            'smooth_weight_inverse_fourth_moment':j, 'quadrature_spread':spread,
            'epsilon_length_one':eps, 'actual_delta_for_example_weight':delta,
            'q_interval_for_example_weight':[(1-delta)**.25,(1+delta)**.25],
            'lorentz_invariant_residual_max':max(invariant_residuals),
            'wightman_common_reference_residual_max':max(wm_residuals),
            'changed_reference_wightman_differences':altered},
        'scope': {'new_scientific_rounds':0, 'actual_readout_certified':False,
                  'finite_backreaction_certified':False, 'geometric_entropy_identified':False,
                  'temperature_exactly_selected_at_positive_error':False,
                  'new_axiom_adopted':False, 'goal_completed':False}}

def close(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-12)
    return a==b

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=compute()
    path=D/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        mode='created_once'
    else:
        assert close(result,json.loads(path.read_text(encoding='utf-8')))
        mode='read_only_match'
    print(json.dumps({'mode':mode,'status':result['status'],
                     'exact':result['exact'],'diagnostics':result['diagnostics']},ensure_ascii=False))

if __name__=='__main__':
    main()
