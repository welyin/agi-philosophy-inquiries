"""Read-only by default. Finite arithmetic for a corollary, not a parent simulation."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

D = Path(__file__).resolve().parent
ROOT = D.parents[3]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def atan_bounds(x, n=40):
    total = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)), F(0))
    other = total + (-1)**n*x**(2*n+1)/F(2*n+1)
    return min(total, other), max(total, other)

def sin_bounds(x, n=18):
    assert 0 <= x <= 1
    total = sum(((-1)**j*x**(2*j+1)/math.factorial(2*j+1) for j in range(n)), F(0))
    other = total + (-1)**n*x**(2*n+1)/math.factorial(2*n+1)
    return min(total, other), max(total, other)

def check():
    source = json.loads((D/'sources.json').read_text(encoding='utf-8'))
    for p, digest in source['local_sha256'].items():
        assert sha(ROOT/p) == digest, p
    assert sha(ROOT/'research_cognition_physics/archive_1234_/_shared/goal_prompt.txt') == source['goal_sha256']
    a, b = atan_bounds(F(1,5)), atan_bounds(F(1,239))
    pi_lo, pi_hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    assert F(314,100) < pi_lo < pi_hi < F(22,7)
    tl, tu = pi_lo/2, pi_hi/2
    kl = sin_bounds(2*pi_lo/65)[0]*F(65,16)
    ku = sin_bounds(2*pi_hi/65)[1]*F(65,16)
    assert F(39,100) < kl < ku < F(2,5)
    delta, sigma = F(1,12500), F(1,10000)
    e = delta+sigma
    phase_hi = F(31,20)*ku*tu
    assert phase_hi < F(341,350) < F(39,40) < pi_lo/2
    # Entire ideal-to-estimate line segment stays inside a fixed inversion domain.
    d_min = sin_bounds(F(17,100)*tl)[0] - 2*e
    d_max = sin_bounds(F(23,100)*tu)[1] + 2*e
    x_min = sin_bounds(F(9,20)*kl*tl)[0] - 4*e
    x_max = sin_bounds(phase_hi)[1] + 4*e
    assert 0 < d_min < d_max < F(19,50)
    assert 0 < x_min < x_max < F(17,20)
    assert F(28,3)**2 < 11**2*(1-F(19,50)**2)
    assert 4**2 < 8**2*(1-F(17,20)**2)
    ez, ey = F(11,4)*e+delta, F(11,2)*e+2*e+delta
    assert ez == F(23,40000) and ey == F(143,100000)
    # The old cycle witness: scope comparison, not new actual equal-calibration claim.
    assert 1-(ku*tu)**2/2 > F(4,5)
    gap = F(39,100)*F(1,5)*F(4,5)/4
    actual_gap = gap-2*delta
    assert gap == F(39,2500) and actual_gap == F(193,12500)
    assert actual_gap > 2*ez
    n = 300_000_000
    assert 2*n*sigma**2 == 6
    assert sum((F(6)**j/math.factorial(j) for j in range(40)), F(0)) > 400

    # Algebra diagnostics over a superset of realizable integral pairs.
    # No trajectory or actual device is claimed for arbitrary box points.
    t = math.pi/2
    k = math.sin(2*math.pi/65)/(16/65)
    max_residual = 0.0
    smallest_slack = 1.0
    for is_ in np.linspace(-t, t, 17):
        for it in np.linspace(-t, t, 17):
            tau = t+.15*is_
            s = k*(t+.35*is_)
            dd = k*.2*it
            pd = .5+.5*math.sin(tau/5)
            px = .5+.25*math.sin(s+dd)
            py = .5+.25*math.sin(s-dd)
            pz = .5+.25*math.sin(s)
            def predict(d,r):
                tauhat = 5*math.asin(2*d-1)
                shat = k*(t+7/3*(tauhat-t))
                xhat = math.asin(4*r-2)
                return .5+.25*math.sin(2*shat-xhat), .5+.25*math.sin(shat)
            yh,zh = predict(pd,px)
            max_residual=max(max_residual,abs(yh-py),abs(zh-pz))
            for sd in [-1,1]:
                for sx in [-1,1]:
                    yh,zh=predict(pd+sd*float(e),px+sx*float(e))
                    smallest_slack=min(smallest_slack,
                        float(ey)-abs(yh-py)-float(delta),
                        float(ez)-abs(zh-pz)-float(delta))
    assert max_residual < 1e-14 and smallest_slack > 0
    # Exact response coefficients: S = k[T+(nu/epsilon)(tau-T)].
    assert F(7,20)/F(3,20)==F(7,3)
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    z=np.diag([1.,-1.]); x=np.array([[0.,1.],[1.,0.]])
    assert max(abs(np.vdot(bell,np.kron(o,np.eye(2))@bell)) for o in [z,x]) < 1e-14
    assert max(abs(np.vdot(bell,np.kron(np.eye(2),o)@bell)) for o in [z,x]) < 1e-14
    s0=k*t; aa=k*.2
    witness={'clock':.5+.5*math.sin(t/5),'x':.5+.25*math.sin(s0),
             'base_z':.5+.25*math.sin(s0),'cycle_z':.5+.25*math.sin(s0-aa)}
    return {
      'status':'PASS_SHARED_RESPONSE_PREDICTION',
      'exact_certificate':{
        'parent_and_instrument_error':str(delta),'calibration_sampling_error':str(sigma),
        'clock_to_z_actual_probability_error':str(ez),
        'clock_x_to_y_actual_probability_error':str(ey),
        'phase_upper_bound':'341/350','clock_domain_upper':'19/50','signal_domain_upper':'17/20',
        'derivative_bounds':[11,8],
        'old_cycle_comparator_gap_lower':str(gap),
        'old_cycle_actual_gap_lower_conditional':str(actual_gap),
        'independent_repeats_per_calibration':n,'joint_sampling_failure_upper':'1/100'},
      'diagnostics':{'integral_box_points':17**2,'max_formula_residual':max_residual,
        'smallest_error_bound_slack':smallest_slack,'old_cycle_comparator_probabilities':witness},
      'scope':{'new_scientific_rounds':0,'new_axioms':0,'new_parent_simulation':False,
        'actual_experiment':False,'same_sample_intermediate_measurement':False,
        'full_geometry_reconstructed':False,'response_law_selected':False,
        'unknown_full_output_transfer':False,
        'physical_permissions':'Qualified reproducible known backgrounds and original terminal tasks only.'}}

if __name__=='__main__':
    result=check()
    if sys.argv[1:]==['--write']:
        with (D/'results.json').open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert not sys.argv[1:], 'Only --write or no argument'
        saved=json.loads((D/'results.json').read_text(encoding='utf-8'))
        assert result==saved, 'Saved result differs'
    print(json.dumps(result,ensure_ascii=False))
