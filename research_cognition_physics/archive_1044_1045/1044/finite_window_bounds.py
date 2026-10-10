"""Rational verification of the explicit finite-window bounds in round 1044.

Default is read-only. No time-grid propagation or cutoff is used as a proof.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import math


def hyperbolic_bounds(x, parity, last_n):
    terms = [x**n / math.factorial(n) for n in range(parity, last_n + 1, 2)]
    first_omitted = x**(last_n + 2) / math.factorial(last_n + 2)
    ratio = x*x / ((last_n + 3)*(last_n + 4))
    return sum(terms), sum(terms) + first_omitted/(1-ratio)


def verify():
    checks = {}
    checks['cosh_174_below_3'] = hyperbolic_bounds(F(174,100),0,10)[1] < 3
    checks['sinh_101_below_1191'] = hyperbolic_bounds(F(101,100),1,9)[1] < F(1191,1000)
    checks['sinh_172_above_27'] = hyperbolic_bounds(F(172,100),1,9)[0] > F(27,10)
    klo,khi,qlo,qhi = F(9,10),F(11,10),F(88,100),F(174,100)
    qp = khi/qlo
    A = khi/(2*qlo)+qhi/(2*klo)
    Ap = 1/(2*qlo)+khi*khi/(2*qlo**3)+qp/(2*klo)+qhi/(2*klo*klo)
    checks['q_envelope'] = qlo*qlo < 2-F(101,100)**2 and qhi*qhi > 4-F(99,100)**2
    checks['A_envelope'] = A < F(16,10) and Ap < F(33,10)
    Dp = 3*qp+F(33,10)*3+F(16,10)*3*qp
    checks['D_derivative_below_20'] = Dp < 20
    pref=4/(2*klo*qlo)
    numerator_deriv = 3*pref*(1/klo+qp/qlo)+pref*3*qp
    checks['r_derivative_below_50'] = numerator_deriv+20 < 50
    checks['t_derivative_below_22'] = 1+20 < 22
    N = 3+(khi/qlo)*3
    Np = 3*qp+(1/qlo+khi*khi/qlo**3)*3+(khi/qlo)*3*qp
    checks['interior_amplitude'] = N < 7 and Np < 18 and 18+7*20 < 160
    h=F(1,50)
    checks['f_L1'] = 2*h/3 < F(12,100)**2
    checks['f_derivative_L1'] = 4*8/(3*h) < 24**2
    beta=F(198,100)
    tail2=(900/beta+F(12,100)**2/(6*F(99,100)**3))/3
    checks['half_axis_13'] = tail2 < 13**2  # pi > 3
    pointC=188/beta+2*F(84,100)/beta**2
    pointIn=F(2412,100)/beta+2*F(12,100)/beta**2
    pointOut=F(3012,100)/beta+2*F(12,100)/beta**2
    checks['interior_39'] = pointC**2 < 6*39**2  # 2*pi > 6
    checks['free_interior_6'] = pointIn**2 < 6*6**2
    checks['out_interior_7'] = pointOut**2 < 6*7**2
    xlo,xhi=F(99,100)**2,F(101,100)**2
    checks['g2_denominator'] = min(xlo*(2-xlo),xhi*(2-xhi)) > F(9995,10000)
    checks['g2_q_upper'] = 2-xlo < F(101,100)**2
    checks['g2_F_upper'] = F(1191,1000)**2/F(9995,10000) < F(142,100)
    checks['g4_denominator'] = xhi*(4-xhi) < F(304,100)
    checks['g4_q_lower'] = 4-xhi > F(172,100)**2
    checks['g4_F_lower'] = 4*F(27,10)**2/F(304,100) > F(95,10)
    T=10**12
    delta=F(100,10**6)+F(120,T)
    gap=F(19,21)-F(71,121)
    finite_gap=gap-4*delta
    checks['finite_probability_gap'] = finite_gap > F(317,1000)
    checks['delta'] = delta < F(101,10**6)
    common_record_lower=F(13456,23657)-2*delta
    checks['common_record_lower'] = F(9799,10000)/6+1 > F(29,25) and common_record_lower > F(55,100)
    E=F(1,10)**2/4+F(101,100)**2
    # sqrt(2(E+8))+1.06 < 6, checked after moving the positive 1.06.
    checks['momentum_first_moment_error'] = 2*(E+8) < (6-F(106,100))**2
    checks['momentum_second_moment_error'] = 2*(E+8)+F(106,100)**2 < 20
    assert all(checks.values()), checks
    return {
        'scope':'Rational envelope verification; the analytic proof supplies interval/time/reference quantifiers.',
        'checks':checks,'all_passed':all(checks.values()),
        'half_time_T':T,'total_evolution_time':2*T,
        'isometry_error_bound':{'rational':str(delta),'decimal':float(delta)},
        'asymptotic_probability_gap_lower':{'rational':str(gap),'decimal':float(gap)},
        'finite_probability_gap_lower':{'rational':str(finite_gap),'decimal':float(finite_gap)},
        'cq_half_diamond_bound':float(2*delta),
        'common_record_probability_lower':{'rational':str(common_record_lower),'decimal':float(common_record_lower)},
        'target_first_moment_error_bound':float(6*delta),
        'target_second_moment_error_bound':float(20*delta),
        'free_energy_support_upper':float(E),
        'pre_read_H0_norm_upper':float(E+8),
        'post_read_direct_sum_H_norm_upper':float(E+12),
        'continued_post_read_H0_norm_upper':float(E+16),
        'per_particle_collision_action_budget_upper':float(2*T*(E+12)),
        'not_claimed':['normalized rare-branch uniform error','force-square bound','GR','SM matching','all cognitive axioms','arbitrarily close parameters distinguished at this fixed T']
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    out=verify()
    body=json.dumps(out,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        Path(__file__).with_name('finite_window_bounds_results.json').write_text(body,encoding='utf-8')
    print(body)
