"""Exact algebra checks; no scalar spectrum, detector, or geometry is simulated."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

D = Path(__file__).resolve().parent
ROOT = next(p for p in D.parents if (p/'research_cognition_physics').is_dir())


def calc():
    baseline = json.loads((D/'source_baseline.json').read_text(encoding='utf-8'))
    for name, expected in baseline['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    # sin(A_t) sin(A): product of exponential coefficients, with -st*w exponent.
    gaussian_coefficients = {1:F(0),-1:F(0)}
    for s in [-1,1]:
        for t in [-1,1]:
            gaussian_coefficients[-s*t] += -F(s*t,4)
    assert gaussian_coefficients == {1:F(1,2),-1:F(-1,2)}
    # The worst angular contribution is -(1-Fc). cosh(H|t|) >= 2(2-Fc)/Fc.
    kernel_checks = 0
    minimum_margin = None
    for fc in [F(1,8),F(1,2),F(3,4),F(15,16)]:
        threshold = 2*(2-fc)/fc
        for factor in [F(1),F(3,2),F(2),F(4)]:
            ch = threshold*factor
            z_minus_one_lower = fc*ch-(2-fc)
            half_cosh_bound = fc*ch/2
            assert z_minus_one_lower >= half_cosh_bound > 0
            margin = z_minus_one_lower-half_cosh_bound
            minimum_margin = margin if minimum_margin is None else min(minimum_margin,margin)
            kernel_checks += 1
    # cosh(log z) >= z/2: combines with prior bound to give 4 e^-H|t|/Fc.
    for z in [F(1),F(2),F(7),F(19,3)]:
        assert (z+1/z)/2 >= z/2
    # Reuse 1232: m*=2c/T -> A=2*c0*c*T and u=8/(c0*c*r^4*T^3).
    c0 = F(23,24)**4
    # These are algebra substitutions, NOT spectral data or an instrument certificate.
    substitutions = 0
    for c, T, r in [(F(1,10),F(8),F(1,2)),(F(3,7),F(13),F(2,5))]:
        A = c0*T*T*(2*c/T)
        u = (16/(r**4*T*T))/A
        assert A == 2*c0*c*T
        assert u == 8/(c0*c*r**4*T**3)
        rho = F(1,10)
        lam2 = min(1/T**2,rho*A/T**4)/2
        assert 0 < lam2*T*T <= 1
        dyson_upper = lam2*lam2*T**4/2
        op_budget = rho*lam2*A/2
        assert op_budget > 0 and dyson_upper+op_budget <= rho*lam2*A
        substitutions += 1
    return {'status':'PASS_WITH_SCOPE',
            'exact':{'gaussian_exp_w_coefficient':str(gaussian_coefficients[1]),
                     'gaussian_exp_minus_w_coefficient':str(gaussian_coefficients[-1]),
                     'kernel_guard_checks':kernel_checks,
                     'minimum_kernel_guard_margin':str(minimum_margin),
                     'generic_window_substitutions':substitutions,
                     'c0':str(c0)},
            'not_computed':['specific BD spectrum or its numeric density',
                            'actual apparatus/source completion',
                            'finite-G geometry remainder'],
            'proof_scope':'same-material positive-frequency spectral existence and generic finite comparison window'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--write',action='store_true')
    args=ap.parse_args()
    result=calc()
    if args.write:
        with (D/'results.json').open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        mode='created_once'
    else:
        assert json.loads((D/'results.json').read_text(encoding='utf-8')) == result
        mode='read_only_match'
    print(json.dumps({'mode':mode,**result},ensure_ascii=False))


if __name__=='__main__':
    main()
