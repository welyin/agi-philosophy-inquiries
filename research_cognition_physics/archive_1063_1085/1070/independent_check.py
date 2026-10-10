"""Independent finite-window translation bound: standard library only.
Fraction certifies the piecewise-linear relaxation exactly. Decimal evaluates
an amplitude-damping witness. Samples check formulas, not the universal proof.
No author imports; no file writes at runtime.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib
import json

checks = 0

def assert_exact(condition):
    global checks
    checks += 1
    if not condition:
        raise AssertionError('exact certificate failed')

def rational_text(q):
    return str(q.numerator) if q.denominator == 1 else str(q)

def decimal(q):
    return Decimal(q.numerator) / Decimal(q.denominator)

# Rational unit vectors include 1D axis, 2D circle, and genuine 3D directions.
directions = [
    (F(0), F(0), F(1)), (F(0), F(0), F(-1)),
    (F(1), F(0), F(0)), (F(0), F(1), F(0)),
    (F(3,5), F(0), F(4,5)), (F(0), F(4,5), F(-3,5)),
    (F(2,3), F(2,3), F(1,3)),
    (F(-2,7), F(3,7), F(6,7)),
]
for vector in directions:
    assert_exact(sum(t*t for t in vector) == 1)

cases = [
    (F(1,2), F(1,4)),
    (F(2,5), F(1,3)),
    (F(1,10), F(7,8)),
    (F(3,4), F(1,4)),
    (F(1,4), F(3,4)),
    (F(0), F(1,4)),
    (F(2,3), F(0)),
    (F(0), F(1)),
]
case_results = []
max_decimal_kraus_residual = Decimal(0)
max_decimal_formula_residual = Decimal(0)

with localcontext() as context:
    context.prec = 80
    tolerance = Decimal('1e-70')
    for r, a in cases:
        assert_exact(0 <= r <= 1 and 0 <= a <= 1 and r+a <= 1)
        # Positivity alone constrains z-component parameters to a diamond.
        # Cutting at b=a makes the objective affine on both resulting polytopes.
        vertices = sorted(set([
            (F(-1),F(0)), (F(0),F(-1)), (F(0),F(1)), (F(1),F(0)),
            (a,1-a), (a,-(1-a)),
        ]))
        values = []
        for b,s in vertices:
            assert_exact(abs(b)+abs(s) <= 1)
            # Exact maximum over the two protected inputs +rz and -rz.
            plus = b+r*s-(a+r)
            minus = b-r*s-(a-r)
            objective = (abs(b-a)+r*abs(1-s))/2
            assert_exact(max(abs(plus),abs(minus))/2 == objective)
            values.append(objective)
        exact_minimum = min(values)
        assert_exact(exact_minimum == r*a/2)
        # Check the analytic nonnegative-slack identities independently over
        # rational b and allowed s. The identities printed below prove the
        # bound for every real feasible b,s, not merely these check points.
        bvalues = sorted(set([F(-1),F(-3,4),F(-1,5),F(0),a/2,a,(a+1)/2,F(1)]))
        identity_count = 0
        for b in bvalues:
            cap = 1-abs(b)
            for eta in [F(0),F(1,3),F(1,2),F(1)]:
                s = -cap+2*eta*cap
                slack = cap-s
                lhs = abs(b-a)+r*(1-s)-r*a
                if b <= 0:
                    rhs = (1-r)*a-(1+r)*b+r*slack
                elif b <= a:
                    rhs = (1-r)*(a-b)+r*slack
                else:
                    rhs = (1+r)*(b-a)+r*slack
                assert_exact(lhs == rhs and rhs >= 0)
                identity_count += 1

        rd, ad = decimal(r), decimal(a)
        root = (Decimal(1)-ad).sqrt()
        root_a = ad.sqrt()
        k = Decimal(1)-root
        # K0=diag(1,sqrt(1-a)); K1=[[0,sqrt(a)],[0,0]].
        # These establish trace preservation and the affine Bloch action.
        kraus_residual = abs(root*root+root_a*root_a-1)
        max_decimal_kraus_residual = max(max_decimal_kraus_residual,kraus_residual)
        assert_exact(kraus_residual < tolerance)
        assert_exact(-tolerance <= k <= ad+tolerance)
        # The general coefficient comparison is algebraic:
        # sqrt(1-a)>=1-a follows from t>=t^2 for t=1-a in [0,1].
        t=1-a
        assert_exact(t >= t*t)
        sample_errors=[]
        for vector in directions:
            x,y,z = [rd*decimal(q) for q in vector]
            error_squared = (k*k*(x*x+y*y)+ad*ad*z*z)/4
            error = error_squared.sqrt()
            # Exact difference of squared-error bound has nonnegative factor
            # r^2*(a^2-k^2)*(1-n_z^2)/4.
            theoretical_gap = rd*rd*(ad*ad-k*k)*(1-decimal(vector[2])**2)/4
            bound_gap = (rd*ad/2)**2-error_squared
            residual = abs(theoretical_gap-bound_gap)
            max_decimal_formula_residual = max(max_decimal_formula_residual,residual)
            assert_exact(residual < tolerance and error <= rd*ad/2+tolerance)
            sample_errors.append(str(error))
        assert_exact(abs(Decimal(sample_errors[0])-rd*ad/2) < tolerance)
        assert_exact(abs(Decimal(sample_errors[1])-rd*ad/2) < tolerance)
        case_results.append({
            'r':rational_text(r), 'a':rational_text(a),
            'r_plus_a':rational_text(r+a),
            'exact_minimum_trace_error':rational_text(exact_minimum),
            'positivity_relaxation_vertices':len(vertices),
            'nonnegative_slack_identities_checked':identity_count,
            'amplitude_damping_gamma':rational_text(a),
            'amplitude_Txy_decimal':str(root),
            'amplitude_Tz':rational_text(1-a),
            'amplitude_tz':rational_text(a),
            'attaining_inputs':['+r z','-r z'],
            'direction_sample_errors':sample_errors,
        })

result={
    'round':1070,
    'scope':'fixed deterministic qubit CPTP approximation on |x|<=r, target x+a z',
    'author_imported':False,
    'libraries':'Python standard library Fraction and Decimal only',
    'cases':case_results,
    'checks':checks,
    'max_decimal_kraus_residual':str(max_decimal_kraus_residual),
    'max_decimal_formula_residual':str(max_decimal_formula_residual),
    'analytic_lower_bound':'max D >= (|b-a|+r|1-s|)/2 >= r*a/2; |b|+|s|<=1',
    'analytic_upper_bound':'amplitude damping gamma=a; ||diag(-k,-k,-a)||=a, k=1-sqrt(1-a)<=a',
    'slack_definition':'u=1-|b|-s >=0',
    'slack_certificates':{
        'b<=0':'2 lower_error-r*a = (1-r)*a-(1+r)*b+r*u',
        '0<=b<=a':'2 lower_error-r*a = (1-r)*(a-b)+r*u',
        'a<=b':'2 lower_error-r*a = (1+r)*(b-a)+r*u',
    },
    'sample_errors_are_not_optimization_proof':True,
    'diamond_norm_or_reference_guarantee_claimed':False,
    'actual_spatial_displacement_derived':False,
    'passed':True,
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
print(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True))
