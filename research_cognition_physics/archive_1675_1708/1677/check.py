"""Uniform a-priori bound on CKM gamma transport for the SM one-loop contract."""
import argparse
import json
import math
from pathlib import Path
from fractions import Fraction

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def gauge_integral(mass2,g1_coefficient,steps):
    # Integrate a piecewise constant *upper* bound on the gauge driving term.
    # y(s) <= y(0) exp(integral(2G)); right endpoints bound the quadrature
    # for exp(integral(2G)). All exported bounds are rounded out again.
    frac=Fraction
    dt=frac(205,1000*steps)
    precision=10**14
    free_growth=precision
    integral=frac(0)
    for k in range(1,steps+1):
        remaining=frac(204*(steps-k),1000*steps)
        g2=1/(1/frac('0.652')**2+frac(19,3)*remaining)
        g3=1/(1/frac('1.222')**2+14*remaining)
        upper=g1_coefficient*frac('0.59')**2+frac(9,4)*g2+8*g3
        step=2*upper*dt
        assert 0<step<1
        growth=1+step+step**2/(2*(1-step/3))
        free_growth=math.ceil(free_growth*growth)
        integral+=dt*frac(free_growth,precision)
    return mass2*integral


def solve():
    quarks=read('../1397/inputs.json')['high_scale']
    ew=read('../1369/inputs.json')
    T=math.log(1e16/91.1876)/(16*math.pi**2)
    assert .204<T<.205
    def above(row): return row[0]+2*row[1]
    def below(row): return row[0]-2*row[1]
    gauges=[above(ew[k]) for k in ('g1_GUT','g2','g3')]
    assert (1/gauges[0]**2-2*4.1*.205)>1/.59**2
    assert gauges[1]<.652 and gauges[2]<1.222
    u,d=quarks['up'],quarks['down']
    # The -9*y_t^4 term supplies the integrating denominator in the top comparison.
    high_squared=lambda row: (Fraction(str(row[0]))+2*Fraction(str(row[1])))**2
    int_u_free=gauge_integral(high_squared(u[2]),Fraction(17,20),512)
    int_d_free=gauge_integral(high_squared(d[2]),Fraction(1,4),512)
    # Rational lower Taylor bound on the exponent proves the chosen log bound.
    iu,idown=.0867,.000023
    x=9*Fraction(str(iu))
    assert sum(x**n/math.factorial(n) for n in range(18))>1+9*int_u_free
    assert Fraction(str(idown))>int_d_free
    inflation=math.exp(3*(iu+idown))*(1+1e-10)
    def ratio(rows,i):
        return math.ceil(1e8*inflation*(above(rows[i])/below(rows[i+1]))**2+1e-6)/1e8
    ru,rd=[ratio(u,i) for i in (0,1)],[ratio(d,i) for i in (0,1)]
    Fu=(1+max(ru))/(1-max(ru))
    Fd=(1+max(rd))/(1-max(rd))
    epsb=2*rd[1]/(1-rd[1])
    delta=1.5*(Fu*idown+Fd*iu)+1e-9
    ang=quarks['angles_theta12_theta13_theta23_delta'][:3]
    assert .225<below(ang[0])<above(ang[0])<.23
    assert .00400<below(ang[1])<above(ang[1])<.00438
    assert .0464<below(ang[2])<above(ang[2])<.0484
    # Coarse elementary bounds valid for every initial Dirac phase by triangle inequality.
    ub0=(.00399,.00438)
    cb0=(.0463,.0484)
    udlo=.972-delta
    cdlo=.222-delta
    ushi=.23+delta
    # Scalar differential inequality |x'| <= 1.5 U x + remainder.
    E=math.exp(1.5*iu)*(1+1e-12)
    rem=1.5*(Fu*idown+epsb*iu)+1e-10
    def col_b(x):
        return [x[0]/E-rem,E*(x[1]+rem)]
    ub,cb=col_b(ub0),col_b(cb0)
    assert min(udlo,cdlo,ub[0],cb[0])>0
    # For a given complex CKM entry, use a quartet sharing that entry to
    # bound |J|/|V_ia|^2. Both quark frame generators are off-diagonal.
    j_over_weights=[[ub[1]*cb[1]/udlo,ushi*cb[1]/ub[0]],
                    [ub[1]*cb[1]/cdlo,ushi*ub[1]/cb[0]]]
    gaps_u=[2*ru[0]/(1-ru[0]),2*Fu]
    gaps_d=[2*rd[0]/(1-rd[0]),epsb]
    contributions={}
    for i,row in enumerate(('u','c')):
        for j,col in enumerate(('d','b')):
            area=j_over_weights[i][j]
            dint=idown*(1 if j==0 else rd[1])
            rad=1.5*area*(gaps_u[i]*dint+gaps_d[j]*iu)
            contributions[row+col]=math.degrees(rad)+1e-10
    bound=math.ceil(sum(contributions.values())*1000+1e-7)/1000
    assert bound<.2
    saved=read('../1675/phase_box_results.json')
    excluded=saved['excluded_gamma_mod180_deg']
    measured=saved['experimental_one_width_deg']
    assert excluded[0]+bound<measured[0]<measured[1]<excluded[1]-bound
    examples=read('../1675/independent_results.json')
    assert max(abs(row['MZ_minus_high_in_deg_mod180']) for row in examples.values())<bound
    return {'all_input_phases':True,'int_largest_up_squared_upper':iu,
            'int_largest_down_squared_upper':idown,'adjacent_up_squared_ratio_bounds':ru,
            'adjacent_down_squared_ratio_bounds':rd,'running_CKM_magnitude_bounds':
            {'ub':ub,'cb':cb,'ud_min':udlo,'cd_min':cdlo,'us_max':ushi},
            'individual_angle_contributions_deg':contributions,'drift_absolute_bound_deg':bound,
            'weak_scale_excluded_open_interval_mod180_deg':[excluded[0]+bound,excluded[1]-bound],
            'direct_one_width_deg':measured}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('results.json')
        print('1677 PASS: all phases and all scale trajectories; complete two-width input box')
    else:
        with (HERE/'results.json').open('x',encoding='utf8') as output:
            output.write(json.dumps(result,indent=2)+'\n')
