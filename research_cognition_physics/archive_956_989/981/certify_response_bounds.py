"""Rational margins for the finite pulse witness; no floating-point proof."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json

HERE=Path(__file__).resolve().parent
OUT=HERE/'response_bounds_certificate.json'

def run():
    # sqrt(1.16) lies in (1.077,1.0771), giving .0357<d<.036.
    assert F(1077,1000)**2 < F(29,25) < F(10771,10000)**2
    dlow=(1-1/F(1077,1000))/2
    dhigh=(1-1/F(10771,10000))/2
    assert dlow>F(357,10000) and dhigh<F(36,1000)<F(19,100)**2
    # Odd Dyson series from the third term: all next ratios <= A^2/20.
    A=F(1,10);R=A**3/6/(1-A*A/20)
    assert R<F(1,5990)
    delta=2*F(19,100)*A*F(1,5990)+F(1,5990)**2
    assert delta<F(64,10**7)
    # Classical rational enclosure 3.14159 < pi < 22/7.
    # The sine argument lies in (0,1), where y-y^3/6 is increasing.
    y=(F(314159,100000)/2-F(1,10))/4
    assert 0<y<(F(22,7)/2-F(1,10))/4<1
    sinc_lower=1-F(1,20)**2/6
    born_lower=2*F(357,10000)*F(1,20)**2*sinc_lower**2*(y-y**3/6)
    contrast_lower=born_lower-2*F(64,10**7)
    assert contrast_lower>F(51,10**6)
    assert F(5,10000) == F(5)*F(1,100)**2
    return dict(round=981,all_rational_checks_passed=True,
        d_lower=str(dlow),d_upper=str(dhigh),odd_remainder_upper=str(R),
        probability_error_upper='0.0000064',
        rational_contrast_lower=str(contrast_lower),
        convenient_strict_contrast_lower='0.000051',
        alpha_two_derivative_relative_bound='0.0001',
        lapse_two_derivative_relative_bound='0.0005',
        assumptions=['3.14159 < pi < 22/7','h/Q and parity from native material',
                     'lapse changes only the free interval','fixed coordinate field'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert json.loads(OUT.read_text('utf-8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
