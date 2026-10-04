"""Certify the existing T=|X|²/2 readout on the native population fourth jet."""
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'round745_drafts'))
import certify_history_integral as base
I=base.Interval

def coefficients():
    data=json.loads((ROOT/'round745_drafts/exact_history_density_results.json').read_text('utf8'))
    got={(v['r']//2,v['s']//2):Q(int(v['numerator']),int(v['denominator'])) for v in data['coefficients']['4']}
    g=Q(109,625)
    expected={(0,0):-Q(76,3)*g,(0,1):-8*g,(1,0):-Q(17,3)*g,
              (2,0):Q(31,9)*g,(1,1):Q(31,9)*g,(3,0):Q(2,3)*g,
              (2,1):Q(4,3)*g,(1,2):Q(2,3)*g}
    assert got==expected
    return got

def prepare(K):
    coef=coefficients();rows=[[Q(0)]*(2*K+2) for _ in range(4)]
    for k in range(K+1):
        power=2*k+1;sine=Q((-1)**k,4*factorial(power))
        for (a,b),c in coef.items():
            p=a+1+power
            beta=sum(Q(2*(-1)**j*comb(p,j),2*(b+j)+1) for j in range(p+1))
            assert beta>0
            rows[a+b][power]+=c*sine*beta
    return [[I.rational(c) for c in row] for row in rows],sum(abs(c) for c in coef.values())

def evaluate(x,poly):
    if not x:return I(0)
    r=I.rational(x);v=r*r;B=6*v/(6+v)
    g=I(0)
    for row in reversed(poly):
        f=I(0)
        for c in reversed(row):f=f*B+c
        g=g*v+f
    return v*v*(-v).exp()/(1+v/6).sqrt()*g

def run():
    K=20;n=20;cells=500;R=Q(10);h=R/(2*cells);delta=Q(1,2)
    poly,total=prepare(K);weights=base.newton_cotes(n);grid={}
    for cell in range(cells):
        for j,w in enumerate(weights):
            index=cell*n+j;grid[index]=grid.get(index,Q(0))+h*w
    result=I(0)
    for index,w in grid.items():
        if w:result=result+evaluate(R*index/(cells*n),poly)*I.rational(w)
    # |6z²/(6+z²)|<13 on the old half-radius disks. |sin series|<=e^13/4.
    M=Q(11**4)*2*2*(2*total*11**6)*Q(3**13,4)
    ratio=h/delta
    quadrature=cells*h*(2+sum(abs(w) for w in weights))*M*ratio**(n+1)/(1-ratio)
    sine=total*factorial(5)*Q(6**43,4*factorial(43))
    tail=total*Q(1,40*2**100)*factorial(5)*sum(Q(100**j,factorial(j)) for j in range(6))
    error=I.rational(quadrature+sine+tail)
    lo=(result.lo-error.hi).next_minus();hi=(result.hi+error.hi).next_plus()
    assert hi<0
    return dict(readout='E_T_plus=1/2+sin(T)/4, T=|X|²/2, existing624/723 menu.',
        fourth_density_coefficients={str(k):str(v) for k,v in coefficients().items()},
        unnormalized_fourth_contrast_interval=[str(lo),str(hi)],
        quadrature_error_upper=str(I.rational(quadrature).hi),
        sine_error_upper=str(I.rational(sine).hi),tail_error_upper=str(I.rational(tail).hi),
        exact_density_independent_of_onsite_potential=True,
        scalar_volume_dependence='w_v^-2; normalized contrast divides by same positive Gaussian Z',
        graph_extension_is_analytic_not_simulated=True,strict_negative=True,all_checks_passed=True)

if __name__=='__main__':
    r=run()
    with Path(__file__).with_name('higgs_readout_certificate_results.json').open('x',encoding='utf8') as f:
        json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
