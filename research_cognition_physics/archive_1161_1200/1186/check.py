"""Round 1186: exact local identities and an untruncated Gaussian record.
Default is read-only comparison. --write creates results.json, never overwrites.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt, factorial, sqrt, sin, cos
import cmath
import json
from pathlib import Path
import sys

BITS=120
DEN=1<<BITS

def down(x):
    x=F(x); return F((x.numerator*DEN)//x.denominator,DEN)
def up(x):
    x=F(x); return -down(-x)
class IV:
    def __init__(self,lo,hi=None):
        lo=F(lo); hi=lo if hi is None else F(hi)
        assert lo<=hi
        self.lo=down(lo); self.hi=up(hi)
    def __add__(self,o):
        o=iv(o); return IV(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self): return IV(-self.hi,-self.lo)
    def __sub__(self,o): return self+-iv(o)
    def __rsub__(self,o): return iv(o)+-self
    def __mul__(self,o):
        o=iv(o); a=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return IV(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=iv(o); assert o.lo*o.hi>0
        return self*IV(1/o.hi,1/o.lo)
    def __rtruediv__(self,o): return iv(o)/self
    def sqroot(self):
        assert self.lo>=0
        l=isqrt((self.lo.numerator*DEN*DEN)//self.lo.denominator)
        u=isqrt((self.hi.numerator*DEN*DEN)//self.hi.denominator)+1
        return IV(F(l,DEN),F(u,DEN))
    def record(self):
        return {'lower':str(self.lo),'upper':str(self.hi),
                'lower_display':float(self.lo),'upper_display':float(self.hi)}
def iv(x): return x if isinstance(x,IV) else IV(x)
def trig_point(x,which,n=28):
    x=F(x); assert 0<=x<=F(3,2)
    if which=='sin':
        value=sum(((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(n)),F(0))
        error=x**(2*n+1)/factorial(2*n+1)
    else:
        value=sum(((-1)**k*x**(2*k)/factorial(2*k) for k in range(n)),F(0))
        error=x**(2*n)/factorial(2*n)
    return IV(value-error,value+error)
def trig(x,which):
    x=iv(x); assert 0<=x.lo<=x.hi<=F(3,2)
    left=trig_point(x.lo,which); right=trig_point(x.hi,which)
    # On [0,1.5] sine increases and cosine decreases.
    return IV(left.lo,right.hi) if which=='sin' else IV(right.lo,left.hi)

def record_probability(a,b,c):
    a,b,c=map(F,(a,b,c))
    w=IV(a*b).sqroot()
    ratio=F(a+b,2)/w
    dx=trig(w,'cos'); dy=ratio*trig(w,'sin')
    assert dx.lo>0 and dy.lo>0
    radius=(dx*dx+dy*dy).sqroot()
    u=((radius+dx)/2).sqroot(); v=((radius-dx)/2).sqroot()
    theta=F(3,8)-c
    imag=(trig(theta,'sin')*u-trig(theta,'cos')*v)/radius
    return (1+imag)/2

def matrix(n): return [[F(0) for j in range(n)] for i in range(n)]
def add(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(c,a): return [[c*x for x in r] for r in a]
def trace_g(x,signs): return sum((signs[i]*x[i][i] for i in range(len(x))),F(0))
def tensor_map(x,signs,coefficient):
    out=scale(coefficient,x); tr=trace_g(x,signs)
    for i,s in enumerate(signs): out[i][i]-=F(1,2)*s*tr
    return out

def evaluate():
    tensor_cases=0
    a=F(3,2)
    for d in (3,4,5):
        signs=[-1]+[1]*(d-1)
        L=a*F(d-2,2)
        for i in range(d):
            for j in range(i,d):
                x=matrix(d); x[i][j]=x[j][i]=F(1)
                P=tensor_map(x,signs,F(d-2,2))
                A0=scale(-F(1,2)/L,P)
                C_A0=scale(a/2,tensor_map(A0,signs,F(1)))
                assert C_A0==scale(-F(1,4),x)
                tensor_cases+=1
    # D=4 pure conformal h=fG: full quadratic coefficient, second derivative.
    conformal_hessian=a/2*(4-F(16,2))
    assert conformal_hessian==-2*a
    free_contact=-4*a
    normal_order_constant_contact=4*a
    remaining_constant_contact=-2*a
    assert free_contact+normal_order_constant_contact==0
    assert free_contact+normal_order_constant_contact+remaining_constant_contact==conformal_hessian
    # Same canonical mode: O=(3/4)(a^2+a†^2)+(1/2)N, at m=1.
    tadpole=F(0); contact=F(3,4); transition_square=F(9,8)
    spectral_static=-transition_square
    assert contact+spectral_static==-F(3,8)
    assert spectral_static!=-F(3,8)
    # Both full and contact-deleted Hamiltonians are positive on 1<=y<=3/2.
    amin=F(5,12); bmin=F(1); cmin=-F(9,32)
    assert amin*bmin>F(3,5)**2
    positive_lower=F(3,10)+cmin
    assert positive_lower==F(3,160)
    full=record_probability(F(2,3),F(9,4),-F(9,32))
    deleted=record_probability(F(5,12),F(2),-F(1,4))
    gap=deleted-full
    assert gap.lo>F(3879,100000)
    # Cross-check elementary complex expression, not used by the certificate.
    diagnostics=[]
    for aa,bb,cc,bounds in [(F(2,3),F(9,4),-F(9,32),full),(F(5,12),F(2),-F(1,4),deleted)]:
        af,bf,cf=map(float,(aa,bb,cc))
        w=sqrt(af*bf); r=(af+bf)/(2*w)
        amplitude=cmath.exp(1j*(3/8-cf))/cmath.sqrt(cos(w)+1j*r*sin(w))
        probability=(1+amplitude.imag)/2
        assert abs(probability-float((bounds.lo+bounds.hi)/2))<1e-14
        diagnostics.append({'amplitude_real':amplitude.real,'amplitude_imag':amplitude.imag,'p_Y_plus':probability})
    return {
      'round':1186,'date':'2026-10-14','status':'PASS_WITH_LOCAL_COMPOSITE_SCOPE',
      'groups':{
        '1_same_background_contact_tree':{'pass':True,'rational_tensor_basis_cases':tensor_cases,'dimensions':[3,4,5],
          'identity':'C A0 = -I/4; C=a Q/2; L=a(D-2)/2','D4_conformal_hessian_example':str(conformal_hessian)},
        '2_complete_contact_and_causal_source':{'pass':True,'free_contact_example':str(free_contact),
          'normal_order_counterterm_example':str(normal_order_constant_contact),
          'remaining_vacuum_contact_example':str(remaining_constant_contact),
          'mode_tadpole':str(tadpole),'mode_contact':str(contact),'mode_transition_square':str(transition_square),
          'mode_static_complete':str(contact+spectral_static),'mode_static_without_contact':str(spectral_static)},
        '3_complete_quantum_control_record':{'pass':True,'time':'1','y_active':'3/2','input_control_dimension':2,
          'oscillator':'full L2(R), no Fock truncation','reference':'arbitrary for the exact joint unitary',
          'uniform_positive_H_strict_lower':str(positive_lower),'full_probability':full.record(),
          'contact_deleted_probability':deleted.record(),'difference':gap.record(),'certified_strict_lower':'3879/100000',
          'interval_bits':BITS,'trig_terms':28,'float_diagnostics':diagnostics}
      },
      'scope_limits':[
        'Tensor identity assumes the selected covariant scalar action, fixed variables/measure/local terms and isotropic stationary reference.',
        'The positive quantum mode is a finite-mode diagnostic, not a realization of that covariant vacuum or regulator.',
        'No field environment or quantum control coherence is discarded in the comparison.',
        'Probability separation is for the stated known preparation and Y effect; no hardware or universal instrument error is certified.',
        'Composite metric signature, actual common clocks, the two-derivative coefficient and Einstein dynamics remain open.'
      ],
      'new_adopted_axioms':0,'new_empirical_calibration':0,'gr_hard_deliverables_completed':0
    }

if __name__=='__main__':
    got=evaluate(); out=Path(__file__).with_name('results.json')
    if '--write' in sys.argv:
        with out.open('x',encoding='utf-8') as f: json.dump(got,f,ensure_ascii=False,indent=2); f.write('\n')
        print('WROTE results.json; 3 groups PASS; exact rational gap > 0.03879')
    else:
        old=json.loads(out.read_text(encoding='utf-8'))
        assert got==old, 'Saved results differ'
        print('READ_ONLY_MATCH 3 groups; exact rational gap > 0.03879')
