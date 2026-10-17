"""Exact independent finite-symbol check. No apparatus or quantum claim."""
from fractions import Fraction as F
import json

# Q(sqrt(3)); independent signs cancel off-diagonal/linear metric terms.
def q(a=0, b=0): return (F(a), F(b))
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b): return (a[0]*b[0]+3*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,w): return (a[0]*w,a[1]*w)
points=[(q(),q(s),F(3,2)) for s in [-1,1]]
points += [(q(F(s,2)),q(0,F(t,2)),F(-1)) for s in [-1,1] for t in [-1,1]]
points += [(q(s),q(),F(1,2)) for s in [-1,1]]
assert len(points)==8 and sum(abs(w) for _,_,w in points)==8
for x,z,w in points: assert add(mul(x,x),mul(z,z))==q(1)
basis=[lambda x,z:q(1),lambda x,z:x,lambda x,z:z,
       lambda x,z:mul(x,x),lambda x,z:mul(z,z),lambda x,z:mul(x,z)]
for fn in basis:
    total=q()
    for x,z,w in points: total=add(total,scale(fn(x,z),w))
    assert total==q()

mu,P,E,k,r,G=map(F,[12,5,13,2,7,'0.01'])
assert E*E==mu*mu+P*P and P/E<F(1,2)
A=G*k*P*P/(r*E)
weighted=F(0)
for x,z,w in points:
    c=x[0]
    assert x[1]==0
    y=abs(c)
    # Original contracted numerator before the special-angle factorization.
    aa=P*P; a=P; d=P*c
    Q=aa-3*a*a+8*a*c*d+y*y*(a*a-3*aa)+2*y*(2*d*d-aa+a*a)
    C_source_gamma=k/E*Q/(1+y)**2
    assert C_source_gamma==k*P*P/E*(4*y-2)
    # Independent substitution into the reverse massive-source bracket, y_m=1.
    pa2=k*k; pb2=P*P; pan=k*c; pbn=P; dot=k*P*c
    bracket=(2*(2*dot*dot*pbn*pbn-2*pan*pbn*dot*pb2+pan*pan*pb2*pb2-dot*dot*pb2)/(E*E)
             +2*(-pa2*pbn*pbn+pan*pan*pbn*pbn+2*pan*pbn*dot+dot*dot-pan*pan*pb2)
             +(-3*pa2*pbn*pbn+pan*pan*pbn*pbn+8*pan*pbn*dot+pa2*pb2-3*pan*pan*pb2))
    C_reverse=bracket/(4*k*E)
    assert C_reverse==k*P*P/E*(3*c*c-1)
    pair=-G*E*k/r*(2+P*P/(E*E))+4*G*k*P*c/r-G/(4*r)*(C_source_gamma+C_reverse)
    weighted += w*pair
assert weighted==A==F(1,182)
assert weighted/8==F(1,1456)
print(json.dumps(dict(status='PASS_EXACT_RATIONAL_EIGHT_DIRECTION_SCREEN',points=8,
    annihilated_basis_functions=len(basis),weight_l1='8',weighted_pair=str(weighted),
    bare_fixed_metric_max_residual_lower=str(weighted/8),
    full_source_geometry_rejected=False,AW_operator_lower_bound=False,
    actual_instrument_or_new_round=False),indent=2))
