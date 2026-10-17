"""Unnumbered 1261 native metric screen; read-only checks and stdout only."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import json
import math
import numpy as np

root=Path(__file__).resolve().parents[4]
path=root/"research_cognition_physics/archive_1234_/1261/parent_check.py"
spec=importlib.util.spec_from_file_location("frozen1261parent",path)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# Exact Q(sqrt(2)) pairs. This tests only a fixed-source test-shell replacement,
# not same-history geometry after changing the full massless stress source.
zero=(F(0),F(0))
one=(F(1),F(0))
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b): return (a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,q): return (a[0]*q,a[1]*q)
def weighted(seq):
    out=zero
    for j,a in enumerate(seq): out=add(out,scale(a,F((-1)**j,8)))
    return out

cs=[one,(F(0),F(1,2)),zero,(F(0),F(-1,2)),(F(-1),F(0)),
    (F(0),F(-1,2)),zero,(F(0),F(1,2))]
sn=[zero,(F(0),F(1,2)),one,(F(0),F(1,2)),zero,
    (F(0),F(-1,2)),(F(-1),F(0)),(F(0),F(-1,2))]
ab=[one,(F(0),F(1,2)),zero,(F(0),F(1,2)),one,
    (F(0),F(1,2)),zero,(F(0),F(1,2))]
basis=[[one]*8,cs,sn,[mul(c,c) for c in cs],[mul(s,s) for s in sn],
       [mul(c,s) for c,s in zip(cs,sn)]]
assert all(weighted(b)==zero for b in basis)
assert weighted(ab)==(F(1,4),F(-1,4))
assert F(7,5)**2<2<F(3,2)**2
# Therefore (sqrt(2)-1)/4 is strictly between 1/10 and 1/8.

mu,P,E,k,r,G=12.,5.,13.,2.,7.,.01
assert E*E==mu*mu+P*P
n=np.array([1.,0.,0.]); pb=P*n
errs=[]; values=[]
for j in range(8):
    phi=j*math.pi/4
    c=math.cos(phi); s=math.sin(phi)
    pc=k*np.array([c,0.,s])
    cbg=module.ordered_contraction(pb,pc,mu,0.,n)
    cgb=module.ordered_contraction(pc,pb,0.,mu,-n)
    expect_b=k*P*P/E*(4*abs(c)-2)
    expect_g=k*P*P/E*(3*c*c-1)
    pair=module.pair(pb,pc,mu,0.,n,r,G)
    expected=k*(-G*E/r*(2+(P/E)**2)+4*G*P*c/r
                -G*P*P/(4*r*E)*(4*abs(c)+3*c*c-3))
    errs.extend([abs(cbg-expect_b),abs(cgb-expect_g),abs(pair-expected)])
    values.append(E+k+pair)  # any constant free part has zero weighted sum
W=sum((-1)**j*x/8 for j,x in enumerate(values))
A=G*k*P*P/(r*E)
assert max(errs)<1e-12
assert abs(W-A*(math.sqrt(2)-1)/4)<1e-12
print(json.dumps({
 "status":"PASS_UNNUMBERED_NATIVE_METRIC_REPLACEMENT_SCREEN",
 "exact_quadratic_basis_cancelled":len(basis),
 "exact_abs_weight":["1/4","-sqrt(2)/4"],
 "W_over_A_strict_bounds":["1/10","1/8"],
 "max_pair_diagnostic_error":max(errs),
 "finite_source_E":E,
 "W_diagnostic":W,
 "scope":"bare final canonical fixed-B test-shell replacement only",
 "full_source_geometry_rejected":False,
 "anti_Wick_operator_obstruction_proved":False,
 "actual_instrument_certified":False
},ensure_ascii=False,indent=2))
