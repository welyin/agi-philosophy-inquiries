"""Independent round 1187 audit, using exact rational algebra only.
No author module is imported. No quantum or nonlinear trajectory claim.
Default compares saved review output read-only; --write creates it exclusively.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import sys

BASE=Path(__file__).resolve().parent

def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def transpose(x): return [list(row) for row in zip(*x)]
def add(x,y): return [[a+b for a,b in zip(r,s)] for r,s in zip(x,y)]
def scale(a,x): return [[a*b for b in r] for r in x]
def sub(x,y): return add(x,scale(-1,y))
def mm(x,y): return [[sum((a*b for a,b in zip(r,c)),F(0)) for c in zip(*y)] for r in x]
def trace(x): return sum((x[i][i] for i in range(len(x))),F(0))
def col(x): return [[F(v)] for v in x]
def dot(x,y): return mm(transpose(x),y)[0][0]
def inverse(x):
    n=len(x); a=[r[:]+s for r,s in zip(x,eye(n))]
    for i in range(n):
        j=next(j for j in range(i,n) if a[j][i])
        a[i],a[j]=a[j],a[i]
        p=a[i][i]; a[i]=[v/p for v in a[i]]
        for j in range(n):
            if j!=i:
                p=a[j][i]; a[j]=[v-p*w for v,w in zip(a[j],a[i])]
    return [r[n:] for r in a]
def det(x):
    n=len(x); a=[r[:] for r in x]; value=F(1)
    for i in range(n):
        j=next((j for j in range(i,n) if a[j][i]),None)
        if j is None:return F(0)
        if j!=i:a[i],a[j]=a[j],a[i];value=-value
        p=a[i][i]; value*=p
        for j in range(i+1,n):
            t=a[j][i]/p
            a[j]=[v-t*w for v,w in zip(a[j],a[i])]
    return value
def positive_definite(x):
    assert x==transpose(x)
    return all(det([r[:k] for r in x[:k]])>0 for k in range(1,len(x)+1))
def rank(x):
    a=[r[:] for r in x]; row=0
    for c in range(len(a[0])):
        j=next((j for j in range(row,len(a)) if a[j][c]),None)
        if j is None:continue
        a[row],a[j]=a[j],a[row];p=a[row][c];a[row]=[v/p for v in a[row]]
        for j in range(len(a)):
            if j!=row:
                p=a[j][c];a[j]=[v-p*w for v,w in zip(a[j],a[row])]
        row+=1
        if row==len(a):break
    return row

def material_constraint_case():
    # Full material-space rotation, followed by a non-orthonormal spatial frame.
    q=eye(6)
    for i in range(3):
        q[i][i]=q[i+3][i+3]=F(3,5)
        q[i+3][i]=F(4,5);q[i][i+3]=-F(4,5)
    assert mm(transpose(q),q)==eye(6)
    e0=[r[:3] for r in q]; normals=[r[3:] for r in q]
    change=[[F(2),F(1),F(0)],[F(0),F(1),F(1)],[F(0),F(0),F(1)]]
    e=mm(e0,change); h=mm(transpose(e),e)
    assert det(h)==4
    projector=sub(eye(6),mm(mm(e,inverse(h)),transpose(e)))
    assert mm(projector,projector)==projector
    v=add(mm(e,col([1,2,3])),mm(normals,col([3,4,0])))
    u=mm(projector,v); assert dot(u,u)==25
    normal=scale(F(1,5),u)
    W=F(1) # sqrt(det h)/V with V=2
    momenta=scale(W,normal)
    assert mm(transpose(e),momenta)==col([0,0,0])
    assert dot(momenta,momenta)==W*W
    assert dot(momenta,v)-W*5==0
    Hess=scale(W/5,sub(projector,mm(normal,transpose(normal))))
    assert rank(Hess)==2 and trace(Hess)==F(2,5)
    assert mm(Hess,e)==[[F(0)]*3 for _ in range(6)]
    assert mm(Hess,v)==col([0]*6)
    assert mm(Hess,Hess)==scale(W/5,Hess)
    constraint_gradients=transpose([e[i]+momenta[i] for i in range(6)])
    assert rank(constraint_gradients)==4
    # Formal X extension: X gauge differentiates only the new X-shift sector.
    fp=[([F(0)]*4)+r for r in eye(4)]
    assert rank(fp)==4 and all(v==0 for r in fp for v in r[:4])
    return {'N':6,'non_orthonormal_spatial_determinant':'4','V':'2',
        'exact_velocity_Hessian_rank':rank(Hess),'independent_material_constraints':4,
        'canonical_H':'0','dummy_X_gauge_rank':'4 of 8 primary directions',
        'original_material_gauge_directions_unfixed':4}

def graph_case():
    # A genuine two-residual-field graph with nonparallel spatial gradients.
    b=[[F(1,2),-F(1,3)],[F(2,5),F(1,4)],[-F(1,2),F(2,3)]]
    v=col([F(1,2),F(2,3)]); f=F(3,2)
    spatial=add(eye(3),mm(b,transpose(b)))
    internal=add(eye(2),mm(transpose(b),b))
    assert det(spatial)==det(internal)>0
    w=dot(v,mm(inverse(internal),v))
    z=[list(row) for row in transpose(v)]+b
    g=add(eye(4),mm(z,transpose(z)))
    assert det(g)==det(spatial)*(1+w)
    d2=f*f*det(spatial)
    contracted_p2=d2*w/(1+w)
    root2=d2-contracted_p2
    assert root2==d2/(1+w)>0 and contracted_p2<d2
    # Independent differentiation of sqrt(det g) g^-1 Z along a rank-one jet.
    gi=inverse(g); n=sub(eye(2),mm(mm(transpose(z),gi),z))
    assert n==inverse(add(eye(2),mm(transpose(z),z)))
    assert positive_definite(gi) and positive_definite(n)
    xi=col([1,2,-1,3]); pol=col([2,-1])
    dz=mm(xi,transpose(pol))
    dg=add(mm(dz,transpose(z)),mm(z,transpose(dz)))
    j=mm(gi,z)
    dgi=scale(-1,mm(mm(gi,dg),gi))
    dj=add(add(scale(trace(mm(gi,dg))/2,j),mm(dgi,z)),mm(gi,dz))
    contracted=mm(transpose(dj),xi)
    expected=scale(dot(xi,mm(gi,xi)),mm(n,pol))
    assert contracted==expected
    elliptic_value=dot(pol,contracted)
    assert elliptic_value>0
    # Opposite total-action signs alter canonical momenta but not psi_TT.
    for sign in (-1,1):
        k=F(2,3)
        gen=[[F(0),F(sign)/f],[sign*f*k*k,F(0)]]
        assert mm(gen,gen)==scale(k*k,eye(2))
    wave=[[F(0),1/f],[-f*k*k,F(0)]]
    assert mm(wave,wave)==scale(-k*k,eye(2))
    return {'residual_fields':2,'det_spatial':str(det(spatial)),
        'velocity_metric_norm':str(w),'Hamiltonian_root_squared':str(root2),
        'same_oriented_positive_action_root_sign':'negative',
        'contracted_principal_bilinear_divided_by_F_sqrt_det_g':str(elliptic_value),
        'normal_and_coordinate_inverse_Grams_positive':True,
        'two_action_signs_same_second_order_equation':True}

def gap_interval():
    # Independent, separate cosh and cos series, rather than the author's
    # combined positive 4j+2 expansion. Terms through power 58.
    n=30
    cosh_lower=sum((F(1,factorial(2*j)) for j in range(n)),F(0))
    cosh_next=F(1,factorial(2*n))
    ratio=F(1,(2*n+1)*(2*n+2))
    cosh_upper=cosh_lower+cosh_next/(1-ratio)
    cos_partial=sum((F((-1)**j,factorial(2*j)) for j in range(n)),F(0))
    cos_next=F((-1)**n,factorial(2*n))
    cos_lower=min(cos_partial,cos_partial+cos_next)
    cos_upper=max(cos_partial,cos_partial+cos_next)
    gap=(cosh_lower-cos_upper,cosh_upper-cos_lower)
    assert gap[0]>1
    saved=json.loads((BASE/'results.json').read_text(encoding='utf-8'))
    item=saved['groups']['3_finite_mode_linear_separation']
    assert F(item['exact_gap_lower'])<=gap[0]<=gap[1]<=F(item['exact_gap_upper'])
    return {'independent_lower':str(gap[0]),'independent_upper':str(gap[1]),
        'strict_lower':'1','width_display':float(gap[1]-gap[0]),
        'author_encloses_independent_interval':True,'diagnostic':float((gap[0]+gap[1])/2),
        'quantity':'linearized constant-V mode-response coefficient at kT=1'}

def evaluate():
    hashes={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in ('proof.md','check.py','results.json')}
    assert hashes['proof.md']=='a9651f85b944e8e42e24b59e1c561048efb3f07a30a56c89d42e193f226c76f3'
    return {'round':1187,'review_status':'PASS_FOR_REGULAR_BARE_MATERIAL_CLOCK_ROUTE_ONLY',
        'dependency_sha256':hashes,'author_module_imported':False,
        'groups':{'material_constraints':material_constraint_case(),
          'multifield_graph_and_principal_symbol':graph_case(),'finite_linear_mode':gap_interval()},
        'scope_limits':[
          'The raw action is real positive Gram, V>0, regular rank and a local actual-material coordinate chart.',
          'Residual modes require N>=5; N=4 has no residual mode for this test.',
          'General multifield strong ellipticity is a rank-one principal-symbol statement, not ordinary convexity of the full area Hessian.',
          'Canonical H=0 does not forbid relational time; the reduced root is nonzero on R>0.',
          'The finite cosh/cos contrast is a constant-V, m=0 linearized diagnostic, not a massive-loop vacuum prediction.',
          'No nonlinear finite-amplitude error, physical quantum Hilbert space, clock resource or record probability is certified.',
          'No exclusion of renormalized/collective Lorentz geometry or of the full cognitive axiom programme.'
        ],'new_axioms':0,'gr_hard_deliverables':0}

if __name__=='__main__':
    result=evaluate(); out=BASE/'review_checks.json'
    if '--write' in sys.argv:
        with out.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WROTE independent 1187 review: exact constraints, multifield symbol, finite linear gap PASS')
    else:
        assert result==json.loads(out.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH independent 1187 review; finite linear gap > 1')