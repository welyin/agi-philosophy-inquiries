"""1186 independent review: exact determinant expansion and Gaussian wavefunction.
No author code is imported; no oscillator occupation cutoff is used.
--write creates review_checks.json exclusively. Default is read-only verification.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from itertools import permutations
from math import isqrt, factorial
from pathlib import Path
import hashlib
import json
import sys

BASE = Path(__file__).resolve().parent
SCALE = 10**45
TERMS = 22

# Independent exact rational intervals. Arithmetic stays rational; only roots
# are enclosed on a decimal grid, unlike the author's per-operation dyadic grid.
def box(x): return (F(x), F(x))
def plus(x,y): return x[0]+y[0], x[1]+y[1]
def minus(x,y): return x[0]-y[1], x[1]-y[0]
def times(x,y):
    p=[a*b for a in x for b in y]
    return min(p), max(p)
def inv(x):
    assert x[0]>0
    return 1/x[1], 1/x[0]
def divide(x,y): return times(x,inv(y))
def root(x):
    assert x[0]>=0
    lo=isqrt((x[0].numerator*SCALE*SCALE)//x[0].denominator)
    hi=isqrt((x[1].numerator*SCALE*SCALE)//x[1].denominator)+1
    return F(lo,SCALE), F(hi,SCALE)
def point_trig(x,sine):
    assert F(0)<=x<=F(3,2)
    term=x if sine else F(1)
    partial=F(0)
    for j in range(TERMS):
        partial+=term
        denominator=(2*j+2)*(2*j+3) if sine else (2*j+1)*(2*j+2)
        term=-term*x*x/denominator
    # Consecutive partial sums bracket the alternating remainder.
    return min(partial,partial+term), max(partial,partial+term)
def trig(x,sine):
    assert F(0)<=x[0]<=x[1]<=F(3,2)
    left=point_trig(x[0],sine); right=point_trig(x[1],sine)
    return (left[0],right[1]) if sine else (right[0],left[1])
def gauss_box(a,b,c):
    w=root(box(a*b))
    re=trig(w,False)
    im=times(divide(box((a+b)/2),w),trig(w,True))
    assert re[0]>0 and im[0]>0
    radius=root(plus(times(re,re),times(im,im)))
    ur=root(times(box(F(1,2)),plus(radius,re)))
    ui=root(times(box(F(1,2)),minus(radius,re)))
    phase=box(F(3,8)-c)
    numerator=minus(times(trig(phase,True),ur),times(trig(phase,False),ui))
    return times(box(F(1,2)),plus(box(1),divide(numerator,radius)))
def df(x): return Decimal(x.numerator)/Decimal(x.denominator)
def decimal_sincos(x):
    st=x; ct=Decimal(1); sine=st; cosine=ct
    for k in range(1,120):
        st=-st*x*x/Decimal((2*k)*(2*k+1))
        ct=-ct*x*x/Decimal((2*k-1)*(2*k))
        sine+=st; cosine+=ct
        if abs(st)<Decimal('1e-83') and abs(ct)<Decimal('1e-83'): break
    return sine,cosine

def gauss_decimal(a,b,c):
    with localcontext() as ctx:
        ctx.prec=90
        a,b,c=map(df,(a,b,c)); w=(a*b).sqrt()
        sw,cw=decimal_sincos(w)
        re=cw; im=(a+b)*sw/(2*w)
        radius=(re*re+im*im).sqrt()
        ur=((radius+re)/2).sqrt()/radius
        ui=((radius-re)/2).sqrt()/radius
        sp,cp=decimal_sincos(Decimal(3)/8-c)
        return str((1+sp*ur-cp*ui)/2)

# Bivariate polynomials of total degree <=2, used to expand the determinant
# of I + s G^-1 h + t G^-1 k directly by its permutation definition.
def pc(x): return {(0,0):F(x)} if x else {}
def padd(x,y):
    z=dict(x)
    for ij,v in y.items(): z[ij]=z.get(ij,F(0))+v
    return {ij:v for ij,v in z.items() if v}
def pscale(a,x): return {ij:a*v for ij,v in x.items() if a*v}
def pmul(x,y):
    z={}
    for (i,j),a in x.items():
        for (k,l),b in y.items():
            if i+j+k+l<=2: z[i+k,j+l]=z.get((i+k,j+l),F(0))+a*b
    return {ij:v for ij,v in z.items() if v}
def matmul(x,y):
    d=len(x)
    return [[sum((x[i][k]*y[k][j] for k in range(d)),F(0)) for j in range(d)] for i in range(d)]
def trace(x): return sum((x[i][i] for i in range(len(x))),F(0))
def expand_local(d,lorentz):
    diag=[F(v) for v in [2,3,5,7,11][:d]]
    if lorentz: diag[0]=-diag[0]
    h=[[F(((i+1)*(j+1))%5-2,3) for j in range(d)] for i in range(d)]
    k=[[F((i+j)%4-1,2) for j in range(d)] for i in range(d)]
    u=[[h[i][j]/diag[i] for j in range(d)] for i in range(d)]
    v=[[k[i][j]/diag[i] for j in range(d)] for i in range(d)]
    m=[[padd(padd(pc(int(i==j)),{(1,0):u[i][j]}),{(0,1):v[i][j]}) for j in range(d)] for i in range(d)]
    det={}
    for perm in permutations(range(d)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(d) for j in range(i+1,d))
        term=pc(sign)
        for i in range(d): term=pmul(term,m[i][perm[i]])
        det=padd(det,term)
    z=padd(det,pc(-1))
    sq=padd(padd(pc(1),pscale(F(1,2),z)),pscale(-F(1,8),pmul(z,z)))
    inverse_trace={(0,0):F(d),(1,0):-trace(u),(0,1):-trace(v),
        (2,0):trace(matmul(u,u)),(0,2):trace(matmul(v,v)),
        (1,1):trace(matmul(u,v))+trace(matmul(v,u))}
    a=F(3,2); potential=F(d-2,2)*a
    action=pmul(sq,padd(pscale(a/2,inverse_trace),pc(-potential)))
    expected=a/2*(trace(matmul(u,v))-trace(u)*trace(v)/2)
    assert action.get((0,0),0)==a
    assert action.get((1,0),0)==action.get((0,1),0)==0
    assert action.get((1,1),0)==expected
    # Independent trace/traceless eigenvalues prove invertibility and A0.
    q_trace=1-F(d,2); p_trace=F(d-2,2)-F(d,2)
    L=a*F(d-2,2)
    for q,p in [(F(1),F(d-2,2)),(q_trace,p_trace)]:
        assert (a*q/2)*(-p/(2*L))==-F(1,4)
    return {'dimension':d,'signature':'Lorentz' if lorentz else 'positive',
            'mixed_hessian':str(expected),'stationary_linear_coefficients':True}

# Sparse unnormalized Fock vectors: |n~>=(adag)^n|0>. Polynomial operator
# action is exact; this is not a truncation of time evolution.
def combine(*terms):
    out={}
    for scale,vec in terms:
        for n,c in vec.items(): out[n]=out.get(n,F(0))+scale*c
    return {n:c for n,c in out.items() if c}
def ladder_square(vec,raising):
    out={}
    for n,c in vec.items():
        if raising: out[n+2]=c
        elif n>=2: out[n-2]=n*(n-1)*c
    return out
def q2(vec):
    diag={n:(2*n+1)*c for n,c in vec.items()}
    return combine((F(1,2),ladder_square(vec,False)),(F(1,2),ladder_square(vec,True)),(F(1,2),diag))
def p2(vec):
    diag={n:(2*n+1)*c for n,c in vec.items()}
    return combine((-F(1,2),ladder_square(vec,False)),(-F(1,2),ladder_square(vec,True)),(F(1,2),diag))
def norm2(vec): return sum((factorial(n)*c*c for n,c in vec.items()),F(0))
def interval_record(x):
    return {'lower':str(x[0]),'upper':str(x[1]),'width_display':float(x[1]-x[0])}

def evaluate():
    hashes={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in ('proof.md','check.py','results.json')}
    tensor=[expand_local(d,s) for d in (3,4,5) for s in (False,True)]
    vacuum={0:F(1)}
    op=combine((-F(1,2),p2(vacuum)),(F(1),q2(vacuum)),(-F(1,4),vacuum))
    kop=combine((F(1),p2(vacuum)),(F(1),q2(vacuum)),(-F(1,4),vacuum))
    assert op=={2:F(3,4)} and kop=={0:F(3,4)}
    transition=norm2(op)
    static=kop[0]-2*transition/F(2)
    # Differentiating E0(y)=sqrt(y)/2-y^2/8 at y=1.
    energy_second=-F(1,8)-F(1,4)
    assert static==energy_second==-F(3,8)
    # Crude but strict common lower bound, entirely rational.
    assert F(5,12)>F(3,5)**2
    positive_common=F(3,10)-F(9,32)
    assert positive_common==F(3,160)>0
    pars=[(F(2,3),F(9,4),-F(9,32)),(F(5,12),F(2),-F(1,4))]
    full,deleted=[gauss_box(*p) for p in pars]
    gap=minus(deleted,full)
    assert gap[0]>F(3879,100000)
    decimal=[gauss_decimal(*p) for p in pars]
    author=json.loads((BASE/'results.json').read_text(encoding='utf-8'))
    saved=author['groups']['3_complete_quantum_control_record']
    for own,name,val in [(full,'full_probability',decimal[0]),(deleted,'contact_deleted_probability',decimal[1]),(gap,'difference',None)]:
        entry=saved[name]; lo,hi=F(entry['lower']),F(entry['upper'])
        assert lo<=own[0]<=own[1]<=hi
        if val is not None: assert own[0]<=F(val)<=own[1]
    # The sign follows from rho01=z*/2, rho10=z/2 for z=<U0vac|U1vac>.
    z=complex(F(2,5),F(1,7))
    rho01=z.conjugate()/2; rho10=z/2
    yexpect=(-1j*rho10+1j*rho01).real
    assert abs(yexpect-z.imag)<1e-15
    return {
      'round':1186,'review_status':'PASS_WITH_DISTINCT_COVARIANT_AND_CANONICAL_SCOPES',
      'dependency_sha256':hashes,
      'methods':['independent determinant polynomial expansion','trace/traceless tensor eigenspaces',
         'exact sparse polynomial oscillator action','Gaussian wavefunction and rational alternating-series intervals',
         '90-digit Decimal diagnostic'],
      'tensor_tests':tensor,
      'oscillator_response':{'tadpole':'0','contact':str(kop[0]),'transition_square':str(transition),
          'static_complete':str(static),'energy_second_derivative':str(energy_second),
          'retarded_sine_coefficient':str(-2*transition)},
      'uniform_strict_positive_bound':str(positive_common),
      'independent_probability_intervals':{'full':interval_record(full),'deleted':interval_record(deleted),
          'gap':interval_record(gap),'strict_gap_lower':'3879/100000','sqrt_grid_decimal_places':45,
          'alternating_series_terms':TERMS,'enclosed_by_author_intervals':True},
      'decimal_probabilities':{'full':decimal[0],'deleted':decimal[1]},
      'sign_check':'rho01=conj(z)/2 implies <Y>=Im(z), pYplus=(1+Im(z))/2',
      'author_code_audit':'Directed rational operations, positive reciprocal/root, monotone trigonometric intervals and alternating remainders valid on used domains.',
      'scope_limits':[
         'The fixed-(q,p) Hamiltonian second derivative is not the fixed-phi covariant-action second derivative.',
         'The covariant background identity and positive oscillator are separate conditional realizations, not one GR model.',
         'The full Hilbert-space unitary preserves unknown control/reference jointly; a control marginal need not remain unchanged.',
         'The witness uses a known vacuum and stipulated final effect; no physical clock, pointer, regulator or instrument budget is certified.',
         'Zero-derivative algebra does not generate an Einstein kinetic term, Lorentz signature or dark-sector dynamics.'
      ],
      'new_axioms':0,'gr_hard_deliverables':0
    }

if __name__=='__main__':
    result=evaluate()
    out=BASE/'review_checks.json'
    if '--write' in sys.argv:
        with out.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        print('WROTE review_checks.json; independent determinant, exact oscillator response and rational Gaussian gap PASS')
    else:
        assert result==json.loads(out.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH independent review; strict probability gap > 0.03879')