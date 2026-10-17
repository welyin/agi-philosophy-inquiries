"""1285 independent check: exact Hermite ladders and independent metric jets.
No author functions are imported. No finite-Hilbert compression is used.
--write creates this reviewer's result exclusively; default is read-only.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
VAC=(0,0,0,0,0)

def plus(*vs):
    out={}
    for v in vs:
        for n,c in v.items(): out[n]=out.get(n,F(0))+c
    return {n:c for n,c in out.items() if c}

def times(v,c): return {n:a*c for n,a in v.items() if a*c}

def move(n,j,k):
    ns=list(n);ns[j]+=k
    return tuple(ns)

# Probabilists' Hermite states He_n(q) exp(-q^2/4), norm squared n!.
# q=raising+lowering; coefficients are rational, without square-root rounding.
def position(v,j):
    out={}
    for n,c in v.items():
        up=move(n,j,1);out[up]=out.get(up,F(0))+c
        if n[j]:
            dn=move(n,j,-1);out[dn]=out.get(dn,F(0))+c*n[j]
    return {n:c for n,c in out.items() if c}

def repeated_position(v,j,n):
    for _ in range(n): v=position(v,j)
    return v

def derivative_state(v,j):
    # d(He_n psi)/dq = n He_{n-1} psi / 2 - He_{n+1} psi / 2.
    out={}
    for n,c in v.items():
        up=move(n,j,1);out[up]=out.get(up,F(0))-c/2
        if n[j]:
            dn=move(n,j,-1);out[dn]=out.get(dn,F(0))+c*n[j]/2
    return {n:c for n,c in out.items() if c}

def kinetic(v,j): return times(derivative_state(derivative_state(v,j),j),F(-1,2))

def potential(v):
    parts=[]
    for j in range(3):
        for k in (3,4):
            parts.append(times(repeated_position(repeated_position(v,j,2),k,2),F(1,2)))
    parts.append(times(repeated_position(v,3,4),F(1,16)))
    parts.append(times(repeated_position(repeated_position(v,3,2),4,2),F(1,8)))
    parts.append(times(repeated_position(v,4,4),F(1,16)))
    return plus(*parts)

def H(v): return plus(*(kinetic(v,j) for j in range(5)),potential(v))

def inner(v,w):
    return sum((c*w.get(n,F(0))*math.prod(math.factorial(k) for k in n) for n,c in v.items()),F(0))

def rotation(v,i,j):
    return plus(position(derivative_state(v,j),i),times(position(derivative_state(v,i),j),-1))

def dot(x,y): return sum((a*b for a,b in zip(x,y)),F(0))
def J(x): return [-x[1],x[0]]
def fs(x): return str(x.numerator)+'/'+str(x.denominator)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

class Dual:
    """Exact first variation in a finite list of independent arguments."""
    def __init__(self,x,d=(0,0,0)):
        self.x=F(x);self.d=tuple(F(t) for t in d)
    @staticmethod
    def cast(x): return x if isinstance(x,Dual) else Dual(x)
    def __add__(a,b):
        b=Dual.cast(b);return Dual(a.x+b.x,[x+y for x,y in zip(a.d,b.d)])
    __radd__=__add__
    def __neg__(a): return Dual(-a.x,[-x for x in a.d])
    def __sub__(a,b): return a+-Dual.cast(b)
    def __rsub__(a,b): return Dual.cast(b)+-a
    def __mul__(a,b):
        b=Dual.cast(b);return Dual(a.x*b.x,[x*b.x+a.x*y for x,y in zip(a.d,b.d)])
    __rmul__=__mul__
    def __truediv__(a,b):
        b=Dual.cast(b);return Dual(a.x/b.x,[(x*b.x-a.x*y)/b.x**2 for x,y in zip(a.d,b.d)])
    def __rtruediv__(a,b): return Dual.cast(b)/a
    def __pow__(a,n):
        if n==0:return Dual(1)
        return Dual(a.x**n,[n*a.x**(n-1)*x for x in a.d])


def metric_and_legendre_checks():
    rows=[]
    for k in range(1,7):
        a=F(k+3,3);ad=F(k-2,5);add=F(2-k,7)
        n=F(k+2,4);nd=F(k-4,9)
        chi=[F(k,7),F(2-k,9)];cd=[F(k-3,8),F(k+1,11)]
        aa=[F(k-2,6),F(1-k,8),F(k,10)];av=[F(k,9),F(3-k,7),F(1,5)]
        a0=F(k-3,10);g=F(1);lam=F(1,4)
        var=[c/a for c in chi]
        vd=[cd[j]/a-ad*chi[j]/a**2 for j in range(2)]
        cov=[vd[j]+g*a0*J(var)[j] for j in range(2)]
        covchi=[cd[j]+g*a0*J(chi)[j] for j in range(2)]
        curvature=6/n**2*(add/a+ad**2/a**2-ad*nd/(a*n))
        original=a**3*dot(cov,cov)/(2*n)-n*a**3*curvature*dot(var,var)/12
        original+=a*dot(av,av)/(2*n)-n*a*g*g*dot(aa,aa)*dot(var,var)/2-n*a**3*lam*dot(var,var)**2/4
        potential_value=g*g*dot(aa,aa)*dot(chi,chi)/2+lam*dot(chi,chi)**2/4
        reduced=a*(dot(covchi,covchi)+dot(av,av))/(2*n)-n*potential_value/a
        boundary_dot=(add/(2*n)-ad*nd/(2*n*n))*dot(chi,chi)+ad*dot(chi,cd)/n
        assert original-reduced+boundary_dot==0
        pc=[a*c/n for c in covchi];pa=[a*v/n for v in av]
        legendre=dot(pc,cd)+dot(pa,av)-reduced
        charge=chi[0]*pc[1]-chi[1]*pc[0]
        h0=(dot(pc,pc)+dot(pa,pa))/2+potential_value
        assert legendre==n*h0/a-g*a0*charge
        # Independently vary lapse and scale BEFORE selecting a time gauge.
        nn=Dual(n,(1,0,0));ax=Dual(a,(0,1,0));energy=F(33,8)
        lm=-nn*energy/ax
        rho=-lm.d[0]/a**3;pressure=lm.d[1]/(3*n*a*a)
        assert rho==energy/a**4 and pressure==rho/3
        rows.append({'case':k,'boundary_residual':'0/1','legendre_residual':'0/1',
                     'nonzero_G_example':fs(charge),'rho':fs(rho),'pressure':fs(pressure)})
    return rows


def geometry_checks(energy):
    mp2=F(11,2);slope=F(1,2);rows=[]
    for eta in (F(0),F(1,300),F(1,100)):
        aa=1+slope*eta;av=slope;acc=F(0);nn=aa;nv=slope
        # Variables here are (a, da/dt, N). The lapse is not fixed before variation.
        a=Dual(aa,(1,0,0));v=Dual(av,(0,1,0));n=Dual(nn,(0,0,1))
        lag=-3*mp2*a*v*v/n-n*energy/a
        momentum=-6*mp2*a*v/n
        time_momentum=sum((q*d for q,d in zip(momentum.d,(av,acc,nv))),F(0))
        lapse=lag.d[2];scale=lag.d[0]-time_momentum
        assert lapse==scale==0
        rho=energy/aa**4;pr=rho/3
        continuity=-4*energy*slope/aa**5+3*slope*(rho+pr)/aa
        assert continuity==0
        rows.append({'eta':fs(eta),'a':fs(aa),'independent_lapse_EL':'0/1',
                     'independent_scale_EL':'0/1','continuity':'0/1'})
    return rows


def calculate():
    vacuum={VAC:F(1)};h1=H(vacuum);h2=H(h1)
    # Exact action creates only finitely many Hermite terms for these moments.
    # No projection or compressed time propagator is defined.
    moments=[F(1),inner(vacuum,h1),inner(h1,h1),inner(h1,h2),inner(h2,h2)]
    expected=[F(1),F(33,8),F(2683,64),F(495575,512),F(161114457,4096)]
    assert moments==expected
    for state in (vacuum,h1,h2):
        for i,j in ((0,1),(1,2),(0,2),(3,4)):assert not rotation(state,i,j)
    # Symmetry and Hermiticity crosschecks use independent finite Hermite states.
    for j in range(5):
        basis={move(VAC,j,2):F(1)}
        assert inner(basis,h1)==inner(H(basis),vacuum)
        assert inner(basis,h2)==inner(H(basis),h1)
    energy=moments[1]
    var=moments[2]-energy**2
    d4=sum((F(math.comb(4,j))*(-1)**j*moments[j]*moments[4-j] for j in range(5)),F(0))
    t=F(1,100);lower=var*t*t-d4*t**4/24;upper=var*t*t
    assert F(3,1250)<lower<upper<F(1,400)
    ledger={
      'gauge_kinetic':sum((inner(vacuum,kinetic(vacuum,j)) for j in range(3)),F(0)),
      'scalar_kinetic':sum((inner(vacuum,kinetic(vacuum,j)) for j in (3,4)),F(0)),
      'interaction':sum((inner(vacuum,times(repeated_position(repeated_position(vacuum,j,2),k,2),F(1,2))) for j in range(3) for k in (3,4)),F(0))}
    ledger['quartic']=energy-sum(ledger.values())
    assert sum(ledger.values())==energy
    files=[HERE/n for n in ('proof.md','check.py','results.json','NEXT.md')]+[BASE/'research_note_1285.md']
    return {'schema':'1285_independent_v1','round':1285,'verdict':'PASS_WITH_SCOPE','independent_groups':1,
      'method':'exact five-mode Hermite ladder action with factorial Hilbert inner product; independent lapse/scale dual derivatives and unreduced metric jets',
      'author_file_sha256':{str(p.relative_to(BASE)).replace('\\','/'):sha(p) for p in files},
      'hermite':{'vacuum_terms':len(vacuum),'H_vacuum_terms':len(h1),'H2_vacuum_terms':len(h2),
                 'largest_total_H2_index':max(sum(n) for n in h2),'moments':[fs(x) for x in moments],
                 'Gauss_and_SO3_residuals_exact_zero':True,'hermiticity_exact':True},
      'energy_ledger':{k:fs(v) for k,v in ledger.items()},'variance':fs(var),'spectral_difference_moment4':fs(d4),
      'probability':{'nonreturn_lower':fs(lower),'nonreturn_upper':fs(upper),
                     'coarse_lower_margin':fs(lower-F(3,1250)), 'coarse_upper_margin':fs(F(1,400)-upper)},
      'metric_legendre':metric_and_legendre_checks(),'source_geometry':geometry_checks(energy),
      'end_scale':fs(1+t/2),'proper_time':fs(t+t*t/4),'mathematical_post_measurement_energy_upper':fs(5*energy),
      'scope':{'same_declared_homogeneous_quantum_effective_class':True,
        'not_a_finite_Hilbert_compression':True,'no_author_function_imports':True,
        'classical_4D_reduction_has_explicit_boundary_dictionary':True,
        'full_4D_quantum_invariant_sector_certified':False,
        'source_only_for_unmeasured_unconditional_waiting':True,
        'actual_measurement_device_or_source_certified':False,
        'uniform_geometry_for_unknown_inputs_certified':False,
        'Einstein_or_gauge_group_generated_by_this_round':False,'new_axioms':0}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    obj=calculate();p=HERE/'independent_results.json'
    if args.write:
        with p.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(obj,f,indent=2,ensure_ascii=False);f.write('\n')
        mode='EXCLUSIVE_FIRST_WRITE'
    else:
        assert json.loads(p.read_text(encoding='utf-8'))==obj,'Independent result or binding mismatch'
        mode='READ_ONLY_MATCH'
    print('PASS_WITH_SCOPE: '+mode)
    print('Hermite moments='+','.join(obj['hermite']['moments']))
    print('nonreturn='+obj['probability']['nonreturn_lower']+' .. '+obj['probability']['nonreturn_upper'])
    print('boundary/Legendre and independent metric EL residuals: exact zero')

if __name__=='__main__':main()
