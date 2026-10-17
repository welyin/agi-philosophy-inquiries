"""1286 independent certificate by orthogonal radial-oscillator Jacobi ladders.
No author code is imported. H powers are exact finite expansions, not propagation
in a truncated Hilbert space. Default is read-only; --write creates JSON once.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import argparse, hashlib, json, math

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
DIMS=(3,2,1)

def add(*vectors):
    out={}
    for v in vectors:
        for n,c in v.items(): out[n]=out.get(n,F(0))+c
    return {n:c for n,c in out.items() if c}

def scale(v,c): return {n:x*c for n,x in v.items() if x*c}

def move(n,j,delta):
    m=list(n);m[j]+=delta
    return tuple(m)

# Basis product L_n^(d/2-1)(|q|^2/2) times the normalized Gaussian.
# Its norm is (d/2)_n/n!, so squared norms below are sums of positive terms.
def radial_square(v,j):
    out={};a=F(DIMS[j],2)-1
    for n,c in v.items():
        k=n[j]
        terms=[(n,2*(2*k+a+1)),(move(n,j,1),-2*(k+1))]
        if k: terms.append((move(n,j,-1),-2*(k+a)))
        for m,b in terms: out[m]=out.get(m,F(0))+c*b
    return {n:c for n,c in out.items() if c}

def kinetic(v,j):
    # T+|q|^2/8 is the radial oscillator with eigenvalue n+d/4.
    return add({n:c*(n[j]+F(DIMS[j],4)) for n,c in v.items()},
               scale(radial_square(v,j),F(-1,8)))

def components(v):
    return {
      'gauge_kinetic':kinetic(v,0),
      'scalar_kinetic':kinetic(v,1),
      'record_kinetic':kinetic(v,2),
      'gauge_scalar_interaction':scale(radial_square(radial_square(v,0),1),F(1,2)),
      'scalar_quartic':scale(radial_square(radial_square(v,1),1),F(119,1248)),
      'scalar_record_interaction':scale(radial_square(radial_square(v,1),2),F(1,2))}

def H(v): return add(*components(v).values())

@lru_cache(None)
def weight(n):
    answer=F(1)
    for k,d in zip(n,DIMS):
        for j in range(k): answer*=F(d+2*j,2*(j+1))
    return answer

def inner(v,w): return sum((c*w.get(n,F(0))*weight(n) for n,c in v.items()),F(0))
def norm2(v): return sum((c*c*weight(n) for n,c in v.items()),F(0))
def fs(x):
    x=F(x);return str(x.numerator)+'/'+str(x.denominator)
def matstr(m): return [[fs(x) for x in row] for row in m]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def ceil_sqrt(x):
    n=math.isqrt(x.numerator//x.denominator)
    if F(n*n)<x:n+=1
    assert x<=n*n and (n==0 or (n-1)**2<x)
    return n

class Dual:
    def __init__(self,value,grad=(0,0)):
        self.v=F(value);self.d=tuple(F(x) for x in grad)
    @staticmethod
    def cast(v):return v if isinstance(v,Dual) else Dual(v)
    def __add__(a,b):
        b=Dual.cast(b);return Dual(a.v+b.v,[x+y for x,y in zip(a.d,b.d)])
    __radd__=__add__
    def __neg__(a):return Dual(-a.v,[-x for x in a.d])
    def __sub__(a,b):return a+-Dual.cast(b)
    def __rsub__(a,b):return Dual.cast(b)+-a
    def __mul__(a,b):
        b=Dual.cast(b);return Dual(a.v*b.v,[x*b.v+a.v*y for x,y in zip(a.d,b.d)])
    __rmul__=__mul__
    def __truediv__(a,b):
        b=Dual.cast(b);return Dual(a.v/b.v,[(x*b.v-a.v*y)/b.v**2 for x,y in zip(a.d,b.d)])
    def __rtruediv__(a,b):return Dual.cast(b)/a

def source_checks(energy):
    rows=[]
    # Unreduced neutral conformal term, independently evaluated before changing variables.
    for j in range(1,5):
        a=F(j+3,3);ad=F(j-2,5);add=F(2-j,7)
        n=F(j+2,4);nd=F(j-3,9);r=F(j,7);rd=F(2-j,11)
        s=r/a;sd=rd/a-ad*r/a**2
        curvature=6/n**2*(add/a+ad**2/a**2-ad*nd/(a*n))
        old=a**3*sd**2/(2*n)-n*a**3*curvature*s*s/12
        red=a*rd**2/(2*n)
        boundary_dot=(add/(2*n)-ad*nd/(2*n*n))*r*r+ad*r*rd/n
        assert old-red+boundary_dot==0
        y=F(j+1,5)
        assert n*a**3*(y/a**2)*(r/a)**2/2==n*y*r*r/(2*a)
        aa=Dual(a,(1,0));nn=Dual(n,(0,1))
        leff=-nn*energy/aa
        rho=-leff.d[1]/a**3;pressure=leff.d[0]/(3*n*a*a)
        assert rho==energy/a**4 and pressure==rho/3
        rows.append({'metric_jet':j,'neutral_boundary_residual':'0/1',
                     'mixed_potential_residual':'0/1','rho':fs(rho),'pressure':fs(pressure)})
    mp=F(11,2);slope2=energy/(3*mp)
    assert slope2==F(4561,1287)
    # Both Einstein equations and continuity, after the independent source variation.
    assert 3*mp*slope2-energy==0
    assert -slope2+(energy+3*(energy/3))/(6*mp)==0
    assert -4*energy+3*(energy+energy/3)==0
    return {'independent_metric_jets':rows,'slope_squared':fs(slope2),
            'lapse_acceleration_continuity_residuals':['0/1']*3}

def calculate():
    code=[{(0,3,0):F(1)},{(0,0,0):F(4,5),(0,6,0):F(3,5)}]
    gram=[[inner(v,w) for w in code] for v in code]
    hgram=[[inner(v,H(w)) for w in code] for v in code]
    ygram=[[inner(v,radial_square(w,1)) for w in code] for v in code]
    energy=F(4561,78)
    assert gram==[[1,0],[0,1]]
    assert hgram==[[energy,0],[0,energy]]
    assert ygram==[[14,0],[0,F(266,25)]]
    # Orthogonal Jacobi supports prove bandwidth directly, without integrating monomials.
    for n in range(9):
        e={(0,n,0):F(1)};he=H(e);ye=radial_square(e,1)
        assert all(abs(m[1]-n)<=2 for m in he)
        assert all(abs(m[1]-n)<=1 for m in ye)
        assert inner(e,he)==F(357*n*n+1683*n+860,156)
    component_grams={name:[[inner(v,components(w)[name]) for w in code] for v in code]
                     for name in components(code[0])}
    for m in component_grams.values():assert m[0][1]==m[1][0]==0
    for j in range(2):assert sum(m[j][j] for m in component_grams.values())==energy
    cert=[]
    for j,v in enumerate(code):
        powers=[v]
        for k in range(4):powers.append(H(powers[-1]))
        norms=[norm2(w) for w in powers]
        bounds=[ceil_sqrt(q) for q in norms]
        fourth=sum(math.comb(4,k)*bounds[k]*bounds[4-k] for k in range(5))
        # Self-adjointness is checked across distinct H powers in the weighted basis.
        for k in range(3):assert inner(powers[k],powers[k+2])==norms[k+1]
        for w in powers:
            for axis in range(3):
                assert inner(w,radial_square(H(w),axis))==inner(radial_square(w,axis),H(w))
        cert.append({'code':j,'norm_squared':[fs(q) for q in norms],
                     'norm_integer_upper':bounds,'Jacobi_expansion_terms':[len(q) for q in powers],
                     'fourth_probability_derivative_upper':fourth})
    total=sum(c['fourth_probability_derivative_upper'] for c in cert)
    assert total==61450034258
    dy=ygram[0][0]-ygram[1][1]
    assert dy==F(84,25)
    pref=dy/4
    # e^(1/2)<sum_(n>=0)(1/2)^n=2, strictly, implies exp(-1/2)>1/2.
    t=F(1,10**6)
    coefficient=pref/2-F(total,24)*(2*t)**2
    assert coefficient==F(1229274982871,3000000000000)>F(2,5)
    uniform=F(2,5)*t*t
    assert uniform==F(1,2500000000000)
    point=pref*t*t/2-F(total,24)*t**4
    assert point>uniform
    deleted=[energy-component_grams['scalar_record_interaction'][j][j] for j in range(2)]
    assert deleted[0]!=deleted[1]
    # Saved author data is used ONLY as a final comparison, never to generate any result.
    author=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    for j,c in enumerate(cert):
        ac=author['exact_H_norm_certificates'][j]
        for key in ('norm_squared','norm_integer_upper','fourth_probability_derivative_upper'):
            assert c[key]==ac[key],('independent mismatch',j,key)
    assert author['normal_code']['compressed_H']==matstr(hgram)
    assert author['normal_code']['compressed_Y']==matstr(ygram)
    assert author['finite_record_certificate']['window_coefficient_lower']==fs(coefficient)
    files=[HERE/x for x in ('proof.md','check.py','results.json','NEXT.md')]+[BASE/'research_note_1286.md']
    return {'schema':'1286_independent_v1','round':1286,'verdict':'PASS_WITH_SCOPE',
      'independent_groups':1,
      'method':'orthogonal generalized-Laguerre radial oscillator Jacobi action; positive weighted norm sums, independent neutral metric jets; no author imports',
      'author_file_sha256':{str(p.relative_to(BASE)).replace('\\','/'):sha(p) for p in files},
      'independent_code_sha256':sha(Path(__file__).resolve()),
      'normal_code':{'Gram':matstr(gram),'compressed_H':matstr(hgram),'compressed_Y':matstr(ygram),
                     'H_bandwidth_at_most_2':True,'Y_bandwidth_at_most_1':True},
      'component_grams':{n:matstr(m) for n,m in component_grams.items()},
      'H_power_certificates':cert,
      'finite_response':{'difference':'p0-p1','leading_coefficient':fs(pref)+' * exp(-1/2)',
         'fourth_derivative_sum_upper':total,'window':[fs(t),fs(2*t)],
         'factorized_window_coefficient_lower':fs(coefficient),
         'strict_uniform_record_gap':fs(uniform),'strict_point_lower':fs(point),
         'deleted_coupling_record_independence':True,'deleted_mean_energies':[fs(x) for x in deleted]},
      'source_geometry':source_checks(energy),
      'scope':{'same_mean_energy_for_all_unknown_code_and_R':True,
        'joint_isometry_not_invariant_code_dynamics':True,'complete_homogeneous_unread_source':True,
        'declared_effective_quantum_model':True,'full_4D_quantum_invariant_sector':False,
        'actual_terminal_instrument_certified':False,'stable_classical_record_certified':False,
        'branch_semiclassical_geometry_certified':False,'new_axioms':0}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=calculate();target=HERE/'independent_results.json'
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='EXCLUSIVE_FIRST_WRITE'
    else:
        assert result==json.loads(target.read_text(encoding='utf-8')),'independent result mismatch'
        mode='READ_ONLY_MATCH'
    print('PASS_WITH_SCOPE: '+mode)
    print('M4 total = '+str(result['finite_response']['fourth_derivative_sum_upper']))
    print('uniform record gap > '+result['finite_response']['strict_uniform_record_gap'])
    print('independent_results_sha256 = '+sha(target))

if __name__=='__main__':main()
