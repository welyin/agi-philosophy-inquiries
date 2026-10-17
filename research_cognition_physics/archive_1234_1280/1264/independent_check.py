"""1264 independent finite audit: no author code/results imported.
Two groups: same-source representation and common finite menu.
Default is read-only. --write creates the result once, never overwrites it.
Uniform operator claims are established in the reviewed proofs; rational
mass-shell diagnostics below are finite checks, not substitutes for proofs.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json
import math

HERE=Path(__file__).resolve().parent
OUT=HERE/'independent_results.json'
OLD={
 '../1253/proof.md':'422c736a1c3d65671fcef1795cc889e44514516344510048be7c8240acc29161',
 '../1262/proof.md':'657a7d8ba0521bf105cf6ec95d1bf806571e482cf993c739164cadba16f7516',
}
# The exact frozen 1262 hash is checked below after fixing this transcription.
OLD['../1262/proof.md']='657a7d8ba0521bf105cbf6ec95d1bf806571e482cf993c739164cadba16f7516'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def fstr(x): return str(x)
def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def rootsq(x):
    a,b=math.isqrt(x.numerator),math.isqrt(x.denominator)
    assert a*a==x.numerator and b*b==x.denominator
    return F(a,b)
def ceilroot(x,k):
    lo,hi=0,1
    while hi**k*x.denominator<x.numerator: hi*=2
    while hi-lo>1:
        mid=(hi+lo)//2
        if mid**k*x.denominator>=x.numerator: hi=mid
        else: lo=mid
    assert hi**k>=x and (hi-1)**k<x
    return hi

def multiply(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,F(0))+x*y
    return {i:x for i,x in out.items() if x}
def power(a,n):
    out={0:F(1)}
    for _ in range(n): out=multiply(out,a)
    return out

def profile_bounds():
    # Independent closed multivariate chain formula (rather than recursively
    # differentiating a three-variable expression).
    eta_polys=[{0:F(1)}]
    eta_bounds=[F(1)]
    for n in range(1,7):
        old=eta_polys[-1]; new={}
        for k,c in old.items():
            new[k+2]=new.get(k+2,F(0))+c
            if k:new[k+1]=new.get(k+1,F(0))-k*c
        new={k:c for k,c in new.items() if c}
        eta_polys.append(new)
        eta_bounds.append(sum((abs(c)*F(k,2)**k for k,c in new.items()),F(0)))
    # At least one of 4-s,s-1 is >=3/2. e^(2/3)<2 proves denominator>1/2.
    x=F(2,3); terms=[x**n/F(math.factorial(n)) for n in range(11)]
    next_term=x**11/F(math.factorial(11))
    exp_upper=sum(terms,F(0))+next_term/(1-x/12)
    assert exp_upper<2
    rho=[F(1)]
    for n in range(1,7):
        rho.append(2*(eta_bounds[n]+sum((F(math.comb(n,k))*2*eta_bounds[k]*rho[n-k]
                                       for k in range(1,n+1)),F(0))))
    norms_by_order={r:F(0) for r in range(7)}
    count=0
    for a0 in range(7):
      for a1 in range(7-a0):
       for a2 in range(7-a0-a1):
        a=(a0,a1,a2); order=sum(a)
        sup=F(0)
        for beta in itertools.product(*(range(ai//2+1) for ai in a)):
            coeff=F(1)
            for ai,bi in zip(a,beta):
                coeff*=F(math.factorial(ai),math.factorial(bi)*math.factorial(ai-2*bi))
            coeff*=F(4)**(order-2*sum(beta))
            sup+=coeff*rho[order-sum(beta)]
        multi=F(math.factorial(order),math.prod(math.factorial(ai) for ai in a))
        norms_by_order[order]+=multi*64*sup*sup
        count+=1
    assert count==84
    S={j:sum((F(math.comb(j+2,r))*norms_by_order[r] for r in range(j+3)),F(0)) for j in (2,4)}
    # C_j^2 <= S/(8*pi) < S/24; x+1 >= sqrt(x).
    C={j:S[j]/24+1 for j in (2,4)}
    for j in C:assert C[j]*C[j]>S[j]/24
    L=max(10**4,ceilroot(625*C[2]/52,2),ceilroot(125*C[4]/281216,4))
    R=2*L
    d2=5*C[2]/(208*R**2);d4=C[4]/(281216*R**4)
    assert d2<=F(1,2000) and d4<=F(1,2000)
    assert C[2]/R**2<=F(13,625)
    return dict(eta=eta_bounds,rho=rho,S=S,C=C,L=L,R=R,D=d2+d4,
                order_norms=norms_by_order,profile_exp_upper=exp_upper)

def source_group(pr):
    # Rational on-shell momenta with exact square-root pair factors. This checks
    # all tensor entries and the proposed representation correction independently.
    base=[(F(5,4),(F(3,4),F(0),F(0))),
          (F(5),(F(2),F(2),F(4))),
          (F(45),(F(42),F(16),F(2)))]
    modes=[]
    for e,p in base:
      for shift in range(3):
       for sign in (-1,1): modes.append((e,tuple(sign*p[(j+shift)%3] for j in range(3))))
    checks=0;witness=None
    for m in (F(26),F(55,2),F(30),F(32)):
      for (ep0,p0),(eq0,q0) in itertools.product(modes,repeat=2):
        ep,eq=ep0*m,eq0*m;p=tuple(m*x for x in p0);q=tuple(m*x for x in q0)
        assert ep*ep==m*m+dot(p,p) and eq*eq==m*m+dot(q,q)
        s=rootsq(ep*eq);b=(ep+eq)/2;omega=ep-eq
        d=tuple(x-y for x,y in zip(p,q));d2=dot(d,d);V=2*m*m-d2
        mink=-ep*eq+dot(p,q); pv=(ep,)+p;qv=(eq,)+q
        tensor=[]
        for i in range(4):
            tensor.append([(pv[i]*qv[j]+qv[i]*pv[j]-(F(-1) if i==j==0 else F(1) if i==j else F(0))*(mink+m*m))/(2*s) for j in range(4)])
        cov=(ep-eq,)+(tuple(qi-pi for pi,qi in zip(p,q)))
        assert all(sum((cov[i]*tensor[i][j] for i in range(4)),F(0))==0 for j in range(4))
        Q=tensor[0][0]-sum((tensor[i][i] for i in range(1,4)),F(0))
        A=V/(2*b); B=V/(2*s)
        Cf=F(1)/(2*s)+V/(8*b*s*(b+s))
        assert Q==(V+omega*omega)/(2*s)
        assert Q-B==omega*omega/(2*s)
        assert B-A==V*omega*omega/(8*b*s*(b+s))
        assert Q-A==omega*omega*Cf
        assert abs(Q-A)<=5*d2/(8*m)+d2*d2/(16*m**3)
        assert abs(Cf)<=F(5)/(4*m)
        if witness is None and Q!=A: witness={'m':str(m),'A':str(A),'Q':str(Q),'difference':str(Q-A)}
        checks+=1
    assert checks==1296 and witness
    # Independent Lyapunov remainder using the *same* profile scale condition.
    lap=F(13,625);grad2=2*lap;elam2=F(101,100)**2
    first=32*elam2+lap/52+F(1,100)*grad2/26
    upper_second=64*elam2+grad2/26
    lower_magnitude=first*first/26
    assert first<33 and upper_second<66 and lower_magnitude<66
    # Positive inverse L_H^-1 preserves the two-sided bound. Integrating twice
    # gives 33 lambda^2, not an uncomputed formal O(lambda^2).
    return {
      'group':'same-source positive-frequency transport and finite remainder',
      'rational_shell_diagnostics':checks,
      'rational_shell_masses_are_diagnostics_not_actual_six_spectrum':True,
      'nonzero_representation_difference':witness,
      'fourier_convention':'(2*pi)^(-3)',
      'profile_derivatives_through':6,
      'eta_derivative_upper':[str(x) for x in pr['eta']],
      'rho_derivative_upper':[str(x) for x in pr['rho']],
      'S':{str(j):str(v) for j,v in pr['S'].items()},
      'Cbar':{str(j):str(v) for j,v in pr['C'].items()},
      'finite_integer_L':pr['L'],'R':pr['R'],'D_upper_exact':str(pr['D']),
      'delta_f_sup_upper':str(lap),'gradient_square_upper':str(grad2),
      'Hprime_upper':str(first),'Hsecond_positive_upper':str(upper_second),
      'Hsecond_negative_magnitude_upper':str(lower_magnitude),
      'finite_Taylor_coefficient':33,
      'static_family_not_time_switched_field':True,
      'fixed_L2_preparation_is_transported_Cauchy_data':True,
    }

def record_group(pr):
    c=power({-1:F(1,2),1:F(1,2)},4)
    N=sum((x*x for x in c.values()),F(0))
    first=sum((F(j*j)*x*x for j,x in c.items()),F(0))/N
    second=sum((F(j)**4*x*x for j,x in c.items()),F(0))/N
    assert N==F(35,128) and first==F(16,7) and second==F(512,35)
    lap_sq=3*second+6*first*first
    assert lap_sq==F(18432,245)<81
    for order in range(4):
        assert sum((F(j)**order*x*(-1 if (j//2)%2 else 1) for j,x in c.items()),F(0))==0
    assert F(64,35)*2*N==1  # normalized 1D position packet
    assert 3<4 and 9>4 # cube radius sqrt(3)L <2L; outside center10L
    # sin(3pi/8)>9/10 and cos(5pi/12)>1/4, verified by squared radicals.
    assert F(31,25)**2<2 and F(3,2)**2>2
    min_clock_slope=3*F(9,10)**5/F(4)
    assert min_clock_slope>F(3,8)
    assert F(1,4)+F(101,100)/2<F(5,6)
    r=F(1,10000);L=pr['L'];D=pr['D']
    packet_each=F(18,13)/L**2
    native=9*r/16-2*packet_each
    state_source_each=2*(r*D+33*r*r)
    source_contrast=native-2*state_source_each
    generic=9*r/16-F(36,13*10**8)-4*(r/F(1000)+33*r*r)
    assert generic==F(70853,1300000000)>r/2
    assert source_contrast>=generic and native>r/2
    # Both task probabilities equal under any one common free time by translation
    # covariance. A contrast gives a half-contrast lower bound for one task.
    assert source_contrast/2>F(1,40000)
    # Native internal projection, not a spatially local field instrument.
    energy_read_bound=F(5,2)*F(101,100)**2
    assert energy_read_bound==F(10201,4000)
    return {
      'group':'normal packets and common finite two-preparation menu',
      'cos4_norm':str(N),'first_derivative_ratio':str(first),'second_derivative_ratio':str(second),
      'three_dimensional_laplacian_square_coefficient':str(lap_sq),
      'clock_slope_lower':str(min_clock_slope),'initial_clock_age':'pi/4','T':'pi/2',
      'r':str(r),'lambda':'log(1+r), strictly smaller than r',
      'normal_state_error_each':str(packet_each),
      'source_comparison_error_each':str(state_source_each),
      'native_two_task_contrast_lower':str(native),
      'source_two_task_contrast_lower':str(source_contrast),
      'scale_independent_contrast_lower':str(generic),
      'scale_independent_display_only':float(generic),
      'any_common_free_time_at_least_one_task_error_lower':str(source_contrast/2),
      'native_nonselective_read_energy_change_upper':str(energy_read_bound),
      'two_preparations_not_single_joint_probability':True,
      'compact_positive_frequency_L2_not_compact_local_field':True,
      'no_actual_local_instrument_or_self_gravity_certification':True,
    }

def run():
    for rel,digest in OLD.items():assert sha((HERE/rel).resolve())==digest,rel
    pr=profile_bounds()
    return {'round':1264,'verdict':'PASS_INDEPENDENT_FINITE_STATIC_SOURCE_RECORD_BRIDGE',
      'independent_groups':2,'code_sha256':sha(Path(__file__).resolve()),
      'old_source_sha256':OLD,'groups':[source_group(pr),record_group(pr)],
      'scope':'static native one-particle same-source finite bridge; no dynamic metric, local field instrument, full interacting source, or GR claim'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as handle:json.dump(result,handle,ensure_ascii=False,indent=2);handle.write('\n')
        print('WRITE_ONCE_PASS_1264_INDEPENDENT',result['independent_groups'])
    else:
        assert json.loads(OUT.read_text(encoding='utf-8'))==result,'Saved independent result mismatch'
        print('READ_ONLY_MATCH_1264_INDEPENDENT',result['independent_groups'])
