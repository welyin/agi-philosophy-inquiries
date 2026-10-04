"""698: exact full-Haar upper bound for the original static centre entry.

Analytic inequalities dominate the actual source before integration. Complete
sphere moments and Weyl constant terms of the majorant are computed exactly.
This is not a certificate of negative offdiagonal kernel or physical RP.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'round698_drafts'))
import polynomial_majorant_probe as moment
import block_majorant_probe as blocks
TARGET=HERE/'joint_static_haar_majorant_results.json'


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def inequality_certificate():
    # beta^2 >= (1+c)/10+(1-c^2)/40; reduce square-root inequality to
    # a polynomial, c=2t-1, then exhibit positive Bernstein coefficients.
    f=mul(mul([15,-4,1],[15,-4,1]),[17,8]);f[0]-=3600
    transformed=[sum(f[j]*math.comb(j,k)*2**k*(-1)**(j-k) for j in range(k,len(f))) for k in range(len(f))]
    rem=transformed[:];quot=[0]*4
    for j in range(5,1,-1):quot[j-2]=-rem[j];rem[j]=0;rem[j-1]-=quot[j-2]
    assert rem==[0]*6
    bern=[sum(F(math.comb(i,j),math.comb(3,j))*quot[j] for j in range(i+1)) for i in range(4)]
    assert all(x>0 for x in bern)
    # Physical upper slope51/1000: g(1)=(9+4sqrt5)/20 <449/500.
    upper_sqrt=F(56,25);assert upper_sqrt**2>5
    assert (9+4*upper_sqrt)/20==F(449,500)
    return dict(beta_square_polynomial=f,unit_interval_factor_quotient=quot,
        positive_Bernstein_coefficients=[str(x) for x in bern],
        physical_endpoint_sqrt5_upper=str(upper_sqrt),physical_odd_upper='(949-51 cos(phi))/1000',
        color_diagonal_upper='(3696-8 sum_odd(r+r^-1)+sum_odd(r^2+r^-2))/10240',
        color_cross_square_upper='(66+25(z+z^-1))/5120',weak_diagonal='1/4',weak_cross_square='1/80')


def matrix_checks():
    old=blocks.old;rng=np.random.default_rng(69851)
    h,_=old.probe.torus(rng.uniform(0,2*np.pi,(4096,4)))
    v,physical=blocks.polar_batch(h);planes=blocks.plane_coefficients(v)
    odd=np.diag(old.b.prior.rep(np.eye(3),-np.eye(2),1)).real<0
    co=h[:,odd].real
    diagonal=.25+np.sum(.9-.1*co-(1-co*co)/40,axis=1)/64
    bounds=[];margins=[]
    for j,(aa,bb,cc,dd) in enumerate(planes):
        t0=old.b.internal.T[2*j];t1=old.b.internal.T[2*j+1]
        ii,jj=np.nonzero(t0)
        assert len(ii)==16 and np.array_equal(abs(t0),abs(t1))
        assert np.all(odd[ii]==odd[jj]) if j<3 else np.all(odd[ii]!=odd[jj])
        cosine=np.mean((h[:,ii]*h[:,jj]).real,axis=1)
        a=diagonal if j<3 else np.full(len(h),.25)
        c2=(33+25*cosine)/2560 if j<3 else np.full(len(h),1/80)
        margin=min(float((a-aa).min()),float((a-bb).min()),float((c2-cc*cc-dd*dd).min()))
        assert margin>-2e-14
        margins.append(margin);bounds.append((a,a,np.sqrt(c2),np.zeros(len(h))))
    physical_upper=np.prod((949-51*h[:,odd].real)/1000,axis=1)*np.prod((2+2*h[:,~odd].real)/16,axis=1)
    assert np.max(physical-physical_upper)<1e-18
    assert np.max(blocks.sphere_moment(planes)-blocks.sphere_moment(bounds))<1e-16
    polys=moment.polynomials();q=diagonal[:32]
    cs=[p[2][:32]**2 for p in bounds[:3]]
    direct=blocks.sphere_moment([tuple(x[:32] for x in p) for p in bounds])
    compact=moment.evaluate(q,cs,polys)
    err=float(np.max(abs(compact/direct-1)));assert err<1e-12
    return dict(samples=4096,coefficient_bound_margins=margins,
        compact_moment_relative_error=err,actual_five_plane_identity=blocks.check(),
        samples_check_implementation_not_the_inequality_proof=True)


def isprime(n):
    return n>=2 and (n==2 or (n%2 and all(n%d for d in range(3,math.isqrt(n)+1,2))))


def primes():
    n=2000000000-(2000000000-1)%528
    while True:
        if isprime(n):yield n
        n-=528


def root(n,p):
    for base in range(2,1000):
        r=pow(base,(p-1)//n,p)
        factors=(2,11) if n==44 else (2,3)
        if pow(r,n,p)==1 and all(pow(r,n//q,p)!=1 for q in factors):return r
    raise AssertionError('root not found')


def chamber(n):
    rows=[]
    for a,b,c in itertools.combinations(range(n),3):
        for d in range(n):
            e=(-a-b-c-d)%n
            if d<e:rows.append((a,b,c,d))
    return np.array(rows,dtype=np.int64)


def rational_mod(x,p):return (x.numerator%p)*pow(x.denominator,-1,p)%p


def modular_average(p,n,coords,polys):
    r=root(n,p);roots=np.array([pow(r,k,p) for k in range(n)],dtype=np.int64)
    count=len(coords);s1=np.zeros(count,dtype=np.int64);s2=s1.copy();physical=np.ones(count,dtype=np.int64)
    for w,odd in moment.weights():
        exponent=(coords@np.array(w,dtype=np.int64))%n
        val=(roots[exponent]+roots[(-exponent)%n])%p
        if odd:
            s1=(s1+val)%p;s2=(s2+roots[(2*exponent)%n]+roots[(-2*exponent)%n])%p
            factor=(1898-51*val)%p
        else:factor=(2+val)%p
        physical=physical*factor%p
    physical=physical*pow((16**8*2000**8)%p,-1,p)%p
    a=(3696-8*s1+s2)%p;a=a*pow(10240,-1,p)%p
    hs=[np.ones(count,dtype=np.int64)]+[np.zeros(count,dtype=np.int64) for _ in range(8)]
    for j in range(3):
        c=(66+25*(roots[coords[:,j]]+roots[(-coords[:,j])%n]))%p
        c=c*pow(5120,-1,p)%p
        for k in range(1,9):hs[k]=(hs[k]+c*hs[k-1])%p
    value=np.zeros(count,dtype=np.int64)
    for k,coeff in enumerate(polys):
        term=np.zeros(count,dtype=np.int64)
        for c in reversed(coeff):term=(term*a+rational_mod(c,p))%p
        value=(value+term*hs[k])%p
    exponents=[coords[:,j] for j in range(4)]+[(-coords.sum(axis=1))%n]
    density=np.ones(count,dtype=np.int64)
    for i,j in ((0,1),(0,2),(1,2),(3,4)):
        d=(exponents[i]-exponents[j])%n
        density=density*((2-roots[d]-roots[(-d)%n])%p)%p
    # Each nonzero orbit has12 points; that cancels the Weyl denominator12.
    total=int(np.sum(value*physical%p*density%p,dtype=np.int64))%p
    assert count*p<2**63 and 2*p*p<2**63
    return total*pow(n**4,-1,p)%p,r


def exact_average():
    polys=moment.polynomials();den=1
    for k,coeff in enumerate(polys):
        for i,c in enumerate(coeff):den=math.lcm(den,c.denominator*10240**i*5120**k)
    den*=12*16**8*2000**8
    # Positive compact polynomial has M <= its value at the component maxima.
    a=F(3,8);c=F(29,1280)
    mmax=sum(sum(x*a**i for i,x in enumerate(poly))*math.comb(k+2,2)*c**k for k,poly in enumerate(polys))
    bound=mmax*F(53,536870912);integer_bound=(bound*den).numerator//(bound*den).denominator+1
    coords=chamber(44);modulus=1;integer=0;rows=[];generator=primes()
    while modulus<=integer_bound:
        p=next(generator);avg,r=modular_average(p,44,coords,polys);residue=avg*(den%p)%p
        step=(residue-integer)%p*pow(modulus,-1,p)%p
        integer+=modulus*step;modulus*=p
        rows.append(dict(prime=p,root=r,scaled_integral_residue=residue))
    assert 0<=integer<=integer_bound
    exact=F(integer,den)
    # Independent oversampling order and a fresh prime detect degree/index errors.
    p=next(generator);other=chamber(48);avg,r=modular_average(p,48,other,polys)
    assert avg==rational_mod(exact,p)
    assert F(0)<exact<F(246,10**12)
    # Inherited free entry; remaining offdiagonal target sufficient for a
    # negative test is only a conditional implication, not a proved lower bound.
    b00=blocks.old.numerical.reference();lower_target=F(727,10**16);alpha=F(3,10000)
    conditional=b00-2*alpha*lower_target+alpha*alpha*exact
    assert conditional<0
    return dict(common_denominator=str(den),integer_upper_bound=str(integer_bound),
        recovered_integer=str(integer),CRT_modulus=str(modulus),prime_certificates=rows,
        degrees=[43,43,43,42],root_grid=[44]*4,full_grid_points=44**4,chamber_points=len(coords),
        check_grid=[48]*4,check_chamber_points=len(other),check_prime=p,check_root=r,
        exact_haar_majorant=str(exact),decimal_majorant=float(exact),
        prior_uniform_majorant=str(bound),positive_compact_polynomial_coefficients=81,
        moment_polynomials=[[str(x) for x in poly] for poly in polys],
        inherited_B00=str(b00),unproved_offdiagonal_sufficient_lower=str(lower_target),
        conditional_vector=[1,str(-alpha)],conditional_quadratic_upper=str(conditional),
        conditional_upper_decimal=float(conditional),offdiagonal_lower_bound_NOT_proved=True)


def run():
    inequalities=inequality_certificate();checks=matrix_checks();exact=exact_average()
    deps=('research_note_657.md','research_note_675.md','research_note_687.md','research_note_695.md',
        'research_note_697.md','joint_full_holonomy_reduction.py','joint_gauss_support_marginal_results.json',
        'round698_drafts/block_majorant_probe.py','round698_drafts/polynomial_majorant_probe.py',
        'round698_drafts/sphere_majorant_probe.py')
    return dict(date='2026-10-02',round=698,tests_run=2,failures=0,errors=0,
        analytic_inequality_certificate=inequalities,original_matrix_checks=checks,
        exact_integral_upper=exact,
        scope=dict(original_static_centre_B11_upper_certified=True,original_B11_exact_value_not_computed=True,
            whole_original_group_and_two_spheres=True,offdiagonal_integral_sign_certificate_not_complete=True,
            full_physical_RP_undecided=True,original_Hb_and_HF_not_replaced=True,
            inherited_space_interfaces_unchanged=True),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=698,upper=result['exact_integral_upper']['decimal_majorant'],
        primes=len(result['exact_integral_upper']['prime_certificates']),all_checks_passed=True)))
