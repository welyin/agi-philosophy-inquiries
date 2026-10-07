"""P981 free-conformal vacuum: a band-limited induced Einstein-tensor report.

Exact rational filter integral + independent tensor and quadrature diagnostics.
This is not a full SM, detector, or quantum-gravity truncation certificate.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TARGET=HERE/'smeared_source_geometry_results.json'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def frac(x): return {'exact':str(x),'decimal':float(x)}

def exact_integral():
    # int_0^1 dw (1-w^2)^4 int_0^w dk k^6 (1-k^2)^4
    return sum((F((-1)**(i+j)*math.comb(4,i)*math.comb(4,j),
                  (7+2*j)*(8+2*(i+j))) for i in range(5) for j in range(5)),F(0))

def rational_certificate(c,I):
    def atan_bounds(x,n):
        partial=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
        other=partial+(-1)**n*x**(2*n+1)/F(2*n+1)
        return min(partial,other),max(partial,other)
    a,b=atan_bounds(F(1,5),12);d,e=atan_bounds(F(1,239),4)
    pi_lo,pi_hi=16*a-4*e,16*b-4*d
    assert F(314159,100000)<pi_lo<pi_hi<F(1571,500)
    # Positive Taylor remainders certify log(10000)<9.211 and 1/e<.368.
    x=F(9211,1000)
    assert sum((x**j/F(math.factorial(j)) for j in range(61)),F(0))>10000
    assert sum((F(1,math.factorial(j)) for j in range(8)),F(0))>F(125,46)
    eta=F(1,10000)*c/(8*F(314159,100000)**2)*(x+F(46,125)+F(1571,500))
    V=c*I/(24*F(314159,100000)**4)*F(1,100000000)/(1-eta)**2
    assert eta<F(38,1000000)
    assert V<F(6118,10**19)
    assert V/F(1,10**10)<F(6118,10**9)
    return dict(method='Machin alternating bounds and positive exponential Taylor sums',
        response_shift_upper=frac(eta),variance_upper=frac(V),
        Chebyshev_upper=frac(V*10**10),
        stated_probability_upper='6118/1000000000',all_exact_inequalities_passed=True)

def quadrature(c, ratio, n):
    nodes,weights=np.polynomial.legendre.leggauss(n)
    w=(nodes+1)/2; weights=weights/2
    # k=w*v; triangular timelike support, both frequency signs integrated later.
    v=w.copy(); k=w[:,None]*v[None,:]
    s=w[:,None]**2-k**2
    base=(weights[:,None]*weights[None,:]*w[:,None]
          *(1-w[:,None]**2)**4*k**6*(1-k**2)**4)
    q=ratio**2*s
    # mu=M is fixed; zero finite spin-2 local coefficient at that scale is input.
    response=1-c/(8*math.pi**2)*q*(np.log(q)-1j*math.pi)
    return dict(I=float(base.sum()),
        response_I=float((base/abs(response)**2).sum()),
        response_difference_I=float((base*abs(1/response-1)**2).sum()),
        sampled_max_response_shift=float(abs(response-1).max()))

def tensor_checks(c):
    eta=np.diag([-1.,1.,1.,1.])
    rng=np.random.default_rng(984)
    ward=trace=energy=0.
    mineig=0.
    for _ in range(24):
        k=rng.normal(size=3); k*=rng.uniform(.05,.8)/np.linalg.norm(k)
        omega=1.+rng.uniform(.1,.7)
        p=np.r_[omega,k]; pl=eta@p; s=omega**2-k@k
        pi=eta+np.outer(pl,pl)/s
        P=(np.einsum('ac,bd->abcd',pi,pi)+np.einsum('ad,bc->abcd',pi,pi))/2
        P-=np.einsum('ab,cd->abcd',pi,pi)/3
        N=c*s*s/(8*math.pi)*P
        ward=max(ward,float(abs(np.einsum('a,abcd->bcd',p,N)).max()))
        trace=max(trace,float(abs(np.einsum('ab,abcd->cd',eta,N)).max()))
        energy=max(energy,abs(float(N[0,0,0,0])-c*float(k@k)**2/(12*math.pi)))
        # Any real tensor contraction is a positive covariance quadratic form.
        flat=N.reshape(16,16)
        mineig=min(mineig,float(np.linalg.eigvalsh(flat).min()))
    # Scalar normalization: canonical identical-particle pair phase space.
    # <(e_ij n_i n_j)^2>=2/15 for traceless e:e=1.
    scalar_rho=F(1,16)*F(2,15)*F(1,4) # coefficient of omega^4/pi
    # Matrix element squared=4 k^4 (e nn)^2; phase space/(2!)=1/(16 pi), k=w/2.
    assert scalar_rho==F(1,480)
    scalar_N=scalar_rho/2
    assert scalar_N==F(1,960)==F(1,120)/8
    # 631 canonical Dirac cut: rho=omega^4/(80 pi); Weyl has half weight.
    assert F(1,80)/2==F(1,20)/8
    return dict(sampled_timelike_momenta=24,ward_max=ward,trace_max=trace,
        energy_component_max=energy,minimum_covariance_eigenvalue=mineig,
        scalar_noise_coefficient_times_pi=str(scalar_N),
        Dirac_noise_coefficient_times_pi=str(F(1,160)))

def run():
    old=json.loads((HERE.parent/'983/sm_curvature_radiation_results.json').read_text('utf-8'))
    m=old['inherited_matter']
    counts=(m['real_scalars'],m['Weyl_components'],m['gauge_vectors'])
    assert counts==(4,45,12)
    c=counts[0]*F(1,120)+counts[1]*F(1,40)+counts[2]*F(1,10)
    assert c==F(m['c']['exact'])==F(283,120)
    I=exact_integral(); assert I>0
    C=float(c*I/24)/math.pi**4
    checks=tensor_checks(float(c))
    assert max(checks[k] for k in ('ward_max','trace_max','energy_component_max'))<1e-13
    assert checks['minimum_covariance_eigenvalue']>-1e-13
    rows=[]
    for ratio in (.1,.03,.01,.003):
        # Uniform proof uses x|log x|<=1/e on 0<=x<=1, not sampling.
        eps=ratio**2*float(c)/(8*math.pi**2)*(abs(math.log(ratio**2))+1/math.e+math.pi)
        assert eps<1
        q48=quadrature(float(c),ratio,48);q80=quadrature(float(c),ratio,80)
        assert abs(q48['I']-float(I))<2e-18
        assert abs(q80['I']-float(I))<2e-18
        # Endpoint q*log(q) slows quadrature; this is a diagnostic, not the bound.
        assert abs(q48['response_I']-q80['response_I'])<1e-14
        V0=C*ratio**4
        VF=float(c)/(24*math.pi**4)*ratio**4*q80['response_I']
        Vdiff=float(c)/(24*math.pi**4)*ratio**4*q80['response_difference_I']
        upper=V0/(1-eps)**2;lower=V0/(1+eps)**2
        assert lower<=VF<=upper
        assert Vdiff<=(eps/(1-eps))**2*V0
        assert q80['sampled_max_response_shift']<=eps
        rows.append(dict(Lambda_over_M=ratio,response_shift_upper=eps,
            leading_variance=V0,retained_response_variance=VF,
            retained_variance_lower=lower,retained_variance_upper=upper,
            rms_upper=math.sqrt(upper),
            same_source_response_difference_rms=math.sqrt(Vdiff),
            same_source_response_difference_rms_upper=eps/(1-eps)*math.sqrt(V0),
            Chebyshev_error_at_tolerance_1e_minus_5=min(1.,upper/1e-10),
            quadrature_48_80_difference=abs(q48['response_I']-q80['response_I'])))
    # Conservation cannot be preserved by substituting independent component noises.
    k=.5;omega=1.2;s=omega**2-k*k
    eta=np.diag([-1.,1.,1.,1.]);p=np.array([omega,k,0.,0.]);pl=eta@p
    pi=eta+np.outer(pl,pl)/s
    P=(np.einsum('ac,bd->abcd',pi,pi)+np.einsum('ad,bc->abcd',pi,pi))/2
    P-=np.einsum('ab,cd->abcd',pi,pi)/3
    N=float(c)*s*s/(8*math.pi)*P
    # Variance of omega*T00+k*T10, with and without their covariance.
    shared=float(omega**2*N[0,0,0,0]+k*k*N[1,0,1,0]+2*omega*k*N[0,0,1,0])
    wrong=float(omega**2*N[0,0,0,0]+k*k*N[1,0,1,0])
    assert abs(shared)<1e-15 and wrong>.001
    sourcepaths=[HERE.parent/'983/sm_curvature_radiation_results.json',
        HERE.parent/'981/drafts/common_parent_contract_v1.md',
        HERE.parent.parent/'archive_629_652/research_note_631.md',
        HERE.parent/'research_note_899.md',HERE.parent/'research_note_940.md']
    return dict(round=984,all_scientific_checks_passed=True,
        inherited_matter=dict(real_scalars=4,Weyl_components=45,gauge_vectors=12,c=frac(c)),
        analytic=dict(filter_integral=frac(I),
            leading_variance_coefficient=C,
            leading_variance_coefficient_times_pi4=frac(c*I/24),
            variance_scaling='Var(G00[f]/Lambda^2)=C*(Lambda/M)^4',
            no_time_smearing='positive omega-independent integrand at fixed nonzero k; divergent integral'),
        tensor_checks=checks,benchmarks=rows,rational_certificate=rational_certificate(c,I),
        negative_control=dict(shared_conservation_variance=shared,
            independent_component_conservation_variance=wrong),
        scope=dict(free_conformal_vacuum=True,actual_SM_thermal_state=False,
            matter_induced_linear_Einstein_report_only=True,
            finite_band_retained_response_bound=True,
            Gaussian_stress_process_claimed=False,
            finite_spacetime_detector_certified=False,
            full_quantum_geometry_error_certified=False,
            full_physical_remainder_certified=False,UV_completion_assumed=False,
            minimum_physical_scale_assumed=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sourcepaths})

def compare(new,old):
    if isinstance(new,dict):
        assert new.keys()==old.keys()
        for k in new:compare(new[k],old[k])
    elif isinstance(new,list):
        assert len(new)==len(old)
        for a,b in zip(new,old):compare(a,b)
    elif isinstance(new,float):assert math.isclose(new,old,rel_tol=1e-10,abs_tol=1e-14),(new,old)
    else:assert new==old,(new,old)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(result,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(result,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:result[k] for k in ('round','analytic','benchmarks','negative_control','scope')},
                     ensure_ascii=False,indent=2))
