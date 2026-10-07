"""973: audit the ACTUAL 964/969 thermal preparation and a monopole completion.
No full stochastic gravity, coherent quantum geometry or cosmological data fit.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'common_source_distribution_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    c=read(STAGE/'964/material_cosmology_results.json')
    old=read(STAGE/'969/local_field_cosmology_results.json')
    k=c['material']['J'];m=c['cosmology']['photon_modes'];nb=c['cosmology']['mean_occupation']
    M=old['local_sources']['complete_mass'];mu=old['shared_background_parameters']['mu']
    T=old['shared_background_parameters']['T'];meanR=m*k*nb;r=nb/(1+nb)
    varR=m*k*k*nb*(1+nb);Nmax=20000
    # Exact negative-binomial occupation distribution. No Monte Carlo.
    pr=np.empty(Nmax+1);pr[0]=(1-r)**m
    for n in range(Nmax):pr[n+1]=pr[n]*r*(n+m)/(n+1)
    ns=np.arange(Nmax+1);R=k*ns
    # Chernoff bound from E exp(s N)=((1-r)/(1-r exp(s)))^m.
    cutoff=Nmax+1;s=math.log(cutoff/(r*(cutoff+m)))
    logtail=-s*cutoff+m*(math.log1p(-r)-math.log1p(-r*math.exp(s)))
    tail=math.exp(logtail)
    assert tail<1e-20 and abs(sum(pr)-1)<2e-12
    mean_check=float(pr@R);var_check=float(pr@((R-meanR)**2))
    assert abs(mean_check-meanR)<1e-9 and abs(var_check-varR)<1e-8
    def primitive(a,x):return 2*(M*a-2*x)*np.sqrt(M*a+x)/(3*M*M)
    def time(a,x):return math.sqrt(mu)*(primitive(a,x)-primitive(1.,x))
    def scale(x):
        x=np.asarray(x);lo=np.ones_like(x,dtype=float);hi=np.full_like(x,8.,dtype=float)
        assert np.all(time(hi,x)>T)
        for _ in range(70):
            mid=(lo+hi)/2;low=time(mid,x)<T
            lo=np.where(low,mid,lo);hi=np.where(low,hi,mid)
        return (lo+hi)/2
    abar=float(scale(meanR));aa=scale(R);zz=1/aa;zbar=1/abar
    assert abs(abar-old['background']['scale_at_T'])<1e-12
    assert np.all(np.diff(aa)>0)
    weights=pr/sum(pr) # only floating normalization; physical tail below bound.
    mean_a=float(weights@aa);std_a=math.sqrt(float(weights@((aa-mean_a)**2)))
    mean_z=float(weights@zz);var_z=float(weights@((zz-mean_z)**2))
    # A finite Lipschitz report, 0 <= f <= 1, L=100 in scale-factor units.
    width=.02
    f=np.minimum(1.,((aa-abar)/width)**2)
    contrast=float(weights@f)
    # Initial Hamiltonian constraint per branch: p_a=-2 sqrt(mu) sqrt(M+R).
    root=np.sqrt(M+R);rootmean=float(weights@root)
    # Exponentially small tail in the unbounded root moments is independently
    # negligible; Jensen identity is analytical and the finite sums diagnostic.
    constraint_at_means=M+meanR-rootmean**2
    equivalent_variance=float(weights@((root-rootmean)**2))
    assert abs(constraint_at_means-equivalent_variance)<2e-10
    cv_H_squared=math.sqrt(varR)/(M+meanR)
    # 350's mixture requirement applied to a concrete constrained two-shell family.
    shells=np.array([2000.,3200.]);ashell=scale(k*shells)
    amid=float(scale(k*float(np.mean(shells))))
    shell_reports=np.minimum(1.,((ashell-amid)/width)**2)
    mix=float(np.mean(shell_reports));report_after_mean=0.
    assert np.all(abs(ashell-amid)>width) and mix==1.
    # Rigorous finite-probability lower certificate: tails of the original Gibbs law
    # lying at least width away already give f=1. Verify endpoints with independent
    # 64-point Gauss-Legendre integration of the defining positive time integral.
    nlo=int(ns[aa<=abar-width][-1]);nhi=int(ns[aa>=abar+width][0])
    nodes,ww=np.polynomial.legendre.leggauss(64)
    def integ(a,x):
        y=1+(a-1)*(nodes+1)/2
        return math.sqrt(mu)*(a-1)/2*float(ww@(y/np.sqrt(M*y+x)))
    tlo=integ(abar-width,k*nlo);thi=integ(abar+width,k*nhi)
    # positive monotonic source: t(a,R) decreases with R.
    assert tlo>T and thi<T
    saturated=float(sum(pr[:nlo+1])+sum(pr[nhi:]))
    # Margin on all floating finite sums/endpoints is deliberately much larger
    # than observed recurrence error; Decimal independent sums supplied separately.
    conservative_report_lower=saturated-1e-9
    assert conservative_report_lower>.5
    # Independent Decimal recurrence certifies the probability lower witness.
    from decimal import Decimal,localcontext
    with localcontext() as ctx:
        ctx.prec=65
        # Use the exact stored nbar defining this diagnostic.
        nd=Decimal.from_float(nb);rd=nd/(1+nd);p=(1-rd)**m;lowprob=Decimal(0)
        for n in range(Nmax+1):
            if n<=nlo or n>=nhi:lowprob+=p
            p*=rd*Decimal(n+m)/Decimal(n+1)
        decimal_probability=str(lowprob)
        md=Decimal.from_float(M);mud=Decimal.from_float(mu);td=Decimal.from_float(T)
        kd=Decimal.from_float(k);rdmean=Decimal.from_float(meanR)
        def dt(ad,xx):
            def f(z):return 2*(md*z-2*xx)*(md*z+xx).sqrt()/(3*md*md)
            return mud.sqrt()*(f(ad)-f(Decimal(1)))
        al,ah=Decimal(1),Decimal(4)
        for _ in range(200):
            mid=(al+ah)/2
            if dt(mid,rdmean)<td:al=mid
            else:ah=mid
        lower_margin=dt(al-Decimal(1)/50,kd*nlo)-td
        upper_margin=td-dt(ah+Decimal(1)/50,kd*nhi)
        assert lower_margin>6 and upper_margin>1
        assert ah-al<Decimal('1e-50')
        certified_endpoints=dict(mean_scale_lower=str(al),mean_scale_upper=str(ah),
            low_shell_time_margin=str(lower_margin),high_shell_time_margin=str(upper_margin),
            decimal_arithmetic_guard='1e-35')
    assert abs(float(decimal_probability)-saturated)<1e-11
    # A supplied branch-width contract yields a sufficient observable bound,
    # but this actual preparation has not been made narrow by raising nb.
    high_occupation_limit=1/math.sqrt(m)
    source_cv=math.sqrt(varR)/meanR
    files=[STAGE/'964/material_cosmology_results.json',STAGE/'969/local_field_cosmology_results.json',
           STAGE.parent/'archive_342_369/research_note_350.md',STAGE/'research_note_871.md',
           HERE/'drafts/common_realization_decision.md']
    return dict(round=973,all_scientific_checks_passed=True,
        frozen_inputs=dict(k=k,modes=m,mean_occupation=nb,ratio=r,M=M,mu=mu,T=T,mean_R=meanR),
        source_statistics=dict(exact_variance=varR,numerical_mean=mean_check,numerical_variance=var_check,
            coefficient_of_variation=source_cv,high_occupation_cv_limit=high_occupation_limit,
            initial_H_squared_relative_std=cv_H_squared),
        summation=dict(max_N=Nmax,normalization=float(sum(pr)),chernoff_tail=tail,
            lower_event_decimal=decimal_probability,decimal_endpoint_check=certified_endpoints),
        monopole_branch_diagnostic=dict(mean_source_scale=abar,mean_branch_scale=mean_a,
            branch_scale_std=std_a,mean_source_redshift=zbar,mean_branch_redshift=mean_z,
            mean_redshift_difference=mean_z-zbar,
            redshift_variance=var_z,report_width=width,report_lipschitz=2/width,
            report_expectation=contrast,report_under_mean_source=0.,
            saturated_tail_probability=saturated,report_certified_lower=conservative_report_lower,
            low_N_max=nlo,high_N_min=nhi,endpoint_time_margins=[tlo-T,T-thi],
            constraint_residual_at_mean_initial_momentum=constraint_at_means,
            root_energy_variance=equivalent_variance),
        mixture_witness=dict(occupation_shells=shells.tolist(),scales=ashell.tolist(),
            scale_of_mean_source=amid,mixture_report=mix,mean_source_report=report_after_mean,
            initial_geometry_correlated_with_energy=True),
        scope=dict(actual_964_Gibbs_energy_statistics=True,
            homogeneous_energy_resolved_completion_is_an_added_diagnostic=True,
            full_stress_noise_or_Einstein_Langevin_solved=False,
            original_969_mean_field_result_refuted=False,
            all_unknown_quantum_inputs_supported=False,
            actual_geometric_detector_constructed=False,
            full_common_positive_quantum_model=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=['https://arxiv.org/abs/quant-ph/0102125','https://arxiv.org/abs/gr-qc/0307032'])
def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'/{i}')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-9,abs_tol=2e-10),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(OUT))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','references')},ensure_ascii=False,indent=2))
