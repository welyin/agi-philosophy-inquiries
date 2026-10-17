"""1263 independent two-group finite audit. Default read-only; first write exclusive.
No author implementation or results imported. Uniform proofs are textual;
finite exact algebra and rational bounds below diagnose them independently.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math

HERE=Path(__file__).resolve().parent
RESULT=HERE/'independent_results.json'
FROZEN={
    '../1253/proof.md':'422c736a1c3d65671fcef1795cc889e44514516344510048be7c8240acc29161',
    '../1261/parent.md':'1fd650bd75b2c345d599c602d6a13cf35f3263bfe5b14a6c6c75b6d2cc73776d',
    '../1262/proof.md':'657a7d8ba0521bf105cbf6ec95d1bf806571e482cf993c739164cadba16f7516',
}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def product(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,F(0))+x*y
    return {i:x for i,x in out.items() if x}

def power(a,n):
    out={0:F(1)}
    for _ in range(n):out=product(out,a)
    return out

def scatter():
    c=power({-1:F(1,2),1:F(1,2)},4)
    norm=sum(x*x for x in c.values())
    d4=sum(F(j)**8*x*x for j,x in c.items())
    assert norm==F(35,128) and d4==544
    ratio=d4/norm
    assert ratio==F(69632,35) and ratio/F(2)**8==F(272,35)<8
    # Zero-extension has no boundary distributions through the fourth derivative.
    for order in range(4):
        assert sum(F(j)**order*x*(-1 if (j//2)%2 else 1) for j,x in c.items())==0
    assert F(105)<F(2)**8
    clock={j:x/F(8) for j,x in power({-1:F(-1,2),0:F(1),1:F(-1,2)},3).items()}
    assert clock[0]==F(5,16)
    nonzero=sum(abs(x) for j,x in clock.items() if j)
    assert nonzero==F(11,16)
    b=F(10)**100;k=F(10)**40;sigma=F(10)**65
    assert sigma**20==b**13 and k**5==b**2
    h=128*k/b+16/sigma
    pbad=(96*sigma/b)**8+h**8+(4*h/k)**8
    # Simple independent rational envelopes replace the author's root algorithm.
    root2=F(1,10**132);root4=F(1,10**66)
    assert pbad<F(1,10**264) and root2**2>pbad and root4**4>pbad
    terms={
        'exact_massive_free':2*b*h*h/52,
        'exact_massless_free':2*b*h*h/k,
        'full_energy_spread':256*h,
        'full_linear_momentum':20*k*h,
        'all_TT_and_quadratic_remainder':F(40,26)*k*h*h,
        'three_spatial_Taylor':2*512*24**2*k*sigma**2/b**2,
        'bounded_completion_tail':2*(1024+64)*k*root2,
        'unbounded_affine_tail':2*64*24*k*sigma*root4/b,
    }
    eps=sum(terms.values(),F(0))
    assert eps<F(3191,10**20)<F(4,10**17)
    assert k-128*k/b-4/sigma>k/2
    support=2*32*2*4+4*2*2*4+F(8*2*4*4,26)
    assert support==F(7616,13)<1024
    assert 1024*k/b<26 and 2/F(3,4)**3<16
    assert F(3,4)>F(2,3) and F(5,4)<F(4,3)
    # Final photon-only p_Cx menu. Physical variance is 1/(4 sigma^2).
    A=8*sigma*sigma*k*k/b**2
    each_tail=F(1,4)/A/(26-13)**2
    assert each_tail==F(1,26**2)/A
    gap=F(5,16)-nonzero/(1+A)-2*each_tail-2*eps
    assert A==8*10**10 and gap>F(31,100)
    # Native angle menu: good denominator in [3k/4,5k/4].
    # Reference guaranteed-success threshold -65q0/4 and free possible-false
    # threshold -39q0/4 have the same worst-case gap 39q0/4.
    assert 26-F(65,4)==F(39,4) and F(13)*F(3,4)==F(39,4)
    angular_variance_tail=F(1,4)/A/F(39,4)**2
    assert angular_variance_tail==F(4,1521)/A
    photon_bad=(4*h/k)**8
    angle_each=angular_variance_tail+photon_bad
    angle_gap=F(5,16)-nonzero/(1+A)-2*angle_each-2*eps
    assert h/k<F(2,10**98) and angle_gap>F(31,100)
    return {
       'group':'complete finite scattering and photon-only bounded menu',
       'cos4_norm':str(norm),'fourth_derivative_norm':str(d4),
       'fourth_derivative_ratio':str(ratio),'centered_x8_sigma8_upper':str(ratio/2**8),
       'clock_fourier':{str(j):str(x) for j,x in sorted(clock.items())},
       'pbad_simple_upper':'1/10^264','full_support_coefficient':str(support),
       'eight_integrated_terms':{name:str(x) for name,x in terms.items()},
       'state_error_exact':str(eps),'state_error_display_only':float(eps),
       'photon_menu_each_tail_exact':str(each_tail),
       'photon_joint_gap_exact':str(gap),'photon_joint_gap_display_only':float(gap),
       'angle_each_wrong_box_exact':str(angle_each),
       'angle_joint_gap_exact':str(angle_gap),'angle_joint_gap_display_only':float(angle_gap),
       'angle_menu_native_only_not_actual_tetrad':True,
       'retiming_reason':'p_Cx commutes with exact free H at every time; same preparation/effect, no postselection',
       'whole_state_reference_uniform_by_mass_blocks':True,
       'all_original_mass_sectors_and_recoil_retained':True,
       'moments_are_centered_on_reference_path':True,
    },(b,k,sigma,eps)

def source(parameters):
    b,k,sigma,eps=parameters;L=b
    # Independent exact rational light cones, including arbitrarily low positive
    # energy, a high energy scale, antipodal saturation and zero-leg extension.
    directions=[(F(1),(F(1),F(0),F(0))),
                (F(5),(F(3),F(4),F(0))),
                (F(5),(F(-3),F(-4),F(0))),
                (F(13),(F(0),F(5),F(12)))]
    shells=[(F(0),(F(0),F(0),F(0)))]
    for scale in (F(1,256),F(1),F(10)**40):
        for e,p in directions:shells.append((scale*e,tuple(scale*x for x in p)))
    pairs=0;saturated=0;zero_legs=0
    for ep,p in shells:
        for eq,q in shells:
            n=ep*eq-sum(x*y for x,y in zip(p,q))
            d2=sum((x-y)**2 for x,y in zip(p,q))
            assert n>=0 and n<=2*ep*eq
            assert d2-(ep-eq)**2==2*n
            # Trace^2 <= |p-q|^2, proved without numerical square roots.
            assert n*n<=ep*eq*d2
            if ep*eq==0:
                assert n==0;zero_legs+=1
            elif n*n==ep*eq*d2 and d2>0:saturated+=1
            pairs+=1
    assert saturated>0 and zero_legs>0
    # Gaussian L1 kernel = 1; the bounded Weyl trace uses Schur, not a
    # general "bounded symbol => same operator bound" assertion.
    for mass in (F(26),F(27),F(31),F(32)):
        for p2 in (F(0),F(1,256),F(1),F(10)**80):
            assert mass**4<=mass**2*(mass**2+p2)
    assert 4*F(9,4*26)==F(9,26)
    fine=F(9,26)/L**2;null=F(2)/L
    B=32+fine+null
    transported=2*B*eps
    # Conditional A Wigner function is the specified positive Gaussian.
    reference=F(32)*3*sigma**2/(2*L**2)+(128*k/b)**2/52+F(3,208)/sigma**2
    assert reference==48*sigma**2/L**2+(128*k/b)**2/52+F(3,208)/sigma**2
    total=transported+reference+fine+null
    assert transported<F(2042,10**18) and total<F(21,10**16)
    assert total/26<F(1,10**16)
    # Conditional linear-Einstein transport only; pi<22/7, w nonnegative
    # normalized (or replace the central term by M*integral(w)).
    conditional_R=F(176,7)*total
    assert conditional_R<F(6,10**14)
    return {
       'group':'same-state bounded leading material trace insertion',
       'exact_null_trace_pairs':pairs,'antipodal_saturation_cases':saturated,
       'zero_leg_continuous_extension_cases':zero_legs,
       'massless_kernel_bound':'abs(trace K)<=abs(p-q)',
       'massless_Gaussian_operator_bound':'<2/L; no band, no k factor',
       'massive_Weyl_trace_bound':'<=32 by left/right kernel Schur',
       'trace_norm_bound_exact':str(B),
       'transported_trace_difference_exact':str(transported),
       'Gaussian_reference_to_mass_exact':str(reference),
       'actual_trace_to_original_mass_exact':str(total),
       'actual_trace_to_original_mass_display_only':float(total),
       'relative_to_mass26_exact':str(total/26),
       'conditional_linear_R_error_exact':str(conditional_R),
       'conditional_time_weight_nonnegative_normalized':True,
       'source_clock_simultaneous_instrument_proved':False,
       'full_interacting_AW_metric_variation_proved':False,
       'leading_material_trace_not_full_energy_or_geometry':True,
    }

def calculate():
    sources={}
    for path,expected in FROZEN.items():
        actual=sha(HERE/path);assert actual==expected,(path,actual);sources[path]=actual
    group1,parameters=scatter()
    return {'round':1263,'status':'PASS_INDEPENDENT_FINITE_SCATTERING_AND_LEADING_TRACE',
       'independent_groups':2,'independent_code_sha256':sha(Path(__file__)),
       'frozen_source_sha256':sources,'group1_scattering':group1,'group2_source':source(parameters),
       'scope':{'bounded_AW_effective_model':True,'unknown_reference_retained':True,
          'trace_not_dephasing':True,'no_new_TR_CO3_gap':True,
          'actual_instruments_not_certified':True,'EFT_remainder_not_certified':True,
          'complete_interacting_source_not_certified':True,'independent_geometry_not_certified':True}}

def compare(a,b,path='root'):
    if isinstance(b,dict):
        assert isinstance(a,dict) and a.keys()==b.keys(),path
        for key in b:compare(a[key],b[key],path+'.'+key)
    elif isinstance(b,list):
        assert isinstance(a,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'.'+str(i))
    elif isinstance(b,float):assert math.isclose(a,b,rel_tol=1e-13,abs_tol=0),(path,a,b)
    else:assert a==b,(path,a,b)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();data=calculate()
    if args.write_results:
        with RESULT.open('x',encoding='utf-8',newline='\n') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS_INDEPENDENT_1263: independent_groups=2')
    else:
        compare(json.loads(RESULT.read_text(encoding='utf-8')),data)
        print('READ_ONLY_MATCH_INDEPENDENT_1263: independent_groups=2')

if __name__=='__main__':main()
