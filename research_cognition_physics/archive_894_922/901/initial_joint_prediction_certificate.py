"""901: use the continuum defect certificate for common initial geometry/reference.
The exact enclosure is not a time-propagation or full quantum-source certificate.
"""
import argparse,json,random
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import wiener_certificate as w
HERE=Path(__file__).resolve().parent;TARGET=HERE/'initial_joint_prediction_certificate_results.json'
def point(p,derivative=None):
    value=0;imag=0
    for k,(r,i) in p.c.items():
        if derivative is not None:
            f=k[derivative]*(2 if derivative==2 else 1);r,i=-f*i,f*r
        phase=(k[1]+k[2])%4
        r,i=((r,i),(-i,r),(-r,-i),(i,-r))[phase]
        value+=r;imag+=i
    assert imag==0
    return Q(value,w.S)
def interval_data(x):return {'lower':w.bound(x.a),'upper':w.bound(x.b)}
def convolution_checks():
    rng=random.Random(901)
    for trial in range(150):
        big=100 if trial<100 else 2**95
        a={tuple(rng.randint(-2,2) for _ in range(3)):(rng.randint(-big,big),rng.randint(-big,big)) for _ in range(15)}
        b={tuple(rng.randint(-2,2) for _ in range(3)):(rng.randint(-big,big),rng.randint(-big,big)) for _ in range(13)}
        ref={}
        for k,(r,i) in a.items():
            for m,(s,t) in b.items():
                n=tuple(x+y for x,y in zip(k,m));x,y=ref.get(n,(0,0));ref[n]=(x+r*s-i*t,y+r*t+i*s)
        assert w.integer_convolve(a,b)=={k:v for k,v in ref.items() if v!=(0,0)}
    return 150

def source_checks(src):
    import sys
    sys.path.insert(0,str(w.ROOT/'scripts'))
    from research_layout import Layout,ResearchRuntime
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        q,_,At,_,_=old.completed(24,1.)
        indices=np.random.default_rng(901).integers(0,24,size=(60,3));points=np.array([q['grid'][tuple(i)] for i in indices]);actual={k:[] for k in ('A','B','C','Y')}
        for j in indices:
            k=tuple(j);x=q['grid'][k][0]
            actual['A'].append(np.sum(At[k]**2)+q['pKp'][k]);actual['B'].append(q['B'][k]+.0004*np.cos(x)**2)
            actual['C'].append(q['C'][k]-.0004*np.sin(x)**2);actual['Y'].append(q['Y'][k])
        differences={k:float(np.max(np.abs(w.evaluate(src[k],points)-actual[k]))) for k in actual}
        assert max(differences.values())<1e-12
        md=[float(np.max(np.abs(w.evaluate(src['M'][j],points)-np.array([q['mom'][tuple(i)][j] for i in indices])))) for j in range(3)]
        assert max(md)<1e-13
    return dict(original572_753_source_max_differences=differences,momentum_differences=md,sample_count=60,
        diagnostic_only_not_basis_of_enclosure=True)

def run(reproduce=True):
    raw=json.loads(w.TARGET.read_text('utf-8'))
    if reproduce:
        new=w.run()
        assert {k:v for k,v in new.items() if k!='elapsed_seconds'}=={k:v for k,v in raw.items() if k!='elapsed_seconds'}
    r=Q(raw['continuum_residual_upper']['exact']);first_e0=Q(raw['C0_error_upper']['exact'])
    lowhat=Q(raw['witness_global_lower']['exact']);highhat=Q(raw['witness_global_upper']['exact'])
    low=lowhat-first_e0;high=highhat+first_e0
    base=json.loads((HERE.parent/'900/continuous_reference_bounds_results.json').read_text('utf-8'))
    coeff=base['coefficients'];c=lambda key:Q(coeff[key]['exact'])
    amin=5*c('C_min')*low**4-c('B_max')+6*c('Y_min')*high**-4
    amax=5*c('C_max')*high**4+7*c('A_max')*low**-8+6*c('Y_max')*low**-4
    assert amin>0
    c0=1/amin;c1=Q(10,8)*(1+amax/amin);e0=c0*r;e1=c1*r
    assert e0<first_e0
    grads=[Q(row['exact']) for row in raw['each_candidate_gradient_W1']]
    Merror=[8*Q('.008868915')*(low**-9*e1+9*low**-10*g*e0) for g in grads]
    # Full initial conformal metric and all path-length/distance comparisons.
    metric_error=4*high**3*e0
    theta=e0/lowhat;length_lower=(1-theta)**2;length_upper=(1+theta)**2
    witness=json.loads(w.WITNESS.read_text('utf-8'));psi=w.P({tuple(row[:3]):tuple(row[3:]) for row in witness['coefficients']})
    ps=point(psi);dp=[point(psi,j) for j in range(3)]
    mhat=[-8*Q('.008868915')*ps**-9*d for d in dp]
    m=[w.I(mhat[j]-Merror[j],mhat[j]+Merror[j]) for j in range(3)]
    assert m[2].a>0
    params=w.parameters();hs,ss=params['hs'],params['ss'];br=Q('.06')*ss;H=Q('1.05')*hs
    FF=2-(H*H+ss*ss)/6;Ch=-FF*H*ss/12*(Q('.025')-Q('.0056')/br)
    psI=w.I(ps-e0,ps+e0);vh=Ch*(psI**-6)
    detY=br*Q('.02')*m[2];det4=vh*detY
    # Inverse norm is coordinate-unit dependent, retained only as error budget.
    maxabs=lambda i:max(abs(i.a),abs(i.b))
    inv_upper=max(Q(50),1/br.a,(1+maxabs(m[0])/Q('.02')+maxabs(m[1])/br.a)/m[2].a)
    assert det4.a>0 and length_lower>0
    src=w.source();diag=source_checks(src);checks=convolution_checks()
    return dict(round=901,date='2026-10-06',formal_reports=901,cumulative_numbered_groups=3686,fresh_numbered_groups=1,
        all_checks_passed=True,argument_scope='Same859 exact continuous Cauchy solution versus explicit dyadic approximation: whole-torus residual, metric/path-length bounds and certified material Jacobian at original point. Not quantum or time-evolution closure.',
        exact_signed_convolution_crosschecks=checks,source_crosschecks=diag,
        continuous_residual=w.bound(r),first_global_C0_bound=w.bound(first_e0),bootstrap_interval=interval_data(w.I(low,high)),
        refined_reaction_min=w.bound(amin),refined_reaction_max=w.bound(amax),C0_error=w.bound(e0),each_C1_error=w.bound(e1),
        initial_metric_component_error=w.bound(metric_error),all_paths_and_distances_relative_interval=interval_data(w.I(length_lower,length_upper)),
        exact_candidate_point_value=w.bound(ps),each_candidate_point_derivative=[w.bound(a) for a in dp],
        each_actual_point_magnetic_derivative_interval=[interval_data(a) for a in m],actual_point_clock_velocity_interval=interval_data(vh),
        actual_spatial_Jacobian_interval=interval_data(detY),actual_spacetime_Jacobian_interval=interval_data(det4),
        coordinate_dependent_spatial_inverse_infinity_upper=w.bound(inv_upper),
        real_continuum_defect_not_collocation=True,reference_point_is_original=True,
        all_model_parameters_original=True,full_spacetime_patch_or_quantum_source_budget=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--reuse-residual',action='store_true');args=p.parse_args()
    r=run(not args.reuse_residual)
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    # The file preserves all exact fractions; output is deliberately compact.
    print(json.dumps({k:v for k,v in r.items() if k in ('round','all_checks_passed','source_crosschecks','full_goal_completed')},ensure_ascii=False,indent=2))
    for key in ('continuous_residual','C0_error','each_C1_error','initial_metric_component_error','coordinate_dependent_spatial_inverse_infinity_upper'):print(key,r[key]['value'])
    for key in ('all_paths_and_distances_relative_interval','actual_spatial_Jacobian_interval','actual_spacetime_Jacobian_interval'):print(key,r[key]['lower']['value'],r[key]['upper']['value'])
