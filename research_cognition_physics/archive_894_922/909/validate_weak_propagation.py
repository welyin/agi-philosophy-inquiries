"""909 reproduce one actual source response; distinguish checks from physical acceptance."""
from pathlib import Path
import hashlib,json,argparse
import numpy as np
import weak_source_propagation as w
HERE=Path(__file__).resolve().parent;TARGET=HERE/'weak_propagation_validation_results.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b):
    if isinstance(a,dict):
        for k in a:close(a[k],b[k])
    elif isinstance(a,(int,float)) and not isinstance(a,bool):assert np.isclose(a,b,rtol=5e-9,atol=2e-12),(a,b)
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):close(x,y)
    else:assert a==b,(a,b)

def run(reproduce=True):
    saved=json.loads((HERE/'weak_propagation_checks_results.json').read_text('utf-8'))
    charts=json.loads((HERE/'observer_chart_check_results.json').read_text('utf-8'))
    quad=json.loads((HERE/'observer_quadrature_check_results.json').read_text('utf-8'))
    assert len(saved['rows'])==6 and len(quad['rows'])==4 and len(charts['rows'])==2
    for p in saved['projection_checks']:
        assert p['B_pair_error']<2e-12 and p['S_pair_error']<2e-17 and p['independent_trig_value_error']<1e-12
    if reproduce:
        bg=w.mb.Background(16);src=w.loop.LoopSource(bg,anchor_order=2,segments=8);src.kernels()
        check=w.projection_check(src,16);close(check,{k:v for k,v in saved['projection_checks'][0].items() if k!='N'})
        with w.ev.ResearchRuntime(w.ev.Layout()).installed():
            import joint_reference_constraint_strata as old
            model=w.ev.Model(16,old);analytic=w.tangent.AnalyticModel(16,old);y0,_=model.initial()
            nodes,items,_,_=w.build_forcing(src,16,channel=2,time_order=2)
            actual,_,_=w.propagate(model,analytic,y0,nodes,items)
        close(actual,{k:v for k,v in saved['rows'][0].items() if k not in ('metadata','source_free_constraint_fourier_coefficients')})
    def row(N,order,steps=1,omit=False):
        return next(a for a in saved['rows'] if a['metadata']['N']==N and a['metadata']['time_order']==order
                    and a['substeps_per_interval']==steps and a['omitted_temporal_potential']==omit)
    a=row(16,4);b=row(24,4);c=row(16,2);d=row(16,4,2);bad=row(16,4,omit=True)
    delta=lambda x,y:{k:abs(x['observations'][k]-y['observations'][k]) for k in ('relational_probe','relational_magnetic')}
    change=dict(time_interpolation2_to4=delta(a,c),time_step_halving=delta(a,d),space16_to24=delta(a,b),
                omit_temporal_potential=delta(a,bad),
                observer_quadrature_range={k:float(np.ptp([q[k] for q in quad['rows']])) for k in ('relational_probe','relational_magnetic')})
    assert change['time_interpolation2_to4']['relational_probe']<1e-8
    assert change['time_step_halving']['relational_probe']<1e-8
    assert change['space16_to24']['relational_probe']>.1
    assert min(q['relational_probe'] for q in quad['rows'])<0<max(q['relational_probe'] for q in quad['rows'])
    assert b['source_free_constraint_fourier_coefficients']['H']['1']<a['source_free_constraint_fourier_coefficients']['H']['1']/100
    assert change['omit_temporal_potential']['relational_probe']>1e-6
    files=[Path(__file__),HERE/'weak_source_propagation.py',HERE/'weak_propagation_checks.py',HERE/'weak_propagation_checks_results.json',
           HERE/'observer_chart_check.py',HERE/'observer_chart_check_results.json',HERE/'observer_quadrature_check.py',HERE/'observer_quadrature_check_results.json']
    return dict(round=909,date='2026-10-06',cumulative_numbered_groups=3694,
        transport_and_diagnostic_checks_passed=True,one_actual_response_reproduced=reproduce,
        all_six_response_cases_previously_executed=True,all_saved_cases_rerun=False,
        source_file_sha256={str(p.relative_to(HERE)):sha(p) for p in files},differences=change,
        argument_scope='Actual908 adjoint-loop first-jet source mapped to finite Fourier/time Duhamel impulses of the full904 mixed902 linear equations. Original time potentials retained; finite observer fails current spatial/quadrature acceptance. Diagnostics do not certify source support, physical Green constraints, or complete872 feedback.',
        original_source_not_replaced=True,source_kernel_density_confusion_avoided=True,
        physical_observer_accuracy_accepted=False,continuous_error_certified=False,
        physical_retarded_Green_certified=False,complete872_quantum_feedback=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    v=run()
    if a.write:TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:close(v,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in v.items() if k!='source_file_sha256'},ensure_ascii=False,indent=2))
