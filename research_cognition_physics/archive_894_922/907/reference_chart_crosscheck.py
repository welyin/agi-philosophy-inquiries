"""907 working: independently refine a candidate material-chart sign change.
Numerical finite-domain warning, not a certified continuum caustic or theory failure.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import material_reference_jets as j
TARGET=Path(__file__).with_name('reference_chart_crosscheck_results.json')

def snapshot(model,analytic,y):
    n=model.N;point=(0,n//4,n//8);F=model.rhs(y)
    directions=[F,*[{k:np.take(model.grad(a),i,axis=3) for k,a in y.items()} for i in range(3)]]
    derivatives=[j.derivative(analytic,y,v) for v in directions]
    ref=j.references(analytic,y)
    return {k:dict(value=ref[k][point].real.tolist(),Jacobian=np.stack([d[k][point] for d in derivatives],axis=-1).tolist(),determinant=float(np.linalg.det(np.stack([d[k][point] for d in derivatives],axis=-1)))) for k in ref}

def run():
    rows=[]
    with j.ev.ResearchRuntime(j.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        for N,steps in ((16,4),(24,2),(24,4),(32,4)):
            start=time.time();model=j.ev.Model(N,old);analytic=j.tangent.AnalyticModel(N,old);y,_=model.initial()
            for _ in range(steps):y=j.ev.rk4(model,y,.005/steps)
            rows.append(dict(N=N,steps=steps,time=.005,data=snapshot(model,analytic,y)))
            print('reference crosscheck',N,steps,rows[-1]['data']['new']['determinant'],'seconds',round(time.time()-start,2),flush=True)
    initial=json.loads(j.TARGET.read_text('utf-8'))['rows'][0]['data']['new']['determinant']
    return dict(round=907,status='working',formal_previous=906,cumulative_previous=3691,initial_N24_new_determinant=initial,rows=rows,
        new_determinant_sign_consistent_across_refinements=all(q['data']['new']['determinant']<0 for q in rows),
        old_determinant_sign_consistent_across_refinements=all(q['data']['old']['determinant']>0 for q in rows),
        numerical_chart_reversal_candidate_only=True,continuous_zero_certified=False,
        previous_notes_had_claimed_full_time_chart_certificate=False,common_physics_failure_proved=False,
        full_source_implementation_complete=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();v=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert v==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(v,ensure_ascii=False,indent=2))
