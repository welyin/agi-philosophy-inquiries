"""Check the proposed diagnostic observer's reference chart, not microscopic continuity."""
from pathlib import Path
import itertools,json,time,argparse
import numpy as np
import weak_source_propagation as w
HERE=Path(__file__).resolve().parent;TARGET=HERE/'observer_chart_check_results.json'

def run():
    rows=[]
    with w.ev.ResearchRuntime(w.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        for N in (16,24):
            model=w.ev.Model(N,old);analytic=w.tangent.AnalyticModel(N,old);y,_=model.initial()
            y=w.ev.rk4(model,y,.00015);F=model.rhs(y)
            directions=[F,*[{k:np.take(model.grad(a),i,axis=3) for k,a in y.items()} for i in range(3)]]
            Jgrid=np.stack([w.refs.derivative(analytic,y,v)['old'] for v in directions],axis=-1)
            nodes=np.linspace(-.025,.025,9);offset=np.array(list(itertools.product(nodes,repeat=3)))
            pts=w.mb.CENTER[None,1:]+offset;fm=w.FourierMap(N,pts);J=fm.evaluate(Jgrid);det=np.linalg.det(J)
            lo=int(np.argmin(det));hi=int(np.argmax(det))
            row=dict(N=N,time=.00015,observer_radius=.025,sampled_min_det=float(det[lo]),sampled_max_det=float(det[hi]),
                     min_position_offset=offset[lo].tolist(),max_position_offset=offset[hi].tolist(),
                     finite_interpolant_sign_change=bool(det[lo]<0<det[hi]),continuous_physical_chart_certificate=False)
            if row['finite_interpolant_sign_change']:
                a=pts[lo].copy();b=pts[hi].copy()
                for _ in range(45):
                    mid=(a+b)/2;d=float(np.linalg.det(w.FourierMap(N,mid[None]).evaluate(Jgrid)[0]))
                    if d<0:a=mid
                    else:b=mid
                row.update(interpolant_zero_bracket_diameter=float(np.linalg.norm(a-b)),
                           interpolant_zero_midpoint=((a+b)/2).tolist())
            rows.append(row);print(row,flush=True)
    return dict(round=909,date='2026-10-06',rows=rows,
         finding_scope='Sign change, if present, is a numerical certificate of a zero of the continuous finite Fourier interpolant, not a validated zero of the original PDE reference. The full support of the proposed observer is therefore not accepted.',
         original908_source_apparatus_changed=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
