"""Hold the actual sourced trajectory fixed and vary only observer quadrature."""
from pathlib import Path
import json,argparse
import weak_source_propagation as w
HERE=Path(__file__).resolve().parent;TARGET=HERE/'observer_quadrature_check_results.json'

def run():
    bg=w.mb.Background(16);src=w.loop.LoopSource(bg,anchor_order=2,segments=8);src.kernels()
    nodes,items,_,meta=w.build_forcing(src,16,channel=2,time_order=4)
    with w.ev.ResearchRuntime(w.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=w.ev.Model(16,old);analytic=w.tangent.AnalyticModel(16,old);y0,_=model.initial()
        _,y,u=w.propagate(model,analytic,y0,nodes,items,observer=False)
        rows=[]
        for order in (2,3,4,6):
            row=w.readonly_scalar_observations(model,analytic,y,u,order=order)
            rows.append(row);print(row,flush=True)
    return dict(round=909,date='2026-10-06',N=16,source_metadata=meta,rows=rows,
         original_source_and_evolution_fixed=True,observer_profile_fixed=True,
         diagnostic_not_a_continuous_prediction=True,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
