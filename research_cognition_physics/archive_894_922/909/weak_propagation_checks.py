"""Actual-source refinement and omission diagnostics, not a continuum certificate."""
from pathlib import Path
import argparse,json,time
import numpy as np
import weak_source_propagation as w
HERE=Path(__file__).resolve().parent
TARGET=HERE/'weak_propagation_checks_results.json'

def run():
    rows=[];projections=[]
    with w.ev.ResearchRuntime(w.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        for N in (16,24):
            start=time.time();bg=w.mb.Background(N);src=w.loop.LoopSource(bg,anchor_order=2,segments=8);src.kernels()
            model=w.ev.Model(N,old);analytic=w.tangent.AnalyticModel(N,old);y0,_=model.initial()
            projections.append(dict(N=N,**w.projection_check(src,N)))
            scenarios=[(2,1,False),(4,1,False)]
            if N==16:scenarios += [(4,2,False),(4,1,True)]
            for order,steps,omit in scenarios:
                nodes,items,_,meta=w.build_forcing(src,N,channel=2,time_order=order)
                row,y,u=w.propagate(model,analytic,y0,nodes,items,substeps=steps,omit_temporal=omit)
                row.update(metadata=meta,source_free_constraint_fourier_coefficients=w.constraint_modes(analytic,y,u))
                rows.append(row)
                print('case',N,order,steps,omit,'seconds',round(time.time()-start,2),
                      'obs',row['observations'],'lowH',row['source_free_constraint_fourier_coefficients']['H'],flush=True)
    return dict(round=909,date='2026-10-06',projection_checks=projections,rows=rows,
       argument_scope='Finite Fourier/time interpolation of the original908 adjoint-loop weak source, propagated by full904 mixed linear equations; includes negative temporal-potential omission control. No continuous support, discretized Ward, physical Green, or complete872 certificate.',
       continuous_error_certified=False,physical_projected_Green_certified=False,complete872_response=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Complete source-propagation diagnostics saved.' if a.write else 'Diagnostics complete.')
