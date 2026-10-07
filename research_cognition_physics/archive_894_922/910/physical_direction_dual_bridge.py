"""910 working: actual863 constrained direction transported by904, paired with908.
Reuses863 conserved future-test identity. No new inverse/Green theorem, full
response computation or continuous numerical error enclosure is asserted.
"""
from pathlib import Path
import sys,json,time,argparse,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
for n in ('908','904'):sys.path.insert(0,str(STAGE/n))
import coupled_boson_tangent as tangent
import material_background as mb
import relational_loop_source as loop
ev=tangent.base;TARGET=HERE/'physical_direction_dual_bridge_results.json'

class PhysicalDirection(mb.Background):
    def __init__(self,N=16,steps=2,T=.00025,drop=1e-14):
        self.N=N;self.steps=steps;self.T=T;self.drop=drop
        a,b,c=.31,.27,.21
        da=-(a*a+b*b)*c*c/(4*a*(b*b+c*c/4));dc=c
        with ev.ResearchRuntime(ev.Layout()).installed():
            import joint_reference_constraint_strata as old
            model=ev.Model(N,old);analytic=tangent.AnalyticModel(N,old);y0,_=model.initial()
            u0={k:np.zeros_like(v) for k,v in y0.items()};u0['A'][...,0,0]=da;u0['A'][...,2,3]=dc
            self.initial_linear_constraints=tangent.tangent_constraints(analytic,y0,u0)
            ends=[];last=[]
            for sign in (-1,1):
                y={k:a.copy() for k,a in y0.items()};u={k:a.copy() for k,a in u0.items()}
                for _ in range(steps):y,u=tangent.step(model,analytic,y,u,sign*T/steps)
                ends.append((mb.pack(u),mb.pack(analytic.jvp(y,u))))
                last.append(tangent.tangent_constraints(analytic,y,u))
            um,dm=ends[0];up,dp=ends[1]
            coeff=np.stack(((up+um)/2-T*(dp-dm)/4,3*(up-um)/4-T*(dp+dm)/4,T*(dp-dm)/4,(T*(dp+dm)-(up-um))/4),axis=-2)
            ft=np.fft.fftn(coeff,axes=(0,1,2)).reshape((-1,4,58))/N**3
            waves=model.k.reshape((-1,3));keep=np.max(abs(ft),axis=(1,2))>=drop
            self.c=ft[keep];self.k=waves[keep]
            self.tail0=float(np.max(np.sum(abs(ft[~keep]),axis=(0,1))))
            self.endpoint_linear_constraints=last
        self.initial_direction=dict(Ax_color_T1=da,Az_color_T4=dc,
            magnetic_energy_first_derivative=2*a*(b*b+c*c/4)*da+(a*a+b*b)*c*dc/2)
    def as_variation(self,points,unused):return self.jets(points)

def run():
    rows=[]
    old=json.loads((STAGE/'908/physical_family_source_checks_results.json').read_text('utf-8'))
    for N in (16,24):
        start=time.time();bg=mb.Background(N);src=loop.LoopSource(bg,anchor_order=2,segments=8);src.kernels()
        direction=PhysicalDirection(N)
        response=src.response(direction.as_variation)
        row=dict(N=N,steps=direction.steps,background_halfwindow=direction.T,source_profile_unchanged=True,
                 initial_direction=direction.initial_direction,initial_linear_constraints=direction.initial_linear_constraints,
                 endpoint_linear_constraints=direction.endpoint_linear_constraints,
                 finite_interpolant_discarded_tail=direction.tail0,retained_direction_modes=len(direction.k),
                 full_source_pairing=response,conditional_future_read_coefficient=.75*response['total'][2])
        if N==16:
            row['previous_finite_family_comparison']=[dict(epsilon=q['epsilon'],
                finite_family_source_minus_exact_tangent=q['full_source']['total'][2]-response['total'][2],
                rebuilt_record_difference_minus_exact_tangent=q['independently_rebuilt_record_derivative'][2]-response['total'][2]) for q in old['rows']]
        rows.append(row);print(json.dumps(row,ensure_ascii=False),flush=True);print('seconds',round(time.time()-start,2),flush=True)
    return dict(round=910,status='working',date='2026-10-06',formal_previous=909,cumulative_previous=3694,
       rows=rows,actual_original_constrained_tangent_implemented=True,
       future_test_identity_inherited_from_863=True,original_source_unchanged=True,
       future_test_kernel_array_computed=False,continuous_actual_support_certified=False,
       actual_continuous_response_error_bound=False,complete872_response=False,full_goal_completed=False,
       source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                     (Path(__file__),STAGE/'904/coupled_boson_tangent.py',STAGE/'908/material_background.py',STAGE/'908/relational_loop_source.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    v=run()
    if a.write:TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
