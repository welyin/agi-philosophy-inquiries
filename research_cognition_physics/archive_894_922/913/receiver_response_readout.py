"""913 working: same912 forced receiver response paired with908 full sources.
Finite interpolation/quadrature diagnostics, not a physical-Pi/error certificate.
"""
from pathlib import Path
import sys,json,hashlib,argparse,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
for n in ('912','908'):sys.path.insert(0,str(STAGE/n))
import receiver_forced_response as response
import material_background as mb
import relational_loop_source as loop
r=response.r;TARGET=HERE/'receiver_response_readout_results.json'

def interpolant(N,T,minus,plus,waves,drop):
    ym,dm=minus;yp,dp=plus
    co=np.stack(((yp+ym)/2-T*(dp-dm)/4,3*(yp-ym)/4-T*(dp+dm)/4,T*(dp-dm)/4,(T*(dp+dm)-(yp-ym))/4),axis=-2)
    ft=np.fft.fftn(co,axes=(0,1,2)).reshape((-1,4,58))/N**3;k=waves.reshape((-1,3))
    keep=np.max(abs(ft),axis=(1,2))>=drop
    out=object.__new__(mb.Background);out.N=N;out.T=T;out.steps=2;out.drop=drop;out.c=ft[keep];out.k=k[keep]
    out.tail0=float(np.max(np.sum(abs(ft[~keep]),axis=(0,1))))
    out.tail1=[float(np.max(np.sum(abs(ft[~keep])*abs(k[~keep,i,None,None]),axis=(0,1)))) for i in range(3)]
    out.tailt=float(np.max(np.sum(abs(ft[~keep])*np.arange(4)[None,:,None]/T,axis=(0,1))))
    return out

def actual(N=17,steps=2,T=.00025):
    start=time.time()
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(N,old);analytic=response.tangent.AnalyticModel(N,old)
        y0,u0,A0,initial=response.initial(model,analytic,old)
        ends=[];constraints=[]
        for sign in (-1,1):
            y={k:v.copy() for k,v in y0.items()};u=u0.copy();A={k:v.copy() for k,v in A0.items()}
            for _ in range(steps):y,u,A=response.step(model,analytic,y,u,A,sign*T/steps)
            dy,du,dA=response.rhs(model,analytic,y,u,A)
            ends.append(((mb.pack(y),mb.pack(dy)),(mb.pack(A),mb.pack(dA))))
            constraints.append(response.total_constraints(model,analytic,y,u,A))
        bg=interpolant(N,T,ends[0][0],ends[1][0],model.k,1e-11)
        Afield=interpolant(N,T,ends[0][1],ends[1][1],model.k,1e-14)
    bg.J0=bg.jacobian(mb.CENTER[None])[0];bg.X0=bg.reference(mb.CENTER[None])[0]
    src=loop.LoopSource(bg,anchor_order=2,segments=8);src.kernels()
    def variation(points,unused):return Afield.jets(points)
    pair=src.response(variation)
    # Conditional continuous one-mode Dirac response identity, not a finished record.
    diagonal_imag=np.array([6.,6.,.75])*np.array(pair['total'])  # 2m=1 for870 m=.5
    row=dict(N=N,steps=steps,T=T,source_anchor_order=2,source_segments=8,
        original_source_and_preparation_unchanged=True,full_pairing=pair,
        conditional_vA_diagonal_imaginary_coefficients=diagonal_imag.tolist(),
        source_stats=src.stats,initial=initial,endpoint_constraint_diagnostics=constraints,
        finite_background_modes=len(bg.k),finite_response_modes=len(Afield.k),
        discarded_response_interpolant_tail=dict(value=Afield.tail0,space=Afield.tail1,time=Afield.tailt),
        elapsed_seconds=round(time.time()-start,3),
        physical_A_certified=False,actual_continuous_pairing_error_certified=False,
        full_vA_field_computed=False,full872_response_computed=False)
    print(json.dumps(row,ensure_ascii=False),flush=True)
    return row

def run():
    rows=[actual(N) for N in (17,25)]
    return dict(round=913,status='working',date='2026-10-06',formal_previous=912,cumulative_previous=3697,
        rows=rows,mesh_difference_not_error_bound=np.abs(np.array(rows[1]['full_pairing']['total'])-np.array(rows[0]['full_pairing']['total'])).tolist(),
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),STAGE/'912/receiver_forced_response.py',STAGE/'908/material_background.py',STAGE/'908/relational_loop_source.py')},
        full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'mesh_difference_not_error_bound':out['mesh_difference_not_error_bound'],'saved':args.write},indent=2))
