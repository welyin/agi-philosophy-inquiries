"""913: independently rebuild deformed original material loops for actual912 A.
This checks the numerical observable dictionary, not continuum error or full872.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import receiver_response_readout as work
r=work.r;response=work.response;mb=work.mb;loop=work.loop
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'independent_record_rebuild_results.json'

def build_pair(N=17,steps=2,T=.00025):
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(N,old);analytic=response.tangent.AnalyticModel(N,old)
        y0,u0,A0,initial=response.initial(model,analytic,old);ends=[]
        for sign in (-1,1):
            y={k:v.copy() for k,v in y0.items()};u=u0.copy();A={k:v.copy() for k,v in A0.items()}
            for _ in range(steps):y,u,A=response.step(model,analytic,y,u,A,sign*T/steps)
            dy,du,dA=response.rhs(model,analytic,y,u,A)
            ends.append(((mb.pack(y),mb.pack(dy)),(mb.pack(A),mb.pack(dA))))
        bg=work.interpolant(N,T,ends[0][0],ends[1][0],model.k,1e-11)
        field=work.interpolant(N,T,ends[0][1],ends[1][1],model.k,1e-14)
    bg.J0=bg.jacobian(mb.CENTER[None])[0];bg.X0=bg.reference(mb.CENTER[None])[0]
    return bg,field,initial

def add_interpolants(bg,field,epsilon):
    # The actual affine spacetime fields y(x)+epsilon A(x), not a new solution.
    out=object.__new__(mb.Background);out.N=bg.N;out.T=bg.T;out.steps=bg.steps
    out.k=field.k;out.c=epsilon*field.c.copy();ids={tuple(k):i for i,k in enumerate(out.k)}
    for k,c in zip(bg.k,bg.c):out.c[ids[tuple(k)]]+=c
    out.J0=bg.J0;out.X0=bg.X0
    return out

def records_from_path(jets,velocity,na,ns,segments,weights):
    color,full=loop.connection(jets['A']);traces=[]
    for AA in (color,full):
        B=np.einsum('bm,bmij->bij',velocity,AA).reshape((na,ns,3,3))/segments
        E,_=loop.exponential(B);U,*_=loop.transports(E)
        traces.append(np.trace(U,axis1=-2,axis2=-1))
    c,f=traces
    return np.sum(weights[:,None]*np.stack((f.real,f.imag,abs(c)**2-1),axis=-1),axis=0)

def rebuild(bg,field,src,epsilon):
    deformed=add_interpolants(bg,field,epsilon)
    anchor,weights,_=loop.anchors(src.anchor_order);ell,dell,_=loop.design_curve(src.segments)
    labels=(anchor[:,None,:]+ell[None,:,:]).reshape((-1,4));x=src.points.copy()
    norms=[]
    # Chord iteration with the old Jacobian is only a solver preconditioner.
    # Final records use the actual deformed Jacobian and deformed fields.
    for it in range(18):
        residual=deformed.reference(x)-labels
        dx=np.einsum('bma,ba->bm',src.Jinv,residual);x-=dx
        norms.append(float(np.max(abs(dx))))
        if norms[-1]<2e-11:break
    residual=deformed.reference(x)-labels
    coordinate_residual=float(np.max(abs(np.einsum('bma,ba->bm',src.Jinv,residual))))
    assert coordinate_residual<2e-8,(epsilon,coordinate_residual,norms)
    assert np.max(abs(x[:,0]))<bg.T
    J=deformed.jacobian(x);velocity=np.linalg.solve(J,np.tile(dell,(src.na,1))[...,None])[...,0]
    jets=deformed.jets(x)
    rec=records_from_path(jets,velocity,src.na,src.ns,src.segments,weights)
    return rec,dict(epsilon=epsilon,iterations=it+1,coordinate_residual=coordinate_residual,
        maximum_coordinate_shift=float(np.max(abs(x-src.points))),
        minimum_Jacobian_det=float(np.min(np.linalg.det(J))),time_min=float(x[:,0].min()),time_max=float(x[:,0].max()))

def run():
    start=time.time();bg,field,init=build_pair();rows=[]
    def variation(points,p):return field.jets(points)
    for segments,epsilons in ((8,(.02,.01)),(16,(.01,))):
        src=loop.LoopSource(bg,2,segments);src.kernels();pred=src.response(variation)
        # Rebuilding at zero checks the record and source use precisely one path convention.
        base=records_from_path(src.jets,src.velocity,src.na,src.ns,segments,src.weights)
        assert np.max(abs(base-src.records))<1e-14
        for eps in epsilons:
            plus,ps=rebuild(bg,field,src,eps);minus,ms=rebuild(bg,field,src,-eps)
            derivative=(plus-minus)/(2*eps);difference=derivative-np.array(pred['total'])
            row=dict(N=17,segments=segments,anchor_order=2,epsilon=eps,
                source_pairing=pred,rebuilt_record_derivative=derivative.tolist(),
                derivative_minus_source_pairing=difference.tolist(),plus=ps,minus=ms,
                elapsed_seconds=round(time.time()-start,3))
            rows.append(row);print(json.dumps(row,ensure_ascii=False),flush=True)
    old=json.loads(work.TARGET.read_text('utf-8'))['rows'][0]
    assert np.allclose(rows[0]['source_pairing']['total'],old['full_pairing']['total'],rtol=1e-10,atol=1e-18)
    return dict(round=913,date='2026-10-06',status='working',rows=rows,
        original913_pairing_reproduced=True,same_material_labels_and_source_profile=True,
        deformed_material_coordinates_and_tangents_recomputed=True,
        direct_nonlinear_records_used_no_first_jet_source_in_rebuild=True,
        continuum_response_error_certified=False,full872_response_computed=False,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),HERE/'receiver_response_readout.py',STAGE/'908/relational_loop_source.py',STAGE/'912/receiver_forced_response.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('saved',args.write,flush=True)
