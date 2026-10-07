"""913: differentiate actual finite material-path/Wilson record directly.
Reference second jets give the implicit path tangent. No curvature-source
formula is used to compute the direct derivative.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import independent_record_rebuild as rebuilt
work=rebuilt.work;mb=work.mb;loop=work.loop;ms=loop.ms
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'discrete_record_tangent_results.json'

def derivative(field,points,axes=()):
    axes=tuple(axes)
    if len(axes)<2:return field.evaluate(points,None if not axes else axes[0])
    theta=points[:,0]/field.T;nt=axes.count(0)
    powers=np.zeros((len(points),4),dtype=np.result_type(points,float))
    for n in range(nt,4):
        factor=1
        for j in range(nt):factor*=n-j
        powers[:,n]=factor*theta**(n-nt)/field.T**nt
    co=field.c.copy()
    for a in axes:
        if a:co*=1j*field.k[:,a-1,None,None]
    flat=co.reshape((len(co),-1));out=[]
    for start in range(0,len(points),128):
        phase=points[start:start+128,1:]@field.k.T
        v=np.cos(phase)@flat.real-np.sin(phase)@flat.imag
        out.append(np.einsum('bp,bpf->bf',powers[start:start+128],v.reshape((-1,4,58))))
    return np.concatenate(out)

def reference_jacobian(bg,points,variation=None,epsilon=0.):
    def jet(axes):
        val=derivative(bg,points,axes)
        if variation is not None:val=val+epsilon*derivative(variation,points,axes)
        return val
    g,phi,A=mb.unpack(jet(()));dg,dphi,dA=mb.unpack(np.stack([jet((i,)) for i in range(4)],axis=1))
    cache={(i,j):jet((i,j)) for i in range(4) for j in range(i,4)}
    _,ddphi,ddA=mb.unpack(np.stack([np.stack([cache[tuple(sorted((i,j)))] for j in range(4)],axis=1) for i in range(4)],axis=1))
    data=ms.geometry(g,phi,dphi,A,dA)
    columns=[ms.variation(data,phi,dphi,A,dg[:,i],dphi[:,i],ddphi[:,i],dA[:,i],ddA[:,i]) for i in range(4)]
    return np.stack(columns,axis=-1)

def exact_record_derivative(bg,field,src):
    p=src.jets;v=field.jets(src.points)
    dX=ms.variation(src.data,p['phi'],p['dphi'],p['A'],v['g'],v['phi'],v['dphi'],v['A'],v['dA'])
    xi=-np.einsum('bma,ba->bm',src.Jinv,dX)
    J=reference_jacobian(bg,src.points)
    Jcheck=bg.jacobian(src.points)
    assert np.max(abs(J-Jcheck))<2e-13
    h=1e-24
    dJ=reference_jacobian(bg,src.points.astype(complex)+1j*h*xi,field,1j*h).imag/h
    dv=-np.einsum('bma,ban,bn->bm',src.Jinv,dJ,src.velocity)
    movedA=v['A']+np.einsum('br,brma->bma',xi,p['dA'])
    colors=loop.connection(p['A']);dcolors=loop.connection(movedA);derivatives=[]
    for AA,dAA,channel in zip(colors,dcolors,src.channels):
        U,before,after,frechetdata=channel
        dB=(np.einsum('bm,bmij->bij',dv,AA)+np.einsum('bm,bmij->bij',src.velocity,dAA)).reshape((src.na,src.ns,3,3))/src.segments
        dE=loop.frechet(frechetdata,dB)
        dtr=np.sum(np.trace(after@dE@before,axis1=-2,axis2=-1),axis=1)
        derivatives.append(dtr)
    dc,df=derivatives;trc=np.trace(src.channels[0][0],axis1=-2,axis2=-1)
    values=np.sum(src.weights[:,None]*np.stack((df.real,df.imag,2*np.real(trc.conj()*dc)),axis=-1),axis=0)
    return values,dict(reference_Jacobian_two_jet_error=float(np.max(abs(J-Jcheck))),
        maximum_material_coordinate_derivative=float(np.max(abs(xi))),
        maximum_path_velocity_derivative=float(np.max(abs(dv))),
        implicit_reference_tangent_error=float(np.max(abs(np.einsum('bam,bm->ba',Jcheck,xi)+dX))),
        derivative_uses_second_jets_and_full_discrete_links=True)

def run():
    start=time.time();bg,field,initial=rebuilt.build_pair();rows=[]
    for segments in (8,16,32):
        src=loop.LoopSource(bg,2,segments);src.kernels()
        pairing=src.response(lambda x,p:field.jets(x));exact,stats=exact_record_derivative(bg,field,src)
        finite=[]
        if segments==8:
            for eps in (.001,.0005):
                plus,ps=rebuilt.rebuild(bg,field,src,eps);minus,ms_=rebuilt.rebuild(bg,field,src,-eps)
                d=(plus-minus)/(2*eps)
                finite.append(dict(epsilon=eps,rebuilt_derivative=d.tolist(),difference_from_discrete_tangent=(d-exact).tolist(),plus=ps,minus=ms_))
        row=dict(N=17,segments=segments,anchor_order=2,source_pairing=pairing,
            direct_discrete_tangent=exact.tolist(),direct_minus_curvature_source=(exact-np.array(pairing['total'])).tolist(),
            geometry=stats,small_step_independent_rebuilds=finite,elapsed_seconds=round(time.time()-start,3))
        rows.append(row);print(json.dumps(row,ensure_ascii=False),flush=True)
    return dict(round=913,date='2026-10-06',status='working',rows=rows,
        original_material_record_differentiated=True,curvature_source_not_used_in_direct_derivative=True,
        full_continuous_error_certified=False,full872_response_computed=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),HERE/'independent_record_rebuild.py',HERE/'receiver_response_readout.py',STAGE/'908/relational_loop_source.py',STAGE/'908/magnetic_reference_source.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('saved',args.write,flush=True)
