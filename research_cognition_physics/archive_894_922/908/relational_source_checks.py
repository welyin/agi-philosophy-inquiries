"""908 complete finite-domain original-source checks. No rigorous PDE error claim."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from material_background import Background
import relational_loop_source as loop
import magnetic_reference_source as ms
HERE=Path(__file__).resolve().parent;TARGET=HERE/'relational_source_checks_results.json'
def jet_checks(q):
    ids=np.arange(0,len(q.points),max(1,len(q.points)//31))[:31]
    p={k:v[ids] for k,v in q.jets.items()};data=ms.geometry(p['g'],p['phi'],p['dphi'],p['A'],p['dA']);rng=np.random.default_rng(908)
    error=0.;duality=0.;omitted_metric=0.;omitted_time=0.;omitted_clock=0.
    for _ in range(12):
        v={k:rng.normal(size=p[k].shape)*.003 for k in ('g','phi','dphi','A','dA')};v['g']=(v['g']+v['g'].swapaxes(-1,-2))/2
        h=1e-24;z={k:p[k].astype(complex)+1j*h*v[k] for k in v}
        exact=ms.geometry(z['g'],z['phi'],z['dphi'],z['A'],z['dA'])['X'].imag/h
        pred=ms.variation(data,p['phi'],p['dphi'],p['A'],v['g'],v['phi'],v['dphi'],v['A'],v['dA'])
        error=max(error,float(np.max(abs(exact-pred))))
        # Original loop coefficients, not unrelated arbitrary synthetic currents.
        for krec in range(3):
            current=q.currents[krec,ids];k=q.ks[krec][ids]
            source=ms.first_jet_source(data,p['phi'],p['dphi'],p['A'],current,k)
            direct=np.sum(current*v['A'],axis=(-2,-1))+np.sum(k*exact,axis=-1)
            duality=max(duality,float(np.max(abs(ms.pair(source,v)-direct))))
            omitted_metric=max(omitted_metric,float(np.max(abs(np.sum(source['g']*v['g'],axis=(-2,-1))))))
            omitted_time=max(omitted_time,float(np.max(abs(np.sum(source['dphi'][:,0]*v['dphi'][:,0],axis=-1)+np.sum(source['dA'][:,0]*v['dA'][:,0],axis=(-2,-1))))))
            omitted_clock=max(omitted_clock,float(np.max(abs(np.sum(source['dphi']*v['dphi'],axis=(-2,-1))))))
    assert error<1e-10 and duality<1e-10
    assert min(omitted_metric,omitted_time,omitted_clock)>1e-16
    return dict(independent_offshell_jets=12*len(ids),reference_derivative_max_error=error,source_adjoint_pairing_max_error=duality,
        omitted_metric_local_pair_defect=omitted_metric,omitted_time_derivative_local_pair_defect=omitted_time,
        omitted_clock_derivative_local_pair_defect=omitted_clock,independent_time_and_connection_jets_retained=True)

def row(bg,order,segments,dojets=False):
    begin=time.time();q=loop.LoopSource(bg,order,segments);q.kernels()
    variations=dict(metric=loop.metric_variation,clock=loop.clock_variation,probe=loop.probe_variation,
        color_gauge=loop.color_gauge_variation,affine_diffeomorphism=loop.affine_diffeo_variation(bg),time_translation=loop.translation_variation(bg,[.0001,0,0,0]))
    responses={k:q.response(v) for k,v in variations.items()}
    assert max(abs(v) for v in responses['color_gauge']['total'])<1e-9
    assert max(abs(v) for v in responses['affine_diffeomorphism']['total'])<1e-9
    assert max(abs(v) for v in responses['probe']['total'])>1e-4
    assert q.stats['clock_margin']>0 and q.stats['minimum_Jacobian_det']>0 and q.stats['maximum_sampled_receiver_radius']<.75
    result=dict(N=bg.N,background_steps_each_direction=bg.steps,interpolation_half_time=bg.T,Fourier_drop_threshold=bg.drop,
        anchor_order=order,segments_per_edge=segments,retained_spatial_modes=len(bg.k),records=q.records.tolist(),
        inverse_and_support_diagnostics=q.stats,source_responses=responses,finite_Fourier_discard_field_bound_on_slab=bg.tail0,
        finite_Fourier_discard_spatial_derivative_bounds=bg.tail1,
        same_fixed_material_curve_and_profile=True)
    if dojets:result['independent_jet_checks']=jet_checks(q)
    print('source row',bg.N,order,segments,'seconds',round(time.time()-begin,2),'probe',responses['probe']['total'][1],flush=True)
    return result,q

def run():
    rows=[];bg=Background(16)
    for seg in (8,16,32):
        r,q=row(bg,2,seg,seg==16);rows.append(r)
        if seg==16:base=q
    r,q=row(bg,3,16);rows.append(r)
    fd=[];analytic=np.array(base.response(loop.probe_variation)['total'])
    for eps in (1e-4,5e-5,2.5e-5,1.25e-5):
        plus=loop.LoopSource(bg,2,16,loop.probe_variation,eps).records
        minus=loop.LoopSource(bg,2,16,loop.probe_variation,-eps).records
        v=(plus-minus)/(2*eps)
        fd.append(dict(epsilon=eps,derivative=v.tolist(),difference_from_full_source=(v-analytic).tolist()))
        print('finite rebuilt probe',eps,'error',float(np.max(abs(v-analytic))),flush=True)
    errors=[abs(r['difference_from_full_source'][1]) for r in fd]
    assert all(3.5<a/b<4.5 for a,b in zip(errors,errors[1:]))
    assert errors[-1]<1e-7
    e=1e-4;plus=loop.LoopSource(bg,2,16,loop.metric_variation,e).records;minus=loop.LoopSource(bg,2,16,loop.metric_variation,-e).records
    metricfd=(plus-minus)/(2*e);ma=np.array(base.response(loop.metric_variation)['total'])
    assert np.max(abs(metricfd-ma))<1e-8
    bgfine=Background(24);r,q=row(bgfine,2,16);rows.append(r)
    bgtime=Background(16,T=.000125);r,q=row(bgtime,2,16);rows.append(r)
    comparisons={}
    def vals(r):return np.r_[r['records'],*[r['source_responses'][k]['total'] for k in ('metric','clock','probe')]]
    for name,a,b in (('path8_to16',0,1),('path16_to32',1,2),('anchor2_to3',1,3),('space16_to24',1,4),('interpolation_time_halving',1,5)):
        comparisons[name]=float(np.max(abs(vals(rows[a])-vals(rows[b]))))
    assert comparisons['path16_to32']<comparisons['path8_to16']
    names=('magnetic_reference_source.py','material_background.py','relational_loop_source.py','drafts/initial_anchor_failure.json','drafts/source_implementation_working.md')
    return dict(round=908,date='2026-10-06',all_checks_passed=True,
        argument_scope='Original fixed material loop, compact4D anchor profile, three quotient-group records and all off-shell first-jet weak sources on the actual902 short-time numerical background. No full PDE error or strong-source/feedback completion.',
        design=dict(L=loop.L,widths=loop.WIDTH.tolist(),profile='normalized product exp(-1/(1-u^2)) on (-1,1)^4 transported by fixed907 initial Jacobian',
            curve='fixed901 dyadic initial psi witness image of the original physical rectangle',reference_X=loop.DESIGN_X.tolist(),reference_J=loop.DESIGN_J.tolist()),
        rows=rows,refinement_comparisons=comparisons,rebuilt_probe_finite_differences=fd,
        rebuilt_metric_finite_difference=dict(epsilon=e,derivative=metricfd.tolist(),full_source=ma.tolist(),difference=(metricfd-ma).tolist()),
        source_file_sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        full_first_jet_weak_source_implemented=True,original_three_records_and_U1_retained=True,
        dynamic_background_error_certified=False,full_compact_support_certified=False,strong_spacetime_source_array_computed=False,
        actual872_forced_response_computed=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_file_sha256','design')},ensure_ascii=False,indent=2))
