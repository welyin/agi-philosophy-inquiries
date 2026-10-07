"""922: same kinematic background defect and the full original readout derivative.
The background is varied as ONE finite Fourier/Hermite field. All original
material labels, old reference, receiver modes and response direction are kept.
This is a partial joint error contribution, not a full observable certificate.
"""
from pathlib import Path
import sys,argparse,json,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'921'))
import kinematic_defect_readout as previous
c=previous.c;weak=previous.weak;r=previous.r;mb=previous.mb
rebuilt=c.ex.previous.rebuilt;discrete=c.aj.previous
weighted=c.aj.prior.weighted
TARGET=HERE/'background_readout_variation_results.json'


def projected_background_residual(fam,bg,t,M=25):
    p=previous.first_jets(weak.GridField(bg,t,M))
    mc=weak.GridField(fam.mc,t,M).get()
    base=previous.kinetic_residual(p,mc,fam.par);N=fam.model.N
    ik=np.rint(fam.model.k).astype(int);ids=tuple(ik[...,i]%M for i in range(3))
    force={k:np.zeros_like(a) for k,a in fam.y(0.).items()}
    for k,a in base.items():
        co=np.fft.fftn(a.reshape((M,M,M)+a.shape[1:]),axes=(0,1,2))/M**3
        force[k]=-np.fft.ifftn(co[ids]*N**3,axes=(0,1,2)).real
    return force,dict(t=t,M=M,base_residual_sample_max={k:c.maximum(a) for k,a in base.items()},
        projected_background_force_max={k:c.maximum(force[k]) for k in ('phi','A')})


def readout(bg,direction,u,z,baseline_points=None):
    src=c.ex.loop.LoopSource(bg,2,8);src.kernels()
    pair=src.response(lambda x,p:direction.jets(x))
    direct,details=discrete.exact_record_derivative(bg,direction,src)
    w,dw,wstats=weighted.weight_jets(u,z,src.points)
    b,da,jstats=c.aj.extractor(bg,direction,src.points)
    xi,dxi,J,p,v,alpha,rho,residual=b;base=(xi,dxi,p,v,alpha,rho,residual)
    ward,gauge,corr=c.aj.prior.ward_pair(src,base,da,w,dw)
    projected=weighted.weighted_pair(src,weighted.projected_inputs(base,da),w,dw)
    err=c.maximum(ward-projected-gauge);assert err<1e-12,err
    out=dict(ordinary_source=pair,discrete_record_derivative=direct.tolist(),
        weighted_source=weighted.complex_list(ward),weighted_principal_correction=weighted.complex_list(corr),
        gauge_quadrature_defect=weighted.complex_list(gauge),Ward_identity_error=err,
        source_geometry=src.stats,discrete_derivative_checks=details,weight=wstats,extractor=jstats,
        old_reference_J_condition_sample_max=float(np.linalg.cond(J).max()),
        maximum_xi=c.maximum(xi),maximum_dalpha=c.maximum(da))
    if baseline_points is not None:out['maximum_material_point_shift']=c.maximum(src.points-baseline_points)
    return out,src


def vec(row,key):
    val=row[key]
    if key=='ordinary_source':return np.array(val['total'])
    if isinstance(val,dict):return np.array(val['real'])+1j*np.array(val['imag'])
    return np.array(val)


def run():
    start=time.time();bg,A,_=rebuilt.build_pair();T=bg.T
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        fam=previous.Family(T,old);forcing=[];force_stats=[]
        for t in (-T,0.,T):
            f,s=projected_background_residual(fam,bg,t);forcing.append(f);force_stats.append(s)
        b,evolution=previous.integrate(fam,forcing)
    print(json.dumps(dict(stage='background_field',force=force_stats,evolution=evolution,elapsed=time.time()-start)),flush=True)
    u,z,mode=weighted.receiver_modes()
    base,src0=readout(bg,A,u,z)
    old=json.loads((STAGE/'915/analytic_weighted_source_results.json').read_text('utf-8'))['analytic_weighted_source']
    old=np.array(old['real'])+1j*np.array(old['imag']);reproduction=c.maximum(vec(base,'weighted_source')-old)
    assert reproduction<1e-12,reproduction
    # These are differentiation parameters, not physical preparations or cutoffs.
    rows=[];derivatives=[]
    for eps in (.5,.25):
        signed=[]
        for sign in (-1,1):
            changed=rebuilt.add_interpolants(bg,b,sign*eps)
            value,_=readout(changed,A,u,z,src0.points)
            row=dict(epsilon=sign*eps,**value);rows.append(row);signed.append(value)
            print(json.dumps(dict(stage='deformed_readout',epsilon=sign*eps,weighted=value['weighted_source'],shift=value['maximum_material_point_shift'],elapsed=time.time()-start)),flush=True)
        values={}
        for key in ('ordinary_source','discrete_record_derivative','weighted_source','weighted_principal_correction','gauge_quadrature_defect'):
            dv=(vec(signed[1],key)-vec(signed[0],key))/(2*eps)
            values[key]=weighted.complex_list(dv) if np.iscomplexobj(dv) else dv.tolist()
        derivatives.append(dict(epsilon=eps,**values))
    # Independent order of differentiation for the actual FINITE Wilson record:
    # D_b (D_A O) is checked by D_A (D_b O), which differentiates b, not A.
    reverse=[]
    for h in (.001,.0005):
        vals=[]
        for sign in (-1,1):
            changed=rebuilt.add_interpolants(bg,A,sign*h)
            src=c.ex.loop.LoopSource(changed,2,8)
            value,stats=discrete.exact_record_derivative(changed,b,src);vals.append(value)
        value=(vals[1]-vals[0])/(2*h)
        reverse.append(dict(response_direction_epsilon=h,derivative=value.tolist()))
    direct=np.array(derivatives[-1]['discrete_record_derivative']);rev=np.array(reverse[-1]['derivative'])
    # Numerical diagnostics, not rigorous error enclosures.
    steps={k:c.maximum(vec(derivatives[0],k)-vec(derivatives[1],k)) for k in ('discrete_record_derivative','weighted_source','weighted_principal_correction','gauge_quadrature_defect')}
    steps['ordinary_source']=c.maximum(np.array(derivatives[0]['ordinary_source'])-np.array(derivatives[1]['ordinary_source']))
    prior921=json.loads((STAGE/'921/kinematic_defect_readout_results.json').read_text('utf-8'))
    partial=vec(derivatives[-1],'weighted_source')+vec(prior921,'original_weighted_Ward_source')
    return dict(round=922,date='2026-10-06',N=17,projection_quadrature_M=25,time_half_window=T,
      force_is_negative_background_kinematic_residual=True,force_slots=['phi','A'],all104_variables_propagated=True,
      initial_error_component_chosen_zero=True,finite_interpolant_initial_value_not_identified_with_exact_initial_data=True,
      force_samples=force_stats,background_component=evolution,baseline=base,baseline915_reproduction_error=reproduction,
      deformed_readouts=rows,centered_background_derivatives=derivatives,derivative_step_discrepancies_not_bounds=steps,
      independent_reverse_order_derivatives=reverse,mixed_discrete_derivative_order_difference=(direct-rev).tolist(),
      partial_921_tangent_plus_922_background_weighted=weighted.complex_list(partial),
      original_mode_information=mode,
      same_original_preparation_and_readout=True,moving_material_path_source_reference_and_weight_evaluation_included=True,
      one_integrable_background_field_used=True,independent_local_jet_patch_used=False,
      background_derivative_Hessian_readout_included=True,
      background_induced_D2F_propagation_included=False,background_induced_receiver_and_stress_changes_included=False,
      full_background_residual_included=False,initial_geometric_constraint_error_included=False,
      finite_operator_dictionary_error_certified=False,physical_GPi_inverse_implemented=False,
      full_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
      contribution_scope='D_Y j_Y(A) and D_Y K(Y,u,z;A)[b] only; A,u,z fixed as spacetime fields',
      elapsed_seconds=round(time.time()-start,3),source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),STAGE/'921/kinematic_defect_readout.py',STAGE/'913/independent_record_rebuild.py',STAGE/'913/discrete_record_tangent.py',STAGE/'915/analytic_material_jets.py',STAGE/'914/weighted_ward_pairing.py')})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('deformed_readouts','baseline','original_mode_information')},ensure_ascii=False,indent=2))
