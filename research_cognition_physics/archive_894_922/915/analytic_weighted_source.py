"""915: replace the actual alpha finite difference with original analytic jets.
This resolves one algorithmic derivative component within the fixed interpolant;
physical, background, support and full-feedback errors are not certified here.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import analytic_material_jets as aj
import weighted_source_budget as budget
ex=aj.ex;prior=aj.prior;weighted=prior.weighted;HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'analytic_weighted_source_results.json'

def reference_second_check(bg,field,points):
    x=points[::64][:8];bc=aj.Cache(bg,x);vc=aj.Cache(field,x)
    bp=bc.reference_inputs();vp=vc.reference_inputs()
    J=np.stack([ex.oldref.reference_delta(bp,bc.reference_inputs((mu,))) for mu in range(4)],axis=-1)
    dx=ex.oldref.reference_delta(bp,vp);xi=np.linalg.solve(J,dx[...,None])[...,0]
    oldxi,olddxi,*_=ex.diff_jet(bg,field,x)
    rows=[]
    for step in (2e-6,1e-6):
        err=[]
        for mu in (0,1):
            xp=x.copy();xm=x.copy();xp[:,mu]+=step;xm[:,mu]-=step
            dp=ex.diff_jet(bg,field,xp)[1];dm=ex.diff_jet(bg,field,xm)[1]
            expected=np.moveaxis(xi.dd[mu],0,1)
            err.append(dict(axis=mu,maximum_difference=ex.maximum((dp-dm)/(2*step)-expected),analytic_second_derivative_max=ex.maximum(expected)))
        rows.append(dict(step=step,rows=err))
    return dict(old_xi_value_error=ex.maximum(xi.v-oldxi),old_xi_first_jet_error=ex.maximum(np.moveaxis(xi.d,0,1)-olddxi),independent_second_jet_diagnostics=rows)

def run():
    start=time.time();bg,field,_=ex.previous.rebuilt.build_pair();src=ex.loop.LoopSource(bg,2,8);src.kernels()
    u,z,mode=weighted.receiver_modes();w,dw,stats=weighted.weight_jets(u,z,src.points)
    b,da,jetstats=aj.extractor(bg,field,src.points)
    xi,dxi,J,p,v,alpha,rho,residual=b;base=(xi,dxi,p,v,alpha,rho,residual)
    analytic,gauge,corr=prior.ward_pair(src,base,da,w,dw)
    oldbase=weighted.extract_values(bg,field,src.points)
    base_errors=dict(xi=ex.maximum(xi-oldbase[0]),dxi=ex.maximum(dxi-oldbase[1]),alpha=ex.maximum(alpha-oldbase[4]))
    print(json.dumps(dict(stage='analytic_derivative_ready',elapsed_seconds=round(time.time()-start,3),base_errors=base_errors),ensure_ascii=False),flush=True)
    # Original FD algorithm, now compared to the independent analytic derivative
    # of exactly the same finite spacetime interpolants and source samples.
    step=1e-6;fd=[]
    for mu in range(4):
        xp=src.points.copy();xm=xp.copy();xp[:,mu]+=step;xm[:,mu]-=step
        fd.append((weighted.extract_values(bg,field,xp)[4]-weighted.extract_values(bg,field,xm)[4])/(2*step))
    fd=np.stack(fd,axis=1);delta=fd-da
    central,_,_=prior.ward_pair(src,base,fd,w,dw)
    coeff=budget.coefficients(src,w,dw);pred=[];component_bounds=[];uniform_bounds=[]
    for c in coeff:
        term=-c['Deta']*delta
        pred.append(np.sum(term));component_bounds.append(float(np.sum(abs(term))));uniform_bounds.append(float(np.sum(abs(c['Deta']))*ex.maximum(delta)))
    pred=np.array(pred);change=central-analytic
    assert ex.maximum(change-pred)<2e-15
    assert np.all(abs(change)<=np.array(component_bounds)+1e-15)
    old=json.loads((STAGE/'914/weighted_ward_pairing_results.json').read_text('utf-8'))
    oldsource=np.array(old['Ward_reduced_weighted_source']['real'])+1j*np.array(old['Ward_reduced_weighted_source']['imag'])
    assert ex.maximum(central-oldsource)<1e-11
    independent=reference_second_check(bg,field,src.points)
    result=dict(round=915,status='working',date='2026-10-06',N=17,samples=len(src.points),
       finite_interpolant_analytic_derivative_used=True,finite_difference_step_in_new_source=None,
       original_preparation_source_and_pairing_unchanged=True,
       analytic_weighted_source=weighted.complex_list(analytic),old_FD_on_same_analytic_base=weighted.complex_list(central),
       old_FD_reference_reproduction_error=ex.maximum(central-oldsource),
       FD_minus_analytic_derivative_max=ex.maximum(delta),FD_minus_analytic_derivative_by_axis=[ex.maximum(delta[:,i]) for i in range(4)],
       exact_linear_error_transport_prediction=weighted.complex_list(pred),actual_source_change=weighted.complex_list(change),
       error_transport_identity_residual=ex.maximum(change-pred),
       finite_array_componentwise_difference_envelope=component_bounds,finite_array_uniform_difference_envelope=uniform_bounds,
       base_errors=base_errors,jet_statistics=jetstats,independent_reference_checks=independent,
       conditional_non_diagonal_vA_matrix_element=weighted.complex_list(1j*np.array([6.,6.,.75])*analytic),
       pure_gauge_quadrature_defect_still_separate=weighted.complex_list(gauge),
       known_algorthmic_FD_truncation_removed_in_exact_arithmetic=True,
       floating_roundoff_interval_certified=False,physical_A_error_certified=False,actual_finite_observable_error_certified=False,
       full872_response_computed=False,full_goal_completed=False,elapsed_seconds=round(time.time()-start,3),
       source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),HERE/'analytic_material_jets.py',HERE/'weighted_source_budget.py',STAGE/'914/material_extractor_naturality.py',STAGE/'914/weighted_ward_pairing.py')})
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
