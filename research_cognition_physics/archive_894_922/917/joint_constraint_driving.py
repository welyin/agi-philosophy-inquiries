"""917: joint residual Ward drivers, same original background/A/receiver.
The source and background contacts are computed independently of the LHS.
Finite-coordinate differences only validate the differentiated residual identity.
"""
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
import receiver_stress_derivative as stress
c=stress.c
import validate_covariant_residual as validation
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'joint_constraint_driving_results.json'

def covariant_gravity_derivative(E,dE,z):
    return dE-np.einsum('brma,brn->bman',z['Gamma'],E)-np.einsum('brmn,bar->bman',z['Gamma'],E)

def ward(values,derivatives,p,base):
    z=base['z']
    divg=np.einsum('bma,bman->bn',z['ig'],covariant_gravity_derivative(values['gravity'],derivatives['gravity'],z))
    exchange=np.einsum('bA,bnA->bn',values['scalar'],base['D'])+np.einsum('bma,bnma->bn',values['YM'],base['F'])
    divY=np.einsum('bmma->ba',derivatives['YM'])+np.einsum('bmmr,bra->ba',z['Gamma'],values['YM'])
    divY+=sum(c.bracket(p['A'][:,m,:],values['YM'][:,m,:]) for m in range(4))
    charge=np.einsum('bA,aAB,bB->ba',values['scalar'],c.REP,p['phi'])
    return divg+exchange,divY+charge

def linear_euler(bg,field,x,par):
    p=c.aj.Cache(bg,x).jets();v=c.aj.Cache(field,x).jets();h=1e-24
    q=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    return {k:q[k].imag/h for k in ('gravity','scalar','YM')}

def run():
    start=time.time();par=c.parameters();bg,field,_=c.ex.previous.rebuilt.build_pair()
    src=c.ex.loop.LoopSource(bg,2,8);x=src.points[::64]
    cache=c.aj.Cache(bg,x);p=cache.jets();dp=[cache.jets((i,)) for i in range(4)]
    v=c.aj.Cache(field,x).jets();base,dbase=validation.coord_derivatives(p,dp,par)
    h=1e-24;shift=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
    delta={k:shift[k].imag/h for k in ('gravity','scalar','YM','D','F')}
    dmetric=shift['z']['ig'].imag/h;dGamma=shift['z']['Gamma'].imag/h
    E=base['gravity'];z=base['z'];covE=covariant_gravity_derivative(E,dbase['gravity'],z)
    contacts={}
    contacts['inverse_metric']=np.einsum('bma,bman->bn',dmetric,covE)
    contact_connection=-np.einsum('brma,brn->bman',dGamma,E)-np.einsum('brmn,bar->bman',dGamma,E)
    contacts['metric_connection']=np.einsum('bma,bman->bn',z['ig'],contact_connection)
    contacts['scalar_covariant_jet']=np.einsum('bA,bnA->bn',base['scalar'],delta['D'])
    contacts['curvature']=np.einsum('bma,bnma->bn',base['YM'],delta['F'])
    background=sum(contacts.values())
    internal={}
    internal['volume_connection']=np.einsum('bmmr,bra->ba',dGamma,base['YM'])
    internal['gauge_connection']=sum(c.bracket(v['A'][:,m,:],base['YM'][:,m,:]) for m in range(4))
    internal['scalar_representation']=np.einsum('bA,aAB,bB->ba',base['scalar'],c.REP,v['phi'])
    backgroundI=sum(internal.values())
    u,_,_=stress.weighted.receiver_modes();T=stress.stress_jet(bg,u,x)
    rA={k:delta[k].copy() for k in ('gravity','scalar','YM')};rA['gravity']-=T['stress']
    predicted=-T['divergence']-background;predictedI=-backgroundI
    rows=[]
    for step in (2e-6,1e-6):
        dr={k:[] for k in rA}
        for mu in range(4):
            xp=x.copy();xm=x.copy();xp[:,mu]+=step;xm[:,mu]-=step
            plus=linear_euler(bg,field,xp,par);minus=linear_euler(bg,field,xm,par)
            for k in dr:dr[k].append((plus[k]-minus[k])/(2*step))
        dr={k:np.stack(a,axis=1) for k,a in dr.items()};dr['gravity']-=T['derivative']
        actual,actualI=ward(rA,dr,p,base)
        row=dict(coordinate_check_step=step,
            diffeomorphism_identity_error=c.maximum(actual-predicted),internal_identity_error=c.maximum(actualI-predictedI),
            actual_diffeomorphism_driver_max=c.maximum(actual),actual_internal_driver_max=c.maximum(actualI),
            wrong_drop_background_contact_error=c.maximum(actual+T['divergence']),
            wrong_drop_receiver_divergence_error=c.maximum(actual+background),
            wrong_drop_internal_volume_contact_error=c.maximum(actualI+backgroundI-internal['volume_connection']))
        rows.append(row)
    assert rows[-1]['diffeomorphism_identity_error']<2e-8,rows
    assert rows[-1]['internal_identity_error']<2e-8,rows
    assert rows[-1]['wrong_drop_receiver_divergence_error']>1e-6
    old=json.loads((HERE/'receiver_stress_derivative_results.json').read_text('utf-8'))
    assert abs(c.maximum(T['divergence'])-old['actual_stress_divergence_sample_max'])<1e-13
    return dict(round=917,date='2026-10-06',status='working',samples=len(x),N=17,
      same_original_background_response_receiver=True,original_receiver_divergence_reproduced=True,
      diffeomorphism_background_contact_parts_max={k:c.maximum(a) for k,a in contacts.items()},
      diffeomorphism_background_contact_total_max=c.maximum(background),receiver_stress_divergence_max=c.maximum(T['divergence']),
      total_diffeomorphism_driver_max=c.maximum(predicted),total_diffeomorphism_driver_components_max=np.max(abs(predicted),axis=0).tolist(),
      internal_background_contact_parts_max={k:c.maximum(a) for k,a in internal.items()},total_internal_driver_max=c.maximum(predictedI),
      independent_differentiated_residual_checks=rows,
      full_background_contacts_including_volume_retained=True,
      derivative_check_is_not_error_certificate=True,full_support_constraint_driver_certified=False,
      actual_finite_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
      elapsed_seconds=round(time.time()-start,3),
      source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),HERE/'receiver_stress_derivative.py',STAGE/'916/covariant_joint_residual.py',STAGE/'916/validate_covariant_residual.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
