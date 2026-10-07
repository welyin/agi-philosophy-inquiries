"""824 working check: original scalar first-jet preservation versus acceleration.

The sources are the frozen 754 diagnostic; the scalar equation source is set
to zero for this calibration. This is not the actual 819 input source.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'820'))
import color_all_mode_certificate as fixed
old=fixed.previous;geo=old.geo
TARGET=HERE/'reference_normal_jet_probe_results.json'

def at_point(field,point):
    n=field.shape[0];wave=geo.waves(n)
    phase=np.exp(1j*np.einsum('...i,i->...',wave,point))
    return np.einsum('ijk...,ijk->...',np.fft.fftn(field,axes=(0,1,2))/n**3,phase).real

def run():
    data=old.setup();q=data['q'];psi=data['psi'];phi=q['phi'];h=phi[...,1]
    dp,_,_,_=old.old.inverse_gauss(q,data['k'],data['sigma'])
    v=psi[...,None]**-6*old.old.kinverse(q,q['p'])
    dv=psi[...,None]**-6*old.old.kinverse(q,dp)
    angular=[0,2,3]
    da_h=2*np.sum(v[...,angular]*dv[...,angular],axis=-1)/h
    # In this chart Gamma^A_BC=(delta^A_B phi_C+delta^A_C phi_B)/(6F).
    F=2-np.sum(phi*phi,axis=-1)/6
    connection_change=-(dv*np.sum(phi*v,axis=-1)[...,None]
                        +v*np.sum(phi*dv,axis=-1)[...,None])/(3*F[...,None])
    assert np.max(abs(connection_change[..., [1,4]]))<1e-13
    drh_dt=-2*v[...,1]*da_h
    assert np.max(abs(dv[..., [1,4]]))<1e-13
    assert np.max(abs(drh_dt))>1e-6
    point=np.array([0,np.pi/2,np.pi/4])
    pvalue=at_point(psi,point)
    # Original scalar reference spatial jets, sampled by the same periodic
    # spectral calibration. No unknown reference accelerations are set to zero.
    hgrad=np.stack([geo.derivative(h,i) for i in range(3)],axis=-1)
    s=phi[...,4];sgrad=np.stack([geo.derivative(s,i) for i in range(3)],axis=-1)
    rh=psi**-4*np.sum(hgrad*hgrad,axis=-1)-v[...,1]**2
    rs=psi**-4*np.sum(sgrad*sgrad,axis=-1)-v[...,4]**2
    rhgrad=np.array([at_point(geo.derivative(rh,i),point) for i in range(3)])
    rsgrad=np.array([at_point(geo.derivative(rs,i),point) for i in range(3)])
    block=np.array([rhgrad[[0,2]],rsgrad[[0,2]]])
    rhs=np.array([at_point(drh_dt,point),0.])
    coordinate_rate=np.linalg.solve(block,rhs)
    mixed_metric=-pvalue**4*coordinate_rate
    assert np.linalg.norm(block@coordinate_rate-rhs)<1e-12
    assert abs(np.linalg.det(block))>1e-10 and np.linalg.norm(mixed_metric)>1e-6
    # The Jacobian shear conclusion is exact in its algebraic variables.
    # The unspecified old rh_t/rs_t are arbitrary here; they cancel from this
    # proof-point determinant and are not claimed to solve the original PDE.
    vh=at_point(v[...,1],point);vs=at_point(v[...,4],point)
    sy=at_point(sgrad[...,1],point)
    determinants=[];shear_errors=[]
    for old_time in ((0.,0.),(.2,-.7),(1.3,2.1)):
        A=np.array([[vh,0,0,0],[vs,0,sy,0],
                    [old_time[0],*rhgrad],[old_time[1],*rsgrad]])
        delta=np.zeros((4,4));delta[2,0]=rhs[0]
        rate=np.linalg.solve(A,delta[:,0])
        shear_errors.append(float(np.linalg.norm(rate-np.array([0,coordinate_rate[0],0,coordinate_rate[1]]))))
        determinants.append(float(np.linalg.det(A+delta)-np.linalg.det(A)))
    assert max(shear_errors)<1e-12 and max(abs(x) for x in determinants)<1e-12
    return dict(round=824,formal_round_completed=False,all_working_checks_passed=True,
        original_grid=q['N'],original_Gauss_particular_used=True,
        radial_velocity_change_max=float(np.max(abs(dv[..., [1,4]]))),
        radial_connection_acceleration_change_max=float(np.max(abs(connection_change[..., [1,4]]))),
        higgs_modulus_acceleration_change_max=float(np.max(abs(da_h))),
        reference_normal_derivative_change_max=float(np.max(abs(drh_dt))),
        proof_point_spectral_reference_block=block.tolist(),
        proof_point_spectral_reference_normal_change=rhs.tolist(),
        reference_coordinate_rate_xz=coordinate_rate.tolist(),
        relational_metric_normal_xz_change=mixed_metric.tolist(),
        determinant_shear_errors=determinants,coordinate_shear_errors=shear_errors,
        source_is_754_diagnostic_with_zero_scalar_equation_source=True,
        actual819_input_source_not_computed=True,
        full_original_reference_normal_derivatives_not_set_to_zero=True,
        finite_menu_joint_state_not_whole_initial_slice=True,
        formal_test_groups_added=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
