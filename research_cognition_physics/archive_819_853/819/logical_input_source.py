"""819 working: local sources of the logical input family.
Original 64-component auxiliary-past principal and mass matrices are retained.
Two compact separated spatial packets calibrate a source difference; no
future background, constraint solution or preparation device is simulated.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
sys.path.insert(0,str(HERE.parent/'818'))
import preparation_noise_bridge as bridge
TARGET=HERE/'logical_input_source_results.json'

def compact(x,center,width):
    y=(x-center)/width;chi=np.zeros_like(x);derivative=np.zeros_like(x)
    mask=abs(y)<1
    chi[mask]=np.exp(-1/(1-y[mask]**2))
    derivative[mask]=chi[mask]*(-2*y[mask]/(1-y[mask]**2)**2)/width
    return chi,derivative

def covariance(rho,fields):
    return np.array([[np.trace(rho@a@b.conj().T) for b in fields] for a in fields])

def run():
    x,w=np.polynomial.legendre.leggauss(256)
    c1,d1=compact(x,-.55,.3);c2,d2=compact(x,.55,.3);cz,dz=compact(x,0.,.8)
    for c,d in ((c1,d1),(c2,d2),(cz,dz)):
        norm=np.sqrt(np.dot(w,c*c));c/=norm;d/=norm
    assert np.max(abs(c1*c2))+np.max(abs(d1*c2))+np.max(abs(c1*d2))==0
    s=np.zeros(64,complex);spin=original.GAMMA[2][np.ix_([30,31],[30,31])]
    es,vs=np.linalg.eigh(spin);s[[30,31]]=vs[:,-1]
    gamma=float(np.vdot(s,original.GAMMA[2]@s).real)
    mass=float(np.vdot(s,original.MASS@s).real)
    assert abs(gamma-1)<1e-12
    mass_action=float(np.linalg.norm(original.MASS@s))
    assert mass_action>1e-3  # mass coupling is retained, even if its diagonal vanishes.
    a=bridge.prep.previous.old.car(2);fields=a+[c.conj().T for c in a]
    frame=bridge.record_frame(a)
    def lifted(rho):return frame@np.kron(rho,np.eye(2)/2)@frame.conj().T
    rZp=lifted((bridge.I2+bridge.Z)/2);rZm=lifted((bridge.I2-bridge.Z)/2)
    rXp=lifted((bridge.I2+bridge.X)/2);rXm=lifted((bridge.I2-bridge.X)/2)
    delta_z=covariance(rZp,fields)-covariance(rZm,fields)
    delta_x=covariance(rXp,fields)-covariance(rXm,fields)
    expected=np.diag([1.,0.,-1.,0.])
    assert np.max(abs(delta_z-expected))<1e-12
    assert np.linalg.norm(delta_x)>1.
    number=[c.conj().T@c for c in a]
    records=[]
    for kappa in (4.,8.,16.,32.):
        f=cz[:,None]*np.exp(1j*kappa*x[:,None])*s
        df=(dz+1j*kappa*cz)[:,None]*np.exp(1j*kappa*x[:,None])*s
        current=np.einsum('ti,ij,tj->t',f.conj(),original.GAMMA[2],df)
        kinetic=np.dot(w,current.imag)
        full=float((kinetic+mass).real)
        assert abs(full-kappa)<1e-11
        source=full*(number[0]+number[1])
        z_contrast=float(np.trace((rZp-rZm)@source).real)
        x_contrast=float(np.trace((rXp-rXm)@source).real)
        assert abs(z_contrast+kappa)<1e-11 and abs(x_contrast)<1e-12
        # Equal local quadratic source for X+ and X- follows pointwise:
        # their covariance difference has only cross-packet legs.
        diagonal_x=float(np.max(abs(np.diag(delta_x))))
        assert diagonal_x<1e-12
        records.append(dict(kappa=kappa,packet_energy=full,
            Z_plus_minus_energy_source_difference=z_contrast,
            X_plus_minus_energy_source_difference=x_contrast))
    mixture=.23*rZp+.77*rXp
    J=4*(number[0]+number[1])
    affine_error=float(abs(np.trace(mixture@J)-.23*np.trace(rZp@J)-.77*np.trace(rXp@J)))
    assert affine_error<1e-12
    return dict(round=819,status='working_not_formal',all_checks_passed=True,
        compact_packet_cross_support_max=float(np.max(abs(c1*c2))),
        original_64_component_mass_action_norm=mass_action,
        mass_diagonal_expectation=mass,
        logical_Z_covariance_difference_error=float(np.max(abs(delta_z-expected))),
        logical_X_covariance_difference_norm=float(np.linalg.norm(delta_x)),
        source_rows=records,source_affinity_error=affine_error,
        original_curved_Cauchy_source_nonzero_lift_proven=False,
        full_constraint_compensation_computed=False,
        gravity_alone_forced_to_change=False,
        formal_numbered_test_groups_added=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))

