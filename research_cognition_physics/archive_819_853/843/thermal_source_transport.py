"""843: reference, record reset and source derivatives in original coefficients.

The full original Nambu matrices are kept. Time evolution here is the old
730 local-symbol diagnostic, NOT the original inhomogeneous Dirac PDE.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'thermal_source_transport_results.json'
sys.path.insert(0,str(HERE.parent/'842'))
from transported_thermal_reference import original,logistic,spectral
vertex=original.vertex;old=vertex.original

def block(a,b):
    z=np.zeros_like(a);return np.block([[a,z],[z,b]])
def run():
    width=.08;start=-.24;steps=96;dt=-start/steps;k=np.array([1.,-1.,0.])
    past=[];flows=[];ends=[];gammas=[]
    for sign in (1,-1):
        h0=old.matrices(start,width,sign*k)[0];u=np.eye(64,dtype=complex)
        for j in range(steps):
            h=old.matrices(start+(j+.5)*dt,width,sign*k)[0]
            step,_=old.step(h,np.zeros_like(h),dt);u=step@u
        h,g=old.matrices(0.,width,sign*k)
        past.append(h0);flows.append(u);ends.append(h);gammas.append(g)
    h0=block(*past);u=block(*flows);end=block(*ends)
    charge=np.block([[np.zeros((64,64)),original.CHARGE],[original.CHARGE,np.zeros((64,64))]])
    assert np.linalg.norm(charge@h0.conj()@charge+h0)<1e-11
    assert np.linalg.norm(charge@u.conj()@charge-u)<1e-11
    q0=np.diag([1.]*10+[0.]*22+[1.]*10+[0.]*22);q=block(q0,q0);r=np.eye(128)-q
    phi=old.collar(0.,width)[1]
    labels=['energy','spatial_scale']+[f'scalar_{j}' for j in range(5)]
    sources=[end,block(*gammas)]+[block(vertex.dmass(phi,np.eye(5)[j]),vertex.dmass(phi,np.eye(5)[j])) for j in range(5)]
    vacuum=spectral(h0,lambda x:(x>0).astype(float))
    pure_now=u@vacuum@u.conj().T
    def data(beta):
        a=spectral(h0,lambda x:logistic(beta*x))
        db=spectral(h0,lambda x:x*logistic(beta*x)*logistic(-beta*x))
        delta=u@(a-vacuum)@u.conj().T
        reset_delta=r@delta@r;derivative=r@u@db@u.conj().T@r
        values=np.array([-.5*np.trace(g@reset_delta).real for g in sources])
        slopes=np.array([-.5*np.trace(g@derivative).real for g in sources])
        return a,delta,reset_delta,values,slopes
    rows=[]
    for beta in (1.,2.,4.):
        a,delta,d,values,slopes=data(beta)
        thermal_now=u@a@u.conj().T
        pure_post=q/2+r@pure_now@r;thermal_post=q/2+r@thermal_now@r
        residual=float(np.linalg.norm(thermal_post-pure_post-d))
        reality=float(np.linalg.norm(charge@d.conj()@charge+d))
        assert max(residual,reality)<1e-11
        step=2e-4
        numeric=(data(beta+step)[3]-data(beta-step)[3])/(2*step)
        error=float(np.max(abs(numeric-slopes)));assert error<2e-7
        omitted=np.array([-.5*np.trace(g@delta).real for g in sources])
        omission=float(np.max(abs(omitted-values)));assert omission>1e-3
        inst=spectral(end,lambda x:logistic(beta*x))
        instantaneous_error=float(np.linalg.norm(thermal_now-inst));assert instantaneous_error>1e-3
        assert max(abs(values[2:]))>1e-3
        rows.append(dict(beta=beta,organized_source_change=dict(zip(labels,values.tolist())),
            beta_source_derivative=dict(zip(labels,slopes.tolist())),
            source_derivative_residual=error,reset_difference_identity_residual=residual,
            self_dual_state_difference_residual=reality,
            omitted_record_projection_max_source_error=omission,
            wrong_instantaneous_thermal_reference_norm_error=instantaneous_error))
    return dict(round=843,all_checks_passed=True,fresh_test_groups=1,
        original_Nambu_pair_dimension=128,self_dual_diagnostic_reset_rank=int(np.trace(q).real),
        full_original_mass_gauge_kinetic_coefficients_retained=True,rows=rows,
        original_ten_spatial_code_modes_simulated=False,original_inhomogeneous_PDE_solved=False,
        source_Ward_and_first_order_completion_are_analytic=True,
        original_continuum_source_values_computed=False,
        reference_and_record_source_changes_both_retained=True,
        full_quantum_gravity_or_native_preparation_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
