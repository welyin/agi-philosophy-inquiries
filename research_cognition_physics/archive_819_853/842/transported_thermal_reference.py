"""842: original 64-component past mass, thermal sums and reset energy.

This explicitly changes the reference. It does not certify the unresolved
regional modular domain of the original pure reference or solve the future PDE.
"""
from pathlib import Path
import argparse,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'transported_thermal_reference_results.json'
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original

def spectral(a,f):
    w,u=np.linalg.eigh(a);return (u*f(w))@u.conj().T
def logistic(x):
    z=np.exp(-np.abs(x));return np.where(x>=0,1/(1+z),z/(1+z))
def entropy(a):
    w=np.linalg.eigvalsh(a);assert w.min()>-1e-12 and w.max()<1+1e-12
    w=np.clip(w,0,1);positive=w>0;below=w<1
    return float(-.5*(np.sum(w[positive]*np.log(w[positive]))+
                      np.sum((1-w[below])*np.log1p(-w[below]))))
def tails(beta,cut,massmax):
    # Sum n>=N of n^p exp(-beta*n), p=0..3. Exact geometric moments.
    n=cut+1;z=math.exp(-beta);a=1-z
    moments=(1/a,z/a**2,z*(1+z)/a**3,z*(1+4*z+z*z)/a**4)
    sums=[z**n*sum(math.comb(p,j)*n**(p-j)*moments[j] for j in range(p+1)) for p in range(4)]
    logz=32*(24*sums[2]+2*sums[0])
    energy=32*(math.sqrt(3)*(24*sums[3]+2*sums[1])+massmax*(24*sums[2]+2*sums[0]))
    return logz,energy

def run():
    positive=np.linalg.eigvalsh(original.MASS);positive=positive[positive>0]
    assert len(positive)==32
    gap=float(positive.min());massmax=float(positive.max());cut=12
    grid=np.arange(-cut,cut+1)
    ks=np.stack(np.meshgrid(grid,grid,grid,indexing='ij'),axis=-1).reshape(-1,3)
    energies=np.sqrt(np.sum(ks*ks,axis=1)[:,None]+positive[None,:]**2)
    thermal=[]
    for beta in (1.,2.,4.):
        n=logistic(-beta*energies)
        logz=float(np.logaddexp(0,-beta*energies).sum());energy=float((energies*n).sum())
        tz,te=tails(beta,cut,massmax)
        # Verify the closed-form tail dominates a large explicit shell sum.
        ns=np.arange(cut+1,cut+501,dtype=float)
        assert 32*np.sum((24*ns**2+2)*np.exp(-beta*ns))<=tz*(1+1e-13)
        thermal.append(dict(beta=beta,cube_half_width=cut,log_partition_partial=logz,
            log_partition_tail_bound=tz,excitation_energy_partial=energy,excitation_energy_tail_bound=te,
            entropy_partial=beta*energy+logz,entropy_tail_bound=beta*te+tz))
    q0=np.diag([1.]*10+[0.]*22+[1.]*10+[0.]*22)
    assert np.linalg.norm(original.CHARGE@q0@original.CHARGE-q0)==0
    assert np.linalg.norm(q0@original.MASS-original.MASS@q0)>1e-3
    rows=[]
    for k in ([0,0,0],[1,-1,0],[2,1,-1]):
        h=original.hamiltonian(np.array(k,float));q=q0.copy();charge=original.CHARGE
        if any(k):
            # C interchanges k and -k. A nonzero momentum fibre alone is not
            # a self-dual CAR system; retain both fibres before counting 1/2.
            zero=np.zeros((64,64),complex)
            h=np.block([[h,zero],[zero,original.hamiltonian(-np.array(k,float))]])
            q=np.block([[q0,zero],[zero,q0]])
            charge=np.block([[zero,original.CHARGE],[original.CHARGE,zero]])
        dim=len(h);rank=int(np.trace(q).real);r=np.eye(dim)-q
        assert np.linalg.norm(charge@h.conj()@charge+h)<1e-12
        abs_h=spectral(h,np.abs)
        eq=float(np.trace(q@abs_h).real)
        for beta in (1.,2.,4.):
            a=spectral(h,lambda x:logistic(beta*x));g=q/2+r@a@r
            de=float(-.5*np.trace(h@(g-a)).real)
            ds=entropy(g)-entropy(a);b=beta*de-ds
            assert -1e-11<=de<=.75*eq+1e-11
            assert -1e-11<=ds<=rank*np.log(2)+1e-11 and b>=-1e-11
            # Reference log pair is evaluated from h, not tiny covariance eigenvalues.
            loga=spectral(h,lambda x:-np.logaddexp(0,-beta*x))
            logia=spectral(h,lambda x:-np.logaddexp(0,beta*x))
            direct=float(-entropy(g)-.5*np.trace(g@loga+(np.eye(dim)-g)@logia).real)
            assert abs(direct-b)<1e-10
            lpast=spectral(h,lambda x:np.logaddexp(0,beta*x)+np.logaddexp(0,-beta*x))
            assert np.linalg.eigvalsh(beta*abs_h+2*np.log(2)*np.eye(dim)-lpast).min()>-1e-11
            rows.append(dict(momentum=k,includes_negative_momentum_partner=bool(any(k)),
                beta=beta,self_dual_code_dimension=rank,total_self_dual_dimension=dim,
                code_absolute_past_energy=eq,reset_excitation_energy_increment=de,
                reset_entropy_increment=ds,reference_relative_entropy=b,
                direct_relative_entropy_residual=abs(direct-b),global_relative_entropy_upper_bound=.75*beta*eq))
    # At k=0 this is a literal self-dual 32-CAR finite factor. Positive-gap
    # limit demonstrates why finite beta does not prove a uniform vacuum bound.
    h=original.MASS;q=q0;r=np.eye(64)-q
    p=original.I-original.occupied(np.zeros(3));g=q/2+r@p@r
    de0=float(-.5*np.trace(h@(g-p)).real)
    assert de0>=gap*(1-1/1024)-1e-12
    return dict(round=842,all_checks_passed=True,fresh_test_groups=1,
        original_past_positive_excitation_modes_per_momentum=32,original_auxiliary_gap=gap,
        original_auxiliary_largest_mass=massmax,thermal_T3_resource_sums=thermal,
        original_mass_reset_calibrations=rows,zero_temperature_reset_energy_increment_at_zero_momentum=de0,
        diagnostic_code_is_actual_ten_spatial_sterile_modes=False,
        reference_change_explicit=True,original_pure_reference_domain_problem_solved=False,
        thermal_reference_Hadamard_and_finite_entropy_are_analytic=True,
        future_reference_is_instantaneous_KMS=False,
        auxiliary_energy_is_current_physical_H=False,
        native_thermal_preparation_or_heat_bath_proven=False,
        original_background_PDE_solved=False,backreaction_solution_computed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
