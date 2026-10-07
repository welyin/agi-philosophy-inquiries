"""836: equal configuration law plus minimal kinetic energy fixes the state.

Finite positive-kinetic calibration and original configuration metric checks.
The continuum Gauss statement is proved in the note, not by the graph grid.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'configuration_energy_rigidity_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime


def laplacian(size,edges):
    out=np.zeros((size,size))
    for i,j,w in edges:
        out[i,i]+=w;out[j,j]+=w;out[i,j]-=w;out[j,i]-=w
    return out


def graph_checks():
    rng=np.random.default_rng(836);size=7
    psi=np.arange(1,size+1,dtype=float);psi/=np.linalg.norm(psi)
    edges=[(i,i+1,1+.1*i) for i in range(size-1)]+[(0,size-1,.4)]
    t=laplacian(size,edges);potential=np.linspace(.2,2.,size)**2
    h=t+np.diag(potential)
    u=rng.normal(size=(size,3))+1j*rng.normal(size=(size,3))
    u/=np.linalg.norm(u,axis=1)[:,None]
    phi=psi[:,None]*u;sigma=phi@phi.conj().T;rho=np.outer(psi,psi)
    assert min(np.linalg.eigvalsh(sigma))>-1e-14
    assert np.max(abs(np.diag(sigma)-psi**2))<1e-15
    difference=float(np.trace(h@(sigma-rho)).real)
    expected=sum(w*psi[i]*psi[j]*np.linalg.norm(u[i]-u[j])**2 for i,j,w in edges)
    assert abs(difference-expected)<1e-13 and difference>0
    constant=np.tile(np.array([1,1j,2])/np.sqrt(6),(size,1))
    cphi=psi[:,None]*constant;crho=cphi@cphi.conj().T
    assert np.linalg.norm(crho-rho)<1e-14
    # Without connected support, independent component phases/coherences cost
    # no kinetic energy, and equality no longer forces the original pure state.
    de=[(0,1,1.),(1,2,1.),(3,4,1.),(4,5,1.),(5,6,1.)]
    vectors=np.zeros((size,2));vectors[:3,0]=1;vectors[3:,1]=1
    dphi=psi[:,None]*vectors;drho=dphi@dphi.T
    dh=laplacian(size,de)+np.diag(potential)
    disconnected_energy=float(np.trace(dh@(drho-rho)))
    assert abs(disconnected_energy)<1e-14 and np.linalg.norm(drho-rho)>.1
    # Without a phase-free initial state, time-reversed phase patterns have
    # the same configuration law and scalar energy but are distinct states.
    phased=psi*np.exp(1j*np.linspace(0,1.3,size))
    reversed_phase=phased.conj()
    phase_energy=float(np.vdot(phased,h@phased).real)
    reverse_energy=float(np.vdot(reversed_phase,h@reversed_phase).real)
    assert abs(phase_energy-reverse_energy)<1e-14
    phase_distance=float(np.linalg.norm(np.outer(phased,phased.conj())-
                                        np.outer(reversed_phase,reversed_phase.conj())))
    assert phase_distance>.1
    return dict(configuration_diagonal_residual=float(np.max(abs(np.diag(sigma)-psi**2))),
                positive_energy_increment=difference,Dirichlet_identity_residual=abs(difference-expected),
                constant_vector_gives_original_pure_state=True,
                disconnected_equal_energy_counterexample=disconnected_energy,
                disconnected_state_difference=float(np.linalg.norm(drho-rho)),
                phased_initial_equal_energy_state_difference=phase_distance)


def native_metric_checks(old):
    rng=np.random.default_rng(8361);maximum=0.;minimum=1e100
    for _ in range(24):
        phi=rng.normal(size=5)*.2
        inv=old.geom.inverse(phi)
        assert min(np.linalg.eigvalsh(inv))>0
        # Two spectral components of a mixed state, expressed as psi(q) u(q).
        a=rng.normal(size=5);b=rng.normal(size=5);c=rng.normal(size=5)
        theta=.4+a@phi;alpha=b@phi;beta=c@phi
        u=np.array([np.cos(theta)*np.exp(1j*alpha),np.sin(theta)*np.exp(1j*beta)])
        du=np.vstack([np.exp(1j*alpha)*(-np.sin(theta)*a+1j*np.cos(theta)*b),
                      np.exp(1j*beta)*(np.cos(theta)*a+1j*np.sin(theta)*c)])
        psi=np.exp(-float(phi@phi));dpsi=-2*phi*psi
        dcomponents=u[:,None]*dpsi[None,:]+psi*du
        before=float(dpsi@inv@dpsi)
        after=float(sum(np.vdot(v,inv@v).real for v in dcomponents))
        extra=psi**2*float(sum(np.vdot(v,inv@v).real for v in du))
        err=abs(after-before-extra);maximum=max(maximum,err);minimum=min(minimum,extra)
        assert err<2e-12 and extra>=0
    return dict(original_curved_inverse_metric_points=24,
                mixed_component_gradient_identity_max_residual=maximum,
                minimum_positive_extra_term=minimum,
                values_are_local_identity_checks_not_a_native_global_time_simulation=True)


def run():
    with ResearchRuntime(Layout()).installed():
        import joint_record_mass_feedback as old
        native=native_metric_checks(old)
    return dict(round=836,all_checks_passed=True,fresh_test_groups=1,
                finite_configuration_graph=graph_checks(),native_curved_metric=native,
                connected_phase_free_support_is_required=True,
                full_configuration_statistics_is_stronger_than_a_finite_menu=True,
                mixed_states_are_included_in_analytic_rigidity=True,
                other_fermionic_material_changes_are_not_excluded=True,
                original_autonomous_compensation_or_all_physics_complete=False,
                scope='For a fixed scalar Hamiltonian with positive kinetic metric, a phase-free wavefunction positive on connected support minimizes kinetic energy among all density matrices with the same configuration law; equality fixes that pure state. Applying this to the 835 frozen-fermion endpoint rules out a scalar-only compensation preserving every configuration statistic and mean energy. Finite menus, disconnected supports, phased initial states and changed fermionic resources are not excluded.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
