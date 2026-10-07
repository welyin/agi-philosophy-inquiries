"""Working 859: one existing probe and original color magnetic reference.
Reuses 572/753 initial-constraint solver through the current read-only loader.
No new formal round or quantum-background equivalence is claimed here.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'shared_probe_reference_probe_results.json'
def run():
    runtime=ResearchRuntime(Layout())
    with runtime.installed():
        import joint_reference_constraint_strata as old
        geo=old.geo
        rows=[]
        for N,epsilon in ((16,.01),(16,.02),(24,.02)):
            q,psi0,tensor0,_,color=old.completed(N,1.)
            x,y,z=np.moveaxis(q['grid'],-1,0)
            m=1.;probe=epsilon*np.sin(x)
            new=dict(q);new['B']=q['B']+epsilon**2*np.cos(x)**2
            new['U']=q['U']+.5*m*m*probe**2
            new['C']=q['C']-m*m*probe**2
            assert float(np.min(new['C']))>1.99
            psi,tensor,stats=geo.solve_hamiltonian(new,initial=psi0)
            assert stats['original_Hamiltonian_residual']<3e-8
            assert np.max(abs(tensor-tensor0))<1e-14
            assert float(np.min(psi-psi0))>0
            psiz=geo.derivative(psi,2)
            # On x=0,y=pi/2 the original h spatial gradient and its z derivative
            # vanish, so the h-normal projector and its z derivative agree with
            # the initial spatial metric for this particular derivative.
            point=(0,N//4,N//8)
            assert psiz[point]<0
            h=float(q['phi'][point][1]);s=float(q['phi'][point][4]);F0=float(q['F'][point])
            vel=psi[point]**-6*F0*(q['p'][point]-q['phi'][point]*np.dot(q['phi'][point],q['p'][point])/12)
            vh=float(vel[1]);assert vh>0
            _,u,_=geo.old.scalar.parameters();br=.06*np.sqrt(u[1])
            coefficient=sum(old.inner(color['F'][i,j],color['F'][i,j]) for i in range(3) for j in range(i+1,3))
            assert coefficient>0
            mag=coefficient*psi**-8
            magz=-8*coefficient*psi[point]**-9*psiz[point]
            # Time derivative of the last row cancels: h has only a time entry,
            # s only t/y, the existing probe only x at this point.
            J=np.array([[vh,0,0,0],[float(vel[4]),0,-br,0],[0,epsilon,0,0],[.137,.211,-.091,magz]])
            jac=float(np.linalg.det(J));formula=vh*br*epsilon*magz
            assert jac>0 and abs(jac-formula)<1e-18
            magnetic_derivative_error=abs(geo.derivative(mag,2)[point]-magz)
            assert magnetic_derivative_error<2e-8
            direct_rho=(.5*psi**-12*q['pKp']+.5*psi**-4*new['B']+new['U']+psi**-8*q['Y'])
            einstein=(-8*psi**-5*geo.laplace(psi)-psi**-12*np.sum(tensor*tensor,axis=(-1,-2))+2*q['tau2']/3-2*direct_rho)
            assert float(np.max(abs(einstein)))<3e-8
            rows.append(dict(N=N,probe_amplitude=epsilon,probe_mass=m,
                original_constraint_residual=float(np.max(abs(einstein))),
                added_probe_momentum_density='0',unchanged_momentum_solution_error=float(np.max(abs(tensor-tensor0))),
                minimum_conformal_response=float(np.min(psi-psi0)),
                min_response_divided_by_amplitude_squared=float(np.min(psi-psi0)/epsilon**2),
                point_psi=float(psi[point]),point_psi_z=float(psiz[point]),point_clock_velocity=vh,
                color_magnetic_coefficient=float(coefficient),point_magnetic_scalar=float(mag[point]),
                point_magnetic_z_derivative=float(magz),magnetic_derivative_crosscheck_error=float(magnetic_derivative_error),
                new_reference_jacobian=jac,jacobian_divided_by_probe_amplitude=jac/epsilon))
    return dict(kind='round_859_working_probe',formal_reports=858,newly_published_numbered_groups=0,
        all_checks_passed=True,rows=rows,probe_species_added_this_probe=0,
        existing_probe_species_in_852=1,source_and_background_changed=True,
        numerical_collocation_not_analytic_existence_proof=True,
        quantum_state_and_reference_transport_completed=False,
        canonical_configuration_only_property_on_reduced_physical_phase_space_proved=False,
        self_consistent_graph_to_continuum_matching_proved=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
