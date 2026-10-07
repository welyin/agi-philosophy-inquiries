"""Working900: actual859 constraint data and h-slice spatial material Jacobian.
This evaluates the original initial data, not a new oscillator. Collocation and
finite differences are not a continuum/interval certificate or time evolution.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'current_reference_initial_budget_results.json'
def point_value(coeff,k,pos):
    phase=np.exp(1j*(k[...,0]*pos[0]+k[...,1]*pos[1]+k[...,2]*pos[2]))
    return np.real(np.einsum('ijk,ijk...->...',phase,coeff))
def run():
    rows=[]
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        geo=old.geo
        for N in (16,24,32):
            q,psi0,At0,_,color=old.completed(N,1.)
            epsilon=.02;xx=q['grid'][...,0];probe=epsilon*np.sin(xx)
            new=dict(q);new['B']=q['B']+epsilon**2*np.cos(xx)**2
            new['U']=q['U']+.5*probe**2;new['C']=q['C']-probe**2
            psi,At,stats=geo.solve_hamiltonian(new,initial=psi0)
            idx=(0,N//4,N//8);pos=q['grid'][idx];ps=float(psi[idx])
            hstar=np.sqrt(geo.old.PAR['h2']);_,u,_=geo.old.scalar.parameters();sstar=np.sqrt(u[1])
            phi=q['phi'][idx];pi=q['p'][idx];FF=float(q['F'][idx])
            v=ps**-6*FF*(pi-phi*np.dot(phi,pi)/12);vh=float(v[1]);assert vh>0
            E=color['E'];Fc=color['F'];S=sum(old.inner(Fc[i,j],Fc[i,j]) for i in range(3) for j in range(i+1,3))
            flux=np.array([sum(old.inner(E[j],Fc[i,j]) for j in range(3)) for i in range(3)])
            assert np.max(abs(flux))<1e-16
            gradpsi=np.array([geo.derivative(psi,j)[idx] for j in range(3)])
            gradM=-8*S*ps**-9*gradpsi
            J=np.vstack(([0.,-.06*sstar,0.],[epsilon,0.,0.],gradM))
            det=float(np.linalg.det(J));assert det>0
            expected=vh*det
            k=geo.waves(N);psihat=np.fft.fftn(psi)/N**3;pihat=np.fft.fftn(q['p'],axes=(0,1,2))/N**3
            def actual_magnetic(point):
                pp=float(point_value(psihat,k,point));pvec=point_value(pihat,k,point)
                h=hstar*(1+.05*np.sin(point[0]+point[1]));s=sstar*(1+.06*np.cos(point[1]))
                ph=np.array([0.,h,0.,0.,s]);fff=2-np.dot(ph,ph)/6
                hv=float((pp**-6*fff*(pvec-ph*np.dot(ph,pvec)/12))[1])
                dh=np.array([.05*hstar*np.cos(point[0]+point[1])]*2+[0.])
                normal_sq=hv**2-pp**-4*np.dot(dh,dh);assert normal_sq>0
                n=np.r_[hv,-pp**-4*dh]/np.sqrt(normal_sq)
                invg=np.diag([-1.,pp**-4,pp**-4,pp**-4]);P=invg+np.outer(n,n)
                spacetimeF=np.zeros((4,4,3,3),complex);spacetimeF[1:,1:]=Fc
                bc=geo.old.PAR['b'][0];velocity=2*bc*pp**-2*E
                spacetimeF[0,1:]=velocity;spacetimeF[1:,0]=-velocity
                gram=2*np.einsum('mnij,rsji->mnrs',spacetimeF,spacetimeF).real
                return float(.5*np.einsum('mr,ns,mnrs',P,P,gram))
            center=actual_magnetic(pos);assert abs(center-S*ps**-8)<1e-16
            difference_rows=[]
            for step in (1e-3,5e-4,2.5e-4,1e-5,5e-6,2.5e-6):
                measured=np.array([(actual_magnetic(pos+step*np.eye(3)[j])-actual_magnetic(pos-step*np.eye(3)[j]))/(2*step) for j in range(3)])
                difference_rows.append(dict(step=step,full_covariant_spatial_derivative=measured.tolist(),maximum_formula_error=float(np.max(abs(measured-gradM)))))
            assert difference_rows[-1]['maximum_formula_error']<2e-8,(N,gradM.tolist(),difference_rows)
            # Coordinate scales are declared; inverse magnitudes alone are not
            # invariant resource costs or evidence for a physical minimum size.
            scales=np.array([sstar,epsilon,S]);Jn=J/scales[:,None]
            inv=np.linalg.inv(J);invn=np.linalg.inv(Jn)
            rows.append(dict(N=N,epsilon=epsilon,point=pos.tolist(),constraint_residual=float(stats['original_Hamiltonian_residual']),
                psi=ps,clock_velocity=vh,color_flux_max=float(np.max(abs(flux))),magnetic_scalar=center,
                actual_spatial_material_J=J.tolist(),spatial_determinant=det,spacetime_determinant=expected,
                inverse_norm2=float(np.linalg.norm(inv,2)),minimum_singular_value=float(np.linalg.svd(J,compute_uv=False)[-1]),
                declared_material_scales=scales.tolist(),normalized_inverse_norm2=float(np.linalg.norm(invn,2)),
                normalized_minimum_singular_value=float(np.linalg.svd(Jn,compute_uv=False)[-1]),
                identity_error=float(np.max(abs(J@inv-np.eye(3)))),full_covariant_derivative_crosschecks=difference_rows))
    oldrows=json.loads((HERE.parent/'859/probe_reference_cone_compatibility_results.json').read_text('utf-8'))['background_rows']
    for row in rows:
        known=[x for x in oldrows if x['N']==row['N'] and x['probe_amplitude']==.02]
        if known:assert abs(row['spacetime_determinant']-known[0]['new_reference_jacobian'])<1e-16
    return dict(working_round=900,formal_round=899,cumulative_numbered_groups=3684,fresh_numbered_groups=0,
        status='working actual initial-data budget; not a completed round',rows=rows,
        original859_initial_constraint_model=True,full_normal_projector_and_color_electric_field_retained=True,
        full_spatial_Jacobian_at_reference_point_evaluated=True,old_arbitrary_determinant_auxiliary_entries_used=False,
        uniform_patch_or_time_error_certificate=False,original_quantum_modes_or_source_M_T_computed=False,
        floating_collocation_not_interval_continuum_proof=True,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
