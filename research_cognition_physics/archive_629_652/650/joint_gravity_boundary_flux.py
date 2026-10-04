"""650: original scalar modes, constraint-compatible data and boundary flux.

Einstein action and 3+1 dimensions are inputs. The asymptotically-flat
completion is proved by contraction in the note, not solved on a finite grid.
"""
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_geometric_boundary_matching as boundary

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gravity_boundary_flux_results.json'


def modes():
    mat,u,_=original.lattice.scalar.parameters()
    phi=boundary.vacuum();F=float(original.F(phi));K=original.metric(phi)
    idx=[1,4];kr=K[np.ix_(idx,idx)]
    Hess=2*np.diag(np.sqrt(u))@mat@np.diag(np.sqrt(u))/F**2
    d,V=np.linalg.eigh(kr);invsqrt=(V*d**-.5)@V.T
    mass2,v=np.linalg.eigh(invsqrt@Hess@invsqrt)
    basis=np.zeros((5,2));basis[idx,:]=invsqrt@v
    return phi,K,Hess,mass2,basis


def original_modes_and_frame_flux():
    spec=importlib.util.spec_from_file_location('boundary_probe650',HERE/'round650_drafts/original_boundary_flux_probe.py')
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    old=probe.run()
    assert old==json.loads((HERE/'round650_drafts/original_boundary_flux_probe_results.json').read_text('utf8'))
    phi,K,Hess,mass2,basis=modes();idx=[1,4]
    assert np.max(abs(basis.T@K@basis-np.eye(2)))<2e-14
    assert np.max(abs(Hess@basis[idx,:]-K[np.ix_(idx,idx)]@basis[idx,:]*mass2))<2e-14
    rows=[]
    for step in (.004,.002,.001):
        got=np.zeros((2,2))
        for i in range(2):
            di=np.eye(5)[idx[i]]*step
            got[i,i]=(original.node_potential(phi+di)+original.node_potential(phi-di)-2*original.node_potential(phi))/step**2
            for j in range(i):
                dj=np.eye(5)[idx[j]]*step
                got[i,j]=got[j,i]=(original.node_potential(phi+di+dj)-original.node_potential(phi+di-dj)-original.node_potential(phi-di+dj)+original.node_potential(phi-di-dj))/(4*step**2)
        rows.append(dict(step=step,original_potential_Hessian_error=float(np.max(abs(got-Hess)))))
    assert rows[-1]['original_potential_Hessian_error']<rows[0]['original_potential_Hessian_error']/12
    return dict(original_probe_reproduced=True,original_radial_mass_squared=mass2.tolist(),
        independent_potential_differences=rows,frame_flux_rows=old['rows'],
        plane_wave_time_dependence_not_assigned_to_cutoff_AF_solution=True)


def AF_completion_certificate():
    phi,K,Hess,mass2,basis=modes();mat,_,_=original.lattice.scalar.parameters()
    # Conservative strict margins around the inherited finite input values.
    assert np.linalg.norm(phi)<.9 and np.linalg.norm(basis,axis=0).max()<1.5
    assert 0<np.linalg.eigvalsh(mat).min() and np.linalg.eigvalsh(mat).max()<1
    R=Q(1);k=Q(7,10);eta=Q(1,100);anorm=Q(3,2);pstar=Q(9,10)
    # |grad chi| <= 8/(e R) < 3/R, chi=exp[-r^2/(R^2-r^2)].
    displacement=anorm*eta;Fmin=Q(2)-(pstar+displacement)**2/6
    D=k+3/R;Kmax=2/Fmin**2
    Bmax=Kmax*displacement**2*D**2
    Delta=2*pstar*displacement+displacement**2
    Umax=Delta**2/(4*Fmin**2)
    C=R**2/16;upper=Q(3,2)
    image=C*(Bmax*upper+2*Umax*upper**5)
    contraction=C*(Bmax+10*Umax*upper**4)
    assert Fmin>1 and image<Q(1,2) and contraction<1
    # A numerical spot-check of the original nonlinear B,U; the proof of
    # the global bounds above does not use this sample as a supremum.
    rng=np.random.default_rng(650);x=rng.uniform(-1,1,(180,3));x=x[np.sum(x*x,axis=1)<1]
    r2=np.sum(x*x,axis=1);chi=np.exp(-r2/(1-r2))
    gradchi=-2*x*chi[:,None]/(1-r2[:,None])**2
    eps=.006;delta=-.004;a=basis[:,1]
    osc=eps*np.cos(float(k)*x[:,0])+delta*np.sin(float(k)*x[:,0])
    grad=gradchi*osc[:,None]
    grad[:,0]+=chi*float(k)*(-eps*np.sin(float(k)*x[:,0])+delta*np.cos(float(k)*x[:,0]))
    fields=phi+chi[:,None]*osc[:,None]*a
    B=np.sum(grad*grad,axis=1)*np.einsum('i,nij,j->n',a,original.metric(fields),a)
    U=original.node_potential(fields)
    assert np.max(B)<float(Bmax) and np.max(U)<float(Umax) and np.min(original.F(fields))>float(Fmin)
    rational=lambda v:dict(numerator=v.numerator,denominator=v.denominator,decimal=float(v))
    return dict(parameter_diamond_abs_eps_plus_abs_delta=float(eta),support_radius=float(R),
        F_lower_bound=rational(Fmin),B_upper_bound=rational(Bmax),U_upper_bound=rational(Umax),
        Newton_image_bound=rational(image),Newton_Lipschitz_bound=rational(contraction),
        interval_for_conformal_factor=[1.,1.5],
        sample_B_max=float(np.max(B)),sample_U_max=float(np.max(U)),
        sample_count=len(x),sample_not_used_as_proof_of_uniform_bound=True,
        nonlinear_AF_constraint_solution_existence_analytic=True,
        numerical_solution_of_full_PDE_not_claimed=True,
        compact_torus_linearization_stability_not_assumed=True)


def integrated_flux_and_gluing():
    # At t=0 the cutoff data have delta1 phi=a chi and
    # normal(delta2 phi)=a k chi on x=0, hence omega_n=k chi^2.
    # Entire plane integral is compactly supported; the artificial interface
    # has a second side with opposite outward normal, not a reflecting wall.
    k=.7;rows=[]
    for order in (64,128,192):
        x,w=np.polynomial.legendre.leggauss(order);r=(x+1)/2;w=w/2
        chi=np.exp(-r*r/(1-r*r))
        one=float(2*np.pi*k*np.sum(w*r*chi**2))
        rows.append(dict(order=order,left_outward_flux=one,right_outward_flux=-one,
            two_sided_sum=0.,nonintegrability_curl_magnitude=one))
    assert abs(rows[-1]['left_outward_flux']-rows[-2]['left_outward_flux'])<2e-13
    assert rows[-1]['left_outward_flux']>.5
    # Full original boundary one-form at arbitrary off-shell jets, with
    # field-dependent F and metric, checked in both frames. Differentiate
    # it independently rather than replacing it with the vacuum KG formula.
    phi0=boundary.vacuum()+np.array([.02,-.03,.01,.02,.04])
    h0=np.diag([-1.2,1.1,.9]);K0=np.array([[.02,.01,0],[.01,-.03,.01],[0,.01,.01]])
    v0=np.array([.03,-.02,.01,.02,-.04]);rng=np.random.default_rng(6501)
    dirs=[]
    for _ in range(2):
        dh=rng.normal(size=(3,3))*.1;dh=(dh+dh.T)/2
        dK=rng.normal(size=(3,3))*.1;dK=(dK+dK.T)/2
        dirs.append((rng.normal(size=5)*.1,dh,dK,rng.normal(size=5)*.1))
    base=(phi0,h0,K0,v0)
    def oneform(point,tangent,frame,side):
        phi,h,K,v=point;dphi,dh,_,_=tangent
        if frame=='J':
            P,j=boundary.jordan_momenta(phi,h,side*K,side*v)
            return np.sqrt(abs(np.linalg.det(h)))*(.5*np.sum(P*dh)+j@dphi)
        hE,KE,vE,PE,jE=boundary.einstein_data(phi,h,side*K,side*v)
        F,f=boundary.data(phi);dhE=F*dh+(f@dphi)*h
        return np.sqrt(abs(np.linalg.det(hE)))*(.5*np.sum(PE*dhE)+jE@dphi)
    def curl(frame,side,step):
        val=0.
        for i,j,sgn in ((0,1,1),(1,0,-1)):
            plus=tuple(b+step*d for b,d in zip(base,dirs[i]))
            minus=tuple(b-step*d for b,d in zip(base,dirs[i]))
            val+=sgn*(oneform(plus,dirs[j],frame,side)-oneform(minus,dirs[j],frame,side))/(2*step)
        return float(val)
    pairs=[]
    for step in (.002,.001,.0005):
        J=curl('J',1,step);E=curl('E',1,step);other=curl('J',-1,step)
        assert abs(J-E)<2e-10 and abs(J+other)<1e-12 and abs(J)>1e-3
        pairs.append(dict(step=step,Jordan_curl=J,Einstein_curl=E,
            frame_difference=abs(J-E),matched_two_sided_curl=J+other))
    return dict(compact_interface_integrals=rows,full_boundary_pairing=pairs,
        general_boundary_jet_check_is_kinematic_only=True,
        AF_family_supplies_separate_on_shell_witness=True,
        two_sided_flux_cancellation_not_single_region_autonomy=True,
        no_assertion_against_boundary_edge_extensions_or_flux_charges=True)


def run():
    deps=('research_note_573.md','research_note_600.md','research_note_618.md',
          'research_note_619.md','research_note_647.md','research_note_648.md','research_note_649.md',
          'joint_curved_quantum_source.py','joint_geometric_boundary_matching.py',
          'round650_drafts/original_boundary_flux_probe.py',
          'round650_drafts/original_boundary_flux_probe_results.json')
    return dict(round=650,tests_run=3,failures=0,errors=0,
        original_modes_and_frame_flux=original_modes_and_frame_flux(),
        AF_completion=AF_completion_certificate(),boundary_flux=integrated_flux_and_gluing(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(given_3plus1_leading_Einstein_five_scalar_action=True,
            exact_small_AF_initial_data_family_and_short_time_classical_development=True,
            fixed_timelike_interface_only=True,first_variation_area_zero_not_exact_fixed_area_family=True,
            matching_between_two_sides_exact_for_restrictions_of_same_global_solution=True,
            no_full_quantum_gravity_or_relational_instrument_completion=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
