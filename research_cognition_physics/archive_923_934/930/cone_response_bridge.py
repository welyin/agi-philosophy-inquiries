"""930: conditional metric constitutive bridge and finite polarized interface test.
Classification of all skewonless nonbirefringent media is a cited theorem,
not inferred from the samples below. Spacetime/Maxwell kinematics are inputs.
"""
from pathlib import Path
import argparse,hashlib,itertools,json
import numpy as np
HERE=Path(__file__).resolve().parent
I2=np.eye(2);J=np.array([[0.,-1.],[1.,0.]])
PAIRS=((0,1),(0,2),(0,3),(2,3),(3,1),(1,2))
EPS=np.zeros((4,4,4,4))
for p in itertools.permutations(range(4)):
    inv=sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
    EPS[p]=(-1)**inv

def product(a,b):return np.einsum('mr,ns->mnrs',a,b)-np.einsum('ms,nr->mnrs',a,b)

def constitutive(g,lam,alpha):
    gi=np.linalg.inv(g)
    return lam*np.sqrt(-np.linalg.det(g))*product(gi,gi)+alpha*EPS

def six(chi):return np.array([[chi[a,b,c,d] for c,d in PAIRS] for a,b in PAIRS])

def symmetric_entries(m):return m[np.triu_indices(6)]

def rank_response(g,lam):
    gi=np.linalg.inv(g);dens=np.sqrt(-np.linalg.det(g));base=product(gi,gi);cols=[]
    for a in range(4):
        for b in range(a,4):
            h=np.zeros((4,4));h[a,b]=h[b,a]=1
            dgi=-gi@h@gi
            dchi=lam*dens*(.5*np.trace(gi@h)*base+product(dgi,gi)+product(gi,dgi))
            cols.append(symmetric_entries(six(dchi)))
    cols += [symmetric_entries(six(dens*base)),symmetric_entries(six(EPS))]
    jac=np.asarray(cols).T;vals=np.linalg.svd(jac,compute_uv=False)
    assert np.linalg.matrix_rank(jac,tol=1e-10)==11
    return int(np.linalg.matrix_rank(jac,tol=1e-10)),float(vals[-1])

def interface(l1,l2,beta):
    A=(l1+l2)*I2+beta*J
    r=np.linalg.solve(A,(l1-l2)*I2-beta*J)
    t=np.linalg.solve(A,2*l1*I2)
    rr=np.linalg.solve(A,(l2-l1)*I2-beta*J)
    tn=np.linalg.solve(A,2*np.sqrt(l1*l2)*I2)
    S=np.block([[r,tn],[tn,rr]])
    denom=(l1+l2)**2+beta**2
    R=((l1-l2)**2+beta**2)/denom;T=4*l1*l2/denom
    # Independent four-equation boundary solve for incident left basis vectors.
    # c-b=a; H_left=lambda1 J(a-b), H_right=lambda2 J c-beta c.
    boundary=np.block([[-I2,I2],[-l1*J,-l2*J+beta*I2]])
    rhs=np.vstack([I2,-l1*J]);sol=np.linalg.solve(boundary,rhs)
    boundary_error=float(np.linalg.norm(sol[:2]-r)+np.linalg.norm(sol[2:]-t))
    flux_error=float(np.linalg.norm(r.T@r+(l2/l1)*t.T@t-I2))
    unitarity=float(np.linalg.norm(S.T@S-np.eye(4)))
    assert max(boundary_error,flux_error,unitarity)<1e-13
    assert np.linalg.norm(r.T@r-R*I2)<1e-13
    recovered_A=2*l1*np.linalg.inv(t)
    inferred_l2=.5*np.trace(recovered_A)-l1;inferred_beta=-.5*np.trace(J@recovered_A)
    assert max(abs(inferred_l2-l2),abs(inferred_beta-beta))<1e-13
    polarized=np.abs(S[:,0])**2
    force=1+R-T
    assert abs(force-2*R)<1e-13
    return dict(lambda_left=l1,lambda_right=l2,axion_jump=beta,
       reflection_matrix=r.tolist(),field_transmission_matrix=t.tolist(),
       flux_normalized_scattering_matrix=S.tolist(),reflection_power=R,transmission_power=T,
       outgoing_port_polarization_probabilities=polarized.tolist(),interface_momentum_per_incoming_energy=force,
       boundary_solve_error=boundary_error,flux_error=flux_error,full_scattering_unitarity_error=unitarity,
       reconstructed_lambda_right=float(inferred_l2),reconstructed_axion_jump=float(inferred_beta))

def run():
    eta=np.diag([-1.,1.,1.,1.]);coframe=np.array([[1.1,.1,0,0],[0,.9,.15,0],[0,0,1.2,.05],[0,0,0,1.]])
    g=coframe.T@eta@coframe;gi=np.linalg.inv(g);density=np.sqrt(-np.linalg.det(g))
    conformal_error=0.;symbol_error=0.;ranks=[]
    directions=[np.array(v,float)/np.linalg.norm(v) for v in ((1,0,0),(0,1,0),(0,0,1),(1,2,-3))]
    for lam,alpha in ((1.,0.),(2.,0.),(1.,2**-.5),(1.,1.)):
        chi=constitutive(g,lam,alpha)
        assert np.linalg.norm(six(chi)-six(chi).T)<1e-13
        for omega in (.7,1.3,2.):conformal_error=max(conformal_error,float(np.linalg.norm(constitutive(omega**2*g,lam,alpha)-chi)))
        for direction in directions:
            for freq,want_rank in ((1.,1),(.7,3)):
                q=coframe.T@np.r_[-freq,direction]
                actual=np.einsum('manb,a,b->mn',chi,q,q)
                v=gi@q;wanted=lam*density*((q@gi@q)*gi-np.outer(v,v))
                symbol_error=max(symbol_error,float(np.linalg.norm(actual-wanted)))
                r=int(np.linalg.matrix_rank(actual,tol=1e-10));assert r==want_rank;ranks.append(r)
    rank,smallest=rank_response(g,1.2)
    assert conformal_error<1e-12 and symbol_error<1e-12
    examples=[interface(1.,1.,0.),interface(1.,2.,0.),interface(1.,1.,2**-.5),interface(1.,1.,1.),interface(.7,1.6,.2)]
    a,b=examples[1:3]
    same_power=abs(a['reflection_power']-b['reflection_power'])
    same_force=abs(a['interface_momentum_per_incoming_energy']-b['interface_momentum_per_incoming_energy'])
    total_variation=.5*float(np.sum(abs(np.array(a['outgoing_port_polarization_probabilities'])-np.array(b['outgoing_port_polarization_probabilities']))))
    assert same_power<1e-14 and same_force<1e-14 and abs(total_variation-16/81)<1e-14
    # Positive energy and its parameter sources follow from the same quadratic Lagrangian.
    electric=np.array([.3,-.7,.2]);magnetic=np.array([-.4,.1,.9]);lam=1.6;alpha=.8
    lag=.5*lam*(electric@electric-magnetic@magnetic)+alpha*(electric@magnetic)
    d=lam*electric+alpha*magnetic;h=lam*magnetic-alpha*electric
    energy=float(electric@d-lag);want_energy=float(.5*lam*(electric@electric+magnetic@magnetic))
    energy_error=abs(energy-want_energy);poynting_error=float(np.linalg.norm(np.cross(electric,h)-lam*np.cross(electric,magnetic)))
    assert max(energy_error,poynting_error)<1e-14
    eps=.01;root=np.sqrt(eps);ratio_bounds=[(1-root)/(1+root),(1+root)/(1-root)]
    scaled_axion_bound=np.sqrt(eps/(1-eps))
    sat_l=interface(1.,ratio_bounds[1],0.);sat_a=interface(1.,1.,2*scaled_axion_bound)
    assert abs(sat_l['reflection_power']-eps)<1e-14 and abs(sat_a['reflection_power']-eps)<1e-14
    return dict(round=930,date='2026-10-06',all_scientific_checks_passed=True,
        common_cone_symbol_samples=len(ranks),null_symbol_rank=1,off_cone_symbol_rank=3,
        symbol_identity_max_error=symbol_error,conformal_response_max_error=conformal_error,
        metric_dilaton_axion_coefficient_jacobian_rank=rank,jacobian_smallest_singular_value=smallest,
        fixed_common_cone_remaining_coefficients=2,
        interfaces=examples,same_power_pair_difference=same_power,same_force_pair_difference=same_force,
        same_power_pair_full_polarization_total_variation=total_variation,
        common_lagrangian_energy_error=energy_error,common_lagrangian_poynting_error=poynting_error,
        source_lambda=float(.5*(electric@electric-magnetic@magnetic)),source_alpha=float(electric@magnetic),
        transparency_budget=eps,lambda_ratio_bounds=ratio_bounds,normalized_axion_jump_bound=float(scaled_axion_bound),
        transparent_exactly_iff_response_parameters_match_in_declared_class=True,
        full_response_readout_recovers_remaining_relative_parameters=True,
        nonbirefringent_classification_is_cited_not_numerically_proved=True,
        common_cone_matching_and_transparent_replacement_are_extra_hypotheses=True,
        no_independent_surface_layers_assumed=True,
        physical_volume_or_absolute_EM_normalization_generated=False,
        interface_carrier_dynamics_and_closed_backreaction_constructed=False,
        all_six_protocols_realized_in_EM_candidate=False,
        spacetime_dimension_Maxwell_kinematics_or_GR_generated=False,
        whole_stage_completed=False,full_goal_completed=False,
        source_hashes={str(Path('930')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    path=HERE/'cone_response_bridge_results.json'
    if a.write:
        assert not path.exists();path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','interfaces')},ensure_ascii=False,indent=2))
