"""978: both positions move, with the complete native internal energy as mass.
Analytic Duhamel/moment bound; no grid evolution is used as its proof.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'moving_complete_source_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def kron(*args):
    x=np.ones((1,1))
    for a in args:x=np.kron(x,a)
    return x

def run():
    old=read(STAGE/'976/finite_internal_reference_results.json')
    coarse=read(STAGE/'977/coarse_source_record_results.json')
    m,h,q,F,W,identities=load('native968',STAGE/'968/internal_relay.py').material()
    t0=old['parameters']['t0'];time_max=old['parameters']['coordinate_interval'][1]
    delta=math.pi/t0;G=1e-30;radius=1e10;core=radius/100;sigma=1e3
    mu=delta*radius/G;radii=np.sqrt(np.array([radius**2,(radius/2)**2])-core**2)
    eta=G*mu/core
    potentials=-G*mu/np.sqrt(radii*radii+core*core)
    assert max(abs(potentials-np.array(old['parameters']['Phi'])))<1e-19
    assert 0<eta<.003 and m['J']<.04
    # Uniform mathematical bounds, independently justified by bounded material
    # h,Q and the initial photon support. Do not use infinity-norm of full M.
    mass_lower=199.;rms_upper=204.
    assert 200-3*m['J']-.2>mass_lower
    material_norm_upper=3*1.04+.2
    photon_norm_upper=.5
    field_linear_upper=.003*2*math.sqrt(2)
    field_square_upper=.003**2/.5*4
    triangle_mass_upper=200+material_norm_upper+photon_norm_upper+field_linear_upper+field_square_upper
    assert triangle_mass_upper<rms_upper
    n=11;I=np.eye(4);ann=np.diag(np.sqrt(np.arange(1,n)),1)
    Hm=kron(h,I)+kron(I,h)+.2*kron(q,q);S=kron(q,I)+kron(I,q)
    HL=kron(Hm,np.eye(n))+.5*kron(np.eye(16),np.diag(np.arange(n)))
    HL+=.003*kron(S,ann+ann.T)+.000018*kron(S@S,np.eye(n))
    M=200*np.eye(352)+kron(HL,np.eye(2))+kron(np.eye(176),np.diag([-m['J'],0.]))
    photon=np.zeros(n);photon[:2]=1/math.sqrt(2)
    plus=W@np.ones(2)/math.sqrt(2);rplus=np.ones(2)/math.sqrt(2)
    states=np.column_stack([np.kron(np.kron(np.kron(W[:,s],plus),photon),rplus) for s in (0,1)])
    moment_matrix=states.T@M@M@states
    mass_means=np.diag(states.T@M@states)
    assert abs(moment_matrix[0,1])<1e-10
    assert math.sqrt(float(np.linalg.eigvalsh(moment_matrix).max()))<rms_upper
    for s in (0,1):
        expected=old['reference_account']['source_rows'][s]
        assert abs(mass_means[s]-expected['mass_mean_with_reference'])<1e-10
        assert abs(moment_matrix[s,s]-mass_means[s]**2-expected['mass_variance_with_reference'])<1e-9
    # For each 3D Gaussian position packet, ||p^2 psi||=sqrt(15)/(4 sigma^2).
    kinetic_rate=math.sqrt(15)/(8*sigma**2)*(1/mass_lower+1/mu)
    kinetic_error=time_max*kinetic_rate
    # Reuse 945's spatial lemma, not its old mass/source norm or error value.
    exponents=radii*radii/(32*sigma**2)
    assert min(exponents)>100*math.log(10)
    tail_probability_upper=2**1.5*1e-100
    tail_terms=(1/core+1/radii)**2*tail_probability_upper
    B2=24*sigma**2/radii**4+tail_terms
    B=np.sqrt(B2)
    potential_error=time_max*(G*mu)*rms_upper*float(max(B))
    numerical_guard=1e-10
    moving_bound=kinetic_error+potential_error+numerical_guard
    assert moving_bound<.005
    record_lower=old['interval_certificate']['uniform_infinite_contrast_lower']-2*moving_bound
    joint_coarse=coarse['uniform_budget']['maximum_joint_output_distance']+moving_bound
    assert record_lower>.48 and joint_coarse<.011
    # Independent Gaussian quadrature checks the spatial moment lemma. It is
    # not the certificate; analytic B above is authoritative.
    x,w=np.polynomial.hermite.hermgauss(16);w=w/math.sqrt(math.pi)
    xx=np.stack(np.meshgrid(x,x,x,indexing='ij'),axis=-1).reshape(-1,3)
    ww=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()
    fluct=2*sigma*xx
    rows=[]
    for l,r in enumerate(radii):
        relative=fluct+np.array([0.,0.,r]);den=np.sum(relative**2,axis=1)+core**2
        u=1/np.sqrt(den);u0=1/math.sqrt(r*r+core*core)
        spatial_variance=float(ww@((u-u0)**2))
        expected_force_kernel=np.sum(ww[:,None]*relative/den[:,None]**1.5,axis=0)
        force_means=[(G*mu*mass_means[s]*expected_force_kernel).tolist() for s in (0,1)]
        assert spatial_variance<float(B2[l])
        rows.append(dict(path=l,radius=float(r),potential=float(potentials[l]),
            spatial_difference_norm_bound=float(B[l]),quadrature_second_moment=spatial_variance,
            analytic_second_moment_upper=float(B2[l]),initial_force_X_by_sender=force_means,
            force_Y_is_exact_opposite=True))
    # Independent source/force finite differences in dimensionless positions,
    # avoiding subtraction of the extremely large constant rest mass mu.
    def potential_per_mass(x,y):return -delta/math.sqrt(float(np.dot(y-x,y-x))+.01**2)
    force_errors=[];translation_errors=[]
    step=1e-4
    for r in radii:
        xpos=np.array([.13,-.07,.02]);ypos=xpos+np.array([0.,0.,r/radius])
        d=ypos-xpos;expected=delta*d/(float(d@d)+.01**2)**1.5
        for i in range(3):
            e=np.eye(3)[i]*step
            deriv=(potential_per_mass(xpos-2*e,ypos)-8*potential_per_mass(xpos-e,ypos)
                   +8*potential_per_mass(xpos+e,ypos)-potential_per_mass(xpos+2*e,ypos))/(12*step)
            force_errors.append(abs(deriv+expected[i]))
        shift=np.array([.33,-.22,.07])
        translation_errors.append(abs(potential_per_mass(xpos+shift,ypos+shift)-potential_per_mass(xpos,ypos)))
    assert max(force_errors)<1e-12 and max(translation_errors)<1e-18
    # Duhamel needs moments on the known comparison trajectory, not on an
    # unproved frozen approximation to the true moving trajectory.
    maximal_force_norm=(G*mu)*rms_upper*2/(3*math.sqrt(3)*core**2)
    momentum_norm_upper=math.sqrt(3)/(2*sigma)+time_max*maximal_force_norm
    speed_norm_upper=momentum_norm_upper/mass_lower
    assert speed_norm_upper<1e-5
    files=[Path(__file__),HERE/'drafts/moving_source_decision.md',HERE/'drafts/common_recovery_audit_v2.md',
        STAGE/'945/joint_record_newton.py',STAGE/'956/native_material_interface.py',
        STAGE/'968/internal_relay.py',STAGE/'976/finite_internal_reference_results.json',
        STAGE/'977/coarse_source_record_results.json']
    return dict(round=978,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(G=G,mu=mu,reference_radius=radius,soft_core=core,position_sigma=sigma,
            time_max=time_max,radii=radii.tolist(),potentials=potentials.tolist()),
        operator_bounds=dict(weak_potential_global_upper=eta,mass_lower=mass_lower,
            initial_rms_mass_upper=rms_upper,analytic_triangle_mass_upper=triangle_mass_upper,
            full_H_lower_excluding_mu=(1-eta)*mass_lower,
            internal_occupation_not_physically_cutoff=True,
            mass_second_moment_exact_on_initial_support=moment_matrix.tolist(),
            mass_means=mass_means.tolist(),source_moment_matrix_off_diagonal=float(moment_matrix[0,1])),
        uniform_moving_error=dict(kinetic_rate=kinetic_rate,kinetic=kinetic_error,
            potential=potential_error,numerical_guard=numerical_guard,total=moving_bound,
            tail_probability_upper=tail_probability_upper,minimum_tail_exponent=float(min(exponents)),
            interval_starts_at_zero=True,arbitrary_sender_and_passive_reference=True,
            actual_moving_record_contrast_lower=record_lower,
            moving_output_to_977_cq_bound=joint_coarse),
        spatial_checks=rows,
        conservation=dict(total_momentum_conserved_exactly=True,total_energy_conserved_exactly=True,
            full_internal_mass_conserved_exactly=True,both_centres_are_dynamical=True,
            force_identity_difference=max(force_errors),translation_difference=max(translation_errors),
            force_norm_upper=maximal_force_norm,speed_norm_upper=speed_norm_upper,
            source_is_complete_internal_energy=True),
        scope=dict(newton_effective_model_input=True,actual_wavepacket_evolution_computed=False,
            uniform_bound_not_ODE_sampling=True,force_signal_above_measurement_error_claimed=False,
            full_post_Newtonian_budget=False,dynamical_Einstein_metric_derived=False,
            entire_SM_matching_completed=False,actual_recombination_device_completed=False,
            full_goal_completed=False),
        references=['https://arxiv.org/abs/1502.00971','https://arxiv.org/abs/gr-qc/9405057',
                    'https://arxiv.org/abs/hep-th/0409156'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-6,abs_tol=1e-11),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write:assert not OUT.exists()
    out=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as dest:json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(OUT))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
