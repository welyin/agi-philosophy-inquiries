"""655: full-graph positive splitting, common time sources and normalization.

Full boson/Gauss convergence is analytic in the note. Numerics retain the
original 64 CAR modes of a two-node conditional spatial factor; they do not
diagonalize the entire confining bosonic quantum model or a chiral field theory.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_mass_state_time_limit as mass
import joint_fermion_gauss_completion as old
import joint_spinor_subgroup_mass as dictionary
import joint_gauss_fermion_influence as car

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_full_graph_transfer_sources_results.json'


def norm(x):return float(np.linalg.norm(x,2))


def blockdiag(a,b):
    out=np.zeros((len(a)+len(b),len(a)+len(b)),complex)
    out[:len(a),:len(a)]=a;out[len(a):,len(a):]=b
    return out


def nambu(h,d):return np.block([[h,d],[d.conj().T,-h.T]])


def graph_data(neutral=False):
    p=mass.PHIS[:2].copy()
    if neutral:p=np.array([[0.,0.,0.,0.,.62],[0.,0.,0.,0.,-.41]])
    h1,d1=old.mass_matrices(p[0]);h2,d2=old.mass_matrices(p[1])
    hm=blockdiag(h1,h2);dm=blockdiag(d1,d2)
    c=old.gauge.group_exp(np.array([.17,-.11,.08,.05,.09,-.06,.04,.12]),3)
    w=old.gauge.group_exp(np.array([.21,-.16,.28]),2);z=np.exp(.13j)
    r=old.representation(c,w,z)
    hop=np.zeros((64,64),complex);hop[:32,32:]=.23*r;hop[32:,:32]=.23*r.conj().T
    ph=dictionary.ph_matrix();u=blockdiag(ph[:32,:32],ph[:32,:32]);v=blockdiag(ph[:32,32:],ph[:32,32:])
    transform=np.block([[u,v],[v.conj(),u.conj()]])
    bm=transform@nambu(hm,dm)@transform.conj().T
    bn=transform@nambu(hop,np.zeros_like(hop))@transform.conj().T
    assert norm(bm[:64,:64])<1e-13 and norm(bn[:64,64:])<1e-13
    return dict(phi=p,link=(c,w,z),h=bn[:64,:64],d=bm[:64,64:],
                normal_BdG=bn,target_BdG=bm+bn,transform=transform)


def multiply_jets(a,b):
    return [a[0]@b[0],a[1]@b[0]+a[0]@b[1],
            a[2]@b[0]+2*a[1]@b[1]+a[0]@b[2]]


def step_jets(bn,d,a,eta=1.):
    e=car.exp_h(bn,a*eta/2)
    jets=[e,-a/2*bn@e,(a*a/4)*(bn@bn)@e]
    eye=np.eye(len(d));zero=np.zeros_like(d);p=d@d.conj().T
    middle=[mass.bdg_transfer(d,a*eta),
        np.block([[2*a*a*eta*p,-a*d],[-a*d.conj().T,zero]]),
        np.block([[2*a*a*p,zero],[zero,zero]])]
    return multiply_jets(multiply_jets(jets,middle),jets)


def power_jets(s,n):
    size=len(s[0]);zero=np.zeros((size,size),complex)
    out=[np.eye(size,dtype=complex),zero,zero]
    for _ in range(n):out=multiply_jets(s,out)
    return out


def sources(jets):
    s,d,dd=jets;solve=np.linalg.solve(np.eye(len(s))+s,d)
    first=.5*np.trace(solve)
    second=.5*np.trace(np.linalg.solve(np.eye(len(s))+s,dd)-solve@solve)
    assert abs(first.imag)<1e-10 and abs(second.imag)<1e-10
    return float(first.real),float(second.real)


def logz_from_step(s,n):
    vals=np.linalg.eigvalsh(s);assert min(vals)>0
    x=n*np.log(vals)/2
    return float(.5*np.sum(np.logaddexp(x,-x)))


def spatial_and_sources_check():
    data=graph_data();bn=data['normal_BdG'];d=data['d'];b=data['target_BdG'];beta=1.
    exact=car.exp_h(b,beta);target_g=np.linalg.inv(np.eye(128)+exact)
    target_sources=sources([exact,-beta*b@exact,beta*beta*b@b@exact]);rows=[]
    for n in (2,4,8,16):
        a=beta/n;jets=step_jets(bn,d,a);power=power_jets(jets,n)
        src=sources(power);g=np.linalg.inv(np.eye(128)+power[0])
        minimum=float(np.linalg.eigvalsh(jets[0])[0]);assert minimum>0
        errors=[abs(src[i]-target_sources[i]) for i in range(2)]
        rows.append(dict(time_sites=n,time_step=a,positive_step_minimum=minimum,
            full64_mode_covariance_error=norm(g-target_g),
            lapse_first_source=src[0],lapse_second_source=src[1],source_errors=errors))
    assert all(rows[i+1]['full64_mode_covariance_error']<.55*rows[i]['full64_mode_covariance_error'] for i in range(3))
    assert all(rows[-1]['source_errors'][j]<rows[0]['source_errors'][j]/30 for j in range(2))
    # Independent finite differences of the full64 partition.
    a=.25;n=4;eps=2e-4
    f0=logz_from_step(step_jets(bn,d,a)[0],n)
    fp=logz_from_step(step_jets(bn,d,a,1+eps)[0],n)
    fm=logz_from_step(step_jets(bn,d,a,1-eps)[0],n)
    fd=[(fp-fm)/(2*eps),(fp-2*f0+fm)/eps**2]
    analytic=sources(power_jets(step_jets(bn,d,a),n))
    fd_errors=[abs(fd[i]-analytic[i]) for i in range(2)]
    assert max(fd_errors)<2e-6
    # Independently gauge-transform both endpoints and the original link.
    rng=np.random.default_rng(655);rr=[];pp=[]
    for p in data['phi']:
        c=old.gauge.group_exp(.2*rng.normal(size=8),3)
        w=old.gauge.group_exp(.2*rng.normal(size=3),2);z=np.exp(.15j*rng.normal())
        r=old.representation(c,w,z);rr.append(r)
        x=z**3*w@(p[:2]+1j*p[2:4]);pp.append(np.r_[x.real,x.imag,p[4]])
    rh=blockdiag(*rr);rn=blockdiag(rh,rh.conj());rt=data['transform']@rn@data['transform'].conj().T
    original_link=old.representation(*data['link']);newlink=rr[0]@original_link@rr[1].conj().T
    hh=np.zeros((64,64),complex);hh[:32,32:]=.23*newlink;hh[32:,:32]=.23*newlink.conj().T
    nn=data['transform']@nambu(hh,np.zeros_like(hh))@data['transform'].conj().T
    p1,q1=old.mass_matrices(pp[0]);p2,q2=old.mass_matrices(pp[1])
    mm=data['transform']@nambu(blockdiag(p1,p2),blockdiag(q1,q2))@data['transform'].conj().T
    cov=norm(step_jets(nn,mm[:64,64:],a)[0]-rt@step_jets(bn,d,a)[0]@rt.conj().T)
    assert cov<2e-12
    return dict(rows=rows,original_full64_mode_target_sources=list(target_sources),
        independent_source_difference_errors=fd_errors,local_gauge_covariance_error=cov,
        actual_hopping_mass_commutator_norm=norm(bn@(b-bn)-(b-bn)@bn),
        hopping_strength_is_allowed_diagnostic_input=.23,
        bosonic_path_integral_or_full_Gauss_partition_not_numerically_computed=True)


def neutral_fock_check():
    data=graph_data(neutral=True);ids=[30,31,62,63]
    h=data['h'][np.ix_(ids,ids)];d=data['d'][np.ix_(ids,ids)]
    nf=car.fock(h,np.zeros_like(d));c=mass.creation(d);hf=nf+c+c.conj().T
    beta=1.;a=.25;n=4;e=car.exp_h(nf,a/2);l=mass.nil_exp(c,a)
    tf=e@l@l.conj().T@e;heat=np.linalg.matrix_power(tf,n)
    bu=step_jets(nambu(h,np.zeros_like(d)),d,a)[0]
    trace=float(np.trace(heat).real)
    trace_error=abs(np.log(trace)-logz_from_step(bu,n))
    assert trace_error<1e-12
    rho=car.exp_h(hf,1.3);rho/=np.trace(rho)
    # A bounded local number-effect is a diagnostic, not an implemented recorder.
    effect=np.diag([.4+.1*((s&1)>0) for s in range(16)])
    target=np.trace(effect@car.exp_h(hf,beta)@rho@car.exp_h(hf,beta)).real
    rows=[]
    for n in (2,4,8,16):
        a=beta/n;e=car.exp_h(nf,a/2);l=mass.nil_exp(c,a)
        tf=e@l@l.conj().T@e;f=np.linalg.matrix_power(tf,n)
        actual=np.trace(effect@f@rho@f).real
        rows.append(dict(time_sites=n,original_state_weighted_effect_error=float(abs(actual-target))))
    # General inserted state statistics have O(a), not the spectral O(a^2), bias.
    assert all(rows[i+1]['original_state_weighted_effect_error']<.56*rows[i]['original_state_weighted_effect_error']
               for i in range(3))
    return dict(actual_decoupled_original_neutral_modes=4,independent_Fock_log_trace_error=trace_error,
                original_state_weighted_rows=rows,not_a_full_bosonic_record_simulation=True)


def stability_and_auxiliary_scope_check():
    mat,u,_=old.original.lattice.scalar.parameters();lam=np.linalg.det(mat)/np.trace(mat)
    big_d=6*old.original.M-u.sum();aa=32/(lam*big_d**2);bb=144/big_d**2
    linear=[]
    for j in range(5):
        e=np.eye(5)[j]
        numerator=mass.delta(e)*np.sqrt(old.original.F(e))
        linear.append(float(np.sum(abs(np.triu(numerator,1)))))
    kappa=np.sqrt(6*old.original.M)*np.linalg.norm(linear)
    direction=np.array([1.,2.,-.4,.6,.8]);direction/=np.linalg.norm(direction)
    rows=[]
    for z in (1.+0j,.3+.8j,1.2-.4j):
        c=2*abs(z)*kappa*aa**.25;b0=2*abs(z)*kappa*bb**.25
        bound=b0+3*c**(4/3)/(4**(4/3)*z.real**(1/3))
        maximum=-float('inf')
        for f in (.8,.2,.01,1e-4,1e-6):
            p=direction*np.sqrt(6*(old.original.M-f));v=float(old.original.node_potential(p))
            cupper=float(np.sum(abs(np.triu(mass.delta(p),1))))
            assert cupper<=kappa/np.sqrt(old.original.F(p))*(1+1e-12)
            log_upper=-z.real*v+2*abs(z)*cupper
            assert log_upper<=bound
            maximum=max(maximum,log_upper)
        rows.append(dict(complex_lapse=[z.real,z.imag],analytic_global_bound=float(bound),
                         maximum_sampled_log_norm_bound=float(maximum)))
    gap=np.log(17/8);tau=.7;n=25
    a_even=gap*tau/(2*np.pi*n);a_odd=gap*tau/(2*np.pi*n+np.pi)
    v0=np.array([1.,1.])/np.sqrt(2)
    states=[]
    for a in (a_even,a_odd):
        v=v0*np.array([1,np.exp(-1j*tau*gap/a)]);states.append(np.outer(v,v.conj()))
    distance=float(np.sum(abs(np.linalg.eigvalsh(states[0]-states[1])))/2)
    assert abs(distance-1)<1e-13
    return dict(complex_time_stability_rows=rows,full_configuration_bound_is_analytic=True,
                original_pair_creation_bound_coefficient=float(kappa),
                excited_auxiliary_realtime_subsequence_trace_distance=distance,
                Euclidean_projection_does_not_imply_realtime_forgetting=True)


def normalization_geometry_check():
    data=json.loads((HERE/'joint_chiral_character_state_results.json').read_text('utf8'))
    lam0=float(Fraction(data['temporal_gluing']['spherical_harmonic_spectrum'][0]['eigenvalue']))
    k=32*np.log(2)-np.log(lam0);a=.25;eps=.8;psi0=1.1
    w0=eps**3*psi0**6;energy=k/a;density=energy/w0;rows=[]
    for psi in (.9,1.1,1.3):
        volume=eps**3*psi**6;cosmo=density*volume
        step=1e-5
        fd=(density*eps**3*(psi+step)**6-density*eps**3*(psi-step)**6)/(2*step)
        expected=6*density*eps**3*psi**5
        assert abs(fd-expected)<2e-6
        rows.append(dict(psi=psi,original_node_volume=volume,count_energy=energy,
            volume_energy=cosmo,count_spatial_derivative=0.,volume_spatial_derivative=expected,
            derivative_difference_check=float(abs(fd-expected))))
    assert abs(rows[1]['volume_energy']-energy)<1e-12 and rows[1]['volume_spatial_derivative']>1
    # Additive volumes supply a consistent scalar factor; bare per-cell counts do not.
    volumes=np.array([.17,.23,.41]);rate=.37
    scalar_error=abs(np.prod(np.exp(-rate*volumes))-np.exp(-rate*sum(volumes)))
    assert scalar_error<1e-15
    return dict(actual_vacuum_coefficient=float(k),matched_reference_volume=w0,
        rows=rows,volume_additive_scalar_gluing_error=float(scalar_error),
        bare_log_normalization_change_on_bisecting_one_cell=float(-k),
        cosmological_coefficient_not_determined_by_flat_transfer=True,
        assumed_sitewise_extension_is_not_actual_spatial_overlap_measure=True)


def run():
    deps=('research_note_585.md','research_note_598.md','research_note_623.md','research_note_640.md',
          'research_note_643.md','research_note_645.md','research_note_646.md','research_note_654.md',
          'joint_mass_state_time_limit.py','joint_fermion_gauss_completion.py',
          'joint_spinor_subgroup_mass.py','joint_gauss_fermion_influence.py',
          'joint_chiral_character_state_results.json')
    return dict(date='2026-10-02',round=655,tests_run=4,failures=0,errors=0,
        spatial_and_sources=spatial_and_sources_check(),neutral_fock=neutral_fock_check(),
        stability_and_auxiliary_scope=stability_and_auxiliary_scope_check(),
        normalization_geometry=normalization_geometry_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Same original finite-graph confining quantum H admits the declared mass splitting, positive-time lapse derivatives and fixed-state records. Spatial conditional numerics retain64 modes. No actual multidimensional overlap-measure identity, full spatial continuum, arbitrary auxiliary real-time limit or quantum GR.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
