"""693: explicit original Haar small-gap modulus and all-source average error.

Uses original two-time/two-direction finite box. No hard admissibility cutoff,
Monte Carlo sign inference, volume-uniform bound or RP conclusion.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_auxiliary_orbit_quadrature as orbit
import joint_rational_physical_limit as soft

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauge_sublevel_control_results.json'
base=soft.base
N=256;SEAMS=4;COORDS=14*SEAMS;DEGREE=22*SEAMS*N
GSTAR=2-np.sqrt(3)


def adjugate(a):
    n=len(a);out=np.zeros_like(a)
    for i in range(n):
        for j in range(n):
            out[j,i]=(-1)**(i+j)*np.linalg.det(np.delete(np.delete(a,i,axis=0),j,axis=1))
    return out


def cayley_su(a):
    n=len(a);minus=np.eye(n)-1j*a
    d=np.linalg.det(minus);delta=float(abs(d)**2)
    raw=(np.eye(n)+1j*a)@adjugate(minus)
    numerator=raw*d.conjugate();numerator[:,0]=raw[:,0]*d
    unit=np.linalg.solve(minus.T,(np.eye(n)+1j*a).T).T
    special=unit.copy();special[:,0]/=np.linalg.det(unit)
    assert np.max(abs(numerator/delta-special))<3e-13
    assert np.max(abs(special.conj().T@special-np.eye(n)))<3e-13
    assert abs(np.linalg.det(special)-1)<3e-13
    return special,numerator,delta


def rational_group(seed,scale):
    rng=np.random.default_rng(seed)
    def herm(n):
        raw=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        return scale*(raw+raw.conj().T)/2
    a=herm(3);b=herm(2);t=scale*rng.normal()
    c,cn,dc=cayley_su(a);w,wn,dw=cayley_su(b)
    dz=1+t*t;zp=(1+1j*t)**2;zm=zp.conjugate();z=zp/dz
    den=dc*dw*dz**6
    numerator=np.zeros((16,16),complex)
    numerator[:6,:6]=np.kron(cn,wn)*zp*dz**5
    numerator[6:9,6:9]=cn.conj()*zm**4*dw*dz**2
    numerator[9:12,9:12]=cn.conj()*zp**2*dw*dz**4
    numerator[12:14,12:14]=wn*zm**3*dc*dz**3
    numerator[14,14]=zp**6*dc*dw;numerator[15,15]=den
    direct=base.mass.dictionary.left_rep(c,w,z)
    error=float(np.max(abs(numerator/den-direct)));assert error<3e-12
    j=base.mass.dictionary.dictionary()
    assert np.max(abs(j@numerator@j.conj().T/den-base.prior.rep(c,w,z)))<3e-12
    return (c,w,z),error


def identity_seam_and_rational_check():
    g5=np.kron(np.kron(np.eye(4),base.internal.spin.G5),np.eye(16))
    g0=np.kron(np.kron(np.eye(4),base.internal.spin.GAMMA[3]),np.eye(16))
    t0=np.zeros((256,256),complex)
    for i,(t,x) in enumerate(base.prior.SITES):
        j=base.prior.SITES.index((1-t,x))
        t0[64*i:64*(i+1),64*j:64*(j+1)]=(-1 if t else 1)*np.eye(64)
    assert np.max(abs(t0.conj().T@t0-np.eye(256)))==0
    rows=[];errors=[]
    for seed,scale in ((69361,.2),(69362,1.),(69363,3.)):
        links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
        for i in range(4):
            group,err=rational_group(seed+i,scale)
            links[0,i]=base.prior.rep(*group);errors.append(err)
        _,_,_,h,gap=base.kernel(links)
        x=g5@h;xs=x-g0@t0;us=np.eye(256)-xs
        unit_error=float(np.max(abs(us.conj().T@us-np.eye(256))))
        assert unit_error<3e-12 and gap>=GSTAR-3e-12
        ws=(xs+xs.conj().T)/2
        identity_error=float(np.max(abs(xs.conj().T@xs-2*ws)))
        assert identity_error<3e-12
        _,logdet=np.linalg.slogdet(h)
        assert logdet>=N*np.log(GSTAR)-1e-11
        rows.append(dict(seed=seed,Cayley_scale=scale,actual_identity_seam_gap=gap,
            spatial_unitary_error=unit_error,spatial_norm_identity_error=identity_error,
            determinant_log_absolute=float(logdet),uniform_log_lower=float(N*np.log(GSTAR))))
    constants=[2**(n*(n-1))*math.prod(math.factorial(j) for j in range(1,n))/
               np.pi**(n*(n+1)/2) for n in (1,2,3)]
    assert all(c<1 for c in constants)
    return dict(rows=rows,uniform_identity_seam_gap_formula='2-sqrt(3)',
        uniform_identity_seam_gap=float(GSTAR),denominator_degree_per_original_group=22,
        original_left_charge_modules=[1,-4,2,-3,6,0],
        cleared_original_representation_error=max(errors),
        Cayley_Haar_density_maxima_U1_U2_U3=constants,
        small_gap_measure_not_inferred_from_samples=True)


def log_modulus(log_delta):
    alpha=1/(DEGREE*(COORDS+2))
    a=np.exp(((N-1)*np.log(3)-N*np.log(GSTAR))/DEGREE)
    constant=12*SEAMS/np.pi+4*COORDS*2.**COORDS*np.sqrt(19)*a
    return min(0.,np.log(constant)+alpha*log_delta),constant,alpha


def source_modulus_check():
    links,e,phis=soft.fixture();u,v,_,h,_=base.kernel(links)
    mat=orbit.entry.matrices(e,phis)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    exact_eps=v@v.conj().T-u@u.conj().T
    eps,_,_=soft.regulate(h,g5,.025,4096)
    exact=soft.soft(exact_eps,mat,0);approx=soft.soft(eps,mat,0)
    _,ry,_=orbit.entry.reflection.reflection_matrices(mat,np.eye(4)[[2,3,0,1]])
    terms=orbit.entry.source_rows(ry)
    actual=sum(c*orbit.entry.integral_coefficient(exact,z) for c,z,_ in terms)
    changed=sum(c*orbit.entry.integral_coefficient(approx,z) for c,z,_ in terms)
    eta=float(np.linalg.norm(eps-exact_eps,2));delta=float(abs(changed-actual))
    gold=(1+np.sqrt(5))/2;az=np.sqrt(5)/4
    bounds=[]
    for c,z,_ in terms:
        nz=np.linalg.norm(z,2);m=256+z.shape[1]//2
        radius=2+gold*nz;lip=.5+az*nz
        bounds.append(abs(c)*m*lip*radius**(m-1))
    bound=sum(bounds)*eta
    assert delta<=bound and eta>0
    # Independent Pfaffian perturbation bound, including a singular endpoint.
    rng=np.random.default_rng(69375);control=[]
    for singular in (False,True):
        a=rng.normal(size=(12,12))+1j*rng.normal(size=(12,12));a=(a-a.T)/20
        if singular:a[-2:,:]=0;a[:,-2:]=0
        b=rng.normal(size=(12,12))+1j*rng.normal(size=(12,12));b=(b-b.T)/200
        left=abs(base.pf(a+b)-base.pf(a))
        right=6*np.linalg.norm(b,2)*max(np.linalg.norm(a,2),np.linalg.norm(a+b,2))**5
        assert left<=right
        control.append(dict(singular_endpoint=singular,actual=float(left),bound=float(right)))
    table=[]
    for exponent in (-2,-1000,-10_000_000,-30_000_000,-60_000_000):
        value,c,alpha=log_modulus(exponent*np.log(10))
        table.append(dict(delta_base10_exponent=exponent,log10_probability_bound=float(value/np.log(10))))
    return dict(original_B_source_absolute_difference=delta,original_epsilon_error=eta,
        actual_source_Pfaffian_Lipschitz_log10_bound=float(np.log10(bound)),
        independent_Pfaffian_checks=control,Haar_coordinate_dimension=COORDS,
        determinant_polynomial_degree_bound=DEGREE,explicit_Haar_exponent_denominator=DEGREE*(COORDS+2),
        explicit_Haar_prefactor=float(c),modulus_log_examples=table,
        delta_half_bound_log10_threshold=float(np.log(.5/c)/alpha/np.log(10)),
        original_heat_moment_prefactor_not_numerically_evaluated=True,
        explicit_bound_not_a_practical_integral_sign_certificate=True)


def run():
    deps=('research_note_614.md','research_note_669.md','research_note_673.md','research_note_676.md',
          'research_note_677.md','research_note_686.md','research_note_692.md',
          'round693_drafts/finite_character_entry_results.json','round693_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=693,tests_run=2,failures=0,errors=0,
        original_uniform_witness=identity_seam_and_rational_check(),all_source_error_control=source_modulus_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(uniform_actual_Haar_small_gap_power_bound=True,
            all_original_spatial_backgrounds_retained=True,hard_spectral_cutoff_added=False,
            complete_original_source_average_error_in_finite_heat_moments=True,
            practical_full_integral_sign_certificate=False,
            general_dynamic_RP_or_HF_or_continuum_or_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();out=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    else:assert out==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=693,all_checks_passed=True,
        degree=DEGREE,exponent_denominator=DEGREE*(COORDS+2),
        half_bound_log10_delta=out['all_source_error_control']['delta_half_bound_log10_threshold'])))
