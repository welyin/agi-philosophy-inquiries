"""1260: same original positive material and a completed two-body 1PN process.
Three author groups; default read-only replay; --write-once exclusive first save.
Fraction fields are authoritative. Floating values are diagnostics/displays.
No normal preparation is momentum-projected. The fixed cutoff belongs to H.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
SRC=BASE/'archive_1106_1119/1112/check.py'
PROOF=HERE/'proof.md'
RESULT=HERE/'results.json'
INPUT_HASHES={
 SRC:'6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f',
 BASE/'archive_1234_/1253/proof.md':'422c736a1c3d65671fcef1795cc889e44514516344510048be7c8240acc29161'
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def calculate():
    for path,expected in INPUT_HASHES.items():
        assert sha(path)==expected,(str(path),'inherited hash mismatch')
    sp=importlib.util.spec_from_file_location('original1112',SRC)
    old=importlib.util.module_from_spec(sp); sp.loader.exec_module(old)
    mat=old.Material(1)
    shape=(4,8,3,3,1,3); dim=864
    I=np.eye(dim,dtype=complex); basis=I.reshape(shape+(dim,))
    S=mat.swap_external(basis).reshape(dim,dim)
    MD=mat.h_action(basis).reshape(dim,dim)
    M=MD+2*I-S
    Pi=np.diag(np.broadcast_to((np.arange(8)%2)[None,:,None,None,None,None],shape).reshape(dim))
    J3=mat.loaded_joint(np.eye(3),math.pi/4).reshape(dim,3)
    J=J3[:,1:3]
    e,W=np.linalg.eigh(M)
    k=np.rint(e-(mat.baseline+1)).astype(int)
    assert np.max(abs(e-(mat.baseline+1+k)))<1e-10
    C=W.conj().T@J
    spectral_residuals=[]
    for j,w in enumerate([1,3,4,4,3,1]):
        r=float(np.linalg.norm(C[k==j].conj().T@C[k==j]-w*np.eye(2)/16))
        assert r<1e-10
        spectral_residuals.append(r)
    phi=F(1,10**6); R=F(10**30); mu=phi*R; a=R/100; sigma=F(10**9)
    nu=F(1,100); N=80_000_000_000
    T=(2*math.pi*N+math.pi/2)/(1-float(phi)); Tmax=F(6*10**11)
    Delta=T*float(phi*phi)/2; beta=Delta/float(mu)
    # Remove only a state-independent common phase; round exact integral gaps.
    # Never exponentiate 1e12 times floating eigenvalue errors.
    rows=[]
    for label,d,b in [('delete_G2',0.,0.),('complete_1PN_static_comparator',Delta,beta)]:
        phases=np.exp(-1j*((math.pi/2+d)*k+b*(mat.baseline+1+k)**2))
        Q=W@(phases[:,None]*C)
        E=Q.conj().T@Pi@Q
        p=float(np.trace(E).real/2)
        assert np.linalg.norm(E-p*np.eye(2))<1e-10
        rows.append({'model':label,'Mr_probability':p,'all_input_effect_scalar_residual':float(np.linalg.norm(E-p*np.eye(2)))})
    # Same ideal internal age, different coordinate endpoints: no main redshift gap.
    Sage=2*math.pi*N+math.pi/2
    Tplus=Sage/(1-float(phi)+float(phi*phi)/2)
    beta_plus=Tplus*float(phi*phi)/(2*float(mu))
    Qsame=W@(np.exp(-1j*((math.pi/2)*k+beta_plus*(mat.baseline+1+k)**2))[:,None]*C)
    p_same=float(np.trace(Qsame.conj().T@Pi@Qsame).real/2)
    assert abs(p_same-rows[0]['Mr_probability'])<1e-12
    # Rational uniform norm envelopes: m in (26,32), normal Gaussian sigma.
    # sqrt15<4, sqrt3<2, sqrt6<3, sqrt135135<384.
    mmin=F(26); mmax=F(32)
    p2=1/sigma**2; p1=1/sigma; relative_position=3*sigma
    umax=F(3,2)/a; grad=1/a**2
    bounds={}
    bounds['exact_SR_kinetic']=Tmax*p2*(1/(2*mmin)+1/(2*mu))
    bounds['Newton_spatial_variation']=Tmax*mmax*mu*grad*relative_position
    bounds['G2_spatial_variation']=Tmax*mmax*mu*(mmax+mu)*umax*grad*relative_position
    bounds['G2_finite_source_M_squared']=Tmax*phi**2*mmax**2/(2*mu)
    pp=umax*p2+3*grad*p1
    bounds['1PN_self_velocity_all_three_axes']=Tmax*F(3,2)*(mu/mmin+mmax/mu)*pp
    bounds['1PN_scalar_cross_all_three_axes']=Tmax*F(7,2)*3*(umax/(4*sigma**2)+grad/(2*sigma))
    bounds['1PN_radial_cross_all_nine_indices']=Tmax*F(1,2)*9*(1/(4*a*sigma**2)+5/(2*a*a*sigma))
    bounds['total_static_comparator_norm']=sum(bounds.values())
    assert bounds['total_static_comparator_norm']<F(1,10**7)
    # Exact SR completion differs from its 1PN Taylor polynomial by <=p^6/(16m^5).
    sr_beyond_1pn=Tmax*F(6,16)/sigma**6*(1/mmin**5+1/mu**5)
    assert sr_beyond_1pn<F(1,10**49)
    # Global positive completion: velocity sandwich norm <=7 G M mu umax nu^2.
    positive_coefficient=1-mu*umax*(1+7*nu**2)
    assert (2*F(N)+F(1,2))*F(22,7)/(1-phi)<Tmax
    assert F(3,4)*F(9,10)**5>F(3,8)
    assert positive_coefficient>F(999,1000)
    # Phase interval theorem: baseline age=3pi/4 mod 2pi, Delta in (1/4,pi/12).
    Delta_lower=(F(2*N)+F(1,2))*F(25,8)*phi**2/(2*(1-phi))
    assert Delta_lower>F(1,4)
    assert (F(N)+F(1,4))*phi**2<(1-phi)/12
    prob_gap_lower=F(3,32)-2*F(1,10**7)
    assert prob_gap_lower>F(9,100)
    # Finite read/trigger comparison error .001 per model is an extra contract,
    # not manufactured here; its inclusion preserves a .09 lower bound.
    with_instrument=prob_gap_lower-2*F(1,1000)
    assert with_instrument>F(9,100)
    # Only the known 2PN static term, NOT a bound on the full GR remainder.
    known_2pn_static=Tmax*mmax*mu*(mmax*mmax+5*mmax*mu+mu*mu)/(4*R**3)

    def encode(x):return {'rational':str(x),'display':float(x)}
    result={'status':'PASS_1260','author_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'original_source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
     'parameters':{'hbar':1,'c':1,'G':1,'R':str(R),'mu':str(mu),'a':str(a),'sigma_per_axis':str(sigma),'cutoff_velocity':str(nu),'N':N,'T':T,'Delta_age':Delta},
     'physical_Mr_and_six_mass_coherences':rows,
     'same_internal_age_reference_deletion':{'T_minus':T,'T_plus':Tplus,'reference_probability_plus':p_same,'reference_probability_minus':rows[0]['Mr_probability'],'actual_effect_gap_bound':'1/5000000','physical_stopping_instrument_constructed':False},
     'bounds':{key:encode(value) for key,value in bounds.items()},
     'positive_H_coefficient':encode(positive_coefficient),
     'SR_completion_beyond_1PN_on_same_comparator':encode(sr_beyond_1pn),
     'certified_Mr_gap_with_mechanical_error':str(prob_gap_lower),
     'conditional_gap_with_extra_instrument_error_001_each':str(with_instrument),
     'known_2PN_static_piece_only':encode(known_2pn_static),
     'scope':{'actual_two_body_motion_and_reaction':True,'same_1253_mass':True,'unknown_reference_preserved':True,'all_EIH_velocity_cross_terms_retained':True,'source_and_geometric_ruler_dictionary_proved':False,'full_GR_remainder_certified':False,'actual_preparation_or_reader_constructed':False,'source_geometry_not_inferred_from_clock_readout':True}}

    result['proof_sha256']=sha(PROOF)
    result['inherited_source_hashes']={str(p.relative_to(BASE)).replace('\\','/'):sha(p) for p in INPUT_HASHES}
    result['scientific_groups']=3
    result['groups']={
     '1_original_material_actual_Mr':{
      'dimension':dim,'mass_range_display':[float(e[0]),float(e[-1])],
      'six_mass_weights_over_16':[1,3,4,4,3,1],
      'spectral_residuals':spectral_residuals,
      'full_unknown_input_effects':result.pop('physical_Mr_and_six_mass_coherences')},
     '2_rational_positivity_and_full_state_budget':{
      'bounds':result.pop('bounds'),
      'positive_H_coefficient':result.pop('positive_H_coefficient'),
      'phase_derivative_lower':encode(F(3,4)*F(9,10)**5),
      'certified_Mr_gap_with_mechanical_error':result.pop('certified_Mr_gap_with_mechanical_error')},
     '3_endpoint_deletion_and_SR_scope':{
      'same_internal_age_reference_deletion':result.pop('same_internal_age_reference_deletion'),
      'SR_completion_beyond_1PN_on_same_comparator':result.pop('SR_completion_beyond_1PN_on_same_comparator'),
      'known_2PN_static_piece_only':result.pop('known_2PN_static_piece_only'),
      'conditional_gap_with_extra_instrument_error_001_each':result.pop('conditional_gap_with_extra_instrument_error_001_each')}
    }
    return result

def compare(actual,saved,path='root'):
    if isinstance(saved,dict):
        assert actual.keys()==saved.keys(),(path,'keys')
        for key in saved: compare(actual[key],saved[key],path+'.'+key)
    elif isinstance(saved,list):
        assert len(actual)==len(saved),(path,'length')
        for i,(a,b) in enumerate(zip(actual,saved)): compare(a,b,f'{path}[{i}]')
    elif isinstance(saved,float):
        assert math.isclose(actual,saved,rel_tol=2e-8,abs_tol=3e-9),(path,actual,saved)
    else:
        assert actual==saved,(path,actual,saved)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args()
    result=calculate()
    if args.write_once:
        with RESULT.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('WRITE_ONCE_PASS_1260: 3 author groups; exclusive results.json created')
    else:
        compare(result,json.loads(RESULT.read_text(encoding='utf-8-sig')))
        print('READ_ONLY_MATCH_1260: 3 author groups; no files changed')
