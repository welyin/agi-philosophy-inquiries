"""986: one covariant mode Legendre map, native dipoles and TT source.

The common process and infinite-occupation record bound are analytic.
Finite matrices check only initial support moments and source identities.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'metric_dipole_bridge_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    sp=importlib.util.spec_from_file_location(name,p);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod);return mod
def kron(*xs):
    y=np.ones((1,1))
    for x in xs:y=np.kron(y,x)
    return y
def frac(x):return dict(exact=str(x),decimal=float(x))
def stable_exp_difference(v,k):
    # ||exp(k s)-1||_vac^2; use expm1 to retain terms when v<<eps_float.
    return math.sqrt(math.expm1(2*k*k*v)-2*math.expm1(k*k*v/2))

def run():
    m,h,q,_,w,_=load('native968_986',STAGE/'968/internal_relay.py').material()
    previous=read(STAGE/'976/finite_internal_reference_results.json')
    n=6;I=np.eye(4);ann=np.diag(np.sqrt(np.arange(1,n)),1)
    x=(ann+ann.T)/math.sqrt(2);p=(ann-ann.T)/(1j*math.sqrt(2))
    # Project x^2,p^2 from a larger space, retaining the top diagonal correction.
    a_big=np.diag(np.sqrt(np.arange(1,n+2)),1)
    xb=(a_big+a_big.T)/math.sqrt(2);pb=(a_big-a_big.T)/(1j*math.sqrt(2))
    x2=(xb@xb)[:n,:n];p2=(pb@pb)[:n,:n]
    hm=kron(h,I)+kron(I,h)+.2*kron(q,q);s=kron(q,I)+kron(I,q)
    g=.003;omega=.5
    hmfull=kron(hm,np.eye(n),np.eye(2))+kron(np.eye(16*n),np.diag([-m['J'],0.]))
    sx=kron(s,x,np.eye(2));ss=kron(s@s,np.eye(n),np.eye(2))
    xx=kron(np.eye(16),x2,np.eye(2));pp=kron(np.eye(16),p2,np.eye(2))
    unit=np.eye(32*n)
    def K(strain):
        return hmfull+omega/2*(math.exp(2*strain)*xx+math.exp(-2*strain)*pp-unit)+math.sqrt(2)*g*math.exp(strain)*sx+g*g/omega*ss
    k0=K(0)
    k_native=hmfull+omega*kron(np.eye(16),np.diag(np.arange(n)),np.eye(2))+math.sqrt(2)*g*sx+g*g/omega*ss
    native_identity=float(np.linalg.norm(k0-k_native))
    c0=omega*(xx-pp)+math.sqrt(2)*g*sx
    c2=2*omega*(xx+pp)+math.sqrt(2)*g*sx
    z=1e-3
    derivative=(K(-2*z)-8*K(-z)+8*K(z)-K(2*z))/(12*z)
    first_error=float(np.linalg.norm(derivative-c0,2))
    second=(K(z)-2*K(0)+K(-z))/(z*z)
    second_error=float(np.linalg.norm(second-c2,2))
    assert native_identity<1e-13 and first_error<1e-10 and second_error<1e-4
    photon=np.zeros(n);photon[:2]=1/math.sqrt(2)
    plus=w@np.ones(2)/math.sqrt(2);ref=np.ones(2)/math.sqrt(2)
    W=np.column_stack([np.kron(np.kron(np.kron(w[:,i],plus),photon),ref) for i in (0,1)])
    mismatch=math.sqrt(2)*g*sx
    mismatch_moment=(mismatch@W).conj().T@(mismatch@W)
    assert np.linalg.eigvalsh(mismatch_moment).min()>1e-7
    k_initial=float(np.linalg.norm(k0@W,2)); assert k_initial<4
    number_bound=(4+3*1.04+.2+g*g/omega*4+4*g)/(omega-2*g)
    assert number_bound<15
    # Same G as 978/985, unit local TT profile, canonical plus polarization.
    G=1e-30;Omega=.02;lam=math.sqrt(4*math.pi*G);v=lam*lam/(2*Omega)
    c1=stable_exp_difference(v,1);ce2=stable_exp_difference(v,2)
    time=previous['parameters']['coordinate_interval'][1]
    error=time*(15.75*ce2+.102*c1)
    # e^a-1 <= a/(1-a), and c_k^2<=e^(2k^2 v)-1.
    # With v<4e-28: c_2<6e-14, c_1<3e-14, certified rationally.
    v_upper=F(4,10**28);c2_upper=F(6,10**14);c1_upper=F(3,10**14)
    assert 8*v_upper/(1-8*v_upper)<c2_upper*c2_upper
    assert 2*v_upper/(1-2*v_upper)<c1_upper*c1_upper
    assert v<float(v_upper) and time<140000
    bound=140000*(F(63,4)*c2_upper+F(51,500)*c1_upper)
    assert bound<F(134,10**9)
    lower=F(4859481,10**7)-2*F(134,10**9)
    # Geometric vacuum lies within the declared weak-strain domain initially.
    strain_tail_variance_bound=v/1e-6  # Chebyshev for |s|>1e-3, initial only
    # Independent classical Legendre check: pi=e^-2s dot y/(N omega)-eta e^-s S.
    eta=math.sqrt(2)*g/omega; legendre=[]
    for strain,N,S,y,piy in [(.02,.97,1.2,.3,.7),(-.03,1.01,-.8,-.4,.2)]:
        dy=N*omega*math.exp(2*strain)*(piy+eta*math.exp(-strain)*S)
        L=math.exp(-2*strain)*dy*dy/(2*N*omega)-N*omega/2*math.exp(-2*strain)*y*y-eta*math.exp(-strain)*S*dy
        expected=N*omega/2*((math.exp(strain)*piy+eta*S)**2+math.exp(-2*strain)*y*y)
        legendre.append(abs(piy*dy-L-expected))
    assert max(legendre)<1e-15
    files=[Path(__file__),HERE/'drafts/metric_dipole_decision.md',STAGE/'968/internal_relay.py',
        STAGE/'976/finite_internal_reference_results.json',STAGE/'981/drafts/common_parent_contract_v1.md',
        STAGE/'985/two_active_sources_results.json']
    return dict(round=986,date='2026-10-07',all_scientific_checks_passed=True,
        inherited=dict(G=G,g=g,omega=omega,J=m['J'],same_native_generator_at_zero_strain=True),
        geometry_input=dict(TT_frequency=Omega,local_profile_value=1,lambda_strain=lam,
            vacuum_strain_variance=v,initial_weak_strain_tail_bound=strain_tail_variance_bound,
            dynamic_geometry_is_single_physical_TT_mode=True),
        operator_checks=dict(native_identity_error=native_identity,Legendre_errors=legendre,
            source_first_derivative_error=first_error,source_second_difference_error=second_error,
            omitted_mixed_source_second_moment=mismatch_moment.real.tolist(),
            mixed_source_is_fixed_by_existing_g=True),
        infinite_domain_bound=dict(initial_K_norm=k_initial,uniform_photon_number_norm_upper=15,
            derived_photon_number_norm_bound=number_bound,quadrature_square_norm_upper=31.5,
            c1=c1,c2=ce2,unrounded_Duhamel_bound=error,
            variance_rational_upper=frac(v_upper),c1_rational_upper=frac(c1_upper),
            c2_rational_upper=frac(c2_upper),rational_Duhamel_upper=frac(bound),
            reported_state_upper=frac(F(134,10**9)),record_contrast_lower=frac(lower),
            all_times_in_original_interval=True,unknown_input_and_passive_reference=True),
        source_dictionary=dict(lapse_source='m0+K(s)',TT_source='dK/ds',
            EM_source='dK/dg',total_joint_energy_conserved=True,
            geometric_force='-Omega^2 Q_G - lambda*N*dK/ds'),
        scope=dict(joint_quantum_mode_process_exists=True,record_transported=True,
            independent_quadrupole_matching_derived=False,full_U1_completed=False,
            spatial_mode_matching_error_certified=False,full_GR_error_certified=False,
            Newton_motion_985_and_TT_glued=False,all_SM_and_cosmology_completed=False,
            physical_minimum_scale_assumed=False,full_goal_completed=False),
        references=['https://arxiv.org/html/1908.06929v2','https://arxiv.org/html/1305.5231v3'],
        source_hashes={str(f.relative_to(ROOT)):sha(f) for f in files})
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-35),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write:assert not OUT.exists()
    out=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(OUT))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
