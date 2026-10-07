"""949: one parameter family preserves free masses and the actual internal h.

Numerical checks cover action scaling and inherited leading coefficients.
No full-field finite-time error window, empirical fit or autonomous apparatus
completion is inferred from these algebraic checks.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'joint_weak_domain_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))

def run():
    scalar=read(STAGE/'946/portal_common_process_results.json')
    protocol=read(STAGE/'947/protocol_field_transport_results.json')
    bridge=read(STAGE/'948/neutral_matter_bridge_results.json')
    assert all(r['all_scientific_checks_passed'] for r in (scalar,protocol,bridge))
    p=scalar['parameters'];bp=bridge['parameters']
    lam0=p['lambda_H'];v0=p['v'];alpha0=p['portal_alpha'];mu=p['mu']
    K0=np.array(p['mass_matrix']);gs0=bp['g_source'];g20=bp['g2'];yf0=bp['y_fermion']
    g10=.35 # Additional diagnostic U(1) input, not a fitted physical value.
    gravity=p['inherited_gravity_and_moving_parameters']
    G0=gravity['G'];mass=gravity['base_mass'];mr=gravity['receiver_mass']
    t=protocol['same_physical_protocol']['duration']
    n=30;kappa=protocol['same_physical_protocol']['internal_clock_scale']
    c=np.array([.5*np.sqrt((j+1)*(n-j)) for j in range(n)])
    hc=kappa*(15*np.eye(31)+np.diag(c,1)+np.diag(c,-1))
    w,U=np.linalg.eigh(hc);ket=np.eye(31)[:,0];finish=np.eye(31)[:,-1]
    clock0=U@(np.exp(-1j*t*w)*(U.conj().T@ket))
    assert np.linalg.norm(clock0-finish)<2e-12
    I=np.eye(2);X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1.])
    generators=[X/2,Y/2,Z/2,I/2]
    expect_gauge_masses=np.array([0,(g20*v0/2)**2,(g20*v0/2)**2,
                                  (g20*g20+g10*g10)*v0*v0/4])
    rng=np.random.default_rng(949)
    # A generic complex generation matrix verifies the full matrix statement;
    # it is a diagnostic input, not the old E_rec Yukawa table.
    Y0=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    M0=Y0*v0/math.sqrt(2)
    A=np.eye(2)+K0;D=np.linalg.inv(A)
    r1=next(r for r in bridge['rows'] if r['spacelike_Q']==1.)
    radii=np.array(gravity['radii'])
    invr=1/np.sqrt(radii*radii+gravity['soft_core']**2)
    grav_phase0=G0*mass*mr*t*(invr[0]-invr[1])
    potential_error=0.;gauge_error=0.;yukawa_error=0.;scaled_response_error=0.
    rows=[]
    fields=rng.normal(scale=.2,size=(30,5)) # p, eta, three Goldstones
    for eps in (1.,.5,.2,.1,.05,.01):
        v=v0/eps;alpha=alpha0/eps;lam=lam0*eps**2
        g2=eps*g20;g1=eps*g10;yf=eps*yf0;gs=eps*gs0
        K=np.array([[mu**2+2*lam*alpha**2,2*lam*alpha*v],
                    [2*lam*alpha*v,2*lam*v*v]])
        assert norm(K-K0)<1e-13
        phi_v=np.array([0.,v/math.sqrt(2)],complex)
        vectors=[g*(T@phi_v) for g,T in zip((g2,g2,g2,g1),generators)]
        mass_matrix=np.array([[2*np.real(np.vdot(a,b)) for b in vectors] for a in vectors])
        gauge_error=max(gauge_error,float(np.max(abs(np.linalg.eigvalsh(mass_matrix)-expect_gauge_masses))))
        actual_yukawa=(eps*Y0)*v/math.sqrt(2)
        yukawa_error=max(yukawa_error,norm(actual_yukawa-M0))
        for pp,eta,pi1,pi2,pi3 in fields:
            radial2=eta**2+pi1*pi1+pi2*pi2+pi3*pi3
            # Full complex doublet, without restricting to the radial field.
            H=np.array([(pi2+1j*pi1)/math.sqrt(2),
                        (v+eta-1j*pi3)/math.sqrt(2)])
            raw=lam*(np.vdot(H,H).real-v*v/2+alpha*pp)**2+.5*mu**2*pp**2
            L=v0*eta+alpha0*pp
            expanded=lam0*L*L+.5*mu**2*pp**2+eps*lam0*L*radial2+eps**2*lam0*radial2**2/4
            potential_error=max(potential_error,abs(raw-expanded))
        mw=g2*v/2;mf=yf*v/math.sqrt(2)
        cw=2*mw*mw/v;cf=mf/v
        C=np.array([[gs,0,0],[0,cw,cf]])
        response=-C.T@np.linalg.inv(np.eye(2)+K)@C
        base=np.array(r1['response_hessian'])
        scaled_response_error=max(scaled_response_error,norm(response/eps**2-base))
        G=eps**2*G0
        grav_phase=G*mass*mr*t*(invr[0]-invr[1])
        assert abs(grav_phase/eps**2-grav_phase0)<1e-13
        # Actual h and T are held fixed. The rejected shortcut changes h.
        shortcut=U@(np.exp(-1j*t*eps*w)*(U.conj().T@ket))
        shortcut_p=float(abs(np.vdot(finish,shortcut))**2)
        analytic_p=math.sin(eps*math.pi/2)**60
        assert abs(shortcut_p-analytic_p)<1e-12
        rows.append(dict(epsilon=eps,v=v,alpha=alpha,lambda_H=lam,
            scalar_masses_squared=np.linalg.eigvalsh(K).tolist(),
            gauge_masses_squared=np.linalg.eigvalsh(mass_matrix).tolist(),
            fermion_mass=mf,W_vertex=cw,fermion_vertex=cf,
            record_W_coefficient=float(response[0,1]),
            record_fermion_coefficient=float(response[0,2]),
            Newton_mass_path_phase=float(grav_phase),
            unchanged_internal_clock_finish_probability=float(abs(np.vdot(finish,clock0))**2),
            uniformly_scaled_h_finish_probability=shortcut_p,
            uniformly_scaled_h_analytic_probability=analytic_p))
    assert potential_error<2e-11 and gauge_error<1e-13 and yukawa_error<1e-13
    assert scaled_response_error<1e-13
    assert next(r for r in rows if r['epsilon']==.5)['uniformly_scaled_h_finish_probability']<1e-8
    files=[Path(__file__),STAGE/'946/portal_common_process_results.json',
           STAGE/'947/protocol_field_transport_results.json',
           STAGE/'948/neutral_matter_bridge_results.json']
    return dict(round=949,date='2026-10-07',all_scientific_checks_passed=True,
        parameter_family=dict(hbar=1.,time=t,gauge_and_Yukawa_and_record_power=1,
             Higgs_v_and_portal_alpha_power=-1,Higgs_lambda_and_Newton_G_power=2,
             material_rest_mass_and_internal_h_power=0,
             spacetime_coordinates_rescaled=False),
        checks=dict(full_doublet_potential_polynomial_error=potential_error,
             gauge_mass_matrix_spectrum_error=gauge_error,
             complex_Yukawa_mass_matrix_error=yukawa_error,
             joint_matter_response_scaling_error=scaled_response_error,
             unchanged_internal_clock_finish_error=float(np.linalg.norm(clock0-finish)),
             inherited_Newton_phase_coefficient=float(grav_phase0)),
        rows=rows,
        scope=dict(same_scalar_free_spectrum_preserved=True,
             gauge_and_Yukawa_masses_preserved=True,
             actual_947_h_not_scaled_to_zero=True,
             nonzero_cross_sector_leading_coefficients_for_every_positive_epsilon=True,
             conditional_common_window_lemma_only=True,
             actual_parent_common_menu_coefficients_bounded=False,
             fixed_experimental_parameter_point_certified=False,
             full_interacting_protocol_or_joint_observations_certified=False,
             full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-10,(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','rows')},ensure_ascii=False,indent=2))
    print(json.dumps([r for r in out['rows'] if r['epsilon'] in (.5,.1)],ensure_ascii=False,indent=2))
