"""946: a gauge-invariant portal's quadratic sector in the SAME record/Newton H.

This is a finite-mass witness for the specified two-scalar effective Hamiltonian.
The full nonlinear SM/Einstein parent matching is explicitly not certified.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'portal_common_process_results.json'
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def norm(a):return float(np.linalg.norm(a,2))
def static_kernel(mu,r):
    if r==0:return 1/(4*np.pi**1.5)-mu*math.exp(mu*mu)*math.erfc(mu)/(4*np.pi)
    return math.exp(mu*mu)/(8*np.pi*r)*(math.exp(-mu*r)*math.erfc(mu-r/2)
               -math.exp(mu*r)*math.erfc(mu+r/2))
def oscillatory_bound(mu,t):
    # For C_mu=int k² exp(-k²)cos(T omega)/(2pi² omega²) dk:
    # one integration by parts and the single maximum of k exp(-k²)/omega.
    x=(-mu*mu+math.sqrt(mu**4+2*mu*mu))/2
    return math.sqrt(x)*math.exp(-x)/math.sqrt(x+mu*mu)/(np.pi**2*t)

def compute(n=256):
    old=load('field944',STAGE/'944/finite_field_communication.py')
    t=30.;radii=np.array([2.,3.]);gs=gr=.5;s=8.
    lam=.5;v=2.;alpha=.5;mu=1.
    K=np.array([[mu*mu+2*lam*alpha*alpha,2*lam*alpha*v],
                [2*lam*alpha*v,2*lam*v*v]])
    eig,O=np.linalg.eigh(K);m=np.sqrt(eig)
    assert eig.min()>0 and abs(np.linalg.det(K)-2*lam*v*v*mu*mu)<1e-12
    # Sectors use the SAME original source charge, receiver, and path labels.
    sectors=[(l,b,y) for l in (0,1) for b in (1.,-1.) for y in (-1.,1.)]
    k,w=old.nodes(n,0.,12.);q,wq=old.nodes(48,-1.,1.)
    measure=(w*k*k/(4*np.pi**2))[:,None]*wq[None,:]
    amplitudes=[];phases=np.zeros(8);field_read=[]
    energy=np.zeros((8,2));thetas=np.zeros(2);Drr=0.;Ddiff=0.
    radial=w*k*k*np.exp(-k*k)/(2*np.pi**2)
    read_variance=0.;mean_source=0.;mean_receiver=np.zeros(2)
    for a,ma in enumerate(m):
        om=np.sqrt(k*k+ma*ma)[:,None];f=np.exp(-k*k/2)[:,None]
        base=f/np.sqrt(2*om)*np.ones((1,len(q)))
        arr=[]
        for j,(l,b,y) in enumerate(sectors):
            source=O[0,a]*base*(gs*b+gr*y*np.exp(-1j*k[:,None]*radii[l]*q[None,:]))
            disp=source*(np.exp(-1j*om*t)-1)/om
            arr.append(disp.reshape(-1)*np.sqrt(measure.reshape(-1)))
            phases[j]+=float(np.sum(measure*abs(source)**2*(om*t-np.sin(om*t))/om**2))
            energy[j,0]+=float(np.sum(measure*om*abs(disp)**2))
            energy[j,1]+=float(2*np.real(np.sum(measure*source.conj()*disp)))
        amplitudes.append(np.asarray(arr))
        field_read.append((O[1,a]*base*np.sqrt(measure)).reshape(-1))
        omega=om[:,0]
        kt=(omega*t-np.sin(omega*t))/omega**3
        noise=(1-np.cos(omega*t))/omega**3
        rr=(1-np.cos(omega*t))/omega**2
        thetas+=gs*gr*O[0,a]**2*np.array([np.sum(radial*np.sinc(k*r/np.pi)*kt) for r in radii])
        Drr+=gr*gr*O[0,a]**2*np.sum(radial*noise)
        Ddiff+=gr*gr*O[0,a]**2*np.sum(radial*(1-np.sinc(k/np.pi))*noise)
        read_variance+=O[1,a]**2*np.sum(radial/(2*omega))
        mean_source-=gs*O[0,a]*O[1,a]*np.sum(radial*rr)
        mean_receiver-=gr*O[0,a]*O[1,a]*np.array([np.sum(radial*np.sinc(k*r/np.pi)*rr) for r in radii])
    amplitude=np.concatenate(amplitudes,axis=1);d=np.concatenate(field_read)
    overlap=amplitude@amplitude.conj().T;diag=np.diag(overlap).real
    field_gram=np.exp(1j*(phases[:,None]-phases[None,:]+overlap.imag)
                      -.5*(diag[:,None]+diag[None,:]-2*overlap.real))
    old945=json.loads((STAGE/'945/joint_record_newton_results.json').read_text('utf-8'))
    params=old945['parameters'];mass=params['base_mass'];mr=params['receiver_mass']
    G=params['G'];core=params['soft_core'];delta=params['position_std']
    kappas=G*mr*(mass+np.arange(31))
    weights=np.array([math.comb(30,j)/2**30 for j in range(31)])
    u=1/np.sqrt(radii*radii+core*core)
    grav=np.array([[np.sum(weights*np.exp(1j*kappas*t*(u[l]-u[ll])))
                   for ll,bb,yy in sectors] for l,b,y in sectors])
    gram=field_gram*grav
    z=amplitude@d
    def sine_kernel(g):
        return g*math.exp(-s*s*read_variance/2)*np.sin(s*(z[:,None]+z.conj()[None,:]))
    sg=sine_kernel(gram)
    assert abs(float(d@d)-read_variance)<1e-13
    gram_errors=max(norm(gram-gram.conj().T),float(np.max(abs(np.diag(gram)-1))))
    assert gram_errors<1e-12
    instrument_eigs=[float(np.linalg.eigvalsh((gram+sign*sg)/2).min()) for sign in (-1,1)]
    assert min(instrument_eigs)>-1e-12

    yy,vy=np.linalg.eigh(Y);rot=np.kron(vy,vy);plus=np.array([1,1],complex)/np.sqrt(2)
    prep=rot.conj().T@np.kron(plus,plus);rho_ra=np.outer(prep,prep.conj())
    charge=(yy[:,None]+yy[None,:]).reshape(-1)
    rho_ra*=np.isclose(charge[:,None],charge[None,:])
    rel=rot.conj().T@((np.kron(Z,X)-np.kron(X,Z))/2)@rot
    read=(np.eye(4)+rel)/2;ql=(I-Y)/2
    path_rho=np.outer(plus,plus.conj())
    rows=[];tables=[]
    for bi in (0,1):
        rho_s=np.zeros((2,2));rho_s[bi,bi]=1
        rho0=np.kron(path_rho,np.kron(rho_s,rho_ra)).reshape(8,2,8,2)
        def density(kernel):return (rho0*kernel[:,None,:,None]).reshape(16,16)
        rho=density(gram);rho_zero=density(field_gram)
        rho_sin=density(sg)
        row=dict(record=float(np.real(np.trace(rho@np.kron(I,np.kron(I,read))))),
            path=float(np.real(np.trace(rho@np.kron(ql,np.eye(8))))),
            path_G0=float(np.real(np.trace(rho_zero@np.kron(ql,np.eye(8))))),
            higgs=float(.5*(np.trace(rho)+np.trace(rho_sin)).real))
        table=[]
        for il in (0,1):
            for ir in (0,1):
                for ih in (0,1):
                    path_effect=ql if il else I-ql
                    record_effect=read if ir else np.eye(4)-read
                    hg=(gram+sg)/2 if ih else (gram-sg)/2
                    effect=np.kron(path_effect,np.kron(I,record_effect))
                    table.append(float(np.real(np.trace(density(hg)@effect))))
        assert min(table)>-1e-12 and abs(sum(table)-1)<1e-12
        assert abs(sum(table[4:])-row['path'])<1e-12
        assert abs(sum(table[j] for j in (2,3,6,7))-row['record'])<1e-12
        assert abs(sum(table[j] for j in (1,3,5,7))-row['higgs'])<1e-12
        rows.append(row);tables.append(table)
    record_gap=rows[0]['record']-rows[1]['record']
    higgs_gap=rows[0]['higgs']-rows[1]['higgs']
    gravity_contrast=rows[0]['path']-rows[0]['path_G0']
    analytic_record=float(.25*np.exp(-2*Drr)*np.sum(np.sin(2*thetas)))
    analytic_higgs=float(np.exp(-s*s*read_variance/2)*np.sin(s*mean_source)*np.mean(np.cos(s*mean_receiver)))
    visibility=float(np.exp(-Ddiff)*np.cos(thetas[0]-thetas[1]))
    analytic_gravity=float(.5*visibility*np.sum(weights*np.sin(kappas*t*(u[0]-u[1]))))
    formula_error=max(abs(record_gap-analytic_record),abs(higgs_gap-analytic_higgs),abs(gravity_contrast-analytic_gravity))
    assert formula_error<1e-12

    # Analytic finite-mass control: both normal modes are kept, no Fock cutoff.
    I0=1/(8*np.pi**1.5);I3=I0/m.min()**3;C=3/(32*np.pi**1.5*m.min())
    A2=4*(gs*gs+gr*gr)*I3
    kinetic=t*2*np.sqrt(15)/(8*mass*delta*delta)
    position=t*(abs(gs)+abs(gr))*delta*np.sqrt(C)*np.sqrt(1+4*A2)
    error0=kinetic+position
    gravity_error=old945['finite_moving_control']['added_gravity_position_bound']
    error=error0+gravity_error
    # These bounds depend only on elementary integrals / erfc, not oscillatory quadrature.
    mix=abs(O[0,0]*O[1,0])
    stat=[static_kernel(float(ma),0) for ma in m]
    cos_bounds=[oscillatory_bound(float(ma),t) for ma in m]
    As_lower=gs*mix*(stat[0]-stat[1]-sum(cos_bounds))
    As_upper=2*gs*mix*sum(stat);Ar_upper=2*gr*mix*sum(stat)
    variance_upper=I0/(2*m.min())
    assert 0<As_lower<=As_upper and s*max(As_upper,Ar_upper)<np.pi/2
    higgs_lower=math.exp(-s*s*variance_upper/2)*math.sin(s*As_lower)*math.cos(s*Ar_upper)
    kernel_pp=[sum(O[0,a]**2*static_kernel(float(ma),float(r)) for a,ma in enumerate(m)) for r in radii]
    theta_lo=gs*gr*(t*np.array(kernel_pp)-I3)
    theta_hi=gs*gr*(t*np.array(kernel_pp)+I3)
    assert min(theta_lo)>0 and max(theta_hi)<np.pi/4
    record_lower=.25*math.exp(-4*gr*gr*I3)*float(np.sin(2*theta_lo).sum())
    theta_abs=gs*gr*(t*I0/m.min()**2+I3)
    vis_lower=math.exp(-4*gr*gr*I3)*math.cos(2*theta_abs)
    gravity_lower=.5*vis_lower*math.sin(kappas.min()*t*(u[0]-u[1]))
    bounds=dict(state_error_G0=float(error0),state_error_GN=float(error),
        higgs_static_contrast_lower=float(higgs_lower),higgs_moving_contrast_lower=float(higgs_lower-2*error),
        record_static_contrast_lower=float(record_lower),record_moving_contrast_lower=float(record_lower-2*error),
        gravity_moving_contrast_lower=float(gravity_lower-error-error0),
        higgs_source_mean_lower=float(As_lower),higgs_mean_absolute_upper=float(As_upper),
        oscillatory_integral_absolute_bounds=cos_bounds,variance_upper=float(variance_upper))
    assert min(bounds[k] for k in ('higgs_moving_contrast_lower','record_moving_contrast_lower','gravity_moving_contrast_lower'))>0
    # Off-diagonal scalar response is the same inverse Hessian: no second portal fit.
    ktest=.7;direct=np.linalg.inv(ktest*ktest*np.eye(2)+K)
    diagonalized=O@np.diag(1/(ktest*ktest+eig))@O.T
    propagator_error=norm(direct-diagonalized)
    assert propagator_error<1e-14 and abs(direct[0,1])>0
    energy_error=float(np.max(abs(energy.sum(axis=1))))
    assert energy_error<1e-12
    # Audit the old measurement wording, not just a new set of passing numerics.
    incompatible=norm(np.diag([1.,0.])@ql-ql@np.diag([1.,0.]))
    assert abs(incompatible-.5)<1e-15
    return dict(parameters=dict(duration=t,lambda_H=lam,v=v,portal_alpha=alpha,mu=mu,
            mass_matrix=K.tolist(),mass_squared=eig.tolist(),radii=radii.tolist(),higgs_read_strength=s,
            inherited_gravity_and_moving_parameters=params),
        common_process=dict(outputs=rows,joint_outcome_order='path_bit,record_bit,higgs_bit; lexicographic 000..111',
            joint_probability_tables=tables,record_probability_contrast=record_gap,higgs_probability_contrast=higgs_gap,
            gravity_path_probability_contrast=gravity_contrast,gram_minimum_eigenvalue=float(np.linalg.eigvalsh(gram).min()),
            higgs_instrument_gram_minimum_eigenvalues=instrument_eigs,closed_form_comparison_error=float(formula_error),
            field_energy_plus_interaction_error=energy_error,portal_propagator_identity_error=propagator_error,
            off_diagonal_static_propagator=float(direct[0,1])),
        finite_moving_control=bounds,
        old_945_reading_correction=dict(path0_vs_coherent_path_commutator_norm=incompatible,
            old_individual_probabilities_preserved=True,three_current_effects_commute=True))

def run():
    out=compute(256);check=compute(320)
    quad=max(abs(out['common_process'][k]-check['common_process'][k]) for k in
        ('record_probability_contrast','higgs_probability_contrast','gravity_path_probability_contrast'))
    assert quad<1e-11
    files=[Path(__file__),STAGE/'944/finite_field_communication.py',
           STAGE/'945/joint_record_newton_results.json',STAGE/'929/joint_protocol_selection_results.json']
    return dict(round=946,date='2026-10-07',all_scientific_checks_passed=True,**out,
        quadrature_comparison_error=float(quad),
        scope=dict(gauge_invariant_parent_potential_is_physical_input=True,
            full_SM_representations_and_parameters_not_derived=True,
            added_scalar_has_no_chiral_gauge_anomaly_contribution=True,
            exact_results_are_for_quadratic_portal_effective_H=True,
            nonlinear_SM_and_full_Einstein_matching_certified=False,
            ideal_Higgs_effect_access_is_declared_input=True,
            full_six_protocol_organization_completed=False,full_goal_completed=False),
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
    else:compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
