"""955: preserve a joint finite protocol at the selected portal recovery point.

New analytic ingredient: an observable-specific Heisenberg bound for the
Higgs radial effect. Existing Gaussian/protocol structure is explicitly reused.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"joint_effective_window_results.json"
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.])
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def norm(a):return float(np.linalg.norm(a,2))

def compute(n):
    field=load("field944",STAGE/"944/finite_field_communication.py")
    portal=load("portal946",STAGE/"946/portal_common_process.py")
    protocol=load("protocol929",STAGE/"929/joint_protocol_selection.py")
    old=read(STAGE/"947/protocol_field_transport_results.json")
    base=read(STAGE/"946/portal_common_process_results.json")
    recovery=read(STAGE/"954/portal_recovery_window_results.json")
    pars=base["parameters"]["inherited_gravity_and_moving_parameters"]
    zeta=recovery["declared_coefficient_window"]["chosen_zeta"]
    K=np.array([[1+.25*zeta*zeta,zeta],[zeta,4.]])
    eigen,O=np.linalg.eigh(K);m=np.sqrt(eigen);mix=abs(O[0,0]*O[1,0])
    T=30.;strength=8.;gs=gr=.5;half_window=1e-5;hmax=np.pi
    mass=pars["base_mass"];mr=pars["receiver_mass"];G=pars["G"]
    delta=pars["position_std"];core=pars["soft_core"];radii=np.array(pars["radii"])
    k,w=field.nodes(n,0.,12.);angle,wa=field.nodes(48,-1.,1.)
    measure=(w*k*k/(4*np.pi**2))[:,None]*wa[None,:]
    radial=w*k*k*np.exp(-k*k)/(2*np.pi**2)
    sectors=[(l,b,y) for l in (0,1) for b in (1.,-1.) for y in (-1.,1.)]
    amplitudes=[];phases=np.zeros(8);read_vectors=[];energy=np.zeros((8,2))
    theta=np.zeros(2);Drr=0.;Ddiff=0.;variance=0.;meanS=0.;meanR=np.zeros(2)
    for a,ma in enumerate(m):
        om=np.sqrt(k*k+ma*ma)[:,None]
        mode=np.exp(-k*k/2)[:,None]/np.sqrt(2*om)*np.ones((1,len(angle)))
        displacements=[]
        for j,(l,b,y) in enumerate(sectors):
            source=O[0,a]*mode*(gs*b+gr*y*np.exp(-1j*k[:,None]*radii[l]*angle))
            disp=source*(np.exp(-1j*om*T)-1)/om
            displacements.append((disp*np.sqrt(measure)).reshape(-1))
            phases[j]+=float(np.sum(measure*abs(source)**2*(om*T-np.sin(om*T))/om**2))
            energy[j,0]+=float(np.sum(measure*om*abs(disp)**2))
            energy[j,1]+=float(2*np.real(np.sum(measure*source.conj()*disp)))
        amplitudes.append(np.array(displacements))
        read_vectors.append((O[1,a]*mode*np.sqrt(measure)).reshape(-1))
        om=om[:,0]
        theta+=gs*gr*O[0,a]**2*np.array([
            np.sum(radial*np.sinc(k*r/np.pi)*(om*T-np.sin(om*T))/om**3) for r in radii])
        Drr+=gr*gr*O[0,a]**2*np.sum(radial*(1-np.cos(om*T))/om**3)
        Ddiff+=gr*gr*O[0,a]**2*np.sum(radial*(1-np.sinc(k/np.pi))*(1-np.cos(om*T))/om**3)
        variance+=O[1,a]**2*np.sum(radial/(2*om))
        meanS-=gs*O[0,a]*O[1,a]*np.sum(radial*(1-np.cos(om*T))/om**2)
        meanR-=gr*O[0,a]*O[1,a]*np.array([
            np.sum(radial*np.sinc(k*r/np.pi)*(1-np.cos(om*T))/om**2) for r in radii])
    amps=np.concatenate(amplitudes,axis=1);d=np.concatenate(read_vectors)
    inner=amps@amps.conj().T;diag=inner.diagonal().real
    F=np.exp(1j*(phases[:,None]-phases[None,:]+inner.imag)
        -.5*(diag[:,None]+diag[None,:]-2*inner.real))
    invr=1/np.sqrt(radii*radii+core*core);kap=G*mass*mr
    phaseG=np.array([kap*T*invr[l] for l,b,y in sectors])
    gram=F*np.exp(1j*(phaseG[:,None]-phaseG[None,:]))
    z=amps@d
    sine=gram*np.exp(-strength*strength*variance/2)*np.sin(strength*(z[:,None]+z.conj()[None,:]))
    eig_min=min(float(np.linalg.eigvalsh((gram+s*sine)/2).min()) for s in (-1.,1.))
    assert eig_min>-1e-12 and norm(gram-gram.conj().T)<1e-12
    y,vy=np.linalg.eigh(Y);rot=np.kron(vy,vy);plus=np.array([1,1],complex)/np.sqrt(2)
    prep=rot.conj().T@np.kron(plus,plus);rho_ra=np.outer(prep,prep.conj())
    charges=(y[:,None]+y[None,:]).reshape(-1);rho_ra*=np.isclose(charges[:,None],charges[None,:])
    rel=rot.conj().T@((np.kron(Z,X)-np.kron(X,Z))/2)@rot
    effectR=(np.eye(4)+rel)/2;effectG=(I-Y)/2
    tables=[]
    for bi in (0,1):
        rho_s=np.zeros((2,2));rho_s[bi,bi]=1
        rho=np.kron(np.outer(plus,plus),np.kron(rho_s,rho_ra)).reshape(8,2,8,2)
        def density(kernel):return (rho*kernel[:,None,:,None]).reshape(16,16)
        table=[]
        for ig in (0,1):
            for ir in (0,1):
                for ih in (0,1):
                    eg=effectG if ig else I-effectG
                    er=effectR if ir else np.eye(4)-effectR
                    kernel=(gram+sine)/2 if ih else (gram-sine)/2
                    table.append(float(np.trace(density(kernel)@np.kron(eg,np.kron(I,er))).real))
        assert min(table)>0 and abs(sum(table)-1)<1e-12
        tables.append(table)
    tables=np.array(tables)
    record_gap=float((tables[0]-tables[1])[[2,3,6,7]].sum())
    higgs_gap=float((tables[0]-tables[1])[[1,3,5,7]].sum())
    analyticR=.25*np.exp(-2*Drr)*np.sin(2*theta).sum()
    analyticH=np.exp(-strength*strength*variance/2)*np.sin(strength*meanS)*np.cos(strength*meanR).mean()
    analyticG=.5*np.exp(-Ddiff)*np.cos(theta[0]-theta[1])*np.sin(kap*T*(invr[0]-invr[1]))
    assert abs(record_gap-analyticR)<1e-12 and abs(higgs_gap-analyticH)<1e-12
    # Full original two-event instrument, arbitrary Q and passive reference.
    ks=protocol.kraus(np.pi/2);bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    choi_error=0.;input_joint=[];complete=np.zeros((2,2),complex)
    for x in (0,1):
        for ybit in (0,1):
            op=ks[ybit]@ks[x];v= np.kron(op,I)@bell;block=np.outer(v,v.conj())
            complete+=op.conj().T@op
            choi_error=max(choi_error,norm(sum(tables[x,j]*block for j in range(8))-block))
    for psi in (protocol.ry(-np.pi/2)[:,q] for q in (0,1)):
        events=np.array([[np.linalg.norm(ks[ybit]@ks[x]@psi)**2 for ybit in (0,1)] for x in (0,1)])
        input_joint.append(np.einsum("xy,xz->xyz",events,tables).sum(axis=(0,1)))
    input_joint=np.array(input_joint)
    actualR=float((input_joint[0]-input_joint[1])[[2,3,6,7]].sum())
    actualH=float((input_joint[0]-input_joint[1])[[1,3,5,7]].sum())
    assert abs(actualR-.6*record_gap)<1e-12 and abs(actualH-.6*higgs_gap)<1e-12
    assert choi_error<1e-12 and norm(complete-I)<1e-12
    # Re-evaluate the old analytic norm bounds at the NEW K.
    I0=1/(8*np.pi**1.5);I1=1/(4*np.pi**2)
    I3=I0/m[0]**3;grad=3/(32*np.pi**1.5*m[0]);A2=2*I3
    kinetic=T*2*np.sqrt(15)/(8*mass*delta*delta)
    position=T*delta*np.sqrt(grad)*np.sqrt(1+4*A2)
    oldG=read(STAGE/"945/joint_record_newton_results.json")
    gravity_position=oldG["finite_moving_control"]["added_gravity_position_bound"]
    moving0=kinetic+position;moving=moving0+gravity_position
    redshift=T*G*mr*hmax*float(invr.max())
    initial_H_norm=hmax+2*np.sqrt(15)/(8*mass*delta*delta)+math.sqrt(I0/(2*m[0]))+G*mr*(mass+hmax)/core
    alltime=half_window*initial_H_norm
    all_error=moving+redshift+alltime
    static_pp=[sum(O[0,a]**2*portal.static_kernel(float(ma),float(r)) for a,ma in enumerate(m)) for r in radii]
    theta_lo=gs*gr*(T*np.array(static_pp)-I3)
    theta_hi=gs*gr*(T*np.array(static_pp)+I3)
    assert theta_lo.min()>0 and theta_hi.max()<np.pi/4
    record_lower=.25*math.exp(-I3)*float(np.sin(2*theta_lo).sum())
    theta_abs=.25*(T*I0/m[0]**2+I3)
    gravity_lower=.5*math.exp(-I3)*math.cos(2*theta_abs)*math.sin(kap*T*(invr[0]-invr[1]))
    stat=[portal.static_kernel(float(ma),0.) for ma in m]
    cos_bound=[portal.oscillatory_bound(float(ma),T) for ma in m]
    As=gs*mix*(stat[0]-stat[1]-sum(cos_bound))
    Ar=2*gr*mix*sum(stat)
    assert As>0 and strength*Ar<np.pi/2
    higgs_lower=math.exp(-strength*strength*I0/(4*m[0]))*math.sin(strength*As)*math.cos(strength*Ar)
    # NEW: energy bound on actual moving state; no prescribed-trajectory
    # substitution. H >= Tkin + Hfield/2 - C - Vmax and h>=0.
    Eupper=hmax+3/(4*mass*delta*delta)
    Cenergy=I0/m[0]**2
    Vmax=G*mr*(mass+hmax)/core
    L=Eupper+Cenergy+Vmax
    field_energy_upper=2*L
    speed_rms_bound=math.sqrt(2*L/mass)
    mixed_gradient_bound=mix*(1/m[0]+1/m[1])*I1
    position_integral=T*math.sqrt(3)*delta+.5*T*T*speed_rms_bound
    weak_position=strength/2*(abs(gs)+abs(gr))*mixed_gradient_bound*position_integral
    momentum_field_norm=math.sqrt(2*I0)*math.sqrt(field_energy_upper)+math.sqrt((I1+m[1]*I0)/2)
    weak_time=half_window*(strength/2*momentum_field_norm+strength*strength*I0/4)
    weak_higgs=weak_position+weak_time
    # This weak bound controls only E_H, not the full cq instrument.
    lowerR=.6*record_lower-2*all_error
    lowerH=.6*higgs_lower-2*weak_higgs
    lowerG=gravity_lower-2*all_error
    rejected=.6*higgs_lower-2*all_error
    assert lowerR>0 and lowerH>0 and lowerG>0 and rejected<0
    # Independent mixed commutator Lipschitz diagnostic at several times.
    max_gradient=0.
    for time in (0.,.7,4.,17.,30.):
        kernel=np.zeros_like(k)
        for a,ma in enumerate(m):
            om=np.sqrt(k*k+ma*ma)
            kernel+=O[0,a]*O[1,a]*np.sin(time*om)/om
        measured=float(np.sum(radial*k*abs(kernel)))
        assert measured <= mixed_gradient_bound*(1+1e-12)
        max_gradient=max(max_gradient,measured)
    return dict(parameters=dict(zeta=zeta,K=K.tolist(),scalar_masses_squared=eigen.tolist(),
        duration=T,higgs_read_strength=strength,readout_half_window=half_window,
        base_mass=mass,position_std=delta,all_preexisting_resources_kept=True),
        fixed_center=dict(joint_tables=tables.tolist(),protocol_joint_tables=input_joint.tolist(),
            record_contrast=actualR,higgs_contrast=actualH,gravity_contrast=float(analyticG),
            instrument_Gram_minimum=eig_min,Choi_recovery_error=choi_error,
            field_energy_plus_interaction_error=float(np.max(abs(energy.sum(axis=1))))),
        bounds=dict(global_cq_instrument_distance=all_error,global_state_movement=moving,
            global_time_window=alltime,internal_redshift=redshift,
            kinetic_energy_upper=L,field_energy_upper=field_energy_upper,
            mixed_retarded_gradient_upper=mixed_gradient_bound,mixed_gradient_diagnostic=max_gradient,
            weak_higgs_position_error=weak_position,weak_higgs_time_window=weak_time,
            weak_higgs_total_error=weak_higgs,higgs_fixed_center_protocol_lower=.6*higgs_lower,
            record_probability_contrast_lower=lowerR,higgs_probability_contrast_lower=lowerH,
            gravity_probability_contrast_lower=lowerG,using_global_bound_for_Higgs_lower=rejected),
        recovery=dict(matter_relative_kernel_bound=.25*zeta*zeta,
            record_W_Q1=abs(.5*.4225*zeta/(10+.25*zeta*zeta)),
            original_954_coefficient_contract_satisfied=True))

def run():
    out=compute(256);second=compute(320)
    discrepancy=max(abs(out["fixed_center"][k]-second["fixed_center"][k]) for k in
        ("record_contrast","higgs_contrast","gravity_contrast"))
    assert discrepancy<1e-11
    files=[Path(__file__),STAGE/"944/finite_field_communication.py",
        STAGE/"946/portal_common_process.py",STAGE/"929/joint_protocol_selection.py",
        STAGE/"945/joint_record_newton_results.json",STAGE/"946/portal_common_process_results.json",
        STAGE/"947/protocol_field_transport_results.json",STAGE/"954/portal_recovery_window_results.json"]
    return dict(round=955,date="2026-10-07",all_scientific_checks_passed=True,**out,
        quadrature_comparison=discrepancy,
        scope=dict(actual_two_protocol_unknown_input_and_passive_reference_transported=True,
            new_quadratic_field_vacuum_declared=True,finite_motion_and_time_window_included=True,
            same_h_B_scalar_K_Newton_process=True,all_three_finite_contrasts_have_analytic_positive_lower_bounds=True,
            full_instrument_bound_not_confused_with_effect_bound=True,
            matter_recovery_is_tree_spacelike_not_full_SM=True,
            same_H_energy_and_total_momentum_conservation_preserved=True,
            complete_covariant_stress_and_quantum_gravity_certified=False,
            all_access_instruments_autonomously_manufactured=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true");args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        previous=read(TARGET)
        assert previous["source_hashes"]==out["source_hashes"] and previous["scope"]==out["scope"]
        for key,val in out["bounds"].items():assert abs(val-previous["bounds"][key])<1e-11,key
    print(json.dumps({k:v for k,v in out.items() if k not in ("source_hashes","fixed_center")},ensure_ascii=False,indent=2))
