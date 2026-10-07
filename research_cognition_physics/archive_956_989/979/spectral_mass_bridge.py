"""979: source-aware transport of static spectral simulation.
The finite matrix is a contract check, not a generated exchange simulator.
Continuum statements are proved in the note, not inferred from discretization.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent; STAGE=HERE.parent; ROOT=HERE.parents[2]
OUT=HERE/'spectral_mass_bridge_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def op(a):return float(np.linalg.norm(a,2))
def kron(*xs):
    a=np.ones((1,1))
    for x in xs:a=np.kron(a,x)
    return a
def evol(h,t):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T
def run():
    X=np.array([[0.,1.],[1.,0.]])
    Z=np.diag([1.,-1.]); K=.5*np.eye(2)+.2*X+.3*Z
    eps=1e-4; defect=eps*(.6*X+.8*Z)
    angle=2e-4; c=1.7; bare_mass=10.; mu=13.; G=.005; core=.4
    R=np.eye(3);R[np.ix_([0,2],[0,2])]=[[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]]
    W=np.eye(3)[:,:2]; V=R[:,:2]; eta=op(V-W)
    internal=np.zeros((3,3));internal[:2,:2]=K+c*np.eye(2)+defect;internal[2,2]=70+c
    Hs=R@internal@R.T
    Ms=bare_mass*np.eye(3)+Hs
    Mt=(bare_mass+c)*np.eye(2)+K
    Mc=Mt+defect
    assert op(V.T@Hs@V-(K+c*np.eye(2)+defect))<1e-12
    assert op(Ms@V-V@Mc)<1e-12
    assert op(K@defect-defect@K)>1e-6
    # A periodic finite position model checks the algebra only.
    n=5;shift=np.roll(np.eye(n),1,axis=0);L=2*np.eye(n)-shift-shift.T
    ix=np.arange(n);dist=np.abs(ix[:,None]-ix[None,:]);dist=np.minimum(dist,n-dist)
    U=np.diag((dist.astype(float)**2+core**2).ravel()**-.5)
    P2=kron(L,np.eye(n));Y2=kron(np.eye(n),L);Ipos=np.eye(n*n)
    def hamiltonian(M):
        return kron(M,Ipos)+kron(np.linalg.inv(M)/2,P2)+kron(np.eye(len(M)),Y2/(2*mu))-G*mu*kron(M,U)
    hs=hamiltonian(Ms);ht=hamiltonian(Mt);hc=hamiltonian(Mc)
    ve=kron(V,Ipos);we=kron(W,Ipos)
    exact_intertwining=op(hs@ve-ve@hc)
    min_mass=min(np.linalg.eigvalsh(Mt).min(),np.linalg.eigvalsh(Mc).min())
    rate_bound=eps*(1+op(P2)/(2*min_mass**2))
    actual_rate=op(hc-ht)
    assert actual_rate<=rate_bound+1e-12
    rows=[]
    for t in (0.,.5,2.,5.):
        actual=op(evol(hs,t)@we-we@evol(ht,t))
        bound=2*eta+t*rate_bound
        assert actual<=bound+1e-11
        rows.append(dict(time=t,all_input_isometry_difference=actual,analytic_uniform_bound=bound))
    total_shift=kron(shift,shift)
    momentum_symmetry=op(hs@kron(np.eye(3),total_shift)-kron(np.eye(3),total_shift)@hs)
    # Actual Newton source and equal/opposite force use the same full mass.
    x=np.array([.2,-.1,.3]);y=np.array([1.1,.4,-.2]);r=y-x;den=float(r@r)+core**2
    force_coeff=G*mu*r/den**1.5
    force_error=max(op(V.T@(f*Ms)@V-f*Mt) for f in force_coeff)
    force_bound=eps*G*mu*2/(3*math.sqrt(3)*core**2)
    assert force_error<=force_bound+1e-12
    source_error=op(V.T@(-mu*Ms/math.sqrt(den))@V+mu*Mt/math.sqrt(den))
    assert source_error<=mu/core*eps+1e-10
    def potential(xx,yy):return -G*mu*Ms/math.sqrt(float((yy-xx)@(yy-xx))+core**2)
    step=1e-4;fd=[];equal_opposite=[]
    for i in range(3):
        e=np.eye(3)[i]*step
        dx=(potential(x-2*e,y)-8*potential(x-e,y)+8*potential(x+e,y)-potential(x+2*e,y))/(12*step)
        dy=(potential(x,y-2*e)-8*potential(x,y-e)+8*potential(x,y+e)-potential(x,y+2*e))/(12*step)
        fd.append(op(-dx-force_coeff[i]*Ms));equal_opposite.append(op(dx+dy))
    assert max(fd)<1e-10 and max(equal_opposite)<1e-10 and exact_intertwining<1e-11
    assert momentum_symmetry<1e-11
    # Same positive comparison dynamics and actual initial p^2 moment in R^3.
    # K is the scaled full 929 protocol H, 0 <= K <= pi; no simulator couplings generated here.
    T=30.;eta_app=1e-5;epsilon_app=1e-6;mass0=20.;mstar=19.;mmax=24.;E2=25.
    sigma=100.;mu_app=20.;G_app=1e-30;a_app=1.
    kinetic_rate=math.sqrt(15)/(8*sigma**2)*(1/mass0+1/mu_app)
    actual_E2_upper=mass0+math.pi+kinetic_rate
    assert actual_E2_upper<E2 and mass0+math.pi<mmax
    B=1+mmax*(E2+mmax)/mstar**2
    simulation_bound=2*eta_app+T*epsilon_app*B
    flat_protocol_bound=T*(kinetic_rate+G_app*mu_app/a_app*(mass0+math.pi))
    total=simulation_bound+flat_protocol_bound
    assert total<.0003
    # The baseline cannot be omitted when it multiplies a relative lapse.
    physical_base=20.;shift_energy=2.;phase=.5*math.pi
    omitted=(1+math.cos(phase*physical_base))/2
    included=(1+math.cos(phase*(physical_base+shift_energy)))/2
    assert abs(omitted-included)>1-1e-12
    # Concrete specialization of the known norm-versus-energy-moment warning.
    tails=[]
    for theta in (.1,.03,.01):
        Lambda=theta**-4;p=math.sin(theta)**2
        tails.append(dict(angle=theta,high_energy=Lambda,high_energy_probability=p,
            all_time_trace_distance_upper=math.sin(2*theta),
            excess_mass_mean=Lambda*p,excess_mass_second_moment=(40*Lambda+Lambda**2)*p))
    assert tails[-1]['all_time_trace_distance_upper']<.021 and tails[-1]['excess_mass_mean']>9999
    files=[Path(__file__),HERE/'drafts/spectral_mass_decision.md',
        STAGE.parent/'archive_429_466/research_note_436.md',
        STAGE.parent/'archive_585_628/research_note_625.md',
        STAGE/'929/joint_protocol_selection.py',STAGE/'929/joint_protocol_selection_results.json',
        STAGE/'research_note_935.md',STAGE/'978/moving_complete_source.py',
        STAGE/'978/moving_complete_source_results.json']
    return dict(round=979,date='2026-10-07',all_scientific_checks_passed=True,
        finite_contract=dict(position_dimension=n*n,internal_dimensions=[2,3],eta=eta,epsilon=eps,
            actual_internal_error=op(defect),noncommuting_error=op(K@defect-defect@K),
            low_mass_intertwining_error=exact_intertwining,inverse_mass_error=op(np.linalg.inv(Mc)-np.linalg.inv(Mt)),
            minimum_mass=float(min_mass),actual_full_H_error=actual_rate,full_H_error_bound=rate_bound,
            weak_potential=G*mu/core,rows=rows,translation_symmetry_error=momentum_symmetry,
            exact_compressed_force_error=force_error,exact_compressed_force_bound=force_bound,
            source_derivative_error=source_error,force_finite_difference_error=max(fd),
            equal_opposite_force_difference=max(equal_opposite)),
        protocol_transport=dict(target='929 complete H scaled by pi/30, including clock and faults',
            unknown_input_and_passive_reference=True,physical_simulator_generated=False,
            source_rule='complete actual simulator energy with matched rest baseline',
            time=T,eta=eta_app,epsilon=epsilon_app,mass0=mass0,mass_lower=mstar,mass_upper=mmax,
            initial_H_minus_mu_norm_bound=E2,gaussian_initial_bound=actual_E2_upper,
            sigma=sigma,mu=mu_app,G=G_app,soft_core=a_app,continuum_coefficient=B,
            moving_simulation_error=simulation_bound,moving_target_to_flat_protocol=flat_protocol_bound,
            full_protocol_output_bound=total,uniform_time_interval=True),
        resource_boundaries=dict(baseline_omission_probability_gap=abs(omitted-included),tail_examples=tails,
            finite_spectral_simulation_does_not_bound_bare_preparation_energy=True,
            actual_simulator_rest_mass_matching_is_input=True),
        scope=dict(existing_exchange_theorem_reused=True,continuum_positions_in_analytic_proof=True,
            finite_position_check_not_continuum_convergence=True,exchange_coupling_list_generated=False,
            native_976_charge_model_identified_with_exchange_simulator=False,
            all_SM_matching_completed=False,full_GR_completed=False,macro_arrow_completed=False,
            resolved_gravity_signal_claimed_in_protocol_example=False,full_goal_completed=False),
        references=['https://arxiv.org/html/1701.05182v4','https://arxiv.org/abs/1502.00971'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-6,abs_tol=1e-10),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not OUT.exists()
    out=run()
    if a.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(OUT.read_text('utf-8-sig')))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
