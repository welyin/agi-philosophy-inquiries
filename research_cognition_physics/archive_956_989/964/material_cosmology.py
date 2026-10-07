"""964: the SAME 958 material in an explicit homogeneous mean-field branch.
No microscopic gravity derivation, optical apparatus or full SM matching.
The photon Gibbs modes are untruncated; the finite menu is not a UV cutoff.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/"material_cosmology_results.json"

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text("utf-8-sig"))
def entropy(rho):
    ev=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert ev.min()>-2e-12
    ev=ev[ev>1e-14]
    return float(-np.sum(ev*np.log(ev)))

def run():
    spec=importlib.util.spec_from_file_location("old956",STAGE/"956/native_material_interface.py")
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    oldres=read(STAGE/"958/capacitive_material_write_results.json")
    m=old.material(1.,.1);h=m["H"];W=m["W"];J=m["J"];kap=.2
    cs=[old.annihilate(i) for i in range(4)]
    ns=[c.T@c for c in cs]
    sec=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    Q=sec.T@(ns[0]+ns[1]-ns[2]-ns[3])@sec/2
    HA=np.kron(h,np.eye(6));HB=np.kron(np.eye(6),h)
    V=kap*np.kron(Q,Q);H=HA+HB+V
    values,vectors=np.linalg.eigh(H)
    chi=oldres["exact_interface"]["conditional_energy_chi"]
    T=oldres["exact_interface"]["controlled_phase_time"]
    assert abs(-2*J-values[0]-chi)<1e-13
    e16=np.eye(16)
    code=np.stack([(e16[5]-e16[6]-e16[9]+e16[10])/2,
        (2*e16[3]+2*e16[12]-e16[5]-e16[6]-e16[9]-e16[10])/np.sqrt(12)],axis=1)
    local=np.kron(W,np.eye(4))@code
    plus=np.array([1.,1.])/np.sqrt(2)
    inputs=[np.kron(local[:,a],local@plus).astype(complex) for a in (0,1)]
    def act(op,x):
        z=x.reshape(6,4,6,4).transpose(0,2,1,3)
        z=(op@z.reshape(36,16)).reshape(6,6,4,4)
        return z.transpose(0,2,1,3).reshape(576)
    def exp(op,x): return float(np.vdot(x,act(op,x)).real)
    def evolve(t,x):
        u=(vectors*np.exp(-1j*t*values))@vectors.conj().T
        return act(u,x)
    mass=100.
    ebar=sum(exp(H,x) for x in inputs)/2
    M=mass+ebar
    var=sum(exp(H@H,x) for x in inputs)/2-ebar**2
    assert abs(ebar+J)<1e-13 and np.min(values)+mass>0
    material_rows=[];energy_err=0.;nonselective_spectrum_err=0.
    for multiple in (0.,.5,1.,2.):
        time=multiple*T
        states=[evolve(time,x) for x in inputs]
        rb=[x.reshape(24,24).T@x.reshape(24,24).conj() for x in states]
        # rho_B[b,d] = sum_a psi[a,b] conj(psi[a,d]).
        rmean=(rb[0]+rb[1])/2
        holevo=entropy(rmean)-(entropy(rb[0])+entropy(rb[1]))/2
        target=(np.exp(1j*J*time)*local[:,0]+local[:,1])/np.sqrt(2)
        probs=[float(np.vdot(target,r@target).real) for r in rb]
        joint=np.array([[(1-p)/2,p/2] for p in probs])
        cols=joint.sum(axis=0)
        MI=sum(float(joint[a,b]*np.log(joint[a,b]/(.5*cols[b])))
               for a in (0,1) for b in (0,1) if joint[a,b]>1e-14)
        energy_err=max(energy_err,max(abs(exp(H,x)-exp(H,y)) for x,y in zip(states,inputs)))
        # Nonzero eigenvalues of 1/2 sum |psi_s><psi_s| from its 2x2 Gram matrix.
        gram=np.array([[np.vdot(x,y)/2 for y in states] for x in states])
        nonselective_spectrum_err=max(nonselective_spectrum_err,
                                     float(np.max(abs(np.linalg.eigvalsh(gram)-.5))))
        material_rows.append(dict(time_over_T=multiple,receiver_probabilities=probs,
            receiver_entropy=entropy(rmean),receiver_holevo=holevo,
            classical_read_mutual_information=MI,
            mean_full_energy=sum(exp(H,x) for x in states)/2))
    atT=material_rows[2]
    contrast=atT["receiver_probabilities"][1]-atT["receiver_probabilities"][0]
    assert abs(contrast-oldres["exact_interface"]["finite_time_rows"][1]["receiver_contrast"])<1e-10
    assert atT["receiver_holevo"]>.69 and material_rows[-1]["receiver_holevo"]<.001
    assert energy_err<1e-12 and nonselective_spectrum_err<1e-12

    # One fixed set of 12 transverse modes: +-axes, two polarizations.
    # Their equal diagonal occupation makes the *mean* stress isotropic.
    k=J; modes=12;R=100.;nbar=R/(modes*k)
    theta0=k/math.log1p(1/nbar)
    Sgamma=modes*((1+nbar)*math.log1p(nbar)-nbar*math.log(nbar))
    V0=(2*math.pi/k)**3
    def F(a): return 2*(M*a-2*R)*math.sqrt(M*a+R)/(3*M*M)
    sqrtmu=T/(F(2)-F(1));mu=sqrtmu**2
    Hinit=math.sqrt((M+R)/mu)
    Mp2=mu/(3*V0)
    def time(a): return sqrtmu*(F(a)-F(1))
    def conformal(a): return 2*sqrtmu/M*(math.sqrt(M*a+R)-math.sqrt(M+R))
    def scale(t):
        lo,hi=1.,4.
        while time(hi)<t: hi*=2
        for _ in range(65):
            mid=(lo+hi)/2
            if time(mid)<t: lo=mid
            else: hi=mid
        return (lo+hi)/2
    # Independent Raychaudhuri evolution. Constraint is not reset per step.
    def rhs(y):
        a,B,eta,work=y
        return np.array([a*B,
            -3*T*T/(2*mu)*(M/a**3+4*R/(3*a**4)),1/a,R*B/a])
    def solve(n):
        y=np.array([1.,T*Hinit,0.,0.]);step=1/n;mx=0.
        for _ in range(n):
            k1=rhs(y);k2=rhs(y+step*k1/2);k3=rhs(y+step*k2/2);k4=rhs(y+step*k3)
            y+=step*(k1+2*k2+2*k3+k4)/6
            a,B,_,_=y
            mx=max(mx,abs(B*B-T*T/mu*(M/a**3+R/a**4)))
        return y,mx
    coarse,errc=solve(600);fine,errf=solve(1200)
    exact=np.array([2.,T*math.sqrt((2*M+R)/(mu*16)),conformal(2)/T,R/2])
    integerr=float(np.max(abs(fine-exact)))
    stepdiff=float(np.max(abs(fine-coarse)))
    assert integerr<2e-9 and stepdiff<3e-8 and errf<2e-10
    rows=[]
    for a in (1.,1.25,1.5,2.):
        hub=math.sqrt((M*a+R)/(mu*a**4))
        rho=(M/a**3+R/a**4)/V0;P=R/(3*V0*a**4)
        rhodot=-hub*(3*M/a**3+4*R/a**4)/V0
        conservation=abs(rhodot+3*hub*(rho+P))
        pa=-6*Mp2*V0*a*a*hub
        constraint=-pa*pa/(12*Mp2*V0*a)+M+R/a
        theta=theta0/a
        current_n=1/math.expm1((k/a)/theta)
        assert abs(current_n-nbar)<1e-10
        # S = beta E + log Z, evaluated anew at the current a and temperature.
        Scheck=modes*((k/a)/theta*current_n-math.log1p(-math.exp(-(k/a)/theta)))
        assert abs(Scheck-Sgamma)<2e-10
        rows.append(dict(a=a,time_over_T=time(a)/T,hubble=hub,
            photon_to_material_gap_ratio=(k/a)/J,temperature=theta,
            total_comoving_matter_energy=M+R/a,
            accumulated_pressure_work=R*(1-1/a),photon_entropy=Scheck,
            continuity_residual=conservation,Hamiltonian_constraint_residual=abs(constraint)))
    # Deliberately wrong one-rule rescaling, including the bound internal H.
    # It preserves photon / internal-gap ratio, unlike this common model.
    wrong_gap_ratio=(k/2)/(J/2)
    correct_gap_ratio=(k/2)/J
    assert wrong_gap_ratio==1 and correct_gap_ratio==.5
    # An explicit omitted finite-multipole *input*, not a parent error proof.
    # If ||delta H|| <= q Hubble^2, Duhamel gives eps <= q H0 log(a_f).
    q=100.;conditional_epsilon=q*Hinit*math.log(2)
    record_lower=oldres["exact_interface"]["analytic_receiver_contrast_lower"]-2*conditional_epsilon
    assert record_lower>.94
    sourcefiles=[Path(__file__),STAGE/"956/native_material_interface.py",
        STAGE/"958/capacitive_material_write.py",STAGE/"958/capacitive_material_write_results.json",
        STAGE/"941/composite_operation_bridge_results.json",STAGE/"963/common_scale_certificate_results.json"]
    return dict(round=964,date="2026-10-07",all_scientific_checks_passed=True,
        kind="same-material homogeneous semiclassical coexistence; not a full common quantum model",
        material=dict(parameters=dict(U=1.,v=.1,kappa=kap,rest_offset=mass),
            active_dimension=36,total_dimension=576,J=J,chi=chi,T=T,
            balanced_prior_internal_energy=ebar,complete_mean_mass=M,
            mass_variance=var,minimum_mass=mass+float(values.min()),
            energy_conservation_error=energy_err,
            nonselective_spectrum_error=nonselective_spectrum_err,
            conditional_read_contrast=contrast,rows=material_rows),
        cosmology=dict(photon_modes=modes,photon_comoving_energy=R,
            mean_occupation=nbar,initial_temperature=theta0,coordinate_volume=V0,
            reduced_planck_mass_squared=Mp2,mu=mu,initial_hubble=Hinit,
            H0_over_material_gap=Hinit/J,conformal_elapsed_over_T=conformal(2)/T,
            scale_at_two_T=scale(2*T),rows=rows,
            independent_Raychaudhuri=dict(fine_endpoint=fine.tolist(),analytic_endpoint=exact.tolist(),
                max_endpoint_error=integerr,step_difference=stepdiff,
                max_dimensionless_constraint_residual=errf)),
        thermal=dict(photon_entropy=Sgamma,full_material_entropy=math.log(2),
            photon_temperature_ratio_at_two=.5,photon_entropy_production=0.,
            thermalization_or_irreversible_arrow_generated=False),
        counterchecks=dict(same_rule_rescaling_gap_ratio=wrong_gap_ratio,
            actual_ratio=correct_gap_ratio,
            missing_pressure_work_if_comoving_energy_declared_constant=R/2),
        conditional_tidal_ledger=dict(input_q=q,
            assumption="operator norm delta H <= q Hubble^2 along fixed declared background",
            local_trace_distance_upper=conditional_epsilon,
            read_contrast_lower_if_assumption_holds=record_lower,
            coefficient_matched_from_material_parent=False,
            perturbed_background_error_certified=False),
        scope=dict(same_958_material_and_readout=True,Einstein_FLRW_and_mass_coupling_are_inputs=True,
            uniform_mean_source_and_pressure_work_closed=True,
            material_clock_and_radiation_do_not_share_all_rate_rescaling=True,
            thermal_initial_state_and_low_correlation_preparation_are_inputs=True,
            semiclassical_background_is_not_full_linear_quantum_process=True,
            gravitational_source_fluctuations_certified=False,
            microscopic_binding_or_optical_readout_derived=False,
            full_SM_and_GR_common_recovery=False,full_goal_completed=False,
            no_new_device_optimization=True),
        references=["https://arxiv.org/html/0810.2712","https://arxiv.org/html/1502.00971",
                    "https://arxiv.org/abs/0802.0658"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sourcefiles})

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    if args.write: assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write("\n")
    else:
        before=read(TARGET)
        assert before["source_hashes"]==out["source_hashes"] and before["scope"]==out["scope"]
        assert abs(before["material"]["conditional_read_contrast"]-out["material"]["conditional_read_contrast"])<1e-11
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))

