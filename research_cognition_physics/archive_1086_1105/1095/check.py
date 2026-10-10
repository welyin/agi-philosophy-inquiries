"""Finite certificate for round1095. Default recomputes without writing."""
from pathlib import Path
from fractions import Fraction as Fr
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent


def split_branch(n, mu, g, c, u0=0.):
    mass=mu*2**(n+1)
    m,u=mass,u0
    leaves=[]
    steps=[]
    released=0.
    x,t=0.,0.
    for k,a in enumerate([g]*n+[c]):
        duration=.5 if k<n else .25
        x+=u*duration
        t+=duration
        energy=m*a*a/2
        leaves.append(dict(mass=m/2,velocity=u-a,birth_time=t,birth_x=x))
        steps.append(dict(parent_mass=m,parent_velocity=u,kick=a,
                          daughter_velocity=u+a,released_energy=energy,
                          event_time=t,event_x=x))
        released+=energy
        m,u=m/2,u+a
    leaves.append(dict(mass=m,velocity=u,birth_time=t,birth_x=x))
    return dict(initial_mass=mass,initial_velocity=u0,core_count=2**(n+1),
                steps=steps,leaves=leaves,released=released)


def calculate():
    groups=[]
    def add(name,condition,**data):
        assert condition,name
        groups.append(dict(name=name,status='PASS',**data))

    errors=[0.,0.,0.]
    for n in range(9):
        for u0 in [-.5,0.,1.25]:
            b=split_branch(n,1.,.75,1.,u0)
            mass=sum(v['mass'] for v in b['leaves'])
            p=sum(v['mass']*v['velocity'] for v in b['leaves'])
            kin=sum(v['mass']*v['velocity']**2/2 for v in b['leaves'])
            errors=[max(errors[0],abs(mass-b['initial_mass'])),
                    max(errors[1],abs(p-b['initial_mass']*u0)),
                    max(errors[2],abs(kin-b['initial_mass']*u0*u0/2-b['released']))]
            assert math.isclose(b['released'],(2**(n+1)-2)*.75**2+1.)
    example=split_branch(3,1.,.75,1.)
    add('all_recoil_branches_and_actual_members',max(errors)<1e-10,
        finite_preparations=27,max_mass_momentum_energy_errors=errors,
        example=example,source_cores_exclude_receiver_and_control=True)

    # Actual velocity update is a declared new primitive. This operator check
    # does not pretend the old fixed-velocity drift Hamiltonian generates it.
    speeds=np.array([-.5,.25,1.])
    mass,a=8.,.75
    fwd,rear=speeds+a,speeds-a
    W=np.zeros((36,6),complex)  # plus velocity, minus velocity, Q+, blank Q-
    pout=np.empty(36)
    hout=np.empty(36)
    for i in range(3):
        for j in range(3):
            for q in range(2):
                for blank in range(2):
                    idx=np.ravel_multi_index((i,j,q,blank),(3,3,2,2))
                    pout[idx]=mass*(fwd[i]+rear[j])/2
                    hout[idx]=mass*(fwd[i]**2+rear[j]**2)/4
        for q in range(2):
            W[np.ravel_multi_index((i,i,q,0),(3,3,2,2)),2*i+q]=1
    pin=np.repeat(mass*speeds,2)
    hin=np.repeat(mass*speeds**2/2+mass*a*a/2,2)
    iso=float(np.linalg.norm(W.conj().T@W-np.eye(6),ord=2))
    pe=float(np.linalg.norm(pout[:,None]*W-W*pin[None,:],ord=2))
    he=float(np.linalg.norm(hout[:,None]*W-W*hin[None,:],ord=2))
    # Unknown qubit/reference tested through its Choi state, with coherent source.
    chi=np.array([1.,1j,-1.])/math.sqrt(3)
    bell=np.array([1.,0.,0.,1.])/math.sqrt(2)
    psi=np.kron(chi,bell).reshape(6,2)
    output=(W@psi).reshape(3,3,2,2,2)
    qr=np.einsum('ijqbr,ijsbt->qrst',output,output.conj()).reshape(4,4)
    expected=np.outer(bell,bell.conj())
    qe=float(np.linalg.norm(qr-expected,ord=2))
    add('split_isometry_full_input_and_conservation',max(iso,pe,he,qe)<1e-12,
        isometry_residual=iso,momentum_intertwiner_residual=pe,
        energy_intertwiner_residual=he,payload_reference_residual=qe,
        rear_payload_is_prepared_blank=True,full_physical_energy_model_certified=False)

    # The new drift velocities change positions, not just stored names.
    source=example['steps'][3]
    source_start_t=example['steps'][2]['event_time']
    source_start_x=example['steps'][2]['event_x']
    source_u=source['parent_velocity']
    emission_dx=source['event_x']-source_start_x
    expected_dx=source_u*(source['event_time']-source_start_t)
    add('actual_preparation_and_emission_worldlines',abs(emission_dx-expected_dx)<1e-12
        and source_u==2.25 and emission_dx==.5625,
        new_input_absolute_time=source_start_t,new_input_absolute_x=source_start_x,
        source_velocity=source_u,emission_displacement_after_input=emission_dx,
        source_preparation_finished_before_new_input=True)

    V,c,d=Fr(3),Fr(1),Fr(1,4)
    n=max(0,math.floor((V-c)/(c-d))+1)
    w=c+n*(c-d)
    delta,eps,D=Fr(1,4),Fr(1,2),Fr(52)
    u=n*(c-d)
    total=delta+(D-u*delta+eps)/w
    h=delta+eps/w
    bound=Fr(1)+D/w
    threshold=V*Fr(1)*w/(w-V)
    assert total<=h+D/w<=bound
    add('finite_complete_delivery_certificate',n==3 and w==Fr(13,4)
        and total==Fr(211,13) and bound==17 and D>threshold
        and D/bound>V and D-u*delta>eps,
        n=n,w=str(w),actual_time=str(total),conservative_time=str(bound),
        complete_average_speed_lower=str(D/bound),threshold_using_h_one=str(threshold),
        core_count=16,split_release_energy=str(Fr(71,8)),
        all_preparation_costs_recorded_but_not_added_to_message_age=True)

    # Reuse round1091's exact encounter integral and payload SWAP.
    nodes,weights=np.polynomial.legendre.leggauss(3)
    integral=float(np.dot(weights,.75*(1-nodes**2)))
    S=np.zeros((4,4))
    for q in range(2):
        for r in range(2): S[2*r+q,2*q+r]=1
    K=np.pi/2*(np.eye(4)-S)
    ev,U=np.linalg.eigh(K)
    gate=(U*np.exp(-1j*ev*integral))@U.conj().T
    swap_err=float(np.linalg.norm(gate-S,ord=2))
    initial=np.zeros(8,complex)  # payload, receiver blank, reference
    initial[0]=initial[5]=1/math.sqrt(2)
    final=(np.kron(gate,np.eye(2))@initial).reshape(2,2,2)
    br=np.einsum('qbr,qcs->brcs',final,final.conj()).reshape(4,4)
    ref_error=float(np.linalg.norm(br-expected,ord=2))
    add('complete_encounter_channel_reused',abs(integral-1)<1e-12
        and swap_err<1e-12 and ref_error<1e-12,
        encounter_integral=integral,full_gate_swap_residual=swap_err,
        receiver_reference_choi_residual=ref_error,
        exact_task=True,receiver_mechanical_backreaction_generated=False)

    eta=.5
    speeds=[]
    value=0.
    for _ in range(20):
        value=eta*(value+1.)
        speeds.append(value+1.)
    completed=[]
    for flight in [4.,10.,100.]:
        distance=52.
        delay=distance*max(0.,1/3.-1/flight)
        completed.append(dict(flight_speed=flight,input_after_processing=delay,
                              completed_speed=distance/(delay+distance/flight)))
    # Use the actual event time, not the looser h=1 sufficient threshold.
    # Here T(D)=(D+c*delta+eps)/w, so exceeding V needs D>9 exactly.
    exact_threshold=V*(c*delta+eps)/(w-V)
    Dmax=8.
    domain_max_speed=Dmax/float((Fr(8)+c*delta+eps)/w)
    add('removing_uniform_role_and_delivery_permissions',max(speeds)<2.
        and all(z['completed_speed']<=3.+1e-12 for z in completed)
        and exact_threshold==9 and Dmax<float(exact_threshold)
        and domain_max_speed<3.,
        proportional_role_inheritance_speed_upper=2.,
        proportional_role_finite_speeds=speeds,
        distance_proportional_processing=completed,
        finite_domain_Dmax=Dmax,exact_threshold_for_N3=float(exact_threshold),
        completed_speed_at_Dmax=domain_max_speed,
        N3_exact_delivery_violates_bound_in_this_domain=False,
        not_a_claim_about_other_preparations_or_larger_N=True)

    Dmin,Tmax,Vmax=Fr(519,10),Fr(1701,100),Fr(301,100)
    gap=Dmin-Vmax*Tmax
    add('finite_calibration_margin_with_matching_task_scope',gap>0,
        minimum_distance=str(Dmin),maximum_total_time=str(Tmax),
        certified_envelope_upper=str(Vmax),positive_distance_margin=str(gap),
        hypothetical_calibration_only=True,
        must_use_same_transfer_error_class_for_envelope_and_delivery=True,
        empirical_data=False)

    return dict(schema='round1095_finite_role_inheritance_v1',round=1095,status='PASS',
        numpy_version=np.__version__,groups=groups,
        scope=dict(new_split_primitive_is_input=True,
            cognitive_derivation_of_split=False,full_six_protocol_model_certified=False,
            full_mechanical_conservation_of_reception=False,
            whole_conjecture_decided=False,Lorentz_derived=False,
            new_adopted_axioms=0,scientific_count_increment=0,scientific_count_total=3860))


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for key in a: compare(a[key],b[key])
    elif isinstance(a,(list,tuple)):
        assert isinstance(b,(list,tuple)) and len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
    else: assert a==b,(a,b)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    output=calculate()
    path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8') as stream:
            json.dump(output,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
    else: compare(output,json.loads(path.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1095,status='PASS',groups=len(output['groups']),
        saved_result_matches=True,result_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        whole_conjecture_decided=False),ensure_ascii=False))


if __name__=='__main__': main()
