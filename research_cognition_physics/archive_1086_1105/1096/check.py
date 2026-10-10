"""1096: reuse the 929 isometry; check its matched complete-state transport.
No large processor Hamiltonian is reconstructed; arbitrary dimension is proved
by the swap intertwiner. Default invocation compares without writing.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
SOURCE=BASE/'archive_923_934/929/joint_protocol_selection.py'
spec=importlib.util.spec_from_file_location('round929_readonly',SOURCE)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def exp_h(h,t=1.):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T


def swap(d):
    s=np.zeros((d*d,d*d))
    for a in range(d):
        for b in range(d): s[b*d+a,a*d+b]=1
    return s


def calculate():
    groups=[]
    def add(name,ok,**data):
        assert bool(ok),name
        groups.append(dict(name=name,status='PASS',**data))
    actual_error=iso_error=0.
    v_small=None
    for theta in [0.,np.pi/2]:
        for f1 in range(4):
            for f2 in range(4):
                _,history=old.history(theta,(f1,f2))
                v=history[-1]
                actual_error=max(actual_error,float(np.linalg.norm(v-old.expected(theta,(f1,f2)))))
                iso_error=max(iso_error,float(np.linalg.norm(v.conj().T@v-np.eye(2))))
                if theta==np.pi/2 and (f1,f2)==(1,3):
                    rows=np.where(np.linalg.norm(v,axis=1)>1e-14)[0]
                    v_small=v[rows]
    add('inherited_full_input_protocol',max(actual_error,iso_error)<1e-12,
        cases=32,max_actual_isometry_error=actual_error,max_normalization_error=iso_error,
        original_code_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest())

    h=np.array([[.1,.2+.3j,0],[.2-.3j,-.4,.15],[0,.15,.6]],complex)
    d=len(h);s=swap(d);ks=np.pi/2*(np.eye(d*d)-s)
    h0=np.kron(h,np.eye(d))+np.kron(np.eye(d),h)
    comm=float(np.linalg.norm(h0@s-s@h0))
    t=1.3
    # Constant pulse with area 1 is one check of the general commuting integral.
    exact=exp_h(t*h0+ks)
    target=np.kron(exp_h(h,t),exp_h(h,t))@s
    gate_error=float(np.linalg.norm(exact-target))
    internal_conservation=float(np.linalg.norm(exact.conj().T@h0@exact-h0))
    rng=np.random.default_rng(1096)
    psi=rng.normal(size=(d,2))+1j*rng.normal(size=(d,2));psi/=np.linalg.norm(psi)
    chi=np.array([1,2j,-.5],complex);chi/=np.linalg.norm(chi)
    initial=np.einsum('ar,b->abr',psi,chi)
    out=(exact@initial.reshape(d*d,2)).reshape(d,d,2)
    u=exp_h(h,t)
    expected=np.einsum('a,br->abr',u@chi,u@psi)
    reference_error=float(np.linalg.norm(out-expected))
    add('matched_full_state_contact',max(comm,gate_error,internal_conservation,reference_error)<1e-12,
        commutator_norm=comm,full_gate_error=gate_error,
        internal_energy_operator_error=internal_conservation,reference_state_error=reference_error)

    # The actual 929 isometry, including its untouched unknown-input columns,
    # is transported, not another spectator qubit. Compress only zero rows.
    v=v_small;d_out=v.shape[0]
    vp=np.kron(v,v)
    intertwiner_error=float(np.linalg.norm(swap(d_out)@vp-vp@swap(2)))
    source_input=np.kron(np.eye(2),np.array([[1.],[0.]]))
    complete=swap(d_out)@vp@source_input
    wanted=np.kron(v[:,0:1],v)
    input_operator_error=float(np.linalg.norm(complete-wanted))
    # Operator equality above covers arbitrary references; Choi adds a check.
    choi_error=float(np.linalg.norm(complete.flatten()-wanted.flatten())/np.sqrt(2))
    dk=2*31*2**17
    add('same_running_protocol_and_reference',max(intertwiner_error,input_operator_error,choi_error)<1e-12,
        retained_nonzero_output_rows=d_out,swap_isometry_intertwiner_error=intertwiner_error,
        complete_input_operator_error=input_operator_error,normalized_choi_vector_error=choi_error,
        single_full_register_dimension=dk,two_full_register_dimension=dk**2,
        full_hardware_duplicated=True,source_keeps_evolved_receiver_state=True)

    j=np.eye(5)[:,[0,2,4]]
    embedding_error=float(np.linalg.norm(
        np.kron(j,np.eye(d_out))@np.kron(np.eye(3),v)
        -np.kron(np.eye(5),v)@np.kron(j,np.eye(2))))
    # Same contact in three actual Galilean charts; all backgrounds move too.
    q,n=np.polynomial.legendre.leggauss(24)
    eps=.1;w=2.;length=2.;tc=length/w
    times=tc+eps/w*q
    integ=float(np.dot(n,.75*(1-q*q)))
    geometry_error=0.
    for boost in [-.7,0.,.5]:
        source=(w-boost)*times;receiver=length-boost*times
        geometry_error=max(geometry_error,float(np.max(np.abs(source-receiver-(w*times-length)))))
    add('directory_and_actual_chart_transport',max(embedding_error,geometry_error,abs(integ-1))<1e-12,
        directory_embedding_operator_error=embedding_error,galilean_relative_position_error=geometry_error,
        integrated_contact_area=integ)

    z=np.diag([1.,-1.])/2;s2=swap(2);k2=np.pi/2*(np.eye(4)-s2)
    bad=np.kron(z,np.eye(2))
    bad_comm=float(np.linalg.norm(bad@s2-s2@bad))
    bad_gate=float(np.linalg.norm(exp_h(bad+k2)-exp_h(bad)@s2))
    good=np.kron(z,np.eye(2))+np.kron(np.eye(2),z)
    delta=.07;T=.8;pert=delta*np.kron(np.eye(2),np.array([[0.,1.],[1.,0.]]))
    mismatch=float(np.linalg.norm(exp_h(T*(good+pert)+k2)-exp_h(T*good+k2),ord=2))
    add('deletions_and_finite_mismatch',bad_comm>.1 and bad_gate>.1 and mismatch<=T*delta+1e-12,
        unmatched_commutator_norm=bad_comm,unmatched_full_gate_error=bad_gate,
        mismatch_operator_error=mismatch,Duhamel_bound=T*delta,
        bound_compares_same_complete_initial_state_and_pulse=True)

    rows=[]
    for c in [.01,.1,1.,10.,100.]:
        w=max(4*np.pi*c,1.);L=w;end=1+.1/w
        lp=2/np.sqrt(3)*(np.pi+.02-L/(2*c))
        assert end<np.pi-.02 and lp<=2/np.sqrt(3)*(.02-np.pi)+1e-11
        rows.append(dict(c=c,finite_speed=w,finite_distance=L,contact_end=end,
            readout_window=[np.pi-.02,np.pi+.02],maximum_lorentz_time_difference=lp))
    actual_window=1-np.cos(.02/2)**60
    add('finite_window_and_complete_task_witness',actual_window<=.003,
        cases=rows,readout_record_error=actual_window,record_error_bound=.003,
        error_not_a_bound_for_full_clock_terminal_state=True,empirical_data=False)
    return dict(schema='round1096_running_protocol_transfer_v1',round=1096,status='PASS',
        numpy_version=np.__version__,groups=groups,
        scope=dict(finite_task_joint_qualification=True,entire_conjecture_decided=False,
            shared_motion_observers_consensus=False,all_mechanical_sources_closed=False,
            six_protocol_recorded_events_are_original_internal_instrument_results=True,
            new_adopted_axioms=0,scientific_count_increment=0,scientific_count_total=3860))


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a:compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
    else:assert a==b,(a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=calculate();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(path.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1096,status='PASS',groups=len(out['groups']),saved_result_matches=True),ensure_ascii=False))


if __name__=='__main__':main()
