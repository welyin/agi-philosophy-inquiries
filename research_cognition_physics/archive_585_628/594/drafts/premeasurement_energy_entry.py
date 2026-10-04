"""594 entry: original instrument's pointer dilation and its energy boundary.

Checks original curved packet kinetic cost and rotation derivatives. It does
not propagate the full joint Hamiltonian or claim a completed stable detector.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_record_constraint_balance as joint
import joint_geometry_work_noise as prior
TARGET=HERE/'premeasurement_energy_entry_results.json'


def angle(s):
    den=np.sqrt(1-np.sin(s)**2/4)
    return .5*np.arccos(np.sin(s)/2),-np.cos(s)/(4*den),3*np.sin(s)/(16*den**3)


def run():
    s=np.linspace(-np.sqrt(12),np.sqrt(12),701)
    A,dA,ddA=angle(s);L,dL,ddL=joint.entry.read.instrument(s)
    assert np.max(abs(np.cos(A)-L[:,0]))<5e-16
    assert np.max(abs(np.sin(A)-L[:,1]))<5e-16
    assert np.max(abs(dA))<=.25*(1+1e-14)
    assert np.max(abs(ddA))<=1/(2*np.sqrt(3))*(1+1e-14)
    errors=[]
    for step in (.02,.01,.005):
        plus=angle(s+step)[0];minus=angle(s-step)[0]
        errors.append(float(np.max(abs((plus-2*A+minus)/step**2-ddA))))
    assert errors[-1]<errors[0]/12
    packet_rows=[]
    for nodes in (32,48,64):
        theta,s,density,amp,L,D=joint.packet(nodes)
        h,s,p,kss=prior.packet_data(nodes)
        x=(h-prior.packet.HCENTER)/prior.packet.HRADIUS
        y=(s-prior.packet.SCENTER)/prior.packet.SRADIUS
        F=joint.entry.read.M-(h*h+s*s)/6
        coords=np.stack((h,s),axis=-1)
        g=F[...,None,None]*(np.eye(2)-coords[..., :,None]*coords[...,None,:]/(6*joint.entry.read.M))
        grad=np.stack((-2*x/(prior.packet.HRADIUS*(1-x*x)**2),-2*y/(prior.packet.SRADIUS*(1-y*y)**2)),axis=-1)
        grad=grad-coords/(2*F[...,None])
        grad=np.broadcast_to(grad,(nodes,*grad.shape)).copy()
        corr=1+.25*np.sin(theta)[:,None,None]*np.sin(s)[None,:,:]
        grad[...,1]+=.125*np.sin(theta)[:,None,None]*np.cos(s)[None,:,:]/corr
        volume=prior.W*np.exp(.72*np.sin(theta))
        def kinetic(vec):
            dens=np.einsum('thsi,hsij,thsj->ths',vec.conj(),g,vec).real
            return prior.HBAR**2/2*float(np.sum(density*dens/volume[:,None,None]))
        before=kinetic(grad)
        A,dA,ddA=angle(s)
        rotated0=np.cos(A)[None,...,None]*grad
        rotated1=np.sin(A)[None,...,None]*grad
        rotated0[...,1]+=-np.sin(A)[None,:,:]*dA[None,:,:]
        rotated1[...,1]+=np.cos(A)[None,:,:]*dA[None,:,:]
        after=kinetic(rotated0)+kinetic(rotated1)
        predicted=float(np.sum(density*D))
        assert abs(after-before-predicted)<5e-11
        mean_A2=float(np.sum(density*A[None,:,:]**2))
        assert mean_A2>(np.pi/6)**2
        packet_rows.append(dict(nodes=nodes,original_node_kinetic_before=before,
                                ideal_pointer_rotated_node_kinetic=after,
                                ideal_source_energy_increment=after-before,
                                original_instrument_energy_increment=predicted,
                                initial_total_variance_kappa_squared_coefficient=mean_A2))
    wmin=prior.W*np.exp(-.72);hbar=prior.HBAR;M=joint.entry.read.M
    a1=.25;a2=1/(2*np.sqrt(3));b0=np.sqrt(6*M)/3
    C1=hbar*a1*np.sqrt(2*M/wmin)
    C0=hbar*hbar/(2*wmin)*(M*(a2+a1*a1)+b0*a1)
    # Explicit assumed budgets only; not estimates of the complete H spectrum.
    E,B=1000.,2e6;duhamel=np.sqrt(B)+C1*np.sqrt(E)+C0
    mu=1.;a=C1/(2*np.sqrt(mu))
    comm0=hbar*hbar/(2*wmin)*(M*a2+b0*a1)
    b=C1*np.sqrt(mu)/2+comm0
    graph_bound=(np.sqrt(B)+b/a)*np.exp(a)-b/a
    return dict(status='594 premeasurement entry; not completed numbered round',checks_passed=True,
                angle_second_derivative_errors=errors,original_joint_packet_rows=packet_rows,
                bounded_coupling_constants=dict(angle_derivative_bound=a1,angle_second_derivative_bound=a2,
                    commutator_sqrt_energy_coefficient=C1,rotation_commutator_constant=C0,
                    assumed_input_E=E,assumed_input_H2=B,duhamel_error_prefactor=duhamel,
                    actual_H0_graph_norm_uniform_bound=graph_bound,
                    input_budget_nonemptiness_not_numerically_certified=True),
                full_joint_unitary_not_numerically_propagated=True,
                terminal_pointer_dephasing_energy_is_analytic_identity=True,
                extra_pointer_and_coupling_are_declared_inputs=True,
                stable_record_and_internal_switch_not_implemented=True,
                dependency_hashes={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in
                    ('research_note_593.md','joint_record_constraint_balance.py','joint_geometry_work_noise.py')})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
