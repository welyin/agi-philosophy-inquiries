"""594: bounded pointer premeasurement in the original joint model.

Reuses the frozen original-packet entry. The extra matrix check verifies the
terminal copy algebra, not numerical evolution of the full joint Hamiltonian.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_pointer_energy_boundary_results.json'
spec=importlib.util.spec_from_file_location('pointer594entry',HERE/'round594_drafts/premeasurement_energy_entry.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def terminal_copy_check():
    identity=np.eye(2,dtype=complex)
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    ready=np.diag([1.,0.]).astype(complex)
    pointer=(identity+.2*sx+.6*sy+.3*sz)/2
    copy=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],complex)
    after=copy@np.kron(pointer,ready)@copy.conj().T
    reduced=np.trace(after.reshape(2,2,2,2),axis1=1,axis2=3)
    assert np.linalg.norm(reduced-np.diag(np.diag(pointer)))<1e-15
    theta,s,density,amp,L,D=entry.joint.packet(48)
    angles=entry.angle(s)[0]
    mean_A=float(np.sum(density*angles[None,:,:]))
    mean_A2=float(np.sum(density*angles[None,:,:]**2))
    rows=[]
    for kappa in (.4,3.,20.):
        # All scalar-position dependence is integrated with the original
        # correlated packet. Copy does not act on H0, so only V changes.
        v=kappa*(np.pi/3*identity+mean_A*sy)
        before=float(np.trace(np.kron(v,identity)@np.kron(pointer,ready)).real)
        final=float(np.trace(np.kron(v,identity)@after).real)
        expected=-.6*kappa*mean_A
        assert abs(final-before-expected)<2e-14
        rows.append(dict(kappa=kappa,pointer_copy_energy_change=final-before,
                         analytic_copy_energy_change=expected))
    residuals=[]
    for A in angles.ravel()[::13]:
        v=np.pi/3*identity+A*sy
        comm=v@sz-sz@v
        residuals.append(abs(np.linalg.norm(comm,2)-2*A))
        variance=float(np.trace(ready@(v@v)).real-np.trace(ready@v).real**2)
        assert abs(variance-A*A)<1e-15
    return dict(rows=rows,mean_A=mean_A,original_packet_A2=mean_A2,
                record_commutator_identity_error=float(max(residuals)),
                copying_erases_pointer_sigma_y_in_old_marginal=True,
                test_pointer_state_is_algebraic_not_actual_propagated_output=True,
                zero_memory_bare_energy_does_not_cancel_interaction_cost=True)


def run():
    frozen=entry.run()
    assert frozen==json.loads(entry.TARGET.read_text('utf8'))
    evidence=dict(rotation_derivative_identity=frozen['angle_second_derivative_errors'],
                  original_joint_packet_energy=frozen['original_joint_packet_rows'],
                  graph_norm_budget=frozen['bounded_coupling_constants'],
                  terminal_copy_and_stability=terminal_copy_check())
    names=('research_note_593.md','joint_record_constraint_balance.py','joint_geometry_work_noise.py',
           'round594_drafts/premeasurement_energy_entry.md','round594_drafts/premeasurement_energy_entry.py',
           'round594_drafts/premeasurement_energy_entry_results.json','round594_drafts/premeasurement_energy_entry_checks.json')
    return dict(round=594,tests_run=4,failures=0,errors=0,evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original joint H0 with a declared neutral zero-bare-energy pointer and bounded time-independent coupling. Controlled premeasurement, uniform output H0-square budget, and interaction-coherence energy balance. Terminal copy and persistent record are not autonomous implementations. No full-H numerical propagation, GR constraints, continuum, or new matter species derived.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
