"""593: same joint quantum geometry, actual records, and missing resources.

Original wave-packet diagnostics; no physical state of an invented Einstein
constraint is constructed. Conditional constraint statements are analytic.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import numpy as np
import joint_geometry_work_noise as prior
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_record_constraint_balance_results.json'
spec=importlib.util.spec_from_file_location('entry593',HERE/'round593_drafts/joint_constraint_entry.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def inherited_joint_check():
    r=entry.run();assert r==json.loads(entry.TARGET.read_text('utf8'))
    return dict(frozen_entry_reproduced=True,evidence=r)


def packet(nodes=48):
    h,s,p,kss=prior.packet_data(nodes)
    z,w=np.polynomial.legendre.leggauss(nodes);theta=.3+.4*z
    pg=w*np.exp(-2/(1-z*z));pg/=pg.sum()
    density=pg[:,None,None]*p[None,:,:]*(1+.25*np.sin(theta)[:,None,None]*np.sin(s)[None,:,:])
    density/=density.sum()
    amplitude=np.sqrt(density).reshape(nodes,-1)
    L,dL,_=entry.read.instrument(s)
    volume=prior.W*np.exp(.72*np.sin(theta))
    D=prior.HBAR**2/(2*volume[:,None,None])*kss[None,:,:]*np.sum(dL*dL,axis=-1)[None,:,:]
    return theta,s,density,amplitude,L,D


def conditional_geometry_check():
    theta,s,density,amp,L,D=packet()
    before=amp@amp.T;pooled=np.zeros_like(before);rows=[]
    for r in range(2):
        v=amp*L[...,r].reshape(1,-1)
        block=v@v.T;probability=float(np.trace(block));conditional=block/probability
        distance=.5*float(np.sum(abs(np.linalg.eigvalsh(conditional-before))))
        assert distance>1e-7
        rows.append(dict(record=r,probability=probability,conditional_geometry_trace_distance=distance,
                         conditional_theta_mean=float(theta@np.diag(conditional))))
        pooled+=block
    err=float(np.linalg.norm(pooled-before))
    assert err<1e-14
    return dict(rows=rows,unconditional_geometry_error=err,
                conditional_geometry_changes_do_not_cancel_total_energy_injection=True,
                no_gravity_constraint_state_claim=True)


def repeated_balance_check():
    theta,s,density,amp,L,D=packet()
    Q=float(np.sum(density*D));G=float(np.sum(density*(-.72*np.cos(theta))[:,None,None]*D))
    rows=[]
    for steps in (1,2,4):
        # All original Kraus histories, no waits. D and its geometry
        # derivative are multiplication operators, so the source law
        # composes without replacing a quantum state by its mean field.
        expected_Q=0.;probability=0.
        for history in itertools.product((0,1),repeat=steps):
            weight=np.ones_like(s)
            for r in history:weight*=L[...,r]**2
            block=density*weight[None,:,:]
            probability+=float(np.sum(block))
            expected_Q+=float(np.sum(block*D))
        assert abs(expected_Q-Q)<2e-15 and abs(probability-1)<2e-14
        rows.append(dict(reads=steps,history_probability=probability,
                         next_read_energy_injection=expected_Q,
                         accumulated_energy_injection=steps*Q,
                         accumulated_geometry_source_jump=steps*G,
                         minimum_boundary_apparatus_energy_loss=steps*Q))
    return dict(rows=rows,zero_wait_exact_operator_iteration=True,
                apparatus_loss_assumes_closed_implementation_and_zero_endpoint_interaction=True,
                apparatus_implementation_not_constructed=True)


def run():
    tests=(inherited_joint_check,conditional_geometry_check,repeated_balance_check)
    evidence={f.__name__:f() for f in tests}
    names=('research_note_590.md','research_note_592.md','joint_geometry_work_noise.py',
           'round593_drafts/joint_constraint_entry.py','round593_drafts/joint_constraint_entry_results.json',
           'round593_drafts/joint_constraint_entry.md','round593_drafts/joint_constraint_entry_checks.json')
    return dict(round=593,tests_run=3,failures=0,errors=0,evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original 590 compact quantum-geometry control branch with original reads. Same-state energy/force changes and H2 budget; conditional obstruction for C_g+H_m without apparatus, not a constructed GR countermodel. Energy-preserving replacements and repeated records need declared resources. No autonomous apparatus, Einstein constraints, or continuum completion.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
