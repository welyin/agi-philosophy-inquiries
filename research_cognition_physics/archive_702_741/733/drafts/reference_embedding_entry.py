"""733 entry: embed the original record block in a changed-background reference.

Exact CAR positivity construction; finite original-matrix calibration only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE));sys.path.insert(0,str(ARCHIVE/'round732_drafts'))
import joint_relative_source_development as old
import source_feedback_entry as family
TARGET=HERE/'reference_embedding_entry_results.json'


def run():
    P0,delta0,values,modes,C=old.initial_modes()
    e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2);ce=C@e.conj()
    spanning=np.column_stack((e,ce,P0@e,P0@ce))
    U,s,_=np.linalg.svd(spanning,full_matrices=False);basis=U[:,s>1e-12]
    K=basis@basis.conj().T;Z=np.eye(128)-K
    Q=np.outer(e,e.conj())+np.outer(ce,ce.conj());R=np.eye(128)-2*Q
    record=lambda P:(R@P@R-P)/2
    assert basis.shape[1]==4
    gamma,beta=.12,.04
    Pstar,_,Bnew,Ds=family.flow(gamma,beta,derivatives=False)
    Pnew=K@P0@K+Z@Pstar@Z
    change=Pnew-Pstar;delta_new=record(Pnew)
    errors=dict(original_block_invariance=old.maxabs(K@P0-P0@K),
        block_reality=old.maxabs(C@K.conj()@C-K),
        new_covariance_reality=old.maxabs(C@Pnew.conj()@C+Pnew-np.eye(128)),
        target_record_difference=old.maxabs(delta_new-delta0),
        retained_block=old.maxabs(K@Pnew@K-K@P0@K),
        occupation=abs(np.vdot(e,Pnew@e)-np.vdot(e,P0@e)))
    assert max(errors.values())<3e-13
    ev=np.linalg.eigvalsh(Pnew);assert ev.min()>-2e-13 and ev.max()<1+2e-13
    rank=int(np.count_nonzero(abs(np.linalg.eigvalsh(change))>1e-10));assert rank<=8
    old_energy=old.source(delta0,old.Btime(0));new_energy=old.source(delta_new,Bnew)
    prep_energy=old.source(change,Bnew)
    assert abs(new_energy-old_energy)>1e-3 and abs(prep_energy)>1e-6
    mixedness=float(np.trace(Pnew@(np.eye(128)-Pnew)).real)
    assert mixedness>1e-6
    deps=('research_note_730.md','research_note_731.md','research_note_732.md',
          'joint_relative_source_development.py','joint_relative_source_development_results.json',
          'round732_drafts/source_feedback_entry.py','round732_drafts/source_feedback_entry_results.json')
    return dict(entry_round=733,new_formal_round=False,changed_background=dict(gamma=gamma,beta=beta),
        record_block_dimension=4,reference_change_rank=rank,
        algebra_errors={k:float(v) for k,v in errors.items()},
        min_covariance_eigenvalue=float(ev.min()),max_covariance_eigenvalue=float(ev.max()),
        mixedness=mixedness,original_record_energy=old_energy,changed_background_record_energy=new_energy,
        reference_replacement_energy=prep_energy,changed_record_sources=[old.source(delta_new,D) for D in Ds],
        saved_record_probability=float(np.vdot(e,Pnew@e).real),
        actual_preparation_and_absolute_source_not_closed=True,
        dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original pure reference defines C-stable smooth four-mode record block. Compression of any changed-background Hadamard covariance plus the old block gives a valid generally mixed Hadamard covariance with the exact old record difference. Finite matrix diagnostic checks positivity, reality and full source replacement cost; no free preparation, full quantum Gauss or nonlinear semiclassical existence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
