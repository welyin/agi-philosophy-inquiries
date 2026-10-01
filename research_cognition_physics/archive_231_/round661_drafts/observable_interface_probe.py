"""661 actual entry: local Weyl components and support, not a completed round."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_local_mirror_process as old


def run():
    q=old.data(nx=2,nt=4);r=q['r'];n=2*r
    jm=np.kron(np.kron(np.eye(len(q['sites'])),old.old.old.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(len(q['sites'])),old.old.old.VP),np.eye(16))
    pv=q['v']@q['v'].conj().T;di=np.linalg.inv(q['D'])
    original=jm.conj().T@q['v']@np.linalg.inv(q['Kl'])
    dirac=jm.conj().T@di@jp
    mature_error=old.err(original-dirac)
    assert mature_error<3e-12
    # Select physically local spin components from the full candidate.
    zero=np.zeros((r,n),complex)
    select=np.block([[jm.conj().T,zero],[zero,jp.T]])
    bare=select@(-np.linalg.inv(q['N']))@select.T
    target=np.block([[np.zeros((r,r)), -original],
                     [original.T,np.zeros((r,r))]])
    difference=old.norm(bare-target)
    # The exact660 map realizes those target Weyl variables by dressed fields.
    to_eta=q['T']@q['S'].conj().T
    choose=np.zeros((2*r,4*r),complex)
    choose[:r,r:2*r]=jm.conj().T@q['v']
    choose[r:,3*r:]=np.eye(r)
    dressed=choose@to_eta
    mapped=dressed@(-np.linalg.inv(q['N']))@dressed.T
    dressed_error=old.err(mapped-target)
    assert dressed_error<3e-12 and difference>.1
    # Local time coordinates of selected original Weyl fields vs candidate.
    half=np.repeat([t>=2 for t,x in q['sites']],32)
    half=np.concatenate((half,half))
    candidate_half=np.repeat([t>=2 for t,x in q['sites']],64)
    candidate_half=np.concatenate((candidate_half,candidate_half))
    across=old.norm(dressed[np.ix_(half,~candidate_half)])
    assert across>.1
    return dict(date='2026-10-02',status='661 entry only; no completed new round',
        nx=2,nt=4,mature_free_Weyl_component_identity_error=mature_error,
        conditional_bare_local_component_Wick_difference=difference,
        transformed_Weyl_component_Wick_error=dressed_error,
        exact_map_across_time_half_norm=across,
        original_same_D_and_all_internal_matrices=True,
        source_free_S9_average_difference_not_claimed=True,
        excludes_only_bare_identification_and_support_preservation_of_this_map=True,
        next='Assess product of original free Weyl local algebra with original auxiliary RP algebra, preserving common sources; do not infer failure of all positive realizations.',
        dependencies={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ('joint_local_mirror_process.py','research_note_660.md','research_note_657.md')})


if __name__=='__main__':
    result=run()
    with (HERE/'observable_interface_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))
