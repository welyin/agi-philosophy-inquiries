"""701 entry: finite-box zero-set witness and diagonal-limit contract audit.
This does not evaluate a continuum interacting integral or prove its existence.
"""
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_anisotropic_time_contract as previous
TARGET=HERE/'joint_scale_entry_results.json'


def run():
    assert previous.run()==json.loads(previous.TARGET.read_text('utf8'))
    rows=[]
    for ns,nt in ((2,2),(3,4),(4,6),(3,8)):
        sites=list(itertools.product(range(nt),range(ns)));n=len(sites)
        shifts=[]
        for temporal in (False,True):
            tmat=np.zeros((n,n),complex)
            for i,(t,x) in enumerate(sites):
                dest=((t+1)%nt,x) if temporal else (t,(x+1)%ns)
                tmat[i,sites.index(dest)]=-1 if temporal and t+1==nt else 1
            shifts.append(tmat)
        x=-np.eye(4*n,dtype=complex)
        for shift,gamma in zip(shifts,(previous.GAMMA[0],previous.GAMMA[3])):
            x+=np.kron(np.eye(n)-(shift+shift.conj().T)/2,np.eye(4))
            x+=np.kron((shift-shift.conj().T)/2,gamma)
        h=np.kron(np.eye(n),previous.G5)@x
        eig=np.linalg.eigvalsh(h);gap=float(min(abs(eig)));assert gap>=1-2e-14
        expected=[]
        for j,l in itertools.product(range(nt),range(ns)):
            kt=(2*j+1)*np.pi/nt;ks=2*np.pi*l/ns
            r=np.sqrt(1+2*(1-np.cos(ks))*(1-np.cos(kt)))
            expected.extend([-r,-r,r,r])
        error=float(np.max(abs(eig-np.sort(expected))));assert error<3e-14
        rows.append(dict(ns=ns,nt=nt,spin_matrix_dimension=4*n,
            full_16_channel_dimension=64*n,internal_tensor_identity_exact=True,gap=gap,spectrum_error=error))
    # Normalization schedule: exact algebra, not a bound on the actual unknown Z_n.
    schedules=[]
    for n in (1,4,12):
        z=F(1,10**n);m=F(10**(2*n));eps=F(1,2**n)
        delta=z*eps/(2*(1+m))
        assert delta<=z/2
        ratio_bound=delta*(1+m)/(z-delta)
        assert ratio_bound<=eps
        schedules.append(dict(n=n,control_Z=str(z),control_M=str(m),epsilon=str(eps),
            sufficient_delta=str(delta),ratio_bound=str(ratio_bound)))
    deps=('research_note_646.md','research_note_659.md','research_note_673.md','research_note_677.md',
        'research_note_686.md','research_note_699.md','research_note_700.md',
        'joint_anisotropic_time_contract.py','research_round_700_checks.json',
        'research_note_382.md','research_note_383.md','research_note_384.md','research_note_386.md',
        'research_note_425.md','research_note_522.md','research_note_523.md')
    return dict(date='2026-10-02',entry_round=701,latest_formal_round=700,new_formal_round=False,
        original_free_AP_zero_set_witnesses=rows,normalization_algebra_controls=schedules,
        scope=dict(new_many_time_sewing_extension_explicit=True,
            global_Haar_zero_set_argument_not_conditional_fixed_spatial_history=True,
            auxiliary_diagonal_approximation_does_not_prove_physical_continuum=True,
            no_uniform_gap_or_uniform_Z_lower_added_as_existence_requirement=True,
            actual_nonzero_Z_unknown=True,effective_resource_bound_unknown=True,
            inherited_negative_fixed_candidate_not_changed=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry=701,formal=700,original_AP_witnesses=len(result['original_free_AP_zero_set_witnesses']),all_checks_passed=True)))
