"""Executed667 entry: actual614 reference under specified node inclusions.

Finite linear algebra is a check, not a proof of an infinite-volume limit.
The entry note supplies the elementary non-normality proof in empty Fock.
No original Gibbs or continuum chiral state is identified with this reference.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_spinor_subgroup_mass as dictionary
import joint_fermion_gauss_completion as matter


def run():
    w=dictionary.ph_matrix();u=w[:32,:32];v=w[:32,32:]
    covariance=v.conj()@v.T
    occupied=[i for i,x in enumerate(np.diag(covariance)) if x.real==1]
    assert len(occupied)==16 and np.allclose(covariance,np.diag(np.diag(covariance)))
    local_mask=sum(1<<i for i in occupied)
    sizes=(1,2,4,8);traces=[];errors=[];masks=[]
    for m in sizes:
        um=np.kron(np.eye(m),u);vm=np.kron(np.eye(m),v)
        c=vm.conj()@vm.T
        errors.extend([float(np.linalg.norm(um@um.conj().T+vm@vm.conj().T-np.eye(32*m))),
                       float(np.linalg.norm(um@vm.T+vm@um.T)),
                       float(np.linalg.norm(c@c-c)),
                       float(np.linalg.norm(c[:32,:32]-covariance))])
        traces.append(float(np.trace(vm.conj().T@vm).real))
        mask=sum(local_mask<<(32*i) for i in range(m));masks.append(mask)
        assert mask.bit_count()==16*m and mask&((1<<32)-1)==local_mask
    assert traces==[16.*m for m in sizes] and max(errors)==0
    overlaps=[[int(x==y) for y in masks] for x in masks]
    assert overlaps==np.eye(4,dtype=int).tolist()
    colour=matter.gauge.group_exp(np.array([.12,-.14,.08,.04,.09,-.06,.05,.11]),3)
    weak=matter.gauge.group_exp(np.array([.17,-.13,.09]),2)
    r=np.kron(dictionary.left_rep(colour,weak,np.exp(.21j)),np.eye(2))
    gauge_error=float(np.linalg.norm(r@covariance-covariance@r))
    assert gauge_error<1e-12
    deps=('joint_spinor_subgroup_mass.py','joint_fermion_gauss_completion.py',
          'research_note_614.md','research_note_653.md','research_note_666.md')
    return dict(date='2026-10-02',entry_for_round=667,formal_round_complete=False,
                inclusion='Append original32-mode nodes with empty-reference Fock inclusions',
                nodes=list(sizes),offdiagonal_Hilbert_Schmidt_squared=traces,
                original_right_occupation_indices=occupied,
                finite_CAR_and_local_covariance_error=max(errors),
                transformed_vacuum_prefix_overlaps=overlaps,
                common_SM_reference_covariance_error=gauge_error,
                infinite_non_normality_requires_note_proof=True,
                no_Gibbs_or_spatial_refinement_limit_claim=True,
                dependency_hashes={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in deps},
                all_checks_passed=True)


if __name__=='__main__':
    result=run();target=HERE/'reference_scale_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_checks_passed')}))
