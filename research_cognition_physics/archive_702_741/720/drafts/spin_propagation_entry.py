"""720 entry: exact spin protection must use the same propagation dictionary."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_smooth_mode_contract as shared
import joint_continuum_record_sources as physical
TARGET=HERE/'spin_propagation_entry_results.json'
PAULI=physical.PAULI


def rotation(axis,angle):
    axis=np.asarray(axis,dtype=float);axis/=np.linalg.norm(axis)
    u=np.cos(angle/2)*np.eye(2)-1j*np.sin(angle/2)*np.einsum('i,ijk->jk',axis,PAULI)
    return u


def lift(u):
    w=np.kron(np.eye(2),u);z=np.zeros_like(w)
    return np.block([[w,z],[z,w.conj()]])


def run():
    rng=np.random.default_rng(72004);coefficient_errors=[]
    for _ in range(6):
        points=rng.normal(size=(2,5))*.25
        m,d,h=shared.coefficients(points,[shared.group.sample(rng) for _ in range(2)],
                                   np.array([.21,.18]))
        u=rotation(rng.normal(size=3),.73)
        w=np.kron(np.eye(32),u)
        coefficient_errors.append(float(max(np.linalg.norm(w@(m+h)@w.conj().T-m-h),
                                            np.linalg.norm(w@d@w.T-d))))
    assert max(coefficient_errors)<1e-12
    rows=[]
    for k in (np.array([.3,0,0]),np.array([.23,.41,-.29]),np.array([0,0,1.1])):
        axis=np.cross(k,np.eye(3)[int(np.argmin(abs(k)))])
        u=rotation(axis,np.pi);v=lift(u)
        B=physical.bdg(k);B0=physical.B0
        flip=float(np.linalg.norm(v@B@v.conj().T-physical.bdg(-k)))
        separation=float(np.linalg.norm(B-v@B@v.conj().T,ord=2)/2)
        attainable=float(np.linalg.norm(B-B0,ord=2))
        assert max(flip,abs(separation-np.linalg.norm(k)),
                   abs(attainable-np.linalg.norm(k)))<1e-12
        u2=rotation([.21,.33,-.61],.83);v2=lift(u2)
        rot=np.array([[np.trace(PAULI[i]@u2@PAULI[j]@u2.conj().T).real/2
                       for j in range(3)] for i in range(3)])
        covariance=float(np.linalg.norm(v2@B@v2.conj().T-physical.bdg(rot@k)))
        assert covariance<1e-12
        rows.append(dict(momentum=k.tolist(),norm=float(np.linalg.norm(k)),
                         lower_bound_for_fixed_spin_commutant=separation,
                         achieved_by_original_mass_matrix=attainable,
                         half_turn_identity_error=flip,
                         simultaneous_space_spin_covariance_error=covariance))
    names=('research_note_604.md','research_note_614.md','research_note_633.md',
           'research_note_667.md','research_note_719.md','joint_continuum_record_sources.py',
           'joint_smooth_mode_contract.py')
    return dict(entry_round=720,new_formal_round=False,
                original_64_mode_spin_scalar_branch_errors=coefficient_errors,
                original_physical_neutral_BdG_rows=rows,
                hypotheses=dict(same_fixed_momentum_and_spin_observable_dictionary=True,
                                tested_scalar_spin_hopping_branch_only=True,
                                allowed_604_Weyl_hopping_not_excluded=True),
                no_general_unification_or_space_dimension_no_go=True,
                no_doubling_recalculation=True,
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
