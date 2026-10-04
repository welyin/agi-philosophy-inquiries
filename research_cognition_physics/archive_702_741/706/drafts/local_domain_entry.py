"""706 entry: graded local Stinespring and two-sided comparison-energy domain.
Includes explicit705 parity completion; not a new completed scientific round.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARC=HERE.parent
sys.path.insert(0,str(ARC))
import numpy as np
import joint_regional_charge_compression as prior
TARGET=HERE/'local_domain_entry_results.json'


def run():
    blocks=[];offset=0
    for name in ('nu','L','u'):
        dat=prior.boundary.label_data(prior.boundary.MODULES[name]);D=dat['dimension']**2
        mu=1+2*sum(dat['casimirs'])
        blocks.append(dict(name=name,D=D,start=offset,energies=mu+np.array([0.,1.,7.,8.]),
                           parity=np.array([1,1,-1,-1])))
        offset+=4*D
    n=offset
    energy=np.concatenate([np.repeat(b['energies'],b['D']) for b in blocks])
    parity=np.concatenate([np.repeat(b['parity'],b['D']) for b in blocks])
    A=np.diag(energy+1);P=np.diag(parity)
    labels=[]
    for b in blocks:
        for j in range(4):
            sigma=b['parity'][j];g=0 if sigma==1 else 2
            labels.append((b,j,g,float(b['energies'][j]-b['energies'][g])))
    env=np.array([0.]+[v[3] for v in labels])
    plus=np.kron(np.eye(len(env)),A)+np.kron(np.diag(env),np.eye(n))
    Sinfty=np.vstack([np.eye(n)]+[np.zeros((n,n)) for _ in labels])
    rng=np.random.default_rng(706);psi=rng.normal(size=n)+1j*rng.normal(size=n)
    psi/=np.linalg.norm(psi)
    rows=[];bad_parity=[]
    for N in (10.,24.,36.,45.):
        keep=energy<=N;ops=[np.diag(keep.astype(float))]
        for b,j,g,benergy in labels:
            D,s=b['D'],b['start'];k=np.zeros((n,n))
            if b['energies'][j]>N:
                k[s+g*D:s+(g+1)*D,s+j*D:s+(j+1)*D]=np.eye(D)
                if b['parity'][j]<0:
                    wrong=np.zeros((n,n))
                    wrong[s:s+D,s+j*D:s+(j+1)*D]=np.eye(D)
                    bad_parity.append(float(np.linalg.norm(P@wrong-wrong@P,2)))
            ops.append(k)
        S=np.vstack(ops)
        iso=float(np.linalg.norm(S.conj().T@S-np.eye(n)))
        forward=float(np.linalg.norm(plus@S-S@A))
        backward=float(np.linalg.norm(A@S.conj().T-S.conj().T@plus))
        parity_error=max(float(np.linalg.norm(P@k-k@P)) for k in ops)
        observed=float(np.linalg.norm(plus@(S-Sinfty)@psi)**2)
        expected=float(2*np.linalg.norm(A@((~keep)*psi))**2)
        assert max(iso,forward,backward,parity_error)<1e-12
        assert abs(observed-expected)<3e-12
        rows.append(dict(cutoff=N,isometry_error=iso,forward_energy_intertwiner_error=forward,
            adjoint_energy_intertwiner_error=backward,even_Kraus_error=parity_error,
            weighted_difference_squared=observed,twice_original_graph_tail=expected))
    assert max(bad_parity)==2.
    deps=('research_note_623.md','research_note_625.md','research_note_637.md',
        'research_note_705.md','joint_regional_charge_compression.py')
    return dict(date='2026-10-02',entry_round=706,latest_formal_round=705,
        entry_checks=2,all_checks_passed=True,rows=rows,
        using_one_reset_vector_across_parities_has_commutator_norm=max(bad_parity),
        local_dimension=n,common_auxiliary_dimension=len(env),
        dependencies={p:hashlib.sha256((ARC/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(graded705_construction_completed=True,
            fixed_auxiliary_energy_labels_not_physical_free_reset=True,
            two_sided_domain_intertwiner=True,
            infinite_strong_graph_limit_proved_analytically=True,
            no_new_effective_Hamiltonian_or_own_Gibbs_family=True,
            actual_second_source_full_history_theorem_next=True,
            not_completed_new_round=True,active_goal_unchanged=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(entry706=True,rows=result['rows'],bad_parity=result[
        'using_one_reset_vector_across_parities_has_commutator_norm'],passed=True)))

