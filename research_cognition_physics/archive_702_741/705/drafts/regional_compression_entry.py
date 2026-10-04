"""705 entry: regional algebra under the existing global spectral compression.
Corollary and calibration only; not a completed numbered scientific round.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ARC=HERE.parent
sys.path.insert(0,str(ARC))
import numpy as np
import joint_preparation_history_limit as previous
TARGET=HERE/'regional_compression_entry_results.json'

def hs(x):return float(np.linalg.norm(x,'fro'))

def run():
    d=previous.data()
    n=24
    while n<len(d['E']) and abs(d['E'][n]-d['E'][n-1])<1e-9:n+=1
    a=d['A'][:n];ell=d['L'][0][:n,:n]
    assert np.linalg.norm(ell-ell.conj().T)<1e-12
    X=np.kron(ell,np.eye(n));Y=np.kron(np.eye(n),ell)
    energies=(a[:,None]+a[None,:]).ravel()
    beta=.12;raw=np.exp(-beta*(energies-energies.min()));rho=raw/raw.sum()
    order=np.argsort(energies)
    normX=float(np.linalg.norm(X,2));normY=float(np.linalg.norm(Y,2))
    kX=float(np.linalg.norm(energies[:,None]*X/energies[None,:],2))
    kY=float(np.linalg.norm(energies[:,None]*Y/energies[None,:],2))
    m2=float(np.dot(energies**2,rho))
    full_comm=hs(X@Y-Y@X);assert full_comm<1e-12
    rows=[]
    for desired in (80,160,320,len(energies)):
        rank=desired
        while rank<len(energies) and abs(energies[order[rank]]-energies[order[rank-1]])<1e-9:rank+=1
        ix=order[:rank];outside=order[rank:]
        x=X[np.ix_(ix,ix)];y=Y[np.ix_(ix,ix)]
        qxp=X[np.ix_(outside,ix)];qyp=Y[np.ix_(outside,ix)]
        comm=x@y-y@x
        defect=Y[np.ix_(ix,outside)]@qxp-X[np.ix_(ix,outside)]@qyp
        residual=hs(comm-defect);assert residual<2e-12
        r=rho[ix]/rho[ix].sum();root=np.sqrt(r)
        weighted=hs(comm*root[None,:])
        direct_bound=normY*hs(qxp*root[None,:])+normX*hs(qyp*root[None,:])
        assert weighted<=direct_bound+2e-12
        cutoff=float(energies[ix].max())
        p_tail=float(rho[outside].sum())
        energy_bound=(normY*kX+normX*kY)*np.sqrt(m2/(1-p_tail))/cutoff
        assert weighted<=energy_bound+2e-12
        rows.append(dict(requested_rank=desired,rank=rank,cutoff=cutoff,
            commutator_operator_norm=float(np.linalg.norm(comm,2)),
            thermal_weighted_commutator_norm=weighted,exact_tail_upper_bound=direct_bound,
            inherited_energy_upper_bound=energy_bound,thermal_discarded_probability=p_tail,
            algebra_identity_residual=residual))
    assert max(r['commutator_operator_norm'] for r in rows[:-1])>1e-5
    assert rows[-1]['thermal_weighted_commutator_norm']<1e-12
    deps=('research_note_617.md','research_note_625.md','research_note_637.md',
          'research_note_644.md','research_note_704.md','joint_preparation_history_limit.py')
    return dict(date='2026-10-02',formal_round_still=704,next_round=705,
        entry_checks=2,all_checks_passed=True,
        calibration=dict(local_dimension=n,total_dimension=len(energies),beta=beta,
            two_independent_copies_of_declared702_neutral_diagnostic=True,
            changed_temperature_declared=True,not_full_graph_or_connected_interaction_simulation=True),
        full_regional_commutator_residual=full_comm,rows=rows,
        dependency_hashes={f:hashlib.sha256((ARC/f).read_bytes()).hexdigest() for f in deps},
        scope=dict(exact_general_compression_identity=True,
            thermal_tail_estimate_corollary_of_existing625_graph_bounds=True,
            finite_cutoff_need_not_preserve_regional_commutation=True,
            no_no_go_for_local_physics_or_continuum=True,
            not_completed_new_scientific_round=True,active_goal_unchanged=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(formal_round=704,entry705=True,rows=result['rows'],passed=True)))

