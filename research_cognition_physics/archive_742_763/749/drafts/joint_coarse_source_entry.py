"""749 entry: audit a specific record-erasure map on731 original constraints.

The two branches are inherited prescribed smooth source tests, NOT states
produced by an actual quantum instrument. No quantum conditional spacetime,
GR emergence, or unique geometric averaging prescription is claimed.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_source_constraint_response as old
TARGET=HERE/'joint_coarse_source_entry_results.json'
def run():
    data=old.setup();q=data['q'];psi0=data['psi'];tensor0=data['tensor'];rows=[]
    for eps in (.04,.02,.01):
        branches=[]
        for t in (eps,-eps):
            z,p,E,E0=old.updated_data(t)
            psi,tensor,info=old.geo.solve_hamiltonian(z,initial=psi0)
            err=old.independent_constraint(z,psi,tensor,t)
            branches.append((psi,tensor,p,E,E0,err))
        meanpsi=(branches[0][0]+branches[1][0])/2
        meantensor=(branches[0][1]+branches[1][1])/2
        p=(branches[0][2]+branches[1][2])/2
        E=(branches[0][3]+branches[1][3])/2
        E0=(branches[0][4]+branches[1][4])/2
        # Averaged prescribed sigma,J,rho are zero, as are averaged color fields.
        average_fields_error=max(old.maximum(p-q['p']),old.maximum(E-q['f']['E']),
                                 old.maximum(E0-q['f']['E0']),old.maximum(meantensor-tensor0))
        gaussian=old.maximum(old.gauss(q,data['k'],p,E,E0))
        momentum=old.maximum(sum(old.geo.derivative(meantensor[..., :,i],i) for i in range(3))+q['mom'])
        defect=old.independent_constraint(q,meanpsi,meantensor,0.)
        rows.append(dict(epsilon=eps,prescribed_branches_not_quantum_record_branches=True,
            largest_branch_H_constraint_error=max(b[-1] for b in branches),
            mean_canonical_matter_and_tensor_error=average_fields_error,
            mean_Gauss_error=gaussian,mean_momentum_error=momentum,
            mean_psi_minus_original=old.maximum(meanpsi-psi0),
            erased_label_mean_H_constraint_defect=defect,
            defect_divided_by_epsilon_squared=defect/(eps*eps)))
        assert average_fields_error<1e-12 and gaussian<1e-12 and momentum<1e-12
        assert max(b[-1] for b in branches)<3e-8
        assert defect>1000*max(b[-1] for b in branches)
    deps=('research_note_634.md','research_note_706.md','research_note_708.md',
          'research_note_724.md','research_note_725.md','research_note_731.md',
          'research_note_741.md','research_note_748.md','joint_source_constraint_response.py')
    return dict(entry_round=749,latest_completed_round=748,formal_test_count_unchanged=3434,
        inherited_model='731 original non-flat matter and all initial constraints; prescribed source family.',
        additional_candidate_map='Erase equally weighted +/- source label and average psi, all canonical matter momenta, and conformal shear. This is one explicit map, not all geometric coarse grainings.',
        rows=rows,no_rigorous_continuum_error_certificate=True,
        quantum_record_realization_not_claimed=True,arbitrary_cognitive_principle_refutation_not_claimed=True,
        dependency_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},
        next='Derive the exact second-order correction to this map from731 fields, retain all source/field covariances, and separately test whether the same marked source family comes from an actual record.')
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
