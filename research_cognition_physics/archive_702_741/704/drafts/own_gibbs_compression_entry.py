"""704 entry: own finite Gibbs initial state in existing625 real histories.
A corollary/application, not a new completed round. Original625 finite64
diagnostic retained, without claiming a full graph thermal calculation.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_source_preserving_compression as old
import joint_geometry_thermal_limit as prior
TARGET=HERE/'own_gibbs_compression_entry_results.json'


def run():
    assert prior.run()==json.loads(prior.TARGET.read_text('utf8'))
    receipt=json.loads((ARCHIVE/'research_round_625_checks.json').read_text('utf8'))
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt.get(key,{}).items():
            assert hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest()==digest
    d=old.diagnostic();weights=(d['E']-d['E'][0]+1)**2
    original_moment=float(weights@d['prob']);rows=[]
    for rank in (6,12,24,40):
        p=float(d['prob'][rank:].sum());own,_,_=old.prepare(d,rank,True)
        compressed,_,_=old.prepare(d,rank,False)
        difference=own-compressed
        identity=p*own.copy();identity[0]-=p
        error=float(np.max(abs(difference-identity)));assert error<3e-16
        weighted=float(weights[:rank]@abs(difference))
        bound=p*(original_moment/(1-p)+weights[0])
        assert weighted<=bound+1e-13
        exact_norm=float(abs(difference).sum())
        shape=old.histories(d,rank,.012,rethermalize=True)
        previous=old.histories(d,rank,.012,rethermalize=False)
        history_gap=float(np.max(abs(shape-previous)))
        assert history_gap<=exact_norm+1e-13
        rows.append(dict(rank=rank,original_thermal_tail=p,exact_identity_error=error,
            initial_trace_norm_difference=exact_norm,weighted_A2_difference=weighted,
            analytic_weighted_bound=bound,finite_real_two_history_gap=history_gap))
    # The weighted difference need not decrease at every finite cutoff.
    # The certified upper bound decreases with the nested thermal tail.
    assert all(rows[i+1]['analytic_weighted_bound']<rows[i]['analytic_weighted_bound']
               for i in range(len(rows)-1))
    deps=('research_note_603.md','research_note_623.md','research_note_624.md','research_note_625.md',
          'research_note_702.md','research_note_703.md','joint_source_preserving_compression.py',
          'joint_source_preserving_compression_results.json','research_round_703_checks.json')
    return dict(date='2026-10-02',entry_round=704,latest_formal_round=703,new_formal_round=False,
        rows=rows,old_diagnostic_dimension=d['dimension'],original_shifted_second_moment=original_moment,
        scope=dict(full_model_corollary_reuses_625_uniform_coefficient_bounds=True,
            original_initial_state_must_be_Gibbs_not_arbitrary_state=True,
            finite_cutoff_states_not_claimed_identical=True,
            own_Gibbs_for_initial_baseline_not_rethermalized_after_each_source=True,
            alternative625_family_not703_log_transfer_derivative_proof=True,
            preparation_geometry_source_jets_not_yet_combined=True,
            numerical_calibration_is_old64_radial_no_full_CAR=True,
            no_free_physical_tail_reset_or_GR_claim=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry=704,formal=703,rows=result['rows'],all_checks_passed=True)))
