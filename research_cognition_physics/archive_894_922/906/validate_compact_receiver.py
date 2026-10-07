"""906 fast independent reproduction; stored trajectories are not silently rerun."""
from pathlib import Path
import argparse,hashlib,json
import compact_receiver_checks as experiment
import receiver_density_pair_identity as pair
HERE=Path(__file__).resolve().parent;TARGET=HERE/'compact_receiver_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(rerun=False):
    raw=json.loads(experiment.TARGET.read_text('utf-8'));paired=json.loads(pair.TARGET.read_text('utf-8'))
    if rerun:assert experiment.run()==raw
    assert experiment.jets()==raw['jets'] and experiment.offshell_action()==raw['offshell_action'] and experiment.discrete()==raw['discrete']
    assert pair.run()==paired
    assert [(q['N'],q['steps']) for q in raw['rows']]==[(17,4),(25,4),(25,8),(33,4)]
    assert raw['rows'][-1]['stress_density_divergence_max']<raw['rows'][0]['stress_density_divergence_max']
    assert raw['rows'][1]['current_density_divergence_max']>raw['rows'][0]['current_density_divergence_max']
    files=[HERE/n for n in ('compact_receiver_propagation.py','compact_receiver_checks.py','compact_receiver_checks_results.json','receiver_density_pair_identity.py','receiver_density_pair_identity_results.json','drafts/initial_refinement_failure.json','drafts/receiver_source_derivation.md')]
    return dict(round=906,date='2026-10-06',all_checks_passed=True,cumulative_numbered_groups=3691,
        argument_scope=raw['argument_scope'],independent_spinor_action_and_density_pair_reproduced=True,
        full_grid_time_integrations_previously_executed=True,full_grid_time_integrations_rerun_by_default=False,
        full_check_command='python -B -X utf8 validate_compact_receiver.py --rerun-flow',
        source_file_sha256={str(p.relative_to(HERE)):sha(p) for p in files},
        finest_receiver_diagnostics=raw['rows'][-1],paired_density_identity=paired,
        local_Ward_residual_is_not_a_certified_observation_error=True,initial_failure_preserved=True,
        original872_total_source_computed=False,finite_time_error_certified=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-flow',action='store_true');a=p.parse_args();out=run(a.rerun_flow)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
