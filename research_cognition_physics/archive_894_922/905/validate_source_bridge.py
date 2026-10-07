"""905 independent source checks and saved Ward diagnostics."""
from pathlib import Path
import argparse,hashlib,json
import source_bridge_checks as experiment
HERE=Path(__file__).resolve().parent;TARGET=HERE/'source_bridge_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(rerun=False):
    raw=json.loads(experiment.TARGET.read_text('utf-8'))
    if rerun:assert experiment.run()==raw
    assert experiment.independent_checks()==raw['independent']
    assert [r['N'] for r in raw['rows']]==[12,16,24]
    for key in ('diffeomorphism_density_max','internal_density_max'):
        assert all(b[key]<a[key] for a,b in zip(raw['rows'],raw['rows'][1:]))
    return dict(round=905,date='2026-10-06',all_checks_passed=True,cumulative_numbered_groups=3690,
        argument_scope=raw['argument_scope'],independent_action_and_constraint_checks_reproduced=True,
        full_Ward_grid_checks_previously_executed=True,full_Ward_checks_rerun_by_default=False,
        full_check_command='python -B -X utf8 validate_source_bridge.py --rerun-ward',
        source_file_sha256={p.name:sha(p) for p in (HERE/'canonical_source_bridge.py',HERE/'source_bridge_checks.py',experiment.TARGET)},
        finest_joint_Ward_diagnostics=raw['rows'][-1],
        original872_quantum_source_computed=False,retarded_solution_numerically_computed=False,
        rigorous_time_error_certificate=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-ward',action='store_true');a=p.parse_args();r=run(a.rerun_ward)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
