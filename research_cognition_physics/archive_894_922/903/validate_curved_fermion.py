"""903 fast validation; --rerun-propagation repeats full background plus Dirac integration."""
from pathlib import Path
import argparse,hashlib,json
import curved_fermion_checks as experiment
HERE=Path(__file__).resolve().parent;TARGET=HERE/'curved_fermion_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(rerun=False):
    raw=json.loads(experiment.TARGET.read_text('utf-8'))
    if rerun:assert experiment.run()==raw
    assert experiment.frame_checks()==raw['frame']
    assert experiment.discrete_checks()==raw['discrete']
    assert [(r['N'],r['steps']) for r in raw['rows']]==[(9,8),(13,8),(13,16),(17,8)]
    assert all(r['gram_error']<2e-10 for r in raw['rows'])
    assert raw['time_halving_N13_max']<1e-10
    errors=raw['common_band_errors_9_13_and_13_17'];assert errors[1]<errors[0]
    return dict(round=903,date='2026-10-06',all_checks_passed=True,cumulative_numbered_groups=3688,
        argument_scope=raw['argument_scope'],independent_frame_and_CAR_checks_reproduced=True,
        recorded_full_propagation_executed=True,full_propagation_rerun_by_default=False,
        full_propagation_command='python -B -X utf8 validate_curved_fermion.py --rerun-propagation',
        source_file_sha256={p.name:sha(p) for p in (HERE/'curved_fermion_propagation.py',HERE/'curved_fermion_checks.py',experiment.TARGET)},
        max_gram_drift=max(r['gram_error'] for r in raw['rows']),
        time_halving_N13_max=raw['time_halving_N13_max'],common_band_errors=errors,
        actual_prepared_covariance_computed=False,actual_quantum_source_computed=False,
        rigorous_continuous_error_bound_computed=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-propagation',action='store_true');a=p.parse_args();r=run(a.rerun_propagation)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
