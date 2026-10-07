"""904 saved mixed tangent diagnostics plus independently repeated small-grid checks."""
from pathlib import Path
import argparse,json,hashlib
import boson_tangent_checks as experiment
HERE=Path(__file__).resolve().parent;TARGET=HERE/'boson_tangent_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(rerun=False):
    raw=json.loads(experiment.TARGET.read_text('utf-8'))
    if rerun:assert experiment.run()==raw
    assert experiment.checks()==raw['checks']
    assert [(r['N'],r['steps']) for r in raw['rows']]==[(12,8),(16,8),(24,8)]
    for key in ('H','M','harmonic'):
        assert all(b['final_linear_constraints'][key]<a['final_linear_constraints'][key] for a,b in zip(raw['rows'],raw['rows'][1:])),key
    assert raw['rows'][-1]['final_linear_constraints']['Gauss']<1e-10
    assert all(r['induced_H5']>1e-9 and r['final_tangent_max']['A']>1e-8 for r in raw['rows'])
    return dict(round=904,date='2026-10-06',all_checks_passed=True,cumulative_numbered_groups=3689,
        argument_scope=raw['argument_scope'],independent_small_grid_checks_reproduced=True,
        recorded_full_tangent_propagation_executed=True,full_three_grid_propagation_rerun_by_default=False,
        full_propagation_command='python -B -X utf8 validate_boson_tangent.py --rerun-propagation',
        source_file_sha256={p.name:sha(p) for p in (HERE/'coupled_boson_tangent.py',HERE/'boson_tangent_checks.py',experiment.TARGET)},
        finest_grid_linear_constraints=raw['rows'][-1]['final_linear_constraints'],
        source_and_constraint_contract_for_retarded_inverse_complete=False,
        quantum_source_or_backreaction_computed=False,rigorous_continuous_error_bound_computed=False,
        full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-propagation',action='store_true');a=p.parse_args();r=run(a.rerun_propagation)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
