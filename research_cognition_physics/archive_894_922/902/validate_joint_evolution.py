"""Review saved902 evolution and reproduce independent action/geometry checks.
Use --rerun-evolution for the complete approximately three-minute integration.
"""
from pathlib import Path
import argparse,json,math,hashlib
import joint_evolution_checks as experiment
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_background_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strip_times(x):
    if isinstance(x,dict):return {k:strip_times(v) for k,v in x.items() if k!='elapsed_seconds'}
    if isinstance(x,list):return [strip_times(v) for v in x]
    return x
def run(rerun=False):
    raw=json.loads(experiment.TARGET.read_text('utf-8'))
    if rerun:assert strip_times(experiment.run())==strip_times(raw)
    identities=experiment.metric_identity_check();variational=experiment.variational_check()
    assert identities==raw['identities'] and variational==raw['variational']
    rows=raw['rows'];assert [(r['N'],r['steps']) for r in rows]==[(12,8),(16,8),(24,8),(24,16),(32,8)]
    spatial=[r for r in rows if r['steps']==8]
    for key in ('H_max','M_max','Gauss_max','harmonic_max'):
        assert all(b['final'][key]<a['final'][key] for a,b in zip(spatial,spatial[1:])),key
    for row in rows:
        assert row['Tend']==.01 and row['final']['minimum_metric_eigenvalue']>.7 and row['final']['minimum_F']>1.8
        assert row['point_initial']['clock_timelike_margin']>0 and row['point_final']['clock_timelike_margin']>0
    assert max(raw['fixed_N24_time_refinement_differences'].values())<3e-12
    notes=[]
    for row in rows:
        n=row['N'];notes.append(dict(N=n,steps=row['steps'],coordinates=[0,2*math.pi*(n//4)/n,2*math.pi*(n//8)/n],
            same_original_reference_point=(n%8==0)))
    # N12 uses z=pi/6; do not compare that point with the pi/4 reference.
    assert not notes[0]['same_original_reference_point'] and all(a['same_original_reference_point'] for a in notes[1:])
    return dict(round=902,date='2026-10-06',formal_reports=902,cumulative_numbered_groups=3687,fresh_numbered_groups=1,
        all_checks_passed=True,argument_scope='Original859 classical Einstein-Yang-Mills-H5-probe evolution in harmonic/temporal gauges: action, stress and energy consistency plus short-time spatial/time-refinement diagnostics. No rigorous time-error certificate or quantum propagation yet.',
        independent_action_checks_reproduced=True,recorded_evolution_executed=True,
        full_evolution_can_be_reproduced_with='python -B -X utf8 validate_joint_evolution.py --rerun-evolution',
        source_file_sha256={p.name:sha(p) for p in (HERE/'common_background_evolution.py',HERE/'joint_evolution_checks.py',experiment.TARGET)},
        original_reference_point_metadata=notes,latest_diagnostics=rows[-1]['final'],
        time_refinement_differences=raw['fixed_N24_time_refinement_differences'],
        metric_work=variational['energy_work'],energy_exchange_errors=variational['energy_exchange_errors'],
        original_color_weak_abelian_and_six_scalars_retained=True,full_metric_evolved=True,
        floating_diagnostics_not_error_enclosure=True,rigorous_time_error_or_future_reference_patch=False,
        quantum_modes_or_source_budget_computed=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-evolution',action='store_true');a=p.parse_args();r=run(a.rerun_evolution)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
