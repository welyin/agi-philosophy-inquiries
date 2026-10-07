"""907: reproduce the finite gauge dictionary audit and inspect saved material jets."""
from pathlib import Path
import argparse,hashlib,json
import connection_dictionary_audit as audit
HERE=Path(__file__).resolve().parent;TARGET=HERE/'material_dictionary_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(rerun=False):
    loops=audit.run();assert loops==json.loads(audit.TARGET.read_text('utf-8'))
    raw=json.loads((HERE/'material_reference_jets_results.json').read_text('utf-8'))
    cross=json.loads((HERE/'reference_chart_crosscheck_results.json').read_text('utf-8'))
    if rerun:
        import material_reference_jets as jets
        import reference_chart_crosscheck as refine
        assert jets.run()==raw and refine.run()==cross
    assert raw['formal_previous']==cross['formal_previous']==906
    assert [(r['N'],r['steps']) for r in cross['rows']]==[(16,4),(24,2),(24,4),(32,4)]
    assert raw['rows'][0]['data']['new']['determinant']>0
    assert all(r['data']['new']['determinant']<0 and r['data']['old']['determinant']>0 for r in cross['rows'])
    for row in raw['rows']:
        assert row['clock_margin']>0
        for k in ('old','new'):
            e=row['data'][k]['time_direction_central_difference_errors']
            assert all(3.8<a/b<4.2 for a,b in zip(e,e[1:]))
    assert not cross['continuous_zero_certified'] and not cross['common_physics_failure_proved']
    names=('connection_dictionary_audit.py','connection_dictionary_audit_results.json',
        'material_reference_jets.py','material_reference_jets_results.json',
        'reference_chart_crosscheck.py','reference_chart_crosscheck_results.json',
        'drafts/material_source_implementation_audit.md','drafts/working_checks.json')
    return dict(round=907,date='2026-10-06',all_checks_passed=True,cumulative_numbered_groups=3692,
        argument_scope='Original finite relational gauge dictionary; exact local counterexample and coherent repair. Actual material chart diagnostics, not a certified time slab or full source.',
        matrix_gauge_and_tangent_checks_reproduced=True,material_grid_time_checks_previously_executed=True,
        material_grid_time_checks_rerun_by_default=False,
        full_check_command='python -B -X utf8 validate_material_dictionary.py --rerun-material',
        correct_adapter_maximum_pure_gauge_record_change=loops['correct_adapter_maximum_pure_gauge_record_change'],
        wrong_adapter_maximum_bounded_record_change=loops['legacy_literal_adapter_maximum_pure_gauge_record_change']/9,
        exact_wrong_adapter_leading_coefficient=loops['exact_wrong_adapter_adjoint_derivative_leading_coefficient'],
        initial_new_reference_determinant=raw['rows'][0]['data']['new']['determinant'],
        final_refined_new_reference_determinants=[r['data']['new']['determinant'] for r in cross['rows']],
        source_file_sha256={n:sha(HERE/n) for n in names},
        continuous_chart_failure_proved=False,abstract_Wilson_theorem_refuted=False,
        original872_total_source_computed=False,finite_time_error_certified=False,
        effective_description_stage_acceptance_allowed=True,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-material',action='store_true');a=p.parse_args();out=run(a.rerun_material)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
