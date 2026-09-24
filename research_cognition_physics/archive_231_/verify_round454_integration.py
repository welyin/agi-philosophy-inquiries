"""Read-only integration of round 454; earlier scientific evidence is immutable."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round453_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round454_integration_checks.json'

def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/'research_round_454_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round']==453
    assert checked['fresh_tests']==dict(run=6,failures=0,errors=0)
    assert len(checked['new_file_hashes'])==3
    assert checked['preserved_draft_hashes']=={}
    for name,sha in checked['new_file_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'shared_relation_revival_audit_results.json')
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    assert saved['scope']==checked['scope']
    for name in ["enlarged_common_interface_positive_construction","fixed_equal_positive_pair_exchanges_only","exact_unknown_logic_endpoint_entangling_map","finite_window_raw_physical_readout_certified","shared_gauge_preparation_input_explicit"]:
        assert saved['scope'][name],name
    for name in ["arbitrary_unknown_gauge_joining_implemented","autonomous_resource_preparation_or_permanent_halting","complete_SoCA_or_quantum_necessity_derived","physical_spatial_dimension_or_GR_generated","full_GR_goal_completed","phase_closure_triggered"]:
        assert not saved['scope'][name],name
    note=HERE/'research_note_454.md'
    assert core.text_checks(note)==checked['text_checks']
    assert checked['text_checks']['display_formulas']==12
    for name in ['shared_relation_revival_audit.py','verify_shared_relation_revival_round.py',Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    links=0
    for doc in [note,HERE/'spark_economy_cognition_bridge_review.md']:
        for link in core.link_parser()(doc.read_text(encoding='utf-8')):
            dest=(doc.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),link
            links+=1
    assert result.pop('science_hashes_verified_231_453')==670
    assert result.pop('unchanged_prior_science_hashes_231_452')==667
    assert result['stage_saved_tests']==2151
    assert result['total_protected_evidence_hashes']==705
    direction=(HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',direction)
    assert latest and int(latest[1])>=454 and int(latest[2])>=2157
    assert '连续自然演化' in direction and '统一内生演化' in direction
    assert '## 82. 第454轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24',rounds=[454],
        execution_mode='independent positive interface construction based on round 453',
        scientific_base_through_round_by_round={454:453},
        fresh_tests_by_round={454:6},fresh_tests=6,
        stage_saved_tests=2157,science_hashes_verified_231_454=673,
        unchanged_prior_science_hashes_231_453=670,
        total_protected_evidence_hashes=708,
        unchanged_prior_evidence_hashes=705,newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        independent_rounds_454_455_share_science_baseline_453=True,
        positive_interface_scope=saved['scope'],
        spark_source_comparison_recorded_without_new_science_test_count=True,
        spark_project_files_modified=False,
        original_229_theorem_not_rewritten=True,
        physical_positions_or_dimension_generated=False,
        full_cognition_to_gr_refuted=False,phase_closure_triggered=False,
        old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=True,all_reported_checks_passed=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
