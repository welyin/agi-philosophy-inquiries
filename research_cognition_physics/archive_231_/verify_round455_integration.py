"""Read-only integration of round 455; earlier scientific evidence is immutable."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round454_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round455_integration_checks.json'

def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/'research_round_455_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round']==453
    assert checked['fresh_tests']==dict(run=6,failures=0,errors=0)
    assert len(checked['new_file_hashes'])==3
    assert checked['preserved_draft_hashes']=={}
    for name,sha in checked['new_file_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'relational_interface_conversion_audit_results.json')
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    assert saved['scope']==checked['scope']
    for name in ["existing_three_spin_unknown_G_L_input","invariant_singlet_aux_resource_explicit","exact_three_constant_exchange_conversion","every_internal_reference_preserved","finite_window_full_channel_error_proved","minimum_carrier_and_closed_five_spin_purity_conditions_proved","no_new_abstract_encoding_universality_claimed"]:
        assert saved['scope'][name],name
    for name in ["depends_on_round_454","singlet_preparation_source_derived","permanent_autonomous_handoff","macro_subject_entangling_network_implemented","full_soca_or_spatial_dimension_or_GR_generated","full_GR_goal_completed","phase_closure_triggered"]:
        assert not saved['scope'][name],name
    note=HERE/'research_note_455.md'
    assert core.text_checks(note)==checked['text_checks']
    assert checked['text_checks']['display_formulas']==13
    for name in ['relational_interface_conversion_audit.py','verify_relational_interface_conversion_round.py',Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    links=0
    for doc in [note,HERE/'spark_economy_cognition_bridge_review.md']:
        for link in core.link_parser()(doc.read_text(encoding='utf-8')):
            dest=(doc.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),link
            links+=1
    assert result.pop('science_hashes_verified_231_454')==673
    assert result.pop('unchanged_prior_science_hashes_231_453')==670
    assert result['stage_saved_tests']==2157
    assert result['total_protected_evidence_hashes']==708
    direction=(HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',direction)
    assert latest and int(latest[1])>=455 and int(latest[2])>=2163
    assert '连续自然演化' in direction and '统一内生演化' in direction
    assert '## 83. 第455轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24',rounds=[455],
        execution_mode='independent positive interface construction based on round 453',
        scientific_base_through_round_by_round={455:453},
        fresh_tests_by_round={455:6},fresh_tests=6,
        stage_saved_tests=2163,science_hashes_verified_231_455=676,
        unchanged_prior_science_hashes_231_454=673,
        total_protected_evidence_hashes=711,
        unchanged_prior_evidence_hashes=708,newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=13,
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
