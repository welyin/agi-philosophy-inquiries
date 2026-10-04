"""Round 517 reproduction and immutable evidence checks."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(),10000))
import verify_qca_finite_region_scope as previous
import verify_interaction_rounds as core
import all_scale_monitored_graph_source as science

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_517_checks.json'
NAMES=['all_scale_monitored_graph_source.py','all_scale_monitored_graph_source_results.json',
       'research_note_517.md']
DRAFTS=['round517_drafts/research_note_517.txt',
        'round517_drafts/research_note_517_before_phase_and_computability_clarification.txt']


def verify():
    base=previous.verify()
    assert base['stage_saved_tests']==2531
    assert base['science_hashes_verified_231_516']==859
    assert base['total_protected_evidence_hashes']==982
    if TARGET.exists():
        for group in ('new_file_hashes','preserved_draft_hashes'):
            for name,sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name)==sha,name
    assert (HERE/NAMES[2]).read_bytes()==(HERE/DRAFTS[0]).read_bytes()
    assert core.digest(HERE/NAMES[0])=='1676690bb960707f543a83551b644cddd873d68eb488818863747994c39fb694'
    for name in (NAMES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved=core.read(HERE/NAMES[1])
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(5,0,0)
    assert (saved['round'],saved['scientific_base_through_round'])==(517,516)
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    yes=['all_finite_distinct_connected_graph_families_theorem',
         'unchanged_model_with_declared_port_instrument',
         'arbitrary_initial_active_state_and_passive_reference',
         'all_binary_tree_sizes_configuration_connected',
         'no_data_or_graph_reset_between_monitoring_steps']
    no=['uniform_in_size_mixing_budget_proved','full_record_environment_decoupled',
        'coherent_uniform_graph_state_generated','autonomous_measurement_and_clock_generated',
        'stable_macroscopic_position_generated','dimension_three_generated','full_GR_goal_completed']
    assert set(saved['scope'])==set(yes+no)
    assert all(saved['scope'][k] is True for k in yes)
    assert all(saved['scope'][k] is False for k in no)
    checked=core.text_checks(HERE/NAMES[2]);assert checked['display_formulas']==12
    links=0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    return dict(date='2026-09-28',round=517,scientific_base_through_round=516,
        fresh_tests=dict(run=5,failures=0,errors=0),saved_results_reproduced=True,
        scientific_results_rewritten=False,previous_scientific_file_hashes_verified=859,
        previous_protected_evidence_hashes_verified=982,
        new_file_hashes={name:core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name:core.digest(HERE/name) for name in DRAFTS},
        text_checks=checked,local_links_checked=links,broken_links=0,
        independent_final_code_and_draft_review_completed=True,
        review_corrections='Retain SWAP energy constant; restrict certification search to computable coefficients.',
        old_numbered_science_experiments_rerun=False,visual_checks_performed=False,
        scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','fresh_tests','all_reported_checks_passed')}))
