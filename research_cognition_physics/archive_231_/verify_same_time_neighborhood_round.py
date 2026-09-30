"""Read-only round 518 scientific reproduction and evidence verification."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(),10000))
import verify_round517_integration as previous
import verify_interaction_rounds as core
import same_time_neighborhood_source as science

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_518_checks.json'
NAMES=['same_time_neighborhood_source.py','same_time_neighborhood_source_results.json',
       'research_note_518.md']
DRAFTS=['round518_drafts/research_note_518.txt',
        'round518_drafts/research_note_518_before_delayed_event_clarification.txt']


def verify():
    base=previous.verify()
    assert base['stage_saved_tests']==2536
    assert base['science_hashes_verified_231_517']==862
    assert base['total_protected_evidence_hashes']==987
    if TARGET.exists():
        for group in ('new_file_hashes','preserved_draft_hashes'):
            for name,sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name)==sha,name
    assert (HERE/NAMES[2]).read_bytes()==(HERE/DRAFTS[0]).read_bytes()
    for name in (NAMES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved=core.read(HERE/NAMES[1])
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    assert (saved['round'],saved['scientific_base_through_round'])==(518,517)
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    yes=['unchanged_440_H_and_declared_port_instrument','same_final_time_joint_radius_two_source',
         'arbitrary_unknown_graph_and_passive_reference','size_independent_success_and_error_bounds',
         'records_and_failure_instrument_retained']
    no=['per_record_posterior_error_guarantee','automatic_failure_reset',
        'arbitrary_large_buffer_generated','good_domain_predictor_generated',
        'autonomous_readout_clock_or_delivery_generated','physical_three_dimensional_space_generated',
        'full_GR_goal_completed']
    assert set(saved['scope'])==set(yes+no)
    assert all(saved['scope'][k] is True for k in yes)
    assert all(saved['scope'][k] is False for k in no)
    checked=core.text_checks(HERE/NAMES[2]);assert checked['display_formulas']==16
    links=0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    return dict(date='2026-09-28',round=518,scientific_base_through_round=517,
        fresh_tests=dict(run=6,failures=0,errors=0),saved_results_reproduced=True,
        scientific_results_rewritten=False,previous_scientific_file_hashes_verified=862,
        previous_protected_evidence_hashes_verified=987,
        new_file_hashes={name:core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name:core.digest(HERE/name) for name in DRAFTS},
        text_checks=checked,local_links_checked=links,broken_links=0,
        independent_final_code_and_draft_review_completed=True,
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
