"""Read-only scientific and frozen-evidence verification for round 516."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(),10000))
import verify_round515_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_516_checks.json'
NAMES=['memory_assisted_target_reception.py','memory_assisted_target_reception_results.json',
       'research_note_516.md']
DRAFTS=['round516_drafts/research_note_516.txt',
        'round516_drafts/memory_assisted_target_reception_before_budget_fix.txt']


def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:base=previous.verify(pending)
    finally:previous.TARGET=old_target
    assert base['stage_saved_tests']==2525
    assert base['science_hashes_verified_231_515']==856
    assert base['total_protected_evidence_hashes']==970
    if TARGET.exists():
        for group in ['new_file_hashes','preserved_draft_hashes']:
            for name,sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name)==sha,name
    assert (HERE/NAMES[2]).read_bytes()==(HERE/DRAFTS[0]).read_bytes()
    for name in [NAMES[0],Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf8'))
    spec=importlib.util.spec_from_file_location('memory_target_reception_checked',HERE/NAMES[0])
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    saved=core.read(HERE/NAMES[1])
    assert saved==json.loads(json.dumps(module.run()))
    assert (saved['round'],saved['scientific_base_through_round'])==(516,515)
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    true=['same_H_source_and_reception','finite_all_input_target_bound_proved',
          'complete_memory_monitoring_and_target_detection_required','no_graph_or_memory_reset_between_cycles']
    false=['full_quantum_environment_same_reduced_channel_claimed','all_network_sizes_numerically_certified',
           'autonomous_readout_and_clock_generated','stable_macroscopic_position_generated',
           'dimension_three_generated','full_GR_goal_completed']
    assert set(saved['scope'])==set(true+false)
    assert all(saved['scope'][k] is True for k in true)
    assert all(saved['scope'][k] is False for k in false)
    checked=core.text_checks(HERE/NAMES[2]);assert checked['display_formulas']==14
    links=0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    return dict(date='2026-09-28',round=516,scientific_base_through_round=515,
        frozen_integration_base_through_round=515,fresh_tests=dict(run=6,failures=0,errors=0),
        saved_results_reproduced=True,scientific_results_rewritten=False,
        previous_scientific_file_hashes_verified=856,previous_protected_evidence_hashes_verified=970,
        text_checks=checked,local_links_checked=links,broken_links=0,
        new_file_hashes={name:core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name:core.digest(HERE/name) for name in DRAFTS},
        visual_rendering_performed=False,legacy_science_tests_rerun=False,
        independent_review='Independent same-H memory history, no-dark-space, fixed-point error/overflow, source-composition, final code and 14-equation draft audit.',
        scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','fresh_tests','all_reported_checks_passed')},ensure_ascii=False))
