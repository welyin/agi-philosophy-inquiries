"""Reproduce round 520 and preserve all 1004 earlier protected files."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_uniform_tree_continuum_review as previous
import verify_interaction_rounds as core
import selective_record_state_learning as science

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_520_checks.json'
FILES=['selective_record_state_learning.py','selective_record_state_learning_results.json',
       'research_note_520.md']
DRAFTS=['round520_drafts/research_note_520.txt']


def verify():
    base=previous.verify()
    assert base['unchanged_numbered_scientific_tests']==2546
    assert base['unchanged_scientific_file_hashes']==868
    assert base['total_protected_including_this_review']==1004
    if TARGET.exists():
        for group in ('new_file_hashes','preserved_draft_hashes'):
            for name,sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name)==sha,name
    assert (HERE/FILES[2]).read_bytes()==(HERE/DRAFTS[0]).read_bytes()
    for name in (FILES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved=core.read(HERE/FILES[1])
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['round'],saved['scientific_base_through_round'])==(520,519)
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    assert saved['scope']['fixed_finite_size_common_step_no_dark_subspace']
    assert saved['scope']['uniform_mean_full_reference_report_bound']
    assert not saved['scope']['macroscopic_positions_or_three_dimensions_derived']
    assert not saved['scope']['full_GR_goal_completed']
    checked=core.text_checks(HERE/FILES[2]);assert checked['display_formulas']==11
    links=0
    for link in core.link_parser()((HERE/FILES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    return dict(date='2026-09-28',round=520,scientific_base_through_round=519,
        fresh_tests=dict(run=6,failures=0,errors=0),saved_results_reproduced=True,
        scientific_results_rewritten=False,previous_scientific_file_hashes_verified=868,
        previous_protected_evidence_hashes_verified=1004,
        new_file_hashes={name:core.digest(HERE/name) for name in FILES},
        preserved_draft_hashes={name:core.digest(HERE/name) for name in DRAFTS},
        inherited_post519_unnumbered_evidence_files=4,
        text_checks=checked,local_links_checked=links,broken_links=0,
        independent_final_code_and_draft_review_completed=True,
        old_numbered_science_experiments_rerun=False,visual_checks_performed=False,
        scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();answer=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('round','fresh_tests',
        'previous_protected_evidence_hashes_verified','all_reported_checks_passed')}))
