"""Verify the complete numbered reports and non-duplicated evidence registration."""
import argparse
import ast
import json
from pathlib import Path
import verify_thermal_direction_bridge as previous
import verify_interaction_rounds as core
from register_round521_523_reports import CONFIG, numbered_content

HERE = Path(__file__).resolve().parent
TARGET = HERE/'round521_523_integration_checks.json'


def verify():
    base = previous.verify()
    assert base['total_protected_including_this_review']==1052
    assert base['unchanged_numbered_scientific_tests']==2552
    assert base['unchanged_scientific_file_hashes']==871
    links, formulas = 0, 0
    newly_created, promoted = {}, {}
    for index,(n,(stem,code,count,scientific_base)) in enumerate(CONFIG.items()):
        checked = core.read(HERE/f'research_round_{n}_checks.json')
        assert checked['registered_tests']==dict(run=6,failures=0,errors=0)
        assert checked['new_scientific_experiments_for_registration']==0
        assert checked['scientific_base_through_round']==scientific_base
        assert checked['cumulative_numbered_tests']==2552+6*(index+1)
        assert checked['cumulative_unique_protected_evidence_files']==1052+2*(index+1)
        for group in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in checked[group].items():
                assert core.digest(HERE/name)==digest,name
        note = HERE/f'research_note_{n}.md'
        draft = HERE/f'round{n}_drafts/research_note_{n}.txt'
        assert note.read_bytes()==draft.read_bytes()==numbered_content(n).encode('utf8')
        assert core.text_checks(note)==checked['text_checks']
        assert checked['text_checks']['display_formulas']==count
        formulas += count
        assert checked['original_full_report_sha256']==core.digest(HERE/(stem+'_review.md'))
        assert checked['scope']==core.read(HERE/(stem+'_results.json'))['scope']
        for link in core.link_parser()(note.read_text('utf8')):
            assert (HERE/link).resolve().exists(),(n,link)
            links += 1
        for path in (note,draft):
            newly_created[path.relative_to(HERE).as_posix()] = core.digest(path)
        for name in (code,stem+'_results.json'):
            promoted[name] = core.digest(HERE/name)
    assert len(newly_created)==6 and len(promoted)==6
    if TARGET.exists():
        old = core.read(TARGET)
        assert old['new_numbered_note_and_draft_hashes']==newly_created
        assert old['previously_protected_source_result_hashes_registered']==promoted
    for name in ('register_round521_523_reports.py',Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-09-30',rounds=list(CONFIG),latest_round=523,
        scientific_base_through_round_by_round={n:row[3] for n,row in CONFIG.items()},
        registered_tests_by_round={n:6 for n in CONFIG},registered_tests=18,
        stage_saved_tests=2570,science_hashes_verified_231_523=880,
        unchanged_prior_science_hashes_231_520=871,
        total_protected_evidence_hashes=1058,unchanged_prior_evidence_hashes=1052,
        tests_transferred_from_prior_unnumbered_reports=18,
        new_scientific_experiments_for_registration=0,
        new_numbered_note_and_draft_hashes=newly_created,
        previously_protected_source_result_hashes_registered=promoted,
        complete_report_bodies_preserved=True,old_reports_code_results_preserved=True,
        new_scientific_display_formulas_checked=formulas,
        local_note_links_checked=links,broken_links=0,
        independent_scientific_review_inherited=True,
        full_result_reproduction_confirmed_by_previous_verifier=True,
        visual_checks_performed=False,active_goal_unchanged=True,
        stage_complete=False,all_reported_checks_passed=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('rounds','stage_saved_tests',
        'total_protected_evidence_hashes','all_reported_checks_passed')}))
