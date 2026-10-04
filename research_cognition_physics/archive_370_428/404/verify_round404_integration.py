"""Integrate round 404 and preserve all prior scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round403_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round404_integration_checks.json'


def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/'research_round_404_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round']==403
    assert checked['batch_scientific_dependencies']==[]
    assert len(checked['new_file_hashes'])==3
    for name,sha in checked['new_file_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'quiescent_program_sector_audit_results.json')
    assert saved['checks']=={'run':6,'failures':0,'errors':0}
    for flag in ('normal_quiescent_sector_explicit',
                 'all_fixed_normal_programs_have_local_reset_limit',
                 'strong_universality_corollary_requires_state_space_clarification'):
        assert saved[flag],flag
    for flag in ('ordinary_finite_task_universality_rejected',
                 'local_reset_implies_global_information_destruction',
                 'convergence_uniform_over_all_program_states',
                 'infinite_classical_program_belongs_to_original_hilbert_sector',
                 'all_infinite_sector_extensions_excluded',
                 'autonomous_continuous_hamiltonian_derived',
                 'spatial_dimension_generated','full_cognition_to_gr_refuted'):
        assert not saved[flag],flag
    formulas=core.text_checks(HERE/'research_note_404.md')['display_formulas']
    assert formulas==12
    links=0
    for path in (HERE/'research_note_404.md',HERE/'spatial_premise_closure_audit.md'):
        for link in core.link_parser()(path.read_text(encoding='utf-8')):
            dest=(path.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),(path,link)
            links+=1
    for name in ('quiescent_program_sector_audit.py','verify_quiescent_program_round.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_403')==520
    assert result.pop('unchanged_prior_science_hashes_231_402')==517
    assert result['stage_saved_tests']==1807
    assert result['total_protected_evidence_hashes']==543
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',(HERE.parent/'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1])>=404 and int(latest[2])>=1813
    result.update(rounds=[404],execution_mode='state-space qualification of sustained physical universality',
                  scientific_base_through_round_by_round={404:403},
                  fresh_tests_by_round={404:6},fresh_tests=6,stage_saved_tests=1813,
                  science_hashes_verified_231_404=523,unchanged_prior_science_hashes_231_403=520,
                  total_protected_evidence_hashes=546,unchanged_prior_evidence_hashes=543,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  formerly_unread_schaeffer_paper_now_checked=True,
                  original_normal_program_sector_and_extended_states_distinguished=True,
                  ordinary_finite_task_universality_rejected=False,
                  all_infinite_sector_extensions_excluded=False,
                  spatial_dimension_generated=False,full_cognition_to_gr_refuted=False,
                  independent_final_agent_review_completed=False,all_reported_checks_passed=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
