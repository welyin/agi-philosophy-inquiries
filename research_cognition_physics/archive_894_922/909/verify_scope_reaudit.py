"""Evidence-only audit; no new physics test or numerical certification."""
from pathlib import Path
import argparse, ast, hashlib, json, re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'drafts/scope_reaudit_checks.json'
NOTE = HERE / 'drafts/effective_scope_reaudit_908.md'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(writing=False):
    evidence = [STAGE / f'research_note_{n}.md' for n in
                (894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908)]
    evidence += [STAGE / '899/drafts/effective_scope_joint_audit.md',
                 STAGE / '908/relational_source_validation_results.json',
                 STAGE / '908/research_round_908_checks.json',
                 HERE / 'drafts/retarded_source_working.md',
                 HERE / 'drafts/working_checks.json']
    preserved = {}
    receipt = json.loads((STAGE / '908/research_round_908_checks.json').read_text('utf-8'))
    working = json.loads((HERE / 'drafts/working_checks.json').read_text('utf-8'))
    for table in (receipt['frozen_inputs'], receipt['new_scientific_and_entry_files'], working['files']):
        for name, value in table.items():
            assert digest(ROOT / name) == value, name
            preserved[name] = value
    v = json.loads((STAGE / '908/relational_source_validation_results.json').read_text('utf-8'))
    assert v['original_first_jet_weak_source_implemented'] is True
    for key in ('strong_spacetime_source_array_computed',
                'full_support_and_dynamical_error_certified',
                'original872_forced_feedback_computed', 'full_goal_completed'):
        assert v[key] is False, key
    assert working['formal_rounds'] == 908
    assert working['cumulative_numbered_groups'] == 3693
    assert working['actual_shifted_source_solver_implemented'] is False
    assert working['original872_forced_response_computed'] is False
    assert not (STAGE / 'research_note_909.md').exists()
    docs = [NOTE]
    docs += [STAGE.parent / n for n in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')]
    docs += [STAGE / n for n in ('README.md', '文件索引.md', '阶段成果总览.md', '跨阶段主题索引.md')]
    docs += [STAGE / '_shared/notes/unified_physics_condition_ledger_current.md']
    count = 0
    for doc in docs:
        prose = re.sub(r'\$\$.*?\$\$', '', doc.read_text('utf-8-sig'), flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)', prose):
            if re.match(r'^[a-zA-Z]+://', link) or link.startswith('#'):
                continue
            target = (doc.parent / link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target == TARGET.resolve()), (doc, link)
            count += 1
    ast.parse(Path(__file__).read_text('utf-8'))
    files = {str(x.relative_to(ROOT)): digest(x) for x in evidence + [NOTE, Path(__file__)]}
    if not writing:
        old = json.loads(TARGET.read_text('utf-8'))
        assert old['evidence_and_audit_hashes'] == files
        assert old['preserved_files'] == preserved
    return dict(date='2026-10-06', evidence_and_scope_checks_passed=True,
                scientific_experiments_rerun=False, new_physics_certificate=False,
                new_scientific_groups=0, formal_rounds=908, cumulative_numbered_groups=3693,
                full_goal_completed=False, local_links_checked=count,
                preserved_files_count=len(preserved), preserved_files=preserved,
                evidence_and_audit_hashes=files,
                judgement=dict(controlled_effective_description_admissible=True,
                               current_stage_completed=False,
                               microscopic_continuity_assumed=False,
                              _physical_minimum_scale_assumed=False,
                               uncut_completion_required_for_stage=False,
                               finite_predictive_source_and_feedback_errors_still_required=True,
                               strong_source_array_and_weak_adjoint_both_required=False))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run(args.write)
    if args.write:
        assert not TARGET.exists(), 'Preserve the previous audit receipt.'
        TARGET.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items()
                      if k not in ('preserved_files', 'evidence_and_audit_hashes')},
                     ensure_ascii=False, indent=2))
