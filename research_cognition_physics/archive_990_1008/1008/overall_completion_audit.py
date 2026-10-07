"""Reproduce coverage/evidence checks; semantic and physical judgments are separate."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'overall_completion_results.json'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(path):
    return path.read_text(encoding='utf-8-sig')


def run():
    receipt = STAGE / '1007/research_round_1007_checks.json'
    previous = read(receipt)
    assert previous['round'] == 1007 and previous['all_delivery_checks_passed']
    assert previous['cumulative_numbered_test_groups_from_1006'] == 3786
    main = HERE / 'overall_operation_hypothesis_v2_2.md'
    old = STAGE / '1005/overall_operation_hypothesis_v2_1.md'
    prose = text(main)
    matches = list(re.finditer(r'^### (P\d\d) ([^\n]+)$', prose, re.M))
    assert [m[1] for m in matches] == [f'P{i:02d}' for i in range(1, 17)]
    old_titles = re.findall(r'^### (P\d\d) ([^\n]+)$', text(old), re.M)
    assert [(m[1], m[2]) for m in matches] == old_titles
    keys = ('参与者与过程', '中间机制', '资源、记录与反作用', '输入与范围')
    coverage = {}
    for idx, match in enumerate(matches):
        end = matches[idx + 1].start() if idx + 1 < len(matches) else prose.index('## 4.')
        section = prose[match.end():end]
        items = re.findall(r'^- \*\*([^\n]+?)：\*\* (.+)$', section, re.M)
        assert tuple(k for k, _ in items) == keys, match[1]
        assert all(value.strip() for _, value in items)
        coverage[match[1]] = {'title': match[2], 'obligations': dict(items)}

    block = prose.split('## 5. 跨部门共同接口：同一对象不能换账', 1)[1].split('## 6.', 1)[0]
    interfaces = [line.split('|')[1:4] for line in block.splitlines()
                  if line.startswith('|') and not line.startswith(('|同一对象', '|---'))]
    assert len(interfaces) == 9 and all(all(row) for row in interfaces)
    audit = text(HERE / 'completion_audit.md')
    requirements = re.findall(r'^\|(R\d\d)\|([^\n]+)$', audit, re.M)
    assert [key for key, _ in requirements] == [f'R{i:02d}' for i in range(1, 10)]
    assert all(len(row.split('|')) == 4 for _, row in requirements)
    assert '当前解释层交付完成' in audit
    assert '它不证明所有候选联合存在' in audit
    assert '必须修正的科学内容0项' in text(HERE / 'drafts/review_notes.md')
    assert 'get_goal' in text(HERE / 'drafts/current_objective_20261008.md')
    for term in ('材料相', '声', '弹性', '输运', '噪声', 'Meissner'):
        assert term in str(coverage['P08']), term
    for term in ('不对称', '生成', '洗出', '探测', '普通物质'):
        assert term in str(coverage['P13']), term
    for term in ('386或425', '332', '333', '有限资源主链', '不新增认知公理'):
        assert term in prose, term

    evidence = [Path(__file__), main, old, receipt,
        HERE / 'completion_audit.md', HERE / 'drafts/current_objective_20261008.md',
        HERE / 'drafts/review_notes.md', HERE / 'follow_on_questions.md',
        HERE / 'drafts/STATUS.md', STAGE / 'research_note_1008.md',
        STAGE / '1006/conventional_superconductivity_adoption_v1.md',
        STAGE / '1007/neutrino_observation_adoption_v1.md',
        STAGE / '993/common_candidate_v1.md', STAGE / 'research_note_991.md',
        STAGE.parent / 'archive_301_341/research_note_332.md',
        STAGE.parent / 'archive_301_341/research_note_333.md']
    return dict(round=1008, date='2026-10-08', kind='overall_mechanism_completion_audit',
        all_audit_checks_passed=True, fresh_physical_test_groups=0, cumulative_test_groups=3786,
        phenomenon_count=len(coverage),
        explanation_obligation_count=sum(len(v['obligations']) for v in coverage.values()),
        interface_count=len(interfaces), requirement_count=len(requirements),
        coverage=coverage, interfaces=interfaces, requirements=dict(requirements),
        current_explanatory_deliverable_completed=True,
        completion_basis='explicit_current_objective_and_independent_semantic_review',
        automated_audit_scope='document_structure_retained_coverage_and_evidence_hashes_only',
        new_cognitive_axioms=False, physical_unification_certified=False,
        independent_physical_derivation_certified=False,
        full_empirical_recovery_certified=False, common_realization_certified=False,
        visual_checks_performed=False,
        source_hashes={str(path.relative_to(ROOT)): sha(path) for path in evidence})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf-8') as dest:
            json.dump(result, dest, ensure_ascii=False, indent=2)
            dest.write('\n')
    else:
        assert result == read(TARGET)
    print(json.dumps({k: v for k, v in result.items()
        if k not in ('coverage', 'interfaces', 'requirements', 'source_hashes')},
        ensure_ascii=False, indent=2))
