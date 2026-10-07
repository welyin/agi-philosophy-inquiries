"""Read-only adoption audit, except for an explicit first --write receipt."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
ROOT = STAGE.parent.parent
TARGET = HERE / 'common_adoption_audit_checks.json'
DOC = HERE / 'research_note_981_working.md'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text('utf-8-sig'))


def run():
    frozen = {}
    def add(items):
        for rel, digest in items.items():
            assert rel not in frozen or frozen[rel] == digest, rel
            frozen[rel] = digest
    for n in range(776, 981):
        receipt = read(STAGE / f'{n}/research_round_{n}_checks.json')
        for key in ('frozen_inputs', 'new_scientific_and_entry_files'):
            add(receipt[key])
    for rel, digest in frozen.items():
        assert sha(ROOT / rel) == digest, rel

    evidence = [STAGE / f'research_note_{n}.md' for n in
                (854,898,899,935,936,938,939,953,957,959,961,965,969,
                 971,973,974,975,976,977,978,979,980)]
    evidence += [STAGE / r for r in (
        '962/drafts/revised_goal_20261007.txt',
        '962/drafts/app_goal_confirmation.json',
        '957/drafts/unified_operation_hypotheses_v0_2.md',
        '953/drafts/priority_reaudit.md',
        '971/drafts/native_parent_map_v0_1.md',
        '939/drafts/common_model_recovery_map.md',
        '980/research_round_980_checks.json',
        '980/finite_thermal_records_results.json',
        '981/drafts/STATUS.md')]
    thermal = read(STAGE / '980/finite_thermal_records_results.json')
    assert thermal['channel']['contacts'] == 600000
    assert thermal['channel']['full_charge_reset_bound'] < .001029
    receipt = read(STAGE / '980/research_round_980_checks.json')
    assert receipt['formal_reports'] == 980
    assert receipt['cumulative_numbered_test_groups_from_979'] == 3765

    prose = DOC.read_text('utf-8')
    for term in ('不是新定理', '正式轮次仍为980', '3765', '不启动',
                 '目标保持active', '共同父模型采用合同', '误差合同'):
        assert term in prose, term
    links = 0
    for target in re.findall(r'\]\(([^)]+)\)', prose):
        if re.match(r'https?://', target):
            continue
        target = target.split('#', 1)[0]
        if not target:
            continue
        p = (DOC.parent / unquote(target)).resolve()
        assert p.exists() or p == TARGET, str(p)
        links += 1

    nav = [STAGE.parent / n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    nav += [STAGE / n for n in ('README.md','文件索引.md','阶段成果总览.md',
                               '跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    heading = '## 980后采用复核：先合并父模型，停止局部扩建'
    for p in nav:
        s = p.read_text('utf-8-sig')
        assert s.count(heading) == 1, str(p)
        assert '981/drafts/research_note_981_working.md' in s
        assert '正式980／累计3765' in s
    formal = [p for p in STAGE.glob('research_note_*.md')
              if re.fullmatch(r'research_note_\d+\.md',p.name)]
    assert max(int(re.search(r'\d+',p.name).group()) for p in formal) == 980
    return dict(audit='common_parent_adoption_after_980',date='2026-10-07',
        all_audit_checks_passed=True,formal_reports=980,
        cumulative_numbered_test_groups=3765,new_scientific_test_groups=0,
        frozen_history_files_verified=len(frozen),local_document_links_checked=links,
        living_navigation_files_checked=len(nav),
        app_goal_changed=False,full_goal_completed=False,
        visual_checks_performed=False,
        frozen_evidence_hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence},
        audit_file_hashes={str(p.relative_to(ROOT)):sha(p) for p in
                           (DOC,Path(__file__),HERE/'publish_adoption_audit.py')})


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(), 'Frozen receipt already exists'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert read(TARGET)==result
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('hashes')},ensure_ascii=False))
