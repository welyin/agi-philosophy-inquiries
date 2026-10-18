"""One-time conditional-stage acceptance; preserve prior scientific bytes."""
from pathlib import Path
import ast
import hashlib
import json
import re
import sys
from urllib.parse import unquote

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
RESEARCH = STAGE.parent
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import relocation_original_bytes


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    acceptance = STAGE/'_shared/stage_acceptance.json'
    assert not acceptance.exists(), 'Already finalized; do not repeat'
    review_path = STAGE/'_shared/stage_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    paper = RESEARCH/'14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md'
    synthesis = STAGE/'_shared/joint_selection_synthesis.md'
    goal = STAGE/'_shared/goal.md'
    for p in (paper, synthesis, goal):
        assert sha(p.read_bytes()) in review, 'Final review must bind final document: '+str(p)
    binds, receipts = 0, []
    files = [paper, synthesis, goal, review_path]
    for n in range(1308, 1313):
        p = STAGE/str(n)/'acceptance.json'
        obj = json.loads(p.read_text(encoding='utf8'))
        assert obj['status'] == 'ACCEPTED_CONDITIONAL_SCOPE'
        assert obj['new_axioms'] == 0
        for rel, digest in obj['sha256'].items():
            assert sha((ROOT/rel).read_bytes()) == digest, rel
            binds += 1
        for dep, digest in obj.get('depends_on', {}).items():
            assert sha((STAGE/dep/'acceptance.json').read_bytes()) == digest
        receipts.extend(dict(round=n, **r) for r in obj['executed_checks'])
        assert all(r['exit_code'] == 0 for r in obj['executed_checks'])
        files.append(p)
    # Future acceptance is the only permitted forward local artifact link.
    markdown = list(STAGE.glob('research_note_*.md')) + [paper, synthesis, goal, review_path, STAGE/'README.md']
    for n in range(1308, 1313):
        markdown.extend((STAGE/str(n)).glob('*.md'))
        markdown.append(STAGE/f'_admission/round{n}_admission.md')
    links = 0
    for p in markdown:
        text = p.read_text(encoding='utf8')
        assert len(re.findall(r'^\s*\$\$\s*$', text, flags=re.M)) % 2 == 0, p
        assert not re.search(r'^\s*\\[\[\]]\s*$', text, flags=re.M), p
        for address in re.findall(r'\]\(([^\n)]+)\)', text):
            address = unquote(address.strip('<>').split('#')[0])
            if not address or '://' in address or address.startswith('mailto:'):
                continue
            dest = (p.parent/address).resolve()
            assert dest == acceptance or dest.exists(), (p, address)
            links += 1
    scripts = list(STAGE.rglob('*.py'))
    for p in scripts:
        ast.parse(p.read_text(encoding='utf-8-sig'), filename=str(p))
    layer_path = RESEARCH/'_migration/conditional_selection_active_1308/layer.json'
    layer_raw = layer_path.read_bytes()
    layer = json.loads(layer_raw)
    assert layer['latest_accepted_round'] == 1312
    prefix = ('\n## 已完成阶段：单参考质量与有限交换任务下的条件性结构选择\n\n'
              '2026-10-18，1308—1312五轮已按公开范围验收，阶段联合审阅通过。'
              '[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md) · '
              '[阶段索引](archive_1308_/README.md) · '
              '[联合验收](archive_1308_/_shared/stage_acceptance.json)。新公理0。\n\n'
              '单参考、球面传递、右单态及完整左Dirac质量条件选择SU(2)标准活动像；'
              '单位交换、保护、相位任务与有限实际交付给颜色—代际分类及耦合—Higgs参数域。'
              '三色三／四／五代在当前任务合同中的整域相同，继续调同一预算不能独选三代。\n\n'
              '[共同输入与剩余自由](archive_1308_/_shared/joint_selection_synthesis.md)区分新增推论、'
              '严格保留的自由和未分类输入。条件性阶段目标完成，不表示唯一生成现实标准模型、'
              '全部参数或万有理论；完整物质表、超荷及其它物理选择仍待新机制。\n\n'
              '两旧目录已收尾为archive_1288_1289与archive_1290_1307，'
              '[历史复算与迁移](_migration/closure_1288_1307/README.md)保持。'
              '下方为历史时点记录，当前阶段状态以本节及验收为准。\n')
    updates = []
    for entry in layer['entries']:
        path = ROOT/entry['destination']
        old = path.read_bytes()
        original = relocation_original_bytes(old, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = original.decode(enc)
        assert len(entry['edits']) == 1 and entry['edits'][0]['before'] == ''
        at = entry['edits'][0]['start']
        new = (text[:at]+prefix+text[at:]).encode(enc)
        entry['edits'][0]['after'] = prefix
        entry['current_sha256'], entry['current_bytes'] = sha(new), len(new)
        assert relocation_original_bytes(new, entry) == original
        updates.append((path, old, new))
    layer.update(stage_complete=True, in_progress_round=None,
                 current_work='conditional stage completed; final filesystem verification',
                 application_goal_status_last_observed='active')
    evidence = dict(schema='conditional_selection_stage_acceptance_v1', date='2026-10-18',
                    status='COMPLETE_CONDITIONAL_SCOPE', new_axioms=0,
                    scientific_rounds=list(range(1308, 1313)),
                    accepted_artifact_bindings=binds,
                    completion_criteria=dict(declared_candidates_and_inputs=True,
                        G_open_activity_class=True, N_integer_domain_and_degeneracy=True,
                        P_true_kernel_domain_and_remaining_family=True,
                        J_same_parent_effective_nonempty=True,
                        new_cross_line_connection=True, paper_and_independent_review=True),
                    execution_receipts=receipts,
                    validation=dict(markdown_files=len(markdown), local_links=links,
                                    python_ast_files=len(scripts), formula_delimiters=True),
                    scope_limits=['declared single-reference chiral mass and color/charge role classes',
                                  'Born task, one-loop all-active finite window, declared error contracts',
                                  'no unique generation number, full species table, global form or empirical numbers',
                                  'no complete instrument or unlimited-energy construction claimed'],
                    sha256={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in files},
                    unchanged_closure_manifest_sha256=sha((RESEARCH/'_migration/closure_1288_1307/manifest.json').read_bytes()))
    assert layer_path.read_bytes() == layer_raw
    for path, old, new in updates:
        assert path.read_bytes() == old
    with acceptance.open('x', encoding='utf8', newline='\n') as f:
        json.dump(evidence, f, ensure_ascii=False, indent=2)
        f.write('\n')
    for path, old, new in updates:
        path.write_bytes(new)
    layer_path.write_text(json.dumps(layer, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps(dict(status=evidence['status'], rounds=evidence['scientific_rounds'],
                         artifact_bindings=binds, validation=evidence['validation']), ensure_ascii=False))


if __name__ == '__main__':
    main()
