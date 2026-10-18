"""One-time per-round acceptance and reversible root-navigation update."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import relocation_original_bytes


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main(n):
    target = STAGE/str(n)
    assert not (target/'acceptance.json').exists(), 'Already accepted; do not repeat'
    reviews = ['independent_review.md']
    if n == 1312:
        reviews += ['independent_scope_review.md']
    for name in reviews:
        assert 'PASS' in (target/name).read_text(encoding='utf8')
    receipts = []
    for script in ('check.py', 'independent_check.py'):
        p = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(target/script), '--check'],
                           cwd=ROOT, capture_output=True, text=True, encoding='utf8')
        assert p.returncode == 0, (script, p.stdout, p.stderr)
        receipts.append(dict(script=script, exit_code=p.returncode, output=p.stdout.strip()))
    names = ['proof.md', 'check.py', 'results.json', 'independent_check.py',
             'independent_results.json', 'NEXT.md'] + reviews
    files = [target/name for name in names]
    files += [STAGE/f'research_note_{n}.md', STAGE/f'_admission/round{n}_admission.md']
    evidence = dict(round=n, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-18',
                    new_axioms=0, stage_complete=False, executed_checks=receipts,
                    main_result=('fixed true-task parameter domain and exact three/four/five-generation degeneracy'
                                 if n == 1311 else 'single-reference full Dirac mass selects faithful SU(2) activity in the declared class'),
                    sha256={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in files},
                    depends_on={str(k):sha((STAGE/str(k)/'acceptance.json').read_bytes())
                                for k in range(1308, n)},
                    scope_limits=['declared roles and finite effective class only',
                                  'no new axioms and no unique empirical Standard Model claim',
                                  'actual complete instruments and higher-order error certificates not added'])
    path = ROOT/'research_cognition_physics/_migration/conditional_selection_active_1308/layer.json'
    raw = path.read_bytes()
    layer = json.loads(raw)
    assert layer['latest_accepted_round'] == n-1
    tail = ('[1312](archive_1308_/research_note_1312.md)已准入，检验单参考完整Dirac质量能否选择活动群；'
            '尚未验收，不将候选结论算作已证明。' if n == 1311 else
            '[1312](archive_1308_/research_note_1312.md)进一步在公开的单参考质量类中选择SU(2)标准活动像。'
            '右单态、球面传递、完整满秩等仍是额外条件；颜色、全局商和现实数值不随之生成。')
    prefix = ('\n## 当前阶段：群、代际与参数关系的条件性生成\n\n'
              f'2026-10-18，最新已验收[{n}](archive_1308_/research_note_{n}.md)。'
              '[阶段目录](archive_1308_/README.md) · [正式目标](archive_1308_/_shared/goal.md)。新公理0，目标继续进行。\n\n'
              '1308—1310已连接单位交换、保护、相位任务、实际有限交付与运行预算。'
              '[1311](archive_1308_/research_note_1311.md)给固定真实核的完整耦合—Higgs尺度域，'
              '并证明三色三／四／五代在这套合同中的整个域相同；继续调整同一预算不能独选三代。'
              + tail + '\n\n'
              '[共同条件整合稿](archive_1308_/_shared/joint_selection_synthesis.md)与'
              '[阶段论文审阅稿](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)正在核验。'
              '阶段联合验收尚未完成，不把局部结果当成整个目标结项。\n\n'
              '旧阶段已归档为archive_1288_1289及archive_1290_1307；'
              '[迁移与历史复算](_migration/closure_1288_1307/README.md)保持。'
              '下方为历史时点记录，以本节及正式目标为准。\n')
    updates = []
    for entry in layer['entries']:
        nav = ROOT/entry['destination']
        current = nav.read_bytes()
        original = relocation_original_bytes(current, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = original.decode(enc)
        assert len(entry['edits']) == 1 and entry['edits'][0]['before'] == ''
        at = entry['edits'][0]['start']
        changed = (text[:at]+prefix+text[at:]).encode(enc)
        entry['edits'][0]['after'] = prefix
        entry['current_sha256'], entry['current_bytes'] = sha(changed), len(changed)
        assert relocation_original_bytes(changed, entry) == original
        updates.append((nav, current, changed))
    layer.update(latest_accepted_round=n, in_progress_round=1312 if n == 1311 else None,
                 current_work='round1312 mass-selection study' if n == 1311 else 'joint synthesis and stage review',
                 application_goal_status_last_observed='active')
    assert path.read_bytes() == raw
    for nav, old, changed in updates:
        assert nav.read_bytes() == old, 'Concurrent navigation edit'
    with (target/'acceptance.json').open('x', encoding='utf8', newline='\n') as f:
        json.dump(evidence, f, ensure_ascii=False, indent=2)
        f.write('\n')
    for nav, old, changed in updates:
        nav.write_bytes(changed)
    path.write_text(json.dumps(layer, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps(dict(accepted_round=n, stage_complete=False, checks=receipts), ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--round', required=True, type=int, choices=[1311, 1312])
    main(p.parse_args().round)
