"""One-time acceptance and reversible navigation update after manual scope review."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    target = STAGE / '1310'
    assert not (target / 'acceptance.json').exists(), 'Already accepted; do not repeat'
    review = (target / 'independent_review.md').read_text(encoding='utf8')
    assert 'PASS' in review, 'Requires completed independent scope review'
    receipts = []
    for script in ('check.py', 'independent_check.py'):
        run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(target / script), '--check'],
                             cwd=ROOT, capture_output=True, text=True, encoding='utf8')
        assert run.returncode == 0, (script, run.stdout, run.stderr)
        receipts.append(dict(script=script, exit_code=run.returncode, output=run.stdout.strip()))
    files = [target / name for name in ('proof.md', 'check.py', 'results.json',
             'independent_check.py', 'independent_results.json', 'independent_review.md', 'NEXT.md')]
    files += [STAGE / 'research_note_1310.md', STAGE / '_admission/round1310_admission.md']
    evidence = dict(round=1310, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-18',
        new_axioms=0, stage_complete=False,
        main_result='finite unpostselected scattering task connects weak coupling and Higgs scale to the finite color-generation budget',
        executed_checks=receipts,
        scope_limits=['declared massive transverse channel and calibrated effective counting model',
                      'fixed true kernel threshold versus upper-kernel necessary bound distinguished',
                      'matching, higher corrections and complete instruments not certified',
                      'no exact three-generation or unique empirical Standard Model claim'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in files},
        depends_on={str(n): sha((STAGE / str(n) / 'acceptance.json').read_bytes()) for n in (1308, 1309)})

    layer_path = ROOT / 'research_cognition_physics/_migration/conditional_selection_active_1308/layer.json'
    layer_raw = layer_path.read_bytes()
    layer = json.loads(layer_raw)
    assert layer['latest_accepted_round'] == 1309
    prefix = ('\n## 当前阶段：群、代际与参数关系的条件性生成\n\n'
        '2026-10-18，目标继续进行，最新已验收[1310](archive_1308_/research_note_1310.md)，'
        '接[1308](archive_1308_/research_note_1308.md)与[1309](archive_1308_/research_note_1309.md)。'
        '[阶段目录](archive_1308_/README.md)与[正式目标](archive_1308_/_shared/goal.md)。新公理0。\n\n'
        '1308—1309在公开的颜色表示族中给出有限预算的颜色／代际分类。1310补上'
        '有限曝光、正交付裕量和完整失败记录到实际带电流截面、同一弱耦合及Higgs尺度的连接，'
        '再接回整数预算。正宽、有质量横向通道与满秩Dirac相位在领先有效层共同非空。'
        '固定真实核的任务阈值与仅用核上界所得的必要外界分别登记。\n\n'
        '群类型、字段／荷表、任务预算和真实误差仍有输入；尚未选出恰三代、现实参数或完整实际仪器。'
        '[当前输入账](archive_1308_/_shared/selection_ledger_after1310.md)。'
        '停止该器件继续细化，按[下一步](archive_1308_/1310/NEXT.md)转向上游共同接收任务。'
        '[仅Abelian补偿选群的筛查](archive_1308_/_admission/group_scope_screen_after1309.md)'
        '因缺少独立任务依据暂计0；1311尚未准入，不新增公理。\n\n'
        '两旧阶段已收尾为archive_1288_1289、archive_1290_1307；'
        '[迁移与历史复算](_migration/closure_1288_1307/README.md)保持。'
        '下方为历史时点记录，当前方向以本节和正式目标为准。\n')
    updates = []
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        current = path.read_bytes()
        original = relocation_original_bytes(current, entry)
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        prior = original.decode(encoding)
        assert len(entry['edits']) == 1 and entry['edits'][0]['before'] == ''
        at = entry['edits'][0]['start']
        new = (prior[:at] + prefix + prior[at:]).encode(encoding)
        entry['edits'][0]['after'] = prefix
        entry['current_sha256'], entry['current_bytes'] = sha(new), len(new)
        assert relocation_original_bytes(new, entry) == original
        updates.append((path, current, new))
    layer.update(latest_accepted_round=1310, in_progress_round=None,
        current_work='upstream joint-reception task screen; round1311 not admitted',
        application_goal_status_last_observed='active')
    assert layer_path.read_bytes() == layer_raw
    for path, old, new in updates:
        assert path.read_bytes() == old, 'Concurrent navigation edit'
    with (target / 'acceptance.json').open('x', encoding='utf8', newline='\n') as handle:
        json.dump(evidence, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    for path, old, new in updates:
        path.write_bytes(new)
    layer_path.write_text(json.dumps(layer, ensure_ascii=False, indent=2) + '\n', encoding='utf8', newline='\n')
    print(json.dumps(dict(accepted_round=1310, checks=receipts,
                         reversible_navigation_files=len(updates), stage_complete=False), ensure_ascii=False))


if __name__ == '__main__':
    main()
