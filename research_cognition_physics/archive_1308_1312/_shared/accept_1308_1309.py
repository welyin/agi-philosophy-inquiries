"""Accept two independently reviewed rounds and update only the reversible live navigation prefix."""
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
    assert all(not (STAGE / str(n) / 'acceptance.json').exists() for n in (1308, 1309))
    evidence = {}
    for n in (1308, 1309):
        tasks = [('check.py', ['--check']), ('independent_check.py', ['--check'] if n == 1308 else [])]
        evidence[n] = []
        for script, args in tasks:
            command = [sys.executable, '-B', '-X', 'utf8', str(STAGE / str(n) / script)] + args
            run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, encoding='utf-8')
            assert run.returncode == 0, (n, script, run.stdout, run.stderr)
            evidence[n].append({'script': script, 'exit_code': run.returncode, 'output': run.stdout.strip()})

    for n in (1308, 1309):
        path = STAGE / str(n)
        files = [path / name for name in ('proof.md', 'check.py', 'results.json',
                                          'independent_check.py', 'independent_results.json',
                                          'independent_review.md', 'NEXT.md')]
        files += [STAGE / f'research_note_{n}.md', STAGE / '_admission' / f'round{n}_admission.md']
        if n == 1308:
            files += [path / 'scope_review.md', path / 'preacceptance_review_corrections.json']
        assert 'PASS' in (path / 'independent_review.md').read_text(encoding='utf-8')
        result = {
            'round': n, 'status': 'ACCEPTED_CONDITIONAL_SCOPE', 'date': '2026-10-18',
            'new_axioms': 0, 'stage_complete': False,
            'main_result': ('Nc=3, ng in {3,4,5} under the declared E/R/P/W structural and leading-order class'
                            if n == 1308 else
                            'exact finite-budget candidate domain and robust envelope; no unique three-generation result'),
            'executed_checks': evidence[n],
            'scope_limits': ['fixed group family and charge menu', 'no full interacting instrument certificate',
                             'no unique empirical Standard Model', 'no new cognitive axiom adoption'],
            'sha256': {p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in files},
        }
        if n == 1309:
            result['depends_on_round1308_acceptance'] = sha((STAGE / '1308/acceptance.json').read_bytes())
        with (path / 'acceptance.json').open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')

    layer_path = ROOT / 'research_cognition_physics/_migration/conditional_selection_active_1308/layer.json'
    original_layer = layer_path.read_bytes()
    layer = json.loads(original_layer)
    prefix = ('\n## 当前阶段：群、代际与参数关系的条件性生成\n\n'
              '2026-10-18，目标继续进行，已验收[1308](archive_1308_/research_note_1308.md)与'
              '[1309](archive_1308_/research_note_1309.md)。新阶段见[archive_1308_](archive_1308_/README.md)，'
              '[正式目标](archive_1308_/_shared/goal.md)保持不变。新公理0。\n\n'
              '1308在公开的开放颜色表示族内，将受保护单位交换、限定相位任务与有限窗口弱运行联立，'
              '条件选出三色、三至五代；没有输入普遍维六上限。1309进一步用有限初端强度／全窗上限及误差预算替换单调要求，'
              '给完整整数可行域。取消统一初端强度下界后，无限候选重新合格。\n\n'
              '这些结论仍固定群类型、荷表和部分任务，只签同作用与领先有效阶的共同结构；'
              '尚未选恰三代、现实参数或全交互仪器。'
              '[选择与剩余输入账](archive_1308_/_shared/selection_ledger_after1309.md)已更新。'
              '下一项先查[强度预算的操作来源](archive_1308_/_admission/finite_strength_origin_screen.md)，'
              '不为选三代调整阈值，不把旧响应界改写算新轮。1310尚未准入。\n\n'
              '两旧阶段已收尾为archive_1288_1289、archive_1290_1307；'
              '[迁移与历史复算](_migration/closure_1288_1307/README.md)保持。'
              '下方为历史时点记录，当前方向以本节和正式目标为准。\n')
    updates = []
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        current = path.read_bytes()
        original = relocation_original_bytes(current, entry)
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        old_text = original.decode(encoding)
        assert len(entry['edits']) == 1 and entry['edits'][0]['before'] == ''
        at = entry['edits'][0]['start']
        new = (old_text[:at] + prefix + old_text[at:]).encode(encoding)
        entry['edits'][0]['after'] = prefix
        entry['current_sha256'], entry['current_bytes'] = sha(new), len(new)
        assert relocation_original_bytes(new, entry) == original
        updates.append((path, current, new))
    layer['latest_accepted_round'] = 1309
    layer['in_progress_round'] = None
    layer['current_work'] = 'finite-strength operational-origin screening; no round1310 admission yet'
    layer['application_goal_status_last_observed'] = 'active'
    assert layer_path.read_bytes() == original_layer
    for path, old, new in updates:
        assert path.read_bytes() == old
    for path, old, new in updates:
        path.write_bytes(new)
    layer_path.write_text(json.dumps(layer, ensure_ascii=False, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({'accepted_rounds': [1308, 1309], 'checks': evidence,
                      'reversible_navigation_files': len(updates), 'stage_complete': False}, ensure_ascii=False))


if __name__ == '__main__':
    main()
