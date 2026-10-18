"""Publish round 1317 only after frozen independent review; preserve old navigation."""
from pathlib import Path
import ast
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    receipt_path = STAGE / '_shared/progress_publication_after1317.json'
    checks_path = STAGE / '_shared/progress_checks_after1317.json'
    acceptance_path = STAGE / '1317/acceptance.json'
    assert not any(p.exists() for p in (receipt_path, checks_path, acceptance_path))
    InteractionGenerationLayout().verify_interaction_navigation()
    old_layer_raw = LAYER.read_bytes()
    old_layer = json.loads(old_layer_raw)
    assert old_layer['latest_accepted_round'] == 1316
    folder = STAGE / '1317'
    review_path = folder / 'independent_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    package = [folder / n for n in ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                                    'independent_check.py', 'independent_results.json')]
    package += [STAGE / 'research_note_1317.md', STAGE / '_admission/round1317_admission.md']
    for path in package:
        assert sha(path.read_bytes()) in review, str(path)
        if path.suffix == '.py':
            ast.parse(path.read_text(encoding='utf8'))
    outputs = {name: json.loads((folder / name).read_bytes())
               for name in ('results.json', 'independent_results.json')}
    assert all(r['status'].startswith('PASS') for r in outputs.values())
    executions = json.loads((STAGE / '_shared/replay_after1317.json').read_bytes())
    assert executions['status'] == 'PASS' and len(executions['checks']) == 2
    for check in executions['checks']:
        assert check['exit_code'] == 0
        assert sha((ROOT / check['script']).read_bytes()) == check['script_sha256']
    package.append(review_path)
    acceptance = dict(round=1317, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active',
        executed_checks=executions['checks'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
        main_result='Full initial-window refusal error and a registered portal norm bound constrain the heavy response among portals with the same light Y.',
        scope_limits=['Fixed chain and local-only permission control are candidate inputs.',
                      'Finite neutral CAR process shares a matrix/source dictionary, not a certified charged QFT receiver.',
                      'No generation, absolute mass, full A/B/C/J or autonomous withdrawal selection.'])
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    old_readme, old_ledger = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = old_readme.decode('utf8')
    readme = readme.replace('最新验收至1316，共四轮新研究', '最新验收至1317，共五轮新研究')
    readme = readme.replace('(_shared/progress_after1316.md)明确', '(_shared/progress_after1317.md)明确')
    anchor = '|[1316](research_note_1316.md)|同一有限接续核约束轻耦合与重响应的广义谱，给完整质量极点的有限夹界|[独审](1316/independent_math_review.md) · [验收](1316/acceptance.json)|'
    assert readme.count(anchor) == 1
    readme = readme.replace(anchor, anchor + '\n|[1317](research_note_1317.md)|全起始窗拒绝误差与入口资源上界约束同Y重响应；删弱给强屏蔽与总许可反例|[独审](1317/independent_review.md) · [验收](1317/acceptance.json)|')
    first = readme.index('下一科学编号1317')
    last = readme.index('停止一般SWAP', first)
    readme = readme[:first] + '''下一科学编号1318，尚未准入。先审[331候选的共同连接](_admission/family_anomaly_selection_screen_after1316.md)：保全部允许中性质量项，核代际条件、离散保护与关联重流能否在同一作用实现；成熟三代计数本身不新开轮。该候选群、表示与运行合同仍为输入。

1317已关闭本链的有限证书问题，停止入口精度和器件优化。一般共同许可不选择端点；三／四／五代、任意轻Y与尺度自由在本链中仍保留。各轮逐项相容未替代A/B/C/J联合验收，不混用批量合同与单位E，不把新链或新群自动算作旧M_Y的推论。

[后继计0筛查](_admission/joint_mechanism_screen_after1316.md)、[一般许可](_admission/joint_consent_scope_screen_after1316.md)和[拒绝权限](_admission/refusal_chain_screen_after1316.md)保留原时点；其中“尚未启动1317”是历史状态，当前以本导航为准。保持轻Y不等于保持完整质量极点。

''' + readme[last:]
    ledger_add = '''
## 1317后的更新

2026-10-19。原时点记录保留，新增公理0。

有限起始窗拒绝误差与事前总入口上界，在固定局部控制权限下给非末端入口界，再给同Y重响应界；强耦合和总许可删弱分别证明资源与控制范围不可省。非交换正例适用于三／四／五通道。有限中性CAR与场论只共享明确核／源字典，未自动完成原弱带荷仪器或A/B/C/J联合实现。

下一编号1318先审331候选的同父资格，成熟三代计数不重报；群、字段、真空、运行正号与完整全局资格分别登记。停止本链精度和长度优化。[联合进展与输入边界](progress_after1317.md)。
'''
    prefix = '''
## 当前目标：从认知交互机制生成物质分工、代际与质量结构

2026-10-19，[1313起新阶段](archive_1313_/README.md)已验收1313—1317五轮，目标保持active，新公理0。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1317.md)。

[1313](archive_1313_/research_note_1313.md)判原有限交接不独选物质角色；[1314—1315](archive_1313_/_shared/progress_after1316.md)区分单位交换的冲突与批量共同分支；[1316](archive_1313_/research_note_1316.md)建立同一有限核的轻耦合—重响应关系与质量极点界。

[1317](archive_1313_/research_note_1317.md)把全起始窗有限拒绝证书和入口资源上界，接到同Y重响应约束；删弱给强屏蔽与总许可反例。控制范围、链、Y和尺度仍是输入；有限中性CAR核的实现不等于原弱带荷仪器已实现，三／四／五代均未被该机制区分。

下一编号1318尚未准入，先审[331候选的代际—保护—完整质量—重流共同连接](archive_1313_/_admission/family_anomaly_selection_screen_after1316.md)。成熟三代计数本身计0，候选群与运行条件不称认知推论。停止本链优化，继续A/B/C/J联合目标，未唯一生成现实标准模型。

上一阶段已归档为[archive_1308_1312](archive_1308_1312/README.md)，[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)及[历史复算](_migration/closure_1308_1312/README.md)保留。下方旧时点记录不改。
'''
    layer = copy.deepcopy(old_layer)
    updates = []
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        original = relocation_original_bytes(raw, entry)
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        text = original.decode(encoding)
        at = entry['edits'][0]['start']
        new = (text[:at] + prefix + text[at:]).encode(encoding)
        entry['current_sha256'], entry['current_bytes'] = sha(new), len(new)
        entry['edits'][0]['after'] = prefix
        assert relocation_original_bytes(new, entry) == original
        updates.append((path, raw, new))
    layer.update(latest_accepted_round=1317, next_round=1318, in_progress_round=None,
                 stage_complete=False, application_goal_status_last_observed='active')
    layer_raw = encoded(layer)
    updates += [(readme_path, old_readme, readme.encode('utf8')),
                (ledger_path, old_ledger, old_ledger + ledger_add.encode('utf8')),
                (LAYER, old_layer_raw, layer_raw)]
    planned = {p: n for p, _, n in updates}
    planned[acceptance_path] = encoded(acceptance)
    checked_md = [p for p in package if p.suffix == '.md']
    checked_md += [readme_path, ledger_path, STAGE / '_shared/progress_after1317.md']
    checked_md += [STAGE / '_admission' / n for n in (
        'joint_consent_scope_screen_after1316.md', 'refusal_chain_screen_after1316.md',
        'joint_mechanism_screen_after1316.md', 'family_anomaly_selection_screen_after1316.md')]
    links = 0
    for path in checked_md:
        text = planned[path].decode('utf8') if path in planned else path.read_text(encoding='utf8')
        assert text.count('$$') % 2 == 0, str(path)
        assert not re.search(r'^\s*\\[\[\]]\s*$', text, re.M), str(path)
        for target in re.findall(r'\]\(([^)\s]+)\)', text):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (path.parent / unquote(target.split('#')[0])).resolve()
            assert resolved.exists() or resolved in planned, (path, target)
            links += 1
    checks = dict(status='PASS', date='2026-10-19', markdown_files=len(checked_md),
        local_links=links, reviewed_rounds=[1317], new_axioms=0,
        old_navigation_boundary_preserved=True,
        frozen_support=len(layer['unchanged_support_sha256']),
        scope='New round, count-zero screens and reversible navigation; not a repository-wide audit.')
    receipt = dict(status='PUBLISHED', latest_accepted_round=1317, new_scientific_rounds=[1317],
        new_axioms=0, stage_complete=False, goal_status='active',
        previous_layer_sha256=sha(old_layer_raw), current_layer_sha256=sha(layer_raw),
        previous_layer_fields={k:v for k,v in old_layer.items() if k not in ('entries','unchanged_support_sha256')},
        previous_navigation_entries=old_layer['entries'],
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): b.decode('utf8')
            for p,b in ((readme_path,old_readme),(ledger_path,old_ledger))},
        acceptance_sha256=sha(encoded(acceptance)),
        current_screen_sha256={p.relative_to(ROOT).as_posix():sha(p.read_bytes())
                              for p in checked_md if p.parent.name == '_admission'},
        original_1308_closure_navigation_unchanged=True)
    fresh = [(acceptance_path, encoded(acceptance)), (checks_path, encoded(checks)),
             (receipt_path, encoded(receipt))]
    for path, raw, _ in updates:
        assert path.read_bytes() == raw, 'Concurrent modification: ' + str(path)
    assert all(not p.exists() for p,_ in fresh)
    written, created = [], []
    try:
        for path, raw, new in updates:
            atomic(path, new)
            written.append((path,raw,new))
        for path, raw in fresh:
            with path.open('xb') as stream:
                stream.write(raw)
            created.append((path,raw))
    except BaseException:
        for path, raw in reversed(created):
            if path.read_bytes() == raw:
                path.unlink()
        for path, old, new in reversed(written):
            if path.read_bytes() == new:
                atomic(path,old)
        raise
    print(json.dumps(checks,ensure_ascii=False))


if __name__ == '__main__':
    main()
