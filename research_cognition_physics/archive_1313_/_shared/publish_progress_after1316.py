"""Publish reviewed rounds once, retaining reversible historical navigation."""
from pathlib import Path
import ast
import copy
import hashlib
import json
import os
import re
import sys
from urllib.parse import unquote

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout, LAYER


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf8')


def atomic(path, raw):
    temp = path.with_name(path.name + '.publish1316.tmp')
    with temp.open('xb') as stream:
        stream.write(raw)
    try:
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def main():
    receipt_path = STAGE / '_shared/progress_publication_after1316.json'
    checks_path = STAGE / '_shared/progress_checks_after1316.json'
    assert not receipt_path.exists() and not checks_path.exists()
    store = InteractionGenerationLayout()
    old_layer_raw = LAYER.read_bytes()
    old_layer = json.loads(old_layer_raw)
    assert old_layer['latest_accepted_round'] == 1314
    for name, expected in old_layer['unchanged_support_sha256'].items():
        assert sha((ROOT / name).read_bytes()) == expected, name
    acceptance = {}
    checked_md = []
    check_details = {
        1315: [{'script': 'check.py', 'status': 'PASS', 'groups': 5},
               {'script': 'independent_check.py', 'status': 'PASS_WITH_SCOPE',
                'target_cases': 8, 'massive_branches': 3}],
        1316: [{'script': 'check.py', 'status': 'PASS', 'cases': 18,
                'source_parameters': 5, 'multiport_deletions': 9},
               {'script': 'independent_check.py', 'status': 'PASS_WITH_SCOPE',
                'sample_count': 20, 'assertion_count': 1055, 'check_categories': 36}]}
    for number, review_name in ((1315, 'independent_review.md'),
                                 (1316, 'independent_math_review.md')):
        folder = STAGE / str(number)
        target = folder / 'acceptance.json'
        assert not target.exists()
        review = (folder / review_name).read_text(encoding='utf8')
        assert 'PASS_WITH_SCOPE' in review
        package = [folder / name for name in ('proof.md', 'check.py', 'results.json',
                    'NEXT.md', 'independent_check.py', 'independent_results.json')]
        package += [STAGE / f'research_note_{number}.md',
                    STAGE / f'_admission/round{number}_admission.md']
        for path in package:
            assert sha(path.read_bytes()) in review, str(path)
            if path.suffix == '.py':
                ast.parse(path.read_text(encoding='utf8'))
        for name in ('results.json', 'independent_results.json'):
            result = json.loads((folder / name).read_bytes())
            assert result['status'].startswith('PASS')
        package.append(folder / review_name)
        checked_md += [path for path in package if path.suffix == '.md']
        acceptance[number] = dict(round=number, status='ACCEPTED_CONDITIONAL_SCOPE',
            date='2026-10-18', new_axioms=0, stage_complete=False, goal_status='active',
            checks='Author and independent programs were run read-only before publication; this publisher verifies their frozen outputs and review bindings.',
            executed_checks=[dict(exit_code=0, **item) for item in check_details[number]],
            sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
            main_result=('Flavor-covariant nonzero batch exchange saturates D=6n in the declared equal B/L exchange class; massive leading joint branch.'
                         if number == 1315 else
                         'Fixed finite Dirac chains give portal-independent light/heavy generalized spectrum and exact finite-window physical pole bounds.'),
            scope_limits=(['unit E not revoked; Nc=3 input', 'not a generation selector or rate calculation',
                           'leading EFT common action, not full instrument certification'] if number == 1315 else
                          ['tree-level canonical finite-chain category', 'Y, n, scale and endpoint rules remain inputs',
                           'not an adopted cognition axiom or a full common A/B/C implementation']))

    readme_path, ledger_path = STAGE / 'README.md', STAGE / '_shared/selection_ledger.md'
    old_readme, old_ledger = readme_path.read_bytes(), ledger_path.read_bytes()
    header = old_readme.decode('utf8').split('## 当前状态')[0]
    readme = header + '''## 当前状态

最新验收至1316，共四轮新研究，作者与独立复算均通过。阶段目标保持active，新公理0。[联合进展及下一步](_shared/progress_after1316.md)明确各项机制的范围与剩余输入。

|轮次|新增判定|证据|
|---|---|---|
|[1313](research_note_1313.md)|三种实际参考交接给相同完整逻辑仪器，该有限任务不独选物质分工|[独审](1313/independent_math_review.md) · [验收](1313/acceptance.json)|
|[1314](research_note_1314.md)|原Dirac资料单值协变生成与单位交换E不相容|[独审](1314/independent_review.md) · [验收](1314/acceptance.json)|
|[1315](research_note_1315.md)|开放批量交换有全味非零Pfaffian与正质量共同领先分支，指定同向交换类D最小为6n|[独审](1315/independent_review.md) · [验收](1315/acceptance.json)|
|[1316](research_note_1316.md)|同一有限接续核约束轻耦合与重响应的广义谱，给完整质量极点的有限夹界|[独审](1316/independent_math_review.md) · [验收](1316/acceptance.json)|

下一科学编号1317，先筛查接入权限与交换职责的共同来源。三／四／五代、任意轻Y与尺度自由仍保留；不把批量合同与原单位E混用，也不将新链资料自动算作原M_Y的推论。各轮逐项相容未替代A/B/C/J联合验收。

停止一般SWAP、同一预算优化、Pfaffian重算、一般Schur正性和无判据的链长拟合。新命题必须能区别真实竞争方案；不自动采用新公理，不默认最少物种或最少代际。目标设定、归档、筛查及综合均计0。

历史启动验收`_shared/start_checks.json`只证明启动快照；1314时点见[旧联合进展](_shared/progress_after1314.md)。当前以各轮验收和本导航为准，保留所有历史证据。
'''
    ledger_append = '''
## 1315—1316后的更新

2026-10-18。原基线及1314时点记录保留；新增公理0。

|主线|新结果|仍需区分|
|---|---|---|
|B 交换合同|1315对M_Y与较弱批量任务给完整非零饱和构造及正质量领先共同分支|单位E仍禁，三色输入，三／四／五代均保留；非零Wilson系数未生成|
|C 共同参数|1316固定链核把轻Y与重响应广义谱联系，物理极点有有限夹界|任意Y、链数、m及末端规则仍输入，多端口保持Y可改变重响应|
|A/J 共同机制|各轮合同和父来源分别明确；新交换与中性链逐项保护相容|不等链长未自动满足原M_Y；1313未选择单端点；尚无同一A/B/C完整机制|

下一步先审真实职责是否约束接入口和交换，而不是继续调候选参数。[联合范围与停止线](progress_after1316.md)。
'''
    new_readme, new_ledger = readme.encode('utf8'), old_ledger + ledger_append.encode('utf8')
    prefix = '''
## 当前目标：从认知交互机制生成物质分工、代际与质量结构

2026-10-18，[1313起新阶段](archive_1313_/README.md)已验收1313—1316四轮，作者与独立复算通过，目标保持active。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1316.md)。新公理0。

[1313](archive_1313_/research_note_1313.md)表明原有限逻辑交接不独选物质角色；[1314](archive_1313_/research_note_1314.md)判原Dirac资料单值协变生成与单位E冲突。[1315](archive_1313_/research_note_1315.md)开放批量合同，给全味非零Pfaffian和正质量共同领先分支；不撤销单位E，也不继承删除前件后的旧颜色选择结论。

[1316](archive_1313_/research_note_1316.md)在固定有限接续核中，消去入口矩阵，建立轻质量耦合—重响应的广义谱关系及完整质量极点的有限夹界。任意轻Y、链数、尺度及接入规则仍是输入；三／四／五代均保留。未由这些条件结果唯一生成现实标准模型。

下一编号1317，先审实际交接与交换职责是否能约束接入权限及共同更新。停止同任务预算优化、Pfaffian重算和链长拟合；不把新链资料自动当作旧生成律的推论，不把逐项相容当成A/B/C/J联合完成。

上一阶段归档为[archive_1308_1312](archive_1308_1312/README.md)，[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)及[历史复算](_migration/closure_1308_1312/README.md)保留。下方旧时点记录不改，当前状态以本节及阶段验收为准。
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
    layer.update(latest_accepted_round=1316, next_round=1317, in_progress_round=None,
                 stage_complete=False, application_goal_status_last_observed='active')
    layer_raw = encoded(layer)
    updates += [(readme_path, old_readme, new_readme), (ledger_path, old_ledger, new_ledger),
                (LAYER, old_layer_raw, layer_raw)]
    planned = {path: new for path, _, new in updates}
    for number, record in acceptance.items():
        planned[STAGE / str(number) / 'acceptance.json'] = encoded(record)
    checked_md += [readme_path, ledger_path, STAGE / '_shared/progress_after1316.md']
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
    checks = dict(status='PASS', date='2026-10-18', markdown_files=len(checked_md),
                  local_links=links, reviewed_rounds=[1315, 1316], new_axioms=0,
                  old_navigation_boundary_preserved=True, frozen_support=len(layer['unchanged_support_sha256']),
                  scope='New round packages, stage navigation and reversible root insertion; not a full repository audit.')
    receipt = dict(status='PUBLISHED', latest_accepted_round=1316,
                   new_scientific_rounds=[1315, 1316], new_axioms=0, stage_complete=False,
                   goal_status='active', previous_layer_sha256=sha(old_layer_raw),
                   current_layer_sha256=sha(layer_raw),
                   previous_layer_fields={k:v for k,v in old_layer.items()
                                          if k not in ('entries', 'unchanged_support_sha256')},
                   previous_navigation_entries=old_layer['entries'],
                   previous_stage_navigation={p.relative_to(ROOT).as_posix(): b.decode('utf8')
                                              for p,b in ((readme_path, old_readme), (ledger_path, old_ledger))},
                   acceptance_sha256={str(n): sha(encoded(r)) for n,r in acceptance.items()},
                   original_1308_closure_navigation_unchanged=True)
    fresh = [(STAGE / str(n) / 'acceptance.json', encoded(r)) for n,r in acceptance.items()]
    fresh += [(checks_path, encoded(checks)), (receipt_path, encoded(receipt))]
    for path, raw, _ in updates:
        assert path.read_bytes() == raw, 'Concurrent change: ' + str(path)
    assert all(not p.exists() for p,_ in fresh)
    written, created = [], []
    try:
        for path, raw, new in updates:
            assert path.read_bytes() == raw
            atomic(path, new)
            written.append((path, raw, new))
        for path, raw in fresh:
            with path.open('xb') as stream:
                stream.write(raw)
            created.append((path, raw))
    except BaseException:
        for path, raw in reversed(created):
            if path.read_bytes() == raw:
                path.unlink()
        for path, old, new in reversed(written):
            if path.read_bytes() == new:
                atomic(path, old)
        raise
    print(json.dumps(checks, ensure_ascii=False))


if __name__ == '__main__':
    main()
