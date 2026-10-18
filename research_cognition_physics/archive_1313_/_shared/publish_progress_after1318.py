"""Publish the reviewed 1318 package; preserve the prior navigation boundary."""
import ast
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    receipt_path = STAGE / '_shared/progress_publication_after1318.json'
    checks_path = STAGE / '_shared/progress_checks_after1318.json'
    acceptance_path = STAGE / '1318/acceptance.json'
    assert not any(p.exists() for p in (receipt_path, checks_path, acceptance_path))
    InteractionGenerationLayout().verify_interaction_navigation()
    old_layer_raw = LAYER.read_bytes()
    old_layer = json.loads(old_layer_raw)
    assert old_layer['latest_accepted_round'] == 1317
    folder = STAGE / '1318'
    review_path = folder / 'independent_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    package = [folder / n for n in ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                                    'independent_check.py', 'independent_results.json')]
    package += [STAGE / 'research_note_1318.md', STAGE / '_admission/round1318_admission.md']
    for path in package:
        assert sha(path.read_bytes()) in review, str(path)
        if path.suffix == '.py':
            ast.parse(path.read_text(encoding='utf8'))
    assert all(json.loads((folder / n).read_bytes())['status'].startswith('PASS')
               for n in ('results.json', 'independent_results.json'))
    executions = json.loads((STAGE / '_shared/replay_after1318.json').read_bytes())
    assert executions['status'] == 'PASS' and len(executions['checks']) == 2
    for check in executions['checks']:
        assert check['exit_code'] == 0
        assert sha((ROOT / check['script']).read_bytes()) == check['script_sha256']
    package.append(review_path)
    acceptance = dict(round=1318, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active',
        executed_checks=executions['checks'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
        main_result='A declared 331 action jointly supports the R-protected full neutral Dirac block, stable local vacuum and correlated ordinary-quark transverse neutral response, including its complete Cartan-sector sum rule and sources.',
        scope_limits=['The generation count reuses known 331 anomaly plus an extra full-threshold color-running contract.',
                      'Groups, representations, Yukawa data, vacuum scales and the running contract remain inputs.',
                      'The transverse response is not full scattering; non-Cartan neutral modes are explicitly accounted for in the relevant block test.',
                      'No single-Higgs SM, complete global-background certificate, all instruments or joint cognitive A/B/C/J completion.'])
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    old_readme, old_ledger = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = old_readme.decode('utf8')
    assert '最新验收至1317，共五轮新研究' in readme
    readme = readme.replace('最新验收至1317，共五轮新研究', '最新验收至1318，共六轮新研究')
    readme = readme.replace('(_shared/progress_after1317.md)明确', '(_shared/progress_after1318.md)明确')
    anchor = '|[1317](research_note_1317.md)|全起始窗拒绝误差与入口资源上界约束同Y重响应；删弱给强屏蔽与总许可反例|[独审](1317/independent_review.md) · [验收](1317/acceptance.json)|'
    assert readme.count(anchor) == 1
    readme = readme.replace(anchor, anchor + '\n|[1318](research_note_1318.md)|同一331作用保允许交叉、稳定真空与Dirac谱；普通夸克的两部门投影及中性横向响应满足共同和式|[独审](1318/independent_review.md) · [验收](1318/acceptance.json)|')
    first = readme.index('下一科学编号1318')
    last = readme.index('历史启动验收', first)
    readme = readme[:first] + '''下一科学编号1319，先做实际角色来源的机制准入，不自动展开新的物理候选。[认知角色查重](_admission/cognitive_family_roles_screen_after1317.md)表明一般身份／伙伴／记录合同对同型复制不敏感；402已实现秘密共享和恢复，尚缺其与物理家族的实际映射，不重建纠错器件。

1318的三代条件限公开331表示及额外运行合同，不排除其它父作用中的四、五代。群、角色、Yukawa、尺度及运行正号仍是输入；多Higgs和重字段保留，未唯一生成现实标准模型。旧单位E与批量合同不混用，逐项相容不替代A/B/C/J联合验收。

[331事前筛查](_admission/family_anomaly_selection_screen_after1316.md)及[1318准入](_admission/round1318_admission.md)保留原时点；正式范围以验收与证明为准。停止该候选的参数扫描、单链精度优化、一般SWAP和重复秘密共享。新增公理0，目标继续。

''' + readme[last:]
    ledger_add = '''
## 1318后的更新

2026-10-19。新公理0，阶段目标继续。

公开331类别的反常与额外全阈值色运行合同给三代；这一成熟计数不重复认领为新发现。同一作用的残余保护允许完整2n Dirac块，给稳定局部真空、跨夸克部门投影关系及Cartan中性模合计交换权和式。完整规范质量核核对了其它中性场在指定源／阶内的角色，参数变分保留质量、流和混合。

共同增益不能消除群／表示、运行合同、Yukawa、真空尺度等输入；也未签单Higgs、所有全局背景或全部仪器。下一项转向实际角色为何不能用同型副本替代；402已覆盖秘密共享构造，未给场—家族映射前不重报。[联合进展](progress_after1318.md)。
'''
    prefix = '''
## 当前目标：从认知交互机制生成物质分工、代际与质量结构

2026-10-19，[1313起新阶段](archive_1313_/README.md)已验收1313—1318六轮，目标保持active，新公理0。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1318.md)。

[1313—1317](archive_1313_/_shared/progress_after1317.md)已区分交接任务的无选择力、单位与批量交换，以及有限操作证书对同Y响应的约束；未由此独选三代或现实参数。

[1318](archive_1313_/research_note_1318.md)在同一公开331候选中，保全部允许中性质量交叉，给稳定局部真空、完整Dirac谱、跨夸克部门的流关系及中性横向交换和式。三代计数复用成熟反常机制并依赖额外运行合同；群、表示、Yukawa和尺度仍输入，多Higgs及重字段保留，未称整个标准模型已唯一生成。

下一编号1319先审[认知角色来源](archive_1313_/_admission/cognitive_family_roles_screen_after1317.md)。一般保身份／伙伴／记录对同型复制不敏感；402已完成秘密共享与恢复，不重建器件。先找实际场—接口映射与有选择力的后续任务，停止修补单个候选和参数拟合，继续A/B/C/J联合目标。

上一阶段归档为[archive_1308_1312](archive_1308_1312/README.md)，[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)及[历史复算](_migration/closure_1308_1312/README.md)保留。下方旧时点记录不改。
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
    layer.update(latest_accepted_round=1318, next_round=1319, in_progress_round=None,
                 stage_complete=False, application_goal_status_last_observed='active')
    layer_raw = encoded(layer)
    updates += [(readme_path, old_readme, readme.encode('utf8')),
                (ledger_path, old_ledger, old_ledger + ledger_add.encode('utf8')),
                (LAYER, old_layer_raw, layer_raw)]
    planned = {p:n for p,_,n in updates}
    planned[acceptance_path] = encoded(acceptance)
    checked_md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
        STAGE / '_shared/progress_after1318.md',
        STAGE / '_admission/cognitive_family_roles_screen_after1317.md']
    links = 0
    for path in checked_md:
        text = planned[path].decode('utf8') if path in planned else path.read_text(encoding='utf8')
        assert text.count('$$') % 2 == 0, str(path)
        assert not re.search(r'^\s*\\[\[\]]\s*$',text,re.M),str(path)
        for target in re.findall(r'\]\(([^)\s]+)\)',text):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (path.parent / unquote(target.split('#')[0])).resolve()
            assert resolved.exists() or resolved in planned,(path,target)
            links += 1
    checks = dict(status='PASS',date='2026-10-19',markdown_files=len(checked_md),
        local_links=links,reviewed_rounds=[1318],new_axioms=0,
        old_navigation_boundary_preserved=True,frozen_support=len(layer['unchanged_support_sha256']),
        scope='New round and direction screen, stage navigation and reversible root insertion; not all repository history.')
    receipt = dict(status='PUBLISHED',latest_accepted_round=1318,new_scientific_rounds=[1318],
        new_axioms=0,stage_complete=False,goal_status='active',
        previous_layer_sha256=sha(old_layer_raw),current_layer_sha256=sha(layer_raw),
        previous_layer_fields={k:v for k,v in old_layer.items() if k not in ('entries','unchanged_support_sha256')},
        previous_navigation_entries=old_layer['entries'],
        previous_stage_navigation={p.relative_to(ROOT).as_posix():b.decode('utf8')
            for p,b in ((readme_path,old_readme),(ledger_path,old_ledger))},
        acceptance_sha256=sha(encoded(acceptance)),
        cognitive_screen_sha256=sha((STAGE/'_admission/cognitive_family_roles_screen_after1317.md').read_bytes()),
        original_1308_closure_navigation_unchanged=True)
    fresh = [(acceptance_path,encoded(acceptance)),(checks_path,encoded(checks)),(receipt_path,encoded(receipt))]
    for path,raw,_ in updates:
        assert path.read_bytes()==raw,'Concurrent modification: '+str(path)
    assert all(not p.exists() for p,_ in fresh)
    written,created=[],[]
    try:
        for path,raw,new in updates:
            atomic(path,new); written.append((path,raw,new))
        for path,raw in fresh:
            with path.open('xb') as stream: stream.write(raw)
            created.append((path,raw))
    except BaseException:
        for path,raw in reversed(created):
            if path.read_bytes()==raw: path.unlink()
        for path,old,new in reversed(written):
            if path.read_bytes()==new: atomic(path,old)
        raise
    print(json.dumps(checks,ensure_ascii=False))


if __name__ == '__main__':
    main()
