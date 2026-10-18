"""Publish one reviewed round without changing frozen navigation originals."""
import ast
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    paths = {name: STAGE / '_shared' / (name + '_after1319.json')
             for name in ('progress_publication', 'progress_checks')}
    acceptance_path = STAGE / '1319/acceptance.json'
    assert not any(p.exists() for p in [*paths.values(), acceptance_path])
    InteractionGenerationLayout().verify_interaction_navigation()
    old_raw = LAYER.read_bytes()
    old = json.loads(old_raw)
    assert old['latest_accepted_round'] == 1318
    folder = STAGE / '1319'
    review_path = folder / 'independent_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    package = [folder / n for n in ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                                    'independent_check.py', 'independent_results.json')]
    package += [STAGE / 'research_note_1319.md', STAGE / '_admission/round1319_admission.md']
    for path in package:
        assert sha(path.read_bytes()) in review, str(path)
        if path.suffix == '.py':
            ast.parse(path.read_text(encoding='utf8'))
    author_status = json.loads((folder / 'results.json').read_bytes())['status']
    assert author_status == 'author_checks_passed_not_acceptance'
    reviewer_status = json.loads((folder / 'independent_results.json').read_bytes())['status']
    assert reviewer_status.startswith('PASS')
    replay = json.loads((STAGE / '_shared/replay_after1319.json').read_bytes())
    assert replay['status'] == 'PASS' and len(replay['checks']) == 2
    for check in replay['checks']:
        assert check['exit_code'] == 0
        assert sha((ROOT / check['script']).read_bytes()) == check['script_sha256']
    package.append(review_path)
    acceptance = dict(round=1319, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active', executed_checks=replay['checks'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
        main_result='A declared reversible antisymmetric-pair interface contract selects internal capacity three; the same finite fermion-pair/boson vertex realizes its full handoff and response, while SO(3) remains a genuine physical-group competitor.',
        scope_limits=['The complete antisymmetric code, same-capacity dual output, no side channel and cyclic amplitude compatibility are additional hypotheses.',
                      'Capacity three is neither generation count nor spatial dimension.',
                      'The maximum tensor-preserving group is not automatically the physical gauge group; the internal phase character is retained.',
                      'Finite mode implementation does not certify a full space-time interacting theory or the joint A/B/C/J goal.'])
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = readme_raw.decode('utf8')
    assert readme.count('最新验收至1318，共六轮新研究') == 1
    readme = readme.replace('最新验收至1318，共六轮新研究', '最新验收至1319，共七轮新研究')
    readme = readme.replace('(_shared/progress_after1318.md)明确', '(_shared/progress_after1319.md)明确')
    rows = readme.splitlines(keepends=True)
    ix = next(i for i, s in enumerate(rows) if s.startswith('|[1318]'))
    rows.insert(ix + 1, '|[1319](research_note_1319.md)|额外二体对偶闭合合同选内部容量3；同顶点给完整交接与响应，但SO(3)仍可实现，未独选物理SU(3)|[独审](1319/independent_review.md) · [验收](1319/acceptance.json)|\n')
    readme = ''.join(rows)
    first, last = readme.index('下一科学编号1319'), readme.index('历史启动验收')
    readme = readme[:first] + '''下一科学编号1320，须先核上述角色合同到既有物质字段及交换许可的实际映射，尚未自动准入。停止纯接口维数、同一转换器精度和最大群的重复分类。[1319准入](_admission/round1319_admission.md)公开全部强前件；正式范围见报告和验收。

1318的三代条件仍限公开331表示与额外运行合同。1319的内部容量3不改变三／四／五代竞争；同一交接也能在SO(3)下实现，数学最大保持群不自动成为真实规范群。不得把条件性角色机制说成已由一般认知原则推出标准模型。阶段A/B/C/J继续，新公理0。

''' + readme[last:]
    ledger_add = '''
## 1319后的更新

2026-10-19。新公理0，阶段继续。额外的全反对称码、同容量对偶、双向完整恢复和无旁路合同给有限误差下的内部维数选择；三槽循环相容限制转换幅及同源响应。有限费米对—玻色父过程与SO(3)竞争共同核验，最大张量保持群没有被当成已生成的规范群。三／四／五代、质量尺度和规范动力学仍开放。

该合同不是层级认知猜想的既成推论。后继须检查同一物质作用的字段—接口—交换映射，不能只再建抽象器件。[联合进展](progress_after1319.md)。
'''
    prefix = '''
## 当前目标：从认知交互机制生成物质分工、代际与质量结构

2026-10-19，[1313起新阶段](archive_1313_/README.md)已验收1313—1319七轮，目标保持active，新公理0。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1319.md)。

[1313—1318](archive_1313_/_shared/progress_after1318.md)保留交接、交换、质量与同父响应的已有结论及剩余输入。1318的331三代选择仍依赖表示及额外运行合同，不排除原三／四／五代竞争。

[1319](archive_1313_/research_note_1319.md)检验额外的二体对偶闭合机制：完整反对称资料在同容量接口间双向恢复，能限制内部容量；同一费米对—玻色顶点给实际交接及响应。循环相容、无旁路等强前件逐项登记。SO(3)也能完成任务，最大张量保持群不自动等于物理规范群；内部维数3不是三代或空间三维。

下一编号1320先审实际物质字段与交换许可的映射，尚未自动准入；不重复接口维数、器件优化、一般张量重建或秘密共享。当前未唯一生成现实标准模型与参数，A/B/C/J联合目标继续。

上一阶段归档为[archive_1308_1312](archive_1308_1312/README.md)，[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)及[历史复算](_migration/closure_1308_1312/README.md)保留。下方旧时点记录不改。
'''
    layer, updates = copy.deepcopy(old), []
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        original = relocation_original_bytes(raw, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        content, at = original.decode(enc), entry['edits'][0]['start']
        new = (content[:at] + prefix + content[at:]).encode(enc)
        entry.update(current_sha256=sha(new), current_bytes=len(new))
        entry['edits'][0]['after'] = prefix
        assert relocation_original_bytes(new, entry) == original
        updates.append((path, raw, new))
    layer.update(latest_accepted_round=1319, next_round=1320, in_progress_round=None,
                 stage_complete=False, application_goal_status_last_observed='active')
    new_layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger_raw + ledger_add.encode('utf8')),
                (LAYER, old_raw, new_layer_raw)]
    planned = {p: n for p, _, n in updates}
    planned[acceptance_path] = encoded(acceptance)
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
        STAGE / '_shared/progress_after1319.md']
    links = 0
    for path in md:
        content = planned[path].decode('utf8') if path in planned else path.read_text(encoding='utf8')
        assert content.count('$$') % 2 == 0
        assert not re.search(r'^\s*\\[\[\]]\s*$', content, re.M), str(path)
        for target in re.findall(r'\]\(([^)\s]+)\)', content):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (path.parent / unquote(target.split('#')[0])).resolve()
            assert resolved.exists() or resolved in planned, (path, target)
            links += 1
    checks = dict(status='PASS', date='2026-10-19', reviewed_rounds=[1319],
        markdown_files=len(md), local_links=links, new_axioms=0,
        old_navigation_boundary_preserved=True, frozen_support=len(layer['unchanged_support_sha256']),
        scope='Round 1319, current stage navigation and original root boundary; not a repository-wide audit.')
    receipt = dict(status='PUBLISHED', latest_accepted_round=1319, new_scientific_rounds=[1319],
        new_axioms=0, stage_complete=False, goal_status='active',
        previous_layer_sha256=sha(old_raw), current_layer_sha256=sha(new_layer_raw),
        previous_layer_fields={k: v for k, v in old.items() if k not in ('entries', 'unchanged_support_sha256')},
        previous_navigation_entries=old['entries'],
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): b.decode('utf8')
                                  for p, b in ((readme_path, readme_raw), (ledger_path, ledger_raw))},
        acceptance_sha256=sha(encoded(acceptance)), original_1308_closure_navigation_unchanged=True)
    fresh = [(acceptance_path, encoded(acceptance)), (paths['progress_checks'], encoded(checks)),
             (paths['progress_publication'], encoded(receipt))]
    for path, raw, _ in updates:
        assert path.read_bytes() == raw, 'Concurrent modification: ' + str(path)
    assert all(not p.exists() for p, _ in fresh)
    written, created = [], []
    try:
        for path, raw, new in updates:
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
        for path, raw, new in reversed(written):
            if path.read_bytes() == new:
                atomic(path, raw)
        raise
    print(json.dumps(checks, ensure_ascii=False))


if __name__ == '__main__':
    main()
