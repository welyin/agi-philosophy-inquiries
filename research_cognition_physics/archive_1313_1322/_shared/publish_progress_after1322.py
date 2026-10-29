"""Publish only the independently reviewed 1322 package; no scientific edits."""
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    folder = STAGE / '1322'
    acceptance_path = folder / 'acceptance.json'
    receipt_path = STAGE / '_shared/progress_publication_after1322.json'
    checks_path = STAGE / '_shared/progress_checks_after1322.json'
    assert not any(p.exists() for p in (acceptance_path, receipt_path, checks_path))
    InteractionGenerationLayout().verify_interaction_navigation()
    old_raw = LAYER.read_bytes()
    old = json.loads(old_raw)
    assert old['latest_accepted_round'] == 1321 and old['in_progress_round'] == 1322
    for n in range(1313, 1322):
        prior = json.loads((STAGE / str(n) / 'acceptance.json').read_bytes())
        for name, digest in prior['sha256'].items():
            assert sha((ROOT / name).read_bytes()) == digest, name
    package = [folder / n for n in ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                                    'independent_check.py', 'independent_results.json')]
    package += [STAGE / 'research_note_1322.md', STAGE / '_admission/round1322_admission.md']
    review_path = folder / 'independent_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    for p in package:
        assert sha(p.read_bytes()) in review, str(p)
    replay = json.loads((STAGE / '_shared/replay_after1322.json').read_bytes())
    assert replay['status'] == 'PASS' and len(replay['checks']) == 2
    for check in replay['checks']:
        assert check['exit_code'] == 0
        assert sha((ROOT / check['script']).read_bytes()) == check['script_sha256']
    package.append(review_path)
    acceptance = dict(round=1322, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active', executed_checks=replay['checks'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
        main_result='An explicit vectorlike right-weak role transports the complete fixed-background finite parent and strong handoff contract, while right gauge, Higgs-radial and effective-scale sources remain distinguishable.',
        scope_limits=[
            'The comparison keeps R_op, not fixed unbroken R_rep, and opens a dimension-five mass menu.',
            'Exact finite-parent equivalence holds at constant nonzero H and zero gauge background.',
            'A relation-variable charge label does not eliminate actual right weak currents.',
            'Neutral readout requires the H-defined relational projector permission.',
            'No full dynamical gauge/Higgs apparatus, empirical parameters or UV completion is certified.'])
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = readme_raw.decode('utf8').replace('最新验收至1321，共九轮新研究', '最新验收至1322，共十轮新研究')
    readme = readme.replace('(_shared/progress_after1321.md)明确', '(_shared/progress_after1322.md)明确')
    rows = readme.splitlines(keepends=True)
    ix = next(i for i, row in enumerate(rows) if row.startswith('|[1321]'))
    rows.insert(ix + 1, '|[1322](research_note_1322.md)|同一固定背景强交接可由右弱双重态承接；保护菜单及规范／Higgs来源区别明确|[独审](1322/independent_review.md) · [验收](1322/acceptance.json)|\n')
    readme = ''.join(rows)
    start, end = readme.index('[1322已准入]'), readme.index('历史启动验收')
    readme = readme[:start] + '''科学轮次验收至1322，进入[阶段论文](../15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)与A/B/C/J逐项联合验收。当前目标仍active，尚不自动登记阶段结项；下一科学轮未准入。新增公理0。

主链1319—1322已经给出声明范围内的共同非空与可区别剩余族。其它机制分支分别保留自己的前件，不拼成一套标准模型；完整现实参数和强合同下全部超荷分类仍未获得。

''' + readme[end:]
    ledger_add = '''
## 1322后的更新

2026-10-19。保持左门户、H背景和强合同，向量样右弱角色的完整有限父算符可运输至1321字典；未破缺保护表示和维数五菜单公开变更。该分支满足R_op，不反证固定R_rep或维数四菜单的结论。

背景任务等价仍保不同右W/Z、径向H和Λ来源；关系中性读口另需H权限。主链共同范围与其它分支分开，进入阶段论文和逐项验收，不新增科学编号或公理。[联合进展](progress_after1322.md)。
'''
    prefix = '''
## 当前目标：认知交互机制的联合验收

2026-10-19，[1313起阶段](archive_1313_/README.md)已验收1313—1322十轮，新公理0。目标仍active，正在完成[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)及逐项独审。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1322.md)。

主链取得内部容量与角色限制、固定交接下的共同系数关系，以及同一有限父模型的三／四／五代和质量剩余族。1322补强合同中的右弱角色替代：任务可相同，实际弱流、Higgs响应与保护菜单仍可区别。未唯一生成完整标准模型或现实参数。

论文将主链与Dirac-only／批量交换、内部链和331分支分别列明；未分类的群／超荷保为开放项，不包装成已证明不足。下一科学轮尚未准入，阶段是否结项以联合验收为准。

上一阶段已归档为[archive_1308_1312](archive_1308_1312/README.md)。下方旧时点记录不改。
'''
    layer, updates = copy.deepcopy(old), []
    for entry in layer['entries']:
        p = ROOT / entry['destination']
        before = p.read_bytes()
        original = relocation_original_bytes(before, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        content, at = original.decode(enc), entry['edits'][0]['start']
        after = (content[:at] + prefix + content[at:]).encode(enc)
        entry.update(current_sha256=sha(after), current_bytes=len(after))
        entry['edits'][0]['after'] = prefix
        assert relocation_original_bytes(after, entry) == original
        updates.append((p, before, after))
    layer.update(latest_accepted_round=1322, next_round=None, in_progress_round=None,
        stage_complete=False, application_goal_status_last_observed='active',
        startup_state='Ten scientific rounds accepted; joint paper acceptance in progress')
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger_raw + ledger_add.encode('utf8')),
                (LAYER, old_raw, encoded(layer))]
    planned = {p: after for p, _, after in updates}
    planned[acceptance_path] = encoded(acceptance)
    links = 0
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
        STAGE / '_shared/progress_after1322.md']
    for p in md:
        content = planned[p].decode('utf8') if p in planned else p.read_text(encoding='utf8')
        assert content.count('$$') % 2 == 0, str(p)
        assert not re.search(r'^\s*\\[\[\]]\s*$', content, re.M), str(p)
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in content), str(p)
        for target in re.findall(r'\]\(([^)\s]+)\)', content):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (p.parent / unquote(target.split('#')[0])).resolve()
            assert resolved.exists() or resolved in planned, (p, target)
            links += 1
    checks = dict(status='PASS', round=1322, markdown_files=len(md), local_links=links,
        prior_scientific_rounds_unchanged=list(range(1313, 1322)), original_navigation_recoverable=True)
    receipt = dict(status='PUBLISHED', latest_accepted_round=1322, new_scientific_rounds=[1322],
        goal_status='active', stage_complete=False, new_axioms=0,
        previous_layer_sha256=sha(old_raw), current_layer_sha256=sha(encoded(layer)),
        previous_layer=old, acceptance_sha256=sha(encoded(acceptance)),
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): b.decode('utf8')
                                   for p, b in ((readme_path, readme_raw), (ledger_path, ledger_raw))})
    fresh = [(acceptance_path, encoded(acceptance)), (checks_path, encoded(checks)),
             (receipt_path, encoded(receipt))]
    for p, before, _ in updates:
        assert p.read_bytes() == before, str(p)
    written, created = [], []
    try:
        for p, before, after in updates:
            atomic(p, after)
            written.append((p, before, after))
        for p, value in fresh:
            with p.open('xb') as stream:
                stream.write(value)
            created.append((p, value))
    except BaseException:
        for p, value in reversed(created):
            if p.read_bytes() == value:
                p.unlink()
        for p, before, after in reversed(written):
            if p.read_bytes() == after:
                atomic(p, before)
        raise
    print(json.dumps(checks))


if __name__ == '__main__':
    main()
