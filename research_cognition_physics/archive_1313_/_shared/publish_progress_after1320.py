"""Publish reviewed round 1320 once, preserving frozen originals and evidence."""
import ast
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    folder = STAGE / '1320'
    acceptance_path = folder / 'acceptance.json'
    receipt_path = STAGE / '_shared/progress_publication_after1320.json'
    checks_path = STAGE / '_shared/progress_checks_after1320.json'
    assert not any(p.exists() for p in (acceptance_path, receipt_path, checks_path))
    InteractionGenerationLayout().verify_interaction_navigation()
    old_raw = LAYER.read_bytes()
    old = json.loads(old_raw)
    assert old['latest_accepted_round'] == 1319 and old['in_progress_round'] == 1320
    previous_acceptance = json.loads((STAGE / '1319/acceptance.json').read_bytes())
    for name, digest in previous_acceptance['sha256'].items():
        assert sha((ROOT / name).read_bytes()) == digest, name
    package = [folder / n for n in ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                                    'independent_check.py', 'independent_results.json')]
    package += [STAGE / 'research_note_1320.md', STAGE / '_admission/round1320_admission.md']
    review_path = folder / 'independent_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    for p in package:
        assert sha(p.read_bytes()) in review, str(p)
        if p.suffix == '.py':
            ast.parse(p.read_text(encoding='utf8'))
    for filename in ('results.json', 'independent_results.json'):
        status = json.loads((folder / filename).read_bytes())['status']
        assert status in ('PASS', 'PASS_WITH_SCOPE', 'author_checks_passed_not_acceptance'), status
    replay = json.loads((STAGE / '_shared/replay_after1320.json').read_bytes())
    assert replay['status'] == 'PASS' and len(replay['checks']) == 2
    for check in replay['checks']:
        assert check['exit_code'] == 0
        assert sha((ROOT / check['script']).read_bytes()) == check['script_sha256']
    package.append(review_path)
    acceptance = dict(round=1320, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active', executed_checks=replay['checks'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
        main_result='The complete endpoint recovery contract restricts an initially open internal role W to a character; full-rank intertwining and connectedness identify it. A fixed simultaneous two-portal parent links transfer error to its common response coefficients.',
        scope_limits=['Complete internal code and limited receiving capacity are additional task conditions; other weak, spin, flavor and momentum slots are separately declared.',
                      'Full-rank covariance and connectedness are required to identify the character.',
                      'The strength-ratio bound assumes constant simultaneous portals; independent sequential control is a genuine alternative.',
                      'The finite CAR parent and leading local-field matching have an explicit conditional mode bridge, not a certificate for arbitrary space-time apparatus.',
                      'The mediator, restricted vertex menu, overlaps and flavor data remain inputs; no unique physical group, family count or empirical parameters are claimed.'])
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = readme_raw.decode('utf8')
    assert readme.count('最新验收至1319，共七轮新研究') == 1
    readme = readme.replace('最新验收至1319，共七轮新研究', '最新验收至1320，共八轮新研究')
    readme = readme.replace('(_shared/progress_after1319.md)明确', '(_shared/progress_after1320.md)明确')
    rows = readme.splitlines(keepends=True)
    ix = next(i for i, s in enumerate(rows) if s.startswith('|[1319]'))
    rows.insert(ix + 1, '|[1320](research_note_1320.md)|完整双端合同限制开放角色W；固定双门户任务约束同源系数，独立顺序控制可消除强度比限制|[独审](1320/independent_review.md) · [验收](1320/acceptance.json)|\n')
    readme = ''.join(rows)
    first, last = readme.index('[1320已准入]'), readme.index('历史启动验收')
    readme = readme[:first] + '''下一科学编号1321，先审[同合同的联合剩余自由](_admission/joint_residual_screen_after1320.md)，未自动准入。保持角色、权限、码和误差，检查三／四／五代及质量／混合的共同非空，不能拼接1314、1315、1318的不同分支。

1320的强度关系限固定、同时开放的两门户；允许独立顺序控制时，一般交接不要求该比例。局域父作用与有限过程的模式桥、其余物质角色和规范权限仍须按各自范围登记；不把单个候选的全部空间装置作为阶段门槛。停止本中介优化，A/B/C/J继续，新公理0。

''' + readme[last:]
    ledger_add = '''
## 1320后的更新

2026-10-19。开放L内部角色W后，完整双端恢复容量合同给dim W=1；满秩交织与连通群再识别W*=D。该结论不预置L颜色单态，也不替代右弱单态或单参考等其它条件。固定同时双门户的实际有限交付，将误差与同父原始电流Gram系数联系；顺序独立开关解除该强度比约束。S、味耦合、模式重叠、顶点菜单及物理群仍输入。

完整有限CAR过程与局域领先匹配分别核验，模式／遗漏作用桥条件明确；未宣称现实装置与完整SM已经生成。下一步冻结共同合同并作跨代际／参数联合剩余自由检验，不重复中介优化。[联合进展](progress_after1320.md)。
'''
    prefix = '''
## 当前目标：从认知交互机制生成物质分工、代际与质量结构

2026-10-19，[1313起新阶段](archive_1313_/README.md)已验收1313—1320八轮，目标保持active，新公理0。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1320.md)。

[1319](archive_1313_/research_note_1319.md)在额外二体对偶合同下选择内部容量，并连接实际交付与同核响应；SO(3)竞争仍保留，未独选真实规范群或三代。[1320](archive_1313_/research_note_1320.md)开放另一物质角色后，以完整恢复及交织条件限制其内部表示；固定双门户的交付误差约束共同交换系数。独立顺序控制可解除比例限制，不能把此关系扩大到一般交接。

下一编号1321先审[联合剩余自由](archive_1313_/_admission/joint_residual_screen_after1320.md)，未自动准入：固定同一合同，检查代际及质量／混合的共同非空，停止单个中介优化。不同轮次的互斥分支不能拼接；物理群、角色前件、模式条件及现实参数仍按来源登记，A/B/C/J目标继续。

上一阶段已归档为[archive_1308_1312](archive_1308_1312/README.md)，[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)及[历史复算](_migration/closure_1308_1312/README.md)保留。下方旧时点记录不改。
'''
    layer, updates = copy.deepcopy(old), []
    for entry in layer['entries']:
        p = ROOT / entry['destination']
        raw = p.read_bytes()
        original = relocation_original_bytes(raw, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        content, at = original.decode(enc), entry['edits'][0]['start']
        new = (content[:at] + prefix + content[at:]).encode(enc)
        entry.update(current_sha256=sha(new), current_bytes=len(new))
        entry['edits'][0]['after'] = prefix
        assert relocation_original_bytes(new, entry) == original
        updates.append((p, raw, new))
    layer.update(latest_accepted_round=1320, next_round=1321, in_progress_round=None,
                 stage_complete=False, application_goal_status_last_observed='active')
    new_layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger_raw + ledger_add.encode('utf8')),
                (LAYER, old_raw, new_layer_raw)]
    planned = {p: n for p, _, n in updates}
    planned[acceptance_path] = encoded(acceptance)
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
        STAGE / '_shared/progress_after1320.md']
    links = 0
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
    checks = dict(status='PASS', date='2026-10-19', reviewed_rounds=[1320],
        markdown_files=len(md), local_links=links, new_axioms=0,
        prior_1319_frozen_package_unchanged=True, original_navigation_recoverable=True,
        scope='Round 1320 and current stage navigation; not a repository-wide audit.')
    receipt = dict(status='PUBLISHED', latest_accepted_round=1320, new_scientific_rounds=[1320],
        new_axioms=0, stage_complete=False, goal_status='active',
        previous_layer_sha256=sha(old_raw), current_layer_sha256=sha(new_layer_raw),
        previous_layer_fields={k: v for k, v in old.items() if k not in ('entries', 'unchanged_support_sha256')},
        previous_navigation_entries=old['entries'],
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): b.decode('utf8')
                                  for p, b in ((readme_path, readme_raw), (ledger_path, ledger_raw))},
        acceptance_sha256=sha(encoded(acceptance)))
    fresh = [(acceptance_path, encoded(acceptance)), (checks_path, encoded(checks)),
             (receipt_path, encoded(receipt))]
    for p, raw, _ in updates:
        assert p.read_bytes() == raw, 'Concurrent modification: ' + str(p)
    written, created = [], []
    try:
        for p, raw, new in updates:
            atomic(p, new)
            written.append((p, raw, new))
        for p, raw in fresh:
            with p.open('xb') as stream:
                stream.write(raw)
            created.append((p, raw))
    except BaseException:
        for p, raw in reversed(created):
            if p.read_bytes() == raw:
                p.unlink()
        for p, raw, new in reversed(written):
            if p.read_bytes() == new:
                atomic(p, raw)
        raise
    print(json.dumps(checks, ensure_ascii=False))


if __name__ == '__main__':
    main()
