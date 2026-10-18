"""Publish independently reviewed 1321; never rewrite frozen scientific material."""
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    folder = STAGE / '1321'
    acceptance_path = folder / 'acceptance.json'
    receipt_path = STAGE / '_shared/progress_publication_after1321.json'
    checks_path = STAGE / '_shared/progress_checks_after1321.json'
    assert not any(p.exists() for p in (acceptance_path, receipt_path, checks_path))
    InteractionGenerationLayout().verify_interaction_navigation()
    old_raw = LAYER.read_bytes()
    old = json.loads(old_raw)
    assert old['latest_accepted_round'] == 1320 and old['in_progress_round'] == 1321
    prior = json.loads((STAGE / '1320/acceptance.json').read_bytes())
    for name, digest in prior['sha256'].items():
        assert sha((ROOT / name).read_bytes()) == digest, name
    package = [folder / n for n in ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                                    'independent_check.py', 'independent_results.json')]
    package += [STAGE / 'research_note_1321.md', STAGE / '_admission/round1321_admission.md']
    review_path = folder / 'independent_review.md'
    review = review_path.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    for p in package:
        assert sha(p.read_bytes()) in review, str(p)
    replay = json.loads((STAGE / '_shared/replay_after1321.json').read_bytes())
    assert replay['status'] == 'PASS' and len(replay['checks']) == 2
    for check in replay['checks']:
        assert check['exit_code'] == 0
        assert sha((ROOT / check['script']).read_bytes()) == check['script_sha256']
    package.append(review_path)
    acceptance = dict(round=1321, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active', executed_checks=replay['checks'],
        sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in package},
        main_result='A complete declared finite-mode parent retains counterrotating terms and admits three, four and five families with variable full-rank Yukawa data under one fixed transfer protocol and uniform error budget.',
        scope_limits=['The full Fock space is relative to the declared mode and scalar-occupation cutoff, not the full continuum field theory.',
                      'Group, matter table, Higgs background, mediator and portal flavor sequence remain declared inputs.',
                      'The task carries the complete internal code, not all unknown flavor data.',
                      'Finite-model probability differences and quadratic mass data do not certify empirical parameters, an interacting pole spectrum, or the previous long-window phase instrument.',
                      'Other weak-role branches and stronger physical realization claims require their own common-object qualification.'])
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = readme_raw.decode('utf8').replace('最新验收至1320，共八轮新研究', '最新验收至1321，共九轮新研究')
    readme = readme.replace('(_shared/progress_after1320.md)明确', '(_shared/progress_after1321.md)明确')
    rows = readme.splitlines(keepends=True)
    ix = next(i for i, row in enumerate(rows) if row.startswith('|[1320]'))
    rows.insert(ix + 1, '|[1321](research_note_1321.md)|同一有限父模式保全反旋转项，固定任务下三／四／五代与满秩质量资料共同非空；未签完整场仪器|[独审](1321/independent_review.md) · [验收](1321/acceptance.json)|\n')
    readme = ''.join(rows)
    start, end = readme.index('[1321已准入]'), readme.index('历史启动验收')
    readme = readme[:start] + '''下一科学编号1322。先审[共同弱角色替代](_admission/weak_role_joint_screen_after1321.md)：保持左门户与实际任务，核带弱荷右伙伴能否在同一有限范围承接，并公开保护表示等变化；不把1313的较弱交接证书直接当作全部强合同反例。

1321补的是明确有限模式父过程的显式误差与联合非空。完整空间模式、动态规范／Higgs、旧长时相位任务及经验参数不因此自动获得认证。停止精化S器件，A/B/C/J继续，新公理0。

''' + readme[end:]
    ledger_add = '''
## 1321后的更新

2026-10-19。相同有限空间模式、完整费米字段和S粒子／反粒子下，保全部反旋转项，用982工具计算真实共同误差；三／四／五代只截取事前固定味序列，不逐代重调耦合、时钟或阈值。共同小质量域保全秩与可变混合，同父短时占据预测可分离。

这补充了有限模式父作用的共同非空，未消除标准字段／规范群、Higgs背景、额外门户及截止选择。三维任务不承载全部未知味资料；新有限预测不替代旧长窗相位仪器或经验参数。下一步核强合同下的弱角色替代，避免拼接1313与1320不同范围。[联合进展](progress_after1321.md)。
'''
    prefix = '''
## 当前目标：从认知交互机制生成物质分工、代际与质量结构

2026-10-19，[1313起新阶段](archive_1313_/README.md)已验收1313—1321九轮，目标保持active，新公理0。[正式目标](archive_1313_/_shared/goal.md) · [联合进展](archive_1313_/_shared/progress_after1321.md)。

[1321](archive_1313_/research_note_1321.md)在同一明确有限模式父作用中保全部反旋转项，以固定门户、时钟和误差同时承接三／四／五代及可变满秩质量资料；不再只并列一个亮态和另一份质量证书。其有限截止、Higgs背景及物理群仍为输入，未认证完整场仪器或现实参数。

下一编号1322先审[共同弱角色替代](archive_1313_/_admission/weak_role_joint_screen_after1321.md)，不能把1313的弱任务反例自动提升为1319—1320全部强合同反例。角色、保护、质量及任务须在同一新分支核对，A/B/C/J继续。

上一阶段已归档为[archive_1308_1312](archive_1308_1312/README.md)，[论文14](14_单参考质量与有限交换任务下的规范结构选择_阶段论文.md)及[历史复算](_migration/closure_1308_1312/README.md)保留。下方旧时点记录不改。
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
    layer.update(latest_accepted_round=1321, next_round=1322, in_progress_round=None,
                 stage_complete=False, application_goal_status_last_observed='active')
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger_raw + ledger_add.encode('utf8')),
                (LAYER, old_raw, encoded(layer))]
    planned = {p: after for p, _, after in updates}
    planned[acceptance_path] = encoded(acceptance)
    links = 0
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
        STAGE / '_shared/progress_after1321.md', STAGE / '_admission/weak_role_joint_screen_after1321.md']
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
    checks = dict(status='PASS', round=1321, markdown_files=len(md), local_links=links,
        prior_1320_frozen_package_unchanged=True, original_navigation_recoverable=True,
        scope='1321 package and affected navigation only')
    receipt = dict(status='PUBLISHED', latest_accepted_round=1321, new_scientific_rounds=[1321],
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
