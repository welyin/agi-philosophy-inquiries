"""Review-bound round-1323 publication; default --check is read-only."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
from urllib.parse import unquote
from start_stage import ROOT, STAGE, LAYER_PATH, sha, require, atomic_replace, publication_lock
from research_layout import relocation_original_bytes
from run_research_interaction_rules_active import InteractionRulesLayout

FOLDER = STAGE / '1323'
ACCEPTANCE = FOLDER / 'acceptance.json'
CHECKS = STAGE / '_shared/progress_checks_after1323.json'
RECEIPT = STAGE / '_shared/progress_publication_after1323.json'
EXPECTED = {
 '1323/proof.md': '1cb7867b2596b275e094751451b03aecc463d5ea0639c09fdcfc4b8cdbed444e',
 '1323/check.py': 'bc3e5b5855f19cc4403416b03739ad3d52f46ed3565e14f525d2cafd5ce42478',
 '1323/results.json': '821a01fca35e63cb1d92d0cab3226600655b569eb29fb86e13ff8b4eeb61f1ad',
 '1323/NEXT.md': '9e5a2828814ce2d65c75bcf36e2ed5250dcb540a2b250ba76974c457f0fe4567',
 'research_note_1323.md': 'e5f92830adc097a55ffe426cef9355094794877baf7ff71254915735fe89e871',
 '_admission/round1323_admission.md': '87796eaa0616b65afbec413431da998b21f01642169fe2c9c38fd8794a394b3e',
 '1323/independent_check.py': '09f1d6e862db69664942c5018ed65bf7417c8a1ee4cff4060d275c4acfd53984',
 '1323/independent_results.json': '66906c918d0c9ee6e29dff747c74145ce3cc43a0e2e2f520ad5ceed92acfbd51',
 '1323/independent_review.md': '61e03374bf46ab3e0117f9fefc8b81eb4af1e3bfbe3c9a4087a33b7d048819df',
 '1323/scope_review.md': 'eac4b569c6d865174b2ba59a2f9d6f3cea558d668a8975f1d91a9a9ffcd2e033',
}

def encoded(value):
    return (json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf8')

def source_map(paths):
    return {p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in paths}

def prepare():
    require(STAGE.relative_to(ROOT).as_posix() == 'research_cognition_physics/archive_1323_', 'Wrong stage')
    for p in (ACCEPTANCE, CHECKS, RECEIPT):
        require(not p.exists(), 'Already published or incomplete prior publication: ' + str(p))
    store = InteractionRulesLayout()
    prior_navigation = store.verify_interaction_rules_navigation()
    old_raw = LAYER_PATH.read_bytes()
    old = json.loads(old_raw)
    require(old['latest_accepted_round'] == 1322 and old['in_progress_round'] in (None, 1323)
            and old['next_round'] == 1323 and not old['stage_complete'], 'Scientific boundary differs')
    goal = STAGE / '_shared/goal.md'
    goal_receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_bytes())
    require(goal_receipt['status'] == 'active' and
            sha(goal.read_bytes()) == goal_receipt['goal_document_sha256'], 'Reviewed goal changed')
    package = [STAGE / name for name in EXPECTED]
    require(len(package) == 10, 'Unexpected evidence count')
    for name, digest in EXPECTED.items():
        require(sha((STAGE / name).read_bytes()) == digest, 'Frozen source changed: ' + name)
    math_review = (FOLDER / 'independent_review.md').read_text(encoding='utf8')
    scope_review = (FOLDER / 'scope_review.md').read_text(encoding='utf8')
    require('PASS_WITH_SCOPE' in math_review and 'PASS_WITH_SCOPE' in scope_review, 'Reviews must pass')
    for name in list(EXPECTED)[:8]:
        require(EXPECTED[name] in math_review, 'Math binding missing: ' + name)
    for name in list(EXPECTED)[:9]:
        require(EXPECTED[name] in scope_review, 'Scope binding missing: ' + name)
    require(sha(goal.read_bytes()) in scope_review, 'Scope goal binding missing')
    replay_path = STAGE / '_shared/replay_after1323.json'
    replay = json.loads(replay_path.read_bytes())
    require(replay['status'] == 'PASS' and replay['round'] == 1323 and len(replay['checks']) == 2
            and replay['new_axioms'] == 0 and replay['stage_complete'] is False, 'Replay incomplete')
    scripts = {(FOLDER / n).relative_to(ROOT).as_posix() for n in ('check.py', 'independent_check.py')}
    require({c['script'] for c in replay['checks']} == scripts, 'Wrong replay scripts')
    for c in replay['checks']:
        require(c['exit_code'] == 0 and c['checks'] > 0 and bool(c['receipt']), 'Replay failed')
        require(sha((ROOT / c['script']).read_bytes()) == c['script_sha256'], 'Replay source changed')
    acceptance = dict(round=1323, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active',
        next_available_round=1324, next_round_admitted=False, scientific_evidence_files=10,
        executed_checks=replay['checks'], sha256=source_map(package),
        main_result='A full unknown-reference extraction contract selects one right-singlet port and a normalized finite-mode interaction. A fixed cross-balance clock distinguishes it from bare linear annihilation; sequential extraction and source balances share the same finite parent.',
        scope_limits=[
            'The full-unknown-reference and short-time balance rules are candidate roles, not adopted axioms.',
            'The capacity theorem concerns one descending port, not all Standard Model right species.',
            'N=1..4 and n=3,4,5 have finite parent witnesses; arbitrary M need not permit perfect transfer.',
            'Finite-window interpolation does not determine a unique global nonpolynomial rule or local field Wilson coefficients.',
            'Failures, partners, preparation, clocks and sequential access remain explicitly accounted for.',
            'No full Higgs/gauge realization, complete 1321/1322 CAR matching, family uniqueness or empirical parameters are certified.'])
    readme_path, ledger_path = STAGE / 'README.md', STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    old_readme = readme_raw.decode('utf8')
    require(old_readme.count('## 启动状态') == 1, 'Unexpected README boundary')
    readme = old_readme.split('## 启动状态', 1)[0] + '''## 当前进展

2026-10-19。首轮1323已独立审阅并验收，新增公理0；目标保持active，尚未阶段结项。[联合进展](_shared/progress_after1323.md)明确A／B／C／J的已证范围与未决接口。

|轮次|新增结果|证据|
|---|---|---|
|[1323](research_note_1323.md)|全未知参考逐份交接选择单下降端口的右角色；短时余额一致性选择归一提取，并有同一有限父模型的组合与来源|[数学独审](1323/independent_review.md) · [范围独审](1323/scope_review.md) · [验收](1323/acceptance.json)|

该合同为待检机制，不将右单态直接等同全部右物种，也不把模式提取当作完整局域Higgs作用。三／四／五味仍共同非空，阶段目标的物质层和共同场字典继续推进。

[下一步](1323/NEXT.md)先筛同一实际参考能否给有效物质响应和参考后态的共同连接。1324为下一个可用科学编号，尚未准入；不重复容量、单余额交换或无选择力的截止扫描。

[发布前图像与事务记录](_shared/progress_publication_after1323.json)保留本阶段原导航及旧层；[启动记录](_shared/navigation_start_receipt.json)仍是计0历史，不增加科学轮次。
'''
    ledger = ledger_raw.decode('utf8')
    before_intro = '2026-10-19。启动账，计0；当前科学成果仍止于1322。'
    require(before_intro in ledger, 'Unexpected ledger boundary')
    ledger = ledger.replace(before_intro, '2026-10-19。下表保留启动时的输入基线；科学成果现验收至1323，最新判定见文末。', 1)
    ledger += '''
## 1323后的更新

[报告](../research_note_1323.md)与[联合进展](progress_after1323.md)给出更强的全未知参考合同。保整个参考、右端及任意伙伴，且成功后无未登记量子旁路时，单下降端口只能承接一维右角色；含N=1的小误差合同亦排除双重态。该量词超出1322原固定背景合同，不推翻旧结果。

另行登记的短时余额一致性在同一纯转换父类内锁定味生成元，选择归一化提取而非普通线性湮灭；顺序组合系数由原作用和实际余额共同决定。有限窗可多项式插值，未证明全场规则唯一。n=3、4、5及一般未由经验值固定的味矩阵继续保留。

A／B有实质局部进展，C的全部右物种与跨部门参数、J的实际局域场匹配尚未完成。1324尚未准入，目标active，不新增公理。[正式验收](../1323/acceptance.json)。
'''
    prefix = '''
## 当前目标：基本相互作用规则的认知来源

2026-10-19，[1323起阶段](archive_1323_/README.md)首轮1323已验收，新增公理0，目标保持active。[正式目标](archive_1323_/_shared/goal.md) · [1323报告](archive_1323_/research_note_1323.md) · [联合进展](archive_1323_/_shared/progress_after1323.md)。

本轮把原参考本身的全部未知资料及伙伴关联纳入逐份交接。在指定无旁路接收合同下，单个下降端口的右角色只能一维；短时余额一致性选择归一提取，并不选择普通线性Higgs菜单。同一有限父过程给出固定校准后的余额差异、完整来源和两次接续，三／四／五味仍共同非空。

这些是公开额外机制下的条件结论，不等于全部标准模型右物种、局域场作用或现实参数已经生成。下一项先筛同一实际参考与物质响应／后态的共同连接，1324尚未准入；阶段C／J及最终联合验收继续，不因首轮角色选择提前结项。

上一阶段保存在[archive_1313_1322](archive_1313_1322/README.md)及[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)。[迁移与历史复算](_migration/closure_1313_1322/README.md)和下方历史记录保留。当前安排以本节与现行目标为准。
'''
    layer, updates = copy.deepcopy(old), []
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        before = path.read_bytes()
        require(sha(before) == entry['current_sha256'], 'Concurrent root navigation edit')
        original = relocation_original_bytes(before, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        text, at = original.decode(enc), entry['edits'][0]['start']
        require(at == text.index('\n') + 1, 'Wrong insertion position')
        after = (text[:at] + prefix + text[at:]).encode(enc)
        entry['edits'][0]['after'] = prefix
        entry.update(current_sha256=sha(after), current_bytes=len(after))
        require(relocation_original_bytes(after, entry) == original, 'Nonreversible prefix')
        updates.append((path, before, after))
    layer.update(latest_accepted_round=1323, in_progress_round=None, next_round=1324,
        next_round_admitted=False, stage_complete=False, application_goal_status_last_observed='active',
        startup_state='Round 1323 accepted with conditional scope; next round not admitted')
    layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger.encode('utf8')), (LAYER_PATH, old_raw, layer_raw)]
    planned = {p.resolve(): after for p, before, after in updates}
    planned[ACCEPTANCE.resolve()] = encoded(acceptance)
    fresh_paths = {CHECKS.resolve(), RECEIPT.resolve()}
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
          STAGE / '_shared/progress_after1323.md', goal]
    links = 0
    for path in md:
        text = planned.get(path.resolve(), path.read_bytes()).decode('utf8')
        require(len(re.findall(r'^\s*\$\$\s*$', text, re.M)) % 2 == 0, 'Unbalanced block formula: ' + str(path))
        require(not re.search(r'^\s*\\[\[\]]\s*$', text, re.M), 'Old block delimiters')
        require(not any(ord(c) < 32 and c not in '\n\r\t' for c in text), 'Control character')
        for target in re.findall(r'\]\(([^)\s]+)\)', text):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (path.parent / unquote(target.split('#')[0])).resolve()
            require(resolved.is_file() or resolved in planned or resolved in fresh_paths,
                    'Missing local link: ' + str(path) + ' -> ' + target)
            links += 1
    sources = source_map(package + [goal, replay_path, STAGE / '_shared/progress_after1323.md',
        STAGE / '_shared/start_stage.py', STAGE / '_shared/goal_start_receipt.json',
        ROOT / 'scripts/run_research_interaction_rules_active.py', Path(__file__).resolve()])
    checks = dict(status='PASS', round=1323, scientific_evidence_files=10,
        markdown_files=len(md), local_links=links, original_navigation_recoverable=True,
        prior_navigation=prior_navigation, scientific_files_modified=False,
        stage_complete=False, next_round_admitted=False)
    receipt = dict(status='PUBLISHED', latest_accepted_round=1323, new_scientific_rounds=[1323],
        new_axioms=0, goal_status='active', stage_complete=False, next_available_round=1324,
        next_round_admitted=False, previous_layer_sha256=sha(old_raw),
        current_layer_sha256=sha(layer_raw), previous_layer=old,
        acceptance_sha256=sha(encoded(acceptance)), source_sha256=sources,
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): raw.decode('utf8')
            for p, raw in ((readme_path, readme_raw), (ledger_path, ledger_raw))})
    return dict(updates=updates, acceptance=encoded(acceptance), checks=checks,
                receipt=receipt, sources=sources)

def publish_locked():
    plan = prepare()
    for path, before, after in plan['updates']:
        require(path.read_bytes() == before, 'Concurrent navigation edit: ' + str(path))
    for name, digest in plan['sources'].items():
        require(sha((ROOT / name).read_bytes()) == digest, 'Reviewed source changed: ' + name)
    written, created = [], []
    def create_fresh(path, raw):
        with path.open('xb') as stream:
            created.append((path, raw))
            stream.write(raw)
    try:
        for path, before, after in plan['updates']:
            require(path.read_bytes() == before, 'Concurrent navigation edit: ' + str(path))
            atomic_replace(path, after)
            written.append((path, before, after))
        create_fresh(ACCEPTANCE, plan['acceptance'])
        for name, digest in plan['sources'].items():
            require(sha((ROOT / name).read_bytes()) == digest, 'Source changed during publication: ' + name)
        plan['checks']['published_navigation'] = InteractionRulesLayout().verify_interaction_rules_navigation()
        checks_raw = encoded(plan['checks'])
        plan['receipt']['checks_sha256'] = sha(checks_raw)
        create_fresh(CHECKS, checks_raw)
        create_fresh(RECEIPT, encoded(plan['receipt']))
    except BaseException:
        for path, expected in reversed(created):
            if path.exists() and expected.startswith(path.read_bytes()):
                path.unlink()
        for path, before, after in reversed(written):
            if path.read_bytes() == after:
                atomic_replace(path, before)
        raise
    print(json.dumps(plan['checks'], ensure_ascii=False))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true', help='Read-only readiness check (default)')
    group.add_argument('--publish', action='store_true', help='Publish reviewed evidence and navigation once')
    args = parser.parse_args()
    if args.publish:
        with publication_lock():
            publish_locked()
    else:
        plan = prepare()
        print(json.dumps(dict(plan['checks'], status='READY', files_written=0), ensure_ascii=False))

if __name__ == '__main__':
    main()
