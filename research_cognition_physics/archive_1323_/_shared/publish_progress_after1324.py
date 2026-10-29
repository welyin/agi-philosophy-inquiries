"""Review-bound round-1324 publication; default --check is read-only."""
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

FOLDER = STAGE / '1324'
ACCEPTANCE = FOLDER / 'acceptance.json'
CHECKS = STAGE / '_shared/progress_checks_after1324.json'
RECEIPT = STAGE / '_shared/progress_publication_after1324.json'
EXPECTED = {'1324/proof.md': '6c44c221014ae1c26e2ad1a870114cccb0c901b68f859e73a1ee25508b7028a1', '1324/check.py': '28398cf2e82a03eb384ed4d60dcf4d3acc85166c79ad9795abf02148be41520b', '1324/results.json': '12b56fe1e07f642d2c6eeb3605a1df5ea1ad00ebde5f1a6ae842d91959b6acda', '1324/NEXT.md': 'f64741052a6955e654b35b131c76a89e7ba7447e55a3cfd37969bc6f7301ed3b', 'research_note_1324.md': 'b385571ada40ac55f97822a2497922f1ac098d81b45c054513de1720ed590faf', '_admission/round1324_admission.md': 'b6f76e23c740414c2949c9c6f03beb03d6a99b223737bcda9dee1dbae4388d7c', '1324/independent_check.py': '1384a7408f5a1443402c5c7026debcd43bc21210762275eb9b08366ad22aed5b', '1324/independent_results.json': '0f26084a05b059b4af5e61a0445d17ea4011a6cbee2c3088fc5cc6c01b85aa1a', '1324/independent_review.md': 'f7f4510e08312c0c3f54dee0515337afff8feb7c1481c514c3a4d6c41c1900ee', '1324/scope_review.md': 'b1ba6259a7966c2176a3610acfa2182e2976296db51633994e83424a9e62b0b3'}


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
    require(old['latest_accepted_round'] == 1323 and old['in_progress_round'] in (None, 1324)
            and old['next_round'] == 1324 and not old['stage_complete'], 'Scientific boundary differs')
    goal = STAGE / '_shared/goal.md'
    goal_receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_bytes())
    require(goal_receipt['status'] == 'active' and
            sha(goal.read_bytes()) == goal_receipt['goal_document_sha256'], 'Reviewed goal changed')
    package = [STAGE / name for name in EXPECTED]
    require(len(package) == 10, 'Unexpected evidence count')
    for name, digest in EXPECTED.items():
        require(sha((STAGE / name).read_bytes()) == digest, 'Frozen source changed: ' + name)
    prior_acceptance = STAGE / '1323/acceptance.json'
    prior_receipt = STAGE / '_shared/progress_publication_after1323.json'
    prior_checks = STAGE / '_shared/progress_checks_after1323.json'
    prior = json.loads(prior_acceptance.read_bytes())
    require(prior['round'] == 1323 and prior['status'] == 'ACCEPTED_CONDITIONAL_SCOPE',
            'Previous round is not accepted')
    for name, digest in prior['sha256'].items():
        require(sha((ROOT / name).read_bytes()) == digest, 'Previous frozen science changed: ' + name)
    prior_publication = json.loads(prior_receipt.read_bytes())
    require(sha(prior_acceptance.read_bytes()) == prior_publication['acceptance_sha256']
            and sha(prior_checks.read_bytes()) == prior_publication['checks_sha256'],
            'Previous publication evidence changed')
    require(sha(old_raw) == prior_publication['current_layer_sha256'],
            'Current layer is not the exact previous publication boundary')
    math_review = (FOLDER / 'independent_review.md').read_text(encoding='utf8')
    scope_review = (FOLDER / 'scope_review.md').read_text(encoding='utf8')
    require('PASS_WITH_SCOPE' in math_review and 'PASS_WITH_SCOPE' in scope_review, 'Reviews must pass')
    for name in list(EXPECTED)[:8]:
        require(EXPECTED[name] in math_review, 'Math binding missing: ' + name)
    for name in list(EXPECTED)[:9]:
        require(EXPECTED[name] in scope_review, 'Scope binding missing: ' + name)
    require(sha(goal.read_bytes()) in scope_review, 'Scope goal binding missing')
    replay_path = STAGE / '_shared/replay_after1324.json'
    replay = json.loads(replay_path.read_bytes())
    require(replay['status'] == 'PASS' and replay['round'] == 1324 and len(replay['checks']) == 2
            and replay['new_axioms'] == 0 and replay['stage_complete'] is False, 'Replay incomplete')
    scripts = {(FOLDER / n).relative_to(ROOT).as_posix() for n in ('check.py', 'independent_check.py')}
    require({c['script'] for c in replay['checks']} == scripts, 'Wrong replay scripts')
    for c in replay['checks']:
        require(c['exit_code'] == 0 and c['checks'] > 0 and bool(c['receipt']), 'Replay failed')
        require(sha((ROOT / c['script']).read_bytes()) == c['script_sha256'], 'Replay source changed')
    acceptance = dict(round=1324, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, goal_status='active',
        next_available_round=1325, next_round_admitted=False, scientific_evidence_files=10,
        executed_checks=replay['checks'], sha256=source_map(package),
        main_result='The same normalized reference parent yields a finite coherent effective conversion, full matter-reference-partner output error, exact source backreaction and a held-out probability distinction after one short-time calibration.',
        scope_limits=[
            'The enlarged finite balance window, coherent preparation and internal phase are declared inputs.',
            'The initial reference is independent of matter and its passive partner; correlated preparation apparatus is not certified.',
            'The complete K<=L+1 parent retains vacuum, transverse dark states and highest right occupation.',
            'A nearly unchanged comparison reference does not imply zero consumption or free reset.',
            'The bare linear menu is a comparison outside the 1323 balance contract, not an internal countermodel.',
            'Preparation response is not a Higgs functional source, a physical pole mass or a complete local field realization.',
            'All right species, family uniqueness and empirical parameters remain unselected.'])
    readme_path, ledger_path = STAGE / 'README.md', STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    old_readme = readme_raw.decode('utf8')
    require(old_readme.count('## 当前进展') == 1, 'Unexpected README boundary')
    readme = old_readme.split('## 当前进展', 1)[0] + '''## 当前进展

2026-10-19。验收至1324，共两轮新研究；新增公理0，目标保持active，尚未阶段结项。[联合进展](_shared/progress_after1324.md)明确A／B／C／J的已证范围与未决接口。

|轮次|新增结果|证据|
|---|---|---|
|[1323](research_note_1323.md)|全未知参考逐份交接选择单下降端口的右角色；短时余额一致性选择归一提取|[数学独审](1323/independent_review.md) · [范围独审](1323/scope_review.md) · [验收](1323/acceptance.json)|
|[1324](research_note_1324.md)|同一父核的有限相干有效转换、完整联合后态误差、精确资源反作用及一次校准后的留出区别|[数学独审](1324/independent_review.md) · [范围独审](1324/scope_review.md) · [验收](1324/acceptance.json)|

参考准备、相位和更大有限窗口是公开输入；同一参考乘各味的共同因子仍不选择现实矩阵。普通线性菜单是1323余额合同外的对照，模式准备响应不是Higgs场源，近似产品后态不代表参考没有消耗。

[下一步](1324/NEXT.md)先筛同一左端与相反参考流的共同父过程。1325为下一个可用科学编号，尚未准入；不再以更多参考数目扫描代替跨部门进展。

[本轮发布事务与前图像](_shared/progress_publication_after1324.json)保留1323后的导航边界；[1323发布记录](_shared/progress_publication_after1323.json)及[启动记录](_shared/navigation_start_receipt.json)保留原历史，未覆盖旧验收。
'''
    ledger = ledger_raw.decode('utf8')
    before_intro = '2026-10-19。下表保留启动时的输入基线；科学成果现验收至1323，最新判定见文末。'
    require(before_intro in ledger, 'Unexpected ledger boundary')
    ledger = ledger.replace(before_intro, before_intro.replace('至1323', '至1324'), 1)
    ledger += '''
## 1324后的更新

[报告](../research_note_1324.md)与[联合进展](progress_after1324.md)将1323同一父核接到实际有限相干参考。完整物质—参考—任意伙伴输出受统一误差控制，真实参考数目与能源另由精确父守恒核算；不能由近似产品后态领取无消耗或免费复位。

一次短时校准后固定系数与时钟，留出准备的实际概率可辨；普通线性菜单明确为余额合同外对照。所有味的有效转换乘同一准备因子，M与三／四／五味仍自由。这是模式准备响应，不是独立Higgs源或质量极点。

扩大窗口、相位准备及实际访问权限已经登记。B／J有有限共同接口，C的全部右端口及跨部门判定仍未完成；1325尚未准入，目标active，新公理0。[正式验收](../1324/acceptance.json)。
'''
    prefix = '''
## 当前目标：基本相互作用规则的认知来源

2026-10-19，[1323起阶段](archive_1323_/README.md)已验收1323—1324两轮，新增公理0，目标保持active。[正式目标](archive_1323_/_shared/goal.md) · [1324报告](archive_1323_/research_note_1324.md) · [联合进展](archive_1323_/_shared/progress_after1324.md)。

1323给出全未知参考逐份交接下的单端口角色和归一原语；1324由同一父作用得到有限相干有效转换、物质—参考—未知伙伴的完整误差及精确资源反作用。一次短时校准后的留出准备仍给可辨预测，保守概率间隔大于0.70。

相干准备、相位和更大有限窗口属于公开资源；近似产品参考不等于没有消耗。普通线性菜单是余额合同外的对照，准备响应不等于Higgs源或现实质量。全部右端口、跨部门关系及局域场字典仍须共同研究，三／四／五味并存。下一项先筛相反参考流共用左端，1325尚未准入，不据两轮局部连接结束阶段。

上一阶段仍保存在[archive_1313_1322](archive_1313_1322/README.md)与[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)。[迁移与历史复算](_migration/closure_1313_1322/README.md)及下方历史记录保留；当前安排以本节与现行目标为准。
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
    layer.update(latest_accepted_round=1324, in_progress_round=None, next_round=1325,
        next_round_admitted=False, stage_complete=False, application_goal_status_last_observed='active',
        startup_state='Round 1324 accepted with conditional scope; next round not admitted')
    layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger.encode('utf8')), (LAYER_PATH, old_raw, layer_raw)]
    planned = {p.resolve(): after for p, before, after in updates}
    planned[ACCEPTANCE.resolve()] = encoded(acceptance)
    fresh_paths = {CHECKS.resolve(), RECEIPT.resolve()}
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
          STAGE / '_shared/progress_after1324.md', goal]
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
    sources = source_map(package + [goal, replay_path, STAGE / '_shared/progress_after1324.md',
        STAGE / '_shared/start_stage.py', STAGE / '_shared/goal_start_receipt.json',
        ROOT / 'scripts/run_research_interaction_rules_active.py', Path(__file__).resolve(),
        prior_acceptance, prior_receipt, prior_checks,
        STAGE / '_shared/publish_progress_after1323.py'])
    checks = dict(status='PASS', round=1324, scientific_evidence_files=10,
        markdown_files=len(md), local_links=links, previous_round_bindings=len(prior['sha256']), original_navigation_recoverable=True,
        prior_navigation=prior_navigation, scientific_files_modified=False,
        stage_complete=False, next_round_admitted=False)
    receipt = dict(status='PUBLISHED', latest_accepted_round=1324, new_scientific_rounds=[1324],
        new_axioms=0, goal_status='active', stage_complete=False, next_available_round=1325,
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
