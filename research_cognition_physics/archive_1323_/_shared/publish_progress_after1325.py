"""Review-bound round-1325 publication; default --check is read-only."""
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

FOLDER = STAGE / '1325'
ACCEPTANCE = FOLDER / 'acceptance.json'
CHECKS = STAGE / '_shared/progress_checks_after1325.json'
RECEIPT = STAGE / '_shared/progress_publication_after1325.json'
EXPECTED = {'1325/proof.md': '30e291fcd272c92b05cf8786175a1f60aa532485aa96b468c30d3cf22d00634a', '1325/check.py': 'f3df09ca17378cd9c61719ff6ac04c6e595343acaa6fc4082eb27cd2ca58aaf3', '1325/results.json': 'ba2cc3f14d3ed704dd519888098ddd463617bad6786be8cda1ec43b30b5a0a7d', '1325/NEXT.md': 'c59094bf9b498601c84ed90e7b6ba2ad493b67fdb588e25f998560359cc1fc18', 'research_note_1325.md': 'a95e45f9c4774a5c7cbe1b5d75eafd81b4117ea8167f297cb4253c9d46e134a1', '_admission/round1325_admission.md': '5be47f2cd0078c5dd666499a0f2d41211ca1a1354e748134ee4a0a88b75c36a4', '1325/independent_check.py': '4d31ee489fd2344883fab072b1be19491d7033ab210af0a8394d94d9244bcc12', '1325/independent_results.json': '84cbd9a3900fc91e315eb7503175b0c3fa41fa4db4682acc2ed5618b1b97ccd0', '1325/independent_review.md': '764543ee91ab4d790859df50dafc0d436fd91949fccf61c4b2edc6d9065c1311', '1325/scope_review.md': '9a15646d3312029292cdc06def8a8f3901bdbacb70a8e8ee2a92e4a1c762179d'}


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
    require(old['latest_accepted_round'] == 1324 and old['in_progress_round'] in (None, 1325)
            and old['next_round'] == 1325 and not old['stage_complete'], 'Scientific boundary differs')
    goal = STAGE / '_shared/goal.md'
    goal_receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_bytes())
    require(goal_receipt['status'] == 'active' and
            sha(goal.read_bytes()) == goal_receipt['goal_document_sha256'], 'Reviewed goal changed')
    observation_path = STAGE / '_shared/goal_status_observation_after1325.json'
    observation = json.loads(observation_path.read_bytes())
    require(observation['application_goal_status_observed'] == 'blocked'
            and observation['objective_document_changed'] is False
            and observation['stage_scientifically_complete'] is False, 'Goal observation changed')
    # Startup receipt records historical activation, not current application state.
    package = [STAGE / name for name in EXPECTED]
    require(len(package) == 10, 'Unexpected evidence count')
    for name, digest in EXPECTED.items():
        require(sha((STAGE / name).read_bytes()) == digest, 'Frozen source changed: ' + name)
    prior_acceptance = STAGE / '1324/acceptance.json'
    prior_receipt = STAGE / '_shared/progress_publication_after1324.json'
    prior_checks = STAGE / '_shared/progress_checks_after1324.json'
    prior = json.loads(prior_acceptance.read_bytes())
    require(prior['round'] == 1324 and prior['status'] == 'ACCEPTED_CONDITIONAL_SCOPE',
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
    replay_path = STAGE / '_shared/replay_after1325.json'
    replay = json.loads(replay_path.read_bytes())
    require(replay['status'] == 'PASS' and replay['round'] == 1325 and len(replay['checks']) == 2
            and replay['new_axioms'] == 0 and replay['stage_complete'] is False, 'Replay incomplete')
    scripts = {(FOLDER / n).relative_to(ROOT).as_posix() for n in ('check.py', 'independent_check.py')}
    require({c['script'] for c in replay['checks']} == scripts, 'Wrong replay scripts')
    for c in replay['checks']:
        require(c['exit_code'] == 0 and c['checks'] > 0 and bool(c['receipt']), 'Replay failed')
        require(sha((ROOT / c['script']).read_bytes()) == c['script_sha256'], 'Replay source changed')
    acceptance = dict(round=1325, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, research_continues=True,
        application_goal_status_last_observed=observation['application_goal_status_observed'],
        next_available_round=1326, next_round_admitted=False, scientific_evidence_files=10,
        executed_checks=replay['checks'], sha256=source_map(package),
        main_result='The same full dual-port reference parent gives conditional right-role selection, optional multiplicity matching, actual finite cross-port coherence and a nonredundant shared-reference curvature constraint on isospectral relative flavor orientation, with common exact handoff and complete sources.',
        scope_limits=['The K>=4 equivariant role conclusion requires complete unknown reference/port/partner inputs, not separately tested ports.', 'Trivial copies remain family multiplicities; equal multiplicities follow only with the optional full-coverage contract.', 'The coherent product reference task has different initial quantifiers from the full unknown fixed-K task.', 'The complete K<=L+2 parent retains vacuum, low boundaries, transverse branches and highest occupation.', 'The shared Gram budget is necessary and sufficient only for the declared initial reference-leakage curvature task; Duhamel is not reversed.', 'The illustrative budget is an input, not a natural constant or empirical mixing prediction; both isospectral examples meet the original common-clock handoff.', 'Each candidate keeps its own complete sources; identical source values are not assumed.', 'This finite single-fermion reference model is not a local Higgs field, a physical mass-pole reconstruction or a complete standard-model realization.'])
    readme_path, ledger_path = STAGE / 'README.md', STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    old_readme = readme_raw.decode('utf8')
    require(old_readme.count('## 当前进展') == 1, 'Unexpected README boundary')
    readme = old_readme.split('## 当前进展', 1)[0] + '## 当前进展\n\n2026-10-19。验收至1325，共三轮新研究；新增公理0，阶段研究继续，尚未结项。[联合进展](_shared/progress_after1325.md)区分科学范围和应用运行状态；[状态观测](_shared/goal_status_observation_after1325.json)记录最近工具状态为blocked，不据此宣称存在科学阻碍。\n\n|轮次|新增结果|证据|\n|---|---|---|\n|[1323](research_note_1323.md)|全未知参考逐份交接选择单下降端口角色；短时余额一致性选择归一提取|[数学独审](1323/independent_review.md) · [范围独审](1323/scope_review.md) · [验收](1323/acceptance.json)|\n|[1324](research_note_1324.md)|同父有限相干有效转换、完整后态误差、真实资源反作用及留出预测|[数学独审](1324/independent_review.md) · [范围独审](1324/scope_review.md) · [验收](1324/acceptance.json)|\n|[1325](research_note_1325.md)|相反参考流共用左端：联合角色、可选重数匹配、实际相干恢复及非冗余跨部门预算|[数学独审](1325/independent_review.md) · [范围独审](1325/scope_review.md) · [验收](1325/acceptance.json)|\n\n含K=4的全输入合同选择两右平凡角色，平凡拷贝仍保为重数；可选满覆盖匹配部门，不选三代。共享参考预算对同谱相对味方向有实际选择力，但预算与相位准备均是公开输入，不预测现实混合角。\n\n[下一步](1325/NEXT.md)和[独立阶段方向审计](_shared/joint_direction_screen_after1325.md)先判断有限参考与物理来源的身份。[驱动来源候选](_admission/driven_reference_source_screen_after1325.md)仍是计0筛查，1326尚未准入；不继续扫描参考余额，也不把全UV或全部装置制造临时设为门槛。\n\n[本轮发布事务与前图像](_shared/progress_publication_after1325.json)保留1324后的导航边界；[1324发布记录](_shared/progress_publication_after1324.json)、[1323发布记录](_shared/progress_publication_after1323.json)及[启动记录](_shared/navigation_start_receipt.json)不变。\n'
    ledger = ledger_raw.decode('utf8')
    before_intro = '2026-10-19。下表保留启动时的输入基线；科学成果现验收至1324，最新判定见文末。'
    require(before_intro in ledger, 'Unexpected ledger boundary')
    ledger = ledger.replace(before_intro, before_intro.replace('至1324', '至1325'), 1)
    ledger += '\n## 1325后的更新\n\n[报告](../research_note_1325.md)与[联合进展](progress_after1325.md)把两种相反参考流放进同一完整父作用。含K=4的全未知等变合同选择两右平凡角色；额外满覆盖才给三部门重数匹配。全未知固定K任务和已知相干参考任务的初态量词保持区分，三／四／五味继续允许。\n\n真实未后选相干读口κ约0.984375；共同参考初始退化预算对同谱相对味方向有非冗余约束。两比较分支亦满足同一奇数谱时钟的全未知交付；预算β是公开示意输入，不能称自然常数或经验混合预测。完整参考、能源及荷账来自同一个父核，不从近似产品后态领取零消耗。\n\n[独立方向审计](joint_direction_screen_after1325.md)及[驱动来源筛查](../_admission/driven_reference_source_screen_after1325.md)界定下一有限物理身份问题：准备响应尚不是Higgs源或现实质量极点。1326未准入，新增公理0，阶段研究继续且尚未结项。[正式验收](../1325/acceptance.json)。\n\n上述1323／1324段的active是当时记录；[最近应用状态观测](goal_status_observation_after1325.json)为blocked，未据此改变目标内容或宣布科学阻碍。\n'
    prefix = '\n## 当前目标：基本相互作用规则的认知来源\n\n2026-10-19，[1323起阶段](archive_1323_/README.md)已验收1323—1325三轮，新增公理0，阶段研究继续、尚未结项。[正式目标](archive_1323_/_shared/goal.md) · [1325报告](archive_1323_/research_note_1325.md) · [联合进展](archive_1323_/_shared/progress_after1325.md)。\n\n同一归一参考原语现连接全未知交接、有限相干有效转换与相反流双端口。1325在不预设共同味数下判定右角色，可选满覆盖匹配部门；实际共享参考预算还能区别同谱的相对味方向，两分支同时保留同钟全未知交付。\n\n相位准备、覆盖及预算是公开前件，三／四／五味和现实参数仍有自由。准备响应尚不是Higgs源或现实质量。下一步先做[有限物理来源筛查](archive_1323_/_admission/driven_reference_source_screen_after1325.md)，1326未准入；[阶段方向审计](archive_1323_/_shared/joint_direction_screen_after1325.md)避免将全UV或所有设备制造临时加为门槛。应用目标最近工具状态另见[状态观测](archive_1323_/_shared/goal_status_observation_after1325.json)，不沿用旧active标签或据此宣称科学失败。\n\n上一阶段仍保存在[archive_1313_1322](archive_1313_1322/README.md)与[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)。[迁移与历史复算](_migration/closure_1313_1322/README.md)及下方历史记录保留；当前安排以本节与现行目标为准。\n'
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
    layer.update(latest_accepted_round=1325, in_progress_round=None, next_round=1326,
        next_round_admitted=False, stage_complete=False, application_goal_status_last_observed=observation['application_goal_status_observed'],
        research_continues=True,
        startup_state='Round 1325 accepted with conditional scope; next round not admitted')
    layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger.encode('utf8')), (LAYER_PATH, old_raw, layer_raw)]
    planned = {p.resolve(): after for p, before, after in updates}
    planned[ACCEPTANCE.resolve()] = encoded(acceptance)
    fresh_paths = {CHECKS.resolve(), RECEIPT.resolve()}
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
          STAGE / '_shared/progress_after1325.md', goal,
          STAGE / '_shared/joint_direction_screen_after1325.md',
          STAGE / '_admission/driven_reference_source_screen_after1325.md']
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
    sources = source_map(package + [goal, replay_path, STAGE / '_shared/progress_after1325.md',
        observation_path, STAGE / '_shared/joint_direction_screen_after1325.md',
        STAGE / '_admission/driven_reference_source_screen_after1325.md',
        STAGE / '_shared/start_stage.py', STAGE / '_shared/goal_start_receipt.json',
        ROOT / 'scripts/run_research_interaction_rules_active.py', Path(__file__).resolve(),
        prior_acceptance, prior_receipt, prior_checks,
        STAGE / '_shared/publish_progress_after1324.py'])
    checks = dict(status='PASS', round=1325, scientific_evidence_files=10,
        markdown_files=len(md), local_links=links, previous_round_bindings=len(prior['sha256']), original_navigation_recoverable=True,
        prior_navigation=prior_navigation, scientific_files_modified=False,
        stage_complete=False, next_round_admitted=False)
    receipt = dict(status='PUBLISHED', latest_accepted_round=1325, new_scientific_rounds=[1325],
        new_axioms=0, research_continues=True,
        application_goal_status_last_observed=observation['application_goal_status_observed'],
        stage_complete=False, next_available_round=1326,
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
