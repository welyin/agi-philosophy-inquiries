"""Finalize the reviewed finite conditional stage once; default check is read-only."""
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

ACCEPTANCE = STAGE / '_shared/stage_acceptance.json'
CHECKS = STAGE / '_shared/stage_completion_checks.json'
RECEIPT = STAGE / '_shared/stage_completion_publication.json'
PAPER = STAGE.parent / '16_未知参考交接的作用选择与动态来源_阶段论文.md'
GOAL = STAGE / '_shared/goal.md'
LEDGER = STAGE / '_shared/joint_premise_ledger_after1326.md'
MATH = STAGE / '_shared/stage_math_review.md'
SCOPE = STAGE / '_shared/stage_scope_review.md'
JOINT_CODE = STAGE / '_shared/check_stage_joint.py'
JOINT = STAGE / '_shared/stage_joint_results.json'
REPLAY = STAGE / '_shared/stage_joint_replay.json'
EXPECTED = {'research_cognition_physics/16_未知参考交接的作用选择与动态来源_阶段论文.md': '964ee3b135be50e342173687d5b13ac5c8131b2c745511a320260b9ea0f4d492', 'research_cognition_physics/archive_1323_/_shared/goal.md': '95bbeb961464c39cdd4e3867d6f00faec868f38fd5cb2cf43c59d4393b7e3eb9', 'research_cognition_physics/archive_1323_/_shared/joint_premise_ledger_after1326.md': '4f7c3fc031ad0cb9c1409e5ad96f338d3d0c6b253727958604948f92b230a257', 'research_cognition_physics/archive_1323_/1323/acceptance.json': 'cde6448206fe24e18bf26bdc28bb366a8ea491deeae1d01437a3f52ee7776120', 'research_cognition_physics/archive_1323_/1324/acceptance.json': '72cefd1ebfe12265f6019c02bc95fcf9471d4a6d1d5336235190c2a9b221caa7', 'research_cognition_physics/archive_1323_/1325/acceptance.json': '296d65941fbdd365dac405a792068d7a19b350886406d1a538395d943bb34e8e', 'research_cognition_physics/archive_1323_/1326/acceptance.json': 'dcfb30abcba21f455f67ffbb74e746726137aa71c7dbf6e412b6651eebb34304', 'research_cognition_physics/archive_1323_/_shared/check_stage_joint.py': 'e80b2c8a252c0c668f5b625f47dd8bd9d87039f679b30bd4d0add0655aff99dd', 'research_cognition_physics/archive_1323_/_shared/stage_joint_results.json': 'fae12d154efd41c8c7923b60e22d473aef67e3480de285e51e498ddff47dde8d', 'research_cognition_physics/archive_1323_/_shared/stage_joint_replay.json': '267aadcb645b865455ea984dd72cc2f6b71c016c95bcb8b9480bebd01b4efeb1', 'research_cognition_physics/archive_1323_/_shared/stage_math_review.md': '98208c229727301c9f72227f1a8d6b82be6ce0e12777b949a9a0f83a42aef965', 'research_cognition_physics/archive_1323_/_shared/stage_scope_review.md': '521fb310dbf466a70397e0131eca294212c2d2f206b06df27ab734892edf2406'}
CRITERIA = {
    'C1': 'operation definitions and dependency ledger',
    'C2': 'declared interaction class and composition/effective transmission',
    'C3': 'fair matter-role comparison and cross-department parameter consequence',
    'C4': 'same-parent nonempty settings, complete sources and finite predictions',
    'C5': 'paper, reproducible evidence and independent joint acceptance',
}


def encoded(value):
    return (json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf8')


def sources(paths):
    return {p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in paths}


def prepare():
    for path in (ACCEPTANCE, CHECKS, RECEIPT):
        require(not path.exists(), 'Stage already finalized or partial transaction: ' + str(path))
    store = InteractionRulesLayout()
    prior_navigation = store.verify_interaction_rules_navigation()
    old_raw = LAYER_PATH.read_bytes()
    old = json.loads(old_raw)
    require(old['latest_accepted_round'] == 1326 and old['in_progress_round'] is None
            and old['next_round'] == 1327 and not old['next_round_admitted']
            and not old['stage_complete'], 'Unexpected stage boundary')
    latest_pub = STAGE / '_shared/progress_publication_after1326.json'
    latest_check = STAGE / '_shared/progress_checks_after1326.json'
    publication = json.loads(latest_pub.read_bytes())
    require(publication['status'] == 'PUBLISHED' and
            sha(old_raw) == publication['current_layer_sha256'], 'Exact 1326 boundary changed')
    require(sha(latest_check.read_bytes()) == publication['checks_sha256'], '1326 checks changed')
    for name, digest in publication['source_sha256'].items():
        require(sha((ROOT / name).read_bytes()) == digest, '1326 publication source changed: ' + name)
    for name, digest in EXPECTED.items():
        require(sha((ROOT / name).read_bytes()) == digest, 'Frozen stage evidence changed: ' + name)
    goal_receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_bytes())
    require(goal_receipt['goal_document_sha256'] == sha(GOAL.read_bytes()), 'Original goal changed')
    observation_path = STAGE / '_shared/goal_resumed_after1325.json'
    observation = json.loads(observation_path.read_bytes())
    require(observation['status'] == 'active' and observation['objective_changed'] is False,
            'Last recorded application observation changed')

    science = {}
    round_receipts = []
    accept_paths = []
    for n in range(1323, 1327):
        path = STAGE / str(n) / 'acceptance.json'
        obj = json.loads(path.read_bytes())
        require(obj['round'] == n and obj['status'] == 'ACCEPTED_CONDITIONAL_SCOPE'
                and obj['new_axioms'] == 0 and obj['scientific_evidence_files'] == 10,
                'Round not independently accepted')
        expected_paths = {(STAGE / str(n) / x).relative_to(ROOT).as_posix() for x in
                          ('proof.md', 'check.py', 'results.json', 'NEXT.md',
                           'independent_check.py', 'independent_results.json',
                           'independent_review.md', 'scope_review.md')}
        expected_paths |= {(STAGE / f'research_note_{n}.md').relative_to(ROOT).as_posix(),
                           (STAGE / f'_admission/round{n}_admission.md').relative_to(ROOT).as_posix()}
        require(set(obj['sha256']) == expected_paths, 'Unexpected scientific file set')
        for name, digest in obj['sha256'].items():
            require(sha((ROOT / name).read_bytes()) == digest, 'Frozen science changed: ' + name)
            require(name not in science, 'Repeated evidence path')
            science[name] = digest
        require(len(obj['executed_checks']) == 2, 'Missing round replay')
        for r in obj['executed_checks']:
            require(r['exit_code'] == 0 and r['checks'] > 0 and r['receipt'], 'Failed round replay')
            require(sha((ROOT / r['script']).read_bytes()) == r['script_sha256'], 'Replay script changed')
            round_receipts.append(dict(round=n, **r))
        accept_paths.append(path)
    require(len(science) == 40, 'Expected forty frozen scientific files')
    require(sha(accept_paths[-1].read_bytes()) == publication['acceptance_sha256'],
            '1326 acceptance differs from publication')

    math = MATH.read_text(encoding='utf8')
    scope = SCOPE.read_text(encoding='utf8')
    require('PASS_WITH_SCOPE' in math and 'PASS_WITH_SCOPE' in scope, 'Final stage reviews required')
    for path in (PAPER, GOAL, LEDGER, JOINT_CODE, JOINT, *accept_paths):
        digest = sha(path.read_bytes())
        require(digest in math, 'Stage mathematical review lacks final binding: ' + str(path))
        require(digest in scope, 'Stage scope review lacks final binding: ' + str(path))
    require(sha(MATH.read_bytes()) in scope, 'Scope review does not bind final mathematical review')
    for criterion in CRITERIA:
        require('| ' + criterion + ' | PASS |' in scope, 'Missing original criterion: ' + criterion)

    joint = json.loads(JOINT.read_bytes())
    require(joint['status'] == 'PASS_STAGE_AUXILIARY_CHECKS_ONLY'
            and joint['new_scientific_rounds'] == 0 and joint['new_axioms'] == 0
            and joint['scientific_files_verified'] == 40, 'Joint auxiliary evidence incomplete')
    require(joint['goal_sha256'] == sha(GOAL.read_bytes())
            and joint['paper_checks']['sha256'] == sha(PAPER.read_bytes()), 'Joint document mismatch')
    for name, digest in joint['source_sha256'].items():
        require(sha((ROOT / name).read_bytes()) == digest, 'Joint source changed: ' + name)
    require([r['acceptance_sha256'] for r in joint['round_evidence']]
            == [sha(p.read_bytes()) for p in accept_paths], 'Joint round receipts differ')
    replay = json.loads(REPLAY.read_bytes())
    require(replay['status'] == joint['status'] and replay['exit_code'] == 0 and replay['receipt']
            and replay['scientific_files_verified'] == 40
            and replay['analytic_arithmetic_checks'] == len(joint['analytic_arithmetic_checks']),
            'Root joint replay incomplete')
    require(replay['script_sha256'] == sha(JOINT_CODE.read_bytes())
            and replay['result_sha256'] == sha(JOINT.read_bytes()), 'Joint replay binding changed')

    completion = {key: dict(status='PASS', original_requirement=value) for key, value in CRITERIA.items()}
    acceptance = dict(schema='interaction_rules_conditional_stage_acceptance_v1',
        date='2026-10-19', status='COMPLETE_CONDITIONAL_SCOPE', stage_complete=True,
        scientific_rounds=list(range(1323, 1327)), new_axioms=0,
        original_goal_unchanged=True, completion_criteria=completion,
        accepted_artifact_bindings=40, scientific_sha256=science,
        executed_round_checks=round_receipts,
        joint_transport=dict(script=JOINT_CODE.relative_to(ROOT).as_posix(),
            results=JOINT.relative_to(ROOT).as_posix(), root_replay=replay,
            new_scientific_rounds=0, no_scientific_rescan=True),
        application_goal_status_last_observed=observation['status'],
        application_goal_completion_recorded=False,
        application_state_scope='Scientific stage completion only; the application goal tool is handled separately and is not predicted here.',
        next_available_round=1327, next_round_admitted=False,
        scope_limits=[
            'Declared SU(2) finite two-mode reference, single-fermion unit-transfer and pure-handoff class.',
            'Full unknown-reference delivery, coherent product-reference evolution and driven low bands are distinct settings of one parameterized parent family.',
            'Normalized balance response is derived only under the added open-short-time task; it is not a new adopted cognitive axiom.',
            'Optional full coverage matches multiplicities; optional reference-curvature budget constrains relative flavor orientation; neither selects three generations or empirical parameters.',
            'The standard linear menu receives no exception to the balance contract, and the 1322 vectorlike branch is excluded only from the specified new interface.',
            'Finite drive source, energy kernel and bounded readout permissions do not identify a complete local Higgs/gauge field or a real mass-pole theory.',
            'No unlimited-energy limit, whole standard model, autonomous controller or theory of everything is claimed.'
        ],
        sha256=sources([PAPER, GOAL, LEDGER, MATH, SCOPE, JOINT_CODE, JOINT, REPLAY, *accept_paths]))

    readme_path = STAGE / 'README.md'
    selection_path = STAGE / '_shared/selection_ledger.md'
    readme_raw, selection_raw = readme_path.read_bytes(), selection_path.read_bytes()
    readme_text = readme_raw.decode('utf8')
    require(readme_text.count('## 当前进展') == 1, 'Unexpected README structure')
    readme = readme_text.split('## 当前进展', 1)[0] + '## 当前进展\n\n2026-10-19。**1323—1326阶段按原五项标准完成条件性联合验收。** 新增公理0；本阶段科学结项；后续场身份研究尚未启动。目录暂保archive_1323_，1327未准入。应用目标工具的最终状态另行登记，本页不预报其变更。\n\n[论文16](../16_未知参考交接的作用选择与动态来源_阶段论文.md) · [阶段验收](_shared/stage_acceptance.json) · [数学联合审阅](_shared/stage_math_review.md) · [范围联合审阅](_shared/stage_scope_review.md) · [联合复核结果](_shared/stage_joint_results.json) · [前件与剩余自由](_shared/joint_premise_ledger_after1326.md)。\n\n|轮次|新增结果|证据|\n|---|---|---|\n|[1323](research_note_1323.md)|全未知参考交付的角色与归一提取、实际顺序组合|[验收](1323/acceptance.json)|\n|[1324](research_note_1324.md)|有限相干转换、完整后态误差、参考反作用与留出后果|[验收](1324/acceptance.json)|\n|[1325](research_note_1325.md)|双端口角色、可选重数匹配与非冗余跨部门预算|[验收](1325/acceptance.json)|\n|[1326](research_note_1326.md)|动态参考全部低带、实际来源、测量工作及固定校准后的响应|[验收](1326/acceptance.json)|\n\n同一参数化父族连接全未知交付、相干产品任务和驱动低带，分别保各自初态与设置。论文第6节给固定门户点和三味合法缩放，不把不同数值表拼接。作用选择、公平角色竞争、跨部门关系和真实有限后果满足原目标；不是完整标准模型或现实参数的终结。\n\n有限模式SU(2)、驱动f及继承Dirac能量核尚未自动取得完整局域弱场、Higgs或现实质量极点身份。任意有限重数、预算、谱与部分作用类别仍有输入。下一方向应先判断归一规则的局域场身份或设备接口身份；未自动采用新公理、创建新目标或准入新轮。\n\n[阶段收尾检查](_shared/stage_completion_checks.json)与[可逆发布记录](_shared/stage_completion_publication.json)保留1326发布后的导航前图像。四轮科学包40件及旧支持文件保持不变；[1326发布记录](_shared/progress_publication_after1326.json)和各历史轮次均保留。\n'
    selection = selection_raw.decode('utf8')
    old_intro = '2026-10-19。下表保留启动时的输入基线；科学成果现验收至1326，最新判定见文末。'
    require(old_intro in selection and '## 1323—1326阶段联合结项' not in selection,
            'Unexpected selection-ledger boundary')
    selection = selection.replace(old_intro,
        '2026-10-19。下表与逐轮更新保留各历史时点；1323—1326阶段已按原五项标准联合验收，最终范围见文末。', 1)
    selection += '\n## 1323—1326阶段联合结项\n\n原五项要求已按[论文16](../../16_未知参考交接的作用选择与动态来源_阶段论文.md)、[数学审阅](stage_math_review.md)、[范围审阅](stage_scope_review.md)及[阶段验收](stage_acceptance.json)逐项确认，新增公理0。四轮40件科学材料、完整来源与20项联合运输算术均核验，科学阶段完成，后续场身份研究尚未启动。\n\n这是指定未知参考访问／余额合同和有限有效类别下的条件性选择；标准线性菜单也接受同一合同检验。右角色选择只针对新接口，三／四／五及任意有限共同重数保留；β预算及经验矩阵没有被生成。低带f源不是Higgs字段身份，完整局域场连接仍是后续方向。\n\n共同父族通过固定门户和明确缩放承接不同任务设置，未把互斥准备合并为一次实验。[联合前件账](joint_premise_ledger_after1326.md)和原各轮记录保持冻结。1327未准入，未自动创建下一目标；应用工具状态由主线另行记录。本节以前的active/blocked/待结项均是各历史时点记录。\n'
    prefix = '\n## 已完成阶段：未知参考交接的作用选择与动态来源\n\n2026-10-19，1323—1326四轮已按原五项要求完成条件性联合验收，新增公理0。[论文16](16_未知参考交接的作用选择与动态来源_阶段论文.md) · [阶段目录](archive_1323_/README.md) · [联合验收](archive_1323_/_shared/stage_acceptance.json) · [前件与剩余自由](archive_1323_/_shared/joint_premise_ledger_after1326.md)。\n\n声明的全未知参考合同选择活动右角色与归一原语，真实组合及双端口给跨部门约束；同一父族进一步连到有限动态参考的完整低带、真实来源与留出预测。共同参数与设置已核，原四轮科学包保持冻结。\n\n结项范围是有限、条件性的作用与接口选择；不表示完整标准模型、三代或现实参数已推出。有限参考的局域场身份仍需后续检验，不能把f改名Higgs或把SU(2)模式直接认成全规范场。本阶段科学结项，后续场身份研究尚未启动。目录暂保archive_1323_，1327未准入；应用目标工具状态另行登记，没有自动开启新阶段。\n\n上一阶段仍在[archive_1313_1322](archive_1313_1322/README.md)及[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)。[迁移与历史复算](_migration/closure_1313_1322/README.md)和下方历史时点保留；当前科学结论以本节及阶段验收为准。\n'
    layer = copy.deepcopy(old)
    updates = []
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        before = path.read_bytes()
        require(sha(before) == entry['current_sha256'], 'Concurrent root navigation edit')
        original = relocation_original_bytes(before, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        text, at = original.decode(enc), entry['edits'][0]['start']
        require(at == text.index('\n') + 1 and len(entry['edits']) == 1
                and entry['edits'][0]['before'] == '', 'Unexpected reversible prefix')
        after = (text[:at] + prefix + text[at:]).encode(enc)
        entry['edits'][0]['after'] = prefix
        entry.update(current_sha256=sha(after), current_bytes=len(after))
        require(relocation_original_bytes(after, entry) == original, 'Prefix not reversible')
        updates.append((path, before, after))
    layer.update(stage_complete=True, in_progress_round=None, next_round=1327,
        next_round_admitted=False, application_goal_status_last_observed=observation['status'],
        research_continues=False,
        current_work='Conditional finite stage scientifically complete; application goal completion is recorded separately.',
        startup_state='Rounds 1323-1326 and all original stage criteria independently accepted')
    layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (selection_path, selection_raw, selection.encode('utf8')),
                (LAYER_PATH, old_raw, layer_raw)]
    planned = {p.resolve(): after for p, before, after in updates}
    planned[ACCEPTANCE.resolve()] = encoded(acceptance)
    forward = {CHECKS.resolve(), RECEIPT.resolve()}
    markdown = [ROOT / name for name in science if name.endswith('.md')]
    markdown += [PAPER, GOAL, LEDGER, MATH, SCOPE, readme_path, selection_path]
    links = 0
    for path in markdown:
        text = planned.get(path.resolve(), path.read_bytes()).decode('utf8')
        require(len(re.findall(r'^\s*\$\$\s*$', text, re.M)) % 2 == 0, 'Unbalanced block math: ' + str(path))
        require(not any(line.strip() in (chr(92)+'[', chr(92)+']') for line in text.splitlines()), 'Legacy block delimiters')
        require(not any(ord(c) < 32 and c not in '\r\n\t' for c in text), 'Control character')
        for target in re.findall(r'\]\(([^)\s]+)\)', text):
            if '://' in target or target.startswith('#'):
                continue
            linked = (path.parent / unquote(target.split('#')[0])).resolve()
            require(linked.is_file() or linked in planned or linked in forward,
                    'Missing local link: ' + str(path) + ' -> ' + target)
            links += 1

    bound = sources([PAPER, GOAL, LEDGER, MATH, SCOPE, JOINT_CODE, JOINT, REPLAY,
        *accept_paths, observation_path, latest_pub, latest_check, Path(__file__).resolve(),
        STAGE / '_shared/start_stage.py', STAGE / '_shared/goal_start_receipt.json',
        ROOT / 'scripts/run_research_interaction_rules_active.py'])
    bound.update(science)
    checks = dict(status='PASS', scientific_stage_complete=True, new_axioms=0,
        completion_criteria=completion, accepted_scientific_files=40,
        scientific_files_modified=False, original_goal_unchanged=True,
        markdown_files=len(markdown), local_links=links,
        analytic_arithmetic_checks=len(joint['analytic_arithmetic_checks']),
        previous_navigation=prior_navigation, original_navigation_recoverable=True,
        next_round_admitted=False, application_goal_completion_recorded=False)
    receipt = dict(schema='interaction_rules_stage_completion_publication_v1',
        date='2026-10-19', status='PUBLISHED', scientific_stage_complete=True,
        scientific_rounds=list(range(1323, 1327)), new_axioms=0,
        application_goal_status_last_observed=observation['status'],
        application_goal_completion_recorded=False,
        previous_layer_sha256=sha(old_raw), current_layer_sha256=sha(layer_raw),
        previous_layer=old, stage_acceptance_sha256=sha(encoded(acceptance)),
        source_sha256=bound,
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): before.decode('utf8')
            for p, before in ((readme_path, readme_raw), (selection_path, selection_raw))})
    return dict(updates=updates, acceptance=encoded(acceptance), checks=checks,
                receipt=receipt, sources=bound)


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
            require(sha((ROOT / name).read_bytes()) == digest, 'Source changed during finalization: ' + name)
        plan['checks']['completed_navigation'] = InteractionRulesLayout().verify_interaction_rules_navigation()
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
    group.add_argument('--publish', action='store_true', help='Finalize reviewed stage once')
    args = parser.parse_args()
    if args.publish:
        with publication_lock():
            publish_locked()
    else:
        plan = prepare()
        print(json.dumps(dict(plan['checks'], status='READY', files_written=0), ensure_ascii=False))


if __name__ == '__main__':
    main()
