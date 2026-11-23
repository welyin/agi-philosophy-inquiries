"""Review-bound round-1326 publication; default --check is read-only."""
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

FOLDER = STAGE / '1326'
ACCEPTANCE = FOLDER / 'acceptance.json'
CHECKS = STAGE / '_shared/progress_checks_after1326.json'
RECEIPT = STAGE / '_shared/progress_publication_after1326.json'
EXPECTED = {'1326/proof.md': '135283e3dd25d9f0dfb87d4f9e04f2e8c4c146634fd16c73085e05bacc1f3c03', '1326/check.py': '1c9dcbb8adfe622a93dde3c23a26661921456e5c3ab07868f23506ff0a9643eb', '1326/results.json': 'afac159e39c6be579dd82e88201ae888a440ed26c6ee7ccf43683c2f1d471309', '1326/NEXT.md': 'ce4f0097d959f0f61841633fe5f9efa7495701719b0a178a27baf2dc6295f5a6', 'research_note_1326.md': '88f8fb9d1d6773b076cf97ad69b32e7a7a1c631b56d49188c0453f6f3a645329', '_admission/round1326_admission.md': '18bb6035737f8e4b951f4156058b2fee407092cc044d48261c3de58572a4f5ae', '1326/independent_check.py': 'e2d6ee26fbbe9081ea2e253be45ae494becee3e117f253ca5b9e9669be6490fc', '1326/independent_results.json': '74c42e94725fe6067e2665d6b514193d0bd865118900f9cd6f3e95fd71485bfb', '1326/independent_review.md': 'dccb075773ef70cd514c1a6d76ead07e4927ff80ed0a3c99ed95283470ecc10e', '1326/scope_review.md': '3a88c84ac09ce23fc364ae11a7318ecfdb443b66e82a9b42326b08ab08733b4c'}


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
    require(old['latest_accepted_round'] == 1325 and old['in_progress_round'] in (None, 1326)
            and old['next_round'] == 1326 and not old['stage_complete'], 'Scientific boundary differs')
    goal = STAGE / '_shared/goal.md'
    goal_receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_bytes())
    require(goal_receipt['status'] == 'active' and
            sha(goal.read_bytes()) == goal_receipt['goal_document_sha256'], 'Reviewed goal changed')
    observation_path = STAGE / '_shared/goal_resumed_after1325.json'
    observation = json.loads(observation_path.read_bytes())
    require(observation['status'] == 'active'
            and observation['objective_changed'] is False
            and observation['current_scientific_blocker'] is None, 'Goal resumption changed')
    # Startup receipt records historical activation, not current application state.
    package = [STAGE / name for name in EXPECTED]
    require(len(package) == 10, 'Unexpected evidence count')
    for name, digest in EXPECTED.items():
        require(sha((STAGE / name).read_bytes()) == digest, 'Frozen source changed: ' + name)
    prior_acceptance = STAGE / '1325/acceptance.json'
    prior_receipt = STAGE / '_shared/progress_publication_after1325.json'
    prior_checks = STAGE / '_shared/progress_checks_after1325.json'
    prior = json.loads(prior_acceptance.read_bytes())
    require(prior['round'] == 1325 and prior['status'] == 'ACCEPTED_CONDITIONAL_SCOPE',
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
    replay_path = STAGE / '_shared/replay_after1326.json'
    replay = json.loads(replay_path.read_bytes())
    require(replay['status'] == 'PASS' and replay['round'] == 1326 and len(replay['checks']) == 2
            and replay['new_axioms'] == 0 and replay['stage_complete'] is False, 'Replay incomplete')
    scripts = {(FOLDER / n).relative_to(ROOT).as_posix() for n in ('check.py', 'independent_check.py')}
    require({c['script'] for c in replay['checks']} == scripts, 'Wrong replay scripts')
    for c in replay['checks']:
        require(c['exit_code'] == 0 and c['checks'] > 0 and bool(c['receipt']), 'Replay failed')
        require(sha((ROOT / c['script']).read_bytes()) == c['script_sha256'], 'Replay source changed')
    acceptance = dict(round=1326, status='ACCEPTED_CONDITIONAL_SCOPE', date='2026-10-19',
        new_axioms=0, stage_complete=False, research_continues=True,
        application_goal_status_last_observed=observation['status'],
        next_available_round=1327, next_round_admitted=False, scientific_evidence_files=10,
        executed_checks=replay['checks'], sha256=source_map(package),
        main_result='The same normalized interaction with an explicit finite dynamic reference gives the complete isolated low band, partner-preserving encoding, actual fixed-parameter drive sources and backreaction, and separated held-out response intervals after one calibration.',
        scope_limits=['The reference free kernel, drive, Dirac momentum dictionary and finite control/readout permissions are declared inputs, not newly adopted cognitive axioms.', 'Old resonant handoff and driven low-band experiments are different settings of the same parameterized parent family.', 'The finite occupation cutoff is an exact declared model; no infinite-cutoff limit is certified.', 'All ordered low levels and the full reference gap are controlled; product preparation and canonical dressed encoding have separate permissions and error bounds.', 'The true source includes the encoding connection; reference baseline, interaction energy and measurement work remain in the complete ledger.', 'Effective drive exchange is accounted for, but an autonomous quantum controller is not constructed.', 'The linear comparison deletes the balance contract and only matches the first-order block at the calibration point.', 'The finite drive source and Dirac energy block do not identify a local Higgs field, physical mass pole or complete standard-model realization.', 'Stage-wide parameter embeddings and all original completion criteria still require the joint paper and independent acceptance.'])
    readme_path, ledger_path = STAGE / 'README.md', STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    old_readme = readme_raw.decode('utf8')
    require(old_readme.count('## 当前进展') == 1, 'Unexpected README boundary')
    readme = old_readme.split('## 当前进展', 1)[0] + '## 当前进展\n\n2026-10-19。验收至1326，共四轮新研究；新增公理0，阶段尚未结项。[联合进展](_shared/progress_after1326.md)记录本轮结果；[目标恢复证据](_shared/goal_resumed_after1325.json)记录应用目标已恢复active，原目标内容未变。\n\n|轮次|新增结果|证据|\n|---|---|---|\n|[1323](research_note_1323.md)|全未知参考交接的单下降角色与归一原语|[数学独审](1323/independent_review.md) · [范围独审](1323/scope_review.md) · [验收](1323/acceptance.json)|\n|[1324](research_note_1324.md)|同父有限相干过程、完整后态、资源及留出预测|[数学独审](1324/independent_review.md) · [范围独审](1324/scope_review.md) · [验收](1324/acceptance.json)|\n|[1325](research_note_1325.md)|相反参考流共享左端的角色、可选重数匹配及相对味预算|[数学独审](1325/independent_review.md) · [范围独审](1325/scope_review.md) · [验收](1325/acceptance.json)|\n|[1326](research_note_1326.md)|动态参考的全部低带、真实驱动源、反作用及有限留出区间|[数学独审](1326/independent_review.md) · [范围独审](1326/scope_review.md) · [验收](1326/acceptance.json)|\n\n同一归一作用现连接全未知交付、有限相干准备和有隙动态参考。三种任务的初态量词与自由核设置不同；不能将旧共振交付与驱动低带并作一次实验。1326保持固定校准，两个有限响应区间严格分离，但f源尚未被认作Higgs场。\n\n[联合前件账](_shared/joint_premise_ledger_after1326.md)及[下一步](1326/NEXT.md)用于原五项标准的阶段论文与独立验收。不同数值见证须通过明确缩放／直和证明共同性；不将全部现实参数唯一、无截止完成或全部设备制造增作门槛。1327未准入。\n\n[本轮发布事务及前图像](_shared/progress_publication_after1326.json)仅保留一次1325后的导航边界；[1325发布记录](_shared/progress_publication_after1325.json)、[1324发布记录](_shared/progress_publication_after1324.json)、[1323发布记录](_shared/progress_publication_after1323.json)及[启动记录](_shared/navigation_start_receipt.json)保留不变。\n'
    ledger = ledger_raw.decode('utf8')
    before_intro = '2026-10-19。下表保留启动时的输入基线；科学成果现验收至1325，最新判定见文末。'
    require(before_intro in ledger, 'Unexpected ledger boundary')
    ledger = ledger.replace(before_intro, before_intro.replace('至1325', '至1326'), 1)
    ledger += '\n## 1326后的更新\n\n[报告](../research_note_1326.md)、[联合进展](progress_after1326.md)及[正式验收](../1326/acceptance.json)把已选归一原语接到同一有限动态参考的全部低带、任意伙伴、实际驱动源及测量工作。完整参考谱隙保证未显式枚举部门不进入低带，未将谱接近代签完整动力学。\n\n固定f=1第一阶校准后，f=2的有限低支响应区间严格不交；门户没有逐过程重调。线性对照删除余额合同，不是完整合同内反模型。参考自由核、驱动、动量字典与读口权限仍是公开输入，实际f源未被改名Higgs源。\n\n原共振交付、相干产品过程和驱动低带是同一参数化父族的不同设置。[联合前件账](joint_premise_ledger_after1326.md)保全部量词和1322竞争边界；[下一步](../1326/NEXT.md)是共同参数嵌入、阶段论文及原五项独立验收，尚未结项，1327未准入。\n\n[恢复证据](goal_resumed_after1325.json)记录应用目标恢复active、目标内容未变；1325的blocked观测仍作为历史记录保留。\n'
    prefix = '\n## 当前目标：基本相互作用规则的认知来源\n\n2026-10-19，[1323起阶段](archive_1323_/README.md)已验收1323—1326四轮，新增公理0，阶段尚未结项。[正式目标](archive_1323_/_shared/goal.md) · [1326报告](archive_1323_/research_note_1326.md) · [联合进展](archive_1323_/_shared/progress_after1326.md)。\n\n同一归一参考作用已连接全未知交接、有限相干转换、双端口角色／跨部门预算，以及动态参考的全部低带、真实驱动源和反作用。一次校准后的有限响应区间严格分离；全未知交付和驱动低带仍是不同实验设置，f源尚未取得Higgs身份。\n\n接下来按[联合前件账](archive_1323_/_shared/joint_premise_ledger_after1326.md)与[下一步](archive_1323_/1326/NEXT.md)完成共同参数嵌入、阶段论文和原五项独立验收。1327未准入，不再重复参考余额扫描，也不新增全UV或全部仪器门槛。应用目标恢复active的时间证据见[恢复记录](archive_1323_/_shared/goal_resumed_after1325.json)。\n\n上一阶段仍保存在[archive_1313_1322](archive_1313_1322/README.md)与[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)。[迁移与历史复算](_migration/closure_1313_1322/README.md)及下方历史记录保留；当前安排以本节与现行目标为准。\n'
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
    layer.update(latest_accepted_round=1326, in_progress_round=None, next_round=1327,
        next_round_admitted=False, stage_complete=False, application_goal_status_last_observed=observation['status'],
        research_continues=True,
        startup_state='Round 1326 accepted with conditional scope; next round not admitted')
    layer_raw = encoded(layer)
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger.encode('utf8')), (LAYER_PATH, old_raw, layer_raw)]
    planned = {p.resolve(): after for p, before, after in updates}
    planned[ACCEPTANCE.resolve()] = encoded(acceptance)
    fresh_paths = {CHECKS.resolve(), RECEIPT.resolve()}
    md = [p for p in package if p.suffix == '.md'] + [readme_path, ledger_path,
          STAGE / '_shared/progress_after1326.md', goal,
          STAGE / '_shared/joint_premise_ledger_after1326.md']
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
    sources = source_map(package + [goal, replay_path, STAGE / '_shared/progress_after1326.md',
        observation_path, STAGE / '_shared/joint_premise_ledger_after1326.md',
        STAGE / '_shared/start_stage.py', STAGE / '_shared/goal_start_receipt.json',
        ROOT / 'scripts/run_research_interaction_rules_active.py', Path(__file__).resolve(),
        prior_acceptance, prior_receipt, prior_checks,
        STAGE / '_shared/publish_progress_after1325.py'])
    checks = dict(status='PASS', round=1326, scientific_evidence_files=10,
        markdown_files=len(md), local_links=links, previous_round_bindings=len(prior['sha256']), original_navigation_recoverable=True,
        prior_navigation=prior_navigation, scientific_files_modified=False,
        stage_complete=False, next_round_admitted=False)
    receipt = dict(status='PUBLISHED', latest_accepted_round=1326, new_scientific_rounds=[1326],
        new_axioms=0, research_continues=True,
        application_goal_status_last_observed=observation['status'],
        stage_complete=False, next_available_round=1327,
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
