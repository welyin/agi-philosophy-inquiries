"""Finalize the declared stage only after every round and joint review is frozen."""
import copy
import json
import re
from urllib.parse import unquote
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout

PAPER = ROOT / 'research_cognition_physics/15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md'
LEDGER = STAGE / '_shared/joint_scope_ledger_after1322.md'
REVIEW = STAGE / '_shared/stage_paper_review_after1322.md'


def main():
    acceptance_path = STAGE / '_shared/stage_acceptance_after1322.json'
    checks_path = STAGE / '_shared/stage_checks_after1322.json'
    receipt_path = STAGE / '_shared/stage_publication_after1322.json'
    assert not any(p.exists() for p in (acceptance_path, checks_path, receipt_path))
    navigation_check = InteractionGenerationLayout().verify_interaction_navigation()
    old_raw = LAYER.read_bytes()
    old = json.loads(old_raw)
    assert old['latest_accepted_round'] == 1322 and old['in_progress_round'] is None
    assert old['stage_complete'] is False
    round_packages, bindings = [], {}
    for n in range(1313, 1323):
        p = STAGE / str(n) / 'acceptance.json'
        acceptance = json.loads(p.read_bytes())
        assert acceptance['round'] == n
        assert acceptance['status'] == 'ACCEPTED_CONDITIONAL_SCOPE'
        assert acceptance['new_axioms'] == 0
        assert len(acceptance['executed_checks']) == 2
        assert all(c['exit_code'] == 0 for c in acceptance['executed_checks'])
        for name, digest in acceptance['sha256'].items():
            assert sha((ROOT / name).read_bytes()) == digest, name
            if name in bindings:
                assert bindings[name] == digest
            bindings[name] = digest
        bindings[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
        round_packages.append(dict(round=n, acceptance_sha256=sha(p.read_bytes()),
                                   bound_files=len(acceptance['sha256'])))
    review = REVIEW.read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    for p in (PAPER, LEDGER, STAGE / '_shared/goal.md', STAGE / '1322/acceptance.json'):
        assert sha(p.read_bytes()) in review, str(p)
        bindings[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    bindings[REVIEW.relative_to(ROOT).as_posix()] = sha(REVIEW.read_bytes())
    for key in ('A', 'B', 'C', 'J'):
        assert key in review
    paper_text = PAPER.read_text(encoding='utf8')
    assert '本稿进入联合审阅' not in paper_text
    acceptance = dict(
        status='ACCEPTED_DECLARED_STAGE_SCOPE', date='2026-10-19',
        scientific_rounds=list(range(1313, 1323)), round_count=10,
        author_recomputations=10, independent_recomputations=10,
        new_axioms=0, new_scientific_rounds_in_synthesis=0,
        stage_goal='从认知交互机制生成物质分工、代际与质量结构',
        requirements={
            'open_candidates_and_dependencies': 'PASS: layered task contracts and open/fixed/unclassified inputs are separated in the paper and joint ledger.',
            'mechanism_selection_A_B_C': 'PASS_WITH_SCOPE: actual role selection and common counterfamilies; unit-exchange obstruction; parameter relations and remaining mass families.',
            'substantive_selection_gain': 'PASS: internal full-code dimensions/characters and common-source parameter constraints reduce independent inputs under declared mechanisms.',
            'joint_nonempty_model_and_predictions': 'PASS_WITH_SCOPE: 1319--1322 provide explicit field/interface/finite-parent dictionaries, partner-preserving failure channels, sources and bounded predictions.',
            'reproducibility_and_independent_acceptance': 'PASS: ten separately accepted rounds and a final independently reviewed paper; historical results remain unchanged.'},
        main_joint_chain=[1319, 1320, 1321, 1322],
        separate_branches=[[1314, 1315], [1316, 1317], [1318]],
        actual_role_comparison_baseline=[1313],
        limits=[
            'This ends the declared mechanism stage, not unique Standard Model or empirical-parameter generation.',
            'Fixed R_rep and dimension-four mass menus are not refuted by the R_op vectorlike alternative.',
            'Three/four/five-family common feasibility concerns the fixed 1% contract, not all thresholds.',
            'The full hypercharge family, actual gauge group and all matter species are not classified in the strong common model.',
            'Finite-mode certificates do not certify dynamic gauge/Higgs instruments, all unknown flavors or continuum UV completion.',
            'Separate branches are not stitched into a single Standard Model.'],
        independent_review=REVIEW.relative_to(ROOT).as_posix(),
        sha256=bindings,
        round_packages=round_packages)
    readme_path = STAGE / 'README.md'
    ledger_path = STAGE / '_shared/selection_ledger.md'
    readme_raw, ledger_raw = readme_path.read_bytes(), ledger_path.read_bytes()
    readme = readme_raw.decode('utf8')
    before = '最新验收至1322，共十轮新研究，作者与独立复算均通过。阶段目标保持active，新公理0。'
    after = '1313—1322共十轮已完成逐轮及联合验收，本阶段按原目标声明范围结项。作者与独立复算各10组，新公理0。'
    assert before in readme
    readme = readme.replace(before, after, 1)
    start, end = readme.index('科学轮次验收至1322，进入'), readme.index('历史启动验收')
    readme = readme[:start] + '''[阶段论文15](../15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)、[共同范围与输入账](_shared/joint_scope_ledger_after1322.md)、[最终独审](_shared/stage_paper_review_after1322.md)及[逐项验收](_shared/stage_acceptance_after1322.json)已保存。整合计0轮；下一科学轮未启动。

主链1319—1322给角色容量与字符限制、共同参数关系及真实剩余族。其它机制分支分别保留自己的前件，未拼成一套标准模型。完整现实参数、实际群选择和强合同下全部超荷分类仍未获得；长期生成目标继续作为未来方向，不自动重设新目标。

''' + readme[end:]
    ledger_add = '''
## 本阶段联合验收

2026-10-19。十轮科学包及[论文15](../../15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md)按原goal的五项交付验收，见[共同分类](joint_scope_ledger_after1322.md)、[独审](stage_paper_review_after1322.md)与[正式验收](stage_acceptance_after1322.json)。新科学轮次0，新公理0。得到真实机制选择及有范围剩余族；未分类输入仍明确开放。没有宣称唯一生成完整标准模型或经验参数，不自动开启下一目标。
'''
    prefix = '''
## 已完成阶段：认知交接机制的角色选择与参数剩余自由

2026-10-19，[1313—1322阶段](archive_1313_/README.md)十轮已通过逐轮和联合验收。[论文15](15_认知交接机制的物质角色选择与参数剩余自由_阶段论文.md) · [共同范围与剩余输入](archive_1313_/_shared/joint_scope_ledger_after1322.md) · [正式验收](archive_1313_/_shared/stage_acceptance_after1322.json)。作者／独立核验各10组，新公理0。

完整内部交接条件选择角色容量与字符；固定同时门户把交付误差与同源交换系数联系。同一有限父模型仍保三／四／五代及质量开族，右弱角色替代具有相同有限任务而不同物理源。各分支的真实增益、剩余族和未分类输入分别列明。

本阶段按原目标声明范围完成；未唯一生成完整标准模型或现实参数，也未把有限截止证书当作全场仪器。下一阶段应先提出能区分现有竞争族的独立任务／来源要求，不自动开始新编号或采用新公理。上一阶段已归档为[archive_1308_1312](archive_1308_1312/README.md)。下方旧时点记录不改。
'''
    layer, updates = copy.deepcopy(old), []
    for entry in layer['entries']:
        p = ROOT / entry['destination']
        before_raw = p.read_bytes()
        original = relocation_original_bytes(before_raw, entry)
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        content, at = original.decode(enc), entry['edits'][0]['start']
        after_raw = (content[:at] + prefix + content[at:]).encode(enc)
        entry.update(current_sha256=sha(after_raw), current_bytes=len(after_raw))
        entry['edits'][0]['after'] = prefix
        assert relocation_original_bytes(after_raw, entry) == original
        updates.append((p, before_raw, after_raw))
    layer.update(stage_complete=True, goal_completion_ready=True, next_round=None,
                 startup_state='Ten accepted rounds and final paper independently accepted within original stage scope')
    updates += [(readme_path, readme_raw, readme.encode('utf8')),
                (ledger_path, ledger_raw, ledger_raw + ledger_add.encode('utf8')),
                (LAYER, old_raw, encoded(layer))]
    planned = {p: after_raw for p, _, after_raw in updates}
    planned[acceptance_path] = encoded(acceptance)
    md = [PAPER, LEDGER, REVIEW, STAGE / '_shared/goal.md', readme_path, ledger_path]
    links = 0
    for p in md:
        text = planned[p].decode('utf8') if p in planned else p.read_text(encoding='utf8')
        assert text.count('$$') % 2 == 0, str(p)
        assert not re.search(r'^\s*\\[\[\]]\s*$', text, re.M), str(p)
        for target in re.findall(r'\]\(([^)\s]+)\)', text):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (p.parent / unquote(target.split('#')[0])).resolve()
            assert resolved.exists() or resolved in planned, (p, target)
            links += 1
    checks = dict(status='PASS', scientific_packages=10, frozen_bindings=len(bindings),
                  author_checks=10, independent_checks=10, markdown_files=len(md), local_links=links,
                  navigation_before=navigation_check, original_navigation_recoverable=True,
                  new_axioms=0, scope='Evidence integrity and publication checks; mathematical judgment is the independent review.')
    receipt = dict(status='STAGE_PUBLISHED', new_scientific_rounds=0, new_axioms=0,
        previous_layer=old, current_layer_sha256=sha(encoded(layer)),
        acceptance_sha256=sha(encoded(acceptance)),
        application_goal_status_last_observed='active', application_completion_call_pending=True,
        previous_stage_navigation={p.relative_to(ROOT).as_posix(): b.decode('utf8')
                                   for p, b in ((readme_path, readme_raw), (ledger_path, ledger_raw))})
    fresh = [(acceptance_path, encoded(acceptance)), (checks_path, encoded(checks)), (receipt_path, encoded(receipt))]
    for p, before_raw, _ in updates:
        assert p.read_bytes() == before_raw, str(p)
    written, created = [], []
    try:
        for p, before_raw, after_raw in updates:
            atomic(p, after_raw)
            written.append((p, before_raw, after_raw))
        for p, raw in fresh:
            with p.open('xb') as stream:
                stream.write(raw)
            created.append((p, raw))
    except BaseException:
        for p, raw in reversed(created):
            if p.read_bytes() == raw:
                p.unlink()
        for p, before_raw, after_raw in reversed(written):
            if p.read_bytes() == after_raw:
                atomic(p, before_raw)
        raise
    print(json.dumps({k: v for k, v in checks.items() if k != 'navigation_before'}))


if __name__ == '__main__':
    main()
