"""Synchronize the user-authorized four-laws goal while preserving science."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1201_'
HERE = STAGE / '_shared'
OLD = BASE / 'archive_1161_'
PREV = HERE / 'after1209_foundation_audit_checks.json'
OUT = HERE / 'four_laws_goal_start_checks.json'
spec = importlib.util.spec_from_file_location('thermo_goal_tools', HERE / 'publish_1206.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links


def verify():
    d = read(OUT)
    ns, nb = m.frozen(d)
    m.m.validate_map(d['management_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    app = read(HERE / 'goal_four_laws_tool_receipt.json')['goal']
    objective = (HERE / 'goal_four_laws_prompt.txt').read_text(encoding='utf-8').strip()
    assert app['objective'] == objective and app['status'] == 'active'
    digest = hashlib.sha256(objective.encode('utf-8')).hexdigest()
    ledger = read(HERE / 'dependency_ledger.json')
    assert digest == ledger['goal_objective_sha256'] == d['goal_objective_sha256']
    assert ledger['goal_document'] == 'goal_four_laws.md'
    assert ledger['latest_management_receipt'] == OUT.name
    assert ledger['latest_accepted_round'] == 1209 and ledger['stage_check_groups'] == 18
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    print(json.dumps(dict(status='PASS', application_goal_status='active', goal_changed=True,
        target='thermodynamic_laws_0_1_2_3', latest_round=1209, stage_rounds=9,
        stage_check_groups=18, new_scientific_rounds=0, new_adopted_axioms=0,
        old_GR_science_files_verified=ns, old_GR_acceptance_bindings_verified=nb,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        prior_goal_marked_complete=False, stage_complete=False)))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists(), 'Goal revision already published'
    prev = read(PREV)
    ns, nb = m.frozen(prev)
    m.m.validate_map(prev['management_sha256'])
    app = read(HERE / 'goal_four_laws_tool_receipt.json')['goal']
    objective = (HERE / 'goal_four_laws_prompt.txt').read_text(encoding='utf-8').strip()
    assert objective == app['objective'] and app['status'] == 'active'
    digest = hashlib.sha256(objective.encode('utf-8')).hexdigest()
    status = ('当前目标已按用户要求改为从认知操作推导热力学第零、第一、第二、第三定律，应用状态active。'
        '已验收1201—1209共9轮18组作为基础；有限实例与参考选择不等于四律已完成。'
        '下一项先核四律命题、认知依赖与真实权限，第三定律重点审低温资源及纯态准备。')
    ledger = read(HERE / 'dependency_ledger.json')
    ledger['finite_process_requirement_status_before_four_laws'] = ledger['stage_requirement_status']
    req = dict(zeroth_law='finite comparisons exist; equilibrium equivalence and common scale proof pending',
        first_law='full finite process energy accounting exists; common work/heat definitions and premise chain to consolidate',
        second_law='relative entropy and reference selection results exist; compatible cyclic/feedback scope to prove',
        third_law='low-temperature and unattainability statements require new scoped derivation and resource audit',
        common_process='1201-1205 finite witness retained within its original scope',
        geometry_handoff='same clock/region/flow and independent area bridges remain explicit',
        stage_synthesis='four-laws paper and independent per-law audit pending',
        stage_deliverables_all_completed=False)
    ledger.update(goal_document='goal_four_laws.md', goal_start_receipt='goal_four_laws_tool_receipt.json',
        goal_objective_sha256=digest, prior_goal_objective_sha256=prev['goal_objective_sha256'],
        prior_goal_document='goal.md', prior_goal_completed=False,
        current_work=status, current_direction_document='goal_four_laws.md',
        latest_management_receipt=OUT.name, stage_requirement_status=req,
        stage_deliverables=['四大定律的明确命题、解析证明及相容适用范围',
            '认知要求到操作条件到定理的依赖与剩余输入表',
            '同一过程有限预测及完整资源/误差证书',
            '阶段论文、逐律独审及几何桥交接'])
    write(HERE / 'dependency_ledger.json', ledger)
    p = STAGE / 'README.md'
    value = p.read_text(encoding='utf-8')
    parts = value.split('\n\n')
    parts[0] = '# 1201起：认知操作热力学与四大定律'
    parts[1] = status
    parts[2] = ('[当前四律目标](_shared/goal_four_laws.md) · '
        '[应用回执](_shared/goal_four_laws_tool_receipt.json) · [旧目标（历史）](_shared/goal.md)。')
    parts[3] = ('推进顺序：认知操作 → 热力学四大定律 → 热力学—几何对应 → 有效GR。'
        '新目标按四律逐项验收，既有有限过程、参考分类及误差证据保留。')
    value = '\n\n'.join(parts).replace('_shared/after1209_foundation_audit_checks.json', '_shared/' + OUT.name)
    value = value.replace('[当前下一项](1209/NEXT.md)', '[1209原安排](1209/NEXT.md)')
    value += '\n\n[当前下一项与四律验收](_shared/goal_four_laws.md)：取代旧的有限交付结项安排，不追溯标记旧目标完成。\n'
    p.write_text(value, encoding='utf-8')
    p = STAGE / '文件索引.md'
    value = p.read_text(encoding='utf-8')
    value = value.replace('[完整目标](_shared/goal.md)', '[当前四律目标](_shared/goal_four_laws.md)')
    value = value.replace('[启动回执](_shared/goal_start_receipt.json)', '[当前启动回执](_shared/goal_four_laws_tool_receipt.json)')
    value = value.replace('[当前下一项](1209/NEXT.md)', '[当前下一项](_shared/goal_four_laws.md)；[1209原安排](1209/NEXT.md)')
    value = value.replace('_shared/after1209_foundation_audit_checks.json', '_shared/' + OUT.name)
    value += '\n- 2026-10-14用户调整目标：[四大定律完整正文](_shared/goal_four_laws_prompt.txt)。未新增科学轮次；[旧目标](_shared/goal.md)保留。\n'
    p.write_text(value, encoding='utf-8')
    p = OLD / '_shared/navigation_notice.md'
    _, history = p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    prefix = ('## 当前目标（1201阶段）：从认知操作推导热力学四大定律\n\n'
        '2026-10-14，用户明确调整目标并删除旧应用目标后，新的四律目标已成功启用，状态 **active**。\n\n'
        '[阶段目录](archive_1201_/README.md) · [现行目标](archive_1201_/_shared/goal_four_laws.md) · '
        '[应用回执](archive_1201_/_shared/goal_four_laws_tool_receipt.json) · [阶段索引](archive_1201_/文件索引.md)。\n\n'
        '分别推导热平衡与温标、完整能源与功热、熵/循环限制、低温熵与零温不可达性；'
        '四律须在相容的共同操作与资源范围内成立。额外条件研究其认知来源，不将待证定律换名公理化。'
        '普通热力学不自动给Unruh、面积熵或GR；几何身份继续明确交接。\n\n'
        '**当前状态：** ' + status + ' [最新报告](archive_1201_/research_note_1209.md) · '
        '[当前核验](archive_1201_/_shared/' + OUT.name + ')。旧GR与旧有限热力学目标未追溯标为完成。\n\n')
    p.write_text(prefix + '以下保留前阶段的历史状态' + history, encoding='utf-8')
    oldledger = read(OLD / '_shared/dependency_ledger.json')
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(OLD / '_shared/dependency_ledger.json', oldledger)
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navverify['status'] == 'PASS'
    write(OUT, dict(status='PUBLISHING'))
    docs = [HERE / 'goal_four_laws.md', STAGE / 'README.md', STAGE / '文件索引.md'] + [BASE / name for name in
        ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']]
    count = links(docs)
    for p in docs:
        value = p.read_text(encoding='utf-8')
        assert value.count('$$') % 2 == 0, p
        assert not any(ord(c) < 32 and c not in '\r\n\t' for c in value), p
    audit = dict(prev['audit_files_sha256'])
    for name in ['goal_four_laws.md', 'goal_four_laws_prompt.txt', 'goal_four_laws_tool_receipt.json']:
        p = HERE / name
        audit[p.relative_to(ROOT).as_posix()] = sha(p)
    management = {name: sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    write(OUT, dict(schema='four_laws_goal_start_v1', status='PASS', date='2026-10-14',
        accepted_rounds=prev['accepted_rounds'], stage_check_groups=18,
        new_scientific_rounds=0, new_adopted_axioms=0, new_empirical_calibration=0,
        GR_hard_deliverables_completed=0, application_goal_status='active',
        goal_changed=True, user_authorized_goal_revision=True, prior_goal_marked_complete=False,
        goal_objective_sha256=digest, prior_goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256=prev['acceptance_sha256'], science_files_sha256=prev['science_files_sha256'],
        screen_files_sha256=prev['screen_files_sha256'], audit_files_sha256=audit,
        management_sha256=management, prior_progress_path=PREV.relative_to(ROOT).as_posix(),
        prior_progress_sha256=sha(PREV), navigation_apply=navapply, navigation_check=navcheck,
        navigation_verify=navverify, local_links_checked=count, missing_links=[],
        old_GR_science_files_verified=ns, old_GR_acceptance_bindings_verified=nb,
        stage_requirement_status=req, stage_complete=False,
        prior_goal_turn_classification='progress: audited and indexed thermodynamic scope without new science',
        current_turn_classification='user-authorized goal replacement; four-laws acceptance and navigation synchronized'))
    verify()


if __name__ == '__main__':
    main()
