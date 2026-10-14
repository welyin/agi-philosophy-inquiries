"""Publish a priority/scope audit without creating or accepting a research round."""
from pathlib import Path
import importlib.util
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1201_'
HERE = STAGE / '_shared'
OLD = BASE / 'archive_1161_'
OUT = HERE / 'conditions_audit_checks.json'
PREV = HERE / 'progress_checks_1205.json'
AUDIT = HERE / 'thermodynamic_conditions_geometry_audit.md'
spec = importlib.util.spec_from_file_location('thermo_publish_tools', HERE / 'publish_1201.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links


def science(receipt):
    old, bindings = m.old_check(read(HERE / 'startup_checks.json'))
    for key in ['science_files_sha256', 'screen_files_sha256']:
        m.validate_map(receipt[key])
    for number, digest in receipt['acceptance_sha256'].items():
        p = STAGE / number / 'acceptance.json'
        assert sha(p) == digest
        for name, h in read(p)['files_sha256'].items():
            assert sha((p.parent / name).resolve()) == h
    return len(old['science_files_sha256']), bindings


def verify():
    d = read(OUT)
    nscience, bindings = science(d)
    m.validate_map(d['management_sha256'])
    m.validate_map(d['audit_files_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    assert not (STAGE / '1206').exists(), 'Unaccepted preflight must remain unnumbered'
    ledger = read(HERE / 'dependency_ledger.json')
    assert ledger['latest_accepted_round'] == 1205
    assert ledger['stage_check_groups'] == 10
    assert ledger['goal_objective_sha256'] == d['goal_objective_sha256']
    print(json.dumps({'status': 'PASS', 'latest_round': 1205, 'stage_rounds': 5,
                      'stage_check_groups': 10, 'new_scientific_rounds': 0,
                      'old_science_files_verified': nscience,
                      'old_acceptance_bindings_verified': bindings,
                      'stage_complete': False, 'goal_changed': False}))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists(), 'Audit already published'
    assert AUDIT.is_file(), 'Reviewed condition audit required'
    prev = read(PREV)
    science(prev)
    m.validate_map(prev['management_sha256'])
    assert not (STAGE / '1206').exists()
    preflight = run(HERE / 'whole_trial_pointer_preflight/independent_check.py')
    assert preflight['all_checks_passed'] is True
    assert preflight == read(HERE / 'whole_trial_pointer_preflight/independent_results.json')
    status = ('最新验收1205；新阶段5轮、10组。固定任务的准备、内部接触、记录及完整资源误差已接续。'
              '按用户提醒，优先审计可扩展平衡、许可操作与实际能源／时间／区域的共同身份，'
              '复用1020、305及317等已有定理；暂缓阶段结项和器件统计细化。')
    ledger = read(HERE / 'dependency_ledger.json')
    ledger.update(current_work=status,
                  current_direction_document='thermodynamic_conditions_geometry_audit.md',
                  latest_management_receipt=OUT.name)
    req = ledger['stage_requirement_status']
    req.update(independently_generated_stable_reference=(
        '1205 certifies independent fixed-window preparation; equilibrium stability under '
        'composed actual cycles and its cognitive justification are not thereby established'),
        stage_synthesis=(
        'draft exists; closure deferred for structural-condition audit requested by user; '
        'separate instance choices, cognitive conditions and geometry identity'),
        thermodynamic_structural_basis=(
        'reuse 1020/305/317 and later solved domain bridges; next specify actual permitted '
        'cycles, retained resources and same reference/energy/clock/region identity'),
        stage_deliverables_all_completed=False)
    write(HERE / 'dependency_ledger.json', ledger)
    p = STAGE / 'README.md'
    paragraphs = p.read_text(encoding='utf-8').split('\n\n')
    paragraphs[1] = status
    text = '\n\n'.join(paragraphs).replace('_shared/progress_checks_1205.json',
                                         '_shared/' + OUT.name)
    text += ('\n\n[热力学条件与几何桥审计](_shared/thermodynamic_conditions_geometry_audit.md)：'
             '本次用户提醒后的当前方向；1205/NEXT保留历史安排。'
             '全部原始指针统计预审暂存于[_shared](_shared/whole_trial_pointer_preflight/README.md)，不计1206。\n')
    p.write_text(text, encoding='utf-8')
    p = STAGE / '文件索引.md'
    text = p.read_text(encoding='utf-8').replace(
        '[下一项](1205/NEXT.md)', '[当前下一项](_shared/thermodynamic_conditions_geometry_audit.md)；[1205原安排](1205/NEXT.md)')
    text = text.replace('_shared/progress_checks_1205.json', '_shared/' + OUT.name)
    text += ('\n- [热力学条件与几何桥审计](_shared/thermodynamic_conditions_geometry_audit.md)：'
             '管理与研究取舍，不增加轮次、公理或完成计数。\n'
             '- [暂缓的全部原始指针预审](_shared/whole_trial_pointer_preflight/README.md)\n')
    p.write_text(text, encoding='utf-8')
    p = OLD / '_shared/dependency_ledger.json'
    oldledger = read(p)
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(p, oldledger)
    p = OLD / '_shared/navigation_notice.md'
    before, after = p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    before = re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)', '**当前状态：** ' + status +
        ' [1205报告](archive_1201_/research_note_1205.md) · '
        '[条件审计与下一项](archive_1201_/_shared/thermodynamic_conditions_geometry_audit.md) · '
        '[当前核验](archive_1201_/_shared/' + OUT.name + ')。旧GR阶段40轮、130组保留；目标继续进行。',
        before, flags=re.S)
    p.write_text(before + '以下保留前阶段的历史状态' + after, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navverify['status'] == 'PASS'
    write(OUT, {'status': 'PUBLISHING'})
    documents = [AUDIT, HERE / 'whole_trial_pointer_preflight/README.md', STAGE / 'README.md',
                 STAGE / '文件索引.md', *[BASE / n for n in
                 ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']]]
    nlinks = links(documents)
    for p in documents:
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in p.read_text(encoding='utf-8'))
    mgmt = {name: sha(ROOT / name) for name in prev['management_sha256']}
    mgmt[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    auditfiles = [AUDIT, *sorted((HERE / 'whole_trial_pointer_preflight').glob('*'))]
    write(OUT, {'schema': 'thermodynamics_structural_conditions_audit_v1',
        'date': '2026-10-14', 'status': 'PASS', 'accepted_rounds': list(range(1201, 1206)),
        'stage_check_groups': 10, 'new_scientific_rounds': 0, 'new_adopted_axioms': 0,
        'new_empirical_calibration': 0, 'GR_hard_deliverables_completed': 0,
        'application_goal_status': 'active', 'stage_complete': False, 'goal_changed': False,
        'goal_objective_sha256': ledger['goal_objective_sha256'],
        'acceptance_sha256': prev['acceptance_sha256'],
        'science_files_sha256': prev['science_files_sha256'],
        'screen_files_sha256': prev['screen_files_sha256'], 'management_sha256': mgmt,
        'audit_files_sha256': {p.relative_to(ROOT).as_posix(): sha(p) for p in auditfiles},
        'prior_progress_path': PREV.relative_to(ROOT).as_posix(), 'prior_progress_sha256': sha(PREV),
        'navigation_apply': navapply, 'navigation_check': navcheck, 'navigation_verify': navverify,
        'local_links_checked': nlinks, 'missing_links': [], 'stage_requirement_status': req,
        'deferred_preflight_reproduced': True,
        'scope': 'priority and structural-condition audit; frozen science untouched; no new theorem acceptance'})
    verify()


if __name__ == '__main__':
    main()
