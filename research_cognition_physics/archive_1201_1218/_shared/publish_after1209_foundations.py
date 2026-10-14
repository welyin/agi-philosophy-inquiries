"""Index the foundations audit without accepting a new scientific round."""
from pathlib import Path
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1201_'
HERE = STAGE / '_shared'
OLD = BASE / 'archive_1161_'
PREV = HERE / 'progress_checks_1209.json'
OUT = HERE / 'after1209_foundation_audit_checks.json'
AUDITS = [HERE / name for name in [
    'after1209_thermodynamic_foundations.md',
    'thermodynamics_goal_scope_audit_after1209.md',
    'thermal_operational_bridge_audit_after1209.md']]
spec = importlib.util.spec_from_file_location('thermo_audit_tools', HERE / 'publish_1206.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links


def verify():
    d = read(OUT)
    ns, nb = m.frozen(d)
    m.m.validate_map(d['management_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    ledger = read(HERE / 'dependency_ledger.json')
    assert ledger['latest_accepted_round'] == 1209 and ledger['stage_check_groups'] == 18
    assert ledger['latest_management_receipt'] == OUT.name
    assert ledger['goal_objective_sha256'] == d['goal_objective_sha256']
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    print(json.dumps(dict(status='PASS', latest_round=1209, stage_rounds=9,
        stage_check_groups=18, new_scientific_rounds=0, new_adopted_axioms=0,
        old_GR_science_files_verified=ns, old_GR_acceptance_bindings_verified=nb,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        goal_changed=False, stage_complete=False)))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists(), 'Already indexed'
    assert all(p.is_file() for p in AUDITS)
    prev = read(PREV)
    ns, nb = m.frozen(prev)
    m.m.validate_map(prev['management_sha256'])
    index_line = ('\n- [1209后热力学认知基础审计](_shared/after1209_thermodynamic_foundations.md)：'
        '区分有限交付、参考角色与几何桥；不新增轮次或采用公理。'
        '[独立目标审计](_shared/thermodynamics_goal_scope_audit_after1209.md)、'
        '[功热及共同身份审计](_shared/thermal_operational_bridge_audit_after1209.md)。\n')
    for name in ['README.md', '文件索引.md']:
        p = STAGE / name
        value = p.read_text(encoding='utf-8')
        assert index_line not in value
        value = value.replace('_shared/progress_checks_1209.json', '_shared/' + OUT.name)
        p.write_text(value + index_line, encoding='utf-8')
    ledger = read(HERE / 'dependency_ledger.json')
    ledger['latest_management_receipt'] = OUT.name
    ledger['current_scope_audit'] = AUDITS[0].name
    req = ledger['stage_requirement_status']
    req['independently_generated_stable_reference'] = (
        '1202-1205 satisfy the independent finite preparation branch of the goal; '
        'all-n stability is a distinct unadopted reference-role candidate, not a retroactive requirement')
    req['thermodynamic_work_heat_scope'] = (
        'Driven material energy and directional exchange are defined in the fixed process; '
        'no universal strong-coupling work/heat split or arbitrary-cycle work extraction is claimed')
    req['stage_synthesis'] = (
        'Needs an integrated scoped deliverable and explicit reference-role/geometry handoff; '
        'the 1209 classification alone is neither cognitive adoption nor stage completion')
    write(HERE / 'dependency_ledger.json', ledger)
    oldledger = read(OLD / '_shared/dependency_ledger.json')
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(OLD / '_shared/dependency_ledger.json', oldledger)
    p = OLD / '_shared/navigation_notice.md'
    value = p.read_text(encoding='utf-8')
    value = value.replace('archive_1201_/_shared/progress_checks_1209.json',
                          'archive_1201_/_shared/' + OUT.name)
    needle = '。旧GR阶段40轮、130组保留；目标继续进行。'
    assert value.count(needle) == 1
    value = value.replace(needle, '。[热力学基础审计](archive_1201_/_shared/' +
        AUDITS[0].name + ')区分有限交付与待审角色，不新增轮次。旧GR阶段40轮、130组保留；目标继续进行。')
    p.write_text(value, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navverify['status'] == 'PASS'
    write(OUT, dict(status='PUBLISHING'))
    docs = AUDITS + [STAGE / 'README.md', STAGE / '文件索引.md'] + [BASE / name for name in
        ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']]
    count = links(docs)
    for p in docs:
        value = p.read_text(encoding='utf-8')
        assert value.count('$$') % 2 == 0, p
        assert not any(ord(c) < 32 and c not in '\r\n\t' for c in value), p
    audit = dict(prev['audit_files_sha256'])
    audit.update({p.relative_to(ROOT).as_posix(): sha(p) for p in AUDITS})
    management = {name: sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    write(OUT, dict(schema='thermodynamic_foundations_audit_after1209_v1', status='PASS',
        date='2026-10-14', accepted_rounds=prev['accepted_rounds'], stage_check_groups=18,
        new_scientific_rounds=0, new_adopted_axioms=0, new_empirical_calibration=0,
        GR_hard_deliverables_completed=0, application_goal_status='active',
        stage_complete=False, goal_changed=False,
        goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256=prev['acceptance_sha256'], science_files_sha256=prev['science_files_sha256'],
        screen_files_sha256=prev['screen_files_sha256'], audit_files_sha256=audit,
        management_sha256=management, prior_progress_path=PREV.relative_to(ROOT).as_posix(),
        prior_progress_sha256=sha(PREV), navigation_apply=navapply, navigation_check=navcheck,
        navigation_verify=navverify, local_links_checked=count, missing_links=[],
        old_GR_science_files_verified=ns, old_GR_acceptance_bindings_verified=nb,
        stage_requirement_status=req,
        scope='User concern and scope review; no new theorem, axiom adoption, or stage closure'))
    verify()


if __name__ == '__main__':
    main()
