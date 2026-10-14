"""Accept the scoped native-reference selection result; preserve prior science."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1201_'
HERE = STAGE / '_shared'
ROUND = STAGE / '1206'
OLD = BASE / 'archive_1161_'
PREV = HERE / 'conditions_audit_checks.json'
OUT = HERE / 'progress_checks_1206.json'
spec = importlib.util.spec_from_file_location('thermo_publish_tools', HERE / 'publish_1201.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links


def frozen(receipt):
    old, bindings = m.old_check(read(HERE / 'startup_checks.json'))
    for key in ['science_files_sha256', 'screen_files_sha256', 'audit_files_sha256']:
        m.validate_map(receipt[key])
    for number, digest in receipt['acceptance_sha256'].items():
        p = STAGE / number / 'acceptance.json'
        assert sha(p) == digest
        for name, h in read(p)['files_sha256'].items():
            assert sha((p.parent / name).resolve()) == h, (number, name)
    return len(old['science_files_sha256']), bindings


def verify():
    d = read(OUT)
    ns, nb = frozen(d)
    m.validate_map(d['management_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    ledger = read(HERE / 'dependency_ledger.json')
    assert ledger['latest_accepted_round'] == 1206 and ledger['stage_check_groups'] == 12
    assert ledger['goal_objective_sha256'] == d['goal_objective_sha256']
    assert ledger['latest_management_receipt'] == OUT.name
    print(json.dumps({'status': 'PASS', 'latest_round': 1206, 'stage_rounds': 6,
        'stage_check_groups': 12, 'old_GR_science_files_verified': ns,
        'old_GR_acceptance_bindings_verified': nb,
        'current_stage_science_files_verified': len(d['science_files_sha256']),
        'goal_status': 'active', 'goal_changed': False, 'stage_complete': False}))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists() and not (ROUND / 'acceptance.json').exists(), 'Already published'
    prev = read(PREV)
    ns, nb = frozen(prev)
    m.validate_map(prev['management_sha256'])
    review = (ROUND / 'independent_review.md').read_text(encoding='utf-8')
    assert '限定签核通过' in review
    for name in ['proof.md', 'check.py', 'results.json', 'independent_check.py', 'independent_results.json']:
        assert sha(ROUND / name) in review, 'Review binding mismatch: ' + name
    author = run(ROUND / 'check.py')
    independent = run(ROUND / 'independent_check.py')
    assert author['status'] == 'AUTHOR_CHECKS_PASS' and author['results_match']
    assert independent['all_checks_passed'] and independent == read(ROUND / 'independent_results.json')
    r = read(ROUND / 'results.json')
    for result in [r, independent]:
        for name, digest in result['source_sha256'].items():
            assert sha(STAGE / name) == digest, name
    ac = r['algebraic_certificate']
    assert ac['discriminant'] == independent['native_spectrum']['discriminant'] == 2736128576000000
    assert ac['mod_11_values'] == independent['native_spectrum']['mod11_values_x0_to10']
    assert ac['number_energy_2x2_determinant'] == independent['fixed_n_resonance']['two_equation_determinant'] == '81/40'
    assert independent['qq_full_product_stationarity']['log_product_constraint_rank'] == 6
    assert r['native_QQ_graph']['edges'] == independent['qq_full_product_stationarity']['qq_undirected_edge_count'] == 18
    cold = independent['qq_full_product_stationarity']['cold_reference_commutator']
    assert F('.0003202089318') < F(cold['interval_lower_exact']) < F(cold['interval_upper_exact']) < F('.0003202089319')
    assert F(r['original_cold_reference_witness']['entry_absolute_lower']) < F(cold['interval_upper_exact'])
    assert F(cold['interval_lower_exact']) < F(r['original_cold_reference_witness']['entry_absolute_upper'])
    gap = independent['near_resonance_diagnostics']['samples'][-1]
    assert gap['n'] == 30 and gap['distinct_occupation_energies'] == r['finite_resolution_scope']['occupation_types'] == 5456
    assert F('.0004237212509') < F(gap['minimum_gap_lower']) < F(gap['minimum_gap_upper']) < F('.0004237212510')
    ranks = independent['permitted_conversion_criterion']['three_subject_qutrit_exact_ranks']
    assert [ranks[k] for k in ['cross_subject_exchange_only', 'with_one_conversion', 'with_conversion_but_delete_last_subject_edge']] == [4, 5, 3]
    files = ['proof.md', 'check.py', 'results.json', 'independent_check.py',
        'independent_results.json', 'independent_review.md', 'NEXT.md',
        '../research_note_1206.md', '../_shared/goal.md']
    hashes = {name: sha((ROUND / name).resolve()) for name in files}
    scope = ('原四能级材料全部有限固定份数的精确保裸能操作均固定任意对角同态产品边缘，'
        '不选择Gibbs；原完整QQ接触的满秩裸平稳独立产品全时不变仅有最大混合。'
        '给原冷参考非零交换子、近共振有限精度边界与另一材料的转换秩删弱例；不新增实际操作权限。')
    write(ROUND / 'acceptance.json', {
        'schema': 'operational_thermodynamics_round_acceptance_v1', 'round': 1206,
        'date': '2026-10-14', 'status': 'ACCEPTED_WITH_SCOPE', 'check_groups': 2,
        'scope': scope, 'certified': {
            'all_finite_fixed_n_integer_resonance_lattice': '{0}',
            'all_bare_energy_preserving_unitaries_fix_diagonal_iid_marginal': True,
            'arbitrary_purifying_reference_joint_state_fixed': False,
            'full_QQ_stationary_full_rank_product_unique': 'I4/4 tensor I4/4',
            'native_cold_commutator_strict_interval': ['.0003202089318', '.0003202089319'],
            'no_n_uniform_positive_occupation_gap': True,
            'finite_n30_occupation_types': 5456,
            'finite_n30_minimum_gap_strict_interval': ['.0004237212509', '.0004237212510'],
            'alternative_qutrit_exchange_conversion_deletion_ranks': [4, 5, 3],
            'finite_log_constraint_bound_requires_actual_preparation_implementation_readout_errors': True},
        'remaining_inputs': ['native h and Q, fixed copy number and selected operation permissions',
            'full-rank independent bare-stationary products for QQ uniqueness only',
            'explicit alternative conversion/readout permissions and finite error parameters',
            'full auxiliary return, work and correlation rules for future equilibrium selection',
            'same actual region/clock/physical-flow identity for geometry'],
        'limits': r['limits'], 'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'GR_hard_deliverables_completed': 0, 'stage_complete': False,
        'root_executed_both_codes': True,
        'review': 'Scoped independent standard-library exact algebra, separate spectral-projector cold witness, rational ranks and interval spectra; root proof/scope audit; no empirical test.',
        'files_sha256': hashes})
    status = ('最新验收1206；新阶段6轮、12组。原材料全部有限精确保裸能自接触不选择热人口；'
        '原完整QQ对独立平稳产品的全时不扰动仅留最大混合。下一项明确带完整辅助资源与关联的共同参考比较，'
        '复用1020；实际区域、钟和物理流的几何身份仍开放，暂缓阶段结项。')
    ledger = read(HERE / 'dependency_ledger.json')
    req = ledger['stage_requirement_status']
    req.update(thermodynamic_structural_basis=(
        '1206 proves all-finite exact bare-energy self-contact is nonselective and native full-QQ '
        'stationary independent-product selection is overstrong; next define full auxiliary/work/correlation '
        'return permissions using 1020; no axiom adopted'),
        stage_synthesis='closure deferred; structural reference selection and same-region physical-flow identity remain open',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1206, accepted_stage_rounds=list(range(1201, 1207)),
        stage_check_groups=12, current_direction_document='../1206/NEXT.md', current_work=status,
        latest_management_receipt=OUT.name, stage_requirement_status=req)
    write(HERE / 'dependency_ledger.json', ledger)
    p = STAGE / 'README.md'
    parts = p.read_text(encoding='utf-8').split('\n\n')
    parts[1] = status
    text = '\n\n'.join(parts).replace('_shared/conditions_audit_checks.json', '_shared/' + OUT.name)
    text = text.replace('## 接续与边界', '- [1206报告](research_note_1206.md)：全部有限自接触的参考选择边界；[验收](1206/acceptance.json)、[当前下一项](1206/NEXT.md)。\n\n## 接续与边界', 1)
    text = text.replace('本次用户提醒后的当前方向；1205/NEXT保留历史安排。',
        '本次用户提醒后的方向审计；[当前下一项](1206/NEXT.md)接续1206，1205/NEXT保留历史安排。')
    p.write_text(text, encoding='utf-8')
    p = STAGE / '文件索引.md'
    text = p.read_text(encoding='utf-8').replace('最新验收1205；5轮、10组', '最新验收1206；6轮、12组')
    text = text.replace('## 管理与接续', '|1206|[原材料的平衡选择与实际接触权限](research_note_1206.md)|[证明](1206/proof.md)、[代码](1206/check.py)、[结果](1206/results.json)、[独立核验](1206/independent_check.py)、[审阅](1206/independent_review.md)、[验收](1206/acceptance.json)|\n\n## 管理与接续', 1)
    text = re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)', r'\1\n', text)
    text = text.replace('[当前下一项](_shared/thermodynamic_conditions_geometry_audit.md)',
        '[当前下一项](1206/NEXT.md)；[结构条件审计](_shared/thermodynamic_conditions_geometry_audit.md)')
    text = text.replace('_shared/conditions_audit_checks.json', '_shared/' + OUT.name)
    text += '\n- [边界、参考与物理流的几何预审](_shared/thermal_geometry_selection_preflight/README.md)：复用旧条件定理，未认证当前几何样本。\n'
    p.write_text(text, encoding='utf-8')
    oldledger = read(OLD / '_shared/dependency_ledger.json')
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(OLD / '_shared/dependency_ledger.json', oldledger)
    p = OLD / '_shared/navigation_notice.md'
    prefix, rest = p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    prefix = re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)', '**当前状态：** ' + status +
        ' [1206报告](archive_1201_/research_note_1206.md) · [下一项](archive_1201_/1206/NEXT.md) · '
        '[当前核验](archive_1201_/_shared/' + OUT.name + ')。旧GR阶段40轮、130组保留；目标继续进行。', prefix, flags=re.S)
    p.write_text(prefix + '以下保留前阶段的历史状态' + rest, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navverify['status'] == 'PASS'
    write(OUT, {'status': 'PUBLISHING'})
    screen = HERE / 'thermal_geometry_selection_preflight/README.md'
    documents = [STAGE / 'README.md', STAGE / '文件索引.md', STAGE / 'research_note_1206.md',
        screen, *ROUND.glob('*.md'), *[BASE / name for name in
        ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']]]
    nlinks = links(documents)
    for p in documents:
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in p.read_text(encoding='utf-8'))
    science = dict(prev['science_files_sha256'])
    science.update({(ROUND / name).resolve().relative_to(ROOT).as_posix(): h for name, h in hashes.items()})
    management = {name: sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    screens = dict(prev['screen_files_sha256'])
    screens[screen.relative_to(ROOT).as_posix()] = sha(screen)
    accs = dict(prev['acceptance_sha256'])
    accs['1206'] = sha(ROUND / 'acceptance.json')
    write(OUT, {'schema': 'thermodynamics_progress_1206_v1', 'date': '2026-10-14', 'status': 'PASS',
        'accepted_rounds': list(range(1201, 1207)), 'stage_check_groups': 12,
        'application_goal_status': 'active', 'stage_complete': False, 'goal_changed': False,
        'goal_objective_sha256': ledger['goal_objective_sha256'], 'acceptance_sha256': accs,
        'science_files_sha256': science, 'management_sha256': management,
        'screen_files_sha256': screens, 'audit_files_sha256': prev['audit_files_sha256'],
        'prior_progress_path': PREV.relative_to(ROOT).as_posix(), 'prior_progress_sha256': sha(PREV),
        'previous_GR_science_files_verified': ns, 'previous_GR_acceptance_bindings_verified': nb,
        'navigation_apply': navapply, 'navigation_check': navcheck, 'navigation_verify': navverify,
        'markdown_files_checked': len(documents), 'local_links_checked': nlinks, 'missing_links': [],
        'stage_requirement_status': req, 'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'GR_hard_deliverables_completed': 0, 'science_program_outputs': {'author': author,
            'independent_all_checks_passed': independent['all_checks_passed'],
            'independent_result_sha256': sha(ROUND / 'independent_results.json')},
        'this_goal_turn_result': 'progress: native all-finite reference-selection classification, full-contact product no-go, exact finite witnesses and permission-sensitive alternative; stage remains active'})
    verify()


if __name__ == '__main__':
    main()
