"""Publish scoped composition-robust reference selection after independent review."""
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
ROUND = STAGE / '1207'
OLD = BASE / 'archive_1161_'
PREV = HERE / 'progress_checks_1206.json'
OUT = HERE / 'progress_checks_1207.json'
spec = importlib.util.spec_from_file_location('previous_thermo_publisher', HERE / 'publish_1206.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
sha, read, write, run, links = prior.sha, prior.read, prior.write, prior.run, prior.links


def verify():
    d = read(OUT)
    ns, nb = prior.frozen(d)
    prior.m.validate_map(d['management_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    ledger = read(HERE / 'dependency_ledger.json')
    assert ledger['latest_accepted_round'] == 1207 and ledger['stage_check_groups'] == 14
    assert ledger['goal_objective_sha256'] == d['goal_objective_sha256']
    assert ledger['latest_management_receipt'] == OUT.name
    print(json.dumps({'status': 'PASS', 'latest_round': 1207, 'stage_rounds': 7,
        'stage_check_groups': 14, 'old_GR_science_files_verified': ns,
        'old_GR_acceptance_bindings_verified': nb,
        'current_stage_science_files_verified': len(d['science_files_sha256']),
        'goal_status': 'active', 'goal_changed': False, 'stage_complete': False}))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists() and not (ROUND / 'acceptance.json').exists(), 'Already published'
    prev = read(PREV)
    ns, nb = prior.frozen(prev)
    prior.m.validate_map(prev['management_sha256'])
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
    w = r['native_finite_witness']['first_hit_witness']
    iw = independent['original_cold_amplified_witness']
    assert int(w['N_blocks_exact_decimal_integer']) == int(iw['blocks']) == 4375693525386771
    assert int(w['total_materials_exact_decimal_integer']) == int(iw['total_material_copies']) == 109392338134669275
    assert w['block_length'] == iw['block_length'] == 25
    assert F(w['trace_distance_strict_lower_rational']) > F('.283')
    assert F(iw['trace_distance_strict_lower_exact']) > F('.283')
    assert F(w['energy_operator_norm_bounds'][1]) < F('.0004237212511')
    assert F(w['conditional_trace_distance_bounds'][0]) < F(iw['conditional_asymmetry_exact']) < F(w['conditional_trace_distance_bounds'][1])
    assert independent['operator_norm_accumulation_control'] == {
        'blocks': 4, 'swap_every_block_norm': 4, 'first_hit_norm': 1}
    assert all(c['formula_matches_full_permutation'] and c['is_involution']
               for c in independent['first_hit_full_permutation_checks'])
    assert r['new_cognitive_axioms_adopted'] == 0
    assert not w['current_preparation_readout_or_control_resources_certified']
    files = ['proof.md', 'check.py', 'results.json', 'independent_check.py',
        'independent_results.json', 'independent_review.md', 'cognitive_scope_review.md', 'NEXT.md',
        '../research_note_1207.md', '../_shared/goal.md']
    hashes = {name: sha((ROUND / name).resolve()) for name in files}
    scope = ('给定有限H、满秩平稳参考、独立副本和全部近保裸能联合酉权限，'
        '跨有限组合统一半迹距连续性恰选择全实逆温度Gibbs族；非Gibbs在任意正总能容差下的扰动上确界为1。'
        '首次命中置换使概率放大而总能范数不累加；原冷参考有巨大但有限的数学证书。'
        '未采用新认知公理，未授予原QQ相应实际权限，未认证该证书的准备/读出资源。')
    write(ROUND / 'acceptance.json', {
        'schema': 'operational_thermodynamics_round_acceptance_v1', 'round': 1207,
        'date': '2026-10-14', 'status': 'ACCEPTED_WITH_SCOPE', 'check_groups': 2,
        'scope': scope, 'certified': {
            'full_rank_stationary_uniform_stability_iff_real_beta_Gibbs': True,
            'Gibbs_bound': 'min(1,sqrt(abs(beta)*delta/2))',
            'non_Gibbs_response_supremum_for_each_positive_delta': 1,
            'supremum_one_attained_at_finite_full_rank_size': False,
            'native_exact_response_at_delta_zero': 0,
            'first_hit_trace_formula': '(1-(1-a-b)^M)*abs(a-b)/(a+b)',
            'first_hit_total_energy_operator_error': 'abs(E_x-E_y), independent of M',
            'first_hit_x_binary_effect_attains_the_witness_trace_distance': True,
            'sufficient_comparison_family': 'first-hit permutations and first-x readouts for all required string pairs and finite block counts; not one fixed finite menu',
            'native_cold_block_length': 25,
            'native_cold_blocks': w['N_blocks_exact_decimal_integer'],
            'native_cold_total_materials': w['total_materials_exact_decimal_integer'],
            'native_cold_trace_distance_strict_lower': '.283',
            'native_total_energy_operator_error_upper': '.0004237212511',
            'finite_preparation_transport': 'min(1,2*epsilon_n+sqrt(abs(beta)*delta/2))',
            'full_operation_access_proven_for_native_model': False,
            'finite_witness_physical_resources_certified': False},
        'remaining_inputs': ['independent physical H, full rank and stationary reference role',
            'independently prepared identical copies and cross-all-finite-n uniformity',
            'sufficient permitted first-hit conditional permutations and joint discrimination',
            'full actual control, preparation, timing, memory, energy and readout resource contracts',
            'temperature sign/value, other conserved charges, same physical region/clock/flow identity'],
        'limits': ['material marginal stability does not return arbitrary purifying R',
            'operator energy error is not total implementation cost',
            'all-scale classification does not invalidate a declared finite effective task',
            'not a derivation of passivity, thermalization, geometry or GR'],
        'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'GR_hard_deliverables_completed': 0, 'stage_complete': False,
        'root_executed_both_codes': True,
        'review': 'Root proof audit; independent exact full-basis permutations, original cold-reference arithmetic, another native reference, scope and primary-literature audit; no empirical experiment.',
        'files_sha256': hashes})
    status = ('最新验收1207；新阶段7轮、14组。给定完整近保能操作权限，跨组合统一稳健性恰选择Gibbs形式；'
        '首次命中置换补上稀有配置的整体可辨放大。该条件仍是未采用的认知候选，实际比较权限和资源须另核；'
        '有限参考仍有效，温度正号及同一几何身份未闭合。')
    ledger = read(HERE / 'dependency_ledger.json')
    req = ledger['stage_requirement_status']
    req.update(thermodynamic_structural_basis=(
        '1207: uniform tensor-composition stability under the declared full near-energy-preserving '
        'unitary menu iff real-beta Gibbs; first-hit amplifier proves maximal non-Gibbs response. '
        'Candidate not adopted; actual permitted comparisons and complete resource implementation remain open'),
        stage_synthesis='reference-selection theorem obtained conditionally; cognitive adoption/actual permissions and common region/clock/flow identity remain open',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1207, accepted_stage_rounds=list(range(1201, 1208)),
        stage_check_groups=14, current_direction_document='../1207/NEXT.md', current_work=status,
        latest_management_receipt=OUT.name, stage_requirement_status=req)
    write(HERE / 'dependency_ledger.json', ledger)
    p = STAGE / 'README.md'
    parts = p.read_text(encoding='utf-8').split('\n\n')
    parts[1] = status
    text = '\n\n'.join(parts).replace('_shared/progress_checks_1206.json', '_shared/' + OUT.name)
    text = text.replace('[当前下一项](1206/NEXT.md)', '[1206原安排](1206/NEXT.md)')
    text = text.replace('## 接续与边界', '- [1207报告](research_note_1207.md)：跨组合稳健校准、Gibbs分类和首次命中放大；[验收](1207/acceptance.json)、[当前下一项](1207/NEXT.md)。\n\n## 接续与边界', 1)
    text = text.replace('[1206原安排](1206/NEXT.md)接续1206', '[当前下一项](1207/NEXT.md)接续1207')
    p.write_text(text, encoding='utf-8')
    p = STAGE / '文件索引.md'
    text = p.read_text(encoding='utf-8').replace('最新验收1206；6轮、12组', '最新验收1207；7轮、14组')
    text = text.replace('## 管理与接续', '|1207|[跨组合稳健校准与Gibbs参考选择](research_note_1207.md)|[证明](1207/proof.md)、[代码](1207/check.py)、[结果](1207/results.json)、[独立核验](1207/independent_check.py)、[独审](1207/independent_review.md)、[范围](1207/cognitive_scope_review.md)、[验收](1207/acceptance.json)|\n\n## 管理与接续', 1)
    text = re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)', r'\1\n', text)
    text = text.replace('[当前下一项](1206/NEXT.md)', '[当前下一项](1207/NEXT.md)')
    text = text.replace('_shared/progress_checks_1206.json', '_shared/' + OUT.name)
    p.write_text(text, encoding='utf-8')
    oldledger = read(OLD / '_shared/dependency_ledger.json')
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(OLD / '_shared/dependency_ledger.json', oldledger)
    p = OLD / '_shared/navigation_notice.md'
    prefix, rest = p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    prefix = re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)', '**当前状态：** ' + status +
        ' [1207报告](archive_1201_/research_note_1207.md) · [下一项](archive_1201_/1207/NEXT.md) · '
        '[当前核验](archive_1201_/_shared/' + OUT.name + ')。旧GR阶段40轮、130组保留；目标继续进行。', prefix, flags=re.S)
    p.write_text(prefix + '以下保留前阶段的历史状态' + rest, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navverify['status'] == 'PASS'
    write(OUT, {'status': 'PUBLISHING'})
    documents = [STAGE / 'README.md', STAGE / '文件索引.md', STAGE / 'research_note_1207.md',
        *ROUND.glob('*.md'), *[BASE / name for name in
        ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']]]
    nlinks = links(documents)
    for p in documents:
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in p.read_text(encoding='utf-8'))
    science = dict(prev['science_files_sha256'])
    science.update({(ROUND / name).resolve().relative_to(ROOT).as_posix(): h for name, h in hashes.items()})
    management = {name: sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    accs = dict(prev['acceptance_sha256'])
    accs['1207'] = sha(ROUND / 'acceptance.json')
    write(OUT, {'schema': 'thermodynamics_progress_1207_v1', 'date': '2026-10-14', 'status': 'PASS',
        'accepted_rounds': list(range(1201, 1208)), 'stage_check_groups': 14,
        'application_goal_status': 'active', 'stage_complete': False, 'goal_changed': False,
        'goal_objective_sha256': ledger['goal_objective_sha256'], 'acceptance_sha256': accs,
        'science_files_sha256': science, 'management_sha256': management,
        'screen_files_sha256': prev['screen_files_sha256'], 'audit_files_sha256': prev['audit_files_sha256'],
        'prior_progress_path': PREV.relative_to(ROOT).as_posix(), 'prior_progress_sha256': sha(PREV),
        'previous_GR_science_files_verified': ns, 'previous_GR_acceptance_bindings_verified': nb,
        'navigation_apply': navapply, 'navigation_check': navcheck, 'navigation_verify': navverify,
        'markdown_files_checked': len(documents), 'local_links_checked': nlinks, 'missing_links': [],
        'stage_requirement_status': req, 'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'GR_hard_deliverables_completed': 0, 'science_program_outputs': {'author': author,
            'independent_all_checks_passed': independent['all_checks_passed'],
            'independent_result_sha256': sha(ROUND / 'independent_results.json')},
        'previous_goal_turn_classification': 'progress: 1206 accepted and frozen with scope and next task',
        'this_goal_turn_result': 'progress: conditional Gibbs selection from composition-robust calibration; first-hit amplification closes rare-probability gap; actual permissions and geometry remain open'})
    verify()


if __name__ == '__main__':
    main()
