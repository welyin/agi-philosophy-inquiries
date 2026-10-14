"""Publish the fixed preparation bridge only after scoped independent review."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util, json, re, subprocess, sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1201_'
HERE = STAGE / '_shared'
ROUND = STAGE / '1205'
OLD = BASE / 'archive_1161_'
OUT = HERE / 'progress_checks_1205.json'
s = importlib.util.spec_from_file_location('publish1201_tools', HERE / 'publish_1201.py')
m = importlib.util.module_from_spec(s)
s.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links
old_check, validate_map = m.old_check, m.validate_map


def execute_science(path):
    result = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(path)],
                            capture_output=True, text=True, encoding='utf-8')
    assert result.returncode == 0, result.stderr + result.stdout
    return json.loads(result.stdout)


def verify():
    d = read(OUT)
    old, n = old_check(read(HERE / 'startup_checks.json'))
    for key in ['science_files_sha256', 'management_sha256', 'screen_files_sha256']:
        validate_map(d[key])
    for number, digest in d['acceptance_sha256'].items():
        p = STAGE / number / 'acceptance.json'
        assert sha(p) == digest
        for name, h in read(p)['files_sha256'].items():
            assert sha((p.parent / name).resolve()) == h
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    print(json.dumps({'status': 'PASS', 'latest_round': 1205, 'new_stage_rounds': 5,
        'stage_check_groups': 10, 'old_GR_science_files_verified': len(old['science_files_sha256']),
        'old_GR_bindings_verified': n, 'current_stage_science_files': len(d['science_files_sha256']),
        'goal_status': 'active', 'stage_complete': False}))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists() and not (ROUND / 'acceptance.json').exists()
    pp = HERE / 'progress_checks_1204.json'
    prev = read(pp)
    for key in ['science_files_sha256', 'management_sha256', 'screen_files_sha256']:
        validate_map(prev[key])
    for number, h in prev['acceptance_sha256'].items():
        assert sha(STAGE / number / 'acceptance.json') == h
    old, n = old_check(read(HERE / 'startup_checks.json'))
    review = (ROUND / 'independent_review.md').read_text(encoding='utf-8')
    assert '限定签核' in review, 'Independent final review required'
    files = ['proof.md', 'check.py', 'results.json', 'independent_check.py',
             'independent_results.json', 'independent_review.md', 'NEXT.md',
             '../research_note_1205.md', '../_shared/goal.md']
    for name in files:
        assert (ROUND / name).is_file(), name
    author_output = execute_science(ROUND / 'check.py')
    independent_output = execute_science(ROUND / 'independent_check.py')
    assert author_output['status'] == 'AUTHOR_CHECKS_PASS'
    assert independent_output['all_checks_passed'] is True
    r = read(ROUND / 'results.json')
    independent_saved = read(ROUND / 'independent_results.json')
    assert independent_output == independent_saved
    assert independent_saved['fixed_tables']['dimensions'] == r['fixed_program']['dimensions']
    assert independent_saved['fixed_tables']['type_counts'] == r['fixed_program']['types']
    assert independent_saved['resource_bounds']['reversible_gates_upper'] == r['fixed_program']['reversible_gate_bound']
    assert F(independent_saved['whole_path_clock']['epsilon_clock_exact']) < F(r['normal_clock']['epsilon_clock_upper'])
    assert F(independent_saved['whole_path_clock']['clock_entropy_rational_upper']) < F(r['normal_clock']['conditional_clock_entropy_upper'])
    for t in [600, 300]:
        assert F(independent_saved['native_flow_transports'][str(t)]['actual_Q_radius_exact']) <= F(r['native_flow_transport'][f't{t}_lambda1']['new_actual_energy_radius_upper'])
    assert r['fixed_program']['actual_word_count_two_pools'] <= 278276
    assert r['fixed_program']['reversible_gate_bound'] == 14006089679648
    assert r['fixed_program']['all_seed_garbage_failure_stock_and_reference_kept']
    assert F(r['normal_clock']['epsilon_clock_upper']) < F(4, 10**15)
    assert r['normal_clock']['same_conserved_success_event_at_both_contact_sections']
    assert F(r['normal_clock']['conditional_clock_entropy_upper']) < F(5, 10**16)
    assert F(r['thermodynamic_transport']['finite_entropy_and_energy_extra_upper']) < F(2, 10**12)
    assert F(r['native_flow_transport']['t300_lambda1']['Q_dir_lower']) > F('.02756774')
    assert F(r['native_flow_transport']['t600_lambda1']['Q_dir_lower']) > F('.07314762')
    hashes = {name: sha((ROUND / name).resolve()) for name in files}
    scope = ('固定N=100两能窗的保记录可逆程序、1000材料库存、同一正常正能钟与原QQ/实体K接续；'
             '按守恒种子成功事件S运输两个截面，末计算旗标单独计误差。给冷热全支的有限流、'
             '完整能源和有限熵修正；具体硬件/初态仍输入，未扩到所有旧菜单或共同几何。')
    acc = {'schema': 'operational_thermodynamics_round_acceptance_v1', 'round': 1205,
        'date': '2026-10-14', 'status': 'ACCEPTED_WITH_SCOPE', 'check_groups': 2,
        'scope': scope, 'certified': {
            'fixed_program': r['fixed_program'], 'normal_clock': r['normal_clock'],
            'thermodynamic_transport': r['thermodynamic_transport'],
            'native_flow_transport': r['native_flow_transport'],
            'full_material_energy_width_bound': '1094',
            'separate_final_flag_and_conserved_success_history': True,
            'same_preparation_H_covers_all_older_menus': False},
        'remaining_inputs': ['fixed h/QQ and degenerate logic registers with specified gates',
            'independent seeds, pure work bits, known-level stocks and Gaussian clock',
            'terminal classification/calibration and independent repeated preparations',
            'physical CO/SR region, signal and metric identity'],
        'limits': r['not_claimed'], 'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'GR_hard_deliverables_completed': 0, 'stage_complete': False,
        'root_executed_both_codes': True,
        'review': 'Root proof and evidence audit plus independent circuit, integer, clock and conditional-history checks; no physical experiment.',
        'files_sha256': hashes}
    write(ROUND / 'acceptance.json', acc)
    status = ('最新验收1205；新阶段5轮、10组。固定两能窗准备的可逆门、工作位、库存与正能钟已接原QQ和实体记录；'
              '守恒成功历史与末旗标分开，完整能源和有限熵修正有上界。下一项核全部原始K记录能否消除成功位末读依赖，再完成阶段条件定理和论文审计。')
    req = dict(read(HERE / 'dependency_ledger.json')['stage_requirement_status'])
    req.update(independently_generated_stable_reference='fixed independent energy-window references with task-specific preparation resource certificate in 1205',
        internal_control_clock_record_preparation_resource_closure='1205 fixed cold/hot chain: finite circuit, full stocks/records, autonomous clock and end-to-end error; explicit hardware/initial inputs',
        stage_synthesis='fixed instance complete; test raw-pointer contrast without success-flag readout, then final theorem/paper and requirement audit',
        stage_deliverables_all_completed=False)
    ledger = read(HERE / 'dependency_ledger.json')
    ledger.update(latest_accepted_round=1205, accepted_stage_rounds=list(range(1201, 1206)),
        stage_check_groups=10, current_direction_document='../1205/NEXT.md', current_work=status,
        latest_management_receipt='progress_checks_1205.json', stage_requirement_status=req)
    write(HERE / 'dependency_ledger.json', ledger)
    p = STAGE / 'README.md'
    txt = p.read_text(encoding='utf-8')
    paragraphs = txt.split('\n\n'); paragraphs[1] = status; txt = '\n\n'.join(paragraphs)
    txt = txt.replace('## 接续与边界', '- [1205报告](research_note_1205.md)：固定准备资源、守恒成功历史和全链接续；[验收](1205/acceptance.json)、[下一项](1205/NEXT.md)。\n\n## 接续与边界', 1)
    txt = txt.replace('_shared/progress_checks_1204.json', '_shared/progress_checks_1205.json')
    txt = txt.replace('1201—1204已接参考、内部接触／记录及有限实际流；固定准备资源上界和阶段综合仍待完成',
                      '1201—1205已接固定参考准备、内部接触／记录及有限实际流；阶段条件综合与最终交付审计仍待完成')
    p.write_text(txt, encoding='utf-8')
    p = STAGE / '文件索引.md'; txt = p.read_text(encoding='utf-8')
    txt = txt.replace('最新验收1204；4轮、8组', '最新验收1205；5轮、10组')
    txt = txt.replace('## 管理与接续', '|1205|[固定参考准备、守恒成功历史与自主接续](research_note_1205.md)|[证明](1205/proof.md)、[代码](1205/check.py)、[结果](1205/results.json)、[独立核验](1205/independent_check.py)、[审阅](1205/independent_review.md)、[验收](1205/acceptance.json)|\n\n## 管理与接续', 1)
    txt = re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)', r'\1\n', txt)
    txt = txt.replace('[下一项](1204/NEXT.md)', '[下一项](1205/NEXT.md)').replace('_shared/progress_checks_1204.json', '_shared/progress_checks_1205.json')
    p.write_text(txt, encoding='utf-8')
    d = read(OLD / '_shared/dependency_ledger.json')
    d['latest_management_receipt'] = '../../archive_1201_/_shared/progress_checks_1205.json'
    write(OLD / '_shared/dependency_ledger.json', d)
    p = OLD / '_shared/navigation_notice.md'; txt = p.read_text(encoding='utf-8')
    prefix, rest = txt.split('以下保留前阶段的历史状态', 1)
    prefix = re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)', '**当前状态：** ' + status +
        ' [1205报告](archive_1201_/research_note_1205.md) · [验收](archive_1201_/1205/acceptance.json) · [当前核验](archive_1201_/_shared/progress_checks_1205.json)。旧GR阶段40轮、130组和冻结稿均保留；目标继续进行。', prefix, flags=re.S)
    p.write_text(prefix + '以下保留前阶段的历史状态' + rest, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check'); assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation'); assert navverify['status'] == 'PASS'
    write(OUT, {'status': 'PUBLISHING'})
    paths = [STAGE / 'README.md', STAGE / '文件索引.md', STAGE / 'research_note_1205.md', *ROUND.glob('*.md'),
             *[BASE / name for name in ['README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md']]]
    nlinks = links(paths)
    science = dict(prev['science_files_sha256'])
    science.update({(ROUND / name).resolve().relative_to(ROOT).as_posix(): h for name, h in hashes.items()})
    management = {name: sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    accs = dict(prev['acceptance_sha256']); accs['1205'] = sha(ROUND / 'acceptance.json')
    write(OUT, {'schema': 'thermodynamics_progress_1205_v1', 'date': '2026-10-14', 'status': 'PASS',
        'accepted_rounds': list(range(1201, 1206)), 'stage_check_groups': 10,
        'application_goal_status': 'active', 'stage_complete': False,
        'goal_objective_sha256': ledger['goal_objective_sha256'], 'acceptance_sha256': accs,
        'science_files_sha256': science, 'management_sha256': management,
        'screen_files_sha256': dict(prev['screen_files_sha256']),
        'prior_progress_path': pp.relative_to(ROOT).as_posix(), 'prior_progress_sha256': sha(pp),
        'previous_GR_science_files_verified': len(old['science_files_sha256']), 'previous_GR_bindings': n,
        'navigation_apply': navapply, 'navigation_check': navcheck, 'navigation_verify': navverify,
        'markdown_files_checked': len(paths), 'local_links_checked': nlinks, 'missing_links': [],
        'stage_requirement_status': req, 'new_adopted_axioms': 0, 'new_empirical_calibration': 0,
        'GR_hard_deliverables_completed': 0,
        'science_program_outputs': {'author': author_output,
            'independent_all_checks_passed': independent_output['all_checks_passed'],
            'independent_result_sha256': sha(ROUND / 'independent_results.json')},
        'this_goal_turn_result': 'progress: fixed preparation resource and autonomous full-history bridge; final stage synthesis pending'})
    verify()


if __name__ == '__main__':
    main()
