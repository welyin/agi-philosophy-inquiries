"""Accept two scoped four-law interfaces without changing the active objective."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'research_cognition_physics'
STAGE = BASE/'archive_1201_'
HERE = STAGE/'_shared'
OLD = BASE/'archive_1161_'
PREV = HERE/'progress_checks_1211.json'
OUT = HERE/'progress_checks_1213.json'
spec = importlib.util.spec_from_file_location('prior_thermo_tools', HERE/'publish_1206.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links
NUMBERS = (1212, 1213)
REVIEWED = ['proof.md', 'check.py', 'results.json', 'independent_check.py',
            'independent_results.json', 'NEXT.md']
AUDIT = HERE/'four_laws_dependency_after1211.md'


def verify():
    data = read(OUT); old, bindings = m.frozen(data)
    m.m.validate_map(data['management_sha256'])
    assert sha(PREV) == data['prior_progress_sha256']
    ledger = read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round'] == 1213 and ledger['stage_check_groups'] == 26
    assert ledger['latest_management_receipt'] == OUT.name
    assert ledger['goal_objective_sha256'] == data['goal_objective_sha256']
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert run(OLD/'_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT/'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    print(json.dumps(dict(status='PASS', latest_round=1213, stage_rounds=13,
        stage_check_groups=26, old_GR_science_files_verified=old,
        old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(data['science_files_sha256']),
        goal_changed=False, application_goal_status='active', stage_complete=False)))


def main():
    if '--verify' in sys.argv: return verify()
    assert not OUT.exists(), 'Already published'
    previous = read(PREV)
    old, bindings = m.frozen(previous)
    m.m.validate_map(previous['management_sha256'])
    outputs, additions, new_acceptances = {}, {}, {}
    scopes = {
        1212: '原有限参考的完整机器/交互/关联工作价值账、实际能源窗口的严格有限容量损失，以及全尝试记录未归还和300接触段钟工作价值为负；限定循环必要界，不认证已实现热机。',
        1213: '原同h材料、正能连续谱辅助与纯零H记录，在新增全局保能耦合下可有限端点精确重置且完整保未知SR；只否定有限能源矩/最坏能源转移/完整内账足以推出端点不可达的弱命题。'}
    for n in NUMBERS:
        folder = STAGE/str(n)
        assert not (folder/'acceptance.json').exists()
        names = REVIEWED+[f'../research_note_{n}.md']
        review = (folder/'independent_review.md').read_text(encoding='utf-8')
        assert '限定签核通过' in review
        for name in names: assert sha((folder/name).resolve()) in review, (n, name)
        for name in ['proof.md', 'NEXT.md', f'../research_note_{n}.md']:
            value = (folder/name).read_text(encoding='utf-8')
            assert not any(x in value for x in (r'\[', r'\]', r'\(', r'\)')), (n, name)
        author = run(folder/'check.py')
        independent = run(folder/'independent_check.py')
        for output in (author, independent):
            assert output['status'] == 'PASS'
            assert output['all_checks_passed'] and output['results_match']
        for obj in (read(folder/'results.json'), read(folder/'independent_results.json')):
            for name, digest in obj.get('source_sha256', {}).items():
                options = [STAGE/name, BASE/name, ROOT/name]
                source = next((p for p in options if p.is_file()), None)
                assert source and sha(source) == digest, (n, name)
        ir = read(folder/'independent_results.json')
        assert ir['all_checks_passed']
        assert ir['independent_code_sha256'] == sha(folder/'independent_check.py')
        outputs[str(n)] = dict(author=author, independent=independent)

    a12 = read(STAGE/'1212/results.json')
    for row in a12['finite_flow_bounds']:
        assert F(row['hot_only_work_value_gap_lower']) > 0
    cycles = {x['t']:x for x in a12['actual_factory_cycle_obstruction']}
    assert F(cycles[300]['contact_clock_work_value_change_interval'][1]) < 0
    assert F(cycles[600]['raw_K_nonblank_lower']) > F('.653')
    assert not a12['scope']['original_QQ_heat_engine_implemented']
    assert not a12['scope']['unconditional_four_laws_derived']
    i12 = read(STAGE/'1212/independent_results.json')
    assert i12['generic_complete_balance_cases'] == 24
    assert i12['max_balance_error'] < 1e-11
    a13 = read(STAGE/'1213/results.json')
    assert a13['check_groups'] == 8 and a13['all_checks_passed']
    assert all(g['pass'] for g in a13['groups'].values())
    boundaries = a13['groups']['threshold_domain_and_claim_boundaries']
    assert boundaries['weak_energy_only_endpoint_unattainability_refuted']
    assert not boundaries['old_QQ_or_local_SR_permission_proved']
    assert not boundaries['all_resources_finite_claim']
    assert not boundaries['third_law_or_Masanes_Oppenheim_refuted']
    assert not boundaries['shifted_Gaussian_in_DX_for_positive_gap']
    assert a13['groups']['unknown_SR_and_complete_record']['excited_output_residual'] == 0

    for n in NUMBERS:
        folder = STAGE/str(n)
        names = REVIEWED+['independent_review.md', f'../research_note_{n}.md',
                          '../_shared/goal_four_laws.md']
        hashes = {name:sha((folder/name).resolve()) for name in names}
        certified = (dict(complete_correlation_interaction_and_machine_balance=True,
            native_finite_reservoir_deficit_positive=True, original_program_noncyclic=True,
            actual_300_clock_work_value_negative=True, actual_600_clock_gain_sign_certified=False,
            free_energy_equals_extractable_work=False, all_cycle_completions_excluded=False,
            conditional_Clausius_to_Kelvin=True, reverse_equivalence_proved=False)
            if n == 1212 else dict(positive_continuous_spectrum_auxiliary=True,
            same_h_ground_stock_not_supplied=True, pure_record_consumed=True,
            global_unitary_preserves_unknown_SR=True, bounded_energy_transfer=True,
            finite_energy_first_second_moments=True, Gaussian_final_position_second_moment_finite=False,
            new_global_interaction_permission=True, native_local_clock_implementation=False,
            permanent_or_released_cold_product=False, all_cognitive_axioms_countermodel=False,
            general_third_law_disproved=False))
        acceptance = dict(schema='operational_thermodynamics_round_acceptance_v1',
            round=n, date='2026-10-15', status='ACCEPTED_WITH_SCOPE', check_groups=2,
            scope=scopes[n], certified=certified, root_executed_both_codes=True,
            new_adopted_axioms=0, new_empirical_calibration=0,
            GR_hard_deliverables_completed=0, stage_complete=False,
            files_sha256=hashes,
            limits=['Existing cognitive axioms and the active four-law goal are unchanged.',
                    'The finite quantum checks do not replace analytic general statements.',
                    'Reference roles, actual permissible processes and resource supply remain scoped.'])
        write(folder/'acceptance.json', acceptance)
        new_acceptances[str(n)] = sha(folder/'acceptance.json')
        additions.update({(folder/name).resolve().relative_to(ROOT).as_posix():digest for name,digest in hashes.items()})

    status = ('最新验收1213；本阶段13轮、26组。1212认证原有限参考的循环工作价值损失、'
              '全尝试记录未归还及300接触段钟工作价值为负；1213给正能连续谱与纯记录的有限端点重置，'
              '严格披露新增非局域耦合及Gaussian输出位置矩问题，只否定弱能源条件推出端点不可达。'
              '四律尚未全部推导；下一步统一热参考、工作/控制与可消耗纯度的资源角色，未新增采用公理。')
    ledger = read(HERE/'dependency_ledger.json')
    req = ledger['stage_requirement_status']
    req.update(first_law='D-P already supplies the autonomous time/complete-energy identity; 1212 consolidates actual total conservation and scoped energy/work-value/heat roles, without inventing a new time axiom',
        second_law='1212: complete cyclic-work necessary bounds, finite-reservoir deficits and actual machine non-return; conditional Clausius-to-Kelvin, no reverse equivalence or native heat-engine certification',
        third_law='1211-1213: cooling resource boundaries and a positive-energy continuous-spectrum endpoint reset; weak energy budgets do not enforce unattainability, full thermal-machine roles and time/resource supply remain open',
        stage_synthesis='Four-laws paper and per-law common-scope completion audit remain pending',
        stage_deliverables_all_completed=False)
    ledger.update(date='2026-10-15', latest_accepted_round=1213,
        accepted_stage_rounds=list(range(1201,1214)), stage_check_groups=26,
        current_direction_document='../1213/NEXT.md', current_scope_audit=AUDIT.name,
        current_work=status, latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json', ledger)
    p = STAGE/'README.md'; value = p.read_text(encoding='utf-8'); parts=value.split('\n\n'); parts[1]=status
    value='\n\n'.join(parts).replace('_shared/'+PREV.name, '_shared/'+OUT.name)
    rows = ('- [1212报告](research_note_1212.md)：有限参考循环工作价值与实际工厂归还；[验收](1212/acceptance.json)。\n\n'
            '- [1213报告](research_note_1213.md)：正能连续谱、纯记录与端点重置的资源边界；[验收](1213/acceptance.json)、[下一项](1213/NEXT.md)。\n\n')
    value=value.replace('## 接续与边界', rows+'## 接续与边界', 1)
    value+='\n[四律认知依赖审计](_shared/'+AUDIT.name+')：时间—能源身份已采用；热参考角色与资源供应条件继续研究。\n'
    p.write_text(value, encoding='utf-8')
    p=STAGE/'文件索引.md'; value=p.read_text(encoding='utf-8')
    value=value.replace('2026-10-14。最新验收1211；11轮、22组', '2026-10-15。最新验收1213；13轮、26组')
    value=value.replace('_shared/'+PREV.name, '_shared/'+OUT.name)
    rows=''.join(f'|{n}|[{"循环工作价值与实际归还" if n==1212 else "正能连续谱与纯记录重置"}](research_note_{n}.md)|[证明]({n}/proof.md)、[代码]({n}/check.py)、[结果]({n}/results.json)、[独审]({n}/independent_review.md)、[验收]({n}/acceptance.json)|\n' for n in NUMBERS)
    value=value.replace('\n## 管理与接续', '\n'+rows+'\n## 管理与接续', 1)
    value=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)', r'\1\n', value)
    value+='\n- [四律认知依赖审计](_shared/'+AUDIT.name+')\n'
    p.write_text(value, encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json', oldledger)
    p=OLD/'_shared/navigation_notice.md'
    prefix, history=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)', '**当前状态：** '+status+
        ' [1212报告](archive_1201_/research_note_1212.md) · [1213报告](archive_1201_/research_note_1213.md) · '
        '[当前核验](archive_1201_/_shared/'+OUT.name+')。旧目标均未追溯标记完成。', prefix, flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+history, encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py', '--apply')
    navcheck=run(OLD/'_shared/update_navigation.py', '--check'); assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py', '--verify-navigation'); assert navverify['status']=='PASS'
    write(OUT, dict(status='PUBLISHING'))
    docs=[STAGE/'README.md', STAGE/'文件索引.md', AUDIT]
    for n in NUMBERS: docs += [STAGE/f'research_note_{n}.md', *(STAGE/str(n)).glob('*.md')]
    docs += [BASE/name for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks=links(docs)
    for p in docs:
        value=p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in value)
    management={name:sha(ROOT/name) for name in previous['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    write(OUT, dict(schema='thermodynamics_progress_1213_v1', date='2026-10-15', status='PASS',
        accepted_rounds=list(range(1201,1214)), stage_check_groups=26,
        goal_changed=False, application_goal_status='active', stage_complete=False,
        goal_objective_sha256=previous['goal_objective_sha256'],
        acceptance_sha256={**previous['acceptance_sha256'], **new_acceptances},
        science_files_sha256={**previous['science_files_sha256'], **additions},
        management_sha256=management, screen_files_sha256=previous['screen_files_sha256'],
        audit_files_sha256={**previous['audit_files_sha256'], AUDIT.relative_to(ROOT).as_posix():sha(AUDIT)},
        prior_progress_path=PREV.relative_to(ROOT).as_posix(), prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old, old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply, navigation_check=navcheck, navigation_verify=navverify,
        local_links_checked=nlinks, missing_links=[], stage_requirement_status=req,
        new_adopted_axioms=0, new_empirical_calibration=0, GR_hard_deliverables_completed=0,
        science_program_outputs=outputs,
        prior_goal_turn_classification='progress: conditional equilibrium and exact cooling-resource boundaries',
        this_goal_turn_result='progress: complete cyclic-work bound and scoped positive-energy endpoint reset with independent checks; four-laws goal active'))
    verify()


if __name__ == '__main__': main()
