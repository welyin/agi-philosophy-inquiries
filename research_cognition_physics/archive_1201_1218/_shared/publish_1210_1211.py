"""Publish two independently checked four-law interface rounds; keep goal open."""
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
PREV = HERE/'four_laws_goal_start_checks.json'
OUT = HERE/'progress_checks_1211.json'
spec = importlib.util.spec_from_file_location('frozen_thermo_tools', HERE/'publish_1206.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links
NUMBERS = (1210, 1211)
REVIEWED = ['proof.md', 'check.py', 'results.json', 'independent_check.py',
            'independent_results.json', 'NEXT.md']


def verify():
    data = read(OUT); old, bindings = m.frozen(data)
    m.m.validate_map(data['management_sha256'])
    assert sha(PREV) == data['prior_progress_sha256']
    ledger = read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round'] == 1211 and ledger['stage_check_groups'] == 22
    assert ledger['latest_management_receipt'] == OUT.name
    assert ledger['goal_objective_sha256'] == data['goal_objective_sha256']
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert run(OLD/'_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT/'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    print(json.dumps(dict(status='PASS', latest_round=1211, stage_rounds=11,
        stage_check_groups=22, old_GR_science_files_verified=old,
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
        1210: '复用1209联合三预算参考角色，在非标量物理H、统一能源单位和共同充分比较权限下，联合资格等价于唯一共同beta，形成精确平衡等价关系；伙伴依赖端口和标量中介删除检验，以及原M600有限温标区间。',
        1211: '原非共振四能级同h材料与零能辅助、精确保裸总能的有限比较类中，无基态部门给严格冷却上界；纯记录保信息路由达到1减无基态概率。区分纯记录、已备同h基态和有限混合库存，低温熵与不可达版本分别处理。'
    }
    # Execute and inspect both full research checks before mutating navigation.
    for n in NUMBERS:
        folder = STAGE/str(n)
        assert not (folder/'acceptance.json').exists()
        names = REVIEWED + [f'../research_note_{n}.md']
        review = (folder/'independent_review.md').read_text(encoding='utf-8')
        assert '限定签核通过' in review
        for name in names: assert sha((folder/name).resolve()) in review, (n, name)
        author = run(folder/'check.py')
        independent = run(folder/'independent_check.py')
        assert author['status'].startswith('PASS') and independent['status'] == 'PASS'
        assert independent['all_checks_passed'] and independent['results_match']
        saved = read(folder/'results.json')
        if n == 1210:
            assert saved == author
        else:
            assert author['all_checks_passed'] and author['results_match']
            assert saved['status'] == 'AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE'
        independent_saved = read(folder/'independent_results.json')
        assert independent_saved['all_checks_passed']
        assert independent_saved['independent_code_sha256'] == sha(folder/'independent_check.py')
        for obj in (saved, independent_saved):
            for path, digest in obj.get('source_sha256', {}).items():
                source = STAGE/path if (STAGE/path).exists() else BASE/path
                assert sha(source.resolve()) == digest, (n, path)
        outputs[str(n)] = dict(author_status=author['status'], independent=independent)
    a = read(STAGE/'1210/results.json')
    assert F(a['native_permission_readout_bridge']['single_event_true_gap_lower']) > F('.1161')
    controls = a['actual_readout_temperature_interface']['controls']
    assert F(controls['hot']['all_valid_data_beta_interval_outer'][1]) < F(controls['cold']['all_valid_data_beta_interval_outer'][0])
    i = read(STAGE/'1211/independent_results.json')
    assert len(i['complete_permutation_cases']) == 4
    assert i['correlated_full_rank_marginals_but_zero_no_ground_probability']
    assert not i['qualifications']['general_third_law_derived']
    a11 = read(STAGE/'1211/results.json')
    for pool in a11['old_correlated_window_pools']:
        independent_pool = i['native_correlated_output_pool_bounds'][pool['reference']]
        assert pool['no_ground_word_count'] == independent_pool['no_ground_words']
        assert F(pool['ideal_correlated_pool_failure_exact']) == F(independent_pool['unavoidable_failure_exact'])
        assert not pool['positive_failure_for_actual_clock_prepared_state_claimed']

    for n in NUMBERS:
        folder = STAGE/str(n)
        names = REVIEWED + ['independent_review.md', f'../research_note_{n}.md',
                            '../_shared/goal_four_laws.md']
        hashes = {name: sha((folder/name).resolve()) for name in names}
        acceptance = dict(schema='operational_thermodynamics_round_acceptance_v1',
            round=n, date='2026-10-14', status='ACCEPTED_WITH_SCOPE', check_groups=2,
            scope=scopes[n], root_executed_both_codes=True,
            new_adopted_axioms=0, new_empirical_calibration=0,
            GR_hard_deliverables_completed=0, stage_complete=False,
            files_sha256=hashes,
            limits=['Mathematical comparison permission is not an actual native clock implementation.',
                    'No new cognitive axiom adopted; no claim that all four laws are derived.',
                    'Earlier finite full-resource process keeps its original qualified scope.'])
        if n == 1210:
            acceptance['certified'] = dict(conditional_equilibrium_equivalence=True,
                non_scalar_energy_required=True, common_comparison_permission_required=True,
                finite_data_prove_all_scale_role=False, ordinary_positive_temperature_selected=False,
                fixed_tolerance_relation_transitive=False,
                beta_interval_envelopes={'hot':['.100','.578'],'cold':['.891','1.585']})
        else:
            acceptance['certified'] = dict(optimal_routing_probability='1-Tr(P_no_ground Omega)',
                same_h_materials_and_zero_energy_auxiliaries_required=True,
                pure_record_resource_consumed=True, arbitrary_material_reference_preserved_globally=True,
                full_rank_marginals_alone_suffice_for_unattainability=False,
                native_actual_clock_exact_cooling_certified=False,
                general_finite_time_third_law_proved=False, nernst_and_unattainability_distinguished=True)
        write(folder/'acceptance.json', acceptance)
        new_acceptances[str(n)] = sha(folder/'acceptance.json')
        additions.update({(folder/name).resolve().relative_to(ROOT).as_posix(): digest
                          for name, digest in hashes.items()})

    status = ('最新验收1211；本阶段11轮、22组。1210给共同参考角色下的条件性第零律、权限删除检验及有限温标读出；'
        '1211给原同h库存的精确冷却资源界和保记录最优路由，区分已有基态、空白记录与有限混合库存。'
        '四律尚未全部推导，参考角色、真实循环权限及一般第三律资源条件继续研究；未新增采用公理。')
    ledger = read(HERE/'dependency_ledger.json')
    req = ledger['stage_requirement_status']
    req.update(zeroth_law='1210: conditional equilibrium equivalence and identifiable beta under the unadopted 1209 role and a common separating comparison menu; native finite readout connected, not a certification of all-n role',
        third_law='1211: sharp same-h stock/zero-energy-record cooling boundary, coherent finite routing, and separate finite-spectrum low-T entropy scope; general cooling resources/time conditions remain open',
        second_law='Next: actual clock/work value, full machine cyclicity, finite-reservoir losses and the permitted relation between Clausius and cyclic work bounds; existing entropy identities reused',
        stage_synthesis='Four-laws paper and per-law common-scope completion audit remain pending',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1211, accepted_stage_rounds=list(range(1201,1212)),
        stage_check_groups=22, current_direction_document='../1211/NEXT.md',
        current_work=status, latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json', ledger)
    p = STAGE/'README.md'; text = p.read_text(encoding='utf-8'); parts=text.split('\n\n'); parts[1]=status
    text='\n\n'.join(parts).replace('_shared/four_laws_goal_start_checks.json','_shared/'+OUT.name)
    new = ('- [1210报告](research_note_1210.md)：共同温标、参考资格与比较权限；[验收](1210/acceptance.json)。\n\n'
           '- [1211报告](research_note_1211.md)：有限库存、保记录路由与第三定律边界；[验收](1211/acceptance.json)、[下一项](1211/NEXT.md)。\n\n')
    text=text.replace('## 接续与边界',new+'## 接续与边界',1)
    p.write_text(text,encoding='utf-8')
    p=STAGE/'文件索引.md'; text=p.read_text(encoding='utf-8')
    text=text.replace('最新验收1209；9轮、18组','最新验收1211；11轮、22组')
    text=text.replace('_shared/four_laws_goal_start_checks.json','_shared/'+OUT.name)
    rows=''.join(f'|{n}|[{"共同温标与比较权限" if n==1210 else "有限库存与第三定律资源边界"}](research_note_{n}.md)|[证明]({n}/proof.md)、[代码]({n}/check.py)、[结果]({n}/results.json)、[独审]({n}/independent_review.md)、[验收]({n}/acceptance.json)|\n' for n in NUMBERS)
    text=text.replace('\n## 管理与接续','\n'+rows+'\n## 管理与接续',1)
    text=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',text)
    p.write_text(text,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    p=OLD/'_shared/navigation_notice.md'
    prefix,history=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+
        ' [1210报告](archive_1201_/research_note_1210.md) · [1211报告](archive_1201_/research_note_1211.md) · '
        '[当前核验](archive_1201_/_shared/'+OUT.name+')。旧目标均未追溯标记完成。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+history,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check'); assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation'); assert navverify['status']=='PASS'
    write(OUT,dict(status='PUBLISHING'))
    docs=[STAGE/'README.md',STAGE/'文件索引.md']
    for n in NUMBERS: docs += [STAGE/f'research_note_{n}.md',*(STAGE/str(n)).glob('*.md')]
    docs += [BASE/name for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks=links(docs)
    for p in docs:
        value=p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in value)
    management={name:sha(ROOT/name) for name in previous['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    write(OUT,dict(schema='thermodynamics_progress_1211_v1',date='2026-10-14',status='PASS',
        accepted_rounds=list(range(1201,1212)),stage_check_groups=22,
        goal_changed=False,application_goal_status='active',stage_complete=False,
        goal_objective_sha256=previous['goal_objective_sha256'],
        acceptance_sha256={**previous['acceptance_sha256'],**new_acceptances},
        science_files_sha256={**previous['science_files_sha256'],**additions},
        management_sha256=management,screen_files_sha256=previous['screen_files_sha256'],
        audit_files_sha256=previous['audit_files_sha256'],
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],stage_requirement_status=req,
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        science_program_outputs=outputs,
        prior_goal_turn_classification='progress: new user-authorized four-laws objective and verified navigation',
        this_goal_turn_result='progress: two complete scoped theorem/counterexample rounds with independent checks; four-laws goal remains active'))
    verify()


if __name__=='__main__': main()
