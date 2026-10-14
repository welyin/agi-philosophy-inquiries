"""Register the reviewed paper draft without adopting roles or closing the stage."""
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
PREV = HERE / 'role_synthesis_checks_after1218.json'
OUT = HERE / 'paper_draft_checks_after1218.json'
spec = importlib.util.spec_from_file_location('frozen_tools', HERE / 'publish_1206.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links
DOCS = [
    'four_laws_stage_paper_draft_after1218.md',
    'paper_audit_zeroth_first_after1218.md',
    'paper_audit_second_third_after1218.md',
    'paper_literature_scope_after1218.md',
    'paper_independent_review_after1218.md',
]
STATUS = (
    '最新科学验收1218；本阶段18轮、36组。四律完整论文预审稿、逐律证据审计和文献范围核对已完成独立审阅；'
    '参考分类、完整能源、Clausius到循环限制、静态与动态第三律分别给出证明及共同范围。'
    'TR热参考角色仍待采用，TW仅为可选充分工作实现，TC仅登记资源范围；'
    '全部限时机器的统一冷却界未获证明。原有限过程与成熟工具的权限边界保留。'
    '本次为论文准备，不新增科学轮次或公理，目标继续有效，阶段未结项。'
)


def requirement_coverage():
    """Evidence routing, not a boolean proof of goal completion."""
    return {
        'zeroth_equilibrium_and_temperature': {
            'status': 'proved_in_declared_role_domain; role_adoption_pending',
            'evidence': ['1209/proof.md sections 2-5', '1210/proof.md sections 2-7',
                         '_shared/paper_audit_zeroth_first_after1218.md section 3'],
            'boundary': 'Mathematical full-menu classification does not certify a restricted actual menu or exact Gibbs membership from finite data.'},
        'first_clock_energy_work_heat': {
            'status': 'scoped_theorem_and_model_prediction',
            'evidence': ['_shared/four_laws_dependency_after1211.md section 2',
                         '1203/proof.md', '1205/proof.md', '1212/proof.md sections 2-3',
                         '_shared/paper_audit_zeroth_first_after1218.md section 4'],
            'boundary': 'D-P identity inherited; native contact and preparation-expanded parents distinct; no delivered mechanical-work claim.'},
        'second_clausius_and_cyclic_relation': {
            'status': 'scoped_forward_relation_proved',
            'evidence': ['1212/proof.md sections 3,8',
                         '_shared/four_laws_stage_paper_draft_after1218.md section 5',
                         '_shared/paper_audit_second_third_after1218.md section 2'],
            'boundary': 'All nonbath/nonwork participants including target return; correlations and non-Gibbs defects retained. Reverse equivalence or reversible attainment not claimed.'},
        'third_static': {
            'status': 'fixed_spectrum_and_uniform_gapped_family_proved',
            'evidence': ['1211/proof.md section 7',
                         '_shared/paper_audit_second_third_after1218.md section 3.1'],
            'boundary': 'Positive-temperature limit, fixed degeneracy and stated gap/dimension conditions; no universal zero residual entropy.'},
        'third_dynamic': {
            'status': 'qualified_finite_protocol_and_fixed_resource_bounds_proved; unrestricted_fixed_time_claim_unproved',
            'evidence': ['1217/proof.md sections 2-7', '1215/proof.md', '1218/proof.md',
                         '_shared/paper_audit_second_third_after1218.md sections 3.2-3.4'],
            'boundary': 'q positivity from independent channel qualification, not a new axiom. Pure stocks and feedback retained. Time/resource family control not inferred from SR.'},
        'cognitive_sources_selection_deletions': {
            'status': 'dependency_and_candidate_audit_ready; no_new_axiom_adopted',
            'evidence': ['_shared/cognitive_thermal_roles_proposal_after1218.md',
                         '_shared/four_laws_stage_paper_draft_after1218.md sections 2,7',
                         '1206/proof.md', '1208/proof.md', '1211/proof.md', '1218/proof.md'],
            'boundary': 'TR is a proposed department role, TW an optional sufficient implementation, TC applicability registration; not universal state restrictions.'},
        'compatible_common_scope': {
            'status': 'nonempty_mathematical_class_and_separate_finite_history_evidence',
            'evidence': ['_shared/paper_audit_second_third_after1218.md section 4',
                         '_shared/after1218_common_realization_review.md'],
            'boundary': 'Nonempty class does not supply every geometric task or actual full-menu permissions.'},
        'finite_preparation_control_records_backreaction': {
            'status': 'fixed_native_whole_history_certificate_verified',
            'evidence': ['1205/proof.md', '_shared/after1218_common_prediction_results.json',
                         '_shared/after1218_joint_theorem_independent_review.md section 5'],
            'boundary': 'All trials kept, original non-Gibbs defects and pure preparation stocks retained; not all-future cooling or cyclic-engine proof.'},
        'clock_position_signal_region_geometry_handoff': {
            'status': 'identities_audited_and_remaining_bridges_listed',
            'evidence': ['_shared/current_process_geometry_interfaces.md',
                         '_shared/thermal_operational_bridge_audit_after1209.md',
                         '_shared/four_laws_stage_paper_draft_after1218.md section 9'],
            'boundary': 'No area definition from entropy, no automatic Unruh/Einstein result; inherited spatial theorems reused with their original scope.'},
        'paper_and_per_law_review': {
            'status': 'complete_pre_adoption_draft_reviewed; final_closure_not_granted',
            'evidence': ['_shared/' + n for n in DOCS],
            'boundary': 'This records a reviewed draft and objective audit, not role adoption or completed derivation under the current axioms.'},
        'reproducibility_and_historical_integrity': {
            'status': 'frozen_science_hashes_and_existing_prediction_replay_verified',
            'evidence': ['_shared/paper_draft_checks_after1218.json'],
            'boundary': 'Metadata checks validate provenance and routing; mathematical review bears theorem validity.'},
    }


def verify():
    d = read(OUT)
    old, bindings = m.frozen(d)
    m.m.validate_map(d['management_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    ledger = read(HERE / 'dependency_ledger.json')
    assert ledger['latest_management_receipt'] == OUT.name
    assert ledger['latest_accepted_round'] == 1218
    assert ledger['stage_check_groups'] == 36
    assert ledger['goal_objective_sha256'] == d['goal_objective_sha256']
    assert ledger['new_adopted_axioms'] == []
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert d['objective_requirement_coverage'] == requirement_coverage()
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    print(json.dumps(dict(status='PASS',latest_scientific_round=1218,
        stage_rounds=18,stage_check_groups=36,paper='reviewed_pre_adoption_draft',
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        new_scientific_rounds=0,new_adopted_axioms=0,goal_changed=False,stage_complete=False)))


def main():
    if '--verify' in sys.argv:
        return verify()
    assert not OUT.exists(), 'Already published'
    prev = read(PREV)
    old, bindings = m.frozen(prev)
    m.m.validate_map(prev['management_sha256'])
    review = (HERE / DOCS[-1]).read_text(encoding='utf-8')
    assert '论文预审通过' in review
    for name in DOCS[:-1]:
        assert sha(HERE / name) in review, name
    for name in DOCS:
        value = (HERE / name).read_text(encoding='utf-8')
        assert not any(x in value for x in (r'\[',r'\]',r'\(',r'\)')), name
        assert not re.search(r'^\s*\$\s*$', value, flags=re.M), name
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in value), name
    links([HERE / name for name in DOCS])
    replay = run(HERE / 'after1218_common_prediction_check.py')
    assert replay['status'] == 'PASS' and replay['all_checks_passed'] and replay['results_match']

    ledger = read(HERE / 'dependency_ledger.json')
    req = ledger['stage_requirement_status']
    req.update(
        stage_synthesis='Full pre-adoption paper, per-law objective audits and primary-literature scope review completed; independently reviewed draft. TR adoption and final scope/closure remain separate.',
        final_stage_paper='reviewed_pre_adoption_draft_not_finalized',
        per_law_evidence_audit='completed_for_declared_conditions; not a claim of unconditional cognitive entailment',
        third_law='Fixed-spectrum/uniform-gapped Nernst and qualified finite-protocol/fixed-budget event bounds proved; all-machine fixed-time bound not claimed. Resource and q qualifications must have independent sources.',
        stage_deliverables_all_completed=False)
    ledger.update(current_work=STATUS,latest_management_receipt=OUT.name,
        paper_draft_document=DOCS[0],current_scope_audit=DOCS[-1],
        current_direction_document='cognitive_thermal_roles_proposal_after1218.md')
    write(HERE / 'dependency_ledger.json', ledger)
    oldledger = read(OLD / '_shared/dependency_ledger.json')
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(OLD / '_shared/dependency_ledger.json', oldledger)

    for name in ['README.md','文件索引.md']:
        p = STAGE / name
        v = p.read_text(encoding='utf-8')
        if name == 'README.md':
            parts = v.split('\n\n')
            parts[1] = STATUS
            v = '\n\n'.join(parts)
        v = v.replace('_shared/' + PREV.name, '_shared/' + OUT.name)
        v += '\n## 四律论文预审与逐项验收准备\n\n'
        v += '- [完整阶段论文预审稿](_shared/' + DOCS[0] + ')\n'
        v += '- [第零／第一律证据审计](_shared/' + DOCS[1] + ') · [第二／第三律及共同过程审计](_shared/' + DOCS[2] + ')\n'
        v += '- [成熟文献与归属边界](_shared/' + DOCS[3] + ') · [完整稿独立审阅](_shared/' + DOCS[4] + ')\n'
        v += '\n已完成采用前论文准备；TR采用、实际资格与最终结项分开处理。科学编号仍止于1218，原冻结稿保留。\n'
        p.write_text(v, encoding='utf-8')
    notice = OLD / '_shared/navigation_notice.md'
    prefix, history = notice.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    prefix = re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)',
        '**当前状态：** ' + STATUS +
        ' [论文预审稿](archive_1201_/_shared/' + DOCS[0] + ') · '
        '[逐项独审](archive_1201_/_shared/' + DOCS[-1] + ') · '
        '[当前核验](archive_1201_/_shared/' + OUT.name + ')。旧目标均未追溯标记完成。', prefix, flags=re.S)
    notice.write_text(prefix + '以下保留前阶段的历史状态' + history, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    assert navcheck['changed_files'] == 0
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navverify['status'] == 'PASS'
    write(OUT,dict(status='PUBLISHING'))
    paths = [HERE / n for n in DOCS] + [STAGE / 'README.md', STAGE / '文件索引.md']
    paths += [BASE / n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks = links(paths)
    management = {name:sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    audits = {**prev['audit_files_sha256'],**{
        (HERE / name).relative_to(ROOT).as_posix():sha(HERE / name) for name in DOCS}}
    write(OUT,dict(schema='four_laws_paper_pre_review_after1218_v1',date='2026-10-15',status='PASS',
        accepted_rounds=prev['accepted_rounds'],stage_check_groups=36,goal_changed=False,
        application_goal_status='active',stage_complete=False,new_scientific_rounds=0,
        goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256=prev['acceptance_sha256'],science_files_sha256=prev['science_files_sha256'],
        management_sha256=management,audit_files_sha256=audits,screen_files_sha256=prev['screen_files_sha256'],
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],stage_requirement_status=req,
        objective_requirement_coverage=requirement_coverage(),
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        supporting_existing_certificate_replay=replay,
        prior_goal_turn_classification='no_progress: last turn verified the already active goal and restated scope; no new scientific or documentary state',
        no_progress_revalidation='Prior receipt verified against current files; safe remaining action was full paper drafting and independent per-requirement audit. Adoption not inferred from automatic continuation.',
        this_goal_turn_result='progress: full paper draft with proofs and shared finite evidence, per-law objective audit and checked primary-literature boundaries, all independently reviewed; no goal closure claim'))
    verify()


if __name__ == '__main__':
    main()
