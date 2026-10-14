"""Verify adopted-role closure, then publish only after the goal tool confirms it."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1201_'
HERE = STAGE / '_shared'
OLD = BASE / 'archive_1161_'
PREV = HERE / 'decision_gate_checks_after1218.json'
READY = HERE / 'four_laws_closure_readiness.json'
OUT = HERE / 'stage_completion_checks_after1218.json'
TOOL = HERE / 'goal_four_laws_completion_receipt.json'
PAPER = BASE / '08_认知热学角色与热力学四律_阶段论文.md'
ADOPTION = HERE / 'thermal_roles_adoption_after1218.md'
AUDIT = HERE / 'four_laws_goal_completion_audit.md'
REVIEW = HERE / 'final_four_laws_paper_review.md'
BUILD = HERE / 'prepare_final_four_laws_paper.py'
spec = importlib.util.spec_from_file_location('frozen_tools', HERE / 'publish_1206.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
sha, read, write, run, links = m.sha, m.read, m.write, m.run, m.links

ROLES = {
    'TR': 'adopted ideal thermal-reference department requirement',
    'TW': 'adopted optional sufficient work implementation subclass',
    'TC': 'adopted resource and theorem-domain registration; not a positivity axiom',
}
STATUS = (
    '1201—1218热力学阶段已完成原目标声明范围的逐项验收，应用目标状态complete。'
    '正式采用TR共同热参考部门要求、TW可选充分工作子类和TC资源登记；'
    '第零律、完整能源第一律、Clausius及循环第二律、Nernst与有限合格操作第三律已有共同定理及有限全历史证据。'
    '正式论文08已完成终审。保留具体模型、操作权限和有效范围；全部机器的统一限时强版本、经验真实性及GR不在本次完成声明中。'
    '阶段仍为18轮、36组科学验收，历史稿完整保留，不启动下一目标。'
)


def validate_docs():
    spec2 = importlib.util.spec_from_file_location('paper_builder', BUILD)
    builder = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(builder)
    assert PAPER.read_text(encoding='utf-8') == builder.build(), 'Paper differs from reviewed build'
    adoption = ADOPTION.read_text(encoding='utf-8')
    assert '采用分层定义，进入论文与逐律验收（推荐）' in adoption
    assert 'call_1d01b5e353a64d76be5fd8c5520ff783' in adoption
    audit = AUDIT.read_text(encoding='utf-8')
    review = REVIEW.read_text(encoding='utf-8')
    assert '逐项目标验收通过' in audit
    assert '正式论文终审通过' in review
    for path in [PAPER, ADOPTION]:
        assert sha(path) in audit and sha(path) in review, path.name
    for path in [PAPER, ADOPTION, AUDIT, REVIEW]:
        value = path.read_text(encoding='utf-8')
        assert not any(x in value for x in (r'\[', r'\]', r'\(', r'\)')), path.name
        assert not re.search(r'^\s*\$\s*$', value, flags=re.M), path.name
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in value), path.name
    return links([PAPER, ADOPTION, AUDIT, REVIEW])


def readiness():
    assert not READY.exists() and not OUT.exists(), 'Readiness already recorded'
    prev = read(PREV)
    old, bindings = m.frozen(prev)
    m.m.validate_map(prev['management_sha256'])
    nlinks = validate_docs()
    prompt = (HERE / 'goal_four_laws_prompt.txt').read_text(encoding='utf-8').strip()
    assert hashlib.sha256(prompt.encode('utf-8')).hexdigest() == prev['goal_objective_sha256']
    replay = run(HERE / 'after1218_common_prediction_check.py')
    assert replay['status'] == 'PASS' and replay['all_checks_passed'] and replay['results_match']
    files = [PAPER, ADOPTION, AUDIT, REVIEW, BUILD, Path(__file__)]
    value = dict(schema='four_laws_closure_readiness_v1',date='2026-10-15',status='READY',
        goal_ready_for_completion=True,goal_objective_sha256=prev['goal_objective_sha256'],
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(prev['science_files_sha256']),
        accepted_rounds=prev['accepted_rounds'],stage_check_groups=36,
        source_authority='Explicit human reply to the recorded role-adoption question, 2026-10-15',
        adopted_thermal_roles=ROLES,new_universal_axioms=0,new_scientific_rounds=0,
        supporting_existing_certificate_replay=replay,local_links_checked=nlinks,
        scope='Original goal in its declared compatible reference, operation, resource and effective domains; not unrestricted machines, laboratory confirmation or GR',
        required_completion_evidence=[AUDIT.relative_to(ROOT).as_posix(),REVIEW.relative_to(ROOT).as_posix()],
        previous_goal_state='blocked; human role decision now explicitly supplied',
        application_goal_marked_complete=False)
    write(READY, value)
    print(json.dumps({k:value[k] for k in ['status','goal_ready_for_completion',
        'old_GR_science_files_verified','old_GR_acceptance_bindings_verified',
        'current_stage_science_files_verified','stage_check_groups','local_links_checked']}))


def verify():
    d = read(OUT)
    old, bindings = m.frozen(d)
    m.m.validate_map(d['management_sha256'])
    assert sha(PREV) == d['prior_progress_sha256']
    assert sha(READY) == d['closure_readiness_sha256']
    assert sha(TOOL) == d['goal_completion_tool_receipt_sha256']
    assert read(TOOL)['goal']['status'] == 'complete'
    ledger = read(HERE / 'dependency_ledger.json')
    assert ledger['status'] == 'complete'
    assert ledger['latest_accepted_round'] == 1218 and ledger['stage_check_groups'] == 36
    assert ledger['latest_management_receipt'] == OUT.name
    assert ledger['adopted_thermal_roles'] == ROLES
    assert ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert run(OLD / '_shared/update_navigation.py', '--check')['changed_files'] == 0
    assert run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')['status'] == 'PASS'
    validate_docs()
    print(json.dumps(dict(status='PASS',goal_status='complete',stage_complete=True,
        stage_scope='declared compatible thermal domains',latest_scientific_round=1218,
        stage_rounds=18,stage_check_groups=36,new_scientific_rounds=0,
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        adopted_thermal_roles=list(ROLES),new_universal_axioms=0,old_GR_goal_completed=False)))


def publish():
    assert not OUT.exists(), 'Already published'
    ready = read(READY)
    assert ready['goal_ready_for_completion'] and ready['status'] == 'READY'
    m.m.validate_map(ready['files_sha256'])
    prev = read(PREV)
    old, bindings = m.frozen(prev)
    m.m.validate_map(prev['management_sha256'])
    assert sha(PREV) == ready['prior_progress_sha256']
    goal = read(TOOL)['goal']
    assert goal['status'] == 'complete'
    assert hashlib.sha256(goal['objective'].encode('utf-8')).hexdigest() == prev['goal_objective_sha256']
    ledger = read(HERE / 'dependency_ledger.json')
    req = ledger['stage_requirement_status']
    req.update(
        zeroth_law='Completed in adopted TR non-scalar, same-object/unit domain: equilibrium equivalence, ordered common beta and conditional finite thermometry; actual sub-menu scope retained',
        first_law='Completed under adopted D-P: same autonomous generator, complete energy conservation, operational energy/work/heat roles and preparation/control/record/backreaction evidence',
        second_law='Completed: exact full Clausius identity with defects/correlations and its cyclic Kelvin implication; feedback/reset retained, no universal material-unital or reverse engine-permission claim',
        third_law='Completed in declared finite resource/operation domains: fixed and uniformly gapped Nernst, general qualified-event bound and explicit safe-work class, fixed-budget bound; stronger all-machine fixed-time result retains independent hypotheses',
        common_process='Nonempty compatible class and certified native full preparation/contact/record history; all finite-source non-Gibbs and pure-stock distinctions preserved',
        geometry_handoff='Required identity audit and remaining same-object clock/position/signal/region, Unruh and independent-area/GR bridges completed as handoff, not derived here',
        stage_synthesis='Formal adopted-role paper and independent objective audit completed',
        final_stage_paper='published 08_认知热学角色与热力学四律_阶段论文.md',
        per_law_evidence_audit='passed against original goal and approved layered definitions',
        stage_deliverables_all_completed=True)
    ledger.update(status='complete',current_work=STATUS,latest_management_receipt=OUT.name,
        latest_goal_status_receipt=TOOL.name,adopted_thermal_roles=ROLES,
        new_adopted_department_requirements=['TR'],new_universal_axioms=0,
        adoption_document=ADOPTION.name,final_paper_document='../../'+PAPER.name,
        current_scope_audit=AUDIT.name,current_direction_document='../../'+PAPER.name,
        next_direction='Same-object thermodynamics-to-geometry handoff; no new goal started')
    write(HERE / 'dependency_ledger.json', ledger)
    oldledger = read(OLD / '_shared/dependency_ledger.json')
    oldledger['latest_management_receipt'] = '../../archive_1201_/_shared/' + OUT.name
    write(OLD / '_shared/dependency_ledger.json', oldledger)
    for name in ['README.md', '文件索引.md']:
        p = STAGE / name
        v = p.read_text(encoding='utf-8')
        if name == 'README.md':
            parts = v.split('\n\n')
            parts[0] = '# 1201—1218：认知操作热力学四律（阶段完成）'
            parts[1] = STATUS
            v = '\n\n'.join(parts)
        v = v.replace('_shared/' + PREV.name, '_shared/' + OUT.name)
        # Preserve earlier decision history while marking its supersession explicitly.
        v = v.replace('## 待决定的理论定位', '## 采用前的理论定位决定记录（已由下节采用与验收接续）')
        v += '\n## 正式采用、论文与结项\n\n'
        v += '- [08阶段论文](../' + PAPER.name + ')\n'
        v += '- [分层热学角色采用](_shared/' + ADOPTION.name + ') · [逐项目标验收](_shared/' + AUDIT.name + ') · [正式论文终审](_shared/' + REVIEW.name + ')\n'
        v += '- [应用完成回执](_shared/' + TOOL.name + ') · [当前结项核验](_shared/' + OUT.name + ')\n\n'
        v += '旧目标、候选报告与预审稿中的状态保留原时点；本节为1201—1218的最终状态。未新增科学轮，未改变1161—1200的GR未完成状态，未启动下一目标。\n'
        p.write_text(v, encoding='utf-8')
    notice = OLD / '_shared/navigation_notice.md'
    _, history = notice.read_text(encoding='utf-8').split('以下保留前阶段的历史状态', 1)
    prefix = '## 已完成阶段（1201—1218）：认知操作热力学四律\n\n'
    prefix += '2026-10-15，用户明确采用分层热学定义后，正式论文与逐项目标验收完成，应用目标状态 **complete**。\n\n'
    prefix += '[08阶段论文](' + PAPER.name + ') · [阶段目录](archive_1201_/README.md) · '
    prefix += '[正式采用](archive_1201_/_shared/' + ADOPTION.name + ') · [逐项验收](archive_1201_/_shared/' + AUDIT.name + ') · '
    prefix += '[完成回执](archive_1201_/_shared/' + TOOL.name + ') · [阶段索引](archive_1201_/文件索引.md)。\n\n'
    prefix += '**当前状态：** ' + STATUS + ' [结项核验](archive_1201_/_shared/' + OUT.name + ')。\n\n'
    prefix += '后继方向是同一实际过程的热力学—几何连接：实际钟、三维位置、信号、区域能流与参考身份，再到Unruh及独立面积响应。保留旧空间接口，不把普通四律或内部控制坐标直接当作Einstein动力学。\n\n'
    notice.write_text(prefix + '以下保留前阶段的历史状态' + history, encoding='utf-8')
    navapply = run(OLD / '_shared/update_navigation.py', '--apply')
    navcheck = run(OLD / '_shared/update_navigation.py', '--check')
    navverify = run(ROOT / 'scripts/run_research_gr_active.py', '--verify-navigation')
    assert navcheck['changed_files'] == 0 and navverify['status'] == 'PASS'
    write(OUT, dict(status='PUBLISHING'))
    paths = [PAPER, ADOPTION, AUDIT, REVIEW, STAGE / 'README.md', STAGE / '文件索引.md']
    paths += [BASE / n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks = links(paths)
    audits = dict(prev['audit_files_sha256'])
    for p in [PAPER, ADOPTION, AUDIT, REVIEW, BUILD, READY, TOOL]:
        audits[p.relative_to(ROOT).as_posix()] = sha(p)
    management = {name:sha(ROOT / name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    result = dict(schema='four_laws_stage_completion_v1',date='2026-10-15',status='PASS',
        application_goal_status='complete',stage_complete=True,goal_changed=False,
        goal_objective_sha256=prev['goal_objective_sha256'],accepted_rounds=prev['accepted_rounds'],
        stage_check_groups=36,new_scientific_rounds=0,new_empirical_calibration=0,
        new_universal_axioms=0,new_adopted_department_requirements=['TR'],adopted_thermal_roles=ROLES,
        acceptance_sha256=prev['acceptance_sha256'],science_files_sha256=prev['science_files_sha256'],
        screen_files_sha256=prev['screen_files_sha256'],audit_files_sha256=audits,
        management_sha256=management,prior_progress_path=PREV.relative_to(ROOT).as_posix(),
        prior_progress_sha256=sha(PREV),closure_readiness_sha256=sha(READY),
        goal_completion_tool_receipt_sha256=sha(TOOL),stage_requirement_status=req,
        final_paper=PAPER.relative_to(ROOT).as_posix(),
        requirement_audit=AUDIT.relative_to(ROOT).as_posix(),
        independent_paper_review=REVIEW.relative_to(ROOT).as_posix(),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        old_GR_goal_completed=False,GR_hard_deliverables_completed=0,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],
        supporting_existing_certificate_replay=ready['supporting_existing_certificate_replay'],
        prior_goal_turn_classification='no progress: genuine role decision pending; goal was blocked',
        this_goal_turn_result='human adopted layered definitions; final paper, explicit safe-work third-law corollary and original-goal audit verified; goal tool confirmed completion in declared scope',
        next_goal_created=False)
    write(OUT, result)
    verify()


if __name__ == '__main__':
    if '--verify' in sys.argv:
        verify()
    elif '--publish' in sys.argv:
        publish()
    else:
        readiness()
