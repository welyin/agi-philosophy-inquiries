"""941 delivery and evidence checks, distinct from scientific recalculation."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_941_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,941):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'composite_operation_bridge_results.json')
    assert r['round']==941 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert r['inherited_internal_dimension']==60 and r['retained_pair_dimension']==6
    assert r['redshift']['all_unknown_input_transfer_error']<2e-12
    assert r['redshift']['finite_record_probability_gap']>.002
    assert r['redshift']['energy_error']<2e-12
    assert r['fixed_momentum']['every_input_and_ancilla_trace_distance_bound']<.1
    assert r['fixed_momentum']['full_unitary_error']<=r['fixed_momentum']['unitary_bound']
    assert r['common_velocity']['all_unknown_input_transfer_error']<3e-12
    assert r['source']['frechet_vs_independent_five_point_error']<1e-6
    assert r['source']['scalar_chain_rule_nonhermitian_norm']>.003
    assert max(r['exact_regrouping'].values())<2e-12
    assert r['mass_energy_and_quantum_equivalence_coupling_are_physical_inputs']
    for k in ('original_field_material_binding_or_full_worldline_matching_proved','dynamical_gravitational_backreaction_or_full_SM_proved','new_clock_precision_optimization_performed','full_goal_completed'):assert not r[k]
    note=STAGE/'research_note_941.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('整体目标未完成','这不是从认知推出等效原理','没有证明原场物种已制造这份复合组织','不再做复合钟精度','不是从原标准模型场论到点粒子的误差界'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/composite_process_decision.md',HERE/'drafts/common_model_with_material_roles.md',STAGE/'942/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至941）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in newdocs+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=941)==list(range(1,942))
    assert '001—941轮共941份' in nav[0].read_text('utf-8-sig')
    assert '231—941的711份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '941：组织更新、复合钟与质量来源的共同扩展' in s and '942/drafts/STATUS.md' in s
    oldfiles=[STAGE/'940/research_round_940_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'composite_operation_bridge.py',HERE/'composite_operation_bridge_results.json',HERE/'drafts/update_navigation941.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=941,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=941,
       cumulative_numbered_test_groups_from_940=previous['cumulative_numbered_test_groups_from_939']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       existing_nontrivial_internal_process_reused=True,
       record_clock_mass_and_full_source_jointly_checked=True,
       finite_momentum_error_and_noncommuting_source_explicit=True,
       equivalence_coupling_and_material_matching_declared_inputs=True,
       candidate_refinement_stopped_for_common_material_selection=True,
       full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
       frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},
       new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run(a.write)
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
