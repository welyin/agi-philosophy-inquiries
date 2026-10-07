"""934 delivery: operation group reconstruction and internal reference boundary."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_934_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,934):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'comparison_group_reference_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==934 and result['all_scientific_checks_passed']
    assert result['group_size']==48 and result['closure_error']<1e-12
    rows=result['moment_rows']
    assert [x['finite_commutant_dimension'] for x in rows]==[1,2,5,15]
    assert [x['continuous_commutant_dimension'] for x in rows]==[1,2,5,14]
    assert max(x['full_channel_frobenius_difference'] for x in rows[:3])<1e-12
    assert abs(rows[3]['full_channel_frobenius_difference']-1)<1e-12
    for row in result['internal_reference_rows']:
        assert abs(row['probability']-row['full_pointer_probability'])<1e-12
        for k in ('common_frame_commutator_error','trace_error','source_marginal_error','reference_marginal_error'):assert row[k]<1e-12
        assert row['smallest_eigenvalue']>-1e-12
    assert abs(result['internal_probability_gap']-1/240)<1e-12
    assert result['continuous_joint_factorization_error']<1e-12
    for k in ('physical_comparison_permissions_and_preparations_remain_inputs','mature_three_design_result_reused','candidate_detail_optimization_stopped'):assert result[k]
    for k in ('all_six_protocols_realized','space_or_gauge_dynamics_derived','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_934.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('成熟','完整能源账','实际比较权限','未结项','三份','1/240'):assert phrase in prose,phrase
    working=HERE/'drafts/operation_hypotheses_screening.md'
    assert '不是七条已采用公理' in working.read_text('utf-8')
    newdocs=[note,working,HERE/'drafts/role_algebra_reuse_review.md',STAGE/'935/drafts/STATUS.md',STAGE/'935/drafts/unified_operation_hypotheses_v0_1.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至934）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    assert '当前认知假说补充：共同指认的六条协议' in nav[1].read_text('utf-8-sig')
    assert '已采用六条共同认知协议' in nav[2].read_text('utf-8-sig')
    links=0
    for doc in newdocs+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=934)==list(range(1,935))
    assert '001—934轮共934份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—934的704份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '934：组合比较的群重建接口' in s and '935/drafts/STATUS.md' in s
    oldfiles=[STAGE/'933/research_round_933_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'comparison_group_reference.py',HERE/'comparison_group_reference_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'933/research_round_933_checks.json')
    out=dict(round=934,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=934,
      cumulative_numbered_test_groups_from_933=previous['cumulative_numbered_test_groups_from_932']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      old_relation_carrier_and_actual_reference_results_reused=True,
      mature_group_reconstruction_and_unitary_design_results_attributed=True,
      internal_joint_readout_not_replaced_by_external_absolute_axis=True,
      finite_group_and_continuous_group_not_inferred_from_same_low_order_data=True,
      no_unique_H_or_group_requirement_added_to_existence_stage=True,
      seven_proposals_registered_not_all_adopted_as_axioms=True,
      fourth_order_instrument_optimization_stopped=True,
      full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
      frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
