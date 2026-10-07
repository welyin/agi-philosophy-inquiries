"""927 delivery: joint task/resource process and scope of future encoding invariance."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_927_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,927):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'joint_encoding_resource_selection_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['all_scientific_checks_passed'] and result['strong_future_contract_mixing_dimension']==1
    uniform,distinct=result['joint_candidate_family']
    for row in (uniform,distinct):
        assert row['physical_dimension']==40 and row['code_dimension']==32
        assert row['positive_total_H_minimum']>1.99
        for key in ('code_intertwining_error','momentum_code_intertwining_error','all_unknown_input_transfer_error','bare_energy_conservation_error','resource_current_identity_error','endpoint_current_identity_error','role_exchange_error_on_code','rotation_error_on_code','instant_encoding_commutator_norm'):assert row[key]<5e-13,key
        assert abs(row['encoding_bare_energy_commutators'][0]-4)<1e-13
        assert abs(row['mean_energy_changes']['hr']-.35)<1e-12
        assert abs(row['mean_energy_changes']['he']+.35)<1e-12
        assert abs(row['mean_energy_changes']['H'])<1e-12
    assert uniform['first_future_derivative_encoding_defect']==0 and uniform['finite_future_encoding_commutator_norm']<1e-12
    assert uniform['finite_probability_difference']<1e-12
    assert distinct['first_future_derivative_encoding_defect']==8
    assert abs(distinct['finite_probability_difference']-2**-.5)<1e-12
    assert distinct['finite_future_encoding_commutator_norm']>1
    for key in ('one_process_carries_records_direction_and_resource','unknown_record_kept_in_joint_reference_code','reference_preparation_common_to_unknown_inputs','resource_shift_not_claimed_new','uniform_coupling_not_assumed_as_cognitive_requirement','strong_future_invariance_is_extra_input','active_encodings_not_claimed_costless'):assert result[key]
    for key in ('autonomous_reencoding_controller_constructed','spacetime_dimension_generated','gauge_group_or_GR_generated','all_physics_inputs_selected','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_927.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==18 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,10)]
    assert '本轮到此停止这类分类' in prose and '925曾统一设g' in prose
    newdocs=[note,HERE/'drafts/STATUS.md',STAGE/'928/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
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
    assert sorted(n for n in nums if n<=927)==list(range(1,928))
    assert '001—927轮共927份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—927的697份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '927：同一过程可相容' in s and '928/drafts/STATUS.md' in s
    oldfiles=[STAGE/'926/research_round_926_checks.json',STAGE/'_shared/notes/candidate_priority_audit_20261006.md']
    newfiles=newdocs+[Path(__file__),HERE/'joint_encoding_resource_selection.py',HERE/'joint_encoding_resource_selection_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'926/research_round_926_checks.json')
    out=dict(round=927,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=927,
      cumulative_numbered_test_groups_from_926=previous['cumulative_numbered_test_groups_from_925']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      single_process_joint_selection_checked=True,full_unknown_input_isometry_checked=True,
      operator_resource_current_checked=True,finite_future_probability_counterexample=True,
      instantaneous_and_future_contracts_distinguished=True,uniform_g_input_not_smuggled_into_conclusion=True,
      autonomous_control_not_claimed=True,prior_record_and_reference_theorems_reused=True,
      classification_branch_stopped=True,next_question_is_upstream_source_selection=True,
      physical_gauge_group_or_GR_generated=False,full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
      frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
