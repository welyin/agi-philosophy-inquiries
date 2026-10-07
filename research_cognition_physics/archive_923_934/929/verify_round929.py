"""929 delivery: joint finite contracts and full-spectral update comparison."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_929_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,929):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'joint_protocol_selection_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['round']==929 and result['all_scientific_checks_passed']
    assert result['unknown_input_and_fault_sectors_checked']==32 and result['repeats']==2
    for key in ('complete_output_isometry_max_error','unknown_input_isometry_max_error','actual_record_and_syndrome_relabeling_error','clock_spectrum_error','clock_spectral_weight_error','complete_gate_weighted_intertwining_error','autonomous_finish_error'):assert result[key]<1e-12,key
    rows=result['compared_candidates']
    assert abs(rows[0]['joint_event_probabilities_for_zero'][0]-.64)<1e-13
    assert abs(rows[1]['joint_event_probabilities_for_zero'][0]-.13)<1e-13
    assert abs(result['event_00_probability_difference']-.51)<1e-13
    for row in rows:
        assert abs(sum(row['joint_event_probabilities_for_zero'])-1)<1e-13
        assert row['nonconstant_00_effect_gap']>.38
    assert (result['logical_steps'],result['clock_dimension'],result['work_qubits'],result['internal_fault_qubits'],result['blank_work_qubits'])==(30,31,13,4,12)
    assert result['full_hilbert_dimension']==4063232
    assert (result['full_H_minimum'],result['full_H_maximum'],result['initial_full_energy'])==(0,30,15)
    assert abs(result['initial_full_energy_variance']-7.5)<1e-12
    assert 0<result['terminal_window_failure']<=result['terminal_window_failure_upper_bound']==.003
    for key in ('full_spectral_input_measure_same_for_all_unknown_inputs_and_candidates','same_six_finite_contracts_jointly_realized','all_blanks_fault_registers_syndromes_and_clock_internal','initial_low_entropy_preparation_and_engineered_H_are_inputs','history_compilation_is_inherited_not_new_physics','full_quantum_postevent_state_and_passive_reference_retained'):assert result[key]
    for key in ('arbitrary_record_phase_attacks_or_controller_errors_certified','eternal_repetition_or_fixed_finite_device_all_scale_claimed','spatial_or_GR_interfaces_realized','all_cognitive_principles_inadequate_for_all_physics_proved','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_929.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('停止这条协议编译分支','完整初始能谱分布相同','任意局部相位攻击','低熵空白','至少.504','没有选出θ','本轮没有证明'):
        assert phrase in prose,phrase
    working=HERE/'drafts/joint_selection_working.md'
    assert '至多一个主体的记录位翻转' in working.read_text('utf-8')
    assert 'theta=0与theta=pi/2' in working.read_text('utf-8')
    newdocs=[note,working,STAGE/'930/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至929）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=929)==list(range(1,930))
    assert '001—929轮共929份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—929的699份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '929：六条有限协议共同实现' in s and '930/drafts/STATUS.md' in s
    oldfiles=[STAGE/'928/research_round_928_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'joint_protocol_selection.py',HERE/'joint_protocol_selection_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'928/research_round_928_checks.json')
    out=dict(round=929,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=929,
      cumulative_numbered_test_groups_from_928=previous['cumulative_numbered_test_groups_from_927']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      same_object_six_declared_finite_contracts_checked=True,
      full_unknown_input_and_passive_reference_isometry_checked=True,
      all_declared_fault_choices_checked=True,
      full_spectral_resource_measure_compared_not_only_energy_mean=True,
      finite_readout_error_cannot_hide_probability_difference=True,
      inherited_clock_and_record_tools_not_claimed_new=True,
      physical_source_and_geometry_not_claimed_generated=True,
      low_entropy_initial_state_and_engineered_H_kept_as_inputs=True,
      protocol_processor_branch_stopped=True,
      universal_program_failure_not_claimed=True,
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
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
