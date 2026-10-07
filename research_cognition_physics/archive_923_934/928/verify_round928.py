"""928 delivery: conditional source naturality, boundary witness, and scoped branch stop."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_928_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,928):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'source_grouping_selection_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
    assert result['round']==928 and result['all_scientific_checks_passed']
    w=result['weighted_family_exact_calibration'];u=result['unweighted_family_exact_calibration']
    assert (w['objects'],w['variables'],w['morphisms'],w['exact_rank'],w['affine_dimension'])==(8,56,124,55,1)
    assert w['candidate_null_vector_residual']==0 and w['rank_with_one_isolated_pair']==56
    assert (u['objects'],u['variables'],u['morphisms'],u['exact_rank'],u['affine_dimension'])==(4,30,92,30,0)
    assert result['weighted_naturality_checks']==75 and result['weighted_naturality_max_error']<1e-14
    assert abs(result['forgotten_weight_max_entry_defect']-1/12)<1e-14
    assert result['finite_grid_b_above_one_positive_minimum']>0 and result['smaller_weight_b_above_one_minimum']<0
    local=result['local_counterexample']
    assert local['coarse_total_difference']==0 and local['coarse_source_difference']==[.25,-.25]
    assert local['boundary_interface_error']==0 and local['restricted_pair_grouping_error']==0
    probe=result['autonomous_quantum_source_probe']
    assert probe['dimension']==32 and probe['minimum_total_energy']==0
    assert probe['coarse_initial_energies']==[[1.,0.],[1.,0.]]
    assert abs(probe['full_initial_energies'][0]-1.5)<1e-13 and abs(probe['full_initial_energies'][1]-1.625)<1e-13
    assert abs(probe['final_plus_probabilities'][0]-1)<1e-13 and abs(probe['final_plus_probabilities'][1])<1e-13
    assert max(abs(x) for x in probe['total_energy_changes'])<1e-13 and abs(probe['coherent_probe_purity']-.5)<1e-13
    for key in ('six_protocol_gap_map_covers_C01_to_C27','all_regrouping_summary_is_extra_assumption','primitive_operator_assignment_and_weights_remain_inputs','source_allocation_not_quantum_state_channel'):assert result[key]
    for key in ('all_six_protocols_jointly_realized','tamper_resistance_and_repeated_reset_constructed','locality_derived_from_protocols','unique_geometric_stress_derived','spacetime_or_GR_generated','physical_minimum_scale_or_exact_continuum_assumed','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_928.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==24 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,13)]
    for phrase in ('本轮到此停止这类分类','完整闭合初始能量不同','协议2','协议6','总容量V','b=1.2','不是声称能同时精确测量全部h_i'):
        assert phrase in prose,phrase
    newdocs=[note,STAGE/'929/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至928）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=928)==list(range(1,929))
    assert '001—928轮共928份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—928的698份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '928：强分组合同约束来源' in s and '929/drafts/STATUS.md' in s
    oldfiles=[STAGE/'927/research_round_927_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'source_grouping_selection.py',HERE/'source_grouping_selection_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'927/research_round_927_checks.json')
    out=dict(round=928,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=928,
      cumulative_numbered_test_groups_from_927=previous['cumulative_numbered_test_groups_from_926']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      weighted_naturality_general_proof_and_finite_rank_calibration=True,
      finite_weight_domain_distinguished_from_universal_domain=True,
      boundary_source_finite_quantum_prediction_witness=True,
      full_initial_energy_difference_explicit=True,
      source_allocation_not_misidentified_as_quantum_state_channel=True,
      six_protocol_gap_map_C01_through_C27_verified=True,
      user_adopted_six_protocol_sections_preserved=True,
      all_regrouping_contract_not_inferred_from_protocol_four=True,
      source_classification_branch_stopped=True,
      original_local_sources_weights_and_geometry_not_derived=True,
      all_six_protocols_jointly_realized=False,
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
