"""937 delivery: classical shared dust reference and original constraints."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_937_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,937):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'dust_common_model_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==937 and result['all_scientific_checks_passed']
    exact=result['exact_reference_and_source_check']
    assert exact['exact_mixed_mass_four_leg']=='135/4928'
    assert exact['exact_total_constraints_zero'] and exact['exact_unit_timelike']
    rows=result['original_full_field_backgrounds']
    assert [(r['N'],r['trace_shift']) for r in rows]==[(16,.01),(16,.03),(24,.01),(24,.03)]
    for r in rows:
        assert r['dust_density']>0 and r['minimum_reference_momentum']>0
        for k in ('new_total_density_residual','new_Einstein_constraint_residual','unchanged_momentum_constraint_residual','unchanged_electroweak_Gauss_residual','color_Gauss_residual'):assert r[k]<1e-10
        assert r['exact_shift_identity_residual']<1e-13
        assert r['canonical_trace_shift_error']<1e-13
        assert r['space_and_matter_canonical_data_preserved']
        assert r['at_rest_mixed_mass_four_leg_coefficient']==0
        assert r['color_electric_coefficient']>0 and r['color_magnetic_coefficient']>0
    moving=result['moving_local_source_check']
    assert abs(moving['explicit_metric_variation'])>1e-4
    for k in ('variation_checks','mixed_mass_checks'):
        assert 3.7<moving[k][0]['error']/moving[k][1]['error']<4.3
    for k in ('new_positive_reference_matter_is_explicit_physical_input','classical_two_derivative_and_zero_fermion_background_scope','full_parent_action_retains_fermions_but_not_fully_quantized','original873_clock_is_different_physical_reference','dust_low_density_limit_not_claimed_uniform','classical_constraints_not_quantum_Ward_certificate','no_new_spatial_dimension_claim'):assert result[k]
    assert not result['full_goal_completed']
    note=STAGE/'research_note_937.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==20 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,11)]
    for phrase in ('新增了实际参考物质','不是另一个已解的全局漂移PDE背景','整体未结项','不是量子Ward核验','停止尘埃图册和精确时钟优化'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/dust_candidate_decision.md',STAGE/'938/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至937）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    assert '当前认知假说补充：共同指认的六条协议' in nav[1].read_text('utf-8-sig')
    assert '已采用六条共同认知协议' in nav[2].read_text('utf-8-sig')
    for q in nav[1:3]:
        s=q.read_text('utf-8-sig')
        assert '先判定方向价值，再投入候选细节' in s
        assert '不预设微观连续或物理最小尺度' in s
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
    assert sorted(n for n in nums if n<=937)==list(range(1,938))
    assert '001—937轮共937份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—937的707份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '937：正能源内部参考与原全场约束共同接入' in s and '938/drafts/STATUS.md' in s
    oldfiles=[STAGE/'936/research_round_936_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'dust_common_model.py',HERE/'dust_common_model_results.json',HERE/'drafts/update_navigation937.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'936/research_round_936_checks.json')
    out=dict(round=937,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=937,
      cumulative_numbered_test_groups_from_936=previous['cumulative_numbered_test_groups_from_935']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      actual_original_background_reused=True,new_reference_input_declared=True,
      known_reference_theory_attributed=True,reference_energy_and_root_sign_checked=True,
      source_metric_variation_retained=True,old873_distinct_clock_preserved=True,
      zero_fermion_background_not_claimed_to_quantize_all_fields=True,
      local_moving_algebra_not_claimed_global_solution=True,
      classical_constraints_not_quantum_Ward=True,dust_optimization_stopped=True,
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
