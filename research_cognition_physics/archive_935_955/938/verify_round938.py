"""938 delivery: physical reference modes and the scope of quantum reuse."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_938_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,938):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'reference_physical_mode_audit_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==938 and result['all_scientific_checks_passed']
    assert result['physical_phase_counts']['0']['additional_physical_pairs']==4
    assert result['physical_phase_counts']['1']['additional_physical_pairs']==4
    rows=result['full_field_finite_wavelength_rows']
    assert [(r['N'],r['amplitude']) for r in rows]==[(16,.0003),(16,.00015),(24,.0003),(24,.00015)]
    for r in rows:
        assert r['physical_coordinate_wavenumber']==1
        for q in r['signs']:
            assert q['native_minus_C_integral']>0
            assert q['min_proper_density']>.011 and q['min_radicand']>5e-5
            assert q['old_non_dust_spatial_constraint_defect']>1e-4
            assert q['source_derivative_error']<1e-7
            for key in ('total_Hamiltonian_density_residual','total_spatial_constraint_residual','full_Einstein_initial_residual','density_weighted_reference_stress_error'):assert q[key]<1e-12
            assert q['original_internal_Gauss_preserved']
    for i in (0,2):assert 3.9<rows[i]['coefficient_relative_error']/rows[i+1]['coefficient_relative_error']<4.1
    assert abs(rows[0]['signs'][0]['native_minus_C_integral']-rows[2]['signs'][0]['native_minus_C_integral'])<1e-12
    # Exact arithmetic check of the analytic sufficient window, separately
    # from the sampled initial data. Parameter bounds reuse 572/900.
    from fractions import Fraction as F
    tau0lo,tau0hi=F('1.732'),F('1.733')
    assert tau0lo**2<3 and tau0hi**2>F('3.003')
    rho_lo=(F('.02')*tau0lo+F('.0001'))/3
    assert rho_lo>F('.0115')
    amp=F('.0003'); tauhi=F('1.743')
    clower=F('.0115')-(2*tauhi*amp+amp*amp)/3
    dupper=F(2,3)*amp*4
    assert clower>F('.01115') and dupper==F('.0008')
    assert F('.01115')**2-dupper**2>0
    assert F('2.002')*F(1,2)**5-2*F('.08')*F(1,2)**-3<0
    for key in ('new_modes_are_not_just_removed_gauge_labels','finite_scale_generator_and_source_change_verified','old_quantum_results_cannot_be_transferred_without_new_physical_mapping','not_a_quantum_completion_or_a_no_go_for_dust','dust_candidate_demoted_from_priority_quantum_route','no_new_cognitive_axiom_or_spatial_dimension_proof'):assert result[key]
    assert not result['full_goal_completed']
    note=STAGE/'research_note_938.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('整体目标未完成','不是量子尘埃不可构造的证明','不是完整时间演化或制备方案','不表示753椭圆解精度提高','下一步优先在已有场表和共同量子分支上整合'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/quantum_recovery_reuse_audit.md',STAGE/'939/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至938）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=938)==list(range(1,939))
    assert '001—938轮共938份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—938的708份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '938：新增参考的物理模式与量子复用边界' in s and '939/drafts/STATUS.md' in s
    oldfiles=[STAGE/'937/research_round_937_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'reference_physical_mode_audit.py',HERE/'reference_physical_mode_audit_results.json',HERE/'drafts/update_navigation938.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'937/research_round_937_checks.json')
    out=dict(round=938,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=938,
      cumulative_numbered_test_groups_from_937=previous['cumulative_numbered_test_groups_from_936']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      old_physical_phase_space_and_perturbative_results_audited=True,
      four_extra_physical_pairs_counted=True,
      exact_rational_regular_window_checked=True,
      actual_finite_wavelength_constraints_and_native_source_checked=True,
      inherited_elliptic_precision_not_claimed_improved=True,
      no_quantum_probability_or_universal_no_go_claimed=True,
      new_reference_candidate_demoted_without_erasing_classical_results=True,
      unlimited_candidate_repair_not_adopted_as_goal=True,
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
