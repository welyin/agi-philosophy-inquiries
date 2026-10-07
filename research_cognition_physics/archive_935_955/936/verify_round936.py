"""936 delivery: native effective generator, source and canonical force."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse,math
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_936_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,936):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'native_source_quantization_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==936 and result['all_scientific_checks_passed']
    assert abs(result['original873_normalized_mixed_coefficient']-.25)<1e-14
    assert max(result['quadrature_comparison'].values())<1e-8
    rows=result['finite_process_rows']
    assert [r['dimension'] for r in rows]==[32,48,80,128]
    for r in rows:
        assert r['boundary_corrected_force_error']<1e-10
        assert r['energy_conservation_error']<1e-12
        assert r['matter_purity']<.99999 and abs(r['canonical_momentum_change'])>.01
        assert r['minimum_H_eigenvalue']>-1e-12
    assert rows[0]['projection_boundary_force_norm']>1e-5
    assert rows[-1]['analytic_generator_tail_bound']<1e-8
    ds=result['source_derivative_rows']
    for key in ('first_error','mixed_error'):assert 3.8<ds[0][key]/ds[1][key]<4.2
    assert abs(result['record_source_derivative']-result['record_source_finite_difference'])<1e-7
    for k in ('same_native_H_for_evolution_source_and_force','source_independent_chart_and_cutoffs_required','calibration_symbol_is_not_full_Einstein_SM','full_physical_symbol_and_matching_still_required','quadrature_comparison_is_not_a_rigorous_integration_certificate'):assert result[k]
    for k in ('original_Q_or_E_rejected','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_936.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==22 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,12)]
    for phrase in ('不是原Einstein—SM的模拟','不是空间维数','尚不是完整共同物理模型','整体未结项','不是严格求积证书'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/native_model_reuse_decision.md',STAGE/'937/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至936）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=936)==list(range(1,937))
    assert '001—936轮共936份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—936的706份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '936：原生有效生成元共同保来源与几何力' in s and '937/drafts/STATUS.md' in s
    oldfiles=[STAGE/'935/research_round_935_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'native_source_quantization.py',HERE/'native_source_quantization_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'935/research_round_935_checks.json')
    out=dict(round=936,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=936,
      cumulative_numbered_test_groups_from_935=previous['cumulative_numbered_test_groups_from_934']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      native_generator_not_endpoint_program_used=True,
      same_quantization_for_force_energy_and_source=True,
      known_quantization_tools_attributed=True,
      finite_mode_symbol_class_not_claimed_to_be_full_physical_model=True,
      canonical_projection_boundary_retained=True,
      integration_convergence_not_claimed_to_be_rigorous_error_certificate=True,
      old_spatial_and_source_results_reused=True,
      complete_physical_symbol_and_effective_matching_remain_open=True,
      compiler_and_calibration_optimization_stopped=True,
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
