"""935 delivery: bounded autonomous source identity, without candidate repair."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse,math
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_935_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,935):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'autonomous_source_identity_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==935 and result['all_scientific_checks_passed']
    assert max(result['endpoint_isometry_errors'].values())<1e-12
    for row in result['two_branch_matrix_rows']:
        assert max(row['identity_errors'])<1e-12
        assert row['program_actual_rate_difference']>.1
    for row in result['full_input_energy_moments']:
        n=row['order']
        assert row['compiled_diagonal']==[.5,.5]
        assert row['direct_diagonal']==[.5,(.5**n+1.5**n)/2]
    for row in result['source_rows']:
        assert abs(row['program']-row['direct_rate'])<1e-12
        assert abs(row['program']-row['target_probability'])<1e-12
    witness=next(r for r in result['source_rows'] if r['source']==.1)
    assert .0904<witness['actual_rate']-witness['program']<.0905
    ds=result['derivative_checks']
    for row in ds:
        assert abs(row['actual_rate'])<1e-10
        assert abs(row['direct_rate']+math.pi/4)<1e-6
    assert all(3.9<a['program_error']/b['program_error']<4.1 for a,b in zip(ds,ds[1:]))
    assert max(r['conservation_error'] for r in result['energy_rows'])<1e-12
    for k in ('clock_baseline_accounted','exact_uniform_universe_time_rescaling_not_used_as_observable','pure_gate_dressed_independent_clock_class_only','clock_optimization_stopped'):assert result[k]
    for k in ('original_Qeff_rejected','all_autonomous_implementations_rejected','gravity_or_standard_model_restored','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_935.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==14 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,8)]
    for phrase in ('B_rate','不是GR或标准模型','同维正H','未结项','不要求H唯一','停止时钟扩展'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/joint_model_identity_review.md',STAGE/'936/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至935）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=935)==list(range(1,936))
    assert '001—935轮共935份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—935的705份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '935：自主实现的程序源与物理来源不能混用' in s and '936/drafts/STATUS.md' in s
    oldfiles=[STAGE/'934/research_round_934_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/unified_operation_hypotheses_v0_1.md']
    newfiles=newdocs+[Path(__file__),HERE/'autonomous_source_identity.py',HERE/'autonomous_source_identity_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'934/research_round_934_checks.json')
    out=dict(round=935,date='2026-10-07',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=935,
      cumulative_numbered_test_groups_from_934=previous['cumulative_numbered_test_groups_from_933']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      established_clock_and_finite_menu_results_reused=True,
      program_source_not_identified_with_actual_rate_source=True,
      reference_external_to_local_source_support_not_external_to_whole=True,
      internal_clock_energy_retained_in_positive_control=True,
      Brate_scoped_trial_contract_not_derived_universal_axiom=True,
      candidate_failure_not_promoted_to_program_failure=True,
      no_unique_H_or_parameter_requirement_added_to_existence_stage=True,
      autonomous_compiler_optimization_stopped=True,
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
