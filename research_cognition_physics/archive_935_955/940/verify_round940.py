"""940 delivery and evidence checks, separate from scientific recalculation."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_940_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,940):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'native_exchange_bridge_results.json')
    assert result['round']==940 and result['all_scientific_checks_passed']
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['exact_rational_decomposition_error']=='0'
    assert result['maximum_Dirac_kernel_or_Ward_error']<1e-12
    assert result['unconserved_source_negative_control_gap']>.1
    assert result['static_unit_source']==dict(gravity_full=.5,gravity_TT_only=0.,electromagnetic_full=-1.,electromagnetic_transverse_only=0.)
    d=result['native_joint_diagnostic']
    assert d['dimension']==1536 and d['physical_bosonic_modes']==8
    assert d['norm_error']<1e-12 and d['energy_conservation_error']<1e-12
    assert d['exact_projected_force_identity_error']<1e-12
    assert d['actual_source_derivative_error']<1e-10
    assert d['gravity_record_difference']>4e-5
    assert d['largest_expected_occupation_boundary_force']>.01
    assert not d['independent_physical_continuum_error_bound_claimed']
    assert result['canonical_weak_field_methods_are_inherited_physics']
    assert result['finite_positive_process_and_own_sources_constructed']
    for k in ('full_original_curved_branch_or_SM_recovery_proved','cognition_generates_Einstein_or_dimension_proved','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_940.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==20 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,11)]
    for phrase in ('整体目标未完成','成熟弱场物理','不是原E_rec全场分支的数值解','没有证明原标量势在该点具有自洽平直真空','停止占据精度优化','原生H是本轮选用的工具'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/native_weak_field_decision.md',STAGE/'941/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至940）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in newdocs+nav:
        s=doc.read_text('utf-8-sig');content=re.sub(r'\$\$.*?\$\$','',s,flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=940)==list(range(1,941))
    assert '001—940轮共940份' in nav[0].read_text('utf-8-sig')
    assert '231—940的710份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig')
        assert '940：同一物质的弱场交换与正量子联合接口' in s and '941/drafts/STATUS.md' in s
    oldfiles=[STAGE/'939/research_round_939_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'native_exchange_bridge.py',HERE/'native_exchange_bridge_results.json',HERE/'drafts/update_navigation940.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=940,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=940,
       cumulative_numbered_test_groups_from_939=previous['cumulative_numbered_test_groups_from_938']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       standard_weak_field_physics_declared_as_input=True,
       full_stress_and_current_covariant_kernels_matched=True,
       positive_native_joint_H_and_own_sources_verified=True,
       finite_occupation_force_boundary_explicit=True,
       finite_model_not_claimed_full_physical_error_certificate=True,
       candidate_refinement_stopped_for_whole_model_selection=True,
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
