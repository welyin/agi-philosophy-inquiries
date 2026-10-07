"""950 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_950_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,950):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'material_vacuum_matching_results.json')
    assert r['round']==950 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    c=r['checks'];loop=r['one_loop'];scope=r['scope']
    assert r['parameters']['literal_Dirac_flavors']==4063232
    assert c['actual_B_core_eigenvalue_error']<1e-12
    assert c['Dirac_projector_trace_error']<1e-13
    assert c['positive_quadrature_relative_agreement']<1e-12
    assert 1.28e-25 < loop['L4'] < 1.30e-25
    assert loop['entire_Q_window_common_response_relative_upper']<3.562e-23
    assert abs(loop['true_to_wrong_species_count_ratio']-4063232)<1e-7
    assert .0465 < r['local_matching_counterexample']['common_relative_response_change'] < .0466
    assert r['weak_family_diagnostics'][0]['collective_loop_factor_diagnostic']>6000
    assert scope['one_closed_fermion_loop_only'] and scope['subtraction_constants_are_matching_inputs']
    assert not scope['all_orders_certified']
    assert not scope['full_metric_source_or_parent_matching_certified']
    assert not scope['literal_flavor_candidate_selected_as_final']
    assert not scope['full_goal_completed']
    note=STAGE/'research_note_950.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==14 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,8)]
    for phrase in ('整体目标未完成','固定局部匹配','停止','4.65%','没有认证'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/material_loop_decision.md',STAGE/'951/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至950）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=950)==list(range(1,951))
    assert '001—950轮共950份' in nav[0].read_text('utf-8-sig')
    assert '231—950的720份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '950：记录材料的真空计数与有限匹配' in s and '951/drafts/STATUS.md' in s
    oldfiles=[STAGE/'949/research_round_949_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'material_vacuum_matching.py',HERE/'material_vacuum_matching_results.json',HERE/'drafts/update_navigation950.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=950,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=950,
       cumulative_numbered_test_groups_from_949=previous['cumulative_numbered_test_groups_from_948']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       actual_material_vacuum_flavor_count_verified=True,
       same_kernel_one_loop_spacelike_bound_verified=True,
       finite_local_matching_counterexample_verified=True,
       literal_flavor_candidate_demoted_without_program_no_go=True,
       no_all_order_or_parent_probability_certificate_claimed=True,
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
