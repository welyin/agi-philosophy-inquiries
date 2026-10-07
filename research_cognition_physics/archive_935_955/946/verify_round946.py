"""946 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_946_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,946):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'portal_common_process_results.json')
    assert r['round']==946 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    q=r['common_process'];b=r['finite_moving_control']
    assert q['gram_minimum_eigenvalue']>0
    assert min(q['higgs_instrument_gram_minimum_eigenvalues'])>0
    assert q['closed_form_comparison_error']<1e-12
    assert q['field_energy_plus_interaction_error']<1e-12
    assert b['higgs_moving_contrast_lower']>.0034
    assert b['record_moving_contrast_lower']>.021
    assert b['gravity_moving_contrast_lower']>.021
    assert b['state_error_GN']<.000501
    for table in q['joint_probability_tables']:
        assert min(table)>0 and abs(sum(table)-1)<1e-12
    assert r['old_945_reading_correction']['path0_vs_coherent_path_commutator_norm']==.5
    assert not r['scope']['nonlinear_SM_and_full_Einstein_matching_certified']
    assert not r['scope']['full_goal_completed']
    note=STAGE/'research_note_946.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==18 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,10)]
    for phrase in ('整体目标未完成','停止门户','0.0034096','不是已证明实际匹配误差','945第5节需限定'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/portal_decision.md',STAGE/'947/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至946）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=946)==list(range(1,947))
    assert '001—946轮共946份' in nav[0].read_text('utf-8-sig')
    assert '231—946的716份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '946：规范不变门户中的记录、径向物质与弱场引力' in s and '947/drafts/STATUS.md' in s
    oldfiles=[STAGE/'945/research_round_945_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'portal_common_process.py',HERE/'portal_common_process_results.json',HERE/'drafts/update_navigation946.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=946,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=946,
       cumulative_numbered_test_groups_from_945=previous['cumulative_numbered_test_groups_from_944']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       same_record_Newton_H_extended_with_two_scalar_portal=True,
       three_commuting_bounded_effects_jointly_verified=True,
       finite_mass_positive_lower_bounds_for_all_three_tasks=True,
       old_945_noncommuting_reading_wording_explicitly_qualified=True,
       full_nonlinear_parent_matching_not_claimed=True,
       physical_portal_and_representation_inputs_explicit=True,
       portal_refinement_stopped_after_decisive_test=True,
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
