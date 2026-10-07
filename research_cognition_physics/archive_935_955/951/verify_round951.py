"""951 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_951_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,951):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'body_sector_bridge_results.json')
    assert r['round']==951 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    c=r['checks'];scope=r['scope']
    assert r['parameters']['actual_947_internal_dimension']==4063232
    assert c['H_intertwining_error']<1e-12
    assert c['empty_sector_generator_error']<1e-12
    assert c['full_unitary_intertwining_error']<2e-12
    assert max(c['physical_source_intertwining_errors'].values())<2e-10
    assert c['finite_material_field_exchange_current_norm']>1
    assert c['number_changing_negative_control_vacuum_leak']>1e-7
    assert r['inherited_by_exact_sector_identity']['bounds_are_frozen_947_not_new_numerical_estimates']
    assert scope['exact_fixed_sector_representation_theorem']
    assert scope['high_energy_matching_coefficients_must_be_kept']
    assert not scope['microscopic_binding_or_SM_material_construction_proved']
    assert not scope['full_SM_Einstein_parent_matching_certified']
    assert not scope['full_goal_completed']
    note=STAGE/'research_note_951.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for phrase in ('整体目标未完成','源无关等距','停止','高能匹配没有消失','不是完整父场到有效物体的匹配误差'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/joint_minimum_menu.md',STAGE/'952/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至951）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=951)==list(range(1,952))
    assert '001—951轮共951份' in nav[0].read_text('utf-8-sig')
    assert '231—951的721份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '951：有限物体内部状态的保源表示' in s and '952/drafts/STATUS.md' in s
    oldfiles=[STAGE/'950/research_round_950_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'body_sector_bridge.py',HERE/'body_sector_bridge_results.json',HERE/'drafts/update_navigation951.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=951,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=951,
       cumulative_numbered_test_groups_from_950=previous['cumulative_numbered_test_groups_from_949']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       exact_material_sector_and_source_identity_verified=True,
       empty_body_sector_not_full_vacuum_subtraction=True,
       full_947_errors_inherited_not_reestimated=True,
       shared_minimum_menu_preserves_all_target_sectors=True,
       microscopic_material_and_full_parent_not_certified=True,
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
