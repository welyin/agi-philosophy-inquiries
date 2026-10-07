"""942 delivery audit, separate from recalculation of physics."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_942_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,942):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'dirac_material_embedding_results.json')
    assert r['round']==942 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert r['new_neutral_Dirac_flavours']==60 and r['Dirac_symbol_dimension']==240
    for k in ('positive_isometry_error','actual_Dirac_intertwiner_error','old_unknown_record_rest_isometry_error','geometry_source_projection_error','lapse_source_projection_error'):assert r[k]<2e-11,k
    assert r['changing_mass_source']['transported_source_error']<1e-6
    assert r['changing_mass_source']['source_gap_if_frame_omitted']>1e-7
    assert r['logical_record']['opposite_momentum_multiplier_difference']>.1
    assert not r['logical_record']['point_local_instrument_or_spatial_role_partition_proved']
    assert r['general_929_embedding']['required_neutral_Dirac_flavours']==4063232
    assert not r['general_929_embedding']['physical_spatial_agents_or_local_access_automatically_realized']
    assert r['actual_material_stress']['relative_covariant_vs_constraint_kernel_error']<1e-12
    for k in ('original_E_rec_action_unchanged','full_new_action_finite_gravity_evolution_or_matching_verified','full_goal_completed'):assert not r[k]
    note=STAGE/'research_note_942.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==20 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,11)]
    for phrase in ('整体目标未完成','原E_rec作用已改变','不是已经实现了多个空间分离的主体','4,063,232','未自动实现','不沿其种类数'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/material_action_decision.md',STAGE/'943/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至942）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=942)==list(range(1,943))
    assert '001—942轮共942份' in nav[0].read_text('utf-8-sig')
    assert '231—942的712份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '942：明确Dirac材料中的内部更新与几何来源' in s and '943/drafts/STATUS.md' in s
    oldfiles=[STAGE/'941/research_round_941_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'dirac_material_embedding.py',HERE/'dirac_material_embedding_results.json',HERE/'drafts/update_navigation942.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=942,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=942,
       cumulative_numbered_test_groups_from_941=previous['cumulative_numbered_test_groups_from_940']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       actual_local_matter_action_and_internal_process_connected=True,
       existing_FW_and_source_transport_not_reclaimed_as_new=True,
       new_species_and_mass_mixing_declared_inputs=True,
       logical_roles_distinguished_from_spatial_participants=True,
       full_929_material_encoding_cost_explicit=True,
       refinement_stopped_for_actual_organization_selection=True,
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
