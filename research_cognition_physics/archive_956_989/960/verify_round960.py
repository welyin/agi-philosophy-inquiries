"""960 delivery audit; science recalculation has a separate readonly entry."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_960_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,960):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'953/drafts/priority_reaudit_checks.json');add(au['frozen_evidence_hashes']);add(au['audit_file_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'material_dispersion_bridge_results.json')
    assert r['round']==960 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    a=r['actual_material_response'];scope=r['scope']
    assert a['maximum_spectral_error']<1e-13
    assert abs(a['chi_from_36_state_spectral_sum']-a['chi_second_order_958'])<1e-15
    assert a['old_read_time_times_difference']>.029
    assert r['finite_geometry_comparison']['relative_energy_mismatch']>.17
    assert r['analytic_effective_domain']['force_relative_error_upper']<.00015
    assert scope['same_material_leading_stationary_Maxwell_matching']
    assert scope['record_phase_response_and_distance_force_share_coefficients']
    assert scope['Maxwell_3plus1_and_dipole_coupling_are_physical_inputs']
    assert not scope['full_real_time_958_bound_transported']
    assert not scope['full_parent_SM_Einstein_certified']
    assert not scope['six_protocols_jointly_certified']
    assert not scope['full_goal_completed']
    assert scope['stop_dispersion_branch_after_this_interface']
    note=STAGE/'research_note_960.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('整体目标未完成','Maxwell','停止','保留','17.20%','物理输入'):
        assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/material_propagation_decision.md',STAGE/'961/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至960）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=960)==list(range(1,961))
    assert '001—960轮共960份' in nav[0].read_text('utf-8-sig')
    assert '231—960的730份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig')
        assert '960：原生材料的电磁色散接口与有限有效域' in s and '961/drafts/STATUS.md' in s
    oldfiles=[STAGE/'959/research_round_959_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'material_dispersion_bridge.py',
        HERE/'material_dispersion_bridge_results.json',HERE/'drafts/update_navigation960.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=960,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=960,
       cumulative_numbered_test_groups_from_959=previous['cumulative_numbered_test_groups_from_958']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       actual_material_polarizability_and_958_leading_energy_matched=True,
       same_Maxwell_energy_and_distance_source_checked=True,
       finite_effective_and_high_frequency_tail_bounds_recorded=True,
       original_real_time_bound_not_claimed=True,
       existing_geometry_mismatch_explicit=True,
       stop_dispersion_branch_and_return_to_global_adoption=True,
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
        with TARGET.open('x',encoding='utf-8') as fobj:
            json.dump(out,fobj,ensure_ascii=False,indent=2);fobj.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
