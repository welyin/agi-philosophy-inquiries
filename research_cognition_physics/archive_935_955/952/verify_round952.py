"""952 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_952_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,952):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'portal_lapse_source_results.json')
    assert r['round']==952 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    c=r['checks'];scope=r['scope'];finite=r['finite_lapse']
    assert c['uniform_lapse_identity_error']<1e-13
    assert c['independent_real_space_bilinear_error']<1e-13
    assert c['independent_resolvent_variation_error']<1e-10
    assert c['exact_positive_remainder_identity_error']<1e-13
    assert c['stationary_energy_lapse_variation_error']<1e-11
    assert .000771 < finite['correct_record_W'] < .000773
    assert finite['wrong_minus_correct_record_W']>.0005
    assert finite['kernel_linearization_remainder_norm']<finite['kernel_remainder_bound']
    assert scope['same_948_record_W_and_fermion_kernel']
    assert scope['constant_lapse_is_only_normalization_check']
    assert not scope['full_Einstein_solution_or_full_stress_certified']
    assert not scope['new_947_finite_time_probability_or_preparation_bound']
    assert not scope['full_goal_completed']
    note=STAGE/'research_note_952.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for phrase in ('整体目标未完成','同阶','停止','不是探测概率','均匀N可以是整体时间单位变化'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/shared_source_decision.md',STAGE/'953/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至952）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=952)==list(range(1,953))
    assert '001—952轮共952份' in nav[0].read_text('utf-8-sig')
    assert '231—952的722份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '952：共同交换核的局部时间来源' in s and '953/drafts/STATUS.md' in s
    oldfiles=[STAGE/'951/research_round_951_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'portal_lapse_source.py',HERE/'portal_lapse_source_results.json',HERE/'drafts/update_navigation952.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=952,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=952,
       cumulative_numbered_test_groups_from_951=previous['cumulative_numbered_test_groups_from_950']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       same_portal_material_and_lapse_kernel_verified=True,
       endpoint_only_same_order_counterexample_verified=True,
       positive_finite_lapse_remainder_verified=True,
       static_state_not_confused_with_947_preparation=True,
       no_full_Einstein_or_detector_certificate_claimed=True,
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
