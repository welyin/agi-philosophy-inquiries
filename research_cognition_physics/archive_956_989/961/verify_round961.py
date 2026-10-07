"""961 delivery audit, not a proof of the mechanism map or of unification."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_961_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,961):
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
    layout=Layout().verify();r=read(HERE/'operational_scale_screen_results.json')
    assert r['round']==961 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert all(r['decisions'].values())
    assert [round(x['echo_clock_ticks']) for x in r['frozen_transfer_cases']]==[5,5,10,10,5]
    assert max(x['all_unknown_payload_echo_isometry_error'] for x in r['frozen_transfer_cases'])<1e-13
    assert r['internal_relation_source_example']['total_current_balance_error']<1e-14
    scope=r['scope']
    assert scope['no_autonomous_growth_process_certified']
    assert scope['no_redshift_or_3D_metric_or_Friedmann_equation_derived']
    assert scope['concept_map_is_not_a_mathematical_proof']
    assert scope['stop_local_optimization'] and not scope['full_goal_completed']
    note=STAGE/'research_note_961.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==6 and re.findall(r'\\tag\{(\d+)\}',prose)==['1','2','3']
    departments=('量子：','时空：','物质：','相互作用与规范：','引力：','热现象与时间箭头：','宇宙演化：')
    for item in departments:assert item in prose,item
    for phrase in ('没有机制','机制相互排斥','机制明确但反日常经验','整体目标未完成',
                   '给定物理怎样实现认知协议','认知操作怎样解释物理结构',
                   '整体机制草图 → 概念冲突检查 → 最小数学检验 → 修订假说 → 系统验证'):
        assert phrase in prose,phrase
    assert prose.count('**最小检验。**')==7
    newdocs=[note,HERE/'drafts/mechanism_priority_entry.md',
             HERE/'drafts/deferred_material_gravity_interface.md',STAGE/'962/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至961）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=961)==list(range(1,962))
    assert '001—961轮共961份' in nav[0].read_text('utf-8-sig')
    assert '231—961的731份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig')
        assert '961：整体机制图与跨部门筛选优先' in s
        assert '961/drafts/mechanism_priority_entry.md' in s and '962/drafts/STATUS.md' in s
        assert '暂缓自动接续局部器件、控制、来源和误差优化' in s
    oldfiles=[STAGE/'960/research_round_960_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'operational_scale_screen.py',
        HERE/'operational_scale_screen_results.json',HERE/'drafts/update_navigation961.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=961,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=961,
       cumulative_numbered_test_groups_from_960=previous['cumulative_numbered_test_groups_from_959']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       all_seven_departments_have_roles_chain_resources_consequences_gaps_and_tests=True,
       mechanisms_and_physical_implementations_separated=True,
       three_kinds_of_mechanism_status_distinguished=True,
       expansion_calibration_screen_reproducible=True,
       former_local_gravity_task_deferred=True,
       concept_map_is_not_a_proof=True,
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
