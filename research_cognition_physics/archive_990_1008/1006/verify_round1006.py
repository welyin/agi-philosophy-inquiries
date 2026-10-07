"""Verify superconductivity adoption 1006; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_1006_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,1006):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for key in ('frozen_inputs','new_scientific_and_entry_files'):add(old[key])
    for path in ('875/drafts/clock_transport_working_checks.json',
        '878/drafts/working_checks.json','884/drafts/working_checks.json',
        '887/drafts/working_checks.json','894/drafts/working_checks.json',
        '896/drafts/working_checks.json','897/drafts/working_checks.json',
        '899/drafts/working_checks.json','900/drafts/working_checks.json',
        '904/drafts/working_checks.json','907/drafts/working_checks.json',
        '909/drafts/working_checks.json','912/drafts/working_checks.json',
        '913/drafts/working_checks.json','915/drafts/working_checks.json',
        '917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    extra={
        '909/drafts/scope_reaudit_checks.json':('preserved_files','evidence_and_audit_hashes'),
        '911/drafts/finite_scope_checkpoint_checks.json':('evidence_and_audit_hashes',),
        '912/drafts/effective_scope_reaudit_checks.json':('evidence_and_audit_hashes',),
        '914/drafts/finite_scope_decision_checks.json':('evidence_hashes','new_document_and_verifier_hashes'),
        '919/drafts/effective_scope_after_918_checks.json':('frozen_inputs_verified','new_document_and_verifier_hashes'),
        '953/drafts/priority_reaudit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '981/drafts/common_adoption_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '985/drafts/common_model_adoption_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '987/drafts/common_overlap_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '998/drafts/boundary_adoption_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes')}
    extra['1004/drafts/resource_adoption_checks.json']=('frozen_evidence_hashes','audit_file_hashes')
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()









    import importlib.util
    spec=importlib.util.spec_from_file_location('mechanism1006',HERE/'superconductivity_adoption_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'superconductivity_adoption_results.json')
    assert core.run()==saved
    assert saved['round']==1006 and saved['kind']=='mechanism_adoption_audit'
    assert saved['all_audit_checks_passed'] and saved['fresh_physical_test_groups']==0
    assert saved['cumulative_test_groups']==3786
    assert tuple(saved['evidence_obligations'])==(
        '参与者与过程','中间机制','资源、记录与反作用','输入与范围')
    assert saved['evidence_obligation_count']==len(saved['evidence_obligations'])==4
    assert all(saved['evidence_obligations'].values())
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('new_cognitive_axioms','microscopic_pairing_derived_from_cognition',
        'real_material_Tc_predicted','full_quantum_record_lifecycle_certified',
        'physical_unification_certified','full_goal_completed','visual_checks_performed'):
        assert saved[key] is False,key
    note=STAGE/'research_note_1006.md'
    prose=note.read_text('utf-8-sig')
    assert prose.count('$$')>0 and prose.count('$$')%2==0
    for term in ('整体目标未完成','累计3786','新增物理试验组0','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'conventional_superconductivity_adoption_v1.md',
        HERE/'superconductivity_adoption_audit.py',
        HERE/'superconductivity_adoption_results.json',Path(__file__),
        HERE/'drafts/adoption_decision.md',HERE/'drafts/review_notes.md',
        HERE/'drafts/publish1006.py',STAGE/'1007/drafts/STATUS.md']
    assert len(newfiles)==len(set(newfiles))==9
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8-sig'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1006）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for pth in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',pth.name)
        if m and pth.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=1006)==list(range(1,1007))
    assert '001—1006轮共1006份' in nav[0].read_text('utf-8-sig')
    assert '231—1006的776份' in nav[5].read_text('utf-8-sig')
    for pth in nav:
        navigation=pth.read_text('utf-8-sig')
        assert navigation.count('## 1006：常规超导与集体相干')==1,pth
        assert navigation.index('## 1006：常规超导与集体相干')<navigation.index(
            '## 1005：辐射反馈与整体假说v2.1'),pth
    oldfiles=[STAGE/'1005/research_round_1005_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    assert prev['all_delivery_checks_passed']
    assert prev['cumulative_numbered_test_groups_from_1004']==3786
    out=dict(round=1006,date='2026-10-08',all_delivery_checks_passed=True,
        formal_reports=1006,fresh_test_groups=0,
        cumulative_numbered_test_groups_from_1005=prev['cumulative_numbered_test_groups_from_1004'],
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,mechanism_evidence_obligations=len(saved['evidence_obligations']),
        base_mechanism_version='v2.1',mechanism_adoption='P08_conventional_superconductivity',
        superconductivity_solver_run=False,new_cognitive_axioms=False,
        microscopic_pairing_derived_from_cognition=False,real_material_Tc_predicted=False,
        full_quantum_record_lifecycle_certified=False,
        physical_unification_certified=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for key in ('frozen_inputs','new_scientific_and_entry_files'):assert before[key]==out[key],key
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run(args.write)
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
