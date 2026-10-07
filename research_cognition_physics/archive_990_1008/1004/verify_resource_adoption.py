"""Verify working resource adoption; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'drafts/resource_adoption_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,1004):
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
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()







    import importlib.util
    spec=importlib.util.spec_from_file_location('resourceaudit',HERE/'resource_preparation_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    result=read(HERE/'resource_preparation_audit_results.json');assert core.run()==result
    assert result['all_audit_checks_passed'] and result['new_scientific_test_groups']==0
    assert result['formal_reports']==1003 and result['cumulative_numbered_test_groups']==3785
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert len(result['resource_classification'])==3
    assert not any(row['natural_supply_certified'] for row in result['resource_classification'])
    assert not result['full_energy_conservation_alone_implies_bare_covariance']
    assert result['covariance_requires_additive_energy_and_invariant_joint_input']
    assert not result['spatial_reference_522_523_generates_time_asymmetry']
    assert not result['natural_preparation_of_1002_auxiliary_certified']
    assert not result['full_goal_completed']
    note=STAGE/'research_note_1004_working.md'
    for term in ('正式仍1003／累计3785','新增物理试验0','整体目标未完成','H₀+V'):
        assert term in note.read_text('utf-8'),term
    files=[note,HERE/'resource_preparation_adoption_v1.md',HERE/'resource_preparation_audit.py',
        HERE/'resource_preparation_audit_results.json',HERE/'drafts/adoption_decision.md',
        HERE/'drafts/NEXT.md',HERE/'drafts/publish_resource_adoption.py',Path(__file__)]
    for pth in files:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in [p for p in files if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    assert '001—1003轮共1003份' in nav[0].read_text('utf-8-sig')
    assert '231—1003的773份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1004工作审计：任务资源与关系参考' in pth.read_text('utf-8-sig')
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1003）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    evidence=[STAGE/'1003/research_round_1003_checks.json',HERE/'drafts/STATUS.md']
    out=dict(date='2026-10-07',kind='working_resource_adoption_audit',
        all_audit_checks_passed=True,formal_reports=1003,cumulative_numbered_test_groups=3785,
        new_scientific_test_groups=0,historical_unique_files_verified=len(frozen),
        historical_manifest_evidence=layout,local_links_checked=links,
        app_goal_changed=False,joint_lifecycle_certified=False,full_goal_completed=False,
        visual_checks_performed=False,
        frozen_evidence_hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence},
        audit_file_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
    if not writing:
        before=read(TARGET)
        for key in ('frozen_evidence_hashes','audit_file_hashes'):assert before[key]==out[key],key
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
        ('frozen_evidence_hashes','audit_file_hashes')},ensure_ascii=False,indent=2))
