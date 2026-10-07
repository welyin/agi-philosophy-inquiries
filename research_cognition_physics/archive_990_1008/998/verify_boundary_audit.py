"""Working boundary-adoption audit checks. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"drafts/boundary_adoption_audit_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,998):
        old=read(STAGE/f"{n}/research_round_{n}_checks.json")
        for key in ("frozen_inputs","new_scientific_and_entry_files"):add(old[key])
    for path in ("875/drafts/clock_transport_working_checks.json",
        "878/drafts/working_checks.json","884/drafts/working_checks.json",
        "887/drafts/working_checks.json","894/drafts/working_checks.json",
        "896/drafts/working_checks.json","897/drafts/working_checks.json",
        "899/drafts/working_checks.json","900/drafts/working_checks.json",
        "904/drafts/working_checks.json","907/drafts/working_checks.json",
        "909/drafts/working_checks.json","912/drafts/working_checks.json",
        "913/drafts/working_checks.json","915/drafts/working_checks.json",
        "917/drafts/working_checks.json"):
        add(read(STAGE/path)["files"])
    extra={
      "909/drafts/scope_reaudit_checks.json":("preserved_files","evidence_and_audit_hashes"),
      "911/drafts/finite_scope_checkpoint_checks.json":("evidence_and_audit_hashes",),
      "912/drafts/effective_scope_reaudit_checks.json":("evidence_and_audit_hashes",),
      "914/drafts/finite_scope_decision_checks.json":("evidence_hashes","new_document_and_verifier_hashes"),
      "919/drafts/effective_scope_after_918_checks.json":("frozen_inputs_verified","new_document_and_verifier_hashes"),
      "953/drafts/priority_reaudit_checks.json":("frozen_evidence_hashes","audit_file_hashes")}
    extra["981/drafts/common_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    extra["985/drafts/common_model_adoption_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    extra["987/drafts/common_overlap_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()









    import importlib.util
    from fractions import Fraction as F
    spec=importlib.util.spec_from_file_location('audit998',HERE/'boundary_scope_audit.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    data=read(HERE/'boundary_scope_audit_results.json')
    assert core.run()==data
    for rel,digest in data['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    goal=read(HERE/'drafts/goal_alignment.json')
    assert goal['status']=='active' and not goal['application_goal_changed']
    assert len(goal['aligned_clauses'])==4
    # Independent source fractions from the old input values.
    late=read(STAGE/'997/vacuum_accessibility_results.json')['inputs']
    M,R,L=[F(late[k]['exact']) for k in ('material_mass','radiation_R','vacuum_cell_L')]
    for row in data['literal_same_window_diagnostic']:
        a=F(row['a']);m=M/a**3;r=R/a**4;total=m+r+L
        assert F(row['radiation_fraction']['exact'])==r/total
        assert F(row['matter_fraction']['exact'])==m/total
        assert F(row['vacuum_fraction']['exact'])==L/total
        assert float(r/total)<.021
    assert not data['joint_cosmological_history_certified']
    assert not data['whole_program_refuted']
    assert data['new_scientific_test_groups']==0 and data['formal_reports']==997
    note=STAGE/'research_note_998_working.md'
    prose=note.read_text('utf-8-sig')
    for term in ('正式997／累计3781保持','不新增正式科学轮次','不反证996或997',
                 '386或425','整体目标未完成','不重复改写应用目标'):
        assert term in prose,term
    assert prose.count('$$')==4
    newfiles=[note,HERE/'boundary_scope_audit.py',HERE/'boundary_scope_audit_results.json',
        HERE/'drafts/boundary_adoption_decision.md',HERE/'drafts/NEXT.md',
        HERE/'drafts/goal_alignment.json',HERE/'drafts/publish_boundary_audit.py',
        HERE/'drafts/publication_record.json',Path(__file__)]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8-sig'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
        '_shared/notes/unified_physics_condition_ledger_current.md')]
    title='## 998工作审计：共同边界与候选投入'
    for p in nav:assert p.read_text('utf-8-sig').count(title)==1,p
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至997）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix=='.md']+nav:
        text=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            dest=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert dest.exists() or (writing and dest==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(nums)==list(range(1,998))
    assert not (STAGE/'research_note_998.md').exists()
    assert '001—997轮共997份' in nav[0].read_text('utf-8-sig')
    assert '231—997的767份' in nav[5].read_text('utf-8-sig')
    preserved=[STAGE/'997/research_round_997_checks.json',HERE/'drafts/STATUS.md']
    out=dict(date='2026-10-07',kind='working_boundary_adoption_audit',
        all_audit_checks_passed=True,formal_reports=997,
        cumulative_numbered_test_groups=3781,new_scientific_test_groups=0,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,app_goal_already_aligned=True,app_goal_changed=False,
        joint_cosmological_history_certified=False,whole_program_refuted=False,
        full_goal_completed=False,visual_checks_performed=False,
        frozen_evidence_hashes={str(p.relative_to(ROOT)):sha(p) for p in preserved},
        audit_file_hashes={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:assert read(TARGET)==out
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
