"""Delivery checks for cooling mechanism 992. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_992_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,992):
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
    import numpy as np
    spec=importlib.util.spec_from_file_location('core992',HERE/'cooling_availability.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'cooling_availability_results.json');core.compare(core.run(),saved)
    assert saved['round']==992 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('physical_reset_implemented','work_extraction_implemented',
                'initial_cosmological_boundary_explained','full_quantum_geometry_verified'):
        assert not saved[key],key
    assert saved['additional_blank_registers_created']==0
    # Original rational/radical four-level matrix, independent of the material importer.
    h=np.array([[0.,-.2,0.,-math.sqrt(3)/80],[-.2,1.,0.,0.],
                [0.,0.,1.,0.],[-math.sqrt(3)/80,0.,0.,-.025]])
    ev,u=np.linalg.eigh(h);beta0=2*math.log(2)
    p=np.exp(-beta0*ev);p/=sum(p)
    actual_entropy=float(-p@np.log(p));e0=float(p@ev)
    errors=[]
    for row in saved['rows']:
        beta=beta0*row['a'];q=np.exp(-beta*ev);z=sum(q);q/=z
        # F(tau) = -T log Z, with unshifted physical h.
        A=e0-actual_entropy/beta+math.log(z)/beta
        errors.append(abs(A-row['availability']))
        assert abs(float(np.sum(abs(p-q))/2)-row['equilibrium_state_trace_distance'])<1e-12
    assert max(errors)<1e-12
    assert abs(saved['expansion_availability_limit']-(e0-min(ev)))<1e-12
    note=STAGE/'research_note_992.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for term in ('整体目标未完成','不是材料已经达到的状态','不自动提取工作',
                 '两者还不是已经共同实现的一次宇宙记录过程','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'cooling_availability.py',HERE/'cooling_availability_results.json',
        Path(__file__),HERE/'drafts/selection.md',HERE/'drafts/mechanism_map_v0_7.md',
        HERE/'drafts/publish992.py',STAGE/'993/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至992）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [pth for pth in newfiles if pth.suffix=='.md']+nav:
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
    assert sorted(n for n in nums if n<=992)==list(range(1,993))
    assert '001—992轮共992份' in nav[0].read_text('utf-8-sig')
    assert '231—992的762份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '992：差异冷却与非平衡资源' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'991/research_round_991_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=992,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=992,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_991=prev['cumulative_numbered_test_groups_from_990']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_free_energy_max_residual=max(errors),
        mechanism_decision='retain differential scaling as athermality formation; keep actual use and boundary selection separate',
        cosmological_arrow_derived=False,full_quantum_geometry_verified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(pth.relative_to(ROOT)):sha(pth) for pth in oldfiles},
        new_scientific_and_entry_files={str(pth.relative_to(ROOT)):sha(pth) for pth in newfiles})
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
