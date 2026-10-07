"""Delivery checks for mechanism screen 991. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_991_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,991):
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
    import numpy as np
    spec=importlib.util.spec_from_file_location('core991',HERE/'dark_mode_screen.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'dark_mode_screen_results.json');core.compare(core.run(),saved)
    assert saved['round']==991 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert saved['scalar_is_new_explicit_physical_input']
    assert not saved['scalar_model_is_original']
    assert not saved['current_experimental_viability_claimed']
    assert not saved['relic_density_computed']
    assert F(saved['masses_squared']['dark']['exact'])==F(1,4)
    assert F(saved['analytic_source_energy_ratio_upper']['exact'])==F(1,500)
    # Independent larger representation verifies that moments do not use the edge.
    a=np.diag(np.sqrt(np.arange(1,12,dtype=float)),1)
    q=(a+a.T)/np.sqrt(1000.)
    residuals=[]
    for row in saved['rows']:
        n=row['n']
        residuals.append(abs((q@q)[n,n]-row['s2']['value']))
        residuals.append(abs(np.linalg.matrix_power(q,4)[n,n]-row['s4']['value']))
    assert max(residuals)<1e-13
    note=STAGE/'research_note_991.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','不是两种测量概率之比','没有识别现实暗物质',
                 '不是用8维截断认证长期量子场演化','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'dark_mode_screen.py',HERE/'dark_mode_screen_results.json',
        Path(__file__),HERE/'drafts/selection.md',HERE/'drafts/mechanism_map_v0_6.md',
        HERE/'drafts/publish991.py',STAGE/'992/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至991）',1)[1].split('### 当前取舍',1)[0]
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
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=991)==list(range(1,992))
    assert '001—991轮共991份' in nav[0].read_text('utf-8-sig')
    assert '231—991的761份' in nav[5].read_text('utf-8-sig')
    for p in nav:assert '991：共同几何与弱可见模式' in p.read_text('utf-8-sig')
    oldfiles=[STAGE/'990/research_round_990_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=991,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=991,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_990=prev['cumulative_numbered_test_groups_from_989']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_fock_moment_max_residual=max(residuals),
        mechanism_decision='retain weak visible access with a shared full geometric source as optional D991',
        original_dark_matter_model_claimed=False,physical_dark_matter_identified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
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
