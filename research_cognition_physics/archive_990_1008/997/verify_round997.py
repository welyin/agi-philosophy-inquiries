"""Delivery checks for vacuum accessibility 997. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_997_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,997):
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
    spec=importlib.util.spec_from_file_location('core997',HERE/'vacuum_accessibility.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'vacuum_accessibility_results.json');core.compare(core.run(),saved)
    assert saved['round']==997 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('physical_vacuum_value_derived','early_boundary_generated',
                'actual_cosmological_fit','full_quantum_geometry_verified','full_goal_completed'):
        assert not saved[key],key
    assert not saved['ideal_causal_tasks']['actual_signal_device_verified']
    assert saved['conditional_future_window']['depends_on_future_continuation']
    # Independent derivative/source calculation from the common cell energy.
    ip=saved['inputs'];M=F(ip['material_mass']['exact']);R=F(ip['radiation_R']['exact'])
    V=F(ip['V0']['exact']);L=F(ip['vacuum_cell_L']['exact']);r=R/L
    assert M==L and F(ip['q']['exact'])==1 and F(ip['r']['exact'])==r
    for row in saved['source_rows']:
        a=F(row['a']);rho=(M/a**3+R/a**4+L)/V;p=(R/(3*a**4)-L)/V
        accel=-(rho+3*p)/(2*L/V)
        assert abs(float(accel)-row['acceleration_over_HLambda_squared'])<1e-15
        assert accel>0 and row['physical_cell_volume_ratio']==int(a**3)
    win=saved['finite_conformal_window'];lo=F(win['lower']['exact']);hi=F(win['upper']['exact'])
    assert F(6577,10000)<lo<hi<F(6582,10000)
    # Dense independent trapezoid is not the certificate; exact monotone sums are.
    x=np.linspace(.25,1,32769)
    estimate=float(np.trapezoid(1/np.sqrt(1+x**3+float(r)*x**4),x))
    assert float(lo)<estimate<float(hi)
    assert 2*F(1,5)<lo and F(11,20)<lo and 2*F(11,20)>1
    fut=saved['conditional_future_window']
    assert F(fut['upper']['exact'])==hi+F(1,4)<1
    assert F(21,10)>2*F(fut['general_upper'])
    # Source matching does not turn a material rest mass into a vacuum density.
    assert F(2)-F(2)*2**3==-14
    for row in saved['availability_same_scale']:assert row['vacuum_shift_residual']<8e-13
    note=STAGE/'research_note_997.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','未来假设','固定几何','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'cosmological_adoption_v1.md',HERE/'vacuum_accessibility.py',
        HERE/'vacuum_accessibility_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish997.py',STAGE/'998/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至997）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=997)==list(range(1,998))
    assert '001—997轮共997份' in nav[0].read_text('utf-8-sig')
    assert '231—997的767份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '997：有效Λ与实际可达范围' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'996/research_round_996_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=997,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=997,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_996=prev['cumulative_numbered_test_groups_from_995']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,finite_causal_window_certified_by_rational_bounds=True,
        adopted_interface='effective positive vacuum source with explicitly limited causal access and future assumptions',
        physical_vacuum_value_derived=False,actual_cosmological_future_predicted=False,
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
