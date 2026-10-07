"""Delivery checks for joint candidate 993. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_993_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,993):
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
    spec=importlib.util.spec_from_file_location('core993',HERE/'joint_matching.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'joint_matching_results.json');core.compare(core.run(),saved)
    assert saved['round']==993 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('native_higgs_matching_computed','prospective_budgets_physically_certified',
                'interacting_resonance_claimed','full_common_process_verified'):
        assert not saved[key],key
    old=read(STAGE/'991/dark_mode_screen_results.json')
    mh2=F(old['masses_squared']['higgs']['exact'])
    md2=F(old['masses_squared']['dark']['exact'])
    assert mh2==4*md2 and saved['free_pair_frequency_equals_higgs_mass']
    # Independent finite-difference mixed derivative of the explicit potential.
    mixed=[]
    for row in saved['static_source_rows']:
        j=F(row['ordinary_calibration_source']['exact'])
        d=F(row['dark_source_contrast']['exact'])
        z1=-j/mh2;z2=-(j+d)/mh2;z3=-d/mh2
        E=lambda z,source: mh2*z*z/2+source*z
        delta=E(z2,j+d)-E(z1,j)-E(z3,d)
        assert delta==F(row['cross_energy_density']['exact'])
        mixed.append(str(delta))
    oldthermal=read(STAGE/'992/cooling_availability_results.json')
    a2=next(x for x in oldthermal['rows'] if x['a']==2)
    th=saved['thermal_interface'];delta=th['prospective_state_trace_error']
    ee=read(STAGE/'980/finite_thermal_records_results.json')['parameters']['energies']
    fdelta=-delta*math.log(delta)-(1-delta)*math.log(1-delta)+delta*math.log(3)
    expected=(2*th['prospective_operator_error']+(max(ee)-min(ee))*delta
        +(a2['temperature']+th['prospective_temperature_error'])*fdelta
        +th['prospective_temperature_error']*math.log(4))
    assert abs(expected-th['error_bound'])<1e-12
    assert a2['availability']-expected>.03525724
    note=STAGE/'research_note_993.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==8
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,5)]
    for term in ('整体目标未完成','选择预算不是证明达标','复用','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'common_candidate_v1.md',HERE/'joint_matching.py',
        HERE/'joint_matching_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish993.py',STAGE/'994/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至993）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=993)==list(range(1,994))
    assert '001—993轮共993份' in nav[0].read_text('utf-8-sig')
    assert '231—993的763份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '993：共同候选与跨部门匹配' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'992/research_round_992_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=993,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=993,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_992=prev['cumulative_numbered_test_groups_from_991']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        conditional_availability_lower_bound=a2['availability']-expected,
        candidate_version='B993 v1: P981 + D991; retain common Higgs dynamics',
        matching_budgets_certified=False,full_common_process_verified=False,
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
