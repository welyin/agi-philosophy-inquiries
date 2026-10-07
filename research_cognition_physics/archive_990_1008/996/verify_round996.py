"""Delivery checks for cosmic CP window 996. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_996_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,996):
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
    spec=importlib.util.spec_from_file_location('core996',HERE/'cosmic_cp_window.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'cosmic_cp_window_results.json');core.compare(core.run(),saved)
    assert saved['round']==996 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('thermal_source_coefficient_evaluated','finite_physical_net_charge_certified',
                'actual_thermal_hierarchy_verified','cosmological_abundance_fitted','full_goal_completed'):
        assert not saved[key],key
    state=saved['phase_averaged_quadratic_preparation']
    assert state['constant_coefficient_shifted_together']
    assert not state['quartic_and_thermal_corrections_certified']
    assert F(state['occupation']['exact'])==F(9,16)
    assert F(state['f_excitation_initial']['exact'])==F(9,800000)
    assert F(state['initial_vacuum_reference_shift']['exact'])==F(1,100000)
    # Exact initial full-reference identity and trace polynomial.
    fi=F(9,800000);fv=F(1,100000);r2=F(233,2500)
    assert fi+fv==F(17,800000)
    shared=saved['shared_coupling'];A=F(shared['norm_constant']['exact'])
    B=F(shared['norm_a_minus3']['exact']);D=F(shared['norm_a_minus6']['exact'])
    gi=F(shared['norm_initial']['exact'])
    assert A==F(7,50)+fv/F(100)+fv*fv*r2
    assert B==fi/F(100)+2*fi*fv*r2 and D==fi*fi*r2 and gi==A+B+D
    window=saved['finite_window'];W=F(window['optical_depth']['exact'])
    # Independent coefficient integration on [1,2], exact arithmetic.
    assert W==(A/2+B*F(15,64)+D*F(127,896))/(10*gi)
    assert F(window['unwashed_shape']['exact'])==F(31,160)
    lower=math.exp(-float(W))*float(F(31,160));y=window['normalized_response']
    assert 0<lower<=y<=float(F(31,160))
    expected=3*float(F(17,200)*fi/gi/F(10))*y
    assert abs(window['yield_per_abs_alpha_over_beta_times_H_over_T']-expected)<1e-18
    assert window['independent_ODE_residual']<1e-12
    assert window['quadrature_comparison']<1e-13
    # Derivative of cumulative washout equals minus the rate/H/a, independently.
    for a in (1.1,1.4,1.8):
        h=1e-5
        def integral(x):return float((A*(F(1)/F(str(x))-F(1,2))
            +B*(F(str(x))**-4-F(1,16))/4+D*(F(str(x))**-7-F(1,128))/7)/(10*gi))
        deriv=(integral(a+h)-integral(a-h))/(2*h)
        rhs=-.1*(float(A)+float(B)*a**-3+float(D)*a**-6)/(float(gi)*a*a)
        assert abs(deriv-rhs)<1e-10
    assert F(saved['counter_kernel_first_moment']['exact'])==0
    note=STAGE/'research_note_996.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','一维筛选方程','同步取C_ref','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'mechanism_map_v0_8.md',HERE/'cosmic_cp_window.py',
        HERE/'cosmic_cp_window_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish996.py',STAGE/'997/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至996）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=996)==list(range(1,997))
    assert '001—996轮共996份' in nav[0].read_text('utf-8-sig')
    assert '231—996的766份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '996：膨胀包络与生成保存约束' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'995/research_round_995_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=996,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=996,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_995=prev['cumulative_numbered_test_groups_from_994']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,response_verified_with_shared_coupling=True,
        adopted_interface='conditional nonperiodic envelope with shared source-washout parameters and H/T suppression',
        thermal_source_coefficient_evaluated=False,finite_physical_net_charge_certified=False,
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
