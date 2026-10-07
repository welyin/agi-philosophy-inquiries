"""Delivery checks for internal CP bridge 995. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_995_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,995):
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
    spec=importlib.util.spec_from_file_location('core995',HERE/'internal_cp_bridge.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'internal_cp_bridge_results.json');core.compare(core.run(),saved)
    assert saved['round']==995 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('actual_thermal_kernel_computed','net_lepton_charge_generated',
                'expanding_cosmological_history_verified','full_goal_completed'):
        assert not saved[key],key
    assert saved['added_fundamental_species']==0 and saved['added_effective_operator_dimension']==7
    # Independent exact diagonal contraction, rather than reuse of the trace function.
    invariant=F(1,10)*F(1,5)+F(1,5)*F(1,10)+F(3,10)*F(3,20)
    assert invariant==F(saved['invariant']['exact'])==F(17,200)
    assert abs(saved['CP_conjugate_invariant']+float(invariant))<1e-15
    # Exact squeezed-state moments and quartic correction; no Fock cutoff.
    state=saved['parity_even_initial_state']
    q2=F(state['canonical_Q2']['exact']);p2=F(state['canonical_P2']['exact'])
    assert q2*p2==F(1,4) and F(state['canonical_Q4']['exact'])==3*q2*q2
    assert F(state['energy_including_quartic']['exact'])==F(3400003,6400000)
    assert F(state['f_second_derivative_including_quartic']['exact'])==F(1499997,80000000000)>0
    assert state['mean_sigma']==0
    rows=saved['quadratic_control']['rows']
    assert rows[0]['f']==rows[-1]['f'] and rows[2]['f']>rows[0]['f']
    assert rows[2]['two_time_cp_factor']>0
    # Independent multi-harmonic periodic-shift check of the cycle identity.
    angles=2*np.pi*np.arange(2048)/2048
    def periodic(phi):return np.sin(phi)+.37*np.cos(3*phi)-.21*np.sin(5*phi)
    periodic_shift_residual=max(abs(float(np.mean(periodic(angles-tau)-periodic(angles))))
                                for tau in (.137,1.23,7.11))
    assert periodic_shift_residual<1e-14
    assert abs(saved['causal_periodic_control']['cycle_integral'])<1e-16
    assert saved['causal_periodic_control']['pointwise_max']>1e-6
    note=STAGE/'research_note_995.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==14
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,8)]
    for term in ('整体目标未完成','不是实际Weinberg热碰撞核','不含四次项','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'internal_weinberg_adoption.md',HERE/'internal_cp_bridge.py',
        HERE/'internal_cp_bridge_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish995.py',STAGE/'996/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至995）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=995)==list(range(1,996))
    assert '001—995轮共995份' in nav[0].read_text('utf-8-sig')
    assert '231—995的765份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '995：内部变化与CP偏置' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'994/research_round_994_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=995,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=995,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_994=prev['cumulative_numbered_test_groups_from_993']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,periodic_shift_max_residual=periodic_shift_residual,
        adopted_interface='Z2-even internal Weinberg modulation with explicit dimension-seven input',
        actual_thermal_kernel_computed=False,net_lepton_charge_generated=False,
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
