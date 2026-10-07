"""Delivery checks for 989. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_989_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,989):
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
    result=read(HERE/'shared_record_quantum_action_results.json')
    assert result['round']==989 and result['all_scientific_checks_passed']
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    import numpy as np
    from fractions import Fraction as F
    def load(name,path):
        spec=importlib.util.spec_from_file_location(name,path)
        m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
    core=load('core989',HERE/'shared_record_quantum_action.py')
    core.compare(core.run(),result)
    # Independent representation: 36D active CAR sector and 24D local relational code.
    old=load('core956',STAGE/'956/native_material_interface.py')
    mat=old.material(1.,.1)
    cs=[old.annihilate(j) for j in range(4)];ns=[c.T@c for c in cs]
    sector=np.eye(16)[:,[j for j in range(16) if j.bit_count()==2]]
    q=sector.T@(ns[0]+ns[1]-ns[2]-ns[3])@sector/2
    H=np.kron(mat['H'],np.eye(6))+np.kron(np.eye(6),mat['H'])+.2*np.kron(q,q)
    ev,vec=np.linalg.eigh(H)
    e=np.eye(16)
    code=np.stack([(e[5]-e[6]-e[9]+e[10])/2,
        (2*e[3]+2*e[12]-e[5]-e[6]-e[9]-e[10])/math.sqrt(12)],axis=1)
    local=np.kron(mat['W'],np.eye(4))@code
    phase=load('phase965',STAGE/'965/material_field_window.py').phases
    T=result['inputs']['time_from_958']
    targetlocal=local@np.diag([phase(np.array([-mat['J']]),T)[0],1])
    target=np.kron(targetlocal,targetlocal)@np.array([-1.,1.,1.,1.])/2
    initial=np.kron(local,local)@np.ones(4)/2
    def apply(u,state):
        x=state.reshape(6,4,6,4).transpose(0,2,1,3)
        x=(u@x.reshape(36,16)).reshape(6,6,4,4)
        return x.transpose(0,2,1,3).reshape(576)
    residuals=[]
    for row in result['rows']:
        u=(vec*phase(ev,row['time']))@vec.conj().T
        actual=apply(u,initial)
        fid=float(abs(np.vdot(target,actual))**2)
        residuals.append(abs(fid-row['native_fidelity']))
    assert max(residuals)<1e-8
    assert 29*97441**2<524737**2
    analytic=result['analytic']
    assert F(analytic['native_witness_lower']['exact'])==F(48147,50000)
    assert F(analytic['finite_prediction_gap_lower']['exact'])==F(23147,50000)
    assert F(analytic['unread_classical_flag_witness_lower']['exact'])==F(60647,100000)
    assert result['flag_commutator_norm']==0
    assert result['decision']['every_relation_label_must_be_coherent'] is False
    assert result['decision']['gravitational_quantization_proved'] is False
    note=STAGE/'research_note_989.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','不是设备无关Bell实验','不升级为宇宙最大化预测收益',
        '没有证明空间、几何、引力中介均必须量子化','整个区间下界来自解析误差'):
        assert term in prose,term
    newfiles=[note,HERE/'shared_record_quantum_action.py',HERE/'shared_record_quantum_action_results.json',
        Path(__file__),HERE/'drafts/mechanism_adoption_decision.md',HERE/'drafts/mechanism_map_v0_5.md',
        HERE/'drafts/publish989.py',STAGE/'990/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至989）',1)[1].split('### 当前取舍',1)[0]
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
        match=re.fullmatch(r'research_note_(\d+).md',p.name)
        if match and p.parent.name.startswith('archive_'):nums.append(int(match[1]))
    assert sorted(n for n in nums if n<=989)==list(range(1,990))
    assert '001—989轮共989份' in nav[0].read_text('utf-8-sig')
    assert '231—989的759份' in nav[5].read_text('utf-8-sig')
    for p in nav:assert '989：共享证据与不可替代的量子作用' in p.read_text('utf-8-sig')
    oldfiles=[STAGE/'988/research_round_988_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=989,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=989,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_988=prev['cumulative_numbered_test_groups_from_987']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_36D_CAR_fidelity_residuals=residuals,
        analytic_rational_window_bounds_verified=True,
        mechanism_decision='classical evidence may select a quantum contact but cannot replace its quantum resources',
        universal_prediction_optimization_adopted=False,full_goal_completed=False,
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
