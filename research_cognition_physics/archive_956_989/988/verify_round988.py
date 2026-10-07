"""Delivery checks for 988. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_988_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,988):
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
    result=read(HERE/"relation_withdrawal_results.json")
    assert result["round"]==988 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    import numpy as np
    spec=importlib.util.spec_from_file_location("core988",HERE/"relation_withdrawal.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    # Independent small-time propagation: scaled Taylor rather than eigh.
    I=np.eye(4)
    h=np.array([[0,-.2,0,0],[-.2,1,0,0],[0,0,1,0],[0,0,0,0.]])
    q=np.zeros((4,4));q[1,2]=q[2,1]=1
    f=np.zeros((4,4));f[0,3]=f[3,0]=-math.sqrt(3)/4;f[3,3]=-.5
    k=lambda a,b,c:np.kron(np.kron(a,b),c)
    H0=k(h,I,I)+k(I,h,I)+k(I,I,h)+.2*(k(q,q,I)+k(I,q,q))
    H1=H0+.05*k(I,f,I)
    def taylor(H,t):
        s=max(0,math.ceil(math.log2(max(1.,2*abs(t)*np.linalg.norm(H,1)))))
        A=-1j*t*H/(2**s);out=np.eye(len(H),dtype=complex);term=out.copy()
        for j in range(1,33):
            term=term@A/j;out+=term
        trunc=(2**s)*math.exp(.5)*(.5**33)/math.factorial(33)
        for _ in range(s):out=out@out
        return out,trunc
    U,b1=taylor(H1,20);V,b2=taylor(H0,5)
    base=np.zeros(64);base[16*3+4*3+1]=1
    formed=U@base
    X=I.copy();X[:,[1,3]]=X[:,[3,1]];XA=k(X,I,I)
    EC=k(I,I,np.diag([0,0,1,0]))
    probs=[float(np.vdot(z,EC@z).real) for z in (V@formed,V@XA@formed)]
    residuals=[abs(p-r['receiver_probability']) for p,r in zip(probs,result['stop_update_only'])]
    assert max(residuals)<1e-12
    # Exact Python-integer computation avoids relying on fixed-width arithmetic.
    K=np.rint(H0*5).astype(object)
    K=np.array([[int(x) for x in row] for row in K],dtype=object)
    E=EC.astype(int).astype(object);XI=XA.astype(int).astype(object)
    for order in range(1,4):
        E=K@E-E@K;C=E@XI-XI@E
        if order<3:assert np.count_nonzero(C)==0
        else:assert C[17,38]==1 and np.count_nonzero(C)==16
    assert result['positive_contrast_lower']>0.00078574
    assert result['park_exact_future_contrast']==0
    assert result['all_input_free_transport_defect']<1e-12
    assert result['resources']['includes_autonomous_control_and_battery'] is False
    note=STAGE/'research_note_988.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','没有把转存器件升级为整个研究纲领的门槛',
        '新增操作输入','未交付一个能量完全闭合的自主装置','转存不是唯一必要机制'):
        assert term in prose,term
    newfiles=[note,HERE/'relation_withdrawal.py',HERE/'relation_withdrawal_results.json',
        Path(__file__),HERE/'drafts/withdrawal_decision.md',HERE/'drafts/mechanism_map_v0_4.md',
        HERE/'drafts/publish988.py',STAGE/'989/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至988）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=988)==list(range(1,989))
    assert '001—988轮共988份' in nav[0].read_text('utf-8-sig')
    assert '231—988的758份' in nav[5].read_text('utf-8-sig')
    for p in nav:assert '988：关系撤销与内部资料保存' in p.read_text('utf-8-sig')
    oldfiles=[STAGE/'987/research_round_987_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=988,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=988,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_987=prev['cumulative_numbered_test_groups_from_986']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_taylor_probability_residuals=residuals,
        exact_python_integer_check=True,taylor_analytic_truncation=[b1,b2],
        mechanism_decision='stopping formation is insufficient; internal relocation works under an explicit contract',
        autonomous_control_and_battery_implemented=False,full_goal_completed=False,
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
