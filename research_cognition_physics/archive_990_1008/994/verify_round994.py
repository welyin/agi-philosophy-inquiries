"""Delivery checks for charge history 994. Passing is not proof of complete physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout
TARGET=HERE/"research_round_994_checks.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,994):
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
    spec=importlib.util.spec_from_file_location('core994',HERE/'charge_history.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'charge_history_results.json');core.compare(core.run(),saved)
    assert saved['round']==994 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('actual_cosmological_rates_computed','initial_asymmetry_generated',
                'real_singlet_carries_B_minus_L','thermal_flavor_cases_are_actual_history'):
        assert not saved[key],key
    # Independent reduced scalar equations from Yukawa and charge definitions.
    q=F(7,237);h=-12*q/7;ls=-9*q;b=12*q;l=3*ls-3*h
    assert b-l==1 and b==F(saved['conversion_B_over_B_minus_L'])
    # Check integer/rational nullspaces with an independent numerical rank.
    rr=np.array([[float(F(x)) for x in row] for row in saved['reaction_rows']])
    wy=np.array([float(F(w)*F(y)) for w,y in zip(saved['susceptibility_weights'],saved['hypercharges'])])
    for branch in saved['neutral_equilibrium_branches']:
        rows=[r for r in rr]+[wy]
        for i,j in branch['fast_weinberg_pairs']:
            r=np.zeros(10);r[3+i]+=1;r[3+j]+=1;r[9]+=2;rows.append(r)
        matrix=np.array(rows)
        assert 10-int(np.linalg.matrix_rank(matrix))==branch['neutral_equilibrium_dimension']
        for vec in branch['basis_chemical_potentials']:
            assert np.max(abs(matrix@np.array([float(F(x)) for x in vec])))<1e-12
    # Partial-flavor result without importing the exact solver.
    q=F(7,257);h=-12*q/7;le=-h;lm=(4*q+h-1)/3;lt=(4*q+h)/3
    assert 9*q+le+lm+lt==0
    assert 12*q==F(saved['one_fast_flavor_example']['B']['exact'])
    assert saved['relaxation'][0]['rows'][-1]['B']>3.54e-4
    assert abs(saved['relaxation'][1]['rows'][-1]['B'])<1e-12
    note=STAGE/'research_note_994.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','转换关系','不是宇宙时间','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'boundary_adoption_v1.md',HERE/'charge_history.py',
        HERE/'charge_history_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish994.py',STAGE/'995/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至994）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=994)==list(range(1,995))
    assert '001—994轮共994份' in nav[0].read_text('utf-8-sig')
    assert '231—994的764份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '994：共同热史与净荷保存' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'993/research_round_993_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=994,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=994,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_993=prev['cumulative_numbered_test_groups_from_992']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,neutral_equilibrium_dimensions=[3,2,0],
        adopted_boundary='explicit flavor charges and protected reaction window; generation remains open',
        cosmological_rates_verified=False,initial_asymmetry_generated=False,
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
