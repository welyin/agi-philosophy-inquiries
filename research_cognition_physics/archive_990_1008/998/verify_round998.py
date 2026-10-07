"""Verify finite receiver 998; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_998_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,998):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for key in ('frozen_inputs','new_scientific_and_entry_files'):add(old[key])
    for path in ('875/drafts/clock_transport_working_checks.json',
        '878/drafts/working_checks.json','884/drafts/working_checks.json',
        '887/drafts/working_checks.json','894/drafts/working_checks.json',
        '896/drafts/working_checks.json','897/drafts/working_checks.json',
        '899/drafts/working_checks.json','900/drafts/working_checks.json',
        '904/drafts/working_checks.json','907/drafts/working_checks.json',
        '909/drafts/working_checks.json','912/drafts/working_checks.json',
        '913/drafts/working_checks.json','915/drafts/working_checks.json',
        '917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    extra={
        '909/drafts/scope_reaudit_checks.json':('preserved_files','evidence_and_audit_hashes'),
        '911/drafts/finite_scope_checkpoint_checks.json':('evidence_and_audit_hashes',),
        '912/drafts/effective_scope_reaudit_checks.json':('evidence_and_audit_hashes',),
        '914/drafts/finite_scope_decision_checks.json':('evidence_hashes','new_document_and_verifier_hashes'),
        '919/drafts/effective_scope_after_918_checks.json':('frozen_inputs_verified','new_document_and_verifier_hashes'),
        '953/drafts/priority_reaudit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '981/drafts/common_adoption_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '985/drafts/common_model_adoption_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '987/drafts/common_overlap_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes'),
        '998/drafts/boundary_adoption_audit_checks.json':('frozen_evidence_hashes','audit_file_hashes')}
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()

    import importlib.util
    from fractions import Fraction as F
    from decimal import Decimal, localcontext
    spec=importlib.util.spec_from_file_location('core998',HERE/'finite_receiver_formation.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'finite_receiver_formation_results.json');core.compare(core.run(),saved)
    assert saved['round']==998 and saved['all_scientific_checks_passed']
    assert saved['thermal_stability_only']
    for key in ('mechanical_stability_certified','unbound_to_bound_formation_simulated',
                'actual_heat_transfer_rate_calculated','actual_work_extracted',
                'record_device_constructed','cosmological_history_generated',
                'quantum_or_gravity_derived','full_goal_completed'):
        assert not saved[key],key
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    exact=lambda obj,key:F(obj[key]['exact'])
    ip=saved['inputs'];s=saved['stable_case']
    c=exact(ip,'c_gravity');A=exact(ip,'radiation_coefficient')
    g0=exact(ip,'initial_T_gravity');r0=exact(ip,'initial_T_radiation')
    T=exact(s,'final_T');Q=exact(s,'released_energy');E=exact(s,'total_energy')
    assert c==3 and A==F(3,4) and E==-c*g0+A*r0**4==-c*T+A*T**4
    assert Q==c*(T-g0)==A*(T**4-r0**4)==F(45,1024)
    assert exact(s,'kinetic_change')==Q and exact(s,'potential_change')==-2*Q
    assert exact(s,'source_energy_change')+exact(s,'receiver_energy_change')==0
    assert exact(s,'radius_ratio')==g0/T==F(497,512)
    Cr=4*A*T**3
    assert exact(s,'radiation_heat_response')==Cr==F(3,8)<c
    assert exact(s,'linear_drift_over_conductance')==1/c-1/Cr==F(-7,3)
    assert exact(s,'entropy_curvature_wrt_transferred_energy')==(1/c-1/Cr)/T**2==F(-28,3)
    # Independent high-precision logarithm checks the rational entropy enclosure.
    lo=exact(s,'total_entropy_change_lower');hi=exact(s,'total_entropy_change_upper')
    with localcontext() as ctx:
        ctx.prec=90
        dec=lambda f:Decimal(f.numerator)/Decimal(f.denominator)
        entropy=-dec(c)*(dec(T)/dec(g0)).ln()+dec(F(4,3)*A*(T**3-r0**3))
        assert Decimal(0)<dec(lo)<entropy<dec(hi)
    assert hi-lo<F(1,10**35)
    # Independent fixed-energy temperature relation and derivative on the whole interval.
    assert c==3 and A/c==F(1,4) and -E/c==F(31,64)
    assert r0**3-1<T**3-1<0
    u=saved['same_energy_unstable_case'];a=exact(u,'T_lower');b=exact(u,'T_upper')
    polynomial=lambda x:A*x**4-c*x-E
    assert 1<a<b<2 and polynomial(a)<=0<=polynomial(b)
    assert 4*A*a**3>c and exact(u,'linear_drift_over_conductance_lower')>0
    boundary=saved['domain_boundary'];tc=exact(boundary,'critical_T')
    assert 4*A*tc**3==c
    assert exact(boundary,'minimum_equilibrium_energy')==-F(3,4)*c*tc
    assert exact(boundary,'illustrative_no_equilibrium_energy')<exact(boundary,'minimum_equilibrium_energy')<E<0
    assert exact(boundary,'infinite_fixed_bath_linear_drift_over_conductance')==1/c>0
    curve=s['curve'];assert len(curve)==17
    for row in curve:
        q=F(row['Q_fraction'])*Q
        assert math.isclose(row['T_gravity'],float(g0+q/c),abs_tol=1e-15)
        assert abs(-float(c)*row['T_gravity']+float(A)*row['T_radiation']**4-float(E))<1e-14
    assert all(x['entropy_change']<y['entropy_change'] for x,y in zip(curve,curve[1:]))

    note=STAGE/'research_note_998.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','固定总能量','固定体积C_V','没有从初始等离子体制造它','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'formation_resource_adoption_v1.md',HERE/'finite_receiver_formation.py',
        HERE/'finite_receiver_formation_results.json',Path(__file__),HERE/'drafts/formation_selection.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish998.py',STAGE/'999/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至998）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=998)==list(range(1,999))
    assert '001—998轮共998份' in nav[0].read_text('utf-8-sig')
    assert '231—998的768份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '998：组织形成与有限接收者' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'997/research_round_997_checks.json',HERE/'drafts/STATUS.md',
              HERE/'drafts/boundary_adoption_audit_checks.json']
    prev=read(oldfiles[0])
    out=dict(round=998,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=998,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_997=prev['cumulative_numbered_test_groups_from_996']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,finite_energy_and_entropy_enclosure_verified=True,
        adopted_interface='classical conditional structural contraction with a changing finite internal receiver',
        thermal_stability_only=True,mechanical_stability_certified=False,
        actual_work_extracted=False,record_device_constructed=False,
        common_cosmological_history_generated=False,quantum_or_gravity_derived=False,
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
