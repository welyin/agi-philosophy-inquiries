"""Verify dynamic thermoelastic interface 1000; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_1000_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,1000):
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
    import numpy as np
    spec=importlib.util.spec_from_file_location('core1000',HERE/'thermoelastic_common_account.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'thermoelastic_common_account_results.json');core.compare(core.run(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1000
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('long_time_ode_simulated','microscopic_phase_derived','transport_coefficients_derived',
                'noise_or_memory_lifetime_certified','all_discrete_nonuniform_modes_damped',
                'einstein_geometry_derived','quantum_channel_constructed','full_goal_completed'):
        assert not saved[key],key
    assert saved['discrete_checkerboard_null_mode_preserved']
    p={key:F(value['exact']) for key,value in saved['parameters'].items()}
    rho,K,b,c,T0,kappa=[p[key] for key in ('rho','K','b','c','T0','kappa')]
    state={key:[F(x['exact']) for x in value] for key,value in saved['finite_state'].items()}
    eps,v,T=state['strain'],state['velocity'],state['temperature']
    assert len(T)==4 and all(t>0 for t in T)
    # Independently differentiate E, S and A with the delivered actual rates.
    def exact(row,key):return F(row[key]['exact'])
    for key,expected in (('matched',F(0)),('omitted_thermal_feedback',F(1,25))):
        row=saved[key]
        de,dv,dt=[[F(x['exact']) for x in row[n]] for n in
                  ('strain_rate','velocity_rate','temperature_rate')]
        dE=sum(rho*vel*acc+(K*ep+b*T0)*dx+c*dtemp for ep,vel,acc,dx,dtemp in zip(eps,v,dv,de,dt))
        dS=sum(b*dx+c*dtemp/t for dx,dtemp,t in zip(de,dt,T))
        dA=sum(rho*vel*acc+K*ep*dx+c*(1-T0/t)*dtemp for ep,vel,acc,dx,dtemp,t in zip(eps,v,dv,de,dt,T))
        assert dE==expected==exact(row,'energy_rate')
        assert dS==F(1,12)==exact(row,'entropy_rate')
        assert dA==dE-T0*dS==exact(row,'availability_rate')
    edge=kappa*sum((T[i]-T[(i+1)%4])**2/(T[i]*T[(i+1)%4]) for i in range(4))
    assert edge==F(1,12)==F(saved['edge_entropy_production']['exact'])
    mode=saved['continuous_linear_symbol'];k=F(mode['k']['exact'])
    expected=[rho*c,rho*kappa*k**2,(K*c+b*b*T0)*k**2,K*kappa*k**4]
    assert [F(x['exact']) for x in mode['polynomial_coefficients']]==expected
    assert expected[1]*expected[2]-expected[0]*expected[3]==F(1,2)==F(mode['routh_margin']['exact'])
    for z in mode['roots']:
        root=complex(z['real'],z['imag'])
        assert root.real<0 and abs(np.polyval([float(x) for x in expected],root))<1e-12
    assert F(mode['isothermal_speed_squared']['exact'])==4
    assert F(mode['adiabatic_speed_squared']['exact'])==F(9,2)
    note=STAGE/'research_note_1000.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','累计3783','棋盘格','未构造全部微观量子通道','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'material_phase_transport_adoption_v1.md',HERE/'thermoelastic_common_account.py',
        HERE/'thermoelastic_common_account_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1000.py',STAGE/'1001/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1000）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1000)==list(range(1,1001))
    assert '001—1000轮共1000份' in nav[0].read_text('utf-8-sig')
    assert '231—1000的770份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1000：宏观相与热—力输运' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'999/research_round_999_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=1000,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1000,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_999=prev['cumulative_numbered_test_groups_from_998']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,energy_entropy_availability_separately_verified=True,
        positive_entropy_alone_insufficient_in_stated_countermodel=True,
        microscopic_phase_derived=False,joint_record_lifecycle_certified=False,
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
