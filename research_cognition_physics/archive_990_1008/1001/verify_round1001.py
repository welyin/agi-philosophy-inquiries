"""Verify nuclear formation interface 1001; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_1001_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,1001):
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
    spec=importlib.util.spec_from_file_location('nuclear1001',HERE/'nuclear_formation_screen.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'nuclear_formation_results.json');core.compare(core.calculate(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1001
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('capture_rate_derived','cosmic_history_simulated','physical_abundances_predicted',
        'mass_uncertainty_propagated','work_extraction_certified','joint_record_lifecycle_certified',
        'cognitive_derivation_of_nuclear_force','full_goal_completed'):assert not saved[key],key
    masses={key:F(value['exact']) for key,value in saved['mass_inputs_MeV'].items()}
    mp,mn,md=[masses[key] for key in ('proton','neutron','deuteron')]
    W=mp+mn;B=W-md
    eg=F(saved['capture_photon_MeV']['exact']);kd=F(saved['deuteron_recoil_MeV']['exact'])
    threshold=F(saved['photodissociation_threshold_MeV']['exact'])
    assert eg+kd==B and (W-eg)**2-eg**2==md*md
    assert md*md+2*md*threshold==W*W
    assert F(saved['matched_energy_residual_MeV']['exact'])==0
    assert F(saved['double_count_residual_MeV']['exact'])==B
    assert math.isclose(math.hypot(float(md),float(B))+float(B)-float(W),
        saved['photon_equals_binding_on_shell_energy_excess_MeV'],rel_tol=2e-9)
    # Independently reconstruct densities and solve the equilibrium condition
    # by monotone bisection, not the generating program's quadratic formula.
    inp=saved['thermal_inputs'];a=inp['neutron_constituent_fraction'];b=inp['proton_constituent_fraction']
    for row in saved['equilibrium_rows']:
        T=row['temperature_MeV'];y=row['Yd']
        nb=inp['eta']*2*inp['zeta3']/math.pi**2*T**3
        K=.75*(2*math.pi/T)**1.5*float(md/(mp*mn))**1.5*math.exp(float(B)/T)
        assert math.isclose(nb*K,row['R'],rel_tol=2e-14)
        lo,hi=0.0,min(a,b)
        for _ in range(180):
            mid=(lo+hi)/2
            if mid-nb*K*(a-mid)*(b-mid)>0:hi=mid
            else:lo=mid
        assert math.isclose((lo+hi)/2,y,rel_tol=5e-13)
        assert math.isclose(row['Yn']+row['Yp']+2*y,1,abs_tol=2e-15)
        assert math.isclose(row['Yp']+y,b,abs_tol=2e-15)
        assert math.isclose(row['Yn']+y,a,abs_tol=2e-15)
        assert row['max_thermal_fugacity']<1e-8
    note=STAGE/'research_note_1001.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for term in ('整体目标未完成','累计3784','四份受限平衡状态','直接复用971','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'nuclear_material_resource_adoption_v1.md',HERE/'nuclear_formation_screen.py',
        HERE/'nuclear_formation_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1001.py',STAGE/'1002/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1001）',1)[1].split('### 当前取舍',1)[0]
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
    for pth in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',pth.name)
        if m and pth.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=1001)==list(range(1,1002))
    assert '001—1001轮共1001份' in nav[0].read_text('utf-8-sig')
    assert '231—1001的771份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1001：核组成与材料资源' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1000/research_round_1000_checks.json',HERE/'drafts/STATUS.md',STAGE/'research_note_971.md']
    prev=read(oldfiles[0])
    out=dict(round=1001,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1001,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_1000=prev['cumulative_numbered_test_groups_from_999']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,old_recoil_interface_reused=True,
        restricted_equilibrium_independently_solved=True,
        full_nuclear_rates_or_yields_verified=False,joint_record_lifecycle_certified=False,
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
