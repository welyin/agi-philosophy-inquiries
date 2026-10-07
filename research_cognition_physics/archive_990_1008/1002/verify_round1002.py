"""Verify thermal record reuse interface 1002; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_1002_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,1002):
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
    spec=importlib.util.spec_from_file_location('reuse1002',HERE/'thermal_record_reuse.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'thermal_record_reuse_results.json');core.compare(core.calculate(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1002
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('original_charge_interaction_implements_writer','all_resources_prepared_from_thermal_states',
        'population_effect_is_qnd_under_full_active_H','infinite_retention_certified',
        'cyclic_auxiliary_restoration_certified','physical_readout_instrument_constructed',
        'joint_clock_communication_lifecycle_certified','full_goal_completed'):assert not saved[key],key
    old=read(STAGE/'980/finite_thermal_records_results.json')
    p=np.array(saved['parameters']['thermal_populations']);e=np.array(saved['parameters']['energies'])
    assert saved['parameters']['thermal_populations']==old['parameters']['thermal_populations']
    assert saved['parameters']['energies']==old['parameters']['energies']
    q=float(p[2]+p[3]);C=float(p[0]+p[1]-p[2]-p[3]);g=saved['parameters']['g']
    assert abs(q-saved['peak_equal_prior_error'])<1e-14
    assert abs(C-saved['population_contrast'])<1e-14
    for row in saved['samples']:
        err=.5*(1-C*math.sin(g*row['time'])**2)
        assert abs(err-row['equal_prior_error'])<2e-12
    intervals=old['eigensystem_certificate']['coarse_analytic_certificate']['thermal_intervals']
    q_hi=sum(F(str(intervals[i][1])) for i in (2,3))
    assert q_hi<F(188903,1000000)
    bound=F(99,100)*F(188903,1000000)+F(1,200)+F(13,14000)+F(1,10000)
    assert bound==F(saved['certified_after_reset_error_bound']['exact'])<F(193043,1000000)
    aux=np.array(saved['parameters']['auxiliary_energies'])
    # Enumerate endpoint joint probabilities, independently of the H exponential.
    joint=np.zeros((4,5))
    for i in range(4):joint[3-i,i+1]=p[i]
    marginal_m=joint.sum(axis=1);marginal_b=joint.sum(axis=0)
    dM=float(e@(marginal_m-p));dB=float(aux@marginal_b-aux[0])
    assert abs(dM-saved['conditional_1_memory_energy_gain'])<1e-13
    assert abs(dB-saved['conditional_1_auxiliary_energy_change'])<1e-12
    assert abs(dM+dB)<1e-12
    def ent(arr):
        arr=np.asarray(arr);arr=arr[arr>0]
        return float(-np.dot(arr,np.log(arr)))
    mutual=ent(marginal_m)+ent(marginal_b)-ent(joint.ravel())
    assert abs(mutual-saved['conditional_1_memory_auxiliary_mutual_information'])<1e-12
    assert marginal_b[0]==0 and math.isclose(sum(marginal_b[1:]),1,abs_tol=1e-15)
    assert saved['unknown_source_dephased_locally']
    assert saved['auxiliary_initial_pure_energy_state_is_input']
    assert saved['engineered_energy_matched_interaction_is_input']
    note=STAGE/'research_note_1002.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for term in ('整体目标未完成','累计3785','不是持续全过程的QND读取','完整记录了来源标签','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'organization_lifecycle_adoption_v1.md',HERE/'thermal_record_reuse.py',
        HERE/'thermal_record_reuse_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1002.py',STAGE/'1003/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1002）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1002)==list(range(1,1003))
    assert '001—1002轮共1002份' in nav[0].read_text('utf-8-sig')
    assert '231—1002的772份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1002：热材料再写与功能分工' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1001/research_round_1001_checks.json',HERE/'drafts/STATUS.md',
              STAGE/'980/finite_thermal_records_results.json']
    prev=read(oldfiles[0])
    out=dict(round=1002,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1002,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_1001=prev['cumulative_numbered_test_groups_from_1000']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,old_thermal_state_reused=True,
        finite_window_rational_bound_independently_verified=True,
        auxiliary_energy_entropy_and_record_all_accounted=True,
        joint_clock_communication_lifecycle_certified=False,
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
