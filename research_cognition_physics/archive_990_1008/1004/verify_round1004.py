"""Verify finite strong-field response 1004; delivery checks are not physical unification."""
from pathlib import Path
import argparse, ast, hashlib, json, math, re, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_1004_checks.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,1004):
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
    extra['1004/drafts/resource_adoption_checks.json']=('frozen_evidence_hashes','audit_file_hashes')
    for path,keys in extra.items():
        value=read(STAGE/path)
        for key in keys:add(value[key])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()







    import importlib.util
    import numpy as np
    spec=importlib.util.spec_from_file_location('horizon1004',HERE/'finite_horizon_response.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'finite_horizon_response_results.json');core.compare(core.calculate(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1004
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('rigorous_quadrature_error_certified','numerical_calibration_is_four_dimensional',
        'finite_window_certifies_finite_HH_reservoir','exact_finite_window_KMS_ratio_claimed',
        'same_background_is_two_self_consistent_backreaction_solutions',
        'autonomous_physical_detector_constructed','hawking_spectrum_derived_from_cognitive_principles',
        'full_goal_completed'):assert not saved[key],key
    p=saved['parameters'];L=p['dimensionless_half_width'];Q=p['dimensionless_frequency_cut']
    assert abs(p['alpha']/p['detector_proper_acceleration']-16)<1e-13
    assert abs(p['alpha']*p['support_half_width']-8)<1e-14
    b=2*math.pi
    tail=4*L*L/math.pi*math.exp(-b*Q)/(1-math.exp(-b*Q))*(Q/b+1/b**2)
    assert math.isclose(tail,saved['spectral_H_minus_B_tail_upper_bound'],rel_tol=1e-14)
    assert 0<tail<6.80e-26
    assert math.isclose(saved['smooth_difference_at_zero'],p['alpha']**2/(24*math.pi),rel_tol=1e-14)
    # A separate scalar power-series evaluation checks the coincidence branch.
    for z in (0.,.001,.01):
        series=(1/12-z*z/240+z**4/6048-z**6/172800)/(2*math.pi)
        assert abs(float(core.smooth_kernel(np.array(z)))-series)<1e-17
    for row in saved['rows']:
        latest=row['evaluations'][-1]
        assert latest['absolute_algorithm_difference']<2.3e-13
        assert row['refinement_difference']<3.5e-13
        assert abs(2*row['H_minus_U']-latest['time_H_minus_B'])<1e-14
        assert row['H_minus_U']>0
    note=STAGE/'research_note_1004.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','累计3786','不是一般R处的固有加速度','测试场比较','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'strong_field_adoption_v1.md',HERE/'finite_horizon_response.py',
        HERE/'finite_horizon_response_results.json',Path(__file__),HERE/'drafts/horizon_adoption_decision.md',
        HERE/'drafts/horizon_review_notes.md',HERE/'drafts/publish1004.py',STAGE/'1005/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1004）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1004)==list(range(1,1005))
    assert '001—1004轮共1004份' in nav[0].read_text('utf-8-sig')
    assert '231—1004的774份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1004：几何、场态与有限强场记录' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1003/research_round_1003_checks.json',HERE/'drafts/STATUS.md',
              HERE/'drafts/resource_adoption_checks.json',HERE/'drafts/NEXT.md']
    prev=read(oldfiles[0])
    out=dict(round=1004,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1004,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_1003=prev['cumulative_numbered_test_groups_from_1002']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,one_plus_one_algorithms_agree=True,
        analytic_frequency_tail_bound_verified=True,
        four_dimensional_radial_spectrum_computed=False,
        physical_backreaction_certified=False,full_goal_completed=False,
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
