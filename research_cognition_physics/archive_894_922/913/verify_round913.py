"""913: independent original-record derivative and archived evidence checks."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.dont_write_bytecode=True;sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import discrete_record_tangent as science
TARGET=HERE/'research_round_913_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,913):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for key in ('frozen_inputs','new_scientific_and_entry_files'):add(old[key])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json',
        '887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json',
        '897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json',
        '904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json',
        '912/drafts/working_checks.json','913/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    audit=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(audit['preserved_files']);add(audit['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):
        add(read(STAGE/path)['evidence_and_audit_hashes'])
    for path,digest in frozen.items():assert sha(ROOT/path)==digest,path
    layout=Layout().verify()
    direct=read(science.TARGET);rebuild=read(science.rebuilt.TARGET);original=read(science.work.TARGET)
    for result in (direct,rebuild,original):
        for path,digest in result['source_hashes'].items():assert sha(STAGE/path)==digest,path
    assert [row['segments'] for row in direct['rows']]==[8,16,32]
    errors=[abs(row['direct_minus_curvature_source'][2]) for row in direct['rows']]
    ratios=[errors[i]/errors[i+1] for i in (0,1)]
    assert all(3.9<x<4.1 for x in ratios)
    finite=direct['rows'][0]['small_step_independent_rebuilds']
    fratio=abs(finite[0]['difference_from_discrete_tangent'][2]/finite[1]['difference_from_discrete_tangent'][2])
    assert abs(finite[1]['difference_from_discrete_tangent'][2])<abs(finite[0]['difference_from_discrete_tangent'][2])
    assert rebuild['rows'][0]['plus']['maximum_coordinate_shift']>.4
    assert rebuild['rows'][1]['plus']['maximum_coordinate_shift']>.2
    assert max(abs(x) for x in rebuild['rows'][1]['derivative_minus_source_pairing'])>1e-7
    for row in direct['rows']:
        assert row['geometry']['reference_Jacobian_two_jet_error']<2e-13
        assert row['geometry']['implicit_reference_tangent_error']<1e-14
        assert row['geometry']['maximum_material_coordinate_derivative']>19
    # Fresh original17^3 background/response and direct8-segment path derivative.
    bg,field,_=science.rebuilt.build_pair()
    src=science.loop.LoopSource(bg,2,8);src.kernels()
    actual,geom=science.exact_record_derivative(bg,field,src)
    expected=direct['rows'][0]['direct_discrete_tangent']
    assert np.allclose(actual,expected,rtol=1e-9,atol=1e-17)
    pairing=src.response(lambda x,p:field.jets(x))
    assert np.allclose(pairing['total'],original['rows'][0]['full_pairing']['total'],rtol=1e-9,atol=1e-17)
    note=STAGE/'research_note_913.md';body=note.read_text('utf-8')
    assert body.count('$$')==10 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(i) for i in range(1,6)]
    docs=[note,STAGE/'914/drafts/STATUS.md']+list((HERE/'drafts').glob('*.md'))
    docs += [STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for file in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',file.name)
        if m and file.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=913)==list(range(1,914))
    oldfiles=[STAGE/'912/research_round_912_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/working_checks.json']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (803,805,869,870,872,899,908,911,912)]
    newfiles=[note,Path(__file__),HERE/'independent_record_rebuild.py',science.rebuilt.TARGET,HERE/'discrete_record_tangent.py',science.TARGET,STAGE/'914/drafts/STATUS.md']
    for file in newfiles:
        if file.suffix=='.py':ast.parse(file.read_text('utf-8'))
    result=dict(round=913,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_912=3698,formal_reports=913,display_equations=5,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
        original17_grid_direct_record_derivative_reproduced=True,
        source_pairing_reproduced=True,path_discretization_difference_ratios=ratios,
        small_step_finite_difference_ratio_diagnostic=fratio,
        frozen_inputs={str(file.relative_to(ROOT)):sha(file) for file in oldfiles},
        new_scientific_and_entry_files={str(file.relative_to(ROOT)):sha(file) for file in newfiles},
        physical_A_error_certified=False,actual_finite_observable_error_certified=False,
        full_vA_field_computed=False,full872_response_computed=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False)
    if not writing:
        old=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert old[k]==result[k],k
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run(args.write)
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
