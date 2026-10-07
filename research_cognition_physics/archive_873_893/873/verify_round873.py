"""Reproduce873 and audit its source matching scope and frozen history."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import clock_mass_source_bridge as experiment
RECEIPT=HERE/'research_round_873_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    previous=json.loads((STAGE/'872/read_vertex_source_closure_results.json').read_text('utf-8'))
    assert science['cumulative_numbered_groups']==previous['cumulative_numbered_groups']+1==3658
    assert science['fresh_numbered_groups']==1
    for n in range(776,873):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_873.md';body=note.read_text('utf-8')
    assert body.count('$$')==18 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,10)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'874/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += list((HERE/'drafts').glob('*.md'))
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$', '', doc.read_text('utf-8-sig'), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=873)==list(range(1,874))
    frozen=[STAGE/'872/research_round_872_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (802,803,853,861,869,871,872)]
    frozen += [STAGE/'861/magnetic_reduced_hamiltonian_results.json',
        STAGE.parent/'archive_629_652/research_note_643.md',
        STAGE/'854/drafts/共同模型阶段报告_截至853.md']
    fresh=[note,STAGE/'874/drafts/STATUS.md',HERE/'clock_mass_source_bridge.py',
        experiment.TARGET,Path(__file__),HERE/'drafts/clock_source_interface_working.md',
        HERE/'drafts/共同模型阶段报告_截至872.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=873,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_872=science['cumulative_numbered_groups'],
        saved_result_reproduced=True,historical_manifest_evidence=history,
        previous_776_through_872_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=873,local_links_checked=links,display_equations=9,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='In the original two-derivative parent sector and regular material-clock chart, two existing neutral receiver mass sources induce the nonzero mixed four-leg physical Hamiltonian coefficient 2alpha U1 U2/v_h^3. The full implicit derivative retains mass-independent fermion constraint terms. It excludes direct source-compatible identification with an entirely quadratic conditional CAR family; it does not refute old fixed-history determinants or prove complete reduced quantization.',
        original_material_mass_source_map=True,
        full_implicit_mass_derivative_at_four_legs=True,
        exact_Grassmann_original_constraint_check=True,
        finite_CAR_calibration_is_full_quantization=False,
        old_conditional_determinants_refuted=False,
        original_graph_continuum_bridge_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope','cumulative_numbered_test_groups_from_872'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
