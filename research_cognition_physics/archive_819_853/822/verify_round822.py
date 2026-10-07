"""Reproduce 822, inspect links and retain all previously frozen science."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import reference_preserving_compensation as experiment
RECEIPT=HERE/'research_round_822_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    for ex in (experiment,):assert ex.run()==json.loads(ex.TARGET.read_text('utf-8'))
    for n in range(776,822):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_822.md';body=note.read_text('utf-8')
    assert body.count('$$')==20 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,11)]
    docs=[note,STAGE/'823/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text('utf-8-sig')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            path=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert path.exists() or (writing and path==RECEIPT.resolve()),(str(doc),target)
            links+=1
    numbers=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if p.parent.name.startswith('archive_') and m:numbers.append(int(m[1]))
    assert sorted(n for n in numbers if n<=822)==list(range(1,823))
    frozen=[STAGE/'821/research_round_821_checks.json',STAGE/'research_note_821.md',
            STAGE/'research_note_796.md',STAGE/'research_note_810.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'823/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
           *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=822,date='2026-10-05',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_821=3600,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_821_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=822,local_links_checked=links,display_equations=10,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='Three Gauss-free weak-electric counterflows from the original nonconstant magnetic invariants have an analytically positive diagonal total-momentum Gram. They replace the radial scalar counterflows in the fixed-geometry all-source construction. The EW Gauss particular remains tangent to the Higgs orbit, so the Higgs modulus and singlet normal velocities, their full initial first jets, and the four original scalar reference values and spatial jets remain fixed. The old special-point reference determinant is preserved. This does not fix all reference normal derivatives, the full relational metric or the actual absolute quantum mean.',
        original_gauge_only_three_direction_counterflows_proven=True,
        reference_preserving_full_linear_material_right_inverse_proven=True,
        h_s_initial_first_jets_and_four_reference_spatial_jets_preserved=True,
        original_special_point_reference_rank_preserved=True,
        full_relational_metric_or_future_geometry_fixed=False,
        actual_quantum_absolute_mean_matched=False,
        finite_strength_or_autonomous_preparation_proven=False,
        graph_to_continuum_map_proven=False,spacetime_group_action_remain_inputs=True,
        original_all_physics_unification_complete=False,visual_checks_performed=False,app_goal_changed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    r=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==r[k],k
    print(json.dumps({k:v for k,v in r.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
