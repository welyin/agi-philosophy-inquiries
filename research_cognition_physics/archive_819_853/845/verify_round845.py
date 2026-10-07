"""Reproduce 845, inspect links and retain all previously frozen science."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import record_source_hbar_filtration as experiment
RECEIPT=HERE/'research_round_845_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    for ex in (experiment,):assert json.loads(json.dumps(ex.run()))==json.loads(ex.TARGET.read_text('utf-8'))
    for n in range(776,845):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_845.md';body=note.read_text('utf-8')
    assert body.count('$$')==24 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    docs=[note,STAGE/'846/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=845)==list(range(1,846))
    frozen=[STAGE/'844/research_round_844_checks.json',STAGE/'research_note_844.md',
            STAGE/'research_note_796.md',STAGE/'research_note_810.md',
            STAGE/'research_note_785.md',STAGE/'research_note_799.md',HERE/'drafts/STATUS.md']
    fresh=[note,STAGE/'846/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
           *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=845,date='2026-10-05',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_844=3629,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_844_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=845,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='For the original finite-menu perturbative outputs with no negative hbar normalization and fixed common background, code, complement weight and remaining preparation, different logical contents agree through hbar squared. The first allowed difference at hbar cubed is the classical six-fermion output kernel paired with the smooth sixth-moment difference. Its original-code projection has 80 signed directions for the selected logical axis. This locates the coefficient but does not prove that the original continuum response is nonzero. The separate scalar feedback calculation is an auxiliary diagnostic, not the original gravitational action.',
        original_unnormalized_source_content_difference_not_before_hbar_cubed=True,
        leading_difference_reduced_to_classical_six_fermion_kernel=True,
        fixed_background_leading_difference_independent_of_common_thermal_rest=True,
        original_code_sixth_degree_selected_axis_support=80,
        normalized_record_observables_excluded_from_hbar_bound=True,
        original_actual_six_leg_response_nonzero_proven=False,
        two_order_self_consistent_background_solution_constructed=False,
        hbar_order_identified_with_real_time_or_coupling_order=False,
        auxiliary_scalar_equation_is_original_action=False,
        finite_coupling_convergence_proven=False,
        autonomous_internal_preparation_proven=False,
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
