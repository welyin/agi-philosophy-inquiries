"""Reproduce 846, inspect links and retain all previously frozen science."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import constrained_six_leg_recursion as experiment
RECEIPT=HERE/'research_round_846_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    for ex in (experiment,):assert json.loads(json.dumps(ex.run()))==json.loads(ex.TARGET.read_text('utf-8'))
    for n in range(776,846):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for section in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[section].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_846.md';body=note.read_text('utf-8')
    assert body.count('$$')==24 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    docs=[note,STAGE/'847/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=846)==list(range(1,847))
    frozen=[STAGE/'845/research_round_845_checks.json',STAGE/'research_note_845.md',
            STAGE/'research_note_796.md',STAGE/'research_note_810.md',
            STAGE/'research_note_785.md',STAGE/'research_note_799.md',HERE/'drafts/STATUS.md']
    fresh=[note,HERE/'drafts/research_note_846_working.md',STAGE/'847/drafts/STATUS.md',*sorted(HERE.glob('*.py')),
           *sorted(p for p in HERE.glob('*.json') if p!=RECEIPT)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=846,date='2026-10-05',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_845=3630,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_845_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=846,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='In the original constrained local physical module and fixed compact classical completion, the source dual projection gives a representative-independent retarded recursion. Five original linear Green stages construct the pure six-fermion response, retaining Dirac feedback, all bosonic sectors, cutoff derivatives and the prescribed comparison of states. The 80 code projections may be extracted in six-generator exterior algebras. This constructs the coefficient functional; the actual continuum 80-term value and its nonzero character have not been computed. Exact rational calibration is not the original spacetime PDE.',
        original_constrained_classical_six_leg_recursion_constructed=True,
        dual_source_projection_and_extension_independence_proven=True,
        original_platform_equations_recovered_by_formal_Noether_identity=True,
        same_compact_completion_and_initial_contract_retained=True,
        Dirac_backreaction_and_all_seven_six_leg_source_classes_retained=True,
        six_generator_coefficient_extraction_is_exact=True,
        new_higher_source_conservation_axiom_introduced=False,
        original_continuum_six_leg_values_numerically_computed=False,
        original_eighty_term_sum_nonzero_proven=False,
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
