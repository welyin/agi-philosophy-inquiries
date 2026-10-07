"""Reproduce 855, check frozen evidence, publication links and stated scope."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import locality_preserving_task_compression as experiment
RECEIPT=HERE/'research_round_855_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,855):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_855.md';body=note.read_text('utf-8')
    assert body.count('$$')==14 and '	' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,8)]
    docs=[note,STAGE/'856/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += list((HERE/'drafts').glob('*.md'))
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
    assert sorted(n for n in numbers if n<=855)==list(range(1,856))
    frozen=[STAGE/'854/research_round_854_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (796,801,844,853,854)]
    frozen += [STAGE.parent/'archive_702_741'/f'research_note_{n}.md' for n in (705,706)]
    frozen += [STAGE.parent/'archive_742_763'/f'research_note_{n}.md' for n in (756,761,763)]
    fresh=[note,STAGE/'856/drafts/STATUS.md',HERE/'locality_preserving_task_compression.py',experiment.TARGET,Path(__file__)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=855,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_854=3640,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_854_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=855,local_links_checked=links,display_equations=7,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='An explicit total-word compression satisfying the 854 record-algebra and finite-jet conditions breaks exact independence of two originally commuting local processes at a cutoff boundary. A tensor product closure preserves the same finite jets and local supports under explicit finite factorization and initial tensor-rank hypotheses. These hypotheses and Gauss compatibility have not been established as a common original graph-continuum dictionary.',
        exact_mixed_record_probability_defect_proved=True,
        finite_order_matching_not_contradicted=True,
        conditional_tensor_repair_preserves_declared_local_supports=True,
        initial_tensor_rank_and_coefficient_factorization_are_extra_conditions=True,
        original_graph_or_all_compressions_disproved=False,
        original_Gauss_and_continuum_tensor_dictionary_completed=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
