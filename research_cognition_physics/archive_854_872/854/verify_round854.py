"""Reproduce 854, check frozen evidence, publication links and stated scope."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import finite_task_jet_bridge as experiment
RECEIPT=HERE/'research_round_854_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,854):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_854.md';body=note.read_text('utf-8')
    assert body.count('$$')==20 and '	' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,11)]
    docs=[note,STAGE/'855/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=854)==list(range(1,855))
    frozen=[STAGE/'853/research_round_853_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/finite_task_bridge_working.md',HERE/'drafts/共同模型阶段报告_截至853.md',HERE/'drafts/joint_model_audit_853_receipt.json']
    frozen += [STAGE/f'research_note_{n}.md' for n in (783,793,795,796,798,801,844,851,852,853)]
    fresh=[note,STAGE/'855/drafts/STATUS.md',HERE/'finite_task_jet_bridge.py',experiment.TARGET,Path(__file__)]
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=854,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_853=3639,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_853_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=854,local_links_checked=links,display_equations=10,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='A fixed finite word space in the inherited physical coefficient representation yields an exactly unitary finite task model matching the original formal process, finite records and bilateral source jets through prescribed order, for all inputs in a fixed finite preparation factor. The entire environment memory is retained. Application to 853 uses its one total actual S and finite final CAR outputs; multiple stages require already compatible physical stages. This does not identify the original graph Hamiltonian, supply a finite-coupling continuum remainder, or prove locality and full gauge algebra preservation.',
        input_and_source_independent_projection=True,
        all_input_amplitudes_and_bilateral_source_jets_matched=True,
        normalized_source_order_budget_N_plus_K=True,
        new_finite_model_exact_unitarity_not_continuum_convergence=True,
        old_graph_Hamiltonian_or_locality_identified=False,
        full_continuum_Gauss_algebra_or_CCR_represented=False,
        uniform_resource_bound_or_all_orders_convergence=False,
        spacetime_group_action_remain_inputs=True,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite a frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
