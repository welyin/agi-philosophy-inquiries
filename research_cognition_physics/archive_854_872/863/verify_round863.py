"""Reproduce 863, verify frozen history, publication scope and references."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import physical_relational_loop_bridge as experiment
RECEIPT=HERE/'research_round_863_checks.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    result=experiment.run()
    assert result==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,863):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_863.md';body=note.read_text('utf-8')
    assert body.count('$$')==24 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    docs=[note,STAGE/'864/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=863)==list(range(1,864))
    frozen=[STAGE/'862/research_round_862_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (765,767,785,796,859,860,861,862)]
    frozen += [STAGE.parent/'archive_554_584'/f'research_note_{n}.md' for n in (569,573)]
    frozen += [STAGE.parent/'archive_742_763'/'research_note_753.md',STAGE.parent/'archive_742_763'/'753/joint_reference_constraint_strata.py',STAGE/'859/shared_probe_reference_probe.py',STAGE/'862/material_link_path_probe.py']
    fresh=[note,STAGE/'864/drafts/STATUS.md',HERE/'physical_relational_loop_bridge.py',experiment.TARGET,Path(__file__),HERE/'smeared_loop_current_probe.py',HERE/'smeared_loop_current_probe_results.json',HERE/'drafts/smeared_relational_loop_working.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=863,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_862=3648,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_862_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=863,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='Finite-resolution relational loop linearization is a smooth compact full mixed conserved test on the original background. An exact constrained color family with unchanged surface references gives a nonzero physical pairing; an explicit conserved companion test proves finite strictly positive variance in the existing physical Hadamard state.',
        original_global_quotient_representation_preserved=True,
        original_constraints_preserved_on_witness_family='analytic constant electric/magnetic energy with zero Gauss and momentum',
        full_mixed_current_nontrivial_on_original_solutions_proved=True,
        numerical_full_quantum_variance_evaluated=False,
        all_order_interacting_loop_defined=False,finite_graph_quantum_matching_proved=False,
        autonomous_internal_preparation_or_record_proved=False,
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
