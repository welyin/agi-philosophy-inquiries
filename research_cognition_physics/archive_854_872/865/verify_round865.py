"""Reproduce 865 and verify frozen evidence, numbering and links."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import relational_loop_fusion_matching as experiment
RECEIPT=HERE/'research_round_865_checks.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,865):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_865.md';body=note.read_text('utf-8')
    assert body.count('$$')==18 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,10)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'866/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=865)==list(range(1,866))
    frozen=[STAGE/'864/research_round_864_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (765,767,785,786,793,796,812,860,861,862,863,864)]
    frozen += [STAGE/'863/physical_relational_loop_bridge.py',STAGE/'863/physical_relational_loop_bridge_results.json']
    frozen += [STAGE.parent/'archive_554_584'/f'research_note_{n}.md' for n in (557,569)]
    fresh=[note,STAGE/'866/drafts/STATUS.md',HERE/'relational_loop_fusion_matching.py',experiment.TARGET,Path(__file__),HERE/'loop_self_contact_probe.py',HERE/'loop_self_contact_probe_results.json',HERE/'drafts/loop_insertion_renormalization_working.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=865,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_864=3650,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_864_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=865,local_links_checked=links,display_equations=9,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='A declared Casimir-perimeter normalization and an unchanged naive same-anchor product cannot both preserve the original quotient-group character fusion. The leading defect is nonzero on the existing constrained physical family. Positive amplification is not positive on the unchanged bounded graph configuration algebra. This is a matching constraint, not a computation of the full interacting loop counterterms.',
        same_anchor_composite_distinguished_from_independent_averages=True,
        full_reference_and_metric_sources_retained=True,
        actual_constraint_background_and_proper_length_used=True,
        declared_finite_candidate_not_computed_counterterm=True,
        all_order_interacting_quantum_loop_defined=False,
        finite_graph_state_and_dynamics_matching_proved=False,
        new_bulk_fields_or_couplings=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists(),'Do not overwrite frozen receipt.'
    result=verify(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for k in ('frozen_inputs','new_scientific_and_entry_files','argument_scope'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
