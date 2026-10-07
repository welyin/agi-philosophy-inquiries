"""Reproduce 866 and verify frozen evidence, numbering and links."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import joint_loop_instrument as experiment
RECEIPT=HERE/'research_round_866_checks.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,866):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_866.md';body=note.read_text('utf-8')
    assert body.count('$$')==18 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,10)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'867/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=866)==list(range(1,867))
    frozen=[STAGE/'865/research_round_865_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (765,767,785,786,793,796,812,860,861,862,863,864,865)]
    frozen += [STAGE/'863/physical_relational_loop_bridge.py',STAGE/'863/physical_relational_loop_bridge_results.json']
    frozen += [STAGE.parent/'archive_554_584'/f'research_note_{n}.md' for n in (557,569)]
    frozen += [STAGE.parent/'archive_585_628'/f'research_note_{n}.md' for n in (598,623)]
    fresh=[note,STAGE/'867/drafts/STATUS.md',HERE/'joint_loop_instrument.py',experiment.TARGET,Path(__file__),HERE/'invariant_joint_loop_probe.py',HERE/'invariant_joint_loop_probe_results.json']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=866,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_865=3651,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_865_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=866,local_links_checked=links,display_equations=9,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='A strictly positive class-function kernel on the original quotient group has an exact finite character-moment cubature. It yields a 1368-label normal CP graph instrument preserving every Gauss branch, the selected joint means, and a finite source electric-energy form. The material classical symbols and physical first sources match; a continuum quantum instrument and autonomous original-action realization are not proved.',
        finite_outcome_labels=1368,exact_rational_Fourier_certificate=True,
        noise_partner_retained=True,electric_off_diagonal_sources_retained=True,
        classical_material_source_bridge=True,continuum_CP_lift_proved=False,
        full_graph_continuum_quantum_dynamics_matching_proved=False,
        finite_control_and_preparation_from_original_action_proved=False,
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
