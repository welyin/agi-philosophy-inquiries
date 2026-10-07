"""Reproduce 867 and verify frozen evidence, numbering and links."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import shared_weyl_readout as experiment
RECEIPT=HERE/'research_round_867_checks.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(writing=False):
    assert experiment.run()==json.loads(experiment.TARGET.read_text('utf-8'))
    for n in range(776,867):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    history=Layout().verify()
    note=STAGE/'research_note_867.md';body=note.read_text('utf-8')
    assert body.count('$$')==24 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(n) for n in range(1,13)]
    for formula in re.findall(r'\$\$(.*?)\$\$',body,re.S):
        assert not re.search(r'(?<!\\)\b(quad|qquad|left|right)\b',formula),formula
    docs=[note,STAGE/'868/drafts/STATUS.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
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
    assert sorted(n for n in numbers if n<=867)==list(range(1,868))
    frozen=[STAGE/'866/research_round_866_checks.json',HERE/'drafts/STATUS.md']
    frozen += [STAGE/f'research_note_{n}.md' for n in (765,767,785,786,793,796,797,809,810,853,860,861,862,863,864,865,866)]
    frozen += [STAGE/'866/joint_loop_instrument.py',STAGE/'866/joint_loop_instrument_results.json',STAGE/'866/invariant_joint_loop_probe.py']
    fresh=[note,STAGE/'868/drafts/STATUS.md',HERE/'shared_weyl_readout.py',experiment.TARGET,Path(__file__),HERE/'drafts/shared_weyl_readout_working.md']
    for p in fresh:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    return dict(round=867,date='2026-10-06',all_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_866=3652,saved_result_reproduced=True,
        historical_manifest_evidence=history,previous_776_through_866_frozen_hashes_verified=True,
        historical_science_rerun=False,formal_reports=867,local_links_checked=links,display_equations=12,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in frozen},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in fresh},
        argument_scope='The original1368-label menu has a common bounded CP lift matching baseline probabilities and every complete first material source in the same physical Hadamard representation. Summable Weyl series preserve the Wick domain and fixed formal S transport. The leading nonselective backreaction has a nonzero original physical witness and a full relative-source Ward identity. Local probe realization, strong causality, nonlinear graph matching and finite bulk coupling remain open.',
        original_menu_labels=1368,all_menu_first_sources_matched=True,
        finite_real_free_readout_strength=True,
        weighted_Weyl_domain_and_Hadamard_smooth_difference='analytic',
        second_order_same_source_backreaction=True,
        original_S_transport='formal bulk coupling',
        quantum_spacetime_propagator_computed_numerically=False,
        matrix_probe_is_not_continuum_CCR=True,
        original_local_probe_implementation_proved=False,
        complete_graph_continuum_dynamics_matching_proved=False,
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
