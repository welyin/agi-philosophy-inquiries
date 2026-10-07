"""Verify working875 without promoting it to a completed research round."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import nonstationary_clock_current_probe as simple
import clock_boundary_record_transport as boundary
RECEIPT=HERE/'drafts/clock_transport_working_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    for experiment in (simple,boundary):
        r=experiment.run()
        assert r==json.loads(experiment.TARGET.read_text('utf-8'))
        assert r['formal_reports']==874 and r['fresh_numbered_groups']==0
        assert r['cumulative_numbered_groups']==3659
    for n in range(776,875):
        old=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in old[key].items():assert sha(ROOT/name)==digest,name
    doc=HERE/'drafts/nonstationary_clock_transport_working.md'
    body=doc.read_text('utf-8')
    assert body.count('$$')==22 and '\t' not in body
    assert re.findall(r'\\tag\{W(\d+)\}',body)==[str(i) for i in range(1,12)]
    docs=[doc,STAGE.parent/'RESEARCH_STATE.md',STAGE/'文件索引.md',
          STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    links=0
    for d in docs:
        text=re.sub(r'\$\$.*?\$\$','',d.read_text('utf-8-sig'),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            p=(d.parent/target.split('#')[0].strip('<>')).resolve()
            assert p.exists() or (writing and p==RECEIPT.resolve()),(d,target)
            links+=1
    files=[doc,Path(__file__),HERE/'nonstationary_clock_current_probe.py',simple.TARGET,
           HERE/'clock_boundary_record_transport.py',boundary.TARGET]
    for f in files:
        if f.suffix=='.py':ast.parse(f.read_text('utf-8'),filename=str(f))
    return dict(date='2026-10-06',working_round=875,formal_reports=874,
        cumulative_numbered_groups=3659,fresh_numbered_groups=0,all_checks_passed=True,
        two_saved_results_reproduced=True,previous_776_through_874_frozen_hashes_verified=True,
        historical_manifest_evidence=Layout().verify(),local_links_checked=links,display_equations=11,
        files={str(p.relative_to(ROOT)):sha(p) for p in files},
        current_and_constant_mass_endpoint_and_ordered_record_transport_checked=True,
        original_field_quantum_ordering_or_state_mapping_proved=False,
        new_second_order_clock_ordering_required_by_original_goal=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not RECEIPT.exists()
    result=run(args.write)
    if args.write:RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(RECEIPT.read_text('utf-8'))
        for key in ('files','formal_reports','cumulative_numbered_groups','fresh_numbered_groups'):
            assert old[key]==result[key],key
    print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
