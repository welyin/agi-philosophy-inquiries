"""Verify working878 without adding a formal round."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
import auxiliary_heat_split_probe as experiment
TARGET=HERE/'drafts/working_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    science=experiment.run()
    assert science==json.loads(experiment.TARGET.read_text('utf-8'))
    assert science['formal_reports']==877 and science['fresh_numbered_groups']==0
    for n in range(776,878):
        receipt=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in receipt[key].items():assert sha(ROOT/name)==digest,name
    note=HERE/'drafts/auxiliary_heat_scope_working.md';body=note.read_text('utf-8')
    assert body.count('$$')==12 and re.findall(r'\\tag\{(W\d+)\}',body)==[f'W{i}' for i in range(1,7)]
    docs=[note,HERE/'drafts/STATUS.md',STAGE.parent/'RESEARCH_STATE.md',
          STAGE/'文件索引.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            p=(doc.parent/target.split('#')[0].strip('<>')).resolve()
            assert p.exists() or (writing and p==TARGET.resolve()),(doc,target)
            links+=1
    files=[note,HERE/'drafts/STATUS.md',Path(__file__),HERE/'auxiliary_heat_split_probe.py',experiment.TARGET]
    for p in files:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    return dict(date='2026-10-06',working_round=878,formal_reports=877,
        cumulative_numbered_groups=3662,fresh_numbered_groups=0,
        all_checks_passed=True,saved_result_reproduced=True,local_links_checked=links,
        historical_776_through_877_frozen_hashes_verified=True,
        files={str(p.relative_to(ROOT)):sha(p) for p in files},
        argument_scope='Only an abstract counterexample to inferring regional split from arbitrary auxiliary heat nuclearity; not a counterexample to the original E theory or cognition principles.',
        full_goal_completed=False,app_goal_changed=False,visual_checks_performed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    if args.write:assert not TARGET.exists()
    result=run(args.write)
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(TARGET.read_text('utf-8'))
        for k in ('files','argument_scope','formal_reports','fresh_numbered_groups'):assert old[k]==result[k],k
    print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
