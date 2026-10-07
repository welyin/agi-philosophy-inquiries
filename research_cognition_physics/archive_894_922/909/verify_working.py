"""909 working checkpoint: preserve908 and record the sourced momentum interface."""
from pathlib import Path
import argparse,ast,hashlib,json,re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'drafts/working_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    old=json.loads((STAGE/'908/research_round_908_checks.json').read_text('utf-8'));prior={}
    for key in ('frozen_inputs','new_scientific_and_entry_files'):
        for name,digest in old[key].items():assert sha(ROOT/name)==digest,name;prior[name]=digest
    assert not (STAGE/'research_note_909.md').exists()
    paths=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/retarded_source_working.md',STAGE/'908/research_round_908_checks.json',
        STAGE/'908/magnetic_reference_source.py',STAGE/'905/canonical_source_bridge.py',STAGE/'904/coupled_boson_tangent.py']
    for p in paths:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'),filename=str(p))
    docs=[HERE/'drafts/STATUS.md',HERE/'drafts/retarded_source_working.md',STAGE/'research_note_908.md',STAGE/'_shared/notes/unified_physics_condition_ledger_current.md']
    docs += [STAGE/p for p in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md')]
    docs += [STAGE.parent/p for p in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    links=0
    for p in docs:
        prose=re.sub(r'\$\$.*?\$\$','',p.read_text('utf-8-sig'),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0].strip('<>')).resolve()
            assert q.exists() or (writing and q==TARGET.resolve()),(p,target)
            links+=1
    return dict(date='2026-10-06',status='working',all_checkpoint_checks_passed=True,
        formal_rounds=908,cumulative_numbered_groups=3693,new_scientific_groups=0,
        prior908_frozen_files_verified=len(prior),local_links_checked=links,
        original_source_momentum_shift_derivation_saved=True,actual_shifted_source_solver_implemented=False,
        original872_forced_response_computed=False,
        files={str(p.relative_to(ROOT)):sha(p) for p in paths},full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run(a.write)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out['files']==json.loads(TARGET.read_text('utf-8'))['files']
    print(json.dumps({k:v for k,v in out.items() if k!='files'},ensure_ascii=False,indent=2))
