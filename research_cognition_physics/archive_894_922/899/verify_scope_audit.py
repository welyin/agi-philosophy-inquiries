"""Verify audit provenance and links; does not rerun or certify mathematical results."""
from pathlib import Path
import argparse, hashlib, json, re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
AUDIT=HERE/'drafts/effective_scope_joint_audit.md'
RECEIPT=HERE/'drafts/working_checks.json'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(writing=False):
    frozen={}
    for n in range(776,899):
        d=json.loads((STAGE/f'{n}/research_round_{n}_checks.json').read_text('utf-8'))
        for key in ('frozen_inputs','new_scientific_and_entry_files'):
            for name,digest in d[key].items():
                assert name not in frozen or frozen[name]==digest, name
                frozen[name]=digest
    for name in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json'):
        d=json.loads((STAGE/name).read_text('utf-8'))
        for name,digest in d['files'].items():
            assert name not in frozen or frozen[name]==digest, name
            frozen[name]=digest
    for name,digest in frozen.items(): assert sha(ROOT/name)==digest, name
    navs=[STAGE.parent/f for f in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    navs += [STAGE/f for f in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for p in [AUDIT,*navs]:
        body=re.sub(r'\$\$.*?\$\$','',p.read_text('utf-8-sig'),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',body):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0].strip('<>')).resolve()
            assert q.exists() or (writing and q==RECEIPT.resolve()),(p,target)
            links+=1
    last=json.loads((STAGE/'898/research_round_898_checks.json').read_text('utf-8'))
    assert last['formal_reports']==898
    assert last['cumulative_numbered_test_groups_from_897']==3683
    assert not (STAGE/'research_note_899.md').exists()
    sources=[STAGE/f'research_note_{n}.md' for n in (860,869,894,895,896,897,898)]
    sources += [STAGE/'873/drafts/共同模型阶段报告_截至872.md',STAGE.parent/'archive_742_763/research_note_761.md',STAGE/'898/research_round_898_checks.json',HERE/'drafts/STATUS.md']
    tracked=[AUDIT,Path(__file__),*sources]
    return dict(date='2026-10-06',status='working audit; not a scientific round',formal_reports=898,
        cumulative_numbered_groups=3683,fresh_numbered_groups=0,
        historical_frozen_unique_files_verified=len(frozen),local_links_checked=links,
        numerical_experiments_rerun=False,mathematical_peer_review=False,
        files={str(p.relative_to(ROOT)):sha(p) for p in tracked})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    if args.write: assert not RECEIPT.exists()
    result=verify(args.write)
    if args.write: RECEIPT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: assert result['files']==json.loads(RECEIPT.read_text('utf-8'))['files']
    print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
