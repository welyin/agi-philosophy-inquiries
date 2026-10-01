"""Non-numbered inheritance audit: provenance/coverage, not new scientific proof."""
import argparse
import hashlib
import json
import re
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
NOTE=HERE/'cognitive_foundation_bridge_605.md'
ENTRY=HERE/'round606_drafts/foundation_inheritance_entry.md'
TARGET=HERE/'cognitive_foundation_bridge_605_results.json'
def run():
    base=core.read(HERE/'research_round_605_checks.json')
    assert base['all_reported_checks_passed'] and base['cumulative_numbered_tests']==3060
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,606):
        r=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in r[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    assert len(historical)==1824
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    rows=re.findall(r'^\|(B\d\d)\|(.+)$',NOTE.read_text('utf8'),re.M)
    assert [name for name,_ in rows]==[f'B{n:02}' for n in range(1,11)]
    paths={NOTE,ENTRY,Path(__file__).resolve()};links=0
    for f in (NOTE,ENTRY):
        for link in core.link_parser()(f.read_text('utf8')):
            path=(f.parent/link).resolve()
            assert path.exists() or path==TARGET.resolve(),(f.name,link)
            if path!=TARGET:paths.add(path)
            links+=1
    sources={str(f.relative_to(HERE)) if f.is_relative_to(HERE) else
             '../'+str(f.relative_to(HERE.parent)):core.digest(f) for f in sorted(paths)}
    return dict(date='2026-10-01',kind='non-numbered early-foundation inheritance audit',
        scientific_baseline=605,cumulative_scientific_checks_unchanged=3060,
        new_scientific_rounds=0,new_scientific_tests=0,inheritance_rows=10,
        evidence_hashes=sources,historical_evidence_files_verified=len(historical),
        local_links_checked=links,broken_links=0,full_230_round_reaudit=False,
        all_FUCP_conditions_in_current_model_proved=False,
        autonomous_design_still_deferred=True,active_goal_unchanged=True,
        document_review_completed=True,independent_review=False,
        all_document_checks_passed=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert core.read(TARGET)==r
    print(json.dumps({k:r[k] for k in ('scientific_baseline','inheritance_rows','historical_evidence_files_verified','all_document_checks_passed')}))
