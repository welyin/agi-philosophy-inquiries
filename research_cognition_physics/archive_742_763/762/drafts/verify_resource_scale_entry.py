"""Check the 762 working entry without increasing completed round counts."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import verify_interaction_rounds as core
import resource_scale_entry as entry
TARGET=HERE/'resource_scale_entry_checks.json'


def verify():
    history=dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for number in range(584,762):
        receipt=core.read(ARCHIVE/f'research_round_{number}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in history or history[name]==digest, name
                history[name]=digest
    history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3637
    for name,digest in history.items():
        assert core.digest(ARCHIVE/name)==digest,name
    result=entry.run()
    assert result==core.read(entry.TARGET)
    for name,digest in result['dependency_hashes'].items():
        assert core.digest(ARCHIVE/name)==digest
    report=HERE/'research_note_762_working.md'
    text=report.read_text('utf8')
    links=list(core.link_parser()(text))
    for link in links:
        path=(HERE/link).resolve()
        assert path.exists() or path==TARGET.resolve(),link
    assert text.count('$$')==10
    assert '正式完成仍为761／3470' in text
    files=('research_note_762_working.md','resource_scale_entry.py',
           'resource_scale_entry_results.json','verify_resource_scale_entry.py',
           'publish_resource_scale_entry.py')
    for name in files:
        if name.endswith('.py'):
            ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-04', kind='round762_working_entry',
        latest_completed_scientific_round=761, cumulative_numbered_tests=3470,
        new_numbered_scientific_tests=0, original_diagnostic_groups=2,
        historical_protected_files=len(history), all_historical_hashes_unchanged=True,
        working_file_hashes={n:core.digest(HERE/n) for n in files},
        frozen_entry_status_sha256=core.digest(HERE/'STATUS.md'),
        saved_diagnostics_reproduced=True, local_links_checked=len(links),
        displayed_formulas=5, broken_links=0, new_goal_or_automation=False,
        full_graph_dynamics_simulated=False, full_continuum_model_proven=False,
        prior_coarse_results_and_525_reused=True,
        primary_source_scope='Hu-Ver daguer 0802.0658 section 3.3.1 used only for hierarchy; '
                             'Egorov and Gaussian results do not certify full quantum gravity.',
        independent_agent_review=False, visual_check=False,
        all_reported_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--write-checks',action='store_true')
    args=p.parse_args()
    result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==core.read(TARGET)
    print(json.dumps({k:result[k] for k in ('kind','latest_completed_scientific_round',
                                          'historical_protected_files','all_reported_checks_passed')}))

