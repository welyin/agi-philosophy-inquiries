"""Integrate three independent complete research rounds based on frozen round 338."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=BASE.parent
TARGET=HERE/'round339_341_integration_checks.json'
ROUNDS=(339,340,341)
STEMS={339:'multispecies_characteristics_audit',
       340:'variable_disformal_coupling_audit',
       341:'disformal_validity_window_audit'}


def verify(pending_report=False):
    protected,old,counts={},{},{}
    for number in range(231,342):
        checked=core.read(HERE/f'research_round_{number}_checks.json')
        assert checked['all_reported_checks_passed']
        for name,sha in checked['new_file_hashes'].items():
            assert core.digest(HERE/name)==sha,name
            if name in protected:
                assert protected[name]==sha
            protected[name]=sha
            if number<=338:
                old[name]=sha
        results=[name for name in checked['new_file_hashes'] if name.endswith('_results.json')]
        assert len(results)==1,(number,results)
        saved=core.read(HERE/results[0])['checks']
        assert saved['failures']==saved['errors']==0
        counts[number]=saved['run']
    assert sum(counts[n] for n in counts if n<=338)==982
    assert len(protected)-6==334 and len(old)-6==325
    proposal=core.read(HERE/'five_node_growth_proposal_checks.json')['new_file_hashes']
    literature=core.read(HERE/'round276_277_integration_checks.json')['new_reference_hashes']
    priority=core.read(HERE/'spacetime_mainline_priority_checks.json')['new_file_hashes']
    for hashes in (proposal,literature,priority):
        for name,sha in hashes.items():
            assert core.digest(HERE/name)==sha,name
    from verify_effective_gravity_rounds import maintenance
    maintenance_hashes=maintenance()
    nav=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md',HERE/'README.md']
    notes=[HERE/f'research_note_{n}.md' for n in ROUNDS]
    review=HERE/'gr_assumption_dependency_review.md'
    parser=core.link_parser(); links=0
    for path in nav+notes+[review]:
        content=path.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in content)
        assert not re.search(rb'\r(?!\n)',path.read_bytes()),path
        for link in parser(content):
            resolved=(path.parent/link).resolve()
            assert resolved.exists() or (pending_report and resolved==TARGET.resolve()),(path.name,link)
            links+=1
    formulas=sum(core.text_checks(p)['display_formulas'] for p in notes)
    # Each new round can be reproduced without importing another new round.
    for number,stem in STEMS.items():
        tree=ast.parse((HERE/(stem+'.py')).read_text(encoding='utf-8'))
        imports=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):
                imports.extend(x.name.split('.')[0] for x in node.names)
            if isinstance(node,ast.ImportFrom) and node.module:
                imports.append(node.module.split('.')[0])
        forbidden={s for n,s in STEMS.items() if n!=number}
        assert not (set(imports)&forbidden),(number,imports)
        text=(HERE/f'research_note_{number}.md').read_text(encoding='utf-8')
        for other in ROUNDS:
            if number!=other:
                assert f'(research_note_{other}.md)' not in text,(number,other)
    for p in (BASE/'research_direction.md',BASE/'RESEARCH_STATE.md'):
        assert '完整研究轮次并行' in p.read_text(encoding='utf-8')
    for name in ('verify_parallel_physics_rounds.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    top=sorted(p.name for p in BASE.iterdir() if p.is_file())
    assert top==sorted(['README.md','research_direction.md','RESEARCH_STATE.md',
                       '可组合认知结构与复量子状态空间_阶段论文.md',
                       '可组合认知结构与有限维量子理论_阶段论文.md'])
    legacy_dir=BASE/'archive_223_230'
    sys.path.insert(0,str(legacy_dir))
    spec=importlib.util.spec_from_file_location('closed_archive',legacy_dir/'verify_archive.py')
    legacy=importlib.util.module_from_spec(spec); spec.loader.exec_module(legacy)
    closed=legacy.verify(False,False)
    assert closed['all_reported_checks_passed']
    diff=subprocess.run(['git','diff','--check'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    assert diff.returncode==0,diff.stdout+diff.stderr
    return {'date':'2026-09-23','rounds':list(ROUNDS),
            'parallel_rounds_share_frozen_base_through':338,
            'cross_round_dependencies_in_batch':False,
            'fresh_tests_by_round':{n:counts[n] for n in ROUNDS},
            'fresh_tests':sum(counts[n] for n in ROUNDS),
            'stage_saved_tests':sum(counts.values()),
            'science_hashes_verified_231_341':len(protected)-6,
            'protected_text_maintenance_files':5,
            'prior_additive_text_maintenance_hashes':maintenance_hashes,
            'frozen_unnumbered_dependency_review_files':1,
            'dependency_review_sha256':core.digest(review),
            'unchanged_prior_science_hashes_231_338':len(old)-6,
            'additional_unnumbered_proposal_files':len(proposal),
            'frozen_literature_files':len(literature),
            'frozen_priority_decision_files':len(priority),
            'total_protected_evidence_hashes':len(protected)+len(proposal)+len(literature)+len(priority),
            'unchanged_prior_evidence_hashes':len(old)+len(proposal)+len(literature)+len(priority),
            'navigation_and_notes_checked':len(nav+notes+[review]),
            'local_links_checked':links,'broken_links':0,
            'new_scientific_display_formulas_checked':formulas,
            'top_level_files':top,'visual_rendering_performed':False,
            'legacy_archive_readonly_check_passed':True,'legacy_science_tests_rerun':False,
            'legacy_stage1_versions_verified':closed['stage1_original_versions_verified'],
            'legacy_stage2_scientific_files_verified':closed['frozen_scientific_files_verified'],
            'git_diff_check_passed':True,'all_reported_checks_passed':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        if TARGET.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
