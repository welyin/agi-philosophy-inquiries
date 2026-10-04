"""Integrate finite-range dictionaries and modular inclusion prerequisites."""
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
TARGET=HERE/'round368_369_integration_checks.json'
ROUNDS=(368,369)
STEMS={368:'quasi_local_encoding_audit',369:'modular_inclusion_limit_audit'}
BASES={368:367,369:367}


def verify(pending_report=False):
    protected,old,counts={},{},{}
    for number in range(231,370):
        checked=core.read(HERE/f'research_round_{number}_checks.json')
        assert checked['all_reported_checks_passed']
        for name,sha in checked['new_file_hashes'].items():
            assert core.digest(HERE/name)==sha,name
            if name in protected:
                assert protected[name]==sha
            protected[name]=sha
            if number<=367:
                old[name]=sha
        results=[name for name in checked['new_file_hashes'] if name.endswith('_results.json')]
        assert len(results)==1,(number,results)
        saved=core.read(HERE/results[0])['checks']
        assert saved['failures']==saved['errors']==0
        counts[number]=saved['run']
        if number in ROUNDS:
            assert checked['parallel_batch']==list(ROUNDS)
            assert checked['scientific_base_through_round']==BASES[number]
            assert checked['batch_scientific_dependencies']==[]
    assert sum(counts[n] for n in counts if n<=367)==1403
    assert len(protected)-6==418 and len(old)-6==412
    proposal=core.read(HERE/'five_node_growth_proposal_checks.json')['new_file_hashes']
    literature=core.read(HERE/'round276_277_integration_checks.json')['new_reference_hashes']
    priority=core.read(HERE/'spacetime_mainline_priority_checks.json')['new_file_hashes']
    for hashes in (proposal,literature,priority):
        for name,sha in hashes.items():
            assert core.digest(HERE/name)==sha,name
    from verify_effective_gravity_rounds import maintenance
    maintenance_hashes=maintenance()
    correction=core.read(HERE/'round346_text_correction.json')
    raw=(HERE/correction['original']).read_bytes()
    corrected=HERE/correction['corrected']
    assert core.digest(HERE/correction['original'])==correction['original_sha256']
    assert core.digest(corrected)==correction['corrected_sha256']
    old_phrase=correction['old_phrase'].encode('utf-8')
    assert raw.count(old_phrase)==correction['occurrences']==1
    assert corrected.read_bytes()==correction['banner'].encode('utf-8')+raw.replace(
        old_phrase,correction['new_phrase'].encode('utf-8'))
    assert core.text_checks(corrected)['display_formulas']==12
    current_maintenance={p.name:core.digest(p) for p in
                         (corrected,HERE/'round346_text_correction.json')}
    assert core.read(HERE/'round345_347_integration_checks.json')['new_text_maintenance_hashes']==current_maintenance
    from verify_relational_geometry_rounds import preserved_draft
    round363_draft=preserved_draft()
    for number in (363,364):
        assert core.read(HERE/f'research_round_{number}_checks.json')['preserved_round363_draft_hashes']==round363_draft
    nav=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md',HERE/'README.md']
    notes=[corrected if n==346 else HERE/f'research_note_{n}.md' for n in ROUNDS]
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
    # The unfrozen round-344 schema fix added only scope metadata; retain both originals.
    draft_dir=HERE/'round344_drafts'
    old_result=core.read(draft_dir/'before_scope_schema_results.json')
    current_result=core.read(HERE/'covariant_gravity_uniqueness_audit_results.json')
    assert {k:v for k,v in current_result.items() if k!='scope'}==old_result
    old_source=(draft_dir/'before_scope_schema.py').read_text(encoding='utf-8')
    new_source=(HERE/'covariant_gravity_uniqueness_audit.py').read_text(encoding='utf-8')
    assert ''.join(line for line in new_source.splitlines(keepends=True)
                   if "'scope':" not in line)==old_source
    draft_hashes={path.relative_to(HERE).as_posix():core.digest(path)
                  for path in sorted(draft_dir.iterdir()) if path.is_file()}
    assert core.read(HERE/'round342_344_integration_checks.json')['preserved_schema_draft_hashes']==draft_hashes
    if TARGET.exists():
        assert core.read(TARGET)['preserved_schema_draft_hashes']==draft_hashes
    for name in ('verify_transport_modular_rounds.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    top=sorted(p.name for p in BASE.iterdir() if p.is_file())
    assert top==sorted(['README.md','research_direction.md','RESEARCH_STATE.md',
                       '可组合认知结构与复量子状态空间_阶段论文.md',
                       '可组合认知结构与有限维量子理论_阶段论文.md'])
    legacy_dir=BASE/'archive_223_230'
    sys.path.insert(0,str(legacy_dir))
    spec=importlib.util.spec_from_file_location('closed_archive',legacy_dir/'verify_archive.py')
    legacy=importlib.util.module_from_spec(spec); spec.loader.exec_module(legacy)
    # Reuse the legacy verifier's explicit pending-report allowance for this report.
    # This changes only the loaded module, never the preserved historical source.
    legacy.REPORT=TARGET
    closed=legacy.verify(False,pending_report)
    assert closed['all_reported_checks_passed']
    diff=subprocess.run(['git','diff','--check'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    assert diff.returncode==0,diff.stdout+diff.stderr
    return {'date':'2026-09-23','rounds':list(ROUNDS),
            'execution_mode':'independent complete rounds in parallel',
            'scientific_base_through_round_by_round':BASES,
            'cross_round_dependencies_in_batch':False,
            'additional_frozen_dependency_rounds':{},
            'fresh_tests_by_round':{n:counts[n] for n in ROUNDS},
            'fresh_tests':sum(counts[n] for n in ROUNDS),
            'stage_saved_tests':sum(counts.values()),
            'science_hashes_verified_231_369':len(protected)-6,
            'protected_text_maintenance_files':7,
            'frozen_text_maintenance_hashes':current_maintenance,
            'prior_additive_text_maintenance_hashes':maintenance_hashes,
            'frozen_unnumbered_dependency_review_files':1,
            'dependency_review_sha256':core.digest(review),
            'unchanged_prior_science_hashes_231_367':len(old)-6,
            'additional_unnumbered_proposal_files':len(proposal),
            'frozen_literature_files':len(literature),
            'frozen_priority_decision_files':len(priority),
            'preserved_schema_draft_hashes':draft_hashes,
            'schema_fix_changed_scientific_values':False,
            'preserved_round363_draft_hashes':round363_draft,
            'draft_tests_counted_as_new_science':False,
            'total_protected_evidence_hashes':len(protected)+len(proposal)+len(literature)+len(priority)+len(draft_hashes)+len(current_maintenance)+len(round363_draft),
            'unchanged_prior_evidence_hashes':len(old)+len(proposal)+len(literature)+len(priority)+len(draft_hashes)+len(current_maintenance)+len(round363_draft),
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
