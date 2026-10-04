"""Integration check for rounds 295-296; previous evidence remains read-only."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
TARGET = HERE/'round295_296_integration_checks.json'


def verify():
    protected, old, counts = {}, {}, {}
    for number in range(231,297):
        checked = core.read(HERE/f'research_round_{number}_checks.json')
        assert checked['all_reported_checks_passed']
        for name,sha in checked['new_file_hashes'].items():
            assert core.digest(HERE/name)==sha,name
            if name in protected:
                assert protected[name]==sha
            protected[name]=sha
            if number<=294:
                old[name]=sha
        results=[name for name in checked['new_file_hashes'] if name.endswith('_results.json')]
        assert len(results)==1,(number,results)
        saved=core.read(HERE/results[0])['checks']
        assert saved['failures']==saved['errors']==0
        counts[number]=saved['run']
    assert sum(counts.values())==576
    proposal=core.read(HERE/'five_node_growth_proposal_checks.json')['new_file_hashes']
    literature=core.read(HERE/'round276_277_integration_checks.json')['new_reference_hashes']
    for hashes in (proposal,literature):
        for name,sha in hashes.items():
            assert core.digest(HERE/name)==sha,name
    nav=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md',HERE/'README.md']
    notes=[HERE/f'research_note_{n}.md' for n in (295,296)]
    parser=core.link_parser()
    links=0
    for path in nav+notes:
        content=path.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\r\n\t' for c in content)
        for link in parser(content):
            assert (path.parent/link).resolve().exists(),(path.name,link)
            links+=1
    formulas=sum(core.text_checks(p)['display_formulas'] for p in notes)
    assert formulas==14
    for p in (BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',HERE/'README.md'):
        assert '下一轮：共同钟尺的因果在线更新' in p.read_text(encoding='utf-8')
    for name in ('finite_clock_witness_audit.py','clock_window_certificate_audit.py',
                 'verify_clock_evidence_rounds.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    top=sorted(p.name for p in BASE.iterdir() if p.is_file())
    assert top==sorted(['README.md','research_direction.md','RESEARCH_STATE.md',
                       '可组合认知结构与复量子状态空间_阶段论文.md',
                       '可组合认知结构与有限维量子理论_阶段论文.md'])
    legacy_dir=BASE/'archive_223_230'
    sys.path.insert(0,str(legacy_dir))
    spec=importlib.util.spec_from_file_location('closed_archive',legacy_dir/'verify_archive.py')
    legacy=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    closed=legacy.verify(False,False)
    assert closed['all_reported_checks_passed']
    diff=subprocess.run(['git','diff','--check'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    assert diff.returncode==0,diff.stdout+diff.stderr
    return {'date':'2026-09-22','rounds':[295,296],
            'fresh_tests':counts[295]+counts[296],'stage_saved_tests':sum(counts.values()),
            'science_hashes_verified_231_296':len(protected)-3,
            'protected_text_maintenance_files':3,
            'unchanged_prior_science_hashes_231_294':len(old)-3,
            'additional_unnumbered_proposal_files':len(proposal),
            'frozen_literature_files':len(literature),
            'total_protected_evidence_hashes':len(protected)+len(proposal)+len(literature),
            'unchanged_prior_evidence_hashes':len(old)+len(proposal)+len(literature),
            'navigation_and_notes_checked':len(nav+notes),'local_links_checked':links,
            'broken_links':0,'new_text_display_formulas_checked':formulas,
            'top_level_files':top,'visual_rendering_performed':False,
            'legacy_archive_readonly_check_passed':True,'legacy_science_tests_rerun':False,
            'legacy_stage1_versions_verified':closed['stage1_original_versions_verified'],
            'legacy_stage2_scientific_files_verified':closed['frozen_scientific_files_verified'],
            'git_diff_check_passed':True,'all_reported_checks_passed':True}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify()
    if args.write_checks:
        if TARGET.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

