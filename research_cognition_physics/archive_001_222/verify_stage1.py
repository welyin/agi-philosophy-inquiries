"""Verify the 222-round archive, original evidence and relocated references.

This checks bytes, coverage, saved regression status and local links. It does
not prove the mathematical claims or rerun the archived numerical tests.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
from urllib.parse import unquote

ARCHIVE = Path(__file__).resolve().parent
HERE = ARCHIVE.parent
ROOT = HERE.parent


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_links(text):
    # Bracketed operators followed by (rho) in LaTeX are not Markdown links.
    text = re.sub(r'\\\[.*?\\\]', '', text, flags=re.S)
    text = re.sub(r'\$\$.*?\$\$', '', text, flags=re.S)
    text = re.sub(r'^```.*?^```\s*$', '', text, flags=re.S | re.M)
    inline = re.findall(r'(?<!!)\[[^\]]*\]\(([^)]+)\)', text)
    references = re.findall(r'^\[[^\]]+\]:\s*(\S+)\s*$', text, flags=re.M)
    for target in inline+references:
        target = unquote(target.split('#')[0].strip().strip('<>'))
        if not target or re.match(r'^[a-zA-Z]+:', target):
            continue
        yield target


def verify(pending_report=False):
    manifest = read_json(ARCHIVE / 'ARCHIVE_MANIFEST.json')
    assert manifest['process_directory'] == 'research_process'
    for removed in ('research_physics_construction', 'research_information_geometry', 'research_cognition_physics'):
        assert not (ARCHIVE / removed).exists(), removed
    originals_by_logical_name = {}
    link_only_updates = []
    for entry in manifest['research_artifacts']:
        p = ARCHIVE / entry['archived_path']
        assert p.is_file() and p.stat().st_size == entry.get('current_bytes', entry['bytes']), str(p)
        assert sha256(p) == entry.get('current_sha256', entry['sha256']), str(p)
        original = p
        if 'original_backup_path' in entry:
            original = ARCHIVE / entry['original_backup_path']
            assert original.stat().st_size == entry['bytes'] and sha256(original) == entry['sha256']
            expected = original.read_bytes()
            allowed = {f'](../{name}/': f'](../../../{name}/' for name in
                       ('research_physics_construction', 'research_information_geometry')}
            for rule in entry['link_replacements']:
                assert allowed[rule['before']] == rule['after']
                before, after = rule['before'].encode(), rule['after'].encode()
                assert expected.count(before) == rule['count']
                expected = expected.replace(before, after)
            assert expected == p.read_bytes(), str(p)
            link_only_updates.append(entry['archived_path'])
        originals_by_logical_name[entry['original_path']] = original
        if p.suffix == '.json':
            read_json(p)
    for entry in manifest['context_files']:
        p = ARCHIVE / entry['path']
        assert sha256(p) == entry['sha256'], str(p)
    library = ARCHIVE / 'research_process'
    rounds = sorted(int(re.search(r'\d+', p.name).group()) for p in library.glob('research_note_*.md'))
    assert rounds == manifest['numbered_rounds'] == list(range(1, 223))
    index = read_json(ARCHIVE / 'ROUND_INDEX.json')
    assert [r['round'] for r in index['rounds']] == rounds
    for r in index['rounds']:
        assert sha256(ARCHIVE / r['note']) == r.get('current_sha256', r['sha256'])
        logical = r['note'].replace('research_process/', 'research_cognition_physics/', 1)
        assert sha256(originals_by_logical_name[logical]) == r['sha256']
    paper = HERE / '可组合认知结构与复量子状态空间_阶段论文.md'
    top_level_files = sorted(p.name for p in HERE.iterdir() if p.is_file())
    assert all(name in {'README.md', 'research_direction.md', 'RESEARCH_STATE.md'}
               or (name.endswith('.md') and '论文' in name) for name in top_level_files), top_level_files
    paper_text = paper.read_text(encoding='utf-8')
    paper_rounds = [int(n) for n in re.findall(r'^\| \[(\d+)\]\[r\d+\] \|', paper_text, re.M)]
    assert paper_rounds == rounds
    definitions = set(re.findall(r'^\[([^\]]+)\]:', paper_text, re.M))
    for label in re.findall(r'\[[^\]]+\]\[(r\d+)\]', paper_text):
        assert label in definitions, label

    old_missing = {(r['from'], r['target']) for r in manifest['preexisting_broken_local_links']}
    originals = [ARCHIVE/r['archived_path'] for r in manifest['research_artifacts']
                 if r['archived_path'].endswith('.md')]
    originals += [ARCHIVE/r['path'] for r in manifest['context_files'] if r['path'].endswith('.md')]
    live = [paper, HERE/'README.md', HERE/'research_direction.md', HERE/'RESEARCH_STATE.md',
            ARCHIVE/'README.md', ARCHIVE/'ROUND_INDEX.md', ROOT/'README.md',
            ROOT/'research_physics_construction/INHERITED_INDEX.md',
            ROOT/'research_physics_construction/NECESSITY_BACKLOG.md']
    link_count = 0
    accepted_old_missing = []
    for p in originals+live:
        text = p.read_text(encoding='utf-8-sig')
        assert not any(ord(ch)<32 and ch not in '\r\n\t' for ch in text), str(p)
        for target in local_links(text):
            link_count += 1
            dest = (p.parent/target).resolve()
            if not dest.exists():
                if pending_report and dest == (ARCHIVE/'stage1_archive_checks.json').resolve():
                    continue
                logical = p.relative_to(ARCHIVE).as_posix() if p.is_relative_to(ARCHIVE) else None
                assert (logical, target) in old_missing, (str(p), target)
                accepted_old_missing.append({'from':logical, 'target':target})

    inherited = read_json(ROOT/'research_physics_construction/inherited_results.json')
    for entry in inherited['artifacts']:
        candidate = originals_by_logical_name.get(entry['path'], ROOT/entry['path'])
        assert sha256(candidate) == entry['sha256'], entry['path']
    tests = read_json(ARCHIVE/'stage1_test_results.json')
    assert tests['success'] and tests['failures'] == tests['errors'] == 0
    compatibility = read_json(ARCHIVE/'stage1_compatibility_checks.json')
    assert compatibility['all_passed']
    layout_checks = read_json(ARCHIVE/'stage1_layout_checks.json')
    assert layout_checks['all_passed']
    produced = [paper, HERE/'README.md', HERE/'research_direction.md', HERE/'RESEARCH_STATE.md',
                ARCHIVE/'run_stage1.py', ARCHIVE/'verify_stage1.py', ARCHIVE/'stage1_test_results.json',
                ARCHIVE/'stage1_compatibility_checks.json', ARCHIVE/'stage1_layout_checks.json', ARCHIVE/'README.md',
                ARCHIVE/'ROUND_INDEX.md', ARCHIVE/'ROUND_INDEX.json', ARCHIVE/'ARCHIVE_MANIFEST.json',
                ARCHIVE/'PATH_MIGRATION.json', ARCHIVE/'paper_math_review_checks.json']
    return {'date':'2026-09-20', 'kind':'Archive support tools and reports moved into archive; verify original evidence, paper coverage and links; not round223.',
            'python':platform.python_version(), 'numbered_rounds':len(rounds),
            'top_level_files':top_level_files,
            'original_research_files_byte_identical':len(manifest['research_artifacts'])-len(link_only_updates),
            'research_files_with_link_only_updates':link_only_updates,
            'original_research_versions_preserved_and_verified':len(manifest['research_artifacts']),
            'deleted_external_snapshot_files':len(manifest['removed_context_files']),
            'external_reading_context_copies_verified':len(manifest['context_files']),
            'original_c0_artifact_hashes_verified':len(inherited['artifacts']),
            'saved_test_report_scope':'Previous full numerical regression; not rerun by this verifier.',
            'archived_unittest_modules':tests['module_count'], 'archived_tests_passed':tests['tests_run'],
            'archived_tests_skipped':tests['skipped'],
            'saved_compatibility_report_scope':'Original stage closure report; current layout validation is recorded separately.',
            'compatibility_checks':compatibility, 'layout_revision_checks':layout_checks,
            'markdown_files_checked':len(originals+live), 'local_links_checked':link_count,
            'new_broken_local_links':0, 'preexisting_broken_links_preserved':accepted_old_missing,
            'paper_round_coverage':len(paper_rounds),
            'review_scope':'Current file integrity, original versions, all-round coverage and links; reads previous regression reports without rerunning tests.',
            'new_file_hashes':{str(p.relative_to(HERE).as_posix()):sha256(p) for p in produced},
            'all_reported_checks_passed':True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    report = verify(pending_report=args.write_checks)
    if args.write_checks:
        (ARCHIVE/'stage1_archive_checks.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in {'new_file_hashes', 'preexisting_broken_links_preserved'}}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
