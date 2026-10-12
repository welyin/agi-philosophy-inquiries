"""Read-only current-stage integrity and closure checks; no image inspection.

Both research scripts have separate default read-only numerical checks. This
verifier checks their signed assets and the scope-bearing archival documents.
"""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
ROOT = BASE.parent
PAPER = BASE / '06_认知组合公理与局域量子过程_阶段论文.md'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-closure', action='store_true')
    args = parser.parse_args()
    ledger = json.loads((HERE / 'dependency_ledger.json').read_text(encoding='utf-8'))
    assert ledger['accepted_rounds'] == [1157, 1158]
    assert ledger['stage_accepted_rounds'] == 2
    assert ledger['stage_numerical_check_groups'] == 8
    assert ledger['scientific_calibration_increment'] == 0
    assert ledger['in_progress_rounds'] == []
    hashes = 0
    for n in (1157, 1158):
        folder = STAGE / str(n)
        acceptance = json.loads((folder / 'acceptance.json').read_text(encoding='utf-8'))
        assert acceptance['status'] == 'PASS'
        for name, expected in acceptance['core_sha256'].items():
            path = (folder / name).resolve()
            assert path.is_relative_to(STAGE)
            assert sha(path) == expected, (n, name)
            hashes += 1
        assert sha(folder / acceptance['review_path']) == acceptance['review_sha256']
    src = json.loads((STAGE / '1157/source_hashes.json').read_text(encoding='utf-8'))
    for entry in src['files']:
        assert sha(ROOT / entry['path']) == entry['sha256'], entry['path']
    src2 = json.loads((STAGE / '1158/source_hashes.json').read_text(encoding='utf-8'))
    for name, expected in src2['sha256'].items():
        assert sha(STAGE / '1158' / name) == expected, name
    docs = list(STAGE.rglob('*.md')) + [PAPER]
    links = 0
    for path in docs:
        body = path.read_text(encoding='utf-8-sig')
        assert body.count('$$') % 2 == 0, ('math delimiters', path)
        for target in re.findall(r'\]\(([^)\r\n]+)\)', body):
            if target.startswith(('http:', 'https:', 'mailto:', '#')):
                continue
            target = target.strip('<>').split('#')[0]
            assert (path.parent / target).resolve().exists(), (path, target)
            links += 1
    scripts = list(STAGE.rglob('*.py'))
    for path in scripts:
        ast.parse(path.read_text(encoding='utf-8'))
    # Keep the single pre-existing reversible navigation layer; original
    # historical navigation must remain recoverable after the new stage.
    sys.path.insert(0, str(ROOT / 'scripts'))
    from research_layout import relocation_original_bytes
    layer = json.loads((BASE / '_migration/field_active_1120/layer.json').read_text(encoding='utf-8'))
    baseline = json.loads((BASE / '_migration/merge_1106_1119/manifest.json').read_text(encoding='utf-8'))
    prior = {e['destination']: e for e in baseline['entries']}
    for entry in layer['entries']:
        current = (ROOT / entry['destination']).read_bytes()
        assert hashlib.sha256(current).hexdigest() == entry['current_sha256']
        original = relocation_original_bytes(current, entry)
        assert hashlib.sha256(original).hexdigest() == prior[entry['original']]['current_sha256']
    closure_hashes = 0
    if args.require_closure:
        assert STAGE.name == 'archive_1157_1158'
        assert not (BASE / 'archive_1157_').exists()
        assert ledger['new_adopted_axioms'] == ['CC1', 'CC2', 'CC3']
        assert ledger['next_round'] is None
        adopted = json.loads((HERE / 'adoption_record.json').read_text(encoding='utf-8'))
        assert adopted['status'] == 'ADOPTED'
        for name, expected in adopted['core_sha256'].items():
            assert sha(HERE / name) == expected, name
            closure_hashes += 1
        manifest = json.loads((HERE / 'directory_migration.json').read_text(encoding='utf-8'))
        assert manifest['status'] == 'VERIFIED_MOVE'
        for entry in manifest['entries']:
            assert (ROOT / entry['new_path']).exists(), entry['new_path']
            if entry['relative'] not in manifest['mutable_after_move']:
                assert sha(ROOT / entry['new_path']) == entry['sha256'], entry['new_path']
        assert 'CC1—CC3正式采用' in PAPER.read_text(encoding='utf-8')
        assert (HERE / 'stage_completion.md').exists()
        assert (HERE / 'completion_scope_review.md').exists()
    print(json.dumps(dict(status='PASS',accepted_rounds=2,numerical_groups=8,
        round_core_hashes=hashes,independent_review_hashes=2,
        frozen_sources=len(src['files'])+len(src2['sha256']),markdown_files=len(docs),
        local_links=links,python_ast=len(scripts),reversible_navigation_files=4,
        closure_hashes=closure_hashes,closure_required=args.require_closure,
        scope='Integrity and scope records only; numeric scripts and analytic reviews are separate evidence.'),
        ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
