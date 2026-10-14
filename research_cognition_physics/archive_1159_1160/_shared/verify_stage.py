"""Read-only evidence/asset check for the finite stage 1159--1160.

This does not repeat scientific experiments or claim to verify their proofs.
--record-once writes the first asset receipt exclusively. Default operation
recomputes the checks and, when present, compares that receipt without writing.
The goal completion receipt is checked if present and may be created later.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
ROOT = BASE.parent
PAPER = BASE / '07_局域相对论过程的区域收支与量子记录_阶段论文.md'
RECEIPT = HERE / 'stage_assets_check.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def ensure(condition, message):
    if not condition:
        raise ValueError(message)


def calculate():
    hashes = {}

    def bind(path):
        path = path.resolve()
        ensure(path.is_file(), f'Missing asset: {path}')
        key = path.relative_to(ROOT).as_posix()
        hashes[key] = digest(path)
        return hashes[key]

    for number in (1159, 1160):
        folder = STAGE / str(number)
        acceptance = read_json(folder / 'acceptance.json')
        ensure(acceptance['status'] == 'PASS', f'Round {number} not accepted')
        ensure(acceptance['check_groups'] == 4, 'Unexpected science count')
        review = folder / 'independent_review.md'
        ensure(bind(review) == acceptance['review_sha256'], 'Review changed')
        review_text = review.read_text(encoding='utf-8-sig')
        for name, expected in acceptance['core_sha256'].items():
            observed = bind(folder / name)
            ensure(observed == expected, f'Accepted core changed: {number}/{name}')
            ensure(expected in review_text, f'Independent binding missing: {name}')
        if number == 1160:
            ensure(bind(folder / 'scope_review.md') == acceptance['scope_review_sha256'],
                   'Internal scope review changed')
            ensure(bind(folder / 'NEXT.md') == acceptance['NEXT_sha256'], 'NEXT changed')
        bind(folder / 'acceptance.json')

    inherited = read_json(STAGE / '1159/source_hashes.json')['sha256']
    for name, expected in inherited.items():
        ensure(digest(BASE / name) == expected, f'Inherited scientific file changed: {name}')
    results = read_json(STAGE / '1160/results.json')
    ensure(results['source_sha256'] == digest(BASE / 'archive_1106_1119/1111/check.py'),
           '1160 inherited runtime changed')
    ensure(results['code_sha256'] == digest(STAGE / '1160/check.py'), '1160 code changed')

    paper_review = HERE / 'stage_paper_review.md'
    scope_review = HERE / 'completion_scope_review.md'
    paper_review_text = paper_review.read_text(encoding='utf-8-sig')
    scope_review_text = scope_review.read_text(encoding='utf-8-sig')
    ensure('PASS' in paper_review_text and 'PASS' in scope_review_text,
           'Stage reviews have not passed')
    for target in (PAPER, HERE / 'common_process_mapping.md', HERE / 'field_bridge_after_1160.md'):
        expected = bind(target)
        ensure(expected in paper_review_text, f'Paper review lacks current binding: {target.name}')
        ensure(expected in scope_review_text, f'Scope review lacks current binding: {target.name}')
    for target in (paper_review, HERE / 'stage_completion.md', HERE / 'NEXT.md',
                   STAGE / '1159/acceptance.json', STAGE / '1160/acceptance.json'):
        ensure(bind(target) in scope_review_text, f'Scope binding changed: {target.name}')
    for number in (1015, 1016, 1017, 1018, 1019, 1027, 1028, 1029):
        target = BASE / f'archive_1009_1043/research_note_{number}.md'
        ensure(digest(target) in paper_review_text, f'Field source binding changed: {number}')

    ledger = read_json(HERE / 'dependency_ledger.json')
    ensure(ledger['latest_accepted_round'] == 1160, 'Stale latest round')
    ensure(ledger['accepted_stage_rounds'] == [1159, 1160], 'Wrong accepted rounds')
    ensure(ledger['stage_check_groups'] == 8, 'Wrong group count')
    ensure(ledger['new_adopted_axioms'] == [], 'Unexpected axiom adoption')
    ensure(ledger['new_scientific_calibration'] == 0, 'Unexpected calibration')

    # Source checks only; no rendered-image inspection. Fenced code is not math.
    markdown = sorted(STAGE.rglob('*.md')) + [PAPER]
    link_count = 0
    for path in markdown:
        raw = path.read_bytes()
        ensure(b'\r' not in raw.replace(b'\r\n', b''), f'Lone CR in {path.name}')
        text = raw.decode('utf-8-sig')
        text = re.sub(r'```[^\n]*\n.*?```', '', text, flags=re.S)
        ensure(text.count('$$') % 2 == 0, f'Unbalanced display math: {path.name}')
        ensure(not re.search(r'\\[\[\]()]', text), f'Legacy math delimiter: {path.name}')
        without_display = re.sub(r'\$\$.*?\$\$', '', text, flags=re.S)
        ensure(len(re.findall(r'(?<!\\)\$', without_display)) % 2 == 0,
               f'Unbalanced inline math: {path.name}')
        for href in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', text):
            href = href.strip('<>')
            if re.match(r'^(?:https?://|mailto:|#)', href):
                continue
            target = unquote(href.split('#', 1)[0])
            if not target:
                continue
            # The single shared notice is rendered from BASE, not HERE.
            parent = BASE if path == HERE / 'navigation_notice.md' else path.parent
            resolved = (parent / target).resolve()
            link_count += 1
            if resolved in (HERE / 'completion_receipt.json', RECEIPT):
                continue  # Output receipts are emitted only after prerequisites pass.
            ensure(resolved.exists(), f'Broken local link: {path.name} -> {href}')

    for path in (STAGE / 'README.md', STAGE / '文件索引.md', HERE / 'goal.md',
                 HERE / 'goal_start_receipt.json', HERE / 'stage_completion.md',
                 HERE / 'NEXT.md', HERE / 'completion_scope_precheck.md',
                 paper_review, scope_review, HERE / 'dependency_ledger.json',
                 HERE / 'progress_checks.json', HERE / 'final_navigation_checks.json',
                 HERE / 'navigation_notice.md', HERE / 'update_navigation.py',
                 HERE / 'verify_stage.py', ROOT / 'scripts/run_research_local_dynamics_active.py'):
        bind(path)

    return dict(schema='local_dynamics_1159_1160_assets_v1', status='PASS',
                accepted_rounds=[1159, 1160], science_check_groups=8,
                new_scientific_calibration=0, new_adopted_axioms=0,
                inherited_source_checks=len(inherited) + 1 + 8,
                markdown_files_checked=len(markdown), local_links_checked=link_count,
                checks='Accepted/reviewed versions, inherited sources, links and math source; not a new experiment',
                sha256=dict(sorted(hashes.items())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.record_once:
        with RECEIPT.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        outcome = 'RECORDED_ONCE'
    elif RECEIPT.exists():
        ensure(read_json(RECEIPT) == result, 'Asset receipt mismatch; investigate, do not overwrite')
        outcome = 'READ_ONLY_MATCH'
    else:
        outcome = 'CHECKED_NOT_RECORDED'
    completion = HERE / 'completion_receipt.json'
    if completion.exists():
        receipt = read_json(completion)
        ensure(receipt['goal']['status'] == 'complete', 'Goal receipt is not complete')
        ensure(receipt['stage_assets_sha256'] == digest(RECEIPT), 'Completed asset binding changed')
    print(json.dumps(dict(status='PASS', mode=outcome,
                          accepted_rounds=result['accepted_rounds'],
                          science_check_groups=result['science_check_groups'],
                          markdown_files=result['markdown_files_checked'],
                          local_links=result['local_links_checked'],
                          bound_assets=len(result['sha256'])), ensure_ascii=False))


if __name__ == '__main__':
    main()
