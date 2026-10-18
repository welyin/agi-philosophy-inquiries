"""Read-only integrity check for the scoped 1288--1289 stage acceptance.

This verifies saved evidence and local navigation; it is not a proof checker.
No scientific programs are rerun or results overwritten.  --write exclusively
creates the closure receipt after final acceptance has been written.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics/archive_1288_'
ACCEPTANCE = BASE / '_shared/stage_acceptance_after1289.json'
OUT = BASE / '_shared/stage_closure_checks_after1289.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def check():
    acceptance = load(ACCEPTANCE)
    assert acceptance['status'] == 'COMPLETE_WITH_SCOPE'
    assert acceptance['stage_complete'] is True
    assert acceptance['scientific_round_count'] == 2
    assert acceptance['author_check_groups'] == 2
    assert acceptance['independent_check_groups'] == 2
    assert acceptance['new_adopted_axioms'] == 0
    assert set(acceptance['required_deliverables']) == {
        'common_contract_and_dictionary', 'new_structural_result',
        'key_common_checks', 'reviewed_paper_and_remaining_freedom'}
    assert all(acceptance['required_deliverables'].values())
    assert sha(BASE / '_shared/goal.md') == acceptance['original_goal_sha256']
    current_bindings = 0
    for rel, digest in acceptance['files_sha256'].items():
        assert sha(ROOT / rel) == digest, rel
        current_bindings += 1
    round_bindings = 0
    source_bindings = 0
    for number in (1288, 1289):
        folder = BASE / str(number)
        record = folder / 'acceptance.json'
        assert sha(record) == acceptance['round_acceptance_sha256'][str(number)]
        data = load(record)
        assert data['status'] == 'ACCEPTED_WITH_SCOPE'
        for rel, digest in data['files_sha256'].items():
            assert sha(folder / rel) == digest, (number, rel)
            round_bindings += 1
        for name in ('results.json', 'independent_results.json'):
            data = load(folder / name)
            sources = data.get('source_hashes', data.get('source_sha256'))
            assert sources is not None, (number, name)
            for rel, digest in sources.items():
                assert sha(ROOT / rel) == digest, (number, name, rel)
                source_bindings += 1
    checked_links = 0
    checked_markdown = 0
    for rel in acceptance['files_sha256']:
        p = ROOT / rel
        if p.suffix != '.md':
            continue
        text = p.read_text(encoding='utf-8-sig')
        assert not any(x in text for x in ('\\[', '\\]', '\\(', '\\)')), rel
        assert sum(line.strip() == '$$' for line in text.splitlines()) % 2 == 0, rel
        for target in re.findall(r'\]\(([^\n]+?)\)', text):
            target = target.strip('<>')
            if re.match(r'^[a-zA-Z]+://', target) or target.startswith('#'):
                continue
            destination = (p.parent / unquote(target.split('#')[0])).resolve()
            assert destination.exists(), (rel, target)
            checked_links += 1
        checked_markdown += 1
    return dict(status='PASS', scope='saved-evidence integrity and local links only',
                scientific_rounds=[1288, 1289], original_goal_unchanged=True,
                stage_acceptance_sha256=sha(ACCEPTANCE),
                current_bindings=current_bindings, round_bindings=round_bindings,
                result_source_bindings=source_bindings,
                checked_markdown=checked_markdown, checked_local_links=checked_links,
                proof_or_experimental_validation=False, image_checks=False,
                scientific_results_overwritten=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = check()
    if args.write:
        with OUT.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    elif OUT.exists():
        assert load(OUT) == result, 'saved closure receipt differs'
    print(json.dumps(result, ensure_ascii=False))
