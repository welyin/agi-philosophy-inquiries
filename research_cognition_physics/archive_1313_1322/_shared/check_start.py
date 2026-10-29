"""Validate new-goal publication and read-only access to frozen evidence."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
import sys
from urllib.parse import unquote

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_research_interaction_generation_active import InteractionGenerationLayout


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    store = InteractionGenerationLayout()
    navigation = store.verify_interaction_navigation()
    closure = store.verify_closure()
    receipt = json.loads((STAGE / '_shared/goal_start_receipt.json').read_text(encoding='utf8'))
    assert receipt['status'] == receipt['tool_result']['goal']['status'] == 'active'
    assert receipt['objective'] == receipt['tool_result']['goal']['objective']
    assert receipt['objective'] == (STAGE / '_shared/goal_prompt.txt').read_text(encoding='utf8').strip()
    review = (STAGE / '_shared/goal_scope_review.md').read_text(encoding='utf8')
    assert 'PASS_WITH_SCOPE' in review
    for name in ('goal.md', 'start_plan.md'):
        assert sha((STAGE / '_shared' / name).read_bytes()) in review
    publication = json.loads((STAGE / '_shared/navigation_start_receipt.json').read_text(encoding='utf8'))
    for name, expected in publication['source_sha256'].items():
        assert sha((ROOT / name).read_bytes()) == expected, name
    links = 0
    markdown = list(STAGE.rglob('*.md'))
    for path in markdown:
        text = path.read_text(encoding='utf8')
        assert not re.search(r'^\s*\\[\[\]]\s*$', text, flags=re.M), path
        assert len(re.findall(r'^\s*\$\$\s*$', text, flags=re.M)) % 2 == 0, path
        for address in re.findall(r'\]\(([^\n)]+)\)', text):
            address = unquote(address.strip().strip('<>').split('#')[0])
            if not address or '://' in address or address.startswith('mailto:'):
                continue
            assert (path.parent / address).resolve().exists(), (path, address)
            links += 1
    scripts = list(STAGE.rglob('*.py')) + [ROOT / 'scripts/run_research_interaction_generation_active.py']
    for path in scripts:
        ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path))
    assert not list(STAGE.glob('research_note_*.md')), 'Startup must not claim a new scientific round'
    result = dict(status='PASS', date='2026-10-18', scope='stage closure and new goal startup only',
                  navigation=navigation, previous_stage_closure=closure,
                  stage_goal_start_record='active', latest_accepted_scientific_round=1312,
                  new_scientific_rounds=0, new_axioms=0, next_scientific_round=1313,
                  goal_scope_review='PASS_WITH_SCOPE', markdown_files=len(markdown),
                  local_links=links, python_ast_files=len(scripts), formula_delimiters=True,
                  no_full_historical_repository_validation_claimed=True,
                  source_sha256={p.relative_to(ROOT).as_posix(): sha(p.read_bytes())
                                 for p in markdown + scripts})
    if args.write:
        with (STAGE / '_shared/start_checks.json').open('x', encoding='utf8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
