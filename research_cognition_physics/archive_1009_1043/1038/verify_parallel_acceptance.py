"""Append-only integration receipt; no new scientific calibration group.

The scientific verifiers and independent reviews are separate evidence. This
script checks their frozen versions, provenance, links, and centralized counts.
Navigation is checked live and is deliberately not frozen for future rounds.
"""
from pathlib import Path
from urllib.parse import unquote
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
ROOT = BASE.parent
OUT = HERE / 'parallel_acceptance_1036_1038.json'
REVIEWS = {
    1035: STAGE/'1035/independent_review_by_1038.md',
    1036: STAGE/'1036/independent_review_by_1037.md',
    1037: STAGE/'1037/independent_review_by_1036.md',
    1038: HERE/'independent_review_by_1037.md',
}
NAV = [BASE/'README.md', BASE/'research_direction.md', BASE/'RESEARCH_STATE.md',
       STAGE/'README.md', STAGE/'文件索引.md']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def checked_links(p, prospective=False):
    s = p.read_text('utf8')
    s = re.sub(r'```.*?```|\$\$.*?\$\$', '', s, flags=re.S)
    count = 0
    for raw in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)', s):
        target = unquote(raw.strip().strip('<>').split('#', 1)[0])
        if not target or re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:', target):
            continue
        q = (p.parent/target).resolve()
        if prospective and q == OUT:
            continue
        assert q.exists(), (str(p), target)
        count += 1
    return count


def evidence():
    frozen, history, receipts = {}, {}, {}
    for n in range(1033, 1039):
        path = STAGE/str(n)/f'research_round_{n}_checks.json'
        d = json.loads(path.read_text('utf8'))
        receipts[str(path.relative_to(BASE))] = sha(path)
        if n == 1036:
            parent, mapping = STAGE, d['artifact_sha256']
            old = d['input_sha256']
        elif n == 1038:
            parent, mapping = BASE, d['frozen_file_sha256']
            old = d['historical_source_sha256']
        else:
            parent, mapping = ROOT, d['source_sha256']
            old = d['historical_source_sha256']
        for rel, digest in mapping.items():
            p = parent/rel
            assert sha(p) == digest, str(p)
            frozen[str(p.relative_to(BASE))] = digest
        for rel, digest in old.items():
            p = BASE/rel
            assert sha(p) == digest, str(p)
            history[str(p.relative_to(BASE))] = digest
    supplement = [HERE/'parallel_integration_1036_1038.md', Path(__file__),
                  HERE/'review_1037_checks.py', HERE/'review_1037_results.json']
    for n, p in REVIEWS.items():
        text = p.read_text('utf8')
        assert '通过' in text and str(n) in text, str(p)
        supplement.append(p)
    # Read statuses as evidence, never rewrite the original author's review flags.
    return dict(
        date='2026-10-08', rounds_accepted=[1036, 1037, 1038],
        cumulative_before=3812, new_calibration_groups=3,
        cumulative_by_round={'1036': 3813, '1037': 3814, '1038': 3815},
        cumulative_after=3815, new_adopted_cognitive_axioms=0,
        independent_agent_reviews={str(n): str(p.relative_to(BASE)) for n,p in REVIEWS.items()},
        human_peer_review=False, all_three_share_one_parent_model=False,
        full_physical_implementation_certified=False, goal_completed=False,
        historical_author_receipts_preserved=True,
        prior_rounds_still_awaiting_independent_agent_review=[1032, 1033, 1034],
        frozen_science_file_count=len(frozen), historical_input_file_count=len(history),
        frozen_science_sha256=frozen, historical_input_sha256=history,
        author_receipt_sha256=receipts,
        supplemental_review_sha256={str(p.relative_to(BASE)): sha(p) for p in supplement})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = evidence()
    live_links = sum(checked_links(p, args.write) for p in NAV)
    for p in [HERE/'parallel_integration_1036_1038.md', *REVIEWS.values()]:
        checked_links(p, args.write)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            json.dump(dict(result, navigation_links_at_acceptance=live_links), stream,
                      ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        old = json.loads(OUT.read_text('utf8'))
        old.pop('navigation_links_at_acceptance')
        assert old == result, 'Frozen integration evidence changed'
    print(json.dumps(dict(passed=True, accepted=result['rounds_accepted'],
                          cumulative=result['cumulative_after'], live_navigation_links=live_links,
                          frozen_science_files=result['frozen_science_file_count'],
                          historical_inputs=result['historical_input_file_count']), ensure_ascii=False))
