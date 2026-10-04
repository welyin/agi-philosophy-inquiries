"""Read-only science and text audit for rounds 270 and 271; no images."""
import argparse
import json
import re
import verify_interaction_rounds as core

CONFIG = {270: ('recent_cell_matching_audit', 8),
          271: ('uniform_edge_growth_audit', 9)}
core.CONFIG.update(CONFIG)


def verify(number, pending=False):
    result = core.verify(number, pending)
    body = (core.HERE/f'research_note_{number}.md').read_text(encoding='utf-8')
    assert not re.search(r'(?<![\\A-Za-z])(qquad|quad)(?![A-Za-z])', body)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        target = core.HERE/f'research_round_{args.round}_checks.json'
        if target.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
