"""Administrative startup checks; no new scientific round or result."""
from pathlib import Path
import hashlib
import json
import re
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OUT = HERE / 'startup_checks.json'
NAVCHECK = BASE / '_migration/thermal_geometry_active_1219/verification.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def links(paths):
    count = 0
    for p in paths:
        text = p.read_text(encoding='utf-8-sig')
        text = re.sub(r'(?ms)^```.*?^```[^\n]*', '', text)
        text = re.sub(r'(?s)\$\$.*?\$\$', '', text)
        parent = BASE if p.name == 'navigation_notice.md' else p.parent
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', text):
            target = unquote(target.strip('<> ').split('#', 1)[0])
            if not target or re.match(r'^[a-zA-Z]+:', target):
                continue
            destination = (parent / target).resolve()
            assert destination.exists() or destination == OUT, (p, target)
            count += 1
    return count


def build():
    prompt = (HERE / 'goal_prompt.txt').read_text(encoding='utf-8').strip()
    receipt = read(HERE / 'goal_start_receipt.json')
    ledger = read(HERE / 'dependency_ledger.json')
    assert receipt['goal']['objective'] == prompt and receipt['goal']['status'] == 'active'
    objective_sha = hashlib.sha256(prompt.encode('utf-8')).hexdigest()
    assert ledger['goal_objective_sha256'] == objective_sha
    assert ledger['status'] == 'active' and ledger['first_round'] == 1219
    assert ledger['accepted_stage_rounds'] == [] and ledger['stage_check_groups'] == 0
    nav = read(NAVCHECK)
    assert nav['status'] == 'PASS'
    for name, expected in nav['file_sha256'].items():
        assert sha(ROOT / name) == expected, name
    assert all(item['exit_code'] == 0 for item in nav['checks'])
    assert nav['checks'][-1]['result']['old_GR_goal_completed'] is False
    assert nav['checks'][-1]['result']['goal_status'] == 'complete'
    policy = read(HERE / 'usage_reset_policy.json')
    assert policy['allowed_uses'] == 1 and policy['trigger_remaining_percent_strictly_less_than'] == 1
    assert not policy['redeemed']
    md = list(STAGE.rglob('*.md')) + [BASE / name for name in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')]
    count = links(md)
    artifacts = [p for p in STAGE.rglob('*') if p.is_file() and p != OUT]
    artifacts += [NAVCHECK, ROOT / 'scripts/run_research_thermal_geometry_active.py']
    return dict(schema='thermal_geometry_startup_checks_v1', date='2026-10-15', status='PASS',
        application_goal_status='active', goal_objective_sha256=objective_sha,
        first_round=1219, latest_historical_accepted_round=1218,
        accepted_scientific_rounds=[], new_scientific_rounds=0, stage_complete=False,
        phase='goal activation and first-question admission; no new theorem claimed',
        source_files_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(set(artifacts))},
        active_navigation_and_old_replay=nav['checks'], old_GR_goal_completed=False,
        four_laws_previous_goal_status='complete_in_declared_scope',
        usage_remaining_percent_at_last_tool_check=policy['latest_observed_remaining_percent'],
        one_conditional_reset_authorized=True, reset_used=False,
        local_markdown_files=len(md), local_links_checked=count, missing_links=[],
        known_historical_limit='Earlier full-layout formula inventory differs from later conjecture edits; already recorded in closure_1161_1218. Not rerun or silently rebound.')


def main():
    if '--write' in sys.argv:
        data = build()
        with OUT.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        data = read(OUT)
        for name, expected in data['source_files_sha256'].items():
            assert sha(ROOT / name) == expected, name
        assert data == build()
    print(json.dumps({k:data[k] for k in ('status', 'application_goal_status', 'first_round',
        'new_scientific_rounds', 'local_links_checked', 'stage_complete', 'reset_used')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
