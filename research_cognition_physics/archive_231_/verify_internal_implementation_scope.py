"""Read-only by default; preserve unnumbered scope-review evidence separately."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import bilateral_receipt_scope_probe as probe
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'internal_implementation_scope_checks.json'
FILES = (
    'internal_implementation_scope_review.md',
    'bilateral_receipt_scope_probe.py',
    'bilateral_receipt_scope_probe_results.json',
    'macroscopic_space_conjecture_review.md',
    'verify_internal_implementation_scope.py',
)


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            compare(actual[key], saved[key])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    elif isinstance(actual, float):
        assert np.isclose(actual, saved, atol=1e-14, rtol=1e-11)
    else:
        assert actual == saved


def verify():
    saved = json.loads((HERE / FILES[2]).read_text(encoding='utf-8'))
    compare(probe.run(), saved)
    assert saved['numbered_round'] is None
    assert saved['numbered_scientific_test_count_increment'] == 0

    links = 0
    for name in (FILES[0], FILES[3]):
        content = (HERE / name).read_text(encoding='utf-8-sig')
        assert content.count('$$') % 2 == 0
        for link in re.findall(r'\]\(([^)]+)\)', content):
            if link.startswith(('https://', 'http://', '#')):
                continue
            assert (HERE / link.split('#')[0]).is_file(), (name, link)
            links += 1
    root_files = sorted(p.name for p in HERE.parent.iterdir() if p.is_file())
    assert len(root_files) == 5, root_files

    completed = subprocess.run(
        [sys.executable, '-B', '-X', 'utf8', str(HERE / 'verify_round503_integration.py')],
        cwd=HERE.parent.parent, capture_output=True, text=True,
        encoding='utf-8', timeout=120, check=True,
    )
    previous = json.loads(completed.stdout.strip().splitlines()[-1])
    assert previous['stage_saved_tests'] == 2450
    assert previous['total_protected_evidence_hashes'] == 892
    assert previous['all_reported_checks_passed']
    hashes = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in FILES}
    if TARGET.exists():
        frozen = json.loads(TARGET.read_text(encoding='utf-8'))
        assert hashes == frozen['additional_unnumbered_evidence_sha256']
    return dict(
        date='2026-09-27', latest_numbered_round=503,
        stage_numbered_scientific_tests=2450,
        numbered_science_files=820, prior_protected_evidence=892,
        new_numbered_scientific_tests=0, reproduced_diagnostic_checks=2,
        local_links_checked=links,
        additional_unnumbered_evidence_sha256=hashes,
        user_macroscopic_space_hypothesis_recorded=True,
        exact_or_limit_operation_semantics_selected=False,
        physical_three_dimensions_derived=False, full_gr_goal_completed=False,
        previous_scientific_evidence_unchanged=True,
        all_checks_passed=True,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in (
        'latest_numbered_round', 'stage_numbered_scientific_tests',
        'reproduced_diagnostic_checks', 'local_links_checked', 'all_checks_passed',
    )}, ensure_ascii=False))
