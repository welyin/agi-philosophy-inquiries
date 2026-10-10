"""Verify frozen 1101 evidence, scope, links and exact certificates."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    a = json.loads((HERE/'acceptance.json').read_text(encoding='utf8'))
    assert a['status'] == 'PASS_CONDITIONAL_LORENTZ_SUBGROUP_AND_SCALE'
    assert a['conditional_full_Lorentz_subgroup']
    assert a['conditional_scale_elimination']
    assert not a['new_S_adopted'] and not a['S_derived_from_cognition']
    assert not a['R_O_FD_derived_from_cognition']
    assert not a['actual_clock_equals_cone_interval']
    assert not a['entire_conjecture_decided']
    assert a['new_adopted_axioms'] == a['scientific_count_increment'] == 0
    for name, digest in a['assets_sha256'].items():
        assert sha(STAGE/name) == digest, name
    for name, digest in a['prior_sources_sha256'].items():
        assert sha(ROOT/name) == digest, name
    reviews = json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))['reviews']
    assert len(reviews) == 2
    for review in reviews:
        assert review['status'] == 'PASS_CONDITIONAL_BRIDGE'
        assert review['actual_recomputation']['groups'] == 4
        assert review['delivery'] == 'read_only_collaboration_message'
        for name, digest in review['source_sha256'].items():
            assert sha(STAGE/name) == digest, name
    links = 0
    for name in a['assets_sha256']:
        p = STAGE/name
        if p.suffix != '.md':
            continue
        body = p.read_text(encoding='utf8')
        assert not re.search(r'\\[\[\]]', body)
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in body)
        assert len(re.findall(r'(?m)^\$\$\s*$', body)) % 2 == 0
        for target in re.findall(r'\]\(([^)]+)\)', body):
            if not target.startswith(('https:', 'http:', '#')):
                assert (p.parent/target.split('#', 1)[0]).exists(), target
                links += 1
    result = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE/'check.py')],
                            capture_output=True, text=True, encoding='utf8')
    assert result.returncode == 0, result.stdout+result.stderr
    check = json.loads(result.stdout)
    assert check['status'] == 'PASS' and check['groups'] == 4
    print(json.dumps(dict(round=1101, status='PASS', assets=len(a['assets_sha256']),
        prior_sources=len(a['prior_sources_sha256']), local_links=links,
        actual_recomputation=check, independent_reviews=2,
        entire_conjecture_decided=False, goal_status='active'), ensure_ascii=False))


if __name__ == '__main__':
    main()
