"""Read-only evidence and recomputation checks for the scoped round1092 result."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    a = json.loads((HERE/'acceptance.json').read_text(encoding='utf8'))
    assert a['status'] == 'PASS_CONDITIONAL_BRIDGE'
    assert not a['entire_conjecture_decided'] and not a['Lorentz_derived']
    assert a['new_adopted_axioms'] == 0
    for path, value in a['assets_sha256'].items():
        assert digest(STAGE/path) == value, path
    for path, value in a['prior_sources_sha256'].items():
        assert digest(ROOT/path) == value, path
    for r in a['independent_reviews']:
        assert r['status'] == 'PASS_CONDITIONAL_BRIDGE'
        assert r['delivery'] == 'read_only_collaboration_message'
        assert r['reviewer_files_created'] is False
        for path, value in r['source_sha256'].items():
            assert digest(STAGE/path) == value, path
    links = 0
    for name in a['assets_sha256']:
        p = STAGE/name
        if p.suffix != '.md':
            continue
        txt = p.read_text(encoding='utf8')
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in txt)
        assert not re.search(r'\\[\[\]]', txt)
        assert len(re.findall(r'(?m)^\$\$\s*$', txt)) % 2 == 0
        for target in re.findall(r'\]\(([^)]+)\)', txt):
            if not target.startswith(('http:', 'https:', '#')):
                assert (p.parent/target.split('#',1)[0]).exists(), (p,target)
                links += 1
    output = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE/'check.py')],
                            capture_output=True, text=True, encoding='utf8')
    assert output.returncode == 0, output.stdout+output.stderr
    recomputation = json.loads(output.stdout)
    assert recomputation['status'] == 'PASS'
    print(json.dumps({'round':1092, 'status':'PASS', 'frozen_assets':len(a['assets_sha256']),
          'prior_sources':len(a['prior_sources_sha256']), 'independent_read_only_reviews':2,
          'local_links':links, 'recomputation':recomputation,
          'entire_conjecture_decided':False, 'goal_status':'active'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
