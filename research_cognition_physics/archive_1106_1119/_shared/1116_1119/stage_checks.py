"""Publication checks for 1116-1119; no scientific experiment is re-executed."""
from pathlib import Path
import ast
import hashlib
import json
import re
from urllib.parse import unquote

STAGE = Path(__file__).resolve().parents[1]
BASE = STAGE.parent
PAPER = BASE / '认知动力学公理与狭义相对论基础动力学_阶段论文.md'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ledger = json.loads((STAGE / '_shared/dependency_ledger.json').read_text(encoding='utf-8'))
    assert ledger['completed_rounds'] == [1116, 1117, 1118, 1119]
    assert not ledger['in_progress_rounds']
    core = {}
    sources = {}
    groups = 0
    for n in ledger['completed_rounds']:
        folder = STAGE / str(n)
        accepted = json.loads((folder/'acceptance.json').read_text(encoding='utf-8'))
        assert accepted['status'] == 'PASS'
        review = folder / 'independent_review.md'
        assert sha(review) == accepted['review_sha256']
        review_text = review.read_text(encoding='utf-8')
        assert 'PASS' in review_text
        for name, expected in accepted['core_sha256'].items():
            p = (folder/name).resolve()
            assert sha(p) == expected, str(p)
            assert expected in review_text
            core[str(p.relative_to(BASE))] = expected
        groups += accepted['numerical_check_groups']
        ast.parse((folder/'check.py').read_text(encoding='utf-8-sig'))
        result = json.loads((folder/'results.json').read_text(encoding='utf-8'))
        for name, expected in result['source_sha256'].items():
            p = BASE/name
            assert sha(p) == expected, ('Changed source', name)
            sources[name] = expected
    assert groups == 22
    adoption_path = STAGE/'_shared/adoption_record.json'
    adoption = json.loads(adoption_path.read_text(encoding='utf-8'))
    assert adoption['status'] == 'ADOPTED' and adoption['axiom_group'] == 'D-P'
    assert ledger['new_adopted_axioms'] == ['D-P']
    for name, expected in adoption['evidence_sha256'].items():
        assert sha(BASE/name) == expected, ('Adoption evidence changed', name)
    final_scope = STAGE/'_shared/final_publication_scope_review.md'
    assert 'PASS' in final_scope.read_text(encoding='utf-8')
    assert sha(PAPER) in final_scope.read_text(encoding='utf-8')
    md_files = sorted(STAGE.rglob('*.md')) + [PAPER]
    local_links = 0
    for p in md_files:
        text = p.read_text(encoding='utf-8-sig')
        assert not any(ord(ch) < 32 and ch not in '\n\r\t' for ch in text), str(p)
        # All research math uses $$ display / $ inline; fenced code is excluded.
        prose = re.sub(r'(?ms)^(```|~~~).*?^\1[^\n]*$', '', text)
        assert not re.search(r'\\[\[\]()]', prose), ('Legacy math delimiter', str(p))
        display = [x for x in prose.splitlines() if x.strip() == '$$']
        assert len(display) % 2 == 0, ('Display delimiter', str(p))
        inline = re.sub(r'\$\$.*?\$\$', '', prose, flags=re.S)
        for i, line in enumerate(inline.splitlines(), 1):
            line = re.sub(r'`[^`]*`', '', line)
            assert len(re.findall(r'(?<!\\)\$', line)) % 2 == 0, ('Inline math', str(p), i)
        for target in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)', text):
            target = target.strip().strip('<>')
            if re.match(r'^(https?://|app://|codex://|#)', target):
                continue
            local = unquote(target.split('#', 1)[0])
            if not local:
                continue
            assert (p.parent/local).resolve().exists(), ('Broken link', str(p), target)
            local_links += 1
    result = dict(
        status='PASS', completed_rounds=ledger['completed_rounds'],
        numerical_check_groups=groups, scientific_calibration_increment=0,
        frozen_core_sha256=core, old_source_sha256=sources,
        markdown_files_checked=len(md_files), local_links_checked=local_links,
        paper_sha256=sha(PAPER), no_scientific_rerun=True,
        adoption_evidence_current=True, final_scope_sha256=sha(final_scope),
        goal_status=ledger['goal_status'], adopted_axioms=ledger['new_adopted_axioms'],
        new_app_task=False, old_files_modified=False,
    )
    out = STAGE/'_shared/stage_publication_checks.json'
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('_sha256')}, ensure_ascii=False))

if __name__ == '__main__':
    main()
