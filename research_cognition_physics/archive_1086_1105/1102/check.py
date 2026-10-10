"""Reproduce inherited certificates for the declared axiom chain, not a new model."""
from pathlib import Path
import hashlib
import json
import re
import runpy
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compute():
    ledger = json.loads((HERE / 'dependency_ledger.json').read_text(encoding='utf8'))
    sources = ledger['sources_sha256']
    for path, expected in sources.items():
        assert digest(ROOT / path) == expected, path
    inherited = {}
    for round_number in (1100, 1101):
        directory = STAGE / str(round_number)
        namespace = runpy.run_path(str(directory / 'check.py'))
        result = namespace['compute']()
        saved = json.loads((directory / 'results.json').read_text(encoding='utf8'))
        assert result == saved, f'Inherited results changed: {round_number}'
        inherited[str(round_number)] = dict(
            status=result['status'], groups=list(result['groups']),
            saved_results_equal=True, source_code_sha256=digest(directory/'check.py'))
    ids = {item['id'] for item in ledger['candidate_extension_clauses']}
    assert ids == {'SR-T', 'SR-A', 'SR-R', 'SR-O', 'SR-F'}
    available = set(ids) | {'SR-S'}
    for claim in ledger['conditional_claims']:
        assert set(claim['depends_on']) <= available
        available.add(claim['id'])
    assert ledger['optional_extensions'] == ['SR-S']
    assert not ledger['full_joint_model_certified']
    assert not ledger['all_clock_matter_Lorentz_certified']
    links = 0
    for file in (HERE/'axioms.md', HERE/'NEXT.md', HERE/'cognitive_review.md', STAGE/'research_note_1102.md'):
        body = file.read_text(encoding='utf8')
        assert len(re.findall(r'(?m)^\$\$\s*$', body)) % 2 == 0
        assert '\\[' not in body and '\\]' not in body
        for target in re.findall(r'\]\(([^)]+)\)', body):
            if target.startswith(('https:', 'http:', '#')):
                continue
            resolved = file.parent/target.split('#', 1)[0]
            # The first run produces this round's result after these checks.
            if resolved.resolve() == (HERE/'results.json').resolve():
                continue
            assert resolved.exists(), (str(file), target)
            links += 1
    return dict(
        round=1102, status='PASS_CANDIDATE_REVIEW_AND_SOURCE_REPLAY',
        inherited_certificates=inherited, inherited_groups_recomputed=9,
        source_files_pinned=len(sources), local_links_checked=links,
        dependency_entries_resolve=True,
        dependency_graph_machine_checked_for_references_only=True,
        proof_validity_established_by='Read-only mathematical reviews and cited analytic proofs, not this schema check',
        candidate_extension_groups=5, newly_adopted_extension_groups=0,
        independence_of_five_groups_proven=False, cognitive_fit_is_not_an_experiment=True,
        SR_S_adopted=False, G_adopted=False,
        old_foundation_definitions_changed=False,
        full_joint_model_certified=False, all_clock_matter_Lorentz_certified=False,
        new_scientific_experiments=0, scientific_count_increment=0,
        scientific_count_total=3860, entire_goal_completed=False)


def main():
    data = compute()
    target = HERE/'results.json'
    if '--write' in sys.argv:
        assert not target.exists(), 'Refuse to overwrite a saved result.'
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    else:
        assert data == json.loads(target.read_text(encoding='utf8'))
    print(json.dumps(dict(round=1102, status=data['status'], inherited_groups=9,
                         new_scientific_experiments=0, joint_model_certified=False),
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
