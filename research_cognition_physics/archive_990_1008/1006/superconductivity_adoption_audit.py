"""Reproduce the superconductivity document audit, not a physical solver."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'superconductivity_adoption_results.json'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    receipt_path = STAGE / '1005/research_round_1005_checks.json'
    previous = read(receipt_path)
    assert previous['round'] == 1005
    assert previous['all_delivery_checks_passed']
    assert previous['fresh_test_groups'] == 0
    assert previous['cumulative_numbered_test_groups_from_1004'] == 3786
    assert previous['physical_unification_certified'] is False
    assert previous['full_goal_completed'] is False

    mechanism = HERE / 'conventional_superconductivity_adoption_v1.md'
    prose = mechanism.read_text(encoding='utf-8-sig')
    sections = re.split(r'^## P08补充说明\s*$', prose, flags=re.M)
    assert len(sections) == 2, 'Expected one final P08 supplement section'
    section = sections[1]
    assert not re.search(r'^## ', section, flags=re.M), 'P08 must be the last section'
    keys = ('参与者与过程', '中间机制', '资源、记录与反作用', '输入与范围')
    entries = re.findall(r'^- \*\*([^\n]+?)：\*\*\s+([^\n]+)$', section, re.M)
    assert tuple(key for key, _ in entries) == keys
    assert all(value.strip() for _, value in entries)
    # Formulae may remain symbolic. These are document evidence obligations;
    # no numerical value or physical truth is inferred from matching prose.
    obligations = dict(entries)
    for term in ('整体目标未完成', '新增物理试验0', '不新增认知公理'):
        assert term in prose, term

    whole = STAGE / '1005/overall_operation_hypothesis_v2_1.md'
    whole_key = str(whole.relative_to(ROOT))
    assert sha(whole) == previous['new_scientific_and_entry_files'][whole_key]
    evidence = [Path(__file__), mechanism,
        STAGE / 'research_note_1006.md',
        HERE / 'drafts/adoption_decision.md',
        HERE / 'drafts/review_notes.md',
        STAGE / '1007/drafts/STATUS.md',
        receipt_path, whole,
        STAGE.parent / 'archive_342_369/research_note_355.md',
        STAGE / '1000/material_phase_transport_adoption_v1.md',
        STAGE / '1004/resource_preparation_adoption_v1.md']
    return dict(round=1006, date='2026-10-08', kind='mechanism_adoption_audit',
        all_audit_checks_passed=True, formal_reports=1006,
        fresh_physical_test_groups=0, cumulative_test_groups=3786,
        section='P08', evidence_obligations=obligations,
        evidence_obligation_count=len(obligations),
        audit_scope='document_structure_evidence_and_declared_adoption_limits',
        adoption_changes=['conventional_pairing_gap_and_stiffness_distinguished',
            'gauge_invariant_electromagnetic_response_and_meissner_scope',
            'weak_link_relational_phase_with_finite_access_resources',
            'collective_order_dissipation_and_selective_records_distinguished'],
        new_cognitive_axioms=False,
        microscopic_pairing_derived_from_cognition=False,
        real_material_Tc_predicted=False,
        full_quantum_record_lifecycle_certified=False,
        physical_unification_certified=False,
        full_goal_completed=False,
        visual_checks_performed=False,
        next_priority='neutrino_production_propagation_detection',
        source_hashes={str(path.relative_to(ROOT)): sha(path) for path in evidence})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf-8') as dest:
            json.dump(result, dest, ensure_ascii=False, indent=2)
            dest.write('\n')
    else:
        assert result == read(TARGET)
    print(json.dumps({key: value for key, value in result.items()
        if key not in ('evidence_obligations', 'source_hashes')},
        ensure_ascii=False, indent=2))
