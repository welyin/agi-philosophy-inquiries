"""Reproduce the neutrino mechanism document audit, not a physical solver."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'neutrino_adoption_results.json'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    receipt_path = STAGE / '1006/research_round_1006_checks.json'
    previous = read(receipt_path)
    assert previous['round'] == 1006
    assert previous['all_delivery_checks_passed']
    assert previous['fresh_test_groups'] == 0
    assert previous['cumulative_numbered_test_groups_from_1005'] == 3786
    assert previous['physical_unification_certified'] is False
    assert previous['full_goal_completed'] is False

    mechanism = HERE / 'neutrino_observation_adoption_v1.md'
    prose = mechanism.read_text(encoding='utf-8-sig')
    sections = re.split(r'^## 7\. P13四项义务\s*$', prose, flags=re.M)
    assert len(sections) == 2, 'Expected one final P13 obligations section'
    section = sections[1]
    assert not re.search(r'^## ', section, flags=re.M), 'P13 must be the last section'
    entries = re.findall(
        r'^([1-4])\. \*\*([^\n]+?)：\*\*\s+([^\n]+)$', section, re.M)
    assert tuple(number for number, _, _ in entries) == ('1', '2', '3', '4')
    keys = ('参与者', '机制', '资源、记录与反作用', '输入与边界')
    assert tuple(key for _, key, _ in entries) == keys
    assert all(value.strip() for _, _, value in entries)
    # Matching the declared obligations audits documentation, not physical truth.
    # No new oscillation calculation or source-detector simulation is performed.
    obligations = {key: value for _, key, value in entries}
    for term in ('无新物理数值试验', '无质量或混合参数预测',
                 '不新增认知公理', '成熟机制的条件性采用'):
        assert term in prose, term

    evidence = [Path(__file__), mechanism,
        STAGE / 'research_note_1007.md',
        HERE / 'drafts/adoption_decision.md',
        HERE / 'drafts/review_notes.md',
        STAGE / '1008/drafts/STATUS.md',
        receipt_path,
        STAGE.parent / 'archive_531_553/research_note_539.md',
        STAGE / 'research_note_995.md',
        STAGE / 'research_note_996.md',
        STAGE.parent / 'archive_702_741/research_note_730.md',
        STAGE / '993/common_candidate_v1.md',
        STAGE / '1005/overall_operation_hypothesis_v2_1.md']
    return dict(round=1007, date='2026-10-08', kind='mechanism_adoption_audit',
        all_audit_checks_passed=True, formal_reports=1007,
        fresh_physical_test_groups=0, cumulative_test_groups=3786,
        section='P13', evidence_obligations=obligations,
        evidence_obligation_count=len(obligations),
        audit_scope='document_structure_evidence_and_declared_adoption_limits',
        adoption_changes=['production_propagation_detection_chain_adopted',
            'field_state_conjugation_and_unitary_three_active_scope_declared',
            'mass_distinguishability_packet_overlap_and_classical_averaging_separated',
            'ordinary_matter_forward_potential_and_adiabatic_conversion_adopted',
            'source_records_resources_and_test_background_backreaction_limits_declared',
            'existing_mass_matching_and_optional_dynamic_W995_interfaces_reused'],
        new_cognitive_axioms=False,
        masses_predicted=False,
        PMNS_fit_performed=False,
        actual_source_detector_simulated=False,
        dynamic_W995_oscillation_certified=False,
        physical_unification_certified=False,
        full_goal_completed=False,
        visual_checks_performed=False,
        next_priority='overall_completeness_and_cross_sector_consistency_audit',
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
