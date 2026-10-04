"""Document integrity and local links only; no new scientific tests or round."""
import argparse
import hashlib
import json
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
REPORT = HERE / 'round755_drafts/cognitive_joint_candidate_working_report.md'
TARGET = HERE / 'round755_drafts/cognitive_joint_candidate_working_report_checks.json'


def verify():
    text = REPORT.read_text('utf8')
    local_links = list(core.link_parser()(text))
    for link in local_links:
        assert (REPORT.parent / link).resolve().exists(), link
    receipt = core.read(HERE / 'research_round_754_checks.json')
    assert receipt['all_reported_checks_passed']
    for key in ('new_file_hashes', 'preserved_draft_hashes'):
        for name, digest in receipt[key].items():
            assert core.digest(HERE / name) == digest, name
    return dict(
        date='2026-10-04', kind='unnumbered_working_report',
        latest_completed_scientific_round=754, new_scientific_tests=0,
        report=REPORT.relative_to(HERE).as_posix(),
        report_sha256=hashlib.sha256(REPORT.read_bytes()).hexdigest(),
        verification_script_sha256=core.digest(Path(__file__)),
        local_links_checked=len(local_links), broken_links=0,
        current_round_frozen_files_unchanged=True,
        all_historical_hashes_covered_separately_by='postcheck_round754.py',
        primary_cognitive_sources=[
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC3371582/',
                 kind='human_dyadic_perceptual_experiment',
                 scope='task_limited; competing information-sharing explanations retained'),
            dict(url='https://web.stanford.edu/~clark/1990s/Clark%2C%20H.H.%20_%20Brennan%2C%20S.E.%20_Grounding%20in%20communication_%201991.pdf',
                 kind='theoretical_chapter',
                 scope='purpose_relative_grounding_and_communication_costs; not energy conservation'),
            dict(url='https://pages.ucsd.edu/~c8johnson/COGS102B/Hutchins95.pdf',
                 kind='distributed_cognition_task_analysis',
                 scope='cockpit_memory; not a universal cosmological principle'),
        ],
        primary_sources_read_online=True, images_inspected=False,
        source_claims_manually_reviewed=True, new_empirical_cognitive_study=False,
        universal_cognitive_necessity_proven=False,
        joint_quantum_gravity_model_completed=False,
        application_goal_changed=False,
        archive_handoff='Preserve this report and receipt separately; inherit their hashes before freezing the next round. They do not increment numbered science totals.',
        all_document_checks_passed=True,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    else:
        assert core.read(TARGET) == result
    print(json.dumps({k: result[k] for k in (
        'kind', 'latest_completed_scientific_round', 'new_scientific_tests',
        'local_links_checked', 'all_document_checks_passed')}))
