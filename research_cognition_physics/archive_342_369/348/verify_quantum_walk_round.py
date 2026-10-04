"""Verify round 348 and retain the frozen terminology correction and drafts."""
import argparse
import json
import verify_autonomous_causality_rounds as previous

CONFIG = {348: ('quantum_walk_relativistic_limit_audit', 17)}
core = previous.core
core.CONFIG.update(CONFIG)


def verify(number, pending=False):
    prior = core.read(core.HERE/'round345_347_integration_checks.json')
    extra = dict(prior['new_text_maintenance_hashes'])
    extra.update(prior['preserved_schema_draft_hashes'])
    for name, sha in extra.items():
        assert core.digest(core.HERE/name) == sha, name
    result = previous.verify(number, pending)
    result['parallel_batch'] = []
    result['execution_mode'] = 'successor of completed independent rounds'
    result['scientific_base_through_round'] = 344
    result['additional_frozen_dependency_rounds'] = [346]
    result['batch_scientific_dependencies'] = []
    result['additional_preserved_files_verified'] = extra
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        path = core.HERE/f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
