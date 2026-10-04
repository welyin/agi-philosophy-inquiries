"""Verify round 403 and the additive round-399 typesetting correction."""
import argparse
import json
import verify_recoverable_region_round as previous

CONFIG = {403: ('passive_control_locality_audit', 8)}
BASES = {403: 402}
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def correction_checks():
    manifest = core.read(core.HERE/'round399_text_correction.json')
    original = core.HERE/manifest['original']
    corrected = core.HERE/manifest['corrected']
    assert core.digest(original) == manifest['original_sha256']
    assert core.digest(corrected) == manifest['corrected_sha256']
    raw = original.read_bytes()
    old = manifest['old_phrase'].encode('utf-8')
    new = manifest['new_phrase'].encode('utf-8')
    assert raw.count(old) == manifest['occurrences'] == 1
    assert corrected.read_bytes() == manifest['banner'].encode('utf-8')+raw.replace(old,new,1)
    assert manifest['only_banner_and_one_missing_plus_changed']
    assert not manifest['numerical_results_changed']
    assert core.text_checks(corrected)['display_formulas'] == 13
    original_checks = core.read(core.HERE/'research_round_399_checks.json')
    for name, sha in original_checks['new_file_hashes'].items():
        assert core.digest(core.HERE/name) == sha, name
    return {name: core.digest(core.HERE/name)
            for name in ('research_note_399_text_v2.md','round399_text_correction.json')}


def verify(number, pending=False):
    result = previous.verify(number,pending)
    extra = correction_checks()
    target = core.HERE/f'research_round_{number}_checks.json'
    if target.exists():
        assert core.read(target)['round399_text_maintenance_hashes'] == extra
    result.update(date='2026-09-23', parallel_batch=[403],
                  execution_mode='passive versus intervened dependence on frozen round-402 baseline',
                  scientific_base_through_round=402,
                  additional_frozen_dependency_rounds=[],batch_scientific_dependencies=[],
                  round399_text_maintenance_hashes=extra,
                  round399_original_science_preserved=True,
                  text_maintenance_is_not_a_scientific_test=True,
                  final_science_review='primary-agent analytic, executable and source review; no independent agent final review')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round',type=int,choices=CONFIG)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify(args.round,args.write_checks)
    if args.write_checks:
        with (core.HERE/f'research_round_{args.round}_checks.json').open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
