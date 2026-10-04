"""Read-only science/text audit for rounds 264 and 265, including text errata."""
import argparse
import json
import re
import verify_interaction_rounds as core

CONFIG = {264: ('ancestral_quotient_audit', 10),
          265: ('coarse_diffusion_audit', 10)}
core.CONFIG.update(CONFIG)
BARE_SPACING = re.compile(r'(?<![\\A-Za-z])(qquad|quad)(?![A-Za-z])')


def corrections():
    manifest_path = core.HERE/'round261_263_text_corrections.json'
    manifest = core.read(manifest_path)
    hashes = {manifest_path.name: core.digest(manifest_path)}
    count = 0
    for entry in manifest['corrections']:
        old = core.HERE/entry['original']
        new = core.HERE/entry['corrected']
        assert core.digest(old) == entry['original_sha256']
        assert core.digest(new) == entry['corrected_sha256']
        text, number = BARE_SPACING.subn(lambda m: '\\'+m.group(), old.read_text(encoding='utf-8'))
        assert number == entry['recovered_backslashes']
        body = new.read_text(encoding='utf-8')
        assert body.split('\n\n', 1)[1] == text
        assert not BARE_SPACING.search(body)
        core.text_checks(new)
        hashes[new.name] = core.digest(new)
        count += number
    assert count == 13
    return hashes


def verify(number, pending=False):
    result = core.verify(number, pending)
    assert not BARE_SPACING.search((core.HERE/f'research_note_{number}.md').read_text(encoding='utf-8'))
    extra = corrections()
    if number == 264:
        result['new_file_hashes'].update(extra)
        result['text_correction_count'] = 13
        result['frozen_original_text_preserved'] = True
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        target = core.HERE/f'research_round_{args.round}_checks.json'
        if target.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
