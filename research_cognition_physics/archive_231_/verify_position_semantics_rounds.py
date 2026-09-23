"""Verify independent position-semantics rounds from frozen baseline 377."""
import argparse
import json
import verify_spatial_dimension_rounds as previous

CONFIG = {
    378: ('frame_position_quotient_audit', 14),
    379: ('position_prediction_closure_audit', 14),
}
BASES = {378: 377, 379: 377}
BATCH = [378, 379]
core = previous.core
core.CONFIG.update(CONFIG)
loaded = previous
while loaded is not None:
    if hasattr(loaded, 'BASES'):
        loaded.BASES.update(BASES)
    loaded = getattr(loaded, 'previous', None)


def preserved_covector_draft():
    name = 'round379_drafts/research_note_379_pre_covector.md'
    original = (core.HERE / name).read_bytes()
    assert original.count(b'p(X_i)') == 2
    newline = b'\r\n' if b'\r\n' in original else b'\n'
    insertion = b'\\ell=p\\cdot dx+\\lambda\\cdot dz,\\qquad' + newline
    expected = original.replace(b'u_i=p(X_i),', insertion + b'u_i=p(X_i),')
    expected = expected.replace(b'p(X_i)', b'\\ell(X_i)')
    assert (core.HERE / 'research_note_379.md').read_bytes() == expected
    hashes = {name: core.digest(core.HERE / name)}
    report = core.HERE / 'research_round_379_checks.json'
    if report.exists():
        assert core.read(report)['preserved_round379_draft_hashes'] == hashes
    return hashes


def verify(number, pending=False):
    result = previous.verify(number, pending)
    result['date'] = '2026-09-23'
    result['parallel_batch'] = BATCH
    result['execution_mode'] = 'independent complete position-semantics rounds in parallel'
    result['scientific_base_through_round'] = BASES[number]
    result['additional_frozen_dependency_rounds'] = []
    result['batch_scientific_dependencies'] = []
    if number >= 379:
        result['preserved_round379_draft_hashes'] = preserved_covector_draft()
        result['covector_notation_correction_changed_science'] = False
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        path = core.HERE / f'research_round_{args.round}_checks.json'
        if path.exists():
            raise RuntimeError('Report already exists; use read-only mode.')
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
