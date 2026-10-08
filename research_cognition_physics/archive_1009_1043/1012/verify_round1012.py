"""Recompute the 1012 finite calibration and verify only its fixed artifacts."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import lepton_mass_operator_classification as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
OUT = HERE / 'research_round_1012_checks.json'
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_231_775 import links

OWN = ('lepton_mass_operator_classification.py',
       'lepton_mass_operator_classification_results.json',
       'review.md', 'NEXT.md', 'verify_round1012.py')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(writing=False):
    fresh = science.run()
    saved = json.loads((HERE / OWN[1]).read_text('utf8'))
    science.compare(fresh, saved)
    for name, digest in saved['historical_source_sha256'].items():
        assert sha(BASE / name) == digest, name
    assert fresh['all_scientific_calibrations_passed'] is True
    assert fresh['operator_count_before_U1'] == 54
    assert fresh['new_cognitive_axioms'] == 0
    for flag in ('complete_general_EFT_basis_claimed', 'physical_spectrum_fit_claimed',
                 'physical_UV_completion_claimed'):
        assert fresh[flag] is False
    evidence = [HERE.parent / 'research_note_1012.md'] + [HERE / n for n in OWN]
    nlinks = 0
    for path in evidence:
        assert path.is_file(), path
        if path.suffix != '.md':
            continue
        for _, _, target, local in links(path.read_text('utf-8-sig')):
            linked = (path.parent / local.replace('\\', '/')).resolve()
            assert linked.exists() or (writing and linked == OUT), (path, target)
            nlinks += 1
    return dict(round=1012, date='2026-10-08', all_delivery_checks_passed=True,
        scientific_result_reproduced=True, frozen_current_files=len(evidence),
        historical_inputs_verified=len(saved['historical_source_sha256']), local_links_checked=nlinks,
        new_calibration_groups=1, cumulative_research_groups=3790,
        analytic_scope='unknown-charge direct lepton mass operators through dimension five; full mass Ward space; specified-spectrum charge selection and isospectral counterexample',
        completeness_from_sampling_claimed=False, primary_proof_and_independent_code_review=True,
        live_navigation_frozen=False, neighboring_work_frozen=False,
        goal_completed=False, new_cognitive_axiom=False, visual_checks_performed=False,
        source_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in evidence})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run(args.write)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        assert result == json.loads(OUT.read_text('utf8'))
    print(json.dumps({k:v for k,v in result.items() if k != 'source_sha256'},
                     ensure_ascii=False, indent=2))
