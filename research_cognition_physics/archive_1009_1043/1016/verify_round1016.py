"""Round 1016 reproduction and delivery checks; current navigation stays live."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import scalar_gauge_vertex_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / 'research_note_1016.md'
OUT = HERE / 'research_round_1016_checks.json'
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_231_775 import links

OWN = ['scalar_gauge_vertex_selection.py', 'scalar_gauge_vertex_selection_results.json',
       'verify_round1016.py', 'review.md', 'input_dependency_update_v0_5.md', 'NEXT.md']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads((HERE / OWN[1]).read_text('utf8')))
    assert fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups'] == 1
    assert fresh['cumulative_test_groups'] == 3794
    assert fresh['new_cognitive_axioms'] == 0
    for key in ['complete_arbitrary_EFT_classification_claimed',
                'group_or_matter_content_generated', 'physical_data_fit_claimed',
                'pointwise_density_used_as_general_action_failure_proof']:
        assert fresh[key] is False, key
    assert len(fresh['independent_noether_jet_witnesses']) == 5
    for witness in fresh['independent_noether_jet_witnesses']:
        assert abs(witness['noether_euler_coefficient']) > 0.01
    assert len(fresh['valid_finite_jet_samples']) == 24
    assert fresh['arbitrary_independent_coefficient_samples'] == 6
    assert fresh['maximum_residuals']['noether'] < 1e-12
    assert fresh['maximum_residuals']['directional_difference'] < 2e-8
    periodic = fresh['periodic_action_checks']
    assert len(periodic) == 6
    assert abs(periodic[0]['integrated_density_variation']) < 1e-12
    assert all(item['integrated_euler_square'] > 1e-6 for item in periodic[1:])
    assert fresh['derivative_improvement_quotient']['nonzero_S_equals_T_permitted']
    for name, digest in fresh['historical_source_sha256'].items():
        assert sha(BASE / name) == digest, name

    assets = [NOTE] + [HERE / name for name in OWN]
    count = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == '.md':
            content = path.read_text('utf-8-sig')
            assert content.count('$$') % 2 == 0, path
            for _, _, target, local in links(content):
                resolved = (path.parent / local.replace('\\', '/')).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                count += 1
    note = NOTE.read_text('utf8')
    for n in range(1, 11):
        assert f'## {n}.' in note
    review = (HERE / 'review.md').read_text('utf8')
    assert '正式报告交叉复核通过' in review
    assert '主代理审阅通过' in review
    assert '独立数学终审通过' in review
    return dict(round=1016, date='2026-10-08', all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3794, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(fresh['historical_source_sha256']),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        independent_noether_obstructions=5, valid_finite_jet_samples=24,
        arbitrary_coefficient_samples=6, periodic_action_families=6,
        necessity_uses_action_modulo_total_derivatives=True,
        gauge_representation_not_assumed_before_necessity_proof=True,
        covariance_rewrite_used_only_after_vertex_constraints=True,
        finite_samples_not_general_proof=True,
        nonzero_equal_raw_improvements_permitted=True,
        independent_mathematical_review=True, independent_reproduction=True,
        complete_arbitrary_EFT_classification_claimed=False,
        group_matter_masses_generated=False, new_cognitive_axiom=False,
        goal_completed=False, visual_checks_performed=False,
        analytic_scope=fresh['scope'],
        historical_source_sha256=fresh['historical_source_sha256'],
        source_sha256={str(path.relative_to(ROOT)).replace('\\', '/'): sha(path) for path in assets})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = verify(prospective=args.write)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        assert result == json.loads(OUT.read_text('utf8'))
    print(json.dumps({k: v for k, v in result.items() if not k.endswith('sha256')}, ensure_ascii=False, indent=2))
