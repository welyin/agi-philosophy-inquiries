"""Reproduce exact source classifications and check round 1017 delivery."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import stress_source_commutant as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / 'research_note_1017.md'
OUT = HERE / 'research_round_1017_checks.json'
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_231_775 import links

OWN = ['stress_source_commutant.py', 'stress_source_commutant_results.json',
       'verify_round1017.py', 'review.md', 'input_dependency_update_v0_6.md', 'NEXT.md']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads((HERE / OWN[1]).read_text('utf8')))
    assert fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups'] == 1
    assert fresh['cumulative_test_groups'] == 3795
    assert fresh['new_cognitive_axioms'] == 0
    assert fresh['full_multigraviton_no_go_claimed'] is False
    assert fresh['gravity_generated'] is False
    assert fresh['target_dimension_or_matter_content_selected'] is False
    families = fresh['exact_polynomial_families']
    assert len(families) == 9
    assert [x['symmetric_commutant_dimension'] for x in families] == [2,1,1,2,1,3,2,2,1]
    assert all(x['all_polynomial_coefficients_checked_exactly'] for x in families)
    assert fresh['rotation_preserves_all_exact_commutant_dimensions']
    assert fresh['on_shell_jet_calibration']['samples'] == 30
    assert fresh['on_shell_jet_calibration']['maximum_divergence'] < 1e-12
    counter = fresh['mixed_vertex_counterexample']
    assert counter['source_span_rank'] == 2
    assert counter['projector_identities_exact']
    assert counter['different_metric_scalar_density_checks'] == 8
    assert counter['maximum_scalar_density_redefinition_residual'] < 2e-12
    assert counter['gravity_action_numerically_rederived'] is False
    bad = fresh['impossible_weighted_source']
    assert bad['exact_curl'] == {'1,1':'-4'}
    assert bad['path_x_then_y'] == '0' and bad['path_y_then_x'] == '4'
    assert bad['direct_on_shell_divergence'][1] == 4.0
    assert 'additive constant in each W' in fresh['retained_freedom']
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
    assert all(f'## {n}.' in note for n in range(1,11))
    review = (HERE / 'review.md').read_text('utf8')
    assert '正式报告交叉复核通过' in review
    assert '主代理审阅通过' in review
    assert '独立数学终审通过' in review
    return dict(round=1017, date='2026-10-08', all_delivery_checks_passed=True,
        scientific_result_reproduced=True, new_calibration_groups=1,
        cumulative_research_groups=3795, local_links_checked=count,
        frozen_current_files=len(assets), historical_input_files=len(fresh['historical_source_sha256']),
        live_navigation_frozen=False, neighboring_round_frozen=False,
        exact_polynomial_families=9, exact_rational_coefficient_nullspaces=True,
        on_shell_jet_calibrations=30, different_metric_density_checks=8,
        finite_samples_not_general_proof=True,
        target_Rn_used_for_global_separation=True,
        additive_vacuum_constants_retained=True,
        mixed_vertex_graph_counterexample_has_complete_classical_action=True,
        preselected_Einstein_action_used_only_for_counterexample=True,
        linear_source_rank_not_nonlinear_metric_uniqueness=True,
        independent_mathematical_review=True, independent_reproduction=True,
        full_multigraviton_no_go_claimed=False,
        gravity_generated=False, new_cognitive_axiom=False,
        goal_completed=False, visual_checks_performed=False,
        analytic_scope=fresh['scope'],
        historical_source_sha256=fresh['historical_source_sha256'],
        source_sha256={str(path.relative_to(ROOT)).replace('\\', '/'):sha(path) for path in assets})


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
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')}, ensure_ascii=False, indent=2))
