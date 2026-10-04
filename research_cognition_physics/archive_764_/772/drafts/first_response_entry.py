"""Reproduce the 772 entry bookkeeping and provenance, not a new experiment.

No original field equation, state, Ward identity, or continuum limit is tested
by this audit. The one-external-leg bound is inherited from round 770.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
TARGET = HERE / 'first_response_entry_results.json'


def vertex_menus(external_legs, loops):
    budget = 2 * loops + external_legs - 2
    # Every k-valent vertex costs k-2 >= 1, hence these bounds are exhaustive.
    arities = list(range(3, budget + 3))
    menus = []
    for counts in itertools.product(range(budget + 1), repeat=len(arities)):
        if sum((k - 2) * count for k, count in zip(arities, counts)) != budget:
            continue
        vertices = sum(counts)
        half_edges = sum(k * count for k, count in zip(arities, counts))
        assert (half_edges - external_legs) % 2 == 0
        internal = (half_edges - external_legs) // 2
        assert internal >= 0 and internal - vertices + 1 == loops
        menus.append(dict(vertices={str(k): count for k, count in zip(arities, counts) if count},
                          internal_lines=internal,
                          physical_hbar_power=internal + external_legs - vertices))
    return menus


def run():
    dependencies = [
        'research_note_601.md', 'research_note_625.md', 'research_note_732.md',
        'research_note_734.md', 'research_note_735.md', 'research_note_754.md',
        'research_note_756.md', 'research_note_763.md', 'research_note_765.md',
        'research_note_767.md', 'research_note_768.md', 'research_note_769.md',
        'research_note_770.md', 'research_note_771.md',
        'joint_local_ward_repair_results.json', 'research_round_771_checks.json',
        'round769_drafts/research_note_769_working.md', 'round772_drafts/STATUS.md',
    ]
    one = vertex_menus(1, 1)
    two = vertex_menus(2, 1)
    assert one == [dict(vertices={'3': 1}, internal_lines=1, physical_hbar_power=1)]
    assert two == [dict(vertices={'4': 1}, internal_lines=1, physical_hbar_power=2),
                   dict(vertices={'3': 2}, internal_lines=2, physical_hbar_power=2)]
    return dict(
        working_round=772, latest_completed_round=771,
        cumulative_scientific_tests=3499, new_scientific_tests=0,
        scope='Finite-order diagram bookkeeping and document provenance only; no new scientific result or numerical field simulation.',
        one_loop_mean_menu=one, one_loop_two_point_menu=two,
        one_leg_bound_reused_from_round_770=True,
        free_two_point_exception='The no-vertex propagator has physical hbar power one and is recorded separately, not forced into the interacting-vertex count.',
        counterterms_recorded_separately=['one-loop linear source contact', 'one-loop quadratic contact', 'background-split and pairing contacts'],
        source_to_initial_constraint_reused_from=[731, 754],
        fixed_source_Cauchy_response_reused_from=[601, 732, 768],
        source_on_original_background_sufficient_for_first_mean=True,
        physical_interacting_state_or_record_map_proven=False,
        finite_strength_self_consistency_proven=False,
        no_original_scientific_test_rerun=True,
        dependency_hashes={name: hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest()
                           for name in dependencies})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps({k: result[k] for k in ('working_round', 'latest_completed_round',
                                            'new_scientific_tests', 'one_loop_two_point_menu')}))
