"""Reuse frozen evidence to audit adoption scope; not a new physical experiment."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'boundary_scope_audit_results.json'


def read(p):
    return json.loads(p.read_text('utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def number(x):
    return {'exact': str(x), 'value': float(x)}


def run():
    late_path = STAGE / '997/vacuum_accessibility_results.json'
    early_path = STAGE / '996/cosmic_cp_window_results.json'
    heat_path = STAGE / '980/finite_thermal_records_results.json'
    late, early, heat = map(read, (late_path, early_path, heat_path))
    assert all(r['all_scientific_checks_passed'] for r in (late, early, heat))
    p = late['inputs']
    M, R, L = (F(p[k]['exact']) for k in
               ('material_mass', 'radiation_R', 'vacuum_cell_L'))
    assert M == L and min(M, R, L) > 0
    assert early['finite_window']['a_initial'] == 1
    assert early['finite_window']['a_final'] == 2
    rows = []
    for a in (F(1), F(2), F(4)):
        # Multiply all physical densities by the same V0*a**4.
        # This avoids floating point cancellation and leaves fractions invariant.
        weighted = (M*a, R, L*a**4)
        total = sum(weighted)
        matter, radiation, vacuum = (x/total for x in weighted)
        omitted_over_radiation = (M*a + L*a**4)/R
        assert matter + radiation + vacuum == 1
        assert radiation < F(21, 1000)
        assert omitted_over_radiation > 47
        rows.append({'a': str(a), 'matter_fraction': number(matter),
                     'radiation_fraction': number(radiation),
                     'vacuum_fraction': number(vacuum),
                     'nonradiation_over_radiation': number(omitted_over_radiation),
                     'H_full_over_H_radiation': math.sqrt(float(1/radiation))})
    assert all(rows[i]['radiation_fraction']['value'] >
               rows[i+1]['radiation_fraction']['value'] for i in range(2))
    # Existing thermal preparation is mixed. This is a reuse diagnostic only.
    populations = heat['parameters']['thermal_populations']
    assert abs(sum(populations)-1) < 1e-14 and min(populations) > 0
    distance_from_any_pure_lower = 1-max(populations)
    assert distance_from_any_pure_lower > .58
    evidence = [late_path, early_path, heat_path,
        STAGE/'research_note_959.md', STAGE/'research_note_973.md',
        STAGE/'993/common_candidate_v1.md',
        STAGE/'994/boundary_adoption_v1.md',
        HERE/'drafts/STATUS.md', HERE/'drafts/boundary_adoption_decision.md',
        Path(__file__)]
    return dict(kind='working_adoption_audit_of_existing_results',
        date='2026-10-07', all_audit_checks_passed=True,
        formal_reports=997, cumulative_numbered_test_groups=3781,
        new_scientific_test_groups=0,
        literal_same_window_diagnostic=rows,
        diagnosis_scope='Rejects literal identification of 997 calibration with '
            '996 radiation-dominated window; neither original result is refuted.',
        matter_radiation_equality_if_formally_extrapolated=number(R/M),
        early_extrapolation_of_fixed_bound_matter_certified=False,
        thermal_reference_distance_from_any_pure_lower=distance_from_any_pure_lower,
        joint_cosmological_history_certified=False,
        whole_program_refuted=False,
        app_goal_already_contains_priority_rule=True,
        app_goal_changed=False, full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)): sha(p) for p in evidence})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    out = run()
    if args.write:
        with TARGET.open('x', encoding='utf-8') as dest:
            json.dump(out, dest, ensure_ascii=False, indent=2)
            dest.write('\n')
    else:
        assert read(TARGET) == out, 'Saved audit differs from recomputation'
    print(json.dumps({k:v for k,v in out.items() if k != 'source_hashes'},
                     ensure_ascii=False, indent=2))
