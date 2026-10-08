"""Local delivery checks; no claim of independent personnel review."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import ast
import hashlib
import json
import re
from urllib.parse import unquote
import discrete_spin_matter_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent/'research_note_1034.md'
OUT = HERE/'research_round_1034_checks.json'
OWN = ['discrete_spin_matter_selection.py', 'discrete_spin_matter_selection_results.json',
       'drafts/discrete_spin_selection_derivation.md', 'selection_audit.md',
       'input_dependency_update_v0_23.md', 'review.md', 'NEXT.md', 'verify_round1034.py']
LINK = re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def local_links(s):
    s = re.sub(r'\$\$.*?\$\$|^~~~.*?^~~~\s*$', '', s, flags=re.S | re.M)
    for m in LINK.finditer(s):
        v = m[1].strip()
        v = v[1:v.index('>')] if v.startswith('<') else v
        path = unquote(v.split('#', 1)[0])
        if path and not re.match(r'^[a-zA-Z]+:', path):
            yield v, path


def exact_block_elimination():
    """Fraction Gaussian elimination of full blocks, independent of numpy solve."""
    y0 = [[F(x, 20) for x in r] for r in ((2, 1, 0), (0, 3, 1), (1, 0, 4))]
    m0 = [[F(x) for x in r] for r in ((4, 1, 0), (1, 5, 1), (0, 1, 6))]
    outputs = []
    # Order N first, LH second, eliminate N by elementary row operations.
    for t in (F(1, 2), F(1), F(2), F(3)):
        y = [[t*x for x in r] for r in y0]
        m = [[t*t*x for x in r] for r in m0]
        k = [m[i]+[y[j][i] for j in range(3)] for i in range(3)]
        k += [y[i]+[F(0)]*3 for i in range(3)]
        for j in range(3):
            pivot = k[j][j]
            assert pivot != 0
            k[j] = [x/pivot for x in k[j]]
            for i in range(6):
                if i != j:
                    a = k[i][j]
                    k[i] = [x-a*z for x, z in zip(k[i], k[j])]
        c = [r[3:] for r in k[3:]]
        assert all(c[i][j] == c[j][i] for i in range(3) for j in range(3))
        # x = M^{-1} Y^T, retained by the same row reduction.
        x = [r[3:] for r in k[:3]]
        assert all(sum(m[i][j]*x[j][a] for j in range(3)) == y[a][i]
                   for i in range(3) for a in range(3))
        assert all(c[a][b] == -sum(y[a][j]*x[j][b] for j in range(3))
                   for a in range(3) for b in range(3))
        outputs.append(c)
    assert all(c == outputs[0] for c in outputs)
    return [[str(v) for v in r] for r in outputs[0]]


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text('utf8')))
    assert fresh['code_sha256'] == sha(HERE/'discrete_spin_matter_selection.py')
    anomaly = fresh['anomalies']
    assert len(anomaly['odd_charge_periodicity']) == 66
    assert len(anomaly['menu_cases']) == 2646
    for r in anomaly['menu_cases']:
        assert r['audited_obstruction_cancels'] == ((r['plus']-r['minus']-r['generations']) % 16 == 0)
    assert anomaly['minimal_menus'] == [(3, 0)]
    assert anomaly['positive_only_solutions_to_64'] == [3, 19, 35, 51]
    # Direct rational SM charge traces; no charge-construction helper.
    d = [6, 3, 3, 2, 1]
    x = [1, 1, -3, -3, 1]
    assert sum(a*b for a, b in zip(d, x)) == -5
    assert sum(a*b**3 for a, b in zip(d, x)) == -125
    assert len(fresh['continuous_lift']) == 28
    for r in fresh['continuous_lift']:
        delta = r['N_with_X5']-r['generations']
        assert F(r['anomaly_traces']['X']) == 5*delta
        assert F(r['anomaly_traces']['X3']) == 125*delta
    words = fresh['mass_operators']['words']
    assert len(words) == 9
    phases = {'q': 1, 'uc': 1, 'dc': 1, 'l': 1, 'ec': 1, 'N': 1,
              'H': -1, 'Hbar': -1, 'S': -1, 'Nminus': -1}
    # Fermions carry i^{odd}; scalars carry +/-1. Direct complex characters.
    for r in words:
        phase = 1
        for field in r['fields']:
            phase *= phases[field] if field in ('H', 'Hbar', 'S') else 1j*phases[field]
        assert r['invariant'] == (phase == 1)
    matrices = fresh['mass_operators']['matrix_cases']
    assert len(matrices) == 6
    assert all(r['parent_covariance_residual'] < 1e-13 and r['Schur_residual'] < 1e-13
               and r['frozen_mass_symmetry_defect'] > 1 for r in matrices)
    freedoms = fresh['remaining_matching_freedom']
    assert len(freedoms) == 4
    for r in freedoms:
        t = F(r['scale'])
        assert F(r['C6_trace'])*t*t == F(949, 90000)
        assert [F(v) for v in r['C5_diagonal']] == [-F(1, 200), -F(1, 75), -F(9, 500)]
    rational_schur = exact_block_elimination()
    assets = [NOTE]+[HERE/name for name in OWN]
    links = 0
    for p in assets:
        assert p.is_file(), p
        if p.suffix == '.py':
            ast.parse(p.read_text('utf8'))
        if p.suffix == '.md':
            s = p.read_text('utf8')
            assert s.count('$$') % 2 == 0, p
            for raw, local in local_links(s):
                target = (p.parent/local).resolve()
                assert target.exists() or (prospective and target == OUT), (p, raw)
                links += 1
    for name, digest in fresh['historical_source_sha256'].items():
        assert sha(BASE/name) == digest
    assert all(f'## {i}.' in NOTE.read_text('utf8') for i in range(1, 11))
    assert '独立代理签审未完成' in (HERE/'review.md').read_text('utf8')
    assert fresh['bordism_classification_adopted_from_literature']
    assert not fresh['full_mixed_global_anomaly_classification_proved_here']
    assert not fresh['existing_P981_has_new_exact_Z4']
    assert not fresh['neutrino_masses_or_generation_number_predicted']
    assert fresh['new_adopted_cognitive_axioms'] == 0 and not fresh['goal_completed']
    return dict(round=1034, date='2026-10-08', local_delivery_checks_passed=True,
        scientific_result_reproduced=True, external_independent_agent_review_completed=False,
        new_calibration_groups=1, cumulative_research_groups=3811,
        new_adopted_cognitive_axioms=0, bordism_classification_adopted_from_literature=True,
        analytic_proof_not_grid_supplies_quantifiers=True,
        odd_charge_cases=66, menu_cases=2646, continuous_lift_cases=28,
        operator_words=9, parent_matrix_cases=6, rational_elimination_cases=4,
        rational_schur_block=rational_schur, physical_implementation_certified=False,
        full_mixed_global_anomaly_classification_completed=False,
        goal_completed=False, visual_checks_performed=False,
        frozen_current_files=len(assets), historical_input_files=len(science.HISTORY),
        local_links_checked=links, live_navigation_frozen=False,
        historical_source_sha256=fresh['historical_source_sha256'],
        source_sha256={str(p.relative_to(ROOT)).replace('\\', '/'):sha(p) for p in assets})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    result = verify(prospective=args.write_receipt)
    if args.write_receipt:
        with OUT.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write('\n')
    else:
        science.compare(result, json.loads(OUT.read_text('utf8')))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')}, ensure_ascii=False))
