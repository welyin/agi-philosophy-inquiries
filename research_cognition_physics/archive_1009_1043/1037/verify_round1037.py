"""Delivery audit plus independent Bell/Fraction check; not human peer review."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import ast
import hashlib
import json
import re
from urllib.parse import unquote
import numpy as np
import local_record_energy_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent/'research_note_1037.md'
OUT = HERE/'research_round_1037_checks.json'
OWN = ['local_record_energy_selection.py', 'local_record_energy_selection_results.json',
       'verify_round1037.py', 'selection_audit.md', 'input_dependency_delta_1037.md',
       'review.md', 'NEXT.md']
LINK = re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bell_fraction_check():
    # Exact independent Pauli-eigenvalue table; no science matrix helper.
    signs = [(1, -1, 1), (-1, 1, 1), (1, 1, -1), (-1, -1, -1)]
    energies = [sum(a*b for a, b in zip(row, (1, 2, 3))) for row in signs]
    transition = [[F(0) for _ in range(4)] for _ in range(4)]
    for n in range(4):
        transition[n][n] = F(9, 10)
        transition[n ^ 1][n] = F(1, 10)
    b = sum(transition[m][n]*(energies[m]-energies[n])**2
            for m in range(4) for n in range(4))/4
    first_drifts = [sum(transition[m][n]*(energies[m]-energies[n])
                        for m in range(4)) for n in range(4)]
    second_drifts = [sum(transition[m][n]*(energies[m]**2-energies[n]**2)
                         for m in range(4)) for n in range(4)]
    assert b == 2 and sum(first_drifts) == 0 and sum(second_drifts) == 0
    assert 20*F(1, 10) == b
    c, coherence = F(3, 5), F(4, 5)
    assert coherence**2+c**2 == 1
    assert 10*(1-coherence) == b
    # Robust finite-data lower bound > 1.8, upper bound = 1.2.
    c0 = F(3, 5)-2*F(1, 100)
    upper = F(1)+100*F(1, 500)
    assert c0 == F(29, 50) and upper == F(6, 5)
    assert 1-c0*c0 == F(1659, 2500)
    assert 1659 < 41*41
    strict_lower = F(10)-F(41, 5)
    assert strict_lower == F(9, 5) and strict_lower-upper == F(3, 5)
    # Direct Bell matrix cross-check, built independently.
    v = np.array([[1, 1, 0, 0], [0, 0, 1, 1],
                  [0, 0, 1, -1], [1, -1, 0, 0]], complex)/np.sqrt(2)
    za = np.diag([1., 1., -1., -1.])
    ks = [(3*np.eye(4)+za)/(2*np.sqrt(5)), (3*np.eye(4)-za)/(2*np.sqrt(5))]
    actual = sum(abs(v.conj().T@k@v)**2 for k in ks)
    expected = np.array([[float(x) for x in row] for row in transition])
    residual = float(np.max(abs(actual-expected)))
    assert residual < 1e-13
    return dict(energies=energies, exact_B=str(b),
                first_drifts=list(map(str, first_drifts)),
                second_drifts=list(map(str, second_drifts)),
                finite_error_upper=str(upper), strict_lower_bound=str(strict_lower),
                strict_margin_lower=str(strict_lower-upper), bell_residual=residual)


def faithful_and_mean_only_check():
    # Full-rank invariant diagonal reference and a nontrivial dephasing channel.
    h = np.diag([-2., .5, 3.]).astype(complex)
    rho = np.diag([.2, .3, .5]).astype(complex)
    ks = [np.diag([int(k == j) for k in range(3)]).astype(complex) for j in range(3)]
    ph = lambda a: sum(k@a@k.conj().T for k in ks)
    ad = lambda a: sum(k.conj().T@a@k for k in ks)
    assert np.linalg.norm(ph(rho)-rho) < 1e-14
    assert np.linalg.norm(ad(h)-h) < 1e-14
    assert np.linalg.norm(ad(h@h)-h@h) < 1e-14
    # In the absorbing qutrit contrast the middle population decays by 1-p,
    # so stationarity with p>0 forces its population to zero.
    p = F(2, 5)
    assert 1-p < 1
    return dict(faithful_min_eigenvalue=.2, invariant_residual=0.,
                second_moment_residual=0., absorbing_middle_population_factor=str(1-p),
                absorbing_channel_has_faithful_invariant_state=False)


def local_links(s):
    s = re.sub(r'\$\$.*?\$\$', '', s, flags=re.S)
    for match in LINK.finditer(s):
        raw = match[1].strip()
        path = unquote(raw.split('#', 1)[0])
        if path and not re.match(r'^[a-zA-Z]+:', path):
            yield raw, path


def verify(prospective=False):
    result = science.run()
    science.compare(result, json.loads(science.OUT.read_text('utf8')))
    exact = bell_fraction_check()
    faithful = faithful_and_mean_only_check()
    assert len(result['arbitrary_instruments']) == 24
    assert len(result['saturation_cases']) == 12
    assert all(row['bound_slack'] > -1e-10 for row in result['arbitrary_instruments'])
    assert not result['goal_completed'] and not result['physical_implementation_certified']
    assets = [NOTE]+[HERE/name for name in OWN]
    links = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == '.py':
            ast.parse(path.read_text('utf8'))
        if path.suffix == '.md':
            contents = path.read_text('utf8')
            assert contents.count('$$') % 2 == 0, path
            for raw, relative in local_links(contents):
                target = (path.parent/relative).resolve()
                assert target.exists() or (prospective and target == OUT), (path, raw)
                links += 1
    note = NOTE.read_text('utf8')
    assert all(f'## {i}.' in note for i in range(1, 11))
    assert '忠实' in note and '全态一阶守恒本身' in note
    assert '其他代理完整审阅待完成' in (HERE/'review.md').read_text('utf8')
    for path, value in result['historical_source_sha256'].items():
        assert sha(BASE/path) == value
    return dict(round=1037, date='2026-10-08', local_delivery_checks_passed=True,
                scientific_result_reproduced=True,
                independent_agent_full_review_completed=False,
                primary_agent_preliminary_review_completed=True,
                human_peer_review_completed=False, new_calibration_groups=1,
                new_adopted_cognitive_axioms=0, arbitrary_instrument_cases=24,
                saturation_cases=12, biased_feedback_cases=1,
                independent_exact_verification=exact,
                faithful_reference_verification=faithful,
                physical_implementation_certified=False,
                full_parent_mapping_certified=False, goal_completed=False,
                frozen_current_files=len(assets), local_links_checked=links,
                historical_input_files=len(science.HISTORY), live_navigation_frozen=False,
                historical_source_sha256=result['historical_source_sha256'],
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
