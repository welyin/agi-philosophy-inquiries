"""1033 local verification; independent human/agent review remains incomplete."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import ast
import hashlib
import json
import re
from urllib.parse import unquote
import numpy as np
import classical_channel_noise_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent/'research_note_1033.md'
OUT = HERE/'research_round_1033_checks.json'
OWN = ['classical_channel_noise_selection.py',
       'classical_channel_noise_selection_results.json',
       'drafts/classical_channel_noise_derivation.md', 'selection_audit.md',
       'input_dependency_update_v0_22.md', 'review.md', 'NEXT.md',
       'verify_round1033.py']
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


def independent_exponential_check():
    """Build Liouville matrices directly; no Schur closed form in this path."""
    z = np.diag([1., -1.]); a = np.kron(z, np.eye(2)); b = np.kron(np.eye(2), z)
    rho0 = np.ones((4, 4), complex)/4
    max_error = 0.
    for g, ga, gb in ((1, .25, .25), (1, 1, 1), (1, .5, 2),
                      (-1, 1, 2), (.3, .4, .3), (0, 0, .5)):
        h = g*a@b
        liouville = np.empty((16, 16), complex)
        for j in range(16):
            r = np.eye(16, dtype=complex)[:, j].reshape(4, 4)
            d = -1j*(h@r-r@h)
            d -= ga/2*(a@(a@r-r@a)-(a@r-r@a)@a)
            d -= gb/2*(b@(b@r-r@b)-(b@r-r@b)@b)
            liouville[:, j] = d.ravel()
        # Every term preserves trace; the chosen effective H is conserved.
        assert np.max(np.abs(np.eye(4).ravel()@liouville)) < 1e-14
        assert np.max(np.abs(h.T.ravel()@liouville)) < 1e-14
        for t in (.0001, .01, .2):
            series, term = np.eye(16, dtype=complex), np.eye(16, dtype=complex)
            for n in range(1, 121):
                term = term@(t*liouville)/n
                series += term
                if np.linalg.norm(term, np.inf) < 1e-16:
                    break
            else:
                raise AssertionError('Taylor convergence not reached')
            rho = (series@rho0.ravel()).reshape(4, 4)
            closed = science.multiplier(g, ga, gb, t)*rho0
            max_error = max(max_error, float(np.max(np.abs(rho-closed))))
            # Full Choi matrix from action on the sixteen matrix units.
            choi = np.zeros((16, 16), complex)
            for i in range(4):
                for j in range(4):
                    choi[4*i:4*i+4, 4*j:4*j+4] = series[:, 4*i+j].reshape(4, 4)/4
            assert np.max(np.abs(choi-choi.conj().T)) < 1e-13
            assert np.linalg.eigvalsh(choi)[0] > -1e-13
    assert max_error < 1e-13
    return max_error


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text('utf8')))
    assert fresh['code_sha256'] == sha(HERE/'classical_channel_noise_selection.py')
    assert len(fresh['grid']) == 80 and len(fresh['feedback']) == 6
    below = above = 0
    for row in fresh['grid']:
        det = F(str(row['gamma_A']))*F(str(row['gamma_B']))-F(str(row['g']))**2
        assert row['classical_threshold_satisfied'] == (det >= 0)
        below += det < 0; above += det >= 0
    assert below + above == 80
    for row in fresh['feedback']:
        assert row['generator_matrix_basis_residual'] < 1e-13
        assert len(row['steps']) == len(row['max_multiplier_error']) == 4
    w = fresh['finite_witness']
    expected = -F(3, 400)+F(9, 10000)/(2*F(97, 100))
    assert expected == F(-273, 38800) == F(w['exact_upper_bound'])
    assert F(w['robust_upper_bound']) == expected+F(1, 1000) < 0
    assert w['direct_witness'] <= float(expected)
    for key, sign in {'II': 1, 'XX': -1, 'YZ': -1, 'ZY': -1}.items():
        assert abs(w['local_pauli_coefficients'][key]-sign/4) < 1e-13
    z = np.diag([1., -1.]); x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    witness = (np.eye(4)-np.kron(x,x)-np.kron(y,z)-np.kron(z,y))/4
    assert np.max(np.abs(np.linalg.eigvalsh(witness)-np.array([-.5,.5,.5,.5]))) < 1e-13
    v = fresh['visibility']['finite_error_example']
    gap = F(989, 1000)**2-1/(1+F(39, 1000))
    assert F(v['exact_positive_gap_lower_bound']) == gap > 0
    assert not v['physical_measurement_data'] and not w['actual_physical_error_certified']
    independent_error = independent_exponential_check()
    assets = [NOTE]+[HERE/name for name in OWN]
    links = 0
    for p in assets:
        assert p.is_file(), p
        if p.suffix == '.py':
            ast.parse(p.read_text('utf8'))
        if p.suffix == '.md':
            s = p.read_text('utf8'); assert s.count('$$') % 2 == 0, p
            for raw, local in local_links(s):
                target = (p.parent/local).resolve()
                assert target.exists() or (prospective and target == OUT), (p, raw)
                links += 1
    for name, digest in fresh['historical_source_sha256'].items():
        assert sha(BASE/name) == digest
    assert all(f'## {i}.' in NOTE.read_text('utf8') for i in range(1,11))
    assert '独立代理签审未完成' in (HERE/'review.md').read_text('utf8')
    assert fresh['new_adopted_cognitive_axioms'] == 0 and not fresh['goal_completed']
    return dict(round=1033, date='2026-10-08', local_delivery_checks_passed=True,
        scientific_result_reproduced=True, external_independent_agent_review_completed=False,
        new_calibration_groups=1, cumulative_research_groups=3810,
        new_adopted_cognitive_axioms=0, analytic_proof_not_grid_supplies_quantifiers=True,
        rate_cases=80, below_threshold_cases=below, at_or_above_threshold_cases=above,
        finite_time_cases=240, feedback_constructions=6, finite_step_cases=24,
        independent_Liouville_exponential_cases=18,
        independent_exponential_residual=independent_error,
        actual_gravitational_channel_certified=False, physical_error_budget_certified=False,
        autonomous_implementation_certified=False, general_classical_gravity_excluded=False,
        goal_completed=False, visual_checks_performed=False,
        frozen_current_files=len(assets), historical_input_files=len(science.HISTORY),
        local_links_checked=links, live_navigation_frozen=False,
        historical_source_sha256=fresh['historical_source_sha256'],
        source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in assets})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args(); result = verify(prospective=args.write_receipt)
    if args.write_receipt:
        with OUT.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False); f.write('\n')
    else:
        saved = json.loads(OUT.read_text('utf8'))
        science.compare(result, saved)
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')}, ensure_ascii=False))
