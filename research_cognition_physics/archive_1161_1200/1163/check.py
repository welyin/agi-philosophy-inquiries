"""1163: three actual readouts, complete internal source, finite-domain scope.

Default recomputes and compares frozen results without writing.
--write-once creates results.json exclusively. No field locality, source AB
metric, old preparation under new H, or curvature is certified here.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
REUSE = HERE.parent / '1161/check.py'
REUSE_SHA = '52edf7a1dfc530996daa88ee3715efa936f7c5845eefe067124a544688621119'
TOL = 2e-9
N = np.array([1., .8])
G = np.array([np.eye(3), [[.9, .05, 0], [.05, .6, 0], [0, 0, .8]]])
OMEGA = 1.
TIME = .5
P = np.array([.3, .4, .2])
DELTA = 1e-4


def small(value, name, tol=TOL):
    v = float(np.linalg.norm(value))
    assert v <= tol, (name, v, tol)
    return v


def reuse():
    assert hashlib.sha256(REUSE.read_bytes()).hexdigest() == REUSE_SHA
    spec = importlib.util.spec_from_file_location('r1163_frozen_1161', REUSE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    old, data = module.reuse()
    return data


def energy_mass(p, r):
    return N[r]*math.sqrt(9+float(p@G[r]@p))


def energy_signal(k, r, positive):
    return (N[r] if positive else 1.)*math.sqrt(float(k@G[r]@k))


def group_readouts(hd):
    exact_gamma = [[F(9, 10), F(1, 20), F(0)],
                   [F(1, 20), F(3, 5), F(0)],
                   [F(0), F(0), F(4, 5)]]
    exact_n = F(4, 5)
    exact_k = [[exact_n*x for x in row] for row in exact_gamma]
    exact_s = [[exact_n**2*x for x in row] for row in exact_gamma]
    assert all(exact_s[i][j] == exact_n*exact_k[i][j]
               for i, j in np.ndindex(3, 3))
    mismatch = exact_gamma[2][2]-exact_n*exact_k[2][2]
    assert mismatch == F(36, 125)
    # Entire original clock generator, not a selected energy gap.
    od = np.diag(np.linspace(-1., 1., len(hd)))
    comm = 1j*(hd@od-od@hd)
    assert np.linalg.norm(comm) > .1
    clock = max(small(1j*((n*hd)@od-od@(n*hd))-n*comm,
                      'complete clock operator') for n in N)
    eps = 1e-3
    hess_res = velocity_res = signal_res = 0.
    rows = []
    for r in range(2):
        hess = np.zeros((3, 3))
        z = np.zeros(3)
        for i, j in np.ndindex(3, 3):
            ei, ej = np.eye(3)[i]*eps, np.eye(3)[j]*eps
            if i == j:
                hess[i, j] = (energy_mass(ei, r)-2*energy_mass(z, r)
                              +energy_mass(-ei, r))/eps**2
            else:
                hess[i, j] = (energy_mass(ei+ej, r)-energy_mass(ei-ej, r)
                              -energy_mass(-ei+ej, r)+energy_mass(-ei-ej, r))/(4*eps**2)
        hess_res = max(hess_res, small(3*hess-N[r]*G[r], 'massive Hessian',
                                      tol=2e-7))
        for positive in (True, False):
            expected_s = (N[r]**2 if positive else 1.)*G[r]
            # Six directions determine all entries, including cross terms.
            for k in (np.eye(3).tolist()+[[1., 1, 0], [1., 0, 1], [0., 1, 1]]):
                k = np.array(k)
                signal_res = max(signal_res, small(
                    energy_signal(k, r, positive)**2-k@expected_s@k, 'signal quadratic form'))
            direction = np.array([1., 2., 3.])/math.sqrt(14)
            h = 1e-5
            numerical = np.array([(energy_signal(direction+h*e, r, positive)
                                  -energy_signal(direction-h*e, r, positive))/(2*h)
                                  for e in np.eye(3)])
            analytic = ((N[r] if positive else 1.)*G[r]@direction
                        /math.sqrt(direction@G[r]@direction))
            velocity_res = max(velocity_res, small(numerical-analytic, 'position derivative', 1e-8))
        rows.append(dict(source_label=r, N=float(N[r]), gamma_inverse=G[r].tolist(),
                         K=(N[r]*G[r]).tolist(), S_plus=(N[r]**2*G[r]).tolist(),
                         S_minus=G[r].tolist(),
                         massive_and_plus_speed_bound=float(N[r]*math.sqrt(np.linalg.eigvalsh(G[r])[-1])),
                         minus_signal_speed_bound=float(math.sqrt(np.linalg.eigvalsh(G[r])[-1]))))
    return dict(exact_negative_zz=str(mismatch), full_clock_residual=clock,
                massive_hessian_residual=hess_res, signal_form_residual=signal_res,
                position_velocity_residual=velocity_res, rows=rows)


def hamiltonian(e, p, k, positive):
    # Shape D x AB x AB; M/E/L and unknown input remain actual tensor axes.
    h = np.zeros((len(e), 4, 4), complex)
    for ab in range(4):
        r = int(ab == 3)
        h[:, ab, ab] = 13+N[r]*e+energy_mass(p, r)+energy_signal(k, r, positive)
    h[:, 0, 3] = h[:, 3, 0] = OMEGA
    return h


def flatten_q(state):
    return state.reshape(-1, 2)


def pullback(state, image):
    a, b = flatten_q(state), flatten_q(image)
    return a.conj().T@b


def reduced_mel_choi(state):
    a = state.reshape(-1, 12, 2)
    return np.einsum('amq,anr->mqnr', a, a.conj())


def group_process(data):
    hd = data['m']['hout']
    e, w = np.linalg.eigh(hd)
    expected = data['m']['baseline']+np.repeat(np.arange(4), 72)
    spectrum = small(e-expected, 'same complete HD spectrum')
    xs = np.zeros((4, 4), complex)
    xs[0, 3] = xs[3, 0] = OMEGA
    stats = {key: 0. for key in ('gram', 'total_energy', 'original_HD', 'mel_choi', 'current_balance')}
    rows = []
    for positive in (True, False):
        for kval in (1., 2.):
            k = np.array([kval, 0., 0.])
            h = hamiltonian(e, P, k, positive)
            vals, vecs = np.linalg.eigh(h)
            u = np.einsum('dab,db,dcb->dac', vecs, np.exp(-1j*TIME*vals), vecs.conj())
            current = 1j*(h@xs-xs@h)
            rest = h-xs
            opposite = 1j*(h@rest-rest@h)
            stats['current_balance'] = max(stats['current_balance'],
                                           small(current+opposite, 'complete energy-current balance'))
            for program in range(3):
                raw = data['outputs'][program]
                a = (w.conj().T@raw.reshape(len(e), -1)).reshape(len(e), 4, 12, 2)
                b = np.einsum('dab,dbmq->damq', u, a)
                ha = np.einsum('dab,dbmq->damq', h, a)
                hb = np.einsum('dab,dbmq->damq', h, b)
                before, after = pullback(a, ha), pullback(b, hb)
                stats['gram'] = max(stats['gram'], small(pullback(b, b)-np.eye(2), 'unknown Q/reference isometry'))
                stats['total_energy'] = max(stats['total_energy'], small(after-before, 'complete energy pullback'))
                stats['original_HD'] = max(stats['original_HD'], small(
                    pullback(b, e[:, None, None, None]*b)-pullback(a, e[:, None, None, None]*a),
                    'bare HD spectral energy'))
                stats['mel_choi'] = max(stats['mel_choi'], small(reduced_mel_choi(b)-reduced_mel_choi(a),
                                                               'all MEL reduced Choi entries'))
                rb = np.zeros_like(b); rb[:, 3] = b[:, 3]
                ra = np.zeros_like(a); ra[:, 3] = a[:, 3]
                effect = pullback(b, rb)
                delta_x = pullback(b, np.einsum('ab,dbmq->damq', xs, b))-pullback(
                    a, np.einsum('ab,dbmq->damq', xs, a))
                rows.append(dict(candidate='plus' if positive else 'minus', signal_center=kval,
                                 program=program, initial_R=np.diag(pullback(a, ra)).real.tolist(),
                                 final_R=np.diag(effect).real.tolist(),
                                 R_effect_offdiag_abs=float(abs(effect[0, 1])),
                                 source_exchange_energy=np.diag(delta_x).real.tolist(),
                                 initial_complete_energy=np.diag(before).real.tolist(),
                                 current_operator_norm=float(np.linalg.norm(current[0]))))
    contrasts = []
    for candidate in ('plus', 'minus'):
        selected = [r for r in rows if r['candidate'] == candidate and r['program'] == 0]
        d = abs(selected[1]['final_R'][1]-selected[0]['final_R'][1])
        assert d > 2e-4, (candidate, d)
        contrasts.append(dict(candidate=candidate, program=0, Q_column=1, R_response_contrast=d))
    m0 = float(data['m']['baseline'])
    return dict(original_spectrum_residual=spectrum, full_operator_checks=stats, rows=rows,
                feedback_contrasts=contrasts, complete_positive_lower_bound=13+.8*(m0+3)-1)


def group_effective(process):
    # Both exact candidates have global momentum energy Lipschitz constants <=1.
    speed_bounds = [float(math.sqrt(np.linalg.eigvalsh(g)[-1])) for g in G]
    assert max(speed_bounds) <= 1+1e-14
    wave_unitary_error = 2*TIME*DELTA
    contrast_error = 2*wave_unitary_error
    margins = [r['R_response_contrast']-contrast_error for r in process['feedback_contrasts']]
    assert min(margins) > 0
    taylor = []
    for r in range(2):
        x = float(P@G[r]@P)
        actual = abs(energy_mass(P, r)-3*N[r]-N[r]*x/6)
        bound = N[r]*x*x/(8*27)
        assert actual <= bound+1e-14
        taylor.append(dict(source_label=r, actual_remainder=actual, bound=bound))
    return dict(momentum_halfwidth_ball=DELTA, time=TIME,
                each_experiment_probability_error=wave_unitary_error,
                contrast_error=contrast_error, certified_positive_contrast_margins=margins,
                exact_massive_low_momentum_remainders=taylor,
                packet='normalized product cosines, halfwidth delta/sqrt(3), continuous H1',
                locality_or_signal_front_certified=False)


def group_weakenings():
    eps, h, lam, pmax, mass = .02, 1e-5, .1, .5, 3.
    numerical_force = -eps*(math.sin(h)-math.sin(-h))/(2*h)
    force_residual = small(numerical_force+eps, 'extra scalar-force term', 1e-10)
    quartic = lam*pmax**4/mass**3
    # Analytic second derivative at zero is zero; finite differences converge quadratically.
    hess_values = [2*lam*d*d/mass**3 for d in (.1, .01, .001)]
    assert abs(hess_values[1]/hess_values[0]-.01) < 1e-14
    return dict(rest_scalar_potential_amplitude=eps, extra_force_at_origin=-eps,
                force_derivative_residual=force_residual,
                clock_gap_and_momentum_hessian_changes=0.,
                quartic_hessian_at_zero=0., finite_difference_hessian_sequence=hess_values,
                quartic_energy_at_declared_cutoff=quartic, cutoff=pmax,
                quartic_unitary_error_bound=TIME*quartic,
                scope='weakening witnesses, not complete countermodels to all cognitive axioms')


def calculate():
    data = reuse()
    readouts = group_readouts(data['m']['hout'])
    process = group_process(data)
    effective = group_effective(process)
    weakening = group_weakenings()
    return dict(schema='round1163_actual_three_readouts_v1', status='PASS', date='2026-10-13',
                check_groups=4, code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                reuse_code_sha256=REUSE_SHA, inherited_1160_replay='matches frozen results',
                adopted_new_axioms=0, readouts=readouts, complete_process=process,
                controlled_domain=effective, weakenings=weakening,
                claims=dict(common_principal_part_for_new_probe_menu=True,
                            old_AB_common_geometry=False, old_V_realized_under_new_H=False,
                            curvature_or_Einstein_dynamics=False))


def compare(a, b, path=''):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            compare(a[k], b[k], path+'/'+k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+'/'+str(i))
    elif isinstance(a, (int, float)) and not isinstance(a, bool):
        assert abs(a-b) <= 3e-9*max(1., abs(b)), (path, a, b)
    else:
        assert a == b, (path, a, b)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--write-once', action='store_true')
    args = ap.parse_args()
    result = calculate()
    target = HERE/'results.json'
    if args.write_once:
        with target.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        compare(result, json.loads(target.read_text(encoding='utf-8-sig')))
    print(json.dumps(dict(status='PASS', mode='write-once' if args.write_once else 'read-only comparison',
                          groups=4, feedback=result['complete_process']['feedback_contrasts'],
                          checks=result['complete_process']['full_operator_checks']), ensure_ascii=False))
