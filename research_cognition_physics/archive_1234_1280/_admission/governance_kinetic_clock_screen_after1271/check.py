"""Read-only replay by default; --write creates results once.

An effective-model diagnostic, not a full cognitive/GR countermodel, not a
new device realization, and not a certificate for actual readout precision.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
E, KB = 1.602176634e-19, 1.380649e-23
STATES = [(0, 0), (0, 1), (1, 0), (1, 1)]


def matrix(rd):
    j, ts, rs = KB * 0.350, 0.070, 580000.0
    energy = np.array([j*(2*n-1)*(2*N-1)/2 for n, N in STATES])
    out, gen = [], np.zeros((4, 4))
    for i, (n, N) in enumerate(STATES):
        for label, resistance in [('L', rs), ('R', rs), ('D', rd)]:
            dst = STATES.index((n, 1-N) if label == 'D' else (1-n, N))
            de = energy[dst]-energy[i]
            rate = de / (E*E*resistance*math.expm1(de/(KB*ts)))
            out.append((i, dst, label, de, rate))
            gen[dst, i] += rate
            gen[i, i] -= rate
    return energy, out, gen


def evaluate():
    q0, q1 = F(290, 333), F(580, 623)
    gap = q1-q0
    # Positive Taylor terms prove exp(7)>1000, with no floating-point claim.
    exp7_lower = sum((F(7)**n/F(math.factorial(n)) for n in range(20)), F(0))
    assert exp7_lower > 1000
    w_lower = F(999, 1000)*gap
    error_per_history = F(1, 1000)
    observed_lower = w_lower - 2*error_per_history
    assert w_lower > F(3, 50) and observed_lower > F(29, 500)
    branches, matrices = [], []
    for rd, exact in [(43000.0, q0), (21500.0, q1)]:
        energy, channels, gen = matrix(rd)
        scale = max(-np.diag(gen))
        gibbs = np.exp(-energy/(KB*0.070)); gibbs /= sum(gibbs)
        assert max(abs(gen @ gibbs))/scale < 1e-14
        assert max(abs(gen.sum(axis=0)))/scale < 1e-14
        reverse = {(i, k, label): rate for i, k, label, de, rate in channels}
        ldb = max(abs(math.log(rate/reverse[k, i, label])+de/(KB*0.070))
                  for i, k, label, de, rate in channels)
        assert ldb < 1e-14
        q = [sum(r for i,k,label,de,r in channels if i==x and label=='D') /
             sum(r for i,k,label,de,r in channels if i==x) for x in range(4)]
        assert max(abs(x-float(exact)) for x in q) < 1e-14
        # Every positive state-dependent speed factor preserves the jump chain.
        for x, speed in enumerate([1/3, 2, 7/4, 3]):
            rates = [speed*r for i,k,label,de,r in channels if i==x]
            assert abs(rates[-1]/sum(rates)-q[x]) < 1e-14
        branches.append({'R_D_ohm':rd, 'D_first_probability_by_state':q,
                         'Gibbs_probability':gibbs.tolist(),
                         'stationary_residual_scaled':float(max(abs(gen@gibbs))/scale),
                         'ldb_residual':ldb})
        matrices.append(gen)
    # Select one actual model state and a finite common observation window.
    # Instrument access and error are NOT supplied by this computation.
    lam0, lam1 = -matrices[0][0,0], -matrices[1][0,0]
    horizon = 7/lam0
    y = -math.expm1(-lam1*horizon)
    p = [y*float(q1), y*(1-float(q1)), 1-y]
    w = (1-float(q0))*p[0]-float(q0)*p[1]
    assert w > float(w_lower)
    return {'status':'PASS_EFFECTIVE_MODEL_SCREEN_ONLY',
            'scope':'Same marked population model and local pure-clock hypothesis; not full metric equivalence.',
            'parameters':{'T_K':0.070,'J_over_kB_K':0.350,'V_V':0.0,'R_S_ohm':580000},
            'exact':{'q0':str(q0),'q1':str(q1),'difference':str(gap),
                     'exp7_Taylor_lower':str(exp7_lower),'finite_W_lower':str(w_lower),
                     'per_history_TV_contract':str(error_per_history),
                     'conditional_record_gap_lower':str(observed_lower)},
            'diagnostic':{'first_event_timeout_seconds':horizon,
                          'lambda0_per_second':float(lam0),'lambda1_per_second':float(lam1),
                          'changed_D_S_noevent':p, 'finite_W':w},
            'material_checks':branches,
            'not_certified':['instrument error','finite reservoir realization','full stress source',
                             'unchanged physical geometry','quantum unknown-input parent',
                             'Einstein or Jacobson selection'],
            'new_numbered_rounds':0,'new_adopted_axioms':0}


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write', action='store_true');args=parser.parse_args()
    result=evaluate(); target=HERE/'results.json'
    if args.write:
        with target.open('x', encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert json.loads(target.read_text(encoding='utf-8'))==result
    print(json.dumps({'status':result['status'],'gap':result['exact']['difference'],
                      'conditional_finite_bound':result['exact']['conditional_record_gap_lower']}))
